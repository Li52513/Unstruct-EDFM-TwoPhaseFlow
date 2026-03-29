import com.comsol.model.Model;
import com.comsol.model.NumericalFeature;
import com.comsol.model.util.ModelUtil;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

public final class ComsolHTCO2ConstPPSingleFracNoWell {
    private static final String DEFAULT_CASE_DIR =
        "Test/Transient/FullCaseTest/H_T_CO2_ConstPP/h_t_co2_constpp_singlefrac_nowell";
    private static final double MAIN_HMAX_M = 4.0;
    private static final double FINE_HMAX_M = 2.0;
    private static final double MAIN_HMIN_FRACTION = 0.25;
    private static final double FINE_HMIN_FRACTION = 0.125;
    private static final double MIN_HMIN_M = 1.0e-4;
    private static final double DP_NORM_MIN = 1.0;
    private static final double DT_NORM_MIN = 1.0;

    private ComsolHTCO2ConstPPSingleFracNoWell() {}

    public static void main(String[] args) {
        final CasePaths paths = resolvePaths(args);
        final boolean skipFineCheck = shouldSkipFineCheck(args);
        int exitCode = 0;

        try {
            ensureDirectories(paths);
            validateInputs(paths);

            final CaseConfig cfg = readCaseConfig(paths.propertyTableCsv);
            final List<ProfileStation> profileStations = readProfileStations(paths.profileStationsCsv);
            final List<MonitorPoint> monitorPoints = readMonitorPoints(paths.monitorPointsCsv);
            final List<TimeRequest> profileTimes = readProfileSchedule(paths.profileScheduleCsv);
            final List<MonitorSample> monitorSamples = readMonitorSchedule(paths.monitorScheduleCsv);
            final double[] studyTimes = buildStudyTimes(profileTimes, monitorSamples);

            ModelUtil.initStandalone(false);
            ModelUtil.loadPreferences();
            configureStandaloneRuntime(paths);
            ModelUtil.showProgress(paths.comsolOutputDir.resolve("comsol_progress.log").toString());

            final double mainHmin = Math.max(cfg.fractureWidthM * MAIN_HMIN_FRACTION, MIN_HMIN_M);
            final Model mainModel = buildAndSolveModel(paths, cfg, studyTimes, "ModelMain", "main", MAIN_HMAX_M, mainHmin);
            writeProfileReferenceCsvs(mainModel, profileStations, profileTimes, paths);
            writeMonitorReferenceCsv(mainModel, monitorPoints, monitorSamples, paths);
            writeRunSummary(paths, cfg, profileStations, monitorPoints, profileTimes, monitorSamples, MAIN_HMAX_M, mainHmin, skipFineCheck);
            mainModel.save(paths.comsolOutputDir.resolve("comsol_model.mph").toString());

            if (!skipFineCheck) {
                final double fineHmin = Math.max(cfg.fractureWidthM * FINE_HMIN_FRACTION, MIN_HMIN_M);
                final Model fineModel = buildAndSolveModel(paths, cfg, studyTimes, "ModelFine", "fine", FINE_HMAX_M, fineHmin);
                fineModel.save(paths.comsolOutputDir.resolve("comsol_model_refined.mph").toString());
                writeMeshCheck(mainModel, fineModel, profileStations, profileTimes, cfg, paths, MAIN_HMAX_M, mainHmin, FINE_HMAX_M, fineHmin);
            } else {
                writeSkippedFineCheck(paths);
            }

            ensureExpectedOutputs(paths, profileTimes);
            System.out.println("[COMSOL] H_T_CO2_ConstPP_SingleFrac_NoWell reference generation completed.");
            System.out.println("[COMSOL] Output directory: " + paths.comsolOutputDir);
        } catch (Exception ex) {
            exitCode = 1;
            ex.printStackTrace(System.err);
            try {
                ensureDirectories(paths);
                try (BufferedWriter out = Files.newBufferedWriter(
                    paths.comsolOutputDir.resolve("comsol_java_error.log"), StandardCharsets.UTF_8)) {
                    out.write(ex.toString());
                    out.newLine();
                    for (StackTraceElement frame : ex.getStackTrace()) {
                        out.write("  at " + frame.toString());
                        out.newLine();
                    }
                }
            } catch (IOException ignored) {
            }
        } finally {
            try {
                ModelUtil.clear();
            } catch (Exception ignored) {
            }
            try {
                ModelUtil.disconnect();
            } catch (Exception ignored) {
            }
            System.exit(exitCode);
        }
    }

    private static boolean shouldSkipFineCheck(String[] args) {
        if (args == null) {
            return false;
        }
        for (String arg : args) {
            final String trimmed = firstNonEmpty(arg);
            if ("--skip-fine-check".equals(trimmed)) {
                return true;
            }
        }
        return false;
    }

    private static CasePaths resolvePaths(String[] args) {
        String caseDirText = null;
        if (args != null && args.length > 0) {
            for (String arg : args) {
                final String trimmed = firstNonEmpty(arg);
                if (trimmed != null && !trimmed.startsWith("--")) {
                    caseDirText = trimmed;
                    break;
                }
            }
        }
        if (caseDirText == null) {
            caseDirText = DEFAULT_CASE_DIR;
        }

        final Path caseDir = Paths.get(caseDirText).normalize();

        final CasePaths paths = new CasePaths();
        paths.caseDir = caseDir;
        paths.engineeringDir = caseDir.resolve("engineering");
        paths.referenceDir = caseDir.resolve("reference");
        paths.comsolInputDir = paths.referenceDir.resolve("comsol_input");
        paths.comsolOutputDir = paths.referenceDir.resolve("comsol");
        paths.profileStationsCsv = paths.engineeringDir.resolve("profile_station_definitions.csv");
        paths.monitorPointsCsv = paths.engineeringDir.resolve("monitor_point_definitions.csv");
        paths.profileScheduleCsv = paths.engineeringDir.resolve("profile_report_schedule.csv");
        paths.monitorScheduleCsv = paths.engineeringDir.resolve("monitor_sample_schedule.csv");
        paths.propertyTableCsv = paths.comsolInputDir.resolve("property_table.csv");
        return paths;
    }

    private static void ensureDirectories(CasePaths paths) throws IOException {
        Files.createDirectories(paths.comsolOutputDir);
    }

    private static void configureStandaloneRuntime(CasePaths paths) throws IOException {
        final String runtimeRootOverride = firstNonEmpty(System.getProperty("codex.comsol.runtimeRoot"));
        final Path runtimeRoot = runtimeRootOverride != null
            ? Paths.get(runtimeRootOverride).toAbsolutePath().normalize()
            : paths.comsolOutputDir.resolve("comsol_runtime");
        final Path recoveryDir = runtimeRoot.resolve("recoveries");
        final Path tempDir = runtimeRoot.resolve("tmp");
        Files.createDirectories(recoveryDir);
        Files.createDirectories(tempDir);

        ModelUtil.setPreference("tempfiles.recovery.checkforrecoveries", "off");
        ModelUtil.setPreference("tempfiles.recovery.autosave", "off");
        ModelUtil.setPreference("tempfiles.recovery.savingeveryactive", "off");
        ModelUtil.setPreference("tempfiles.recovery.recoverydir", recoveryDir.toAbsolutePath().normalize().toString());
        ModelUtil.setPreference("tempfiles.tempfiles.tmpdir", tempDir.toAbsolutePath().normalize().toString());
        ModelUtil.savePreferences();
        System.out.println("[COMSOL] runtime_root=" + runtimeRoot.toString());
        System.out.println("[COMSOL] recovery_dir=" + recoveryDir.toString());
        System.out.println("[COMSOL] tmp_dir=" + tempDir.toString());
        System.out.println("[COMSOL] cs.recoverydir=" + String.valueOf(System.getProperty("cs.recoverydir")));
        System.out.println("[COMSOL] cs.prefsdir=" + String.valueOf(System.getProperty("cs.prefsdir")));
        System.out.println("[COMSOL] cs.tmpdir=" + String.valueOf(System.getProperty("cs.tmpdir")));
    }

    private static void validateInputs(CasePaths paths) {
        final List<Path> required = Arrays.asList(
            paths.profileStationsCsv,
            paths.monitorPointsCsv,
            paths.profileScheduleCsv,
            paths.monitorScheduleCsv,
            paths.propertyTableCsv
        );
        for (Path requiredPath : required) {
            if (!Files.isRegularFile(requiredPath)) {
                throw new IllegalStateException("Missing required input: " + requiredPath);
            }
        }
    }

    private static CaseConfig readCaseConfig(Path csvPath) {
        final CsvTable table = readCsv(csvPath);
        final Map<String, Double> values = new LinkedHashMap<String, Double>();
        for (String[] row : table.rows) {
            values.put(parseRequiredString(table, row, "key"), parseRequiredDouble(table, row, "value"));
        }

        final CaseConfig cfg = new CaseConfig();
        cfg.lxM = requireValue(values, "lx");
        cfg.lyM = requireValue(values, "ly");
        cfg.fracX0Ratio = requireValue(values, "frac_x0_ratio");
        cfg.fracY0Ratio = requireValue(values, "frac_y0_ratio");
        cfg.fracX1Ratio = requireValue(values, "frac_x1_ratio");
        cfg.fracY1Ratio = requireValue(values, "frac_y1_ratio");
        cfg.pInitPa = requireValue(values, "p_init");
        cfg.pLeftPa = requireValue(values, "p_left");
        cfg.pRightPa = requireValue(values, "p_right");
        cfg.tInitK = requireValue(values, "t_init");
        cfg.tLeftK = requireValue(values, "t_left");
        cfg.tRightK = requireValue(values, "t_right");
        cfg.matrixPhi = requireValue(values, "matrix_phi");
        cfg.matrixPerm = requireValue(values, "matrix_perm");
        cfg.matrixCt = requireValue(values, "matrix_ct");
        cfg.matrixRhoR = requireValue(values, "matrix_rho_r");
        cfg.matrixCpR = requireValue(values, "matrix_cp_r");
        cfg.matrixLambdaR = requireValue(values, "matrix_lambda_r");
        cfg.fracturePhi = requireValue(values, "fracture_phi");
        cfg.fractureKt = requireValue(values, "fracture_kt");
        cfg.fractureKn = requireValue(values, "fracture_kn");
        cfg.fractureCt = requireValue(values, "fracture_ct");
        cfg.fractureRhoR = requireValue(values, "fracture_rho_r");
        cfg.fractureCpR = requireValue(values, "fracture_cp_r");
        cfg.fractureLambdaR = requireValue(values, "fracture_lambda_r");
        cfg.co2Rho = requireValue(values, "co2_rho_const");
        cfg.co2Mu = requireValue(values, "co2_mu_const");
        cfg.co2Cp = requireValue(values, "co2_cp_const");
        cfg.co2K = requireValue(values, "co2_k_const");
        cfg.fractureWidthM = requireValue(values, "comsol_thin_band_width");
        cfg.deltaPPa = Math.max(Math.abs(cfg.pLeftPa - cfg.pRightPa), DP_NORM_MIN);
        cfg.deltaTK = Math.max(Math.abs(cfg.tLeftK - cfg.tRightK), DT_NORM_MIN);
        return cfg;
    }

    private static double requireValue(Map<String, Double> values, String key) {
        final Double value = values.get(key);
        if (value == null) {
            throw new IllegalStateException("Missing property-table entry: " + key);
        }
        return value.doubleValue();
    }

    private static Model buildAndSolveModel(
        CasePaths paths,
        CaseConfig cfg,
        double[] studyTimes,
        String modelTag,
        String labelSuffix,
        double hmaxM,
        double hminM
    ) {
        final Model model = ModelUtil.create(modelTag);
        model.modelPath(paths.comsolOutputDir.toString());
        model.label("ComsolHTCO2ConstPPSingleFracNoWell_" + labelSuffix + ".mph");

        configureParameters(model, cfg);
        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 2);
        model.component("comp1").mesh().create("mesh1");

        buildGeometry(model);
        createBoundarySelections(model, cfg);
        createVariables(model);
        createPressurePde(model);
        createTemperaturePde(model);
        configureMesh(model, hmaxM, hminM);
        configureStudy(model, studyTimes);

        model.study("std1").run();
        return model;
    }

    private static void configureParameters(Model model, CaseConfig cfg) {
        model.param().set("Lx", formatQuantity(cfg.lxM, "m"));
        model.param().set("Ly", formatQuantity(cfg.lyM, "m"));
        model.param().set("frac_x0_ratio", formatDouble(cfg.fracX0Ratio));
        model.param().set("frac_y0_ratio", formatDouble(cfg.fracY0Ratio));
        model.param().set("frac_x1_ratio", formatDouble(cfg.fracX1Ratio));
        model.param().set("frac_y1_ratio", formatDouble(cfg.fracY1Ratio));
        model.param().set("x0f", "frac_x0_ratio*Lx");
        model.param().set("y0f", "frac_y0_ratio*Ly");
        model.param().set("x1f", "frac_x1_ratio*Lx");
        model.param().set("y1f", "frac_y1_ratio*Ly");
        model.param().set("pInit", formatQuantity(cfg.pInitPa, "Pa"));
        model.param().set("pLeft", formatQuantity(cfg.pLeftPa, "Pa"));
        model.param().set("pRight", formatQuantity(cfg.pRightPa, "Pa"));
        model.param().set("TInit", formatQuantity(cfg.tInitK, "K"));
        model.param().set("TLeft", formatQuantity(cfg.tLeftK, "K"));
        model.param().set("TRight", formatQuantity(cfg.tRightK, "K"));
        model.param().set("matrix_phi", formatDouble(cfg.matrixPhi));
        model.param().set("matrix_perm", formatQuantity(cfg.matrixPerm, "m^2"));
        model.param().set("matrix_ct", formatQuantity(cfg.matrixCt, "1/Pa"));
        model.param().set("matrix_rho_r", formatQuantity(cfg.matrixRhoR, "kg/m^3"));
        model.param().set("matrix_cp_r", formatQuantity(cfg.matrixCpR, "J/(kg*K)"));
        model.param().set("matrix_lambda_r", formatQuantity(cfg.matrixLambdaR, "W/(m*K)"));
        model.param().set("fracture_phi", formatDouble(cfg.fracturePhi));
        model.param().set("fracture_kt", formatQuantity(cfg.fractureKt, "m^2"));
        model.param().set("fracture_kn", formatQuantity(cfg.fractureKn, "m^2"));
        model.param().set("fracture_ct", formatQuantity(cfg.fractureCt, "1/Pa"));
        model.param().set("fracture_rho_r", formatQuantity(cfg.fractureRhoR, "kg/m^3"));
        model.param().set("fracture_cp_r", formatQuantity(cfg.fractureCpR, "J/(kg*K)"));
        model.param().set("fracture_lambda_r", formatQuantity(cfg.fractureLambdaR, "W/(m*K)"));
        model.param().set("co2_rho_const", formatQuantity(cfg.co2Rho, "kg/m^3"));
        model.param().set("co2_mu_const", formatQuantity(cfg.co2Mu, "Pa*s"));
        model.param().set("co2_cp_const", formatQuantity(cfg.co2Cp, "J/(kg*K)"));
        model.param().set("co2_k_const", formatQuantity(cfg.co2K, "W/(m*K)"));
        model.param().set("wFrac", formatQuantity(cfg.fractureWidthM, "m"));
    }

    private static void buildGeometry(Model model) {
        model.component("comp1").geom("geom1").create("r1", "Rectangle");
        model.component("comp1").geom("geom1").feature("r1").set("base", "corner");
        model.component("comp1").geom("geom1").feature("r1").set("pos", new String[]{"0", "0"});
        model.component("comp1").geom("geom1").feature("r1").set("size", new String[]{"Lx", "Ly"});
        model.component("comp1").geom("geom1").run();
    }

    private static void createBoundarySelections(Model model, CaseConfig cfg) {
        final double tol = Math.max(1.0e-9, 1.0e-6 * Math.max(cfg.lxM, cfg.lyM));

        model.component("comp1").selection().create("leftBnd", "Box");
        model.component("comp1").selection("leftBnd").set("entitydim", "1");
        model.component("comp1").selection("leftBnd").set("condition", "inside");
        model.component("comp1").selection("leftBnd").set("xmin", formatDouble(-tol));
        model.component("comp1").selection("leftBnd").set("xmax", formatDouble(tol));
        model.component("comp1").selection("leftBnd").set("ymin", formatDouble(-tol));
        model.component("comp1").selection("leftBnd").set("ymax", formatDouble(cfg.lyM + tol));
        model.component("comp1").selection("leftBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("leftBnd").set("zmax", formatDouble(tol));

        model.component("comp1").selection().create("rightBnd", "Box");
        model.component("comp1").selection("rightBnd").set("entitydim", "1");
        model.component("comp1").selection("rightBnd").set("condition", "inside");
        model.component("comp1").selection("rightBnd").set("xmin", formatDouble(cfg.lxM - tol));
        model.component("comp1").selection("rightBnd").set("xmax", formatDouble(cfg.lxM + tol));
        model.component("comp1").selection("rightBnd").set("ymin", formatDouble(-tol));
        model.component("comp1").selection("rightBnd").set("ymax", formatDouble(cfg.lyM + tol));
        model.component("comp1").selection("rightBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("rightBnd").set("zmax", formatDouble(tol));
    }

    private static void createVariables(Model model) {
        model.component("comp1").variable().create("var1");
        model.component("comp1").variable("var1").set("Lf", "sqrt((x1f-x0f)^2+(y1f-y0f)^2)");
        model.component("comp1").variable("var1").set("txf", "(x1f-x0f)/Lf");
        model.component("comp1").variable("var1").set("tyf", "(y1f-y0f)/Lf");
        model.component("comp1").variable("var1").set("nxf", "-tyf");
        model.component("comp1").variable("var1").set("nyf", "txf");
        model.component("comp1").variable("var1").set("scoord", "(x-x0f)*txf + (y-y0f)*tyf");
        model.component("comp1").variable("var1").set("ndist", "(x-x0f)*nxf + (y-y0f)*nyf");
        model.component("comp1").variable("var1").set("chi_frac", "if(abs(ndist)<=wFrac/2,if(scoord>=0,if(scoord<=Lf,1,0),0),0)");
        model.component("comp1").variable("var1").set("perm_xx", "matrix_perm*(1-chi_frac) + (fracture_kn + (fracture_kt-fracture_kn)*txf^2)*chi_frac");
        model.component("comp1").variable("var1").set("perm_xy", "((fracture_kt-fracture_kn)*txf*tyf)*chi_frac");
        model.component("comp1").variable("var1").set("perm_yy", "matrix_perm*(1-chi_frac) + (fracture_kn + (fracture_kt-fracture_kn)*tyf^2)*chi_frac");
        model.component("comp1").variable("var1").set("storage_p", "co2_rho_const*((matrix_phi*matrix_ct)*(1-chi_frac) + (fracture_phi*fracture_ct)*chi_frac)");
        model.component("comp1").variable("var1").set("mob_xx", "(co2_rho_const/co2_mu_const)*perm_xx");
        model.component("comp1").variable("var1").set("mob_xy", "(co2_rho_const/co2_mu_const)*perm_xy");
        model.component("comp1").variable("var1").set("mob_yy", "(co2_rho_const/co2_mu_const)*perm_yy");
        model.component("comp1").variable("var1").set("rhoCp_matrix", "matrix_phi*co2_rho_const*co2_cp_const + (1-matrix_phi)*matrix_rho_r*matrix_cp_r");
        model.component("comp1").variable("var1").set("rhoCp_fracture", "fracture_phi*co2_rho_const*co2_cp_const + (1-fracture_phi)*fracture_rho_r*fracture_cp_r");
        model.component("comp1").variable("var1").set("heatcap_eff", "rhoCp_matrix*(1-chi_frac) + rhoCp_fracture*chi_frac");
        model.component("comp1").variable("var1").set("ktherm_matrix", "matrix_phi*co2_k_const + (1-matrix_phi)*matrix_lambda_r");
        model.component("comp1").variable("var1").set("ktherm_fracture", "fracture_phi*co2_k_const + (1-fracture_phi)*fracture_lambda_r");
        model.component("comp1").variable("var1").set("ktherm_eff", "ktherm_matrix*(1-chi_frac) + ktherm_fracture*chi_frac");
        model.component("comp1").variable("var1").set("uDarcyX", "-(perm_xx*px + perm_xy*py)/co2_mu_const");
        model.component("comp1").variable("var1").set("uDarcyY", "-(perm_xy*px + perm_yy*py)/co2_mu_const");
        model.component("comp1").variable("var1").set("betaTx", "co2_rho_const*co2_cp_const*uDarcyX");
        model.component("comp1").variable("var1").set("betaTy", "co2_rho_const*co2_cp_const*uDarcyY");
    }

    private static void createPressurePde(Model model) {
        model.component("comp1").physics().create("pdeP", "CoefficientFormPDE", "geom1");
        model.component("comp1").physics("pdeP").field("dimensionless").field("p");
        model.component("comp1").physics("pdeP").field("dimensionless").component(new String[]{"p"});
        model.component("comp1").physics("pdeP").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("pdeP").feature("cfeq1").set("da", "storage_p");
        model.component("comp1").physics("pdeP").feature("cfeq1").set("c", new String[]{"mob_xx", "mob_xy", "mob_xy", "mob_yy"});
        model.component("comp1").physics("pdeP").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("pdeP").feature("cfeq1").set("f", "0");
        model.component("comp1").physics("pdeP").feature("init1").set("p", "pInit");
        model.component("comp1").physics("pdeP").create("dirPLeft", "DirichletBoundary", 1);
        model.component("comp1").physics("pdeP").feature("dirPLeft").selection().named("leftBnd");
        model.component("comp1").physics("pdeP").feature("dirPLeft").set("r", "pLeft");
        model.component("comp1").physics("pdeP").create("dirPRight", "DirichletBoundary", 1);
        model.component("comp1").physics("pdeP").feature("dirPRight").selection().named("rightBnd");
        model.component("comp1").physics("pdeP").feature("dirPRight").set("r", "pRight");
    }

    private static void createTemperaturePde(Model model) {
        model.component("comp1").physics().create("pdeT", "CoefficientFormPDE", "geom1");
        model.component("comp1").physics("pdeT").field("dimensionless").field("T");
        model.component("comp1").physics("pdeT").field("dimensionless").component(new String[]{"T"});
        model.component("comp1").physics("pdeT").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("pdeT").feature("cfeq1").set("da", "heatcap_eff");
        model.component("comp1").physics("pdeT").feature("cfeq1").set("c", new String[]{"ktherm_eff", "0", "0", "ktherm_eff"});
        model.component("comp1").physics("pdeT").feature("cfeq1").set("be", new String[]{"betaTx", "betaTy"});
        model.component("comp1").physics("pdeT").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("pdeT").feature("cfeq1").set("f", "0");
        model.component("comp1").physics("pdeT").feature("init1").set("T", "TInit");
        model.component("comp1").physics("pdeT").create("dirTLeft", "DirichletBoundary", 1);
        model.component("comp1").physics("pdeT").feature("dirTLeft").selection().named("leftBnd");
        model.component("comp1").physics("pdeT").feature("dirTLeft").set("r", "TLeft");
        model.component("comp1").physics("pdeT").create("dirTRight", "DirichletBoundary", 1);
        model.component("comp1").physics("pdeT").feature("dirTRight").selection().named("rightBnd");
        model.component("comp1").physics("pdeT").feature("dirTRight").set("r", "TRight");
    }

    private static void configureMesh(Model model, double hmaxM, double hminM) {
        model.component("comp1").mesh("mesh1").automatic(false);
        model.component("comp1").mesh("mesh1").feature().create("size1", "Size");
        model.component("comp1").mesh("mesh1").feature("size1").set("custom", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hauto", 1);
        model.component("comp1").mesh("mesh1").feature("size1").set("hmaxactive", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hminactive", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hmax", formatQuantity(hmaxM, "m"));
        model.component("comp1").mesh("mesh1").feature("size1").set("hmin", formatQuantity(hminM, "m"));
        model.component("comp1").mesh("mesh1").feature("size1").set("hgradactive", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hgrad", 1.2);
        ensureFreeTriFeature(model);
        model.component("comp1").mesh("mesh1").run();
    }

    private static void ensureFreeTriFeature(Model model) {
        try {
            model.component("comp1").mesh("mesh1").feature().create("ftri1", "FreeTri");
        } catch (Exception ex) {
            try {
                model.component("comp1").mesh("mesh1").feature("ftri1");
            } catch (Exception missingFeature) {
                throw new RuntimeException("Failed to create mesh feature ftri1", ex);
            }
        }
    }

    private static void configureStudy(Model model, double[] times) {
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").activate("pdeP", true);
        model.study("std1").feature("time").activate("pdeT", true);
        model.study("std1").feature("time").set("tlist", formatTimeList(times));
    }

    private static void writeProfileReferenceCsvs(
        Model model,
        List<ProfileStation> profileStations,
        List<TimeRequest> profileTimes,
        CasePaths paths
    ) throws IOException {
        for (String family : orderedProfileFamilies()) {
            final List<ProfileStation> familyStations = filterProfileStations(profileStations, family);
            final double[][] coords = buildCoordinateMatrix(familyStations);
            final double[] times = new double[profileTimes.size()];
            for (int i = 0; i < profileTimes.size(); ++i) {
                times[i] = profileTimes.get(i).timeS;
            }

            final double[][] pValues = evaluateAtPoints(model, "interpP_" + family, "p", "Pa", coords, times);
            final double[][] tValues = evaluateAtPoints(model, "interpT_" + family, "T", "K", coords, times);
            for (int timeIndex = 0; timeIndex < profileTimes.size(); ++timeIndex) {
                final TimeRequest request = profileTimes.get(timeIndex);
                final Path csvPath = paths.comsolOutputDir.resolve("comsol_profile_" + family + "_" + request.tag + ".csv");
                try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
                    out.write("station_id,label,family,location,target_axis_m,target_x_m,target_y_m,target_time_s,p_ref_pa,t_ref_k");
                    out.newLine();
                    for (int pointIndex = 0; pointIndex < familyStations.size(); ++pointIndex) {
                        final ProfileStation point = familyStations.get(pointIndex);
                        out.write(String.format(
                            Locale.US,
                            "%d,%s,%s,%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                            point.id, point.label, point.family, point.location,
                            point.targetAxis, point.targetX, point.targetY, request.timeS,
                            pValues[timeIndex][pointIndex], tValues[timeIndex][pointIndex]
                        ));
                        out.newLine();
                    }
                }
            }
        }
    }

    private static void writeMonitorReferenceCsv(
        Model model,
        List<MonitorPoint> monitorPoints,
        List<MonitorSample> monitorSamples,
        CasePaths paths
    ) throws IOException {
        final double[][] coords = buildCoordinateMatrix(monitorPoints);
        final double[] times = new double[monitorSamples.size()];
        for (int i = 0; i < monitorSamples.size(); ++i) {
            times[i] = monitorSamples.get(i).targetTimeS;
        }
        final double[][] pValues = evaluateAtPoints(model, "interpMonitorP", "p", "Pa", coords, times);
        final double[][] tValues = evaluateAtPoints(model, "interpMonitorT", "T", "K", coords, times);
        final Path csvPath = paths.comsolOutputDir.resolve("comsol_monitor_timeseries.csv");
        try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
            out.write("sample_id,target_time_s,actual_time_s");
            for (MonitorPoint point : monitorPoints) {
                out.write(",p_ref_" + point.label);
            }
            for (MonitorPoint point : monitorPoints) {
                out.write(",t_ref_" + point.label);
            }
            out.newLine();
            for (int timeIndex = 0; timeIndex < times.length; ++timeIndex) {
                final MonitorSample sample = monitorSamples.get(timeIndex);
                out.write(String.format(Locale.US, "%d,%.12f,%.12f", sample.sampleId, sample.targetTimeS, sample.targetTimeS));
                for (int pointIndex = 0; pointIndex < monitorPoints.size(); ++pointIndex) {
                    out.write(String.format(Locale.US, ",%.12f", pValues[timeIndex][pointIndex]));
                }
                for (int pointIndex = 0; pointIndex < monitorPoints.size(); ++pointIndex) {
                    out.write(String.format(Locale.US, ",%.12f", tValues[timeIndex][pointIndex]));
                }
                out.newLine();
            }
        }
    }

    private static void writeMeshCheck(
        Model mainModel,
        Model fineModel,
        List<ProfileStation> profileStations,
        List<TimeRequest> profileTimes,
        CaseConfig cfg,
        CasePaths paths,
        double mainHmax,
        double mainHmin,
        double fineHmax,
        double fineHmin
    ) throws IOException {
        final TimeRequest finalRequest = profileTimes.get(profileTimes.size() - 1);
        final double[][] coords = buildCoordinateMatrix(profileStations);
        final double[][] mainP = evaluateAtPoints(mainModel, "interpProfileMainCheckP", "p", "Pa", coords, new double[]{finalRequest.timeS});
        final double[][] fineP = evaluateAtPoints(fineModel, "interpProfileFineCheckP", "p", "Pa", coords, new double[]{finalRequest.timeS});
        final double[][] mainT = evaluateAtPoints(mainModel, "interpProfileMainCheckT", "T", "K", coords, new double[]{finalRequest.timeS});
        final double[][] fineT = evaluateAtPoints(fineModel, "interpProfileFineCheckT", "T", "K", coords, new double[]{finalRequest.timeS});

        double pSumSq = 0.0;
        double pMaxAbs = 0.0;
        double tSumSq = 0.0;
        double tMaxAbs = 0.0;
        for (int i = 0; i < profileStations.size(); ++i) {
            final double pAbsErr = Math.abs(mainP[0][i] - fineP[0][i]);
            final double tAbsErr = Math.abs(mainT[0][i] - fineT[0][i]);
            pSumSq += pAbsErr * pAbsErr;
            tSumSq += tAbsErr * tAbsErr;
            pMaxAbs = Math.max(pMaxAbs, pAbsErr);
            tMaxAbs = Math.max(tMaxAbs, tAbsErr);
        }

        final double pL2Norm = Math.sqrt(pSumSq / Math.max(profileStations.size(), 1)) / cfg.deltaPPa;
        final double pLinfNorm = pMaxAbs / cfg.deltaPPa;
        final double tL2Norm = Math.sqrt(tSumSq / Math.max(profileStations.size(), 1)) / cfg.deltaTK;
        final double tLinfNorm = tMaxAbs / cfg.deltaTK;
        final boolean ok = pL2Norm <= 5.0e-3 && tL2Norm <= 5.0e-3;

        final Path reportPath = paths.comsolOutputDir.resolve("comsol_reference_mesh_check.txt");
        try (BufferedWriter out = Files.newBufferedWriter(reportPath, StandardCharsets.UTF_8)) {
            out.write("profile_tag=" + finalRequest.tag);
            out.newLine();
            out.write("time_s=" + formatDouble(finalRequest.timeS));
            out.newLine();
            out.write("main_mesh_hmax_m=" + formatDouble(mainHmax));
            out.newLine();
            out.write("main_mesh_hmin_m=" + formatDouble(mainHmin));
            out.newLine();
            out.write("fine_mesh_hmax_m=" + formatDouble(fineHmax));
            out.newLine();
            out.write("fine_mesh_hmin_m=" + formatDouble(fineHmin));
            out.newLine();
            out.write("pressure_l2_norm=" + formatDouble(pL2Norm));
            out.newLine();
            out.write("pressure_linf_norm=" + formatDouble(pLinfNorm));
            out.newLine();
            out.write("temperature_l2_norm=" + formatDouble(tL2Norm));
            out.newLine();
            out.write("temperature_linf_norm=" + formatDouble(tLinfNorm));
            out.newLine();
            out.write("reference_acceptable=" + (ok ? "true" : "false"));
            out.newLine();
        }
    }

    private static void writeRunSummary(
        CasePaths paths,
        CaseConfig cfg,
        List<ProfileStation> profileStations,
        List<MonitorPoint> monitorPoints,
        List<TimeRequest> profileTimes,
        List<MonitorSample> monitorSamples,
        double hmaxM,
        double hminM,
        boolean skipFineCheck
    ) throws IOException {
        final Path summaryPath = paths.comsolOutputDir.resolve("comsol_run_summary.md");
        try (BufferedWriter out = Files.newBufferedWriter(summaryPath, StandardCharsets.UTF_8)) {
            out.write("# COMSOL Run Summary\n\n");
            out.write("- Case directory: `" + paths.caseDir + "`\n");
            out.write("- Geometry: `" + formatDouble(cfg.lxM) + " m x " + formatDouble(cfg.lyM) + " m`\n");
            out.write("- Fracture representation: thin-band fallback with aperture `" + formatDouble(cfg.fractureWidthM) + " m`\n");
            out.write("- Gravity: disabled (`0, 0, 0`)\n");
            out.write("- Profile station count: `" + profileStations.size() + "`\n");
            out.write("- Monitor count: `" + monitorPoints.size() + "`\n");
            out.write("- Profile times: `" + profileTimes.size() + "`\n");
            out.write("- Monitor samples: `" + monitorSamples.size() + "`\n");
            out.write("- Main mesh target: hmax=`" + formatDouble(hmaxM) + " m`, hmin=`" + formatDouble(hminM) + " m`\n");
            out.write("- Fine check skipped: `" + (skipFineCheck ? "true" : "false") + "`\n");
        }
    }

    private static void ensureExpectedOutputs(CasePaths paths, List<TimeRequest> profileTimes) {
        final List<Path> required = new ArrayList<Path>();
        required.add(paths.comsolOutputDir.resolve("comsol_monitor_timeseries.csv"));
        required.add(paths.comsolOutputDir.resolve("comsol_model.mph"));
        required.add(paths.comsolOutputDir.resolve("comsol_progress.log"));
        for (String family : orderedProfileFamilies()) {
            for (TimeRequest request : profileTimes) {
                required.add(paths.comsolOutputDir.resolve("comsol_profile_" + family + "_" + request.tag + ".csv"));
            }
        }
        for (Path file : required) {
            if (!Files.isRegularFile(file)) {
                throw new IllegalStateException("Missing expected COMSOL output: " + file);
            }
        }
    }

    private static void writeSkippedFineCheck(CasePaths paths) throws IOException {
        final Path reportPath = paths.comsolOutputDir.resolve("comsol_reference_mesh_check.txt");
        try (BufferedWriter out = Files.newBufferedWriter(reportPath, StandardCharsets.UTF_8)) {
            out.write("fine_check_skipped=true");
            out.newLine();
            out.write("reason=--skip-fine-check was set.");
            out.newLine();
        }
    }

    private static double[][] evaluateAtPoints(
        Model model,
        String numericalTag,
        String expression,
        String unit,
        double[][] coords,
        double[] times
    ) {
        model.result().numerical().create(numericalTag, "Interp");
        final NumericalFeature interp = model.result().numerical(numericalTag);
        interp.set("expr", new String[]{expression});
        interp.set("edim", 2);
        interp.set("coorderr", "on");
        interp.set("unit", new String[]{unit});
        interp.set("data", "dset1");
        interp.set("t", times);
        try {
            interp.selection().all();
        } catch (RuntimeException ignored) {
        }

        final List<InterpolationCandidate> candidates = buildInterpolationCandidates(coords);
        final List<String> diagnostics = new ArrayList<String>();
        for (InterpolationCandidate candidate : candidates) {
            interp.set("coord", candidate.coords);
            interp.run();

            final double[][] raw = tryExtractInterpolationData(interp);
            if (hasData(raw)) {
                return orientAsTimePointMatrix(
                    raw, candidate.pointCount, times.length, numericalTag + "[" + candidate.label + "]"
                );
            }

            diagnostics.add(
                candidate.label +
                "{requestedCoord=" + describeMatrix(candidate.coords) +
                ", storedCoord=" + describeMatrix(getCoordinatesSafely(interp)) +
                ", getReal(true)=" + describeMatrix(getRealSafely(interp, true)) +
                ", getReal(false)=" + describeMatrix(getRealSafely(interp, false)) +
                ", getReal()=" + describeMatrix(getRealDefaultSafely(interp)) +
                ", getData(0)=" + describeMatrix(getDataSafely(interp)) +
                "}"
            );
        }

        throw new IllegalStateException(
            "No interpolation data returned for '" + numericalTag + "'." +
            " expectedPointCount=" + coords[0].length +
            ", expectedTimeCount=" + times.length +
            ", attempts=" + String.join(" | ", diagnostics)
        );
    }

    private static double[][] tryExtractInterpolationData(NumericalFeature interp) {
        final double[][] realColumnwise = getRealSafely(interp, true);
        if (hasData(realColumnwise)) {
            return realColumnwise;
        }

        final double[][] realRowwise = getRealSafely(interp, false);
        if (hasData(realRowwise)) {
            return realRowwise;
        }

        final double[][] defaultReal = getRealDefaultSafely(interp);
        if (hasData(defaultReal)) {
            return defaultReal;
        }

        final double[][] exprData = getDataSafely(interp);
        if (hasData(exprData)) {
            return exprData;
        }

        return new double[0][0];
    }

    private static double[][] getRealSafely(NumericalFeature interp, boolean columnwise) {
        try {
            return interp.getReal(columnwise);
        } catch (RuntimeException ex) {
            return new double[0][0];
        }
    }

    private static double[][] getRealDefaultSafely(NumericalFeature interp) {
        try {
            return interp.getReal();
        } catch (RuntimeException ex) {
            return new double[0][0];
        }
    }

    private static double[][] getDataSafely(NumericalFeature interp) {
        try {
            return interp.getData(0);
        } catch (RuntimeException ex) {
            return new double[0][0];
        }
    }

    private static double[][] getCoordinatesSafely(NumericalFeature interp) {
        try {
            return interp.getCoordinates();
        } catch (RuntimeException ex) {
            return new double[0][0];
        }
    }

    private static boolean hasData(double[][] data) {
        return data != null && data.length > 0 && data[0] != null && data[0].length > 0;
    }

    private static String describeMatrix(double[][] data) {
        if (data == null) {
            return "null";
        }
        if (data.length == 0) {
            return "0x0";
        }
        if (data[0] == null) {
            return data.length + "xnull";
        }
        return data.length + "x" + data[0].length;
    }

    private static double[][] orientAsTimePointMatrix(
        double[][] raw,
        int pointCount,
        int timeCount,
        String numericalTag
    ) {
        if (!hasData(raw)) {
            throw new IllegalStateException(
                "No interpolation data returned for '" + numericalTag + "'."
            );
        }

        final int rows = raw.length;
        final int cols = raw[0].length;

        if (rows == timeCount && cols == pointCount) {
            return copyMatrix(raw);
        }

        if (rows == pointCount && cols == timeCount) {
            return transposeMatrix(raw);
        }

        if (timeCount == 1 && cols == pointCount) {
            final double[][] single = new double[1][pointCount];
            System.arraycopy(raw[0], 0, single[0], 0, pointCount);
            return single;
        }

        if (pointCount == 1 && rows == timeCount) {
            final double[][] single = new double[timeCount][1];
            for (int i = 0; i < timeCount; ++i) {
                single[i][0] = raw[i][0];
            }
            return single;
        }

        throw new IllegalStateException(
            "Unexpected interpolation data shape for '" + numericalTag + "': rows=" + rows +
            ", cols=" + cols + ", expected timeCount=" + timeCount +
            ", pointCount=" + pointCount + "."
        );
    }

    private static double[][] copyMatrix(double[][] src) {
        final double[][] copy = new double[src.length][];
        for (int i = 0; i < src.length; ++i) {
            copy[i] = Arrays.copyOf(src[i], src[i].length);
        }
        return copy;
    }

    private static double[][] transposeMatrix(double[][] src) {
        final int rows = src.length;
        final int cols = src[0].length;
        final double[][] dst = new double[cols][rows];
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                dst[j][i] = src[i][j];
            }
        }
        return dst;
    }

    private static double[][] buildCoordinateMatrix(List<? extends SpatialPoint> points) {
        final double[][] coords = new double[2][points.size()];
        for (int i = 0; i < points.size(); ++i) {
            coords[0][i] = points.get(i).targetX;
            coords[1][i] = points.get(i).targetY;
        }
        return coords;
    }

    private static List<InterpolationCandidate> buildInterpolationCandidates(double[][] coords) {
        final List<InterpolationCandidate> candidates = new ArrayList<InterpolationCandidate>();

        final InterpolationCandidate columnMajor = new InterpolationCandidate();
        columnMajor.label = "columns_as_points";
        columnMajor.coords = copyMatrix(coords);
        columnMajor.pointCount = coords[0].length;
        candidates.add(columnMajor);

        final InterpolationCandidate rowMajor = new InterpolationCandidate();
        rowMajor.label = "rows_as_points";
        rowMajor.coords = transposeMatrix(coords);
        rowMajor.pointCount = coords[0].length;
        candidates.add(rowMajor);

        return candidates;
    }

    private static List<ProfileStation> readProfileStations(Path csvPath) {
        final CsvTable table = readCsv(csvPath);
        final List<ProfileStation> points = new ArrayList<ProfileStation>();
        for (String[] row : table.rows) {
            final ProfileStation point = new ProfileStation();
            point.id = (int) parseRequiredDouble(table, row, "station_id");
            point.label = parseRequiredString(table, row, "label");
            point.family = parseRequiredString(table, row, "family");
            point.location = parseRequiredString(table, row, "location");
            point.targetAxis = parseRequiredDouble(table, row, "target_axis_m");
            point.targetX = parseRequiredDouble(table, row, "target_x_m");
            point.targetY = parseRequiredDouble(table, row, "target_y_m");
            points.add(point);
        }
        return points;
    }

    private static List<MonitorPoint> readMonitorPoints(Path csvPath) {
        final CsvTable table = readCsv(csvPath);
        final List<MonitorPoint> points = new ArrayList<MonitorPoint>();
        for (String[] row : table.rows) {
            final MonitorPoint point = new MonitorPoint();
            point.id = (int) parseRequiredDouble(table, row, "point_id");
            point.label = parseRequiredString(table, row, "label");
            point.location = parseRequiredString(table, row, "location");
            point.targetAxis = parseRequiredDouble(table, row, "target_axis_m");
            point.targetX = parseRequiredDouble(table, row, "target_x_m");
            point.targetY = parseRequiredDouble(table, row, "target_y_m");
            points.add(point);
        }
        return points;
    }

    private static List<TimeRequest> readProfileSchedule(Path csvPath) {
        final CsvTable table = readCsv(csvPath);
        final List<TimeRequest> requests = new ArrayList<TimeRequest>();
        for (String[] row : table.rows) {
            final TimeRequest request = new TimeRequest();
            request.tag = parseRequiredString(table, row, "tag");
            request.timeS = parseRequiredDouble(table, row, "target_time_s");
            requests.add(request);
        }
        return requests;
    }

    private static List<MonitorSample> readMonitorSchedule(Path csvPath) {
        final CsvTable table = readCsv(csvPath);
        final List<MonitorSample> samples = new ArrayList<MonitorSample>();
        for (String[] row : table.rows) {
            final MonitorSample sample = new MonitorSample();
            sample.sampleId = (int) parseRequiredDouble(table, row, "sample_id");
            sample.targetTimeS = parseRequiredDouble(table, row, "target_time_s");
            samples.add(sample);
        }
        return samples;
    }

    private static double[] buildStudyTimes(List<TimeRequest> profileTimes, List<MonitorSample> monitorSamples) {
        final Set<Double> values = new LinkedHashSet<Double>();
        for (MonitorSample sample : monitorSamples) {
            values.add(sample.targetTimeS);
        }
        for (TimeRequest request : profileTimes) {
            values.add(request.timeS);
        }
        final List<Double> ordered = new ArrayList<Double>(values);
        Collections.sort(ordered);
        final double[] out = new double[ordered.size()];
        for (int i = 0; i < ordered.size(); ++i) {
            out[i] = ordered.get(i).doubleValue();
        }
        return out;
    }

    private static String formatTimeList(double[] times) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < times.length; ++i) {
            if (i > 0) {
                sb.append(',');
            }
            sb.append(formatDouble(times[i]));
        }
        return sb.toString();
    }

    private static List<String> orderedProfileFamilies() {
        final List<String> families = new ArrayList<String>();
        families.add("matrix_horizontal");
        families.add("fracture_tangent");
        families.add("cross_normal");
        return families;
    }

    private static List<ProfileStation> filterProfileStations(List<ProfileStation> stations, String family) {
        final List<ProfileStation> filtered = new ArrayList<ProfileStation>();
        for (ProfileStation station : stations) {
            if (family.equals(station.family)) {
                filtered.add(station);
            }
        }
        return filtered;
    }

    private static CsvTable readCsv(Path csvPath) {
        try (BufferedReader reader = Files.newBufferedReader(csvPath, StandardCharsets.UTF_8)) {
            final String headerLine = reader.readLine();
            if (headerLine == null) {
                throw new IllegalStateException("CSV file is empty: " + csvPath);
            }

            final CsvTable table = new CsvTable();
            table.headers = splitCsvLine(headerLine);
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.trim().isEmpty()) {
                    continue;
                }
                final String[] row = splitCsvLine(line);
                if (row.length != table.headers.length) {
                    throw new IllegalStateException("CSV column count mismatch in " + csvPath + ": " + line);
                }
                table.rows.add(row);
            }
            return table;
        } catch (IOException ex) {
            throw new IllegalStateException("Failed to read CSV file: " + csvPath, ex);
        }
    }

    private static String[] splitCsvLine(String line) {
        final String[] raw = line.split(",", -1);
        for (int i = 0; i < raw.length; ++i) {
            raw[i] = raw[i].trim();
        }
        return raw;
    }

    private static double parseRequiredDouble(CsvTable table, String[] row, String columnName) {
        final String text = parseRequiredString(table, row, columnName);
        return Double.parseDouble(text);
    }

    private static String parseRequiredString(CsvTable table, String[] row, String columnName) {
        final int index = table.indexOf(columnName);
        if (index < 0 || index >= row.length) {
            throw new IllegalStateException("Missing CSV column: " + columnName);
        }
        return row[index];
    }

    private static String firstNonEmpty(String value) {
        if (value == null) {
            return null;
        }
        final String trimmed = value.trim();
        return trimmed.isEmpty() ? null : trimmed;
    }

    private static String formatQuantity(double value, String unit) {
        return formatDouble(value) + "[" + unit + "]";
    }

    private static String formatDouble(double value) {
        return String.format(Locale.US, "%.12f", value);
    }

}

final class CasePaths {
    Path caseDir;
    Path engineeringDir;
    Path referenceDir;
    Path comsolInputDir;
    Path comsolOutputDir;
    Path profileStationsCsv;
    Path monitorPointsCsv;
    Path profileScheduleCsv;
    Path monitorScheduleCsv;
    Path propertyTableCsv;
}

final class TimeRequest {
    String tag;
    double timeS;
}

class SpatialPoint {
    double targetAxis;
    double targetX;
    double targetY;
}

final class ProfileStation extends SpatialPoint {
    int id;
    String label;
    String family;
    String location;
}

final class MonitorPoint extends SpatialPoint {
    int id;
    String label;
    String location;
}

final class MonitorSample {
    int sampleId;
    double targetTimeS;
}

final class CaseConfig {
    double lxM;
    double lyM;
    double fracX0Ratio;
    double fracY0Ratio;
    double fracX1Ratio;
    double fracY1Ratio;
    double pInitPa;
    double pLeftPa;
    double pRightPa;
    double tInitK;
    double tLeftK;
    double tRightK;
    double matrixPhi;
    double matrixPerm;
    double matrixCt;
    double matrixRhoR;
    double matrixCpR;
    double matrixLambdaR;
    double fracturePhi;
    double fractureKt;
    double fractureKn;
    double fractureCt;
    double fractureRhoR;
    double fractureCpR;
    double fractureLambdaR;
    double co2Rho;
    double co2Mu;
    double co2Cp;
    double co2K;
    double fractureWidthM;
    double deltaPPa;
    double deltaTK;
}

final class InterpolationCandidate {
    String label;
    double[][] coords;
    int pointCount;
}

final class CsvTable {
    String[] headers;
    final List<String[]> rows = new ArrayList<String[]>();

    int indexOf(String columnName) {
        for (int i = 0; i < headers.length; ++i) {
            if (headers[i].equals(columnName)) {
                return i;
            }
        }
        return -1;
    }
}
