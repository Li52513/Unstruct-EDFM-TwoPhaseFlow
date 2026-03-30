import com.comsol.model.Model;
import com.comsol.model.NumericalFeature;
import com.comsol.model.SelectionFeature;
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
    private static final double MAIN_HMAX_M = 2.0;
    private static final double FINE_HMAX_M = 1.0;
    private static final double MAIN_HMIN_M = 0.2;
    private static final double FINE_HMIN_M = 0.1;
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
            final double[] densePressureOnlyTimes = buildUniformStudyTimes(
                studyTimes[studyTimes.length - 1],
                Math.max(1.0, cfg.dtInitS)
            );

            ModelUtil.initStandalone(false);
            ModelUtil.loadPreferences();
            configureStandaloneRuntime(paths);
            ModelUtil.showProgress(paths.comsolOutputDir.resolve("comsol_progress.log").toString());

            final double mainHmin = MAIN_HMIN_M;
            final Model matrixOnlyPressureOnlyModel = buildAndSolveMatrixOnlyPressureOnlyDiagnosticModel(
                paths, cfg, densePressureOnlyTimes, "ModelMatrixOnlyPressureOnly", MAIN_HMAX_M, mainHmin
            );
            writeMatrixOnlyPressureOnlyDiagnostic(
                matrixOnlyPressureOnlyModel, paths, profileStations, profileTimes, cfg, densePressureOnlyTimes
            );
            writeMatrixOnlyPressureOnlyBoundaryDiagnostics(
                matrixOnlyPressureOnlyModel,
                paths,
                cfg,
                profileTimes.get(profileTimes.size() - 1).timeS
            );
            writeMatrixOnlyPressureOnlyHorizontalDiagnostics(
                matrixOnlyPressureOnlyModel,
                paths,
                cfg,
                profileTimes.get(profileTimes.size() - 1).timeS
            );
            matrixOnlyPressureOnlyModel.save(
                paths.comsolOutputDir.resolve("comsol_matrix_only_pressureonly_diagnostic.mph").toString()
            );

            final Model matrixOnlyPressureOnlyTransplantModel = buildAndSolveMatrixOnlyPressureOnlyTransplantModel(
                paths, cfg, densePressureOnlyTimes, "ModelMatrixOnlyPressureOnlyTransplant", MAIN_HMAX_M, mainHmin
            );
            writeMatrixOnlyPressureOnlyDiagnostic(
                matrixOnlyPressureOnlyTransplantModel,
                paths,
                profileStations,
                profileTimes,
                cfg,
                densePressureOnlyTimes,
                "comsol_matrix_only_pressureonly_transplant",
                "matrix_only_pressureonly_transplant_diag",
                "diagMatrixOnlyPressureOnlyTransplant"
            );
            writeMatrixOnlyPressureOnlyBoundaryDiagnostics(
                matrixOnlyPressureOnlyTransplantModel,
                paths,
                cfg,
                profileTimes.get(profileTimes.size() - 1).timeS,
                "comsol_matrix_only_pressureonly_transplant_boundary_diagnostics.csv",
                "matrix_only_pressureonly_transplant_boundary_line",
                "diagBoundaryTransplantLine"
            );
            writeMatrixOnlyPressureOnlyHorizontalDiagnostics(
                matrixOnlyPressureOnlyTransplantModel,
                paths,
                cfg,
                profileTimes.get(profileTimes.size() - 1).timeS,
                "comsol_matrix_only_pressureonly_transplant_horizontal_diagnostics.csv",
                "matrix_only_pressureonly_transplant_horizontal_line",
                "diagHorizontalTransplantLine"
            );
            matrixOnlyPressureOnlyTransplantModel.save(
                paths.comsolOutputDir.resolve("comsol_matrix_only_pressureonly_transplant_diagnostic.mph").toString()
            );

            try {
                final Model matrixOnlyPressureOnlyStationaryModel = buildAndSolveMatrixOnlyPressureOnlyStationaryDiagnosticModel(
                    paths, cfg, "ModelMatrixOnlyPressureOnlyStationary", MAIN_HMAX_M, mainHmin
                );
                writeMatrixOnlyPressureOnlyStationaryDiagnostic(
                    matrixOnlyPressureOnlyStationaryModel,
                    paths,
                    profileStations,
                    cfg,
                    "comsol_matrix_only_pressureonly_stationary",
                    "matrix_only_pressureonly_stationary_diag",
                    "diagMatrixOnlyPressureOnlyStationary"
                );
                matrixOnlyPressureOnlyStationaryModel.save(
                    paths.comsolOutputDir.resolve("comsol_matrix_only_pressureonly_stationary_diagnostic.mph").toString()
                );
            } catch (Exception stationaryEx) {
                System.out.println(
                    "[COMSOL] matrix_only_pressureonly_stationary_diag skipped: " +
                    stationaryEx.getMessage()
                );
            }

            final Model matrixOnlyModel = buildAndSolveMatrixOnlyDiagnosticModel(
                paths, cfg, studyTimes, "ModelMatrixOnly", MAIN_HMAX_M, mainHmin
            );
            writeMatrixOnlyDiagnostic(matrixOnlyModel, paths, profileStations, profileTimes);
            matrixOnlyModel.save(paths.comsolOutputDir.resolve("comsol_matrix_only_diagnostic.mph").toString());

            final Model mainModel = buildAndSolveModel(paths, cfg, studyTimes, "ModelMain", "main", MAIN_HMAX_M, mainHmin);
            writeFractureExchangeDiagnostics(mainModel, paths, cfg, profileTimes.get(profileTimes.size() - 1).timeS);
            writeProfileReferenceCsvs(mainModel, profileStations, profileTimes, paths);
            writeMonitorReferenceCsv(mainModel, monitorPoints, monitorSamples, paths);
            writeRunSummary(paths, cfg, profileStations, monitorPoints, profileTimes, monitorSamples, MAIN_HMAX_M, mainHmin, skipFineCheck);
            mainModel.save(paths.comsolOutputDir.resolve("comsol_model.mph").toString());

            if (!skipFineCheck) {
                final double fineHmin = FINE_HMIN_M;
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
        cfg.dtInitS = requireValue(values, "dt_init");
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
        cfg.fractureWidthM = values.containsKey("fracture_aperture_m")
            ? requireValue(values, "fracture_aperture_m")
            : requireValue(values, "comsol_thin_band_width");
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
        createMatrixTraceCouplings(model);
        createVariables(model);
        createFullPressureCoefficientFunctions(
            model,
            prepareFullPressureCoefficientTables(paths, cfg, "main_pressure_coeffs_" + labelSuffix)
        );
        createMatrixPressurePdes(model);
        createFracturePressurePde(model);
        createMatrixTemperaturePdes(model);
        createFractureTemperaturePde(model);
        logShapeDiagnostics(model);
        configureMesh(model, hmaxM, hminM);
        configureStudy(model, studyTimes);
        configureSegregatedTimeSolver(model, studyTimes, true);
        model.sol("sol1").runAll();
        return model;
    }

    private static Model buildAndSolveMatrixOnlyDiagnosticModel(
        CasePaths paths,
        CaseConfig cfg,
        double[] studyTimes,
        String modelTag,
        double hmaxM,
        double hminM
    ) {
        final Model model = ModelUtil.create(modelTag);
        model.modelPath(paths.comsolOutputDir.toString());
        model.label("ComsolHTCO2ConstPPSingleFracNoWell_matrix_only_diagnostic.mph");

        configureParameters(model, cfg);
        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 2);
        model.component("comp1").mesh().create("mesh1");

        buildMatrixOnlyGeometry(model);
        createExternalBoundarySelections(model, cfg);
        createMatrixOnlyVariables(model);
        createFullPressureCoefficientFunctions(
            model,
            prepareFullPressureCoefficientTables(paths, cfg, "matrix_only_pressure_coeffs")
        );
        createMatrixPressurePdes(model, false);
        createMatrixTemperaturePdes(model, false);
        configureMesh(model, hmaxM, hminM);
        configureStudy(model, studyTimes);
        configureSegregatedTimeSolver(model, studyTimes, false);
        model.sol("sol1").runAll();
        return model;
    }

    private static Model buildAndSolveMatrixOnlyPressureOnlyDiagnosticModel(
        CasePaths paths,
        CaseConfig cfg,
        double[] studyTimes,
        String modelTag,
        double hmaxM,
        double hminM
    ) {
        final Model model = ModelUtil.create(modelTag);
        model.modelPath(paths.comsolOutputDir.toString());
        model.label("ComsolHTCO2ConstPPSingleFracNoWell_matrix_only_pressureonly_diagnostic.mph");

        configureParameters(model, cfg);
        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 2);
        model.component("comp1").mesh().create("mesh1");

        buildMatrixOnlyGeometry(model);
        createExternalBoundarySelections(model, cfg);
        createMatrixOnlyLegacyPressureVariables(model);
        final LegacyPressureOnlyTables diagnosticTables = prepareLegacyPressureOnlyTables(
            paths, cfg, "diag_matrix_only_pressureonly"
        );
        createLegacyPressureOnlyFunctions(model, diagnosticTables);
        createLegacyPressureOnlyPde(model);
        configureLegacyStyleMesh(model, hmaxM, hminM);
        configureLegacyPressureOnlyUniformStudy(model, studyTimes[studyTimes.length - 1], studyTimes[1] - studyTimes[0]);
        model.study("std1").run();
        return model;
    }

    private static Model buildAndSolveMatrixOnlyPressureOnlyTransplantModel(
        CasePaths paths,
        CaseConfig cfg,
        double[] studyTimes,
        String modelTag,
        double hmaxM,
        double hminM
    ) {
        final Model model = ModelUtil.create(modelTag);
        model.modelPath(paths.comsolOutputDir.toString());
        model.label("ComsolHTCO2ConstPPSingleFracNoWell_matrix_only_pressureonly_transplant.mph");

        configureParameters(model, cfg);
        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 2);
        model.component("comp1").mesh().create("mesh1");

        buildMatrixOnlyGeometry(model);
        createTransplantBoundarySelections(model, cfg);
        createMatrixOnlyLegacyPressureVariables(model);
        final LegacyPressureOnlyTables transplantTables = prepareLegacyPressureOnlyTables(
            paths, cfg, "diag_matrix_only_pressureonly_transplant"
        );
        createTransplantPressureOnlyFunctions(model, transplantTables);
        createTransplantPressureOnlyPde(model);
        configureLegacyStyleMesh(model, hmaxM, hminM);
        configureLegacyPressureOnlyUniformStudy(model, studyTimes[studyTimes.length - 1], studyTimes[1] - studyTimes[0]);
        model.study("std1").run();
        return model;
    }

    private static Model buildAndSolveMatrixOnlyPressureOnlyStationaryDiagnosticModel(
        CasePaths paths,
        CaseConfig cfg,
        String modelTag,
        double hmaxM,
        double hminM
    ) {
        final Model model = ModelUtil.create(modelTag);
        model.modelPath(paths.comsolOutputDir.toString());
        model.label("ComsolHTCO2ConstPPSingleFracNoWell_matrix_only_pressureonly_stationary.mph");

        configureParameters(model, cfg);
        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 2);
        model.component("comp1").mesh().create("mesh1");

        buildMatrixOnlyGeometry(model);
        createTransplantBoundarySelections(model, cfg);
        createMatrixOnlyLegacyPressureVariables(model);
        createTransplantPressureOnlyStationaryPde(model);
        configureLegacyStyleMesh(model, hmaxM, hminM);
        configureStationaryPressureOnlyStudy(model);
        model.study("std1").run();
        return model;
    }

    private static void createMatrixTraceCouplings(Model model) {
        createMatrixTraceCoupling(model, "extUp", "x+fracture_map_eps*nxf", "y+fracture_map_eps*nyf");
        createMatrixTraceCoupling(model, "extDown", "x-fracture_map_eps*nxf", "y-fracture_map_eps*nyf");
    }

    private static void createMatrixTraceCoupling(
        Model model,
        String tag,
        String xExpr,
        String yExpr
    ) {
        try {
            model.component("comp1").cpl().remove(tag);
        } catch (RuntimeException ignored) {
        }
        model.component("comp1").cpl().create(tag, "GeneralExtrusion");
        model.component("comp1").cpl(tag).selection().named("matrixAllDom");
        model.component("comp1").cpl(tag).set("method", "closest");
        model.component("comp1").cpl(tag).set("usesrcmap", "off");
        model.component("comp1").cpl(tag).set("dstmap", new String[]{xExpr, yExpr});
        model.component("comp1").cpl(tag).label("MatrixTrace_" + tag);
    }

    private static void logShapeDiagnostics(Model model) {
        try {
            final String[] tags = model.shape().tags();
            final StringBuilder sb = new StringBuilder();
            for (String tag : tags) {
                if (sb.length() > 0) sb.append("; ");
                sb.append(tag).append("->").append(Arrays.toString(model.shape(tag).fieldVariable()));
            }
            System.out.println("[COMSOL] shape_tags=" + sb.toString());
        } catch (RuntimeException ex) {
            System.out.println("[COMSOL] shape_tags_unavailable=" + ex.getMessage());
        }
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
        model.param().set("Lf", "sqrt((x1f-x0f)^2+(y1f-y0f)^2)");
        model.param().set("txf", "(x1f-x0f)/Lf");
        model.param().set("tyf", "(y1f-y0f)/Lf");
        model.param().set("nxf", "-tyf");
        model.param().set("nyf", "txf");
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
        model.param().set("matrix_dx_ref", "Lx/48");
        model.param().set("matrix_dy_ref", "Ly/6");
        model.param().set(
            "mf_transfer_half_distance",
            "max(0.5[m],0.5*(abs(nxf)*matrix_dx_ref + abs(nyf)*matrix_dy_ref))"
        );
        model.param().set("matrix_storage", "co2_rho_const*matrix_phi*matrix_ct");
        model.param().set("fracture_storage", "co2_rho_const*wFrac*fracture_phi*fracture_ct");
        model.param().set("matrix_mob", "co2_rho_const*matrix_perm/co2_mu_const");
        model.param().set("fracture_mob_t", "co2_rho_const*wFrac*fracture_kt/co2_mu_const");
        model.param().set("pmf_coeff", "co2_rho_const*fracture_kn/(co2_mu_const*mf_transfer_half_distance)");
        model.param().set("fracture_map_eps", "max(0.25[m],5*wFrac)");
        model.param().set("rhoCp_matrix", "matrix_phi*co2_rho_const*co2_cp_const + (1-matrix_phi)*matrix_rho_r*matrix_cp_r");
        model.param().set("rhoCp_fracture", "fracture_phi*co2_rho_const*co2_cp_const + (1-fracture_phi)*fracture_rho_r*fracture_cp_r");
        model.param().set("ktherm_matrix", "matrix_phi*co2_k_const + (1-matrix_phi)*matrix_lambda_r");
        model.param().set("ktherm_fracture", "fracture_phi*co2_k_const + (1-fracture_phi)*fracture_lambda_r");
        model.param().set("lambda_int", "2*ktherm_matrix*ktherm_fracture/max(ktherm_matrix+ktherm_fracture,1e-12[W/(m*K)])");
        model.param().set("ht_coeff", "lambda_int/mf_transfer_half_distance");
    }

    private static void buildGeometry(Model model) {
        model.component("comp1").geom("geom1").create("r1", "Rectangle");
        model.component("comp1").geom("geom1").feature("r1").set("base", "corner");
        model.component("comp1").geom("geom1").feature("r1").set("pos", new String[]{"0", "0"});
        model.component("comp1").geom("geom1").feature("r1").set("size", new String[]{"Lx", "Ly"});
        model.component("comp1").geom("geom1").create("lsFrac", "LineSegment");
        model.component("comp1").geom("geom1").feature("lsFrac").set("specify1", "coord");
        model.component("comp1").geom("geom1").feature("lsFrac").set("specify2", "coord");
        model.component("comp1").geom("geom1").feature("lsFrac").set("coord1", new String[]{"x0f", "y0f"});
        model.component("comp1").geom("geom1").feature("lsFrac").set("coord2", new String[]{"x1f", "y1f"});
        model.component("comp1").geom("geom1").feature("lsFrac").set("selresult", "on");
        model.component("comp1").geom("geom1").feature("lsFrac").set("selresultshow", "bnd");
        model.component("comp1").geom("geom1").run();
    }

    private static void buildMatrixOnlyGeometry(Model model) {
        model.component("comp1").geom("geom1").create("r1", "Rectangle");
        model.component("comp1").geom("geom1").feature("r1").set("base", "corner");
        model.component("comp1").geom("geom1").feature("r1").set("pos", new String[]{"0", "0"});
        model.component("comp1").geom("geom1").feature("r1").set("size", new String[]{"Lx", "Ly"});
        model.component("comp1").geom("geom1").run();
    }

    private static void createBoundarySelections(Model model, CaseConfig cfg) {
        final double tol = Math.max(1.0e-9, 1.0e-6 * Math.max(cfg.lxM, cfg.lyM));
        final double x0 = cfg.fracX0Ratio * cfg.lxM;
        final double y0 = cfg.fracY0Ratio * cfg.lyM;
        final double x1 = cfg.fracX1Ratio * cfg.lxM;
        final double y1 = cfg.fracY1Ratio * cfg.lyM;
        final double dx = x1 - x0;
        final double dy = y1 - y0;
        final double len = Math.sqrt(dx * dx + dy * dy);
        final double nx = -dy / len;
        final double ny = dx / len;
        final double xm = 0.5 * (x0 + x1);
        final double ym = 0.5 * (y0 + y1);
        final double probeOffset = 0.15 * cfg.lyM;
        final double probeTol = Math.max(0.5, 0.02 * Math.min(cfg.lxM, cfg.lyM));

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

        model.component("comp1").selection().create("bottomBnd", "Box");
        model.component("comp1").selection("bottomBnd").set("entitydim", "1");
        model.component("comp1").selection("bottomBnd").set("condition", "inside");
        model.component("comp1").selection("bottomBnd").set("xmin", formatDouble(-tol));
        model.component("comp1").selection("bottomBnd").set("xmax", formatDouble(cfg.lxM + tol));
        model.component("comp1").selection("bottomBnd").set("ymin", formatDouble(-tol));
        model.component("comp1").selection("bottomBnd").set("ymax", formatDouble(tol));
        model.component("comp1").selection("bottomBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("bottomBnd").set("zmax", formatDouble(tol));

        model.component("comp1").selection().create("topBnd", "Box");
        model.component("comp1").selection("topBnd").set("entitydim", "1");
        model.component("comp1").selection("topBnd").set("condition", "inside");
        model.component("comp1").selection("topBnd").set("xmin", formatDouble(-tol));
        model.component("comp1").selection("topBnd").set("xmax", formatDouble(cfg.lxM + tol));
        model.component("comp1").selection("topBnd").set("ymin", formatDouble(cfg.lyM - tol));
        model.component("comp1").selection("topBnd").set("ymax", formatDouble(cfg.lyM + tol));
        model.component("comp1").selection("topBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("topBnd").set("zmax", formatDouble(tol));

        model.component("comp1").selection().create("fractureBndBox", "Box");
        model.component("comp1").selection("fractureBndBox").set("entitydim", "1");
        model.component("comp1").selection("fractureBndBox").set("condition", "inside");
        model.component("comp1").selection("fractureBndBox").set("xmin", formatDouble(Math.min(x0, x1) - tol));
        model.component("comp1").selection("fractureBndBox").set("xmax", formatDouble(Math.max(x0, x1) + tol));
        model.component("comp1").selection("fractureBndBox").set("ymin", formatDouble(Math.min(y0, y1) - tol));
        model.component("comp1").selection("fractureBndBox").set("ymax", formatDouble(Math.max(y0, y1) + tol));
        model.component("comp1").selection("fractureBndBox").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("fractureBndBox").set("zmax", formatDouble(tol));

        createExplicitFractureBoundarySelection(model);
        createExplicitMatrixDomainSelection(model);

        model.component("comp1").selection().create("matrixUpDom", "Box");
        model.component("comp1").selection("matrixUpDom").set("entitydim", "2");
        model.component("comp1").selection("matrixUpDom").set("condition", "intersects");
        model.component("comp1").selection("matrixUpDom").set("xmin", formatDouble(xm + probeOffset * nx - probeTol));
        model.component("comp1").selection("matrixUpDom").set("xmax", formatDouble(xm + probeOffset * nx + probeTol));
        model.component("comp1").selection("matrixUpDom").set("ymin", formatDouble(ym + probeOffset * ny - probeTol));
        model.component("comp1").selection("matrixUpDom").set("ymax", formatDouble(ym + probeOffset * ny + probeTol));
        model.component("comp1").selection("matrixUpDom").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("matrixUpDom").set("zmax", formatDouble(tol));

        model.component("comp1").selection().create("matrixDownDom", "Box");
        model.component("comp1").selection("matrixDownDom").set("entitydim", "2");
        model.component("comp1").selection("matrixDownDom").set("condition", "intersects");
        model.component("comp1").selection("matrixDownDom").set("xmin", formatDouble(xm - probeOffset * nx - probeTol));
        model.component("comp1").selection("matrixDownDom").set("xmax", formatDouble(xm - probeOffset * nx + probeTol));
        model.component("comp1").selection("matrixDownDom").set("ymin", formatDouble(ym - probeOffset * ny - probeTol));
        model.component("comp1").selection("matrixDownDom").set("ymax", formatDouble(ym - probeOffset * ny + probeTol));
        model.component("comp1").selection("matrixDownDom").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("matrixDownDom").set("zmax", formatDouble(tol));
    }

    private static void createExternalBoundarySelections(Model model, CaseConfig cfg) {
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

        model.component("comp1").selection().create("bottomBnd", "Box");
        model.component("comp1").selection("bottomBnd").set("entitydim", "1");
        model.component("comp1").selection("bottomBnd").set("condition", "inside");
        model.component("comp1").selection("bottomBnd").set("xmin", formatDouble(-tol));
        model.component("comp1").selection("bottomBnd").set("xmax", formatDouble(cfg.lxM + tol));
        model.component("comp1").selection("bottomBnd").set("ymin", formatDouble(-tol));
        model.component("comp1").selection("bottomBnd").set("ymax", formatDouble(tol));
        model.component("comp1").selection("bottomBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("bottomBnd").set("zmax", formatDouble(tol));

        model.component("comp1").selection().create("topBnd", "Box");
        model.component("comp1").selection("topBnd").set("entitydim", "1");
        model.component("comp1").selection("topBnd").set("condition", "inside");
        model.component("comp1").selection("topBnd").set("xmin", formatDouble(-tol));
        model.component("comp1").selection("topBnd").set("xmax", formatDouble(cfg.lxM + tol));
        model.component("comp1").selection("topBnd").set("ymin", formatDouble(cfg.lyM - tol));
        model.component("comp1").selection("topBnd").set("ymax", formatDouble(cfg.lyM + tol));
        model.component("comp1").selection("topBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("topBnd").set("zmax", formatDouble(tol));
    }

    private static void createTransplantBoundarySelections(Model model, CaseConfig cfg) {
        final double tol = 1.0e-6;

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

    private static void createExplicitFractureBoundarySelection(Model model) {
        final int[] lineEntities = getSelectionEntitiesSafely(model, "geom1_lsFrac_bnd", 1);
        if (lineEntities.length == 0) {
            throw new IllegalStateException(
                "Unable to resolve fracture line selection 'geom1_lsFrac_bnd' from geometry."
            );
        }

        try {
            model.component("comp1").selection().remove("fractureBnd");
        } catch (RuntimeException ignored) {
        }
        final SelectionFeature fractureSelection = model.component("comp1").selection().create("fractureBnd", "Explicit");
        fractureSelection.geom("geom1", 1);
        fractureSelection.set(lineEntities);

        final int resolvedCount = fractureSelection.entities().length;
        if (resolvedCount != 1) {
            throw new IllegalStateException(
                "Expected one explicit fracture line entity, but got " + resolvedCount +
                ". line_count=" + lineEntities.length
            );
        }

        System.out.println(
            "[COMSOL] fracture_selection_source=geom1_lsFrac_bnd" +
            ", entity_count=" + resolvedCount
        );
    }

    private static void createExplicitMatrixDomainSelection(Model model) {
        int[] matrixDomainEntities = getSelectionEntitiesSafely(model, "geom1_r1_dom", 2);
        String sourceTag = "geom1_r1_dom";
        if (matrixDomainEntities.length == 0) {
            matrixDomainEntities = getSelectionEntitiesSafely(model, "geom1_dom", 2);
            sourceTag = "geom1_dom";
        }
        if (matrixDomainEntities.length == 0) {
            matrixDomainEntities = new int[]{1};
            sourceTag = "explicit_domain_1";
        }

        try {
            model.component("comp1").selection().remove("matrixAllDom");
        } catch (RuntimeException ignored) {
        }

        final SelectionFeature matrixAllSelection = model.component("comp1").selection().create("matrixAllDom", "Explicit");
        matrixAllSelection.geom("geom1", 2);
        matrixAllSelection.set(matrixDomainEntities);

        System.out.println(
            "[COMSOL] matrix_domain_selection_source=" + sourceTag +
            ", entity_count=" + matrixDomainEntities.length
        );
    }

    private static int[] getSelectionEntitiesSafely(Model model, String selectionTag, int dim) {
        try {
            return model.selection(selectionTag).entities(dim);
        } catch (RuntimeException first) {
            try {
                return model.selection(selectionTag).entities();
            } catch (RuntimeException second) {
                return new int[0];
            }
        }
    }

    private static void createVariables(Model model) {
        model.component("comp1").variable().create("var1");
        model.component("comp1").variable("var1").set("uMx", "-(matrix_perm/co2_mu_const)*pmx");
        model.component("comp1").variable("var1").set("uMy", "-(matrix_perm/co2_mu_const)*pmy");
        model.component("comp1").variable("var1").set("betaTm_x", "co2_rho_const*co2_cp_const*uMx");
        model.component("comp1").variable("var1").set("betaTm_y", "co2_rho_const*co2_cp_const*uMy");
        model.component("comp1").variable("var1").set("dpfds", "pfTx*txf + pfTy*tyf");
        model.component("comp1").variable("var1").set("uFrac", "-(fracture_kt/co2_mu_const)*dpfds");
        model.component("comp1").variable("var1").set("betaTf_x", "wFrac*co2_rho_const*co2_cp_const*uFrac*txf");
        model.component("comp1").variable("var1").set("betaTf_y", "wFrac*co2_rho_const*co2_cp_const*uFrac*tyf");
        model.component("comp1").variable("var1").set("qmf_up", "pmf_coeff*(extUp(pm)-pf)");
        model.component("comp1").variable("var1").set("qmf_down", "pmf_coeff*(extDown(pm)-pf)");
        model.component("comp1").variable("var1").set("hAdv_up", "qmf_up*co2_cp_const*if(qmf_up>=0,extUp(Tm),Tf)");
        model.component("comp1").variable("var1").set("hAdv_down", "qmf_down*co2_cp_const*if(qmf_down>=0,extDown(Tm),Tf)");
        model.component("comp1").variable("var1").set("Hmf_up", "ht_coeff*(extUp(Tm)-Tf) + hAdv_up");
        model.component("comp1").variable("var1").set("Hmf_down", "ht_coeff*(extDown(Tm)-Tf) + hAdv_down");
        model.component("comp1").variable("var1").set("pMatrixPlot", "pm");
        model.component("comp1").variable("var1").set("TMatrixPlot", "Tm");

        model.component("comp1").variable().create("varFrac");
        model.component("comp1").variable("varFrac").selection().named("fractureBnd");
        model.component("comp1").variable("varFrac").set("xMapUp", "extUp(x)");
        model.component("comp1").variable("varFrac").set("yMapUp", "extUp(y)");
        model.component("comp1").variable("varFrac").set("xMapDown", "extDown(x)");
        model.component("comp1").variable("varFrac").set("yMapDown", "extDown(y)");
        model.component("comp1").variable("varFrac").set("pmUpTrace", "extUp(pm)");
        model.component("comp1").variable("varFrac").set("pmDownTrace", "extDown(pm)");
        model.component("comp1").variable("varFrac").set("TmUpTrace", "extUp(Tm)");
        model.component("comp1").variable("varFrac").set("TmDownTrace", "extDown(Tm)");
        model.component("comp1").variable("varFrac").set("qmfTotal", "qmf_up + qmf_down");
        model.component("comp1").variable("varFrac").set("HmfTotal", "Hmf_up + Hmf_down");
    }

    private static void createMatrixPressurePdes(Model model) {
        createMatrixPressurePdes(model, true);
    }

    private static void createMatrixOnlyVariables(Model model) {
        model.component("comp1").variable().create("var1");
        model.component("comp1").variable("var1").set("uMx", "-(matrix_perm/co2_mu_const)*pmx");
        model.component("comp1").variable("var1").set("uMy", "-(matrix_perm/co2_mu_const)*pmy");
        model.component("comp1").variable("var1").set("betaTm_x", "co2_rho_const*co2_cp_const*uMx");
        model.component("comp1").variable("var1").set("betaTm_y", "co2_rho_const*co2_cp_const*uMy");
        model.component("comp1").variable("var1").set("pMatrixPlot", "pm");
        model.component("comp1").variable("var1").set("TMatrixPlot", "Tm");
    }

    private static void createMatrixOnlyPressureVariables(Model model) {
        model.component("comp1").variable().create("var1");
        model.component("comp1").variable("var1").set("pMatrixPlot", "pm");
    }

    private static void createMatrixOnlyLegacyPressureVariables(Model model) {
        model.component("comp1").variable().create("var1");
        model.component("comp1").variable("var1").set("pMatrixPlot", "p");
    }

    private static void createMatrixPressurePdes(Model model, boolean includeFractureCoupling) {
        model.component("comp1").physics().create("pdePm", "CoefficientFormPDE", "geom1");
        model.component("comp1").physics("pdePm").field("dimensionless").field("pm");
        model.component("comp1").physics("pdePm").field("dimensionless").component(new String[]{"pm"});
        model.component("comp1").physics("pdePm").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("pdePm").feature("cfeq1").set("da", "matrix_storage_fun(pm)");
        model.component("comp1").physics("pdePm").feature("cfeq1").set("c", "matrix_mob_fun(pm)");
        model.component("comp1").physics("pdePm").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("pdePm").feature("cfeq1").set("f", "0");
        model.component("comp1").physics("pdePm").feature("init1").set("pm", "pInit");
        model.component("comp1").physics("pdePm").create("dirPLeft", "DirichletBoundary", 1);
        model.component("comp1").physics("pdePm").feature("dirPLeft").selection().named("leftBnd");
        model.component("comp1").physics("pdePm").feature("dirPLeft").set("r", "pLeft");
        model.component("comp1").physics("pdePm").create("dirPRight", "DirichletBoundary", 1);
        model.component("comp1").physics("pdePm").feature("dirPRight").selection().named("rightBnd");
        model.component("comp1").physics("pdePm").feature("dirPRight").set("r", "pRight");
        if (includeFractureCoupling) {
            createWeakBoundaryContribution(model, "pdePm", "weakFracPm", "-test(pm)*(qmf_up + qmf_down)");
        }
    }

    private static void createFracturePressurePde(Model model) {
        model.component("comp1").physics().create("cbPf", "CoefficientFormBoundaryPDE", "geom1");
        model.component("comp1").physics("cbPf").selection().named("fractureBnd");
        model.component("comp1").physics("cbPf").field("dimensionless").field("pf");
        model.component("comp1").physics("cbPf").field("dimensionless").component(new String[]{"pf"});
        model.component("comp1").physics("cbPf").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("cbPf").feature("cfeq1").set("da", "fracture_storage_fun(pf)");
        model.component("comp1").physics("cbPf").feature("cfeq1").set("c", "fracture_mob_fun(pf)");
        model.component("comp1").physics("cbPf").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("cbPf").feature("cfeq1").set("f", "qmf_up + qmf_down");
        model.component("comp1").physics("cbPf").feature("init1").set("pf", "pInit");
    }

    private static void createLegacyPressureOnlyPde(Model model) {
        model.component("comp1").physics().create("c", "CoefficientFormPDE", "geom1");
        model.component("comp1").physics("c").field("dimensionless").field("p");
        model.component("comp1").physics("c").field("dimensionless").component(new String[]{"p"});
        model.component("comp1").physics("c").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("c").feature("cfeq1").set("da", "S_fun_const(p)");
        model.component("comp1").physics("c").feature("cfeq1").set("c", "lambda_fun_const(p)");
        model.component("comp1").physics("c").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("c").feature("cfeq1").set("f", "0");
        model.component("comp1").physics("c").feature("init1").set("p", "pInit");
        model.component("comp1").physics("c").create("dirLeft", "DirichletBoundary", 1);
        model.component("comp1").physics("c").feature("dirLeft").selection().named("leftBnd");
        model.component("comp1").physics("c").feature("dirLeft").set("r", "pLeft");
        model.component("comp1").physics("c").create("dirRight", "DirichletBoundary", 1);
        model.component("comp1").physics("c").feature("dirRight").selection().named("rightBnd");
        model.component("comp1").physics("c").feature("dirRight").set("r", "pRight");
        createExplicitZeroFluxBoundaryFeature(model, "c", "fluxBottom", "bottomBnd");
        createExplicitZeroFluxBoundaryFeature(model, "c", "fluxTop", "topBnd");
    }

    private static void createExplicitZeroFluxBoundaryFeature(
        Model model,
        String physicsTag,
        String featureTag,
        String selectionTag
    ) {
        final String[] candidateTypes = new String[]{"FluxSource", "FluxBoundary", "ZeroFluxBoundary", "ZeroFlux"};
        RuntimeException lastFailure = null;

        for (String candidateType : candidateTypes) {
            try {
                try {
                    model.component("comp1").physics(physicsTag).feature().remove(featureTag);
                } catch (RuntimeException ignored) {
                }

                model.component("comp1").physics(physicsTag).create(featureTag, candidateType, 1);
                model.component("comp1").physics(physicsTag).feature(featureTag).selection().named(selectionTag);

                final String[] propertyNames = model.component("comp1").physics(physicsTag).feature(featureTag).properties();
                final boolean hasQ = arrayContains(propertyNames, "q");
                final boolean hasG = arrayContains(propertyNames, "g");
                if (!hasQ && !hasG) {
                    throw new IllegalStateException(
                        "Feature type '" + candidateType + "' for '" + featureTag +
                        "' does not expose q/g properties. properties=" + Arrays.toString(propertyNames)
                    );
                }
                if (hasQ) {
                    model.component("comp1").physics(physicsTag).feature(featureTag).set("q", "0");
                }
                if (hasG) {
                    model.component("comp1").physics(physicsTag).feature(featureTag).set("g", "0");
                }

                System.out.println(
                    "[COMSOL] explicit_zero_flux " +
                    "physics=" + physicsTag +
                    ", feature=" + featureTag +
                    ", selection=" + selectionTag +
                    ", type=" + candidateType +
                    ", properties=" + Arrays.toString(propertyNames)
                );
                return;
            } catch (RuntimeException ex) {
                lastFailure = ex;
                try {
                    model.component("comp1").physics(physicsTag).feature().remove(featureTag);
                } catch (RuntimeException ignored) {
                }
            }
        }

        throw new IllegalStateException(
            "Unable to create explicit zero-flux boundary feature for physics '" + physicsTag +
            "', selection '" + selectionTag + "'.",
            lastFailure
        );
    }

    private static boolean arrayContains(String[] values, String target) {
        if (values == null) {
            return false;
        }
        for (String value : values) {
            if (target.equals(value)) {
                return true;
            }
        }
        return false;
    }

    private static void createLegacyPressureOnlyFunctions(Model model, LegacyPressureOnlyTables tables) {
        createInterpolationFunction(model, "lambda_fun_const", tables.lambdaTable, "Pa", "kg/(m*s*Pa)");
        createInterpolationFunction(model, "S_fun_const", tables.storageTable, "Pa", "kg/(m^3*Pa)");
    }

    private static void createFullPressureCoefficientFunctions(Model model, FullPressureCoefficientTables tables) {
        createInterpolationFunction(model, "matrix_mob_fun", tables.matrixMobTable, "Pa", "kg/(m*s*Pa)");
        createInterpolationFunction(model, "matrix_storage_fun", tables.matrixStorageTable, "Pa", "kg/(m^3*Pa)");
        createInterpolationFunction(model, "fracture_mob_fun", tables.fractureMobTable, "Pa", "kg/(s*Pa)");
        createInterpolationFunction(model, "fracture_storage_fun", tables.fractureStorageTable, "Pa", "kg/(m^2*Pa)");
    }

    private static void createInterpolationFunction(
        Model model,
        String tag,
        Path filePath,
        String argUnit,
        String funUnit
    ) {
        model.func().create(tag, "Interpolation");
        model.func(tag).set("funcname", tag);
        model.func(tag).set("source", "file");
        model.func(tag).set("filename", filePath.toString());
        model.func(tag).set("interp", "cubicspline");
        model.func(tag).set("extrap", "linear");
        model.func(tag).set("argunit", new String[]{argUnit});
        model.func(tag).set("fununit", funUnit);
    }

    private static FullPressureCoefficientTables prepareFullPressureCoefficientTables(
        CasePaths paths,
        CaseConfig cfg,
        String prefix
    ) {
        final FullPressureCoefficientTables tables = new FullPressureCoefficientTables();
        tables.matrixMobTable = paths.comsolOutputDir.resolve(prefix + "_matrix_mob_fun.txt");
        tables.matrixStorageTable = paths.comsolOutputDir.resolve(prefix + "_matrix_storage_fun.txt");
        tables.fractureMobTable = paths.comsolOutputDir.resolve(prefix + "_fracture_mob_fun.txt");
        tables.fractureStorageTable = paths.comsolOutputDir.resolve(prefix + "_fracture_storage_fun.txt");

        final double pMin = Math.max(0.0, Math.min(cfg.pRightPa, cfg.pLeftPa) - 5.0e5);
        final double pMax = Math.max(cfg.pRightPa, cfg.pLeftPa) + 5.0e5;
        final double dp = 2.5e3;

        try {
            writeConstantInterpolationTable(tables.matrixMobTable, pMin, pMax, dp, cfg.co2Rho * cfg.matrixPerm / cfg.co2Mu);
            writeConstantInterpolationTable(tables.matrixStorageTable, pMin, pMax, dp, cfg.co2Rho * cfg.matrixPhi * cfg.matrixCt);
            writeConstantInterpolationTable(tables.fractureMobTable, pMin, pMax, dp, cfg.co2Rho * cfg.fractureWidthM * cfg.fractureKt / cfg.co2Mu);
            writeConstantInterpolationTable(tables.fractureStorageTable, pMin, pMax, dp, cfg.co2Rho * cfg.fractureWidthM * cfg.fracturePhi * cfg.fractureCt);
        } catch (IOException ex) {
            throw new IllegalStateException("Failed to write full-model pressure interpolation tables.", ex);
        }

        return tables;
    }

    private static LegacyPressureOnlyTables prepareLegacyPressureOnlyTables(
        CasePaths paths,
        CaseConfig cfg,
        String prefix
    ) {
        final LegacyPressureOnlyTables tables = new LegacyPressureOnlyTables();
        tables.lambdaTable = paths.comsolOutputDir.resolve(prefix + "_lambda_fun.txt");
        tables.storageTable = paths.comsolOutputDir.resolve(prefix + "_S_fun.txt");

        final double pMin = Math.max(0.0, Math.min(cfg.pRightPa, cfg.pLeftPa) - 5.0e5);
        final double pMax = Math.max(cfg.pRightPa, cfg.pLeftPa) + 5.0e5;
        final double dp = 2.5e3;

        try {
            writeConstantInterpolationTable(tables.lambdaTable, pMin, pMax, dp, cfg.co2Rho * cfg.matrixPerm / cfg.co2Mu);
            writeConstantInterpolationTable(tables.storageTable, pMin, pMax, dp, cfg.co2Rho * cfg.matrixPhi * cfg.matrixCt);
        } catch (IOException ex) {
            throw new IllegalStateException("Failed to write legacy pressure-only interpolation tables.", ex);
        }

        return tables;
    }

    private static void writeConstantInterpolationTable(
        Path tablePath,
        double pMinPa,
        double pMaxPa,
        double dpPa,
        double value
    ) throws IOException {
        try (BufferedWriter out = Files.newBufferedWriter(tablePath, StandardCharsets.UTF_8)) {
            double p = pMinPa;
            while (p <= pMaxPa + 0.5 * dpPa) {
                out.write(String.format(Locale.US, "%.0f %.16E", p, value));
                out.newLine();
                p += dpPa;
            }
        }
    }

    private static void createTransplantPressureOnlyFunctions(Model model, LegacyPressureOnlyTables tables) {
        createInterpolationFunction(model, "lambda_fun", tables.lambdaTable, "Pa", "kg/(m*s*Pa)");
        createInterpolationFunction(model, "S_fun", tables.storageTable, "Pa", "kg/(m^3*Pa)");
    }

    private static void createTransplantPressureOnlyPde(Model model) {
        model.component("comp1").physics().create("c", "CoefficientFormPDE", "geom1");
        model.component("comp1").physics("c").field("dimensionless").field("p");
        model.component("comp1").physics("c").field("dimensionless").component(new String[]{"p"});
        model.component("comp1").physics("c").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("c").feature("cfeq1").set("da", "S_fun(p)");
        model.component("comp1").physics("c").feature("cfeq1").set("c", "lambda_fun(p)");
        model.component("comp1").physics("c").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("c").feature("cfeq1").set("f", "0");
        model.component("comp1").physics("c").feature("init1").set("p", "pInit");
        model.component("comp1").physics("c").create("dirLeft", "DirichletBoundary", 1);
        model.component("comp1").physics("c").feature("dirLeft").selection().named("leftBnd");
        model.component("comp1").physics("c").feature("dirLeft").set("r", "pLeft");
        model.component("comp1").physics("c").create("dirRight", "DirichletBoundary", 1);
        model.component("comp1").physics("c").feature("dirRight").selection().named("rightBnd");
        model.component("comp1").physics("c").feature("dirRight").set("r", "pRight");
    }

    private static void createTransplantPressureOnlyStationaryPde(Model model) {
        model.component("comp1").physics().create("c", "CoefficientFormPDE", "geom1");
        model.component("comp1").physics("c").field("dimensionless").field("p");
        model.component("comp1").physics("c").field("dimensionless").component(new String[]{"p"});
        model.component("comp1").physics("c").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("c").feature("cfeq1").set("da", "0");
        model.component("comp1").physics("c").feature("cfeq1").set("c", "matrix_mob");
        model.component("comp1").physics("c").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("c").feature("cfeq1").set("f", "0");
        model.component("comp1").physics("c").feature("init1").set("p", "pInit");
        model.component("comp1").physics("c").create("dirLeft", "DirichletBoundary", 1);
        model.component("comp1").physics("c").feature("dirLeft").selection().named("leftBnd");
        model.component("comp1").physics("c").feature("dirLeft").set("r", "pLeft");
        model.component("comp1").physics("c").create("dirRight", "DirichletBoundary", 1);
        model.component("comp1").physics("c").feature("dirRight").selection().named("rightBnd");
        model.component("comp1").physics("c").feature("dirRight").set("r", "pRight");
    }

    private static void createMatrixTemperaturePdes(Model model) {
        createMatrixTemperaturePdes(model, true);
    }

    private static void createMatrixTemperaturePdes(Model model, boolean includeFractureCoupling) {
        model.component("comp1").physics().create("pdeTm", "CoefficientFormPDE", "geom1");
        model.component("comp1").physics("pdeTm").field("dimensionless").field("Tm");
        model.component("comp1").physics("pdeTm").field("dimensionless").component(new String[]{"Tm"});
        model.component("comp1").physics("pdeTm").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("pdeTm").feature("cfeq1").set("da", "rhoCp_matrix");
        model.component("comp1").physics("pdeTm").feature("cfeq1").set("c", "ktherm_matrix");
        model.component("comp1").physics("pdeTm").feature("cfeq1").set("be", new String[]{"betaTm_x", "betaTm_y"});
        model.component("comp1").physics("pdeTm").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("pdeTm").feature("cfeq1").set("f", "0");
        model.component("comp1").physics("pdeTm").feature("init1").set("Tm", "TInit");
        model.component("comp1").physics("pdeTm").create("dirTLeft", "DirichletBoundary", 1);
        model.component("comp1").physics("pdeTm").feature("dirTLeft").selection().named("leftBnd");
        model.component("comp1").physics("pdeTm").feature("dirTLeft").set("r", "TLeft");
        model.component("comp1").physics("pdeTm").create("dirTRight", "DirichletBoundary", 1);
        model.component("comp1").physics("pdeTm").feature("dirTRight").selection().named("rightBnd");
        model.component("comp1").physics("pdeTm").feature("dirTRight").set("r", "TRight");
        if (includeFractureCoupling) {
            createWeakBoundaryContribution(model, "pdeTm", "weakFracTm", "-test(Tm)*(Hmf_up + Hmf_down)");
        }
    }

    private static void createFractureTemperaturePde(Model model) {
        model.component("comp1").physics().create("cbTf", "CoefficientFormBoundaryPDE", "geom1");
        model.component("comp1").physics("cbTf").selection().named("fractureBnd");
        model.component("comp1").physics("cbTf").field("dimensionless").field("Tf");
        model.component("comp1").physics("cbTf").field("dimensionless").component(new String[]{"Tf"});
        model.component("comp1").physics("cbTf").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("cbTf").feature("cfeq1").set("da", "wFrac*rhoCp_fracture");
        model.component("comp1").physics("cbTf").feature("cfeq1").set("c", "wFrac*ktherm_fracture");
        model.component("comp1").physics("cbTf").feature("cfeq1").set("be", new String[]{"betaTf_x", "betaTf_y"});
        model.component("comp1").physics("cbTf").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("cbTf").feature("cfeq1").set("f", "Hmf_up + Hmf_down");
        model.component("comp1").physics("cbTf").feature("init1").set("Tf", "TInit");
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

    private static void configureLegacyStyleMesh(Model model, double hmaxM, double hminM) {
        model.component("comp1").mesh("mesh1").feature().create("size1", "Size");
        model.component("comp1").mesh("mesh1").feature("size1").set("custom", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hmaxactive", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hminactive", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hmax", formatQuantity(hmaxM, "m"));
        model.component("comp1").mesh("mesh1").feature("size1").set("hmin", formatQuantity(hminM, "m"));
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

    private static void createWeakBoundaryContribution(Model model, String physicsTag, String featureTag, String weakExpr) {
        try {
            try {
                model.component("comp1").common().remove(featureTag);
            } catch (RuntimeException ignored) {
            }
            model.component("comp1").common().create(featureTag, "WeakContribution");
            model.component("comp1").common(featureTag).selection().named("fractureBnd");
            model.component("comp1").common(featureTag).set("weakExpression", weakExpr);
            model.component("comp1").common(featureTag).set("integrationOrder", "4");
            model.component("comp1").common(featureTag).label("WeakCoupling_" + physicsTag);
        } catch (RuntimeException ex) {
            throw new IllegalStateException(
                "Unable to create component-level weak boundary contribution for physics '" + physicsTag +
                "' on selection 'fractureBnd' using property 'weakExpression'.",
                ex
            );
        }
    }

    private static void configureStudy(Model model, double[] times) {
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").activate("pdePm", true);
        model.study("std1").feature("time").activate("pdeTm", true);
        try {
            model.study("std1").feature("time").activate("cbPf", true);
        } catch (RuntimeException ignored) {
        }
        try {
            model.study("std1").feature("time").activate("cbTf", true);
        } catch (RuntimeException ignored) {
        }
        model.study("std1").feature("time").set("tlist", formatTimeList(times));
    }

    private static void configurePressureOnlyStudy(Model model, double[] times, boolean includeFractureFields) {
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").activate("pdePm", true);
        if (includeFractureFields) {
            try {
                model.study("std1").feature("time").activate("cbPf", true);
            } catch (RuntimeException ignored) {
            }
        }
        try {
            model.study("std1").feature("time").activate("pdeTm", false);
        } catch (RuntimeException ignored) {
        }
        try {
            model.study("std1").feature("time").activate("cbTf", false);
        } catch (RuntimeException ignored) {
        }
        model.study("std1").feature("time").set("tlist", formatTimeList(times));
    }

    private static void configureLegacyPressureOnlyStudy(Model model, double[] times) {
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").activate("c", true);
        model.study("std1").feature("time").set("tlist", formatTimeList(times));
    }

    private static void configureLegacyPressureOnlyUniformStudy(Model model, double endTimeS, double dtS) {
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").activate("c", true);
        model.study("std1").feature("time").set(
            "tlist",
            String.format(Locale.US, "range(0,%.12f,%.12f)", Math.max(1.0, dtS), endTimeS)
        );
    }

    private static void configureStationaryPressureOnlyStudy(Model model) {
        model.study().create("std1");
        model.study("std1").create("stat", "Stationary");
        model.study("std1").feature("stat").activate("c", true);
    }

    private static void configureSegregatedTimeSolver(Model model, double[] times, boolean includeFractureFields) {
        model.study("std1").createAutoSequences("sol");
        model.sol("sol1").study("std1");
        model.sol("sol1").attach("std1");
        model.sol("sol1").feature("st1").set("study", "std1");
        model.sol("sol1").feature("st1").set("studystep", "time");
        model.sol("sol1").feature("v1").set("control", "time");
        model.sol("sol1").feature("t1").set("control", "time");
        model.sol("sol1").feature("t1").set("tlist", formatTimeList(times));
        model.sol("sol1").feature("t1").set("plot", "off");
        model.sol("sol1").feature("t1").set("probesel", "all");
        model.sol("sol1").feature("t1").set("probefreq", "tsteps");
        model.sol("sol1").feature("t1").set("rtol", "1e-5");
        model.sol("sol1").feature("t1").set("maxorder", "2");

        try {
            model.sol("sol1").feature("t1").feature().remove("fc1");
        } catch (RuntimeException ignored) {
        }
        try {
            model.sol("sol1").feature("t1").feature().remove("se1");
        } catch (RuntimeException ignored) {
        }

        model.sol("sol1").feature("t1").create("se1", "Segregated");
        model.sol("sol1").feature("t1").feature("se1").set("maxsegiter", 12);
        model.sol("sol1").feature("t1").feature("se1").set("ntolfact", 1.0);

        final String[] pressureVars = includeFractureFields
            ? new String[]{"comp1_pm", "comp1_pf"}
            : new String[]{"comp1_pm"};
        final String[] temperatureVars = includeFractureFields
            ? new String[]{"comp1_Tm", "comp1_Tf"}
            : new String[]{"comp1_Tm"};

        createSegregatedStep(model, "segP", pressureVars);
        createSegregatedStep(model, "segT", temperatureVars);
    }

    private static void createSegregatedStep(Model model, String stepTag, String[] variables) {
        model.sol("sol1").feature("t1").feature("se1").create(stepTag, "SegregatedStep");
        model.sol("sol1").feature("t1").feature("se1").feature(stepTag).set("segvar", variables);
        model.sol("sol1").feature("t1").feature("se1").feature(stepTag).set("linsolver", "dDef");
        model.sol("sol1").feature("t1").feature("se1").feature(stepTag).set("subdtech", "const");
        model.sol("sol1").feature("t1").feature("se1").feature(stepTag).set("subiter", 1);
        model.sol("sol1").feature("t1").feature("se1").feature(stepTag).set("maxsubiter", 8);
    }

    private static void writeFractureExchangeDiagnostics(
        Model model,
        CasePaths paths,
        CaseConfig cfg,
        double finalTimeS
    ) throws IOException {
        final int sampleCount = 11;
        final double[][] coords = new double[2][sampleCount];
        final double x0 = cfg.fracX0Ratio * cfg.lxM;
        final double y0 = cfg.fracY0Ratio * cfg.lyM;
        final double x1 = cfg.fracX1Ratio * cfg.lxM;
        final double y1 = cfg.fracY1Ratio * cfg.lyM;
        for (int i = 0; i < sampleCount; ++i) {
            final double ratio = (double) i / (double) (sampleCount - 1);
            coords[0][i] = x0 + ratio * (x1 - x0);
            coords[1][i] = y0 + ratio * (y1 - y0);
        }

        final double[] times = new double[]{finalTimeS};
        final double[][] xMapUp = evaluateAtPoints(model, "diagXMapUp", "xMapUp", "m", coords, times, 1, "fractureBnd");
        final double[][] yMapUp = evaluateAtPoints(model, "diagYMapUp", "yMapUp", "m", coords, times, 1, "fractureBnd");
        final double[][] xMapDown = evaluateAtPoints(model, "diagXMapDown", "xMapDown", "m", coords, times, 1, "fractureBnd");
        final double[][] yMapDown = evaluateAtPoints(model, "diagYMapDown", "yMapDown", "m", coords, times, 1, "fractureBnd");
        final double[][] pmUp = evaluateAtPoints(model, "diagPmUp", "pmUpTrace", "Pa", coords, times, 1, "fractureBnd");
        final double[][] pmDown = evaluateAtPoints(model, "diagPmDown", "pmDownTrace", "Pa", coords, times, 1, "fractureBnd");
        final double[][] pf = evaluateAtPoints(model, "diagPf", "pf", "Pa", coords, times, 1, "fractureBnd");
        final double[][] qmfUp = evaluateAtPoints(model, "diagQmfUp", "qmf_up", "kg/(m^2*s)", coords, times, 1, "fractureBnd");
        final double[][] qmfDown = evaluateAtPoints(model, "diagQmfDown", "qmf_down", "kg/(m^2*s)", coords, times, 1, "fractureBnd");
        final double[][] qmfTotal = evaluateAtPoints(model, "diagQmfTotal", "qmfTotal", "kg/(m^2*s)", coords, times, 1, "fractureBnd");
        final double[][] tmUp = evaluateAtPoints(model, "diagTmUp", "TmUpTrace", "K", coords, times, 1, "fractureBnd");
        final double[][] tmDown = evaluateAtPoints(model, "diagTmDown", "TmDownTrace", "K", coords, times, 1, "fractureBnd");
        final double[][] tf = evaluateAtPoints(model, "diagTf", "Tf", "K", coords, times, 1, "fractureBnd");
        final double[][] hmfUp = evaluateAtPoints(model, "diagHmfUp", "Hmf_up", "W/m^2", coords, times, 1, "fractureBnd");
        final double[][] hmfDown = evaluateAtPoints(model, "diagHmfDown", "Hmf_down", "W/m^2", coords, times, 1, "fractureBnd");
        final double[][] hmfTotal = evaluateAtPoints(model, "diagHmfTotal", "HmfTotal", "W/m^2", coords, times, 1, "fractureBnd");

        final double fractureLength = Math.sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
        final Path csvPath = paths.comsolOutputDir.resolve("comsol_fracture_exchange_diagnostics.csv");
        try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
            out.write("sample_id,s_m,x_m,y_m,x_map_up_m,y_map_up_m,x_map_down_m,y_map_down_m,pm_up_pa,pm_down_pa,pf_pa,qmf_up,qmf_down,qmf_total,tm_up_k,tm_down_k,tf_k,hmf_up,hmf_down,hmf_total");
            out.newLine();
            for (int i = 0; i < sampleCount; ++i) {
                final double s = fractureLength * (double) i / (double) (sampleCount - 1);
                out.write(String.format(
                    Locale.US,
                    "%d,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12e,%.12e,%.12e,%.12f,%.12f,%.12f,%.12e,%.12e,%.12e",
                    i,
                    s,
                    coords[0][i],
                    coords[1][i],
                    xMapUp[0][i],
                    yMapUp[0][i],
                    xMapDown[0][i],
                    yMapDown[0][i],
                    pmUp[0][i],
                    pmDown[0][i],
                    pf[0][i],
                    qmfUp[0][i],
                    qmfDown[0][i],
                    qmfTotal[0][i],
                    tmUp[0][i],
                    tmDown[0][i],
                    tf[0][i],
                    hmfUp[0][i],
                    hmfDown[0][i],
                    hmfTotal[0][i]
                ));
                out.newLine();
            }
        }

        System.out.println(
            "[COMSOL] fracture_exchange_diag " +
            "pm_up_span=" + formatDouble(span(pmUp[0])) +
            ", pm_down_span=" + formatDouble(span(pmDown[0])) +
            ", pf_span=" + formatDouble(span(pf[0])) +
            ", qmf_total_maxabs=" + formatDouble(maxAbs(qmfTotal[0])) +
            ", tf_span=" + formatDouble(span(tf[0])) +
            ", hmf_total_maxabs=" + formatDouble(maxAbs(hmfTotal[0]))
        );
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

            final double[][] pValues = evaluateProfileFamily(model, familyStations, times, true, "interpP_" + family);
            final double[][] tValues = evaluateProfileFamily(model, familyStations, times, false, "interpT_" + family);
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

    private static void writeMatrixOnlyDiagnostic(
        Model model,
        CasePaths paths,
        List<ProfileStation> profileStations,
        List<TimeRequest> profileTimes
    ) throws IOException {
        final List<ProfileStation> familyStations = filterProfileStations(profileStations, "matrix_horizontal");
        if (familyStations.isEmpty()) {
            return;
        }

        final double[] times = new double[profileTimes.size()];
        for (int i = 0; i < profileTimes.size(); ++i) {
            times[i] = profileTimes.get(i).timeS;
        }

        final double[][] coords = buildCoordinateMatrix(familyStations);
        final double[][] pValues = evaluateAtPoints(model, "diagMatrixOnlyProfileP", "pMatrixPlot", "Pa", coords, times, 2, null);
        final double[][] tValues = evaluateAtPoints(model, "diagMatrixOnlyProfileT", "TMatrixPlot", "K", coords, times, 2, null);

        for (int timeIndex = 0; timeIndex < profileTimes.size(); ++timeIndex) {
            final TimeRequest request = profileTimes.get(timeIndex);
            final Path csvPath = paths.comsolOutputDir.resolve("comsol_matrix_only_profile_matrix_horizontal_" + request.tag + ".csv");
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

        final int finalIndex = profileTimes.size() - 1;
        System.out.println(
            "[COMSOL] matrix_only_diag " +
            "p_span=" + formatDouble(span(pValues[finalIndex])) +
            ", t_span=" + formatDouble(span(tValues[finalIndex]))
        );
    }

    private static void writeMatrixOnlyPressureOnlyDiagnostic(
        Model model,
        CasePaths paths,
        List<ProfileStation> profileStations,
        List<TimeRequest> profileTimes,
        CaseConfig cfg,
        double[] solutionTimes
    ) throws IOException {
        writeMatrixOnlyPressureOnlyDiagnostic(
            model,
            paths,
            profileStations,
            profileTimes,
            cfg,
            solutionTimes,
            "comsol_matrix_only_pressureonly",
            "matrix_only_pressureonly_diag",
            "diagMatrixOnlyPressureOnly"
        );
    }

    private static void writeMatrixOnlyPressureOnlyDiagnostic(
        Model model,
        CasePaths paths,
        List<ProfileStation> profileStations,
        List<TimeRequest> profileTimes,
        CaseConfig cfg,
        double[] solutionTimes,
        String filePrefix,
        String logLabel,
        String tagBasePrefix
    ) throws IOException {
        final List<ProfileStation> familyStations = filterProfileStations(profileStations, "matrix_horizontal");
        if (familyStations.isEmpty()) {
            return;
        }

        final double[] times = new double[profileTimes.size()];
        for (int i = 0; i < profileTimes.size(); ++i) {
            times[i] = profileTimes.get(i).timeS;
        }

        final double[][] coords = buildCoordinateMatrix(familyStations);
        final double offsetY = Math.min(cfg.lyM - 1.0e-6, Math.max(1.0e-6, 0.5 * cfg.lyM + 0.137));
        final double offsetX = Math.min(cfg.lxM - 1.0e-6, Math.max(1.0e-6, 0.137));
        final double[][] offsetCoords = buildYOffsetCoordinateMatrix(familyStations, offsetY);
        final double[][] xyOffsetCoords = buildXYOffsetCoordinateMatrix(familyStations, offsetX, offsetY, cfg.lxM);
        final String explicitDatasetTag = tagBasePrefix + "ExplicitDset";
        ensureExplicitSolutionDataset(model, explicitDatasetTag);
        final double[][] pValues = evaluateAtPoints(model, tagBasePrefix + "ProfileP", "pMatrixPlot", "Pa", coords, times, 2, null);
        final double[][] pValuesExplicitDataset = evaluateAtPoints(
            model,
            tagBasePrefix + "ProfilePExplicitDset",
            "pMatrixPlot",
            "Pa",
            coords,
            times,
            2,
            null,
            explicitDatasetTag
        );
        final int[] solutionNumbers = mapTimesToSolutionNumbers(solutionTimes, times);
        final double[][] pValuesBySolnum = evaluateAtPointsBySolutionNumber(
            model,
            tagBasePrefix + "ProfilePSolnum",
            "pMatrixPlot",
            "Pa",
            coords,
            solutionNumbers,
            2,
            null
        );
        final double[][] pValuesOffset = evaluateAtPoints(model, tagBasePrefix + "ProfilePShifted", "pMatrixPlot", "Pa", offsetCoords, times, 2, null);
        final double[][] pValuesXYOffset = evaluateAtPoints(model, tagBasePrefix + "ProfilePShiftedXY", "pMatrixPlot", "Pa", xyOffsetCoords, times, 2, null);

        for (int timeIndex = 0; timeIndex < profileTimes.size(); ++timeIndex) {
            final TimeRequest request = profileTimes.get(timeIndex);
            final Path csvPath = paths.comsolOutputDir.resolve(filePrefix + "_profile_matrix_horizontal_" + request.tag + ".csv");
            try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
                out.write("station_id,label,family,location,target_axis_m,target_x_m,target_y_m,target_time_s,p_ref_pa,p_linear_pa,p_error_pa");
                out.newLine();
                for (int pointIndex = 0; pointIndex < familyStations.size(); ++pointIndex) {
                    final ProfileStation point = familyStations.get(pointIndex);
                    final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (point.targetX / cfg.lxM);
                    final double pRef = pValues[timeIndex][pointIndex];
                    out.write(String.format(
                        Locale.US,
                        "%d,%s,%s,%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                        point.id, point.label, point.family, point.location,
                        point.targetAxis, point.targetX, point.targetY, request.timeS,
                        pRef, pLinear, pRef - pLinear
                    ));
                    out.newLine();
                }
            }

            final Path explicitDatasetCsvPath = paths.comsolOutputDir.resolve(
                filePrefix + "_profile_matrix_horizontal_explicitdset_" + request.tag + ".csv"
            );
            try (BufferedWriter out = Files.newBufferedWriter(explicitDatasetCsvPath, StandardCharsets.UTF_8)) {
                out.write("station_id,label,family,location,target_axis_m,target_x_m,target_y_m,target_time_s,p_ref_pa,p_linear_pa,p_error_pa");
                out.newLine();
                for (int pointIndex = 0; pointIndex < familyStations.size(); ++pointIndex) {
                    final ProfileStation point = familyStations.get(pointIndex);
                    final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (point.targetX / cfg.lxM);
                    final double pRef = pValuesExplicitDataset[timeIndex][pointIndex];
                    out.write(String.format(
                        Locale.US,
                        "%d,%s,%s,%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                        point.id, point.label, point.family, point.location,
                        point.targetAxis, point.targetX, point.targetY, request.timeS,
                        pRef, pLinear, pRef - pLinear
                    ));
                    out.newLine();
                }
            }

            final Path solnumCsvPath = paths.comsolOutputDir.resolve(filePrefix + "_profile_matrix_horizontal_solnum_" + request.tag + ".csv");
            try (BufferedWriter out = Files.newBufferedWriter(solnumCsvPath, StandardCharsets.UTF_8)) {
                out.write(
                    "station_id,label,family,location,target_axis_m,target_x_m,target_y_m,target_time_s," +
                    "matched_solution_index,matched_solution_time_s,p_ref_pa,p_linear_pa,p_error_pa"
                );
                out.newLine();
                for (int pointIndex = 0; pointIndex < familyStations.size(); ++pointIndex) {
                    final ProfileStation point = familyStations.get(pointIndex);
                    final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (point.targetX / cfg.lxM);
                    final double pRef = pValuesBySolnum[timeIndex][pointIndex];
                    out.write(String.format(
                        Locale.US,
                        "%d,%s,%s,%s,%.12f,%.12f,%.12f,%.12f,%d,%.12f,%.12f,%.12f,%.12f",
                        point.id, point.label, point.family, point.location,
                        point.targetAxis, point.targetX, point.targetY, request.timeS,
                        solutionNumbers[timeIndex], solutionTimes[solutionNumbers[timeIndex] - 1],
                        pRef, pLinear, pRef - pLinear
                    ));
                    out.newLine();
                }
            }

            final Path offsetCsvPath = paths.comsolOutputDir.resolve(filePrefix + "_profile_matrix_horizontal_yoffset_" + request.tag + ".csv");
            try (BufferedWriter out = Files.newBufferedWriter(offsetCsvPath, StandardCharsets.UTF_8)) {
                out.write("station_id,label,family,location,target_axis_m,target_x_m,target_y_m,target_time_s,p_ref_pa,p_linear_pa,p_error_pa");
                out.newLine();
                for (int pointIndex = 0; pointIndex < familyStations.size(); ++pointIndex) {
                    final ProfileStation point = familyStations.get(pointIndex);
                    final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (point.targetX / cfg.lxM);
                    final double pRef = pValuesOffset[timeIndex][pointIndex];
                    out.write(String.format(
                        Locale.US,
                        "%d,%s,%s,%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                        point.id, point.label, point.family, point.location,
                        point.targetAxis, point.targetX, offsetY, request.timeS,
                        pRef, pLinear, pRef - pLinear
                    ));
                    out.newLine();
                }
            }

            final Path xyOffsetCsvPath = paths.comsolOutputDir.resolve(filePrefix + "_profile_matrix_horizontal_xyoffset_" + request.tag + ".csv");
            try (BufferedWriter out = Files.newBufferedWriter(xyOffsetCsvPath, StandardCharsets.UTF_8)) {
                out.write("station_id,label,family,location,target_axis_m,target_x_m,target_y_m,target_time_s,p_ref_pa,p_linear_pa,p_error_pa");
                out.newLine();
                for (int pointIndex = 0; pointIndex < familyStations.size(); ++pointIndex) {
                    final ProfileStation point = familyStations.get(pointIndex);
                    final double shiftedX = Math.min(cfg.lxM - 1.0e-6, point.targetX + offsetX);
                    final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (shiftedX / cfg.lxM);
                    final double pRef = pValuesXYOffset[timeIndex][pointIndex];
                    out.write(String.format(
                        Locale.US,
                        "%d,%s,%s,%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                        point.id, point.label, point.family, point.location,
                        point.targetAxis, shiftedX, offsetY, request.timeS,
                        pRef, pLinear, pRef - pLinear
                    ));
                    out.newLine();
                }
            }
        }

        final int finalIndex = profileTimes.size() - 1;
        final double exactL2 = computeNormalizedLinearProfileL2(pValues[finalIndex], familyStations, cfg);
        final double explicitDatasetL2 = computeNormalizedLinearProfileL2(pValuesExplicitDataset[finalIndex], familyStations, cfg);
        final double solnumL2 = computeNormalizedLinearProfileL2(pValuesBySolnum[finalIndex], familyStations, cfg);
        final double offsetL2 = computeNormalizedLinearProfileL2(pValuesOffset[finalIndex], familyStations, cfg);
        final double xyOffsetL2 = computeNormalizedLinearProfileL2WithXOffset(pValuesXYOffset[finalIndex], familyStations, cfg, offsetX);
        final double timeVsSolnumMaxAbs = computeMaxAbsDifference(pValues[finalIndex], pValuesBySolnum[finalIndex]);
        final double timeVsExplicitDatasetMaxAbs = computeMaxAbsDifference(pValues[finalIndex], pValuesExplicitDataset[finalIndex]);
        System.out.println(
            "[COMSOL] " + logLabel + " " +
            "p_span=" + formatDouble(span(pValues[finalIndex])) +
            ", p_l2_norm=" + formatDouble(exactL2) +
            ", explicit_dset_l2_norm=" + formatDouble(explicitDatasetL2) +
            ", solnum_l2_norm=" + formatDouble(solnumL2) +
            ", time_vs_solnum_maxabs_pa=" + formatDouble(timeVsSolnumMaxAbs) +
            ", time_vs_explicit_dset_maxabs_pa=" + formatDouble(timeVsExplicitDatasetMaxAbs) +
            ", xoffset=" + formatDouble(offsetX) +
            ", yoffset=" + formatDouble(offsetY) +
            ", yoffset_l2_norm=" + formatDouble(offsetL2) +
            ", xyoffset_l2_norm=" + formatDouble(xyOffsetL2)
        );
    }

    private static void writeMatrixOnlyPressureOnlyBoundaryDiagnostics(
        Model model,
        CasePaths paths,
        CaseConfig cfg,
        double finalTimeS
    ) throws IOException {
        writeMatrixOnlyPressureOnlyBoundaryDiagnostics(
            model,
            paths,
            cfg,
            finalTimeS,
            "comsol_matrix_only_pressureonly_boundary_diagnostics.csv",
            "matrix_only_pressureonly_boundary_line",
            "diagBoundaryLine"
        );
    }

    private static void writeMatrixOnlyPressureOnlyStationaryDiagnostic(
        Model model,
        CasePaths paths,
        List<ProfileStation> profileStations,
        CaseConfig cfg,
        String filePrefix,
        String logLabel,
        String tagBasePrefix
    ) throws IOException {
        final List<ProfileStation> familyStations = filterProfileStations(profileStations, "matrix_horizontal");
        if (familyStations.isEmpty()) {
            return;
        }

        final double[][] coords = buildCoordinateMatrix(familyStations);
        final double[][] pValues = evaluateAtPointsBySolutionNumber(
            model,
            tagBasePrefix + "ProfileP",
            "pMatrixPlot",
            "Pa",
            coords,
            new int[]{1},
            2,
            null
        );

        final Path csvPath = paths.comsolOutputDir.resolve(filePrefix + "_profile_matrix_horizontal.csv");
        try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
            out.write("station_id,label,family,location,target_axis_m,target_x_m,target_y_m,p_ref_pa,p_linear_pa,p_error_pa");
            out.newLine();
            for (int pointIndex = 0; pointIndex < familyStations.size(); ++pointIndex) {
                final ProfileStation point = familyStations.get(pointIndex);
                final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (point.targetX / cfg.lxM);
                final double pRef = pValues[0][pointIndex];
                out.write(String.format(
                    Locale.US,
                    "%d,%s,%s,%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                    point.id, point.label, point.family, point.location,
                    point.targetAxis, point.targetX, point.targetY,
                    pRef, pLinear, pRef - pLinear
                ));
                out.newLine();
            }
        }

        final double exactL2 = computeNormalizedLinearProfileL2(pValues[0], familyStations, cfg);
        final double exactLinf = computeNormalizedLinearProfileLinf(pValues[0], familyStations, cfg);
        System.out.println(
            "[COMSOL] " + logLabel + " " +
            "p_span=" + formatDouble(span(pValues[0])) +
            ", p_l2_norm=" + formatDouble(exactL2) +
            ", p_linf_norm=" + formatDouble(exactLinf)
        );
    }

    private static void writeMatrixOnlyPressureOnlyBoundaryDiagnostics(
        Model model,
        CasePaths paths,
        CaseConfig cfg,
        double finalTimeS,
        String fileName,
        String logLabel,
        String tagBase
    ) throws IOException {
        final double[] probeXs = new double[]{
            0.0,
            1.0,
            5.0,
            10.0,
            20.0,
            200.0,
            cfg.lxM - 20.0,
            cfg.lxM - 10.0,
            cfg.lxM - 5.0,
            cfg.lxM - 1.0,
            cfg.lxM
        };
        final int yCount = 41;
        final double[] yValues = new double[yCount];
        for (int i = 0; i < yCount; ++i) {
            yValues[i] = cfg.lyM * (double) i / (double) (yCount - 1);
        }

        final Path csvPath = paths.comsolOutputDir.resolve(fileName);
        try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
            out.write("line_tag,x_m,y_m,target_time_s,p_ref_pa");
            out.newLine();

            for (double probeX : probeXs) {
                final double[][] coords = new double[2][yCount];
                for (int i = 0; i < yCount; ++i) {
                    coords[0][i] = Math.min(cfg.lxM, Math.max(0.0, probeX));
                    coords[1][i] = yValues[i];
                }
                final double[][] pValues = evaluateAtPoints(
                    model,
                    tagBase + "_" + sanitizeTag(String.format(Locale.US, "%.3f", probeX)),
                    "pMatrixPlot",
                    "Pa",
                    coords,
                    new double[]{finalTimeS},
                    2,
                    null
                );

                final String lineTag = "x=" + formatDouble(probeX);
                for (int i = 0; i < yCount; ++i) {
                    out.write(String.format(
                        Locale.US,
                        "%s,%.12f,%.12f,%.12f,%.12f",
                        lineTag,
                        coords[0][i],
                        coords[1][i],
                        finalTimeS,
                        pValues[0][i]
                    ));
                    out.newLine();
                }

                final double mean = average(pValues[0]);
                final double lineSpan = span(pValues[0]);
                final double linearMean = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (probeX / cfg.lxM);
                System.out.println(
                    "[COMSOL] " + logLabel + " " +
                    "x=" + formatDouble(probeX) +
                    ", mean=" + formatDouble(mean) +
                    ", span_y=" + formatDouble(lineSpan) +
                    ", linear_mean=" + formatDouble(linearMean)
                );
            }
        }
    }

    private static void writeMatrixOnlyPressureOnlyHorizontalDiagnostics(
        Model model,
        CasePaths paths,
        CaseConfig cfg,
        double finalTimeS
    ) throws IOException {
        writeMatrixOnlyPressureOnlyHorizontalDiagnostics(
            model,
            paths,
            cfg,
            finalTimeS,
            "comsol_matrix_only_pressureonly_horizontal_diagnostics.csv",
            "matrix_only_pressureonly_horizontal_line",
            "diagHorizontalLine"
        );
    }

    private static void writeMatrixOnlyPressureOnlyHorizontalDiagnostics(
        Model model,
        CasePaths paths,
        CaseConfig cfg,
        double finalTimeS,
        String fileName,
        String logLabel,
        String tagBase
    ) throws IOException {
        final double[] probeYs = new double[]{
            0.0,
            0.01,
            0.1,
            1.0,
            5.0,
            10.0,
            20.0,
            cfg.lyM - 10.0,
            cfg.lyM - 5.0,
            cfg.lyM - 1.0,
            cfg.lyM - 0.1,
            cfg.lyM - 0.01,
            cfg.lyM
        };
        final int xCount = 81;
        final double[] xValues = new double[xCount];
        for (int i = 0; i < xCount; ++i) {
            xValues[i] = cfg.lxM * (double) i / (double) (xCount - 1);
        }

        final Path csvPath = paths.comsolOutputDir.resolve(fileName);
        try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
            out.write("line_tag,x_m,y_m,target_time_s,p_ref_pa,p_linear_pa,p_error_pa");
            out.newLine();

            for (double probeYRaw : probeYs) {
                final double probeY = Math.min(cfg.lyM, Math.max(0.0, probeYRaw));
                final double[][] coords = new double[2][xCount];
                for (int i = 0; i < xCount; ++i) {
                    coords[0][i] = xValues[i];
                    coords[1][i] = probeY;
                }
                final double[][] pValues = evaluateAtPoints(
                    model,
                    tagBase + "_" + sanitizeTag(String.format(Locale.US, "%.3f", probeY)),
                    "pMatrixPlot",
                    "Pa",
                    coords,
                    new double[]{finalTimeS},
                    2,
                    null
                );

                final double[] errors = new double[xCount];
                for (int i = 0; i < xCount; ++i) {
                    final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (xValues[i] / cfg.lxM);
                    errors[i] = pValues[0][i] - pLinear;
                    out.write(String.format(
                        Locale.US,
                        "y=%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                        formatDouble(probeY),
                        xValues[i],
                        probeY,
                        finalTimeS,
                        pValues[0][i],
                        pLinear,
                        errors[i]
                    ));
                    out.newLine();
                }

                final double errorL2 = computeNormalizedL2(errors, cfg.deltaPPa);
                final double errorLinf = computeNormalizedLinf(errors, cfg.deltaPPa);
                System.out.println(
                    "[COMSOL] " + logLabel + " " +
                    "y=" + formatDouble(probeY) +
                    ", l2_norm=" + formatDouble(errorL2) +
                    ", linf_norm=" + formatDouble(errorLinf) +
                    ", span_x=" + formatDouble(span(pValues[0]))
                );
            }
        }
    }

    private static double[][] buildYOffsetCoordinateMatrix(List<ProfileStation> points, double yValue) {
        final double[][] coords = new double[2][points.size()];
        for (int i = 0; i < points.size(); ++i) {
            coords[0][i] = points.get(i).targetX;
            coords[1][i] = yValue;
        }
        return coords;
    }

    private static double[][] buildXYOffsetCoordinateMatrix(
        List<ProfileStation> points,
        double xOffset,
        double yValue,
        double maxX
    ) {
        final double[][] coords = new double[2][points.size()];
        for (int i = 0; i < points.size(); ++i) {
            coords[0][i] = Math.min(maxX - 1.0e-6, Math.max(1.0e-6, points.get(i).targetX + xOffset));
            coords[1][i] = yValue;
        }
        return coords;
    }

    private static double computeNormalizedLinearProfileL2(double[] values, List<ProfileStation> points, CaseConfig cfg) {
        double sumSq = 0.0;
        for (int i = 0; i < points.size(); ++i) {
            final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (points.get(i).targetX / cfg.lxM);
            final double err = values[i] - pLinear;
            sumSq += err * err;
        }
        return Math.sqrt(sumSq / Math.max(1, points.size())) / cfg.deltaPPa;
    }

    private static double computeNormalizedLinearProfileLinf(double[] values, List<ProfileStation> points, CaseConfig cfg) {
        double maxAbs = 0.0;
        for (int i = 0; i < points.size(); ++i) {
            final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (points.get(i).targetX / cfg.lxM);
            maxAbs = Math.max(maxAbs, Math.abs(values[i] - pLinear));
        }
        return maxAbs / cfg.deltaPPa;
    }

    private static double computeMaxAbsDifference(double[] lhs, double[] rhs) {
        double maxAbs = 0.0;
        final int count = Math.min(lhs.length, rhs.length);
        for (int i = 0; i < count; ++i) {
            maxAbs = Math.max(maxAbs, Math.abs(lhs[i] - rhs[i]));
        }
        return maxAbs;
    }

    private static int[] mapTimesToSolutionNumbers(double[] solutionTimes, double[] requestedTimes) {
        final int[] out = new int[requestedTimes.length];
        for (int i = 0; i < requestedTimes.length; ++i) {
            int bestIndex = 0;
            double bestDistance = Double.POSITIVE_INFINITY;
            for (int j = 0; j < solutionTimes.length; ++j) {
                final double distance = Math.abs(solutionTimes[j] - requestedTimes[i]);
                if (distance < bestDistance) {
                    bestDistance = distance;
                    bestIndex = j;
                }
            }
            out[i] = bestIndex + 1;
        }
        return out;
    }

    private static double computeNormalizedL2(double[] errors, double scale) {
        double sumSq = 0.0;
        for (double err : errors) {
            sumSq += err * err;
        }
        return Math.sqrt(sumSq / Math.max(1, errors.length)) / Math.max(scale, DP_NORM_MIN);
    }

    private static double computeNormalizedLinf(double[] errors, double scale) {
        double maxAbs = 0.0;
        for (double err : errors) {
            maxAbs = Math.max(maxAbs, Math.abs(err));
        }
        return maxAbs / Math.max(scale, DP_NORM_MIN);
    }

    private static double computeNormalizedLinearProfileL2WithXOffset(
        double[] values,
        List<ProfileStation> points,
        CaseConfig cfg,
        double xOffset
    ) {
        double sumSq = 0.0;
        for (int i = 0; i < points.size(); ++i) {
            final double shiftedX = Math.min(cfg.lxM - 1.0e-6, points.get(i).targetX + xOffset);
            final double pLinear = cfg.pLeftPa + (cfg.pRightPa - cfg.pLeftPa) * (shiftedX / cfg.lxM);
            final double err = values[i] - pLinear;
            sumSq += err * err;
        }
        return Math.sqrt(sumSq / Math.max(1, points.size())) / cfg.deltaPPa;
    }

    private static double average(double[] values) {
        double sum = 0.0;
        for (double value : values) {
            sum += value;
        }
        return sum / Math.max(1, values.length);
    }

    private static String sanitizeTag(String text) {
        return text.replace('-', 'm').replace('.', 'p');
    }

    private static void writeMonitorReferenceCsv(
        Model model,
        List<MonitorPoint> monitorPoints,
        List<MonitorSample> monitorSamples,
        CasePaths paths
    ) throws IOException {
        final double[] times = new double[monitorSamples.size()];
        for (int i = 0; i < monitorSamples.size(); ++i) {
            times[i] = monitorSamples.get(i).targetTimeS;
        }
        final double[][] pValues = evaluateMonitorPoints(model, monitorPoints, times, true, "interpMonitorP");
        final double[][] tValues = evaluateMonitorPoints(model, monitorPoints, times, false, "interpMonitorT");
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
        final double[][] mainP = evaluateProfileFamily(mainModel, profileStations, new double[]{finalRequest.timeS}, true, "interpProfileMainCheckP");
        final double[][] fineP = evaluateProfileFamily(fineModel, profileStations, new double[]{finalRequest.timeS}, true, "interpProfileFineCheckP");
        final double[][] mainT = evaluateProfileFamily(mainModel, profileStations, new double[]{finalRequest.timeS}, false, "interpProfileMainCheckT");
        final double[][] fineT = evaluateProfileFamily(fineModel, profileStations, new double[]{finalRequest.timeS}, false, "interpProfileFineCheckT");

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
            out.write("- Fracture representation: explicit_lower_dimensional_fracture\n");
            out.write("- Fracture aperture: `" + formatDouble(cfg.fractureWidthM) + " m`\n");
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

    private static double[][] evaluateProfileFamily(
        Model model,
        List<ProfileStation> points,
        double[] times,
        boolean pressure,
        String tagBase
    ) {
        return evaluateMixedLocations(
            model,
            points,
            extractProfileLocations(points),
            times,
            pressure ? "pMatrixPlot" : "TMatrixPlot",
            pressure ? "pf" : "Tf",
            pressure ? "Pa" : "K",
            tagBase
        );
    }

    private static double[][] evaluateMonitorPoints(
        Model model,
        List<MonitorPoint> points,
        double[] times,
        boolean pressure,
        String tagBase
    ) {
        return evaluateMixedLocations(
            model,
            points,
            extractMonitorLocations(points),
            times,
            pressure ? "pMatrixPlot" : "TMatrixPlot",
            pressure ? "pf" : "Tf",
            pressure ? "Pa" : "K",
            tagBase
        );
    }

    private static double[][] evaluateMixedLocations(
        Model model,
        List<? extends SpatialPoint> points,
        List<String> locations,
        double[] times,
        String matrixExpr,
        String fractureExpr,
        String unit,
        String tagBase
    ) {
        final double[][] merged = new double[times.length][points.size()];
        final List<Integer> matrixIdx = new ArrayList<Integer>();
        final List<Integer> fractureIdx = new ArrayList<Integer>();
        for (int i = 0; i < points.size(); ++i) {
            if ("fracture".equalsIgnoreCase(locations.get(i))) {
                fractureIdx.add(i);
            } else {
                matrixIdx.add(i);
            }
        }

        if (!matrixIdx.isEmpty()) {
            final double[][] matrixValues = evaluateAtPoints(
                model,
                tagBase + "_matrix",
                matrixExpr,
                unit,
                buildSubsetCoordinateMatrix(points, matrixIdx),
                times,
                2,
                null
            );
            scatterSubsetValues(merged, matrixValues, matrixIdx);
        }
        if (!fractureIdx.isEmpty()) {
            final double[][] fractureValues = evaluateAtPoints(
                model,
                tagBase + "_fracture",
                fractureExpr,
                unit,
                buildSubsetCoordinateMatrix(points, fractureIdx),
                times,
                1,
                "fractureBnd"
            );
            scatterSubsetValues(merged, fractureValues, fractureIdx);
        }
        return merged;
    }

    private static List<String> extractProfileLocations(List<ProfileStation> points) {
        final List<String> out = new ArrayList<String>(points.size());
        for (ProfileStation point : points) out.add(point.location);
        return out;
    }

    private static List<String> extractMonitorLocations(List<MonitorPoint> points) {
        final List<String> out = new ArrayList<String>(points.size());
        for (MonitorPoint point : points) out.add(point.location);
        return out;
    }

    private static double[][] buildSubsetCoordinateMatrix(List<? extends SpatialPoint> points, List<Integer> indices) {
        final double[][] coords = new double[2][indices.size()];
        for (int i = 0; i < indices.size(); ++i) {
            final SpatialPoint point = points.get(indices.get(i).intValue());
            coords[0][i] = point.targetX;
            coords[1][i] = point.targetY;
        }
        return coords;
    }

    private static void scatterSubsetValues(double[][] merged, double[][] subsetValues, List<Integer> indices) {
        for (int t = 0; t < subsetValues.length; ++t) {
            for (int i = 0; i < indices.size(); ++i) {
                merged[t][indices.get(i).intValue()] = subsetValues[t][i];
            }
        }
    }

    private static double[][] evaluateAtPoints(
        Model model,
        String numericalTag,
        String expression,
        String unit,
        double[][] coords,
        double[] times,
        int entityDim,
        String selectionName
    ) {
        return evaluateAtPoints(model, numericalTag, expression, unit, coords, times, entityDim, selectionName, "dset1");
    }

    private static double[][] evaluateAtPoints(
        Model model,
        String numericalTag,
        String expression,
        String unit,
        double[][] coords,
        double[] times,
        int entityDim,
        String selectionName,
        String datasetTag
    ) {
        model.result().numerical().create(numericalTag, "Interp");
        final NumericalFeature interp = model.result().numerical(numericalTag);
        interp.set("expr", new String[]{expression});
        interp.set("edim", entityDim);
        interp.set("coorderr", "on");
        interp.set("unit", new String[]{unit});
        interp.set("data", datasetTag);
        interp.set("t", times);
        if (selectionName != null && !selectionName.trim().isEmpty()) {
            try {
                interp.selection().named(selectionName);
            } catch (RuntimeException ex) {
                throw new IllegalStateException("Failed to bind interpolation selection '" + selectionName + "'.", ex);
            }
        } else {
            try {
                interp.selection().all();
            } catch (RuntimeException ignored) {
            }
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

    private static double[][] evaluateAtPointsBySolutionNumber(
        Model model,
        String numericalTag,
        String expression,
        String unit,
        double[][] coords,
        int[] solutionNumbers,
        int entityDim,
        String selectionName
    ) {
        return evaluateAtPointsBySolutionNumber(
            model, numericalTag, expression, unit, coords, solutionNumbers, entityDim, selectionName, "dset1"
        );
    }

    private static double[][] evaluateAtPointsBySolutionNumber(
        Model model,
        String numericalTag,
        String expression,
        String unit,
        double[][] coords,
        int[] solutionNumbers,
        int entityDim,
        String selectionName,
        String datasetTag
    ) {
        model.result().numerical().create(numericalTag, "Interp");
        final NumericalFeature interp = model.result().numerical(numericalTag);
        interp.set("expr", new String[]{expression});
        interp.set("edim", entityDim);
        interp.set("coorderr", "on");
        interp.set("unit", new String[]{unit});
        interp.set("data", datasetTag);

        final String[] propertyNames = interp.properties();
        if (!arrayContains(propertyNames, "solnum")) {
            throw new IllegalStateException(
                "Interpolation feature '" + numericalTag + "' does not expose solnum property. properties=" +
                Arrays.toString(propertyNames)
            );
        }
        interp.set("solnum", solutionNumbers);

        if (selectionName != null && !selectionName.trim().isEmpty()) {
            try {
                interp.selection().named(selectionName);
            } catch (RuntimeException ex) {
                throw new IllegalStateException("Failed to bind interpolation selection '" + selectionName + "'.", ex);
            }
        } else {
            try {
                interp.selection().all();
            } catch (RuntimeException ignored) {
            }
        }

        final List<InterpolationCandidate> candidates = buildInterpolationCandidates(coords);
        final List<String> diagnostics = new ArrayList<String>();
        for (InterpolationCandidate candidate : candidates) {
            interp.set("coord", candidate.coords);
            interp.run();

            final double[][] raw = tryExtractInterpolationData(interp);
            if (hasData(raw)) {
                return orientAsTimePointMatrix(
                    raw, candidate.pointCount, solutionNumbers.length, numericalTag + "[" + candidate.label + "]"
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
            "No interpolation data returned for '" + numericalTag + "' using solnum." +
            " expectedPointCount=" + coords[0].length +
            ", expectedTimeCount=" + solutionNumbers.length +
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

    private static void ensureExplicitSolutionDataset(Model model, String datasetTag) {
        try {
            model.result().dataset().remove(datasetTag);
        } catch (RuntimeException ignored) {
        }

        model.result().dataset().create(datasetTag, "Solution");
        try {
            model.result().dataset(datasetTag).set("solution", "sol1");
        } catch (RuntimeException ex) {
            throw new IllegalStateException(
                "Failed to bind explicit solution dataset '" + datasetTag + "' to sol1.",
                ex
            );
        }

        try {
            model.result().dataset(datasetTag).set("studystep", "time");
        } catch (RuntimeException ignored) {
        }
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

    private static double span(double[] values) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (double value : values) {
            min = Math.min(min, value);
            max = Math.max(max, value);
        }
        return max - min;
    }

    private static double maxAbs(double[] values) {
        double out = 0.0;
        for (double value : values) {
            out = Math.max(out, Math.abs(value));
        }
        return out;
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

    private static double[] buildUniformStudyTimes(double endTimeS, double dtS) {
        final double safeDt = Math.max(1.0, dtS);
        final int intervals = Math.max(1, (int) Math.ceil(endTimeS / safeDt));
        final double[] out = new double[intervals + 1];
        for (int i = 0; i < intervals; ++i) {
            out[i] = Math.min(endTimeS, i * safeDt);
        }
        out[intervals] = endTimeS;
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
    double dtInitS;
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

final class LegacyPressureOnlyTables {
    Path lambdaTable;
    Path storageTable;
}

final class FullPressureCoefficientTables {
    Path matrixMobTable;
    Path matrixStorageTable;
    Path fractureMobTable;
    Path fractureStorageTable;
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
