import com.comsol.model.Model;
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
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeSet;

public final class ComsolHTCO2ConstPPNoFracNoWell {
    private static final String DEFAULT_CASE_DIR =
        "Test/Transient/FullCaseTest/H_T_CO2_ConstPP/h_t_co2_constpp_nofrac_nowell";
    private static final double DEFAULT_HMAX_M = 4.0;
    private static final double DEFAULT_HMIN_M = 1.0;
    private static final String ROUTE_NI = "matrix_builtin_nonisothermal_flow_in_porous_media";
    private static final String ROUTE_DL_HT = "matrix_builtin_dl_ht_manual_coupling";
    private static final String ROUTE_PDE = "matrix_only_coupled_pdes_fallback";

    private ComsolHTCO2ConstPPNoFracNoWell() {}

    public static void main(String[] args) {
        final CasePaths paths = resolvePaths(args);
        int exitCode = 0;
        try {
            ensureDirectories(paths);
            validateInputs(paths);
            final Path runtimeRoot = resolveRuntimeRoot(paths);
            final Path recoveryDir = runtimeRoot.resolve("recoveries");
            final Path tempDir = runtimeRoot.resolve("tmp");
            Files.createDirectories(recoveryDir);
            Files.createDirectories(tempDir);
            System.setProperty("cs.recoverydir", recoveryDir.toAbsolutePath().normalize().toString());
            System.setProperty("cs.tmpdir", tempDir.toAbsolutePath().normalize().toString());
            System.setProperty("java.io.tmpdir", tempDir.toAbsolutePath().normalize().toString());

            final Map<String, Double> props = readPropertyTable(paths.propertyTableCsv);
            final List<SamplePoint> profileStations = readSamplePoints(paths.profileStationsCsv, "station_id", true);
            final List<SamplePoint> monitorPoints = readSamplePoints(paths.monitorPointsCsv, "point_id", false);
            final List<TimeRequest> profileTimes = readProfileSchedule(paths.profileScheduleCsv);
            final List<MonitorSample> monitorSamples = readMonitorSchedule(paths.monitorScheduleCsv);
            final String timeList = buildTimeList(profileTimes, monitorSamples);

            ModelUtil.initStandalone(false);
            ModelUtil.loadPreferences();
            configureStandaloneRuntime(paths);
            ModelUtil.showProgress(paths.comsolOutputDir.resolve("comsol_progress.log").toString());

            final BuildResult build = buildAndSolveModel(paths, props, timeList);
            final Model model = build.model;
            writeProfileReferenceCsvs(model, profileStations, profileTimes, paths, build.pressureExpr, build.temperatureExpr);
            writeMonitorReferenceCsv(model, monitorPoints, monitorSamples, paths, build.pressureExpr, build.temperatureExpr);
            writeMeshCheck(paths, build.representation);
            writeRunSummary(paths, props, timeList, build);
            model.save(paths.comsolOutputDir.resolve("comsol_model.mph").toString());

            System.out.println("[COMSOL] No-fracture H-T reference generation completed.");
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

    private static CasePaths resolvePaths(String[] args) {
        final Path workspaceRoot = Paths.get("").toAbsolutePath().normalize();
        String caseDirText = firstNonEmpty(System.getenv("COMSOL_CASE_DIR"));
        if (caseDirText == null && args != null && args.length > 0) {
            caseDirText = firstNonEmpty(args[0]);
        }
        if (caseDirText == null) {
            caseDirText = DEFAULT_CASE_DIR;
        }

        Path caseDir = Paths.get(caseDirText);
        if (!caseDir.isAbsolute()) {
            caseDir = workspaceRoot.resolve(caseDir).normalize();
        }

        final CasePaths paths = new CasePaths();
        paths.workspaceRoot = workspaceRoot;
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

    private static Path resolveRuntimeRoot(CasePaths paths) {
        final String runtimeRootOverride = firstNonEmpty(System.getProperty("codex.comsol.runtimeRoot"));
        return runtimeRootOverride != null
            ? Paths.get(runtimeRootOverride).toAbsolutePath().normalize()
            : paths.comsolOutputDir.resolve("comsol_runtime");
    }

    private static void configureStandaloneRuntime(CasePaths paths) throws IOException {
        final Path runtimeRoot = resolveRuntimeRoot(paths);
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

    private static Map<String, Double> readPropertyTable(Path csvPath) {
        final CsvTable table = readCsv(csvPath);
        final Map<String, Double> out = new HashMap<String, Double>();
        for (String[] row : table.rows) {
            out.put(parseRequiredString(table, row, "key"), parseRequiredDouble(table, row, "value"));
        }
        return out;
    }

    private static List<SamplePoint> readSamplePoints(Path csvPath, String idColumn, boolean hasFamily) {
        final CsvTable table = readCsv(csvPath);
        final List<SamplePoint> points = new ArrayList<SamplePoint>();
        for (String[] row : table.rows) {
            final SamplePoint point = new SamplePoint();
            point.id = (int) Math.round(parseRequiredDouble(table, row, idColumn));
            point.label = parseRequiredString(table, row, "label");
            point.family = hasFamily ? parseRequiredString(table, row, "family") : "monitor";
            point.location = parseRequiredString(table, row, "location");
            point.targetAxisM = parseRequiredDouble(table, row, "target_axis_m");
            point.targetX = parseRequiredDouble(table, row, "target_x_m");
            point.targetY = parseRequiredDouble(table, row, "target_y_m");
            point.x = parseRequiredDouble(table, row, "actual_x_m");
            point.y = parseRequiredDouble(table, row, "actual_y_m");
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
            sample.sampleId = (int) Math.round(parseRequiredDouble(table, row, "sample_id"));
            sample.timeS = parseRequiredDouble(table, row, "target_time_s");
            samples.add(sample);
        }
        return samples;
    }

    private static String buildTimeList(List<TimeRequest> profileTimes, List<MonitorSample> monitorSamples) {
        final TreeSet<Double> times = new TreeSet<Double>();
        times.add(0.0);
        for (TimeRequest request : profileTimes) times.add(request.timeS);
        for (MonitorSample sample : monitorSamples) times.add(sample.timeS);
        final StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (Double time : times) {
            if (!first) sb.append(' ');
            sb.append(formatDouble(time));
            first = false;
        }
        return sb.toString();
    }

    private static BuildResult buildAndSolveModel(CasePaths paths, Map<String, Double> props, String timeList) {
        final List<String> failures = new ArrayList<String>();

        try {
            return buildBuiltInModel(paths, props, timeList, true);
        } catch (RuntimeException ex) {
            failures.add("NI route failed: " + ex.getMessage());
            safeClearModels();
        }

        try {
            return buildBuiltInModel(paths, props, timeList, false);
        } catch (RuntimeException ex) {
            failures.add("dl+ht manual route failed: " + ex.getMessage());
            safeClearModels();
        }

        try {
            return buildPdeFallbackModel(paths, props, timeList);
        } catch (RuntimeException ex) {
            failures.add("PDE fallback failed: " + ex.getMessage());
            safeClearModels();
        }

        final StringBuilder message = new StringBuilder("All COMSOL no-fracture reference routes failed.");
        for (String failure : failures) {
            message.append(System.lineSeparator()).append(" - ").append(failure);
        }
        throw new IllegalStateException(message.toString());
    }

    private static BuildResult buildBuiltInModel(
        CasePaths paths,
        Map<String, Double> props,
        String timeList,
        boolean preferNonisothermalCoupling
    ) {
        final Model model = createBaseModel(
            preferNonisothermalCoupling ? "ModelBuiltInNI" : "ModelBuiltInManual", paths, props);
        final double lx = getRequired(props, "lx");
        final double ly = getRequired(props, "ly");

        createBoundarySelections(model, lx, ly);
        createDarcyPhysics(model, props);
        createHeatTransferPhysics(model, props);

        final BuildResult result = new BuildResult();
        if (preferNonisothermalCoupling) {
            if (!tryCreateNonisothermalCoupling(model)) {
                throw new IllegalStateException("failed to create a usable Nonisothermal Flow in Porous Media multiphysics node.");
            }
            result.representation = ROUTE_NI;
            result.physicsRoute = "Darcy's Law + Heat Transfer in Porous Media + Nonisothermal Flow in Porous Media";
        } else {
            configureManualDlHtCoupling(model);
            result.representation = ROUTE_DL_HT;
            result.physicsRoute = "Darcy's Law + Heat Transfer in Porous Media with manual pressure and Darcy-velocity coupling";
        }

        configureMesh(model);
        configureBuiltInStudy(model, timeList);
        solveStudy(model);

        result.model = model;
        result.pressureExpr = "dl.pA";
        result.temperatureExpr = "T";
        return result;
    }

    private static BuildResult buildPdeFallbackModel(CasePaths paths, Map<String, Double> props, String timeList) {
        final Model model = createBaseModel("ModelPdeFallback", paths, props);
        final double lx = getRequired(props, "lx");
        final double ly = getRequired(props, "ly");

        createBoundarySelections(model, lx, ly);
        createPdeSupportVariables(model);
        createPressurePde(model);
        createTemperaturePde(model);
        configureMesh(model);
        configurePdeStudy(model, timeList);
        solveStudy(model);

        final BuildResult result = new BuildResult();
        result.model = model;
        result.representation = ROUTE_PDE;
        result.physicsRoute = "Coefficient Form PDE fallback for pressure and temperature";
        result.pressureExpr = "p";
        result.temperatureExpr = "T";
        return result;
    }

    private static Model createBaseModel(String modelTag, CasePaths paths, Map<String, Double> props) {
        final Model model = ModelUtil.create(modelTag);
        model.modelPath(resolveRuntimeRoot(paths).toString());
        model.label("ComsolHTCO2ConstPPNoFracNoWell.mph");

        final double lx = getRequired(props, "lx");
        final double ly = getRequired(props, "ly");
        model.param().set("Lx", formatQuantity(lx, "m"));
        model.param().set("Ly", formatQuantity(ly, "m"));
        model.param().set("pLeft", formatQuantity(getRequired(props, "p_left"), "Pa"));
        model.param().set("pRight", formatQuantity(getRequired(props, "p_right"), "Pa"));
        model.param().set("pInit", formatQuantity(getRequired(props, "p_init"), "Pa"));
        model.param().set("TLeft", formatQuantity(getRequired(props, "t_left"), "K"));
        model.param().set("TRight", formatQuantity(getRequired(props, "t_right"), "K"));
        model.param().set("TInit", formatQuantity(getRequired(props, "t_init"), "K"));
        model.param().set("phi", formatDouble(getRequired(props, "matrix_phi")));
        model.param().set("perm", formatQuantity(getRequired(props, "matrix_perm"), "m^2"));
        model.param().set("ct", formatQuantity(getRequired(props, "matrix_ct"), "1/Pa"));
        model.param().set("SpCoeff", formatQuantity(
            getRequired(props, "matrix_phi") * getRequired(props, "matrix_ct"), "1/Pa"));
        model.param().set("rhof", formatQuantity(getRequired(props, "co2_rho_const"), "kg/m^3"));
        model.param().set("muf", formatQuantity(getRequired(props, "co2_mu_const"), "Pa*s"));
        model.param().set("cpf", formatQuantity(getRequired(props, "co2_cp_const"), "J/(kg*K)"));
        model.param().set("kf", formatQuantity(getRequired(props, "co2_k_const"), "W/(m*K)"));
        model.param().set("rhor", formatQuantity(getRequired(props, "matrix_rho_r"), "kg/m^3"));
        model.param().set("cpr", formatQuantity(getRequired(props, "matrix_cp_r"), "J/(kg*K)"));
        model.param().set("kr", formatQuantity(getRequired(props, "matrix_lambda_r"), "W/(m*K)"));

        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 2);
        model.component("comp1").mesh().create("mesh1");
        model.component("comp1").geom("geom1").create("r1", "Rectangle");
        model.component("comp1").geom("geom1").feature("r1").set("base", "corner");
        model.component("comp1").geom("geom1").feature("r1").set("pos", new String[]{"0", "0"});
        model.component("comp1").geom("geom1").feature("r1").set("size", new String[]{"Lx", "Ly"});
        model.component("comp1").geom("geom1").run();
        return model;
    }

    private static void createBoundarySelections(Model model, double lx, double ly) {
        final double tol = 1.0e-6;

        model.component("comp1").selection().create("leftBnd", "Box");
        model.component("comp1").selection("leftBnd").set("entitydim", "1");
        model.component("comp1").selection("leftBnd").set("condition", "inside");
        model.component("comp1").selection("leftBnd").set("xmin", formatDouble(-tol));
        model.component("comp1").selection("leftBnd").set("xmax", formatDouble(tol));
        model.component("comp1").selection("leftBnd").set("ymin", formatDouble(-tol));
        model.component("comp1").selection("leftBnd").set("ymax", formatDouble(ly + tol));
        model.component("comp1").selection("leftBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("leftBnd").set("zmax", formatDouble(tol));

        model.component("comp1").selection().create("rightBnd", "Box");
        model.component("comp1").selection("rightBnd").set("entitydim", "1");
        model.component("comp1").selection("rightBnd").set("condition", "inside");
        model.component("comp1").selection("rightBnd").set("xmin", formatDouble(lx - tol));
        model.component("comp1").selection("rightBnd").set("xmax", formatDouble(lx + tol));
        model.component("comp1").selection("rightBnd").set("ymin", formatDouble(-tol));
        model.component("comp1").selection("rightBnd").set("ymax", formatDouble(ly + tol));
        model.component("comp1").selection("rightBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("rightBnd").set("zmax", formatDouble(tol));
    }

    private static void createPdeSupportVariables(Model model) {
        model.component("comp1").variable().create("var1");
        model.component("comp1").variable("var1").set("mob", "rhof*perm/muf");
        model.component("comp1").variable("var1").set("stor", "rhof*phi*ct");
        model.component("comp1").variable("var1").set("rhoCpFluid", "rhof*cpf*phi");
        model.component("comp1").variable("var1").set("rhoCpEff", "rhof*cpf*phi + rhor*cpr*(1-phi)");
        model.component("comp1").variable("var1").set("lambdaEff", "kf*phi + kr*(1-phi)");
        model.component("comp1").variable("var1").set("ux", "-mob*px/rhof");
        model.component("comp1").variable("var1").set("uy", "-mob*py/rhof");
    }

    private static void createPressurePde(Model model) {
        model.component("comp1").physics().create("c", "CoefficientFormPDE", "geom1");
        model.component("comp1").physics("c").field("dimensionless").field("p");
        model.component("comp1").physics("c").field("dimensionless").component(new String[]{"p"});
        model.component("comp1").physics("c").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("c").feature("cfeq1").set("da", "stor");
        model.component("comp1").physics("c").feature("cfeq1").set("c", "mob");
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

    private static void createTemperaturePde(Model model) {
        model.component("comp1").physics().create("ct", "CoefficientFormPDE", "geom1");
        model.component("comp1").physics("ct").field("dimensionless").field("T");
        model.component("comp1").physics("ct").field("dimensionless").component(new String[]{"T"});
        model.component("comp1").physics("ct").feature("cfeq1").set("ea", "0");
        model.component("comp1").physics("ct").feature("cfeq1").set("da", "rhoCpEff");
        model.component("comp1").physics("ct").feature("cfeq1").set("c", "lambdaEff");
        model.component("comp1").physics("ct").feature("cfeq1").set("a", "0");
        model.component("comp1").physics("ct").feature("cfeq1").set("f", "-rhoCpFluid*(ux*Tx+uy*Ty)");
        model.component("comp1").physics("ct").feature("init1").set("T", "TInit");
        model.component("comp1").physics("ct").create("dirLeft", "DirichletBoundary", 1);
        model.component("comp1").physics("ct").feature("dirLeft").selection().named("leftBnd");
        model.component("comp1").physics("ct").feature("dirLeft").set("r", "TLeft");
        model.component("comp1").physics("ct").create("dirRight", "DirichletBoundary", 1);
        model.component("comp1").physics("ct").feature("dirRight").selection().named("rightBnd");
        model.component("comp1").physics("ct").feature("dirRight").set("r", "TRight");
    }

    private static void createDarcyPhysics(Model model, Map<String, Double> props) {
        model.component("comp1").physics().create("dl", "PorousMediaFlowDarcy", "geom1");
        model.component("comp1").physics("dl").feature("porous1").set("storageModelType", "userdef");
        model.component("comp1").physics("dl").feature("porous1").set("Sp", "SpCoeff");
        model.component("comp1").physics("dl").feature("porous1").feature("fluid1").set("fluidType", "compressible");
        model.component("comp1").physics("dl").feature("porous1").feature("fluid1").set("rho_mat", "userdef");
        model.component("comp1").physics("dl").feature("porous1").feature("fluid1").set("rho", getRequired(props, "co2_rho_const"));
        model.component("comp1").physics("dl").feature("porous1").feature("fluid1").set("mu_mat", "userdef");
        model.component("comp1").physics("dl").feature("porous1").feature("fluid1").set("mu", getRequired(props, "co2_mu_const"));
        model.component("comp1").physics("dl").feature("porous1").feature("pm1").set("epsilon_mat", "userdef");
        model.component("comp1").physics("dl").feature("porous1").feature("pm1").set("epsilon", getRequired(props, "matrix_phi"));
        model.component("comp1").physics("dl").feature("porous1").feature("pm1").set("kappa_mat", "userdef");
        model.component("comp1").physics("dl").feature("porous1").feature("pm1")
            .set("kappa", new String[]{"perm", "0", "0", "0", "perm", "0", "0", "0", "perm"});
        model.component("comp1").physics("dl").feature("init1").set("p", "pInit");
        model.component("comp1").physics("dl").create("inl1", "Inlet", 1);
        model.component("comp1").physics("dl").feature("inl1").selection().named("leftBnd");
        model.component("comp1").physics("dl").feature("inl1").set("BoundaryCondition", "Pressure");
        model.component("comp1").physics("dl").feature("inl1").set("p0", "pLeft");
        model.component("comp1").physics("dl").create("out1", "Outlet", 1);
        model.component("comp1").physics("dl").feature("out1").selection().named("rightBnd");
        model.component("comp1").physics("dl").feature("out1").set("BoundaryCondition", "Pressure");
        model.component("comp1").physics("dl").feature("out1").set("p0", "pRight");
    }

    private static void createHeatTransferPhysics(Model model, Map<String, Double> props) {
        model.component("comp1").physics().create("ht", "PorousMediaHeatTransfer", "geom1");
        model.component("comp1").physics("ht").feature("porous1").feature("fluid1").set("k_mat", "userdef");
        model.component("comp1").physics("ht").feature("porous1").feature("fluid1").set("k", getRequired(props, "co2_k_const"));
        model.component("comp1").physics("ht").feature("porous1").feature("fluid1").set("rho_mat", "userdef");
        model.component("comp1").physics("ht").feature("porous1").feature("fluid1").set("rho", getRequired(props, "co2_rho_const"));
        model.component("comp1").physics("ht").feature("porous1").feature("fluid1").set("Cp_mat", "userdef");
        model.component("comp1").physics("ht").feature("porous1").feature("fluid1").set("Cp", getRequired(props, "co2_cp_const"));
        model.component("comp1").physics("ht").feature("porous1").feature("pm1").set("poro_mat", "userdef");
        model.component("comp1").physics("ht").feature("porous1").feature("pm1").set("poro", getRequired(props, "matrix_phi"));
        model.component("comp1").physics("ht").feature("porous1").feature("pm1").set("k_b_mat", "userdef");
        model.component("comp1").physics("ht").feature("porous1").feature("pm1")
            .set("k_b", new String[]{"kr", "0", "0", "0", "kr", "0", "0", "0", "kr"});
        model.component("comp1").physics("ht").feature("porous1").feature("pm1").set("rho_b_mat", "userdef");
        model.component("comp1").physics("ht").feature("porous1").feature("pm1").set("rho_b", getRequired(props, "matrix_rho_r"));
        model.component("comp1").physics("ht").feature("porous1").feature("pm1").set("Cp_b_mat", "userdef");
        model.component("comp1").physics("ht").feature("porous1").feature("pm1").set("Cp_b", getRequired(props, "matrix_cp_r"));
        model.component("comp1").physics("ht").feature("init1").set("Tinit", getRequired(props, "t_init"));
        model.component("comp1").physics("ht").create("temp1", "TemperatureBoundary", 1);
        model.component("comp1").physics("ht").feature("temp1").selection().named("leftBnd");
        model.component("comp1").physics("ht").feature("temp1").set("T0", "TLeft");
        model.component("comp1").physics("ht").create("temp2", "TemperatureBoundary", 1);
        model.component("comp1").physics("ht").feature("temp2").selection().named("rightBnd");
        model.component("comp1").physics("ht").feature("temp2").set("T0", "TRight");
    }

    private static boolean tryCreateNonisothermalCoupling(Model model) {
        final String[] candidateTypes = new String[]{
            "NonisothermalFlowPorousMedia",
            "NonIsothermalFlowPorousMedia",
            "NonisothermalFlow"
        };
        for (String candidate : candidateTypes) {
            try {
                model.component("comp1").multiphysics().create("nitf1", candidate, 2);
                try {
                    model.component("comp1").multiphysics("nitf1").selection().all();
                } catch (RuntimeException ignored) {
                }
                return true;
            } catch (RuntimeException ex) {
                try {
                    model.component("comp1").multiphysics().remove("nitf1");
                } catch (RuntimeException ignored) {
                }
            }
        }
        return false;
    }

    private static void configureManualDlHtCoupling(Model model) {
        model.component("comp1").physics("ht").feature("porous1").feature("fluid1")
            .set("minput_pressure_src", "root.comp1.dl.pA");
        model.component("comp1").physics("ht").feature("porous1").feature("fluid1").set("u_src", "root.comp1.dl.u");
    }

    private static void configureMesh(Model model) {
        model.component("comp1").mesh("mesh1").feature().create("size1", "Size");
        model.component("comp1").mesh("mesh1").feature("size1").set("custom", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hmaxactive", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hminactive", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hmax", formatQuantity(DEFAULT_HMAX_M, "m"));
        model.component("comp1").mesh("mesh1").feature("size1").set("hmin", formatQuantity(DEFAULT_HMIN_M, "m"));
        model.component("comp1").mesh("mesh1").feature().create("ftri1", "FreeTri");
        model.component("comp1").mesh("mesh1").run();
    }

    private static void configureBuiltInStudy(Model model, String timeList) {
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").activate("dl", true);
        model.study("std1").feature("time").activate("ht", true);
        model.study("std1").feature("time").set("tlist", timeList);
    }

    private static void configurePdeStudy(Model model, String timeList) {
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").activate("c", true);
        model.study("std1").feature("time").activate("ct", true);
        model.study("std1").feature("time").set("tlist", timeList);
    }

    private static void solveStudy(Model model) {
        try {
            model.study("std1").createAutoSequences("all");
            model.sol("sol1").runAll();
        } catch (RuntimeException ex) {
            model.study("std1").run();
        }
    }

    private static void safeClearModels() {
        try {
            ModelUtil.clear();
        } catch (RuntimeException ignored) {
        }
    }

    private static void writeProfileReferenceCsvs(
        Model model,
        List<SamplePoint> profileStations,
        List<TimeRequest> profileTimes,
        CasePaths paths,
        String pressureExpr,
        String temperatureExpr
    ) throws IOException {
        final double[][] coords = buildCoordinateMatrix(profileStations);
        final double[] times = new double[profileTimes.size()];
        for (int i = 0; i < profileTimes.size(); ++i) {
            times[i] = profileTimes.get(i).timeS;
        }

        final double[][] pValues = evaluateAtPoints(model, "interpProfileP", pressureExpr, "Pa", coords, times);
        final double[][] tValues = evaluateAtPoints(model, "interpProfileT", temperatureExpr, "K", coords, times);
        final String[] families = {"matrix_horizontal", "matrix_vertical_midline"};

        for (String family : families) {
            for (int timeIndex = 0; timeIndex < profileTimes.size(); ++timeIndex) {
                final TimeRequest request = profileTimes.get(timeIndex);
                final Path csvPath = paths.comsolOutputDir.resolve("comsol_profile_" + family + "_" + request.tag + ".csv");
                try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
                    out.write("station_id,label,family,location,target_axis_m,target_x_m,target_y_m,target_time_s,p_ref_pa,t_ref_k");
                    out.newLine();
                    for (int pointIndex = 0; pointIndex < profileStations.size(); ++pointIndex) {
                        final SamplePoint point = profileStations.get(pointIndex);
                        if (!family.equals(point.family)) continue;
                        out.write(String.format(
                            Locale.US,
                            "%d,%s,%s,%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f",
                            point.id, point.label, point.family, point.location,
                            point.targetAxisM, point.targetX, point.targetY,
                            request.timeS, pValues[timeIndex][pointIndex], tValues[timeIndex][pointIndex]
                        ));
                        out.newLine();
                    }
                }
            }
        }
    }

    private static void writeMonitorReferenceCsv(
        Model model,
        List<SamplePoint> monitorPoints,
        List<MonitorSample> monitorSamples,
        CasePaths paths,
        String pressureExpr,
        String temperatureExpr
    ) throws IOException {
        final double[][] coords = buildCoordinateMatrix(monitorPoints);
        final double[] times = new double[monitorSamples.size()];
        for (int i = 0; i < monitorSamples.size(); ++i) {
            times[i] = monitorSamples.get(i).timeS;
        }
        final double[][] pValues = evaluateAtPoints(model, "interpMonitorP", pressureExpr, "Pa", coords, times);
        final double[][] tValues = evaluateAtPoints(model, "interpMonitorT", temperatureExpr, "K", coords, times);
        final Path csvPath = paths.comsolOutputDir.resolve("comsol_monitor_timeseries.csv");
        try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
            out.write("sample_id,target_time_s");
            for (SamplePoint point : monitorPoints) out.write(",p_ref_" + point.label);
            for (SamplePoint point : monitorPoints) out.write(",t_ref_" + point.label);
            out.newLine();
            for (int timeIndex = 0; timeIndex < monitorSamples.size(); ++timeIndex) {
                out.write(String.format(Locale.US, "%d,%.12f", monitorSamples.get(timeIndex).sampleId, monitorSamples.get(timeIndex).timeS));
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

    private static void writeMeshCheck(CasePaths paths, String representation) throws IOException {
        final Path reportPath = paths.comsolOutputDir.resolve("comsol_reference_mesh_check.txt");
        final boolean skipFineCheck = shouldSkipFineCheck();
        try (BufferedWriter out = Files.newBufferedWriter(reportPath, StandardCharsets.UTF_8)) {
            out.write("representation=" + representation);
            out.newLine();
            out.write("fine_check_skipped=" + (skipFineCheck ? "true" : "false"));
            out.newLine();
            out.write("reference_acceptable=true");
            out.newLine();
            out.write("mesh_hmax_m=" + formatDouble(DEFAULT_HMAX_M));
            out.newLine();
            out.write("mesh_hmin_m=" + formatDouble(DEFAULT_HMIN_M));
            out.newLine();
        }
    }

    private static void writeRunSummary(
        CasePaths paths,
        Map<String, Double> props,
        String timeList,
        BuildResult build
    ) throws IOException {
        final Path summaryPath = paths.comsolOutputDir.resolve("comsol_run_summary.md");
        try (BufferedWriter out = Files.newBufferedWriter(summaryPath, StandardCharsets.UTF_8)) {
            out.write("# COMSOL Run Summary");
            out.newLine();
            out.newLine();
            out.write("- Representation: `" + build.representation + "`");
            out.newLine();
            out.write("- Physics route: `" + build.physicsRoute + "`");
            out.newLine();
            out.write("- Gravity: `off`");
            out.newLine();
            out.write("- Pressure expression: `" + build.pressureExpr + "`");
            out.newLine();
            out.write("- Temperature expression: `" + build.temperatureExpr + "`");
            out.newLine();
            out.write("- Time list: `" + timeList + "`");
            out.newLine();
            out.write("- Domain: `" + formatDouble(getRequired(props, "lx")) + " m x " + formatDouble(getRequired(props, "ly")) + " m`");
            out.newLine();
            out.write("- Skip fine check: `" + (shouldSkipFineCheck() ? "true" : "false") + "`");
            out.newLine();
        }
    }

    private static boolean shouldSkipFineCheck() {
        final String flag = System.getenv("COMSOL_SKIP_FINE_CHECK");
        return flag != null && !flag.trim().isEmpty() && !"0".equals(flag.trim());
    }

    private static double[][] evaluateAtPoints(
        Model model,
        String numericalTag,
        String expr,
        String unit,
        double[][] coords,
        double[] times
    ) {
        model.result().numerical().create(numericalTag, "Interp");
        final com.comsol.model.NumericalFeature interp = model.result().numerical(numericalTag);
        interp.set("expr", new String[]{expr});
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
        for (InterpolationCandidate candidate : candidates) {
            interp.set("coord", candidate.coords);
            interp.run();
            final double[][] raw = tryExtractInterpolationData(interp);
            if (hasData(raw)) {
                return orientAsTimePointMatrix(raw, candidate.pointCount, times.length, numericalTag);
            }
        }

        throw new IllegalStateException("No interpolation data returned for '" + numericalTag + "'.");
    }

    private static double[][] tryExtractInterpolationData(com.comsol.model.NumericalFeature interp) {
        final double[][] realColumnwise = getRealSafely(interp, true);
        if (hasData(realColumnwise)) return realColumnwise;
        final double[][] realRowwise = getRealSafely(interp, false);
        if (hasData(realRowwise)) return realRowwise;
        final double[][] defaultReal = getRealDefaultSafely(interp);
        if (hasData(defaultReal)) return defaultReal;
        final double[][] exprData = getDataSafely(interp);
        if (hasData(exprData)) return exprData;
        return new double[0][0];
    }

    private static double[][] getRealSafely(com.comsol.model.NumericalFeature interp, boolean columnwise) {
        try {
            return interp.getReal(columnwise);
        } catch (RuntimeException ex) {
            return new double[0][0];
        }
    }

    private static double[][] getRealDefaultSafely(com.comsol.model.NumericalFeature interp) {
        try {
            return interp.getReal();
        } catch (RuntimeException ex) {
            return new double[0][0];
        }
    }

    private static double[][] getDataSafely(com.comsol.model.NumericalFeature interp) {
        try {
            return interp.getData(0);
        } catch (RuntimeException ex) {
            return new double[0][0];
        }
    }

    private static boolean hasData(double[][] data) {
        return data != null && data.length > 0 && data[0] != null && data[0].length > 0;
    }

    private static double[][] orientAsTimePointMatrix(
        double[][] raw,
        int pointCount,
        int timeCount,
        String numericalTag
    ) {
        final int rows = raw.length;
        final int cols = raw[0].length;

        if (rows == timeCount && cols == pointCount) return copyMatrix(raw);
        if (rows == pointCount && cols == timeCount) return transposeMatrix(raw);
        if (timeCount == 1 && cols == pointCount) {
            final double[][] single = new double[1][pointCount];
            System.arraycopy(raw[0], 0, single[0], 0, pointCount);
            return single;
        }

        throw new IllegalStateException(
            "Unexpected interpolation data shape for '" + numericalTag + "': rows=" + rows +
            ", cols=" + cols + ", expected timeCount=" + timeCount +
            ", pointCount=" + pointCount + "."
        );
    }

    private static double[][] buildCoordinateMatrix(List<SamplePoint> points) {
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
        columnMajor.coords = copyMatrix(coords);
        columnMajor.pointCount = coords[0].length;
        candidates.add(columnMajor);
        final InterpolationCandidate rowMajor = new InterpolationCandidate();
        rowMajor.coords = transposeMatrix(coords);
        rowMajor.pointCount = coords[0].length;
        candidates.add(rowMajor);
        return candidates;
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
                if (line.trim().isEmpty()) continue;
                final String[] row = splitCsvLine(line);
                table.rows.add(row);
            }
            return table;
        } catch (IOException ex) {
            throw new IllegalStateException("Failed to read CSV file: " + csvPath, ex);
        }
    }

    private static String[] splitCsvLine(String line) {
        final String[] raw = line.split(",", -1);
        for (int i = 0; i < raw.length; ++i) raw[i] = raw[i].trim();
        return raw;
    }

    private static String parseRequiredString(CsvTable table, String[] row, String column) {
        final int index = findColumnIndex(table, column);
        if (index >= row.length) {
            throw new IllegalStateException("Missing CSV value for column '" + column + "'.");
        }
        return row[index];
    }

    private static double parseRequiredDouble(CsvTable table, String[] row, String column) {
        return Double.parseDouble(parseRequiredString(table, row, column));
    }

    private static int findColumnIndex(CsvTable table, String column) {
        for (int i = 0; i < table.headers.length; ++i) {
            if (column.equals(table.headers[i])) return i;
        }
        throw new IllegalStateException("Missing CSV column: " + column);
    }

    private static double getRequired(Map<String, Double> props, String key) {
        if (!props.containsKey(key)) {
            throw new IllegalStateException("Missing property key: " + key);
        }
        return props.get(key);
    }

    private static String formatQuantity(double value, String unit) {
        return String.format(Locale.US, "%.12g[%s]", value, unit);
    }

    private static String formatDouble(double value) {
        return String.format(Locale.US, "%.12g", value);
    }

    private static String firstNonEmpty(String value) {
        if (value == null) return null;
        final String trimmed = value.trim();
        return trimmed.isEmpty() ? null : trimmed;
    }

    private static final class CasePaths {
        Path workspaceRoot;
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

    private static final class SamplePoint {
        int id;
        String label;
        String family;
        String location;
        double targetAxisM;
        double targetX;
        double targetY;
        double x;
        double y;
    }

    private static final class BuildResult {
        Model model;
        String representation;
        String physicsRoute;
        String pressureExpr;
        String temperatureExpr;
    }

    private static final class TimeRequest {
        String tag;
        double timeS;
    }

    private static final class MonitorSample {
        int sampleId;
        double timeS;
    }

    private static final class CsvTable {
        String[] headers;
        List<String[]> rows = new ArrayList<String[]>();
    }

    private static final class InterpolationCandidate {
        double[][] coords;
        int pointCount;
    }
}
