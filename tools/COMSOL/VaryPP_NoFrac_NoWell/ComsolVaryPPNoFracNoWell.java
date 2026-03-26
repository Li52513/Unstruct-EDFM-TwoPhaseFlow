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
import java.util.List;
import java.util.Locale;

public final class ComsolVaryPPNoFracNoWell {
    private static final String DEFAULT_CASE_DIR =
        "Test/Transient/FullCaseTest/H_CO2_VaryPP/h_co2_varypp_nofrac_nowell";
    private static final double LX_M = 400.0;
    private static final double LY_M = 40.0;
    private static final double P_LEFT_PA = 12.0e6;
    private static final double P_RIGHT_PA = 8.0e6;
    private static final double P_INIT_PA = 10.0e6;
    private static final double T_INIT_K = 360.0;
    private static final double DELTA_P_PA = P_LEFT_PA - P_RIGHT_PA;
    private static final double MAIN_HMAX_M = 1.0;
    private static final double MAIN_HMIN_M = 0.1;
    private static final double FINE_HMAX_M = 0.5;
    private static final double FINE_HMIN_M = 0.05;
    private static final String TIME_LIST = "range(0,500,100000)";

    private ComsolVaryPPNoFracNoWell() {}

    public static void main(String[] args) {
        final CasePaths paths = resolvePaths(args);
        final boolean skipFineCheck = shouldSkipFineCheck();
        int exitCode = 0;

        try {
            ensureDirectories(paths);
            validateInputs(paths);

            final List<SamplePoint> profileStations = readSamplePoints(paths.profileStationsCsv, "station_id");
            final List<SamplePoint> monitorPoints = readSamplePoints(paths.monitorPointsCsv, "point_id");
            final List<TimeRequest> profileTimes = readProfileSchedule(paths.profileScheduleCsv);
            final double[] monitorTimes = readTimeSchedule(paths.monitorScheduleCsv);

            ModelUtil.initStandalone(false);
            ModelUtil.loadPreferences();
            ModelUtil.showProgress(paths.comsolOutputDir.resolve("comsol_progress.log").toString());

            final Model mainModel = buildAndSolveModel(paths, "ModelMain", "main", MAIN_HMAX_M, MAIN_HMIN_M);
            writeProfileReferenceCsvs(mainModel, profileStations, profileTimes, paths);
            writeMonitorReferenceCsv(mainModel, monitorPoints, monitorTimes, paths);
            mainModel.save(paths.comsolOutputDir.resolve("comsol_model.mph").toString());

            if (!skipFineCheck) {
                final Model fineModel = buildAndSolveModel(paths, "ModelFine", "fine", FINE_HMAX_M, FINE_HMIN_M);
                fineModel.save(paths.comsolOutputDir.resolve("comsol_model_refined.mph").toString());
                writeMeshCheck(mainModel, fineModel, profileStations, profileTimes, paths);
            } else {
                writeSkippedFineCheck(paths);
            }

            System.out.println("[COMSOL] VaryPP reference generation completed.");
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

    private static boolean shouldSkipFineCheck() {
        final String flag = System.getenv("COMSOL_SKIP_FINE_CHECK");
        return flag != null && !flag.trim().isEmpty() && !"0".equals(flag.trim());
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
        paths.comsolInputDir = caseDir.resolve("reference").resolve("comsol_input");
        paths.comsolOutputDir = caseDir.resolve("reference").resolve("comsol");
        paths.profileStationsCsv = caseDir.resolve("profile_station_definitions.csv");
        paths.monitorPointsCsv = caseDir.resolve("monitor_point_definitions.csv");
        paths.profileScheduleCsv = caseDir.resolve("profile_report_schedule.csv");
        paths.monitorScheduleCsv = caseDir.resolve("monitor_sample_schedule.csv");
        paths.lambdaTable = paths.comsolInputDir.resolve("lambda_fun.txt");
        paths.storageTable = paths.comsolInputDir.resolve("S_fun.txt");
        paths.rhoTable = paths.comsolInputDir.resolve("rho_fun.txt");
        paths.muTable = paths.comsolInputDir.resolve("mu_fun.txt");
        paths.cfTable = paths.comsolInputDir.resolve("cf_fun.txt");
        return paths;
    }

    private static void ensureDirectories(CasePaths paths) throws IOException {
        Files.createDirectories(paths.comsolOutputDir);
    }

    private static void validateInputs(CasePaths paths) {
        final List<Path> required = Arrays.asList(
            paths.profileStationsCsv,
            paths.monitorPointsCsv,
            paths.profileScheduleCsv,
            paths.monitorScheduleCsv,
            paths.lambdaTable,
            paths.storageTable,
            paths.rhoTable,
            paths.muTable,
            paths.cfTable
        );
        for (Path requiredPath : required) {
            if (!Files.isRegularFile(requiredPath)) {
                throw new IllegalStateException("Missing required input: " + requiredPath);
            }
        }
    }

    private static Model buildAndSolveModel(
        CasePaths paths,
        String modelTag,
        String labelSuffix,
        double hmaxM,
        double hminM
    ) {
        final Model model = ModelUtil.create(modelTag);
        model.modelPath(paths.comsolOutputDir.toString());
        model.label("ComsolVaryPPNoFracNoWell_" + labelSuffix + ".mph");

        model.param().set("Lx", formatQuantity(LX_M, "m"));
        model.param().set("Ly", formatQuantity(LY_M, "m"));
        model.param().set("pLeft", formatQuantity(P_LEFT_PA, "Pa"));
        model.param().set("pRight", formatQuantity(P_RIGHT_PA, "Pa"));
        model.param().set("pInit", formatQuantity(P_INIT_PA, "Pa"));
        model.param().set("tInit", formatQuantity(T_INIT_K, "K"));

        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 2);
        model.component("comp1").mesh().create("mesh1");

        model.component("comp1").geom("geom1").create("r1", "Rectangle");
        model.component("comp1").geom("geom1").feature("r1").set("base", "corner");
        model.component("comp1").geom("geom1").feature("r1").set("pos", new String[]{"0", "0"});
        model.component("comp1").geom("geom1").feature("r1").set("size", new String[]{"Lx", "Ly"});
        model.component("comp1").geom("geom1").run();

        createBoundarySelections(model);
        createInterpolationFunctions(model, paths);
        createCoefficientPde(model);
        configureMesh(model, hmaxM, hminM);
        configureStudy(model);

        model.study("std1").run();
        return model;
    }

    private static void createBoundarySelections(Model model) {
        final double tol = 1.0e-6;

        model.component("comp1").selection().create("leftBnd", "Box");
        model.component("comp1").selection("leftBnd").set("entitydim", "1");
        model.component("comp1").selection("leftBnd").set("condition", "inside");
        model.component("comp1").selection("leftBnd").set("xmin", formatDouble(-tol));
        model.component("comp1").selection("leftBnd").set("xmax", formatDouble(tol));
        model.component("comp1").selection("leftBnd").set("ymin", formatDouble(-tol));
        model.component("comp1").selection("leftBnd").set("ymax", formatDouble(LY_M + tol));
        model.component("comp1").selection("leftBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("leftBnd").set("zmax", formatDouble(tol));

        model.component("comp1").selection().create("rightBnd", "Box");
        model.component("comp1").selection("rightBnd").set("entitydim", "1");
        model.component("comp1").selection("rightBnd").set("condition", "inside");
        model.component("comp1").selection("rightBnd").set("xmin", formatDouble(LX_M - tol));
        model.component("comp1").selection("rightBnd").set("xmax", formatDouble(LX_M + tol));
        model.component("comp1").selection("rightBnd").set("ymin", formatDouble(-tol));
        model.component("comp1").selection("rightBnd").set("ymax", formatDouble(LY_M + tol));
        model.component("comp1").selection("rightBnd").set("zmin", formatDouble(-tol));
        model.component("comp1").selection("rightBnd").set("zmax", formatDouble(tol));
    }

    private static void createInterpolationFunctions(Model model, CasePaths paths) {
        createInterpolationFunction(model, "rho_fun", paths.rhoTable, "Pa", "kg/m^3");
        createInterpolationFunction(model, "mu_fun", paths.muTable, "Pa", "Pa*s");
        createInterpolationFunction(model, "cf_fun", paths.cfTable, "Pa", "1/Pa");
        createInterpolationFunction(model, "lambda_fun", paths.lambdaTable, "Pa", "kg/(m*s*Pa)");
        createInterpolationFunction(model, "S_fun", paths.storageTable, "Pa", "kg/(m^3*Pa)");
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

    private static void createCoefficientPde(Model model) {
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

    private static void configureMesh(Model model, double hmaxM, double hminM) {
        model.component("comp1").mesh("mesh1").feature().create("size1", "Size");
        model.component("comp1").mesh("mesh1").feature("size1").set("custom", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hmaxactive", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hminactive", true);
        model.component("comp1").mesh("mesh1").feature("size1").set("hmax", formatQuantity(hmaxM, "m"));
        model.component("comp1").mesh("mesh1").feature("size1").set("hmin", formatQuantity(hminM, "m"));
        model.component("comp1").mesh("mesh1").feature().create("ftri1", "FreeTri");
        model.component("comp1").mesh("mesh1").run();
    }

    private static void configureStudy(Model model) {
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").activate("c", true);
        model.study("std1").feature("time").set("tlist", TIME_LIST);
    }

    private static void writeProfileReferenceCsvs(
        Model model,
        List<SamplePoint> profileStations,
        List<TimeRequest> profileTimes,
        CasePaths paths
    ) throws IOException {
        final double[][] coords = buildCoordinateMatrix(profileStations);
        final double[] times = new double[profileTimes.size()];
        for (int i = 0; i < profileTimes.size(); ++i) {
            times[i] = profileTimes.get(i).timeS;
        }

        final double[][] values = evaluateAtPoints(model, "interpProfileMain", coords, times);
        for (int timeIndex = 0; timeIndex < profileTimes.size(); ++timeIndex) {
            final TimeRequest request = profileTimes.get(timeIndex);
            final Path csvPath = paths.comsolOutputDir.resolve("comsol_profile_" + request.tag + ".csv");
            try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
                out.write("station_id,label,x_m,y_m,time_s,p_ref_pa");
                out.newLine();
                for (int pointIndex = 0; pointIndex < profileStations.size(); ++pointIndex) {
                    final SamplePoint point = profileStations.get(pointIndex);
                    out.write(String.format(
                        Locale.US,
                        "%d,%s,%.12f,%.12f,%.12f,%.12f",
                        point.id, point.label, point.x, point.y, request.timeS, values[timeIndex][pointIndex]
                    ));
                    out.newLine();
                }
            }
        }
    }

    private static void writeMonitorReferenceCsv(
        Model model,
        List<SamplePoint> monitorPoints,
        double[] monitorTimes,
        CasePaths paths
    ) throws IOException {
        final double[][] coords = buildCoordinateMatrix(monitorPoints);
        final double[][] values = evaluateAtPoints(model, "interpMonitorMain", coords, monitorTimes);
        final Path csvPath = paths.comsolOutputDir.resolve("comsol_monitor_timeseries.csv");
        try (BufferedWriter out = Files.newBufferedWriter(csvPath, StandardCharsets.UTF_8)) {
            out.write("time_s");
            for (SamplePoint point : monitorPoints) {
                out.write(",p_ref_" + point.label);
            }
            out.newLine();
            for (int timeIndex = 0; timeIndex < monitorTimes.length; ++timeIndex) {
                out.write(String.format(Locale.US, "%.12f", monitorTimes[timeIndex]));
                for (int pointIndex = 0; pointIndex < monitorPoints.size(); ++pointIndex) {
                    out.write(String.format(Locale.US, ",%.12f", values[timeIndex][pointIndex]));
                }
                out.newLine();
            }
        }
    }

    private static void writeMeshCheck(
        Model mainModel,
        Model fineModel,
        List<SamplePoint> profileStations,
        List<TimeRequest> profileTimes,
        CasePaths paths
    ) throws IOException {
        final TimeRequest finalRequest = profileTimes.get(profileTimes.size() - 1);
        final double[][] coords = buildCoordinateMatrix(profileStations);
        final double[][] mainValues = evaluateAtPoints(
            mainModel, "interpProfileMainCheck", coords, new double[]{finalRequest.timeS});
        final double[][] fineValues = evaluateAtPoints(
            fineModel, "interpProfileFineCheck", coords, new double[]{finalRequest.timeS});

        double sumSq = 0.0;
        double maxAbs = 0.0;
        for (int i = 0; i < profileStations.size(); ++i) {
            final double absErr = Math.abs(mainValues[0][i] - fineValues[0][i]);
            sumSq += absErr * absErr;
            maxAbs = Math.max(maxAbs, absErr);
        }

        final double l2Norm = Math.sqrt(sumSq / Math.max(profileStations.size(), 1)) / Math.abs(DELTA_P_PA);
        final double linfNorm = maxAbs / Math.abs(DELTA_P_PA);
        final boolean ok = l2Norm < 1.0e-3;

        final Path reportPath = paths.comsolOutputDir.resolve("comsol_reference_mesh_check.txt");
        try (BufferedWriter out = Files.newBufferedWriter(reportPath, StandardCharsets.UTF_8)) {
            out.write("profile_tag=" + finalRequest.tag);
            out.newLine();
            out.write("time_s=" + formatDouble(finalRequest.timeS));
            out.newLine();
            out.write("main_mesh_hmax_m=" + formatDouble(MAIN_HMAX_M));
            out.newLine();
            out.write("main_mesh_hmin_m=" + formatDouble(MAIN_HMIN_M));
            out.newLine();
            out.write("fine_mesh_hmax_m=" + formatDouble(FINE_HMAX_M));
            out.newLine();
            out.write("fine_mesh_hmin_m=" + formatDouble(FINE_HMIN_M));
            out.newLine();
            out.write("profile_l2_norm=" + formatDouble(l2Norm));
            out.newLine();
            out.write("profile_linf_norm=" + formatDouble(linfNorm));
            out.newLine();
            out.write("reference_acceptable=" + (ok ? "true" : "false"));
            out.newLine();
        }
    }

    private static void writeSkippedFineCheck(CasePaths paths) throws IOException {
        final Path reportPath = paths.comsolOutputDir.resolve("comsol_reference_mesh_check.txt");
        try (BufferedWriter out = Files.newBufferedWriter(reportPath, StandardCharsets.UTF_8)) {
            out.write("fine_check_skipped=true");
            out.newLine();
            out.write("reason=COMSOL_SKIP_FINE_CHECK was set.");
            out.newLine();
        }
    }

    private static double[][] evaluateAtPoints(
        Model model,
        String numericalTag,
        double[][] coords,
        double[] times
    ) {
        model.result().numerical().create(numericalTag, "Interp");
        final NumericalFeature interp = model.result().numerical(numericalTag);
        interp.set("expr", new String[]{"p"});
        interp.set("edim", 2);
        interp.set("coorderr", "on");
        interp.set("unit", new String[]{"Pa"});
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

    private static double[][] buildCoordinateMatrix(List<SamplePoint> points) {
        final double[][] coords = new double[2][points.size()];
        for (int i = 0; i < points.size(); ++i) {
            coords[0][i] = points.get(i).x;
            coords[1][i] = points.get(i).y;
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

    private static List<SamplePoint> readSamplePoints(Path csvPath, String idColumn) {
        final CsvTable table = readCsv(csvPath);
        final List<SamplePoint> points = new ArrayList<SamplePoint>();
        for (String[] row : table.rows) {
            final SamplePoint point = new SamplePoint();
            point.id = (int) parseRequiredDouble(table, row, idColumn);
            point.label = parseRequiredString(table, row, "label");
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

    private static double[] readTimeSchedule(Path csvPath) {
        final CsvTable table = readCsv(csvPath);
        final double[] times = new double[table.rows.size()];
        for (int i = 0; i < table.rows.size(); ++i) {
            times[i] = parseRequiredDouble(table, table.rows.get(i), "target_time_s");
        }
        return times;
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

    private static final class CasePaths {
        Path workspaceRoot;
        Path caseDir;
        Path comsolInputDir;
        Path comsolOutputDir;
        Path profileStationsCsv;
        Path monitorPointsCsv;
        Path profileScheduleCsv;
        Path monitorScheduleCsv;
        Path lambdaTable;
        Path storageTable;
        Path rhoTable;
        Path muTable;
        Path cfTable;
    }

    private static final class SamplePoint {
        int id;
        String label;
        double x;
        double y;
    }

    private static final class TimeRequest {
        String tag;
        double timeS;
    }

    private static final class InterpolationCandidate {
        String label;
        double[][] coords;
        int pointCount;
    }

    private static final class CsvTable {
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
}
