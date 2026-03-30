import com.comsol.model.*;
import com.comsol.model.util.*;

public class PressureOnlyProbe {
  private static String describe(double[][] data) {
    if (data == null) return "null";
    int rows = data.length;
    int cols = rows > 0 && data[0] != null ? data[0].length : 0;
    return rows + "x" + cols;
  }
  public static void main(String[] args) {
    int exitCode = 0;
    try {
      ModelUtil.initStandalone(false);
      Model model = ModelUtil.create("m");
      model.param().set("Lx", "400[m]");
      model.param().set("Ly", "40[m]");
      model.param().set("pLeft", "1.2e7[Pa]");
      model.param().set("pRight", "8e6[Pa]");
      model.param().set("pInit", "1e7[Pa]");
      model.param().set("storage", "700*0.1*5e-9[1/Pa]");
      model.param().set("mob", "700*1e-13[m^2]/(6e-5[Pa*s])");
      model.component().create("comp1", true);
      model.component("comp1").geom().create("geom1", 2);
      model.component("comp1").mesh().create("mesh1");
      model.component("comp1").geom("geom1").create("r1", "Rectangle");
      model.component("comp1").geom("geom1").feature("r1").set("size", new String[]{"Lx","Ly"});
      model.component("comp1").geom("geom1").run();
      model.component("comp1").selection().create("leftBnd", "Box");
      model.component("comp1").selection("leftBnd").set("entitydim", "1");
      model.component("comp1").selection("leftBnd").set("condition", "inside");
      model.component("comp1").selection("leftBnd").set("xmin", "-1e-6");
      model.component("comp1").selection("leftBnd").set("xmax", "1e-6");
      model.component("comp1").selection("leftBnd").set("ymin", "-1e-6");
      model.component("comp1").selection("leftBnd").set("ymax", "Ly+1e-6");
      model.component("comp1").selection().create("rightBnd", "Box");
      model.component("comp1").selection("rightBnd").set("entitydim", "1");
      model.component("comp1").selection("rightBnd").set("condition", "inside");
      model.component("comp1").selection("rightBnd").set("xmin", "Lx-1e-6");
      model.component("comp1").selection("rightBnd").set("xmax", "Lx+1e-6");
      model.component("comp1").selection("rightBnd").set("ymin", "-1e-6");
      model.component("comp1").selection("rightBnd").set("ymax", "Ly+1e-6");
      System.out.println("left_count=" + model.component("comp1").selection("leftBnd").entities(1).length);
      System.out.println("right_count=" + model.component("comp1").selection("rightBnd").entities(1).length);
      model.component("comp1").physics().create("pde1", "CoefficientFormPDE", "geom1");
      model.component("comp1").physics("pde1").field("dimensionless").field("p");
      model.component("comp1").physics("pde1").field("dimensionless").component(new String[]{"p"});
      model.component("comp1").physics("pde1").feature("cfeq1").set("ea", "0");
      model.component("comp1").physics("pde1").feature("cfeq1").set("da", "storage");
      model.component("comp1").physics("pde1").feature("cfeq1").set("c", "mob");
      model.component("comp1").physics("pde1").feature("cfeq1").set("f", "0");
      model.component("comp1").physics("pde1").feature("init1").set("p", "pInit");
      model.component("comp1").physics("pde1").create("dirL", "DirichletBoundary", 1);
      model.component("comp1").physics("pde1").feature("dirL").selection().named("leftBnd");
      model.component("comp1").physics("pde1").feature("dirL").set("r", "pLeft");
      model.component("comp1").physics("pde1").create("dirR", "DirichletBoundary", 1);
      model.component("comp1").physics("pde1").feature("dirR").selection().named("rightBnd");
      model.component("comp1").physics("pde1").feature("dirR").set("r", "pRight");
      model.component("comp1").mesh("mesh1").automatic(true);
      model.study().create("std1");
      model.study("std1").create("time", "Transient");
      model.study("std1").feature("time").activate("pde1", true);
      model.study("std1").feature("time").set("tlist", "range(0,2e5,2e7)");
      model.study("std1").createAutoSequences("sol");
      model.study("std1").run();
      model.result().numerical().create("interp1", "Interp");
      model.result().numerical("interp1").set("expr", new String[]{"p"});
      model.result().numerical("interp1").set("edim", 2);
      model.result().numerical("interp1").set("data", "dset1");
      model.result().numerical("interp1").set("coord", new double[][]{{100,200,300},{20,20,20}});
      model.result().numerical("interp1").set("t", new double[]{2e7});
      model.result().numerical("interp1").run();
      double[][] a = model.result().numerical("interp1").getReal(true);
      double[][] b = model.result().numerical("interp1").getReal(false);
      double[][] c = model.result().numerical("interp1").getReal();
      double[][] d = model.result().numerical("interp1").getData(0);
      System.out.println("getReal(true)=" + describe(a));
      System.out.println("getReal(false)=" + describe(b));
      System.out.println("getReal()=" + describe(c));
      System.out.println("getData(0)=" + describe(d));
      double[][] values = a.length > 0 ? a : (b.length > 0 ? b : (c.length > 0 ? c : d));
      for (int i = 0; i < values.length; ++i) {
        for (int j = 0; j < values[i].length; ++j) {
          System.out.println("v[" + i + "][" + j + "]=" + values[i][j]);
        }
      }
    } catch (Exception ex) {
      exitCode = 1;
      ex.printStackTrace(System.err);
    } finally {
      try { ModelUtil.clear(); } catch (Exception ignored) {}
      try { ModelUtil.disconnect(); } catch (Exception ignored) {}
      System.exit(exitCode);
    }
  }
}
