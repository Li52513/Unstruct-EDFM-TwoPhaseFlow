import com.comsol.model.*;
import com.comsol.model.util.*;

public class SolverProbe {
  public static void main(String[] args) {
    int exitCode = 0;
    try {
      ModelUtil.initStandalone(false);
      Model model = ModelUtil.create("m");
      model.component().create("comp1", true);
      model.component("comp1").geom().create("geom1", 2);
      model.component("comp1").geom("geom1").create("r1", "Rectangle");
      model.component("comp1").geom("geom1").run();
      model.component("comp1").mesh().create("mesh1");
      model.component("comp1").mesh("mesh1").automatic(true);
      model.component("comp1").physics().create("pde1", "CoefficientFormPDE", "geom1");
      model.component("comp1").physics("pde1").field("dimensionless").field("u");
      model.component("comp1").physics("pde1").field("dimensionless").component(new String[]{"u"});
      model.component("comp1").physics("pde1").feature("cfeq1").set("da", "1");
      model.component("comp1").physics("pde1").feature("cfeq1").set("c", "1");
      model.component("comp1").physics("pde1").feature("cfeq1").set("f", "0");
      model.component("comp1").physics("pde1").feature("init1").set("u", "0");
      model.study().create("std1");
      model.study("std1").create("time", "Transient");
      model.study("std1").feature("time").activate("pde1", true);
      model.study("std1").feature("time").set("tlist", "range(0,1,2)");
      model.study("std1").createAutoSequences("sol");
      for (String stag : model.sol().tags()) {
        System.out.println("SOL=" + stag);
        for (String ftag : model.sol(stag).feature().tags()) {
          System.out.println("  FEAT=" + ftag);
          try {
            for (String sftag : model.sol(stag).feature(ftag).feature().tags()) {
              System.out.println("    SUB=" + sftag);
              try {
                for (String ssftag : model.sol(stag).feature(ftag).feature(sftag).feature().tags()) {
                  System.out.println("      SUBSUB=" + ssftag);
                }
              } catch (Exception ignored3) {}
            }
          } catch (Exception ignored2) {}
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
