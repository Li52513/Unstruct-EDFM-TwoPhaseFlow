#include   "SinglePhase_CO2_TH_withoutWell.h"
#include  "SinglePhase_CO2_TH_withWell.h"
#include    "TestCapillaryTerm_v2.h"
#include    "run_IMPES_Iteration_TimeTerm_AnalyticalTest.h"
#include    "run_IMPES_Iteration_TwoPhase_BL_Numerical.h"
#include    "run_IMPES_Iteration_TwoPhase_withWell.h"
#include    "run_FC_P_IMPES_I_TwoPhase_withWell.h"
#include    "SinglePhase_TH_CO2_withwell_revised.h"
#include    "EDFM_withFracture_Geomtrycreate.h"
#include    "2D-EDFM-SinglePhase-CO2-HT.h"
#include    "3D_Improved_EDFM_test.h"
#include    "2D_EDFM_test.h"
#include    "2D_EDFM_MeshTest_Benchmark.h"

#include    "3D_TwistFracIntersectionTest.h"
#include	"3D_EDFM_GeomTest.h"
#include    "3D_EDFM_MeshTest_Benchmark.h"
#include    "3D_EDFM_Transmissibility_test.h"
#include    "3D_InitializerTest.h"
#include    "3D_InitializerTest_multiFrac.h"
#include    "3D_PropTest_HD.h"

#include    "3D_BoundarySetup_Test.h"

#include "test_FVM_Grad_Benchmark.h"

int main()
{
    // =========================================================
    // 【Chapter1】EDFM裂缝网格快速生成方法
    // =========================================================
    
    //return EDFM_test_2D();                            //全流程测试包含单根裂缝的2D_EDFM四种加速算法的性能                      结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\MeshTest\2D_EDFM\SingleFrac_Test
     
    //return EDFM_DFN_test_2D();                        //全流程测试包含DFN的2D_EDFM四种加速算法的性能                           结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\MeshTest\2D_EDFM\DFN_Test
    
    //return TwistFracIntersectionTest_3D();            //3D_EDFM中扭曲裂缝相交测试

    //return Improved_EDFM_test_3D();                   //全流程测试3D_Improved_EDFM四种加速算法的性能以及交互多边形生成的准确性  结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\MeshTest\3D_EDFM

    //return RunBenchmark_Step1_Distance_Accuracy();    //测试3D_EDFM中基岩网格单元到裂缝交互单元的距离d的准确性
    
    //RunBenchmark_Transmissibility_Static();           //测试 NNC 和 FF 的静态传导率 (Flow & Heat) 是否符合解析解
    //return 0;

    //return EDFM_DFN_Geomtest_2D();
    //return Improved_3D_EDFM_MeshTest();

    // =========================================================
    // 【Chapter2】EDFM-超临界CO2单相渗流换热模拟研究
    // =========================================================
    //RunTest_Initialization_And_Viz();                 //测试包含单根裂缝的初始化场的可视化验证                               结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\FieldOperator 通过tecplot进行可视化
    
    //RunTest_DFN_Initialization_And_Viz();             //测试包含DFN网络的初始化场的可视化验证                                结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\FieldOperator 通过tecplot进行可视化

    //RunTest_Property_Accuracy_Sweep();                //测试物性计算准确性扫描 (One-by-One Verify)                           结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\PropertyTest 通过Property_Accuracy_Sweep_HD.csv来表征

    //RunTest_BoundaryCondition_Export();					 //测试3D_EDFM中边界条件设置的正确性,保持均匀分布                                   结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\BoundaryConditionSetup    
    //return 0;

    //测试梯度算子
         run_Benchmark_2D_EDFM_Grad();
         run_Benchmark_3D_EDFM_Grad();
     return 0;
    
    //return  EDFM_withFracture_Geomtry(); 
    //return run_IMPES_Iteration_TwoPhase_WellCase();
    //return SinglePhase_CO2_TH_withWell_reviese();
}
