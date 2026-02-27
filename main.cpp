
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
#include "Test_DOFMapper.h"
#include "Test_ADVar.h"

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
    
    RunTest_DFN_Initialization_And_Viz();             //测试包含DFN网络的初始化场的可视化验证                                结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\FieldOperator 通过tecplot进行可视化

    //RunTest_Property_Accuracy_Sweep();                //测试物性计算准确性扫描 (One-by-One Verify)                           结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\PropertyTest 通过Property_Accuracy_Sweep_HD.csv来表征

    //RunTest_BoundaryCondition_Export();					 //测试3D_EDFM中边界条件设置的正确性,保持均匀分布                                   结果与可视化脚本置于：\2D-Unstr-Quadrilateral-EDFM\Test\BoundaryConditionSetup    
    //return 0;

    ////测试梯度算子
    //     run_Benchmark_2D_EDFM_Grad();
    //     run_Benchmark_3D_EDFM_Grad();
    //     // 2. 进阶工况 A: 二次非线性场 (Quadratic)
    //    // P = x^2 + y^2 + z^2
    //    // Grad = (2x, 2y, 2z)
    //     run_Advanced_Accuracy_Test(
    //         "Case A: Quadratic Field (P = x^2 + y^2 + z^2)",
    //         [](const Vector& p) { return p.m_x * p.m_x + p.m_y * p.m_y + p.m_z * p.m_z; },
    //         [](const Vector& p) { return Vector(2.0 * p.m_x, 2.0 * p.m_y, 2.0 * p.m_z); },
    //         true // expect non-zero error
    //     );

    //     // 3. 进阶工况 B: 交叉项场 (Cross-term)
    //     // P = xy + yz
    //     // Grad = (y, x+z, y)
    //     // 这也是二次的，但包含耦合项
    //     run_Advanced_Accuracy_Test(
    //         "Case B: Mixed Quadratic (P = xy + yz)",
    //         [](const Vector& p) { return p.m_x * p.m_y + p.m_y * p.m_z; },
    //         [](const Vector& p) { return Vector(p.m_y, p.m_x + p.m_z, p.m_y); },
    //         true
    //     );

    //     // 4. 进阶工况 C: 三角函数场 (Trigonometric)
    //     // P = sin(x) + cos(y)
    //     // Grad = (cos(x), -sin(y), 0)
    //     run_Advanced_Accuracy_Test(
    //         "Case C: Trigonometric (P = sin(x) + cos(y))",
    //         [](const Vector& p) { return std::sin(p.m_x) + std::cos(p.m_y); },
    //         [](const Vector& p) { return Vector(std::cos(p.m_x), -std::sin(p.m_y), 0.0); },
    //         true
    //     );

    // return 0;

    //测试DOFmap
    //RunDOFMapperVerification();
    //return 0;

	//测试ADVar
   /* Run_ADVar_Comprehensive_Tests();
    return 0;*/
    //return  EDFM_withFracture_Geomtry(); 
    //return run_IMPES_Iteration_TwoPhase_WellCase();
    //return SinglePhase_CO2_TH_withWell_reviese();
}
