//#include <iostream>
//#include <memory>
//#include <utility>
//#include <vector>
//#include "MeshManager.h"
//#include "FieldRegistry.h"
//#include "FaceFieldRegistry.h"
//#include "PhysicalPropertiesManager.h"
//#include "BCAdapter.h"
//#include "TemperatureBCAdapter.h"
//#include "PressureBC.h"
//#include "TemperatureBC.h"
//#include "CapRelPerm.h"
//#include "DiffusionCentral.h"
//#include "ConvectionUpwind.h"
//#include "ConvectionUpwind_Flux.h"
//#include "Timeterm_BDF.h"
//#include "Solver_AssemblerCOO.h"
//#include "LinearSolver_Eigen.h"
//#include "FVM_WellDOF.h"
//#include "FVM_WellCoupling.h"
//#include "FVM_Peaceman.h"
//#include "FVM_SourceTerm_WellHeat.h"
//#include "Solver_TimeLoopUtils.h"
//
//
//
//
//int main() 
//{
//    //定义几何区域参数
//    double lengthX = 100, lengthY = 100, lengthZ = 0;
//    //确定网格划分策略及参数
//    int sectionNumX = 20, sectionNumY = 20, sectionNumZ = 0;
//    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
//    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
//
//    // 1) 构造并预处理网格
//    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
//    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
//
//    // 1) 构造并预处理网格
//    MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
//    mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OrthogonalCorrection);
//    // 2) 网格单元场和网格面场注册
//    FieldRegistry reg;
//    FaceFieldRegistry freg;
//
//    //3) 给定地层初始参数
//    InitFields ic;
//    {
//        //==压力==//
//        ic.p_w0 = 1e7;
//        ic.dp_wdx = 0;
//        ic.dp_wdy = 0;
//        ic.dp_wdz = 0;
//
//        //==温度==//
//        ic.T0 = 593.15;
//        ic.dTdx = 0;
//        ic.dTdy = 0;
//        ic.dTdz = 0;
//
//        //==初始水相饱和度==//
//        double sgr = 0.05;  //CO2的残余饱和度
//        ic.swr = 1 - sgr;   //水相初始饱和度
//    }
//
//    //4) 主变量场p_w,s_w,T的注册
//    {
//        Initializer::createPrimaryFields_TwoPhase_HT_IMPES(mgr.mesh(), reg);
//        Initializer::fillBaseDistributions_TwoPhase_HT_IMPES(mgr.mesh(), reg, ic);
//    }
//
//    // 5) 物性参数注册
//    PhysicalPropertiesManager ppm;
//
//
//
//
//
//
//
//    return 0;
//
//}
