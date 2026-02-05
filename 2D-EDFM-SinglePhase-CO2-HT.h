#pragma once
#include <chrono>
#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include "MeshManager.h"



int EDFM_SinglePhase_CO2_HT_2D()
{
    /*-------------------------几何参数及网格设置与划分-----------------------*/
    double lengthX = 1, lengthY = 1, lengthZ = 0;
    int sectionNumX = 20, sectionNumY = 20, sectionNumZ = 0;
    bool usePrism = true;       ///  3D情况usePrism=true 进行“扫掠 ”；
    bool useQuadBase = false;    /// 在 2D 情况：false 为非结构化三角, true  为非结构化四边形； 在 3D 扫掠模式：仅用来选择底面网格类型
   
    std::cout << "\n ==========Mesh: BUILDING ==========\n";
	MeshManager mgr(lengthX, lengthY, lengthZ, sectionNumX, sectionNumY, sectionNumZ, usePrism, useQuadBase);
	mgr.BuildSolidMatrixGrid(NormalVectorCorrectionMethod::OverRelaxed); //这里输入面法矢量修正方法；其中MinimumCorrection-最小修正值法；OrthogonalCorrection-正交修正法；OverRelaxed-超松弛修正法  当前三维几何还不能计算几何和非正交性
	std::cout << " \n ==========Mesh :BUILDING COMPLETED ==========\n";

    /*--------------------------------------------------------------------------*/

    /*-------------------------裂缝添加与裂缝网格划分--------------------------=*/
	mgr.setDFNRandomSeed(12345); //设置随机数种子
    mgr.generateDFN(/*N=*/10,/*minPoint=*/{ 0.0,0.0,0.0 },/*maxPoint=*/{ 1.0,1.0,0.0 },/*Lmin=*/0.5,/*Lmax=*/1.4,/*alpha=*/0,/*kappa=*/0,/*avoidOverlap=*/true);

    mgr.setDistanceMetric(DistanceMetric::CrossAwareGauss);// 设置距离度量方式AreaWeight CellCenter NodeAverage CrossAwareGauss
	mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GlobalAABB); 
    

	// 导出网格、裂缝信息在matlab中进行可视化
    mgr.exportMesh("mesh");
    mgr.exportFractures("fractures");
	/*--------------------------------------------------------------------------*/



	return 0;

}
