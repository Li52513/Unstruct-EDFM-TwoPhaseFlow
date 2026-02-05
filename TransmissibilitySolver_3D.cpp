#include "TransmissibilitySolver_3D.h"
#include "SolverContrlStrName.h" // 包含 PhysicalProperties_string 定义
#include "2D_FracIndex.h" 

#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_map>

using namespace PhysicalProperties_string;

// =========================================================
// 内部辅助函数
// =========================================================

/**
 * @brief 计算点 P 到线段 AB (由 A, B, length 定义) 的垂直距离
 * @details 用于 F-F 传导率计算中的 d_i
 */
static double _pointToLineDistance(const Vector& P, const Vector& A, const Vector& B, double length)
{
    if (length < 1e-12) return (P - A).Mag(); // 退化为点距离

    Vector AP = P - A;
    Vector AB = B - A;
    
    // 直接使用 Variable3D 重载的 & 运算符进行叉乘，并计算模长
    // 物理意义：平行四边形面积
    double areaParallelogram = (AP & AB).Mag();
    
    return areaParallelogram / length;
}

// =========================================================================
// 静态传导率计算：NNC (Matrix - Fracture)
// =========================================================================
void TransmissibilitySolver_3D::Calculate_Transmissibility_NNC(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
{
    Rock rockStr;
    FractureProperties_String fracStr;

    const double eps = 1e-30;
    const double lam_fluid_ref = 0.6; // [Static Ref]参考流体导热系数 (W/mK)

    // 获取场数据 (使用 const 版本的 meshMgr 接口)
	auto& Kxx = fieldMgr.getOrCreateMatrixScalar(rockStr.k_xx_tag, 1e-15)->data;    // 基岩渗透率场
	auto& Kyy = fieldMgr.getOrCreateMatrixScalar(rockStr.k_yy_tag, 1e-15)->data;    // 基岩渗透率场
	auto& Kzz = fieldMgr.getOrCreateMatrixScalar(rockStr.k_zz_tag, 1e-15)->data;    // 基岩渗透率场
	auto& LamM = fieldMgr.getOrCreateMatrixScalar(rockStr.lambda_tag, 2.0)->data;   // 基岩导热系数场
	auto& PhiM = fieldMgr.getOrCreateMatrixScalar(rockStr.phi_tag, 0.1)->data;      // 基岩孔隙度场

	auto& KnF = fieldMgr.getOrCreateFractureScalar(fracStr.k_n, 1e-12)->data;       // 裂缝法向渗透率场
	auto& LamF = fieldMgr.getOrCreateFractureScalar(fracStr.lambda_tag, 2.0)->data;	// 裂缝导热系数场
	auto& PhiF = fieldMgr.getOrCreateFractureScalar(fracStr.phi_tag, 1.0)->data;    // 裂缝孔隙度场

	const std::string tag_T_NNC_Flow = "T_NNC_Flow";                                // 渗流传导率标签
	const std::string tag_T_NNC_Heat = "T_NNC_Heat";                                // 传热传导率标签

    auto& T_flow = fieldMgr.createNNCScalar(tag_T_NNC_Flow, 0.0)->data;
    auto& T_heat = fieldMgr.createNNCScalar(tag_T_NNC_Heat, 0.0)->data;

    const auto& pairs = meshMgr.getInteractionPairs();
    const auto& frNet = meshMgr.fracture_network(); 

    auto fracMap = frNet.buildFractureIDMap();

    std::cout << "[Transmissibility] Calculating NNC (M-F) for " << pairs.size() << " pairs..." << std::endl;

    for (size_t i = 0; i < pairs.size(); ++i)
    {
        const auto& pair = pairs[i];

		int mGlobalID = pair.matrixCellGlobalID;                                    // 基岩单元全局ID
		int mIdx = meshMgr.mesh().getCellIndex(mGlobalID);                          // 基岩单元索引
		int fElemID = pair.fracElementGlobalID;     			                    // 裂缝单元全局ID      

        double A = pair.intersectionArea;                                           // 该交互多边形的面积          
		double d = pair.distMatrixToFracPlane;									    // 基岩单元中心到裂缝平面的距离
		Vector n = pair.polygonNormal;              			                    // 多边形法向量（已归一化）

        double w = 0.0;
        auto it = fracMap.find(pair.fracMacroID);
		if (it != fracMap.end()) w = it->second->aperture;                          // 裂缝开度
        else continue;

        // --- Flow ---
        double Km_n = Kxx[mIdx] * n.m_x * n.m_x +
            Kyy[mIdx] * n.m_y * n.m_y +
            Kzz[mIdx] * n.m_z * n.m_z;
        double Kf_n = KnF[fElemID];

        // 保持 4.0 系数 经典文献（如 Karimi-Fard et al. 2004, Moinfar 2013）推导出的形状因子 (Shape Factor) 对应的分母项确实是 $4 K_f$
        double term_flow_m = d / (Km_n + eps);
        double term_flow_f = w / (4.0 * Kf_n + eps);
        T_flow[i] = A / (term_flow_m + term_flow_f + eps);

        // --- Heat ---
        double lam_eff_m = PhiM[mIdx] * lam_fluid_ref + (1.0 - PhiM[mIdx]) * LamM[mIdx];
        double lam_eff_f = PhiF[fElemID] * lam_fluid_ref + (1.0 - PhiF[fElemID]) * LamF[fElemID];

        double term_heat_m = d / (lam_eff_m + eps);
        double term_heat_f = w / (4.0 * lam_eff_f + eps);
        T_heat[i] = A / (term_heat_m + term_heat_f + eps);
    }
}

// =========================================================================
// 静态传导率计算：FF (Fracture - Fracture)
// =========================================================================
void TransmissibilitySolver_3D::Calculate_Transmissibility_FF(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
{
    FractureProperties_String fracStr; // [修改]
    const double eps = 1e-30;
	const double lam_fluid_ref = 0.6; // [Static Ref]参考流体导热系数 (W/mK)

    const auto& frNet = meshMgr.fracture_network();
    const auto& ffIntersections = frNet.ffIntersections;
    const auto& indexer = frNet.fracElemIndex;
    const auto& fracturesVec = frNet.getFractures();

    auto& KtF = fieldMgr.getOrCreateFractureScalar(fracStr.k_t, 1e-12)->data;
    auto& LamF = fieldMgr.getOrCreateFractureScalar(fracStr.lambda_tag, 2.0)->data;
    auto& PhiF = fieldMgr.getOrCreateFractureScalar(fracStr.phi_tag, 1.0)->data;

    const std::string tag_T_FF_Flow = "T_FF_Flow";
    const std::string tag_T_FF_Heat = "T_FF_Heat";

    auto& T_ff_flow = fieldMgr.createFFScalar(tag_T_FF_Flow, 0.0)->data;
    auto& T_ff_heat = fieldMgr.createFFScalar(tag_T_FF_Heat, 0.0)->data;

    std::cout << "[Transmissibility] Calculating FF (F-F) connections..." << std::endl;

	// 构建裂缝ID到索引的映射，便于快速查找
    std::unordered_map<int, size_t> idToIndexMap;
    idToIndexMap.reserve(fracturesVec.size());
    for (size_t k = 0; k < fracturesVec.size(); ++k) {
        idToIndexMap[fracturesVec[k].id] = k;
    }

    size_t connectionIdx = 0;

    for (const auto& interaction : ffIntersections)
    {
        int id1 = interaction.fracID_1;
        int id2 = interaction.fracID_2;

		// 查找裂缝索引
        auto it1 = idToIndexMap.find(id1);
        auto it2 = idToIndexMap.find(id2);

        if (it1 == idToIndexMap.end() || it2 == idToIndexMap.end()) {
            connectionIdx += interaction.segments.size();
            continue;
        }

        size_t idx1 = it1->second;
        size_t idx2 = it2->second;
        const auto& f1 = fracturesVec[idx1];
        const auto& f2 = fracturesVec[idx2];

        // 调和平均开度 (Harmonic Aperture)
        double w_harm = 2.0 * f1.aperture * f2.aperture / (f1.aperture + f2.aperture + eps);

        for (const auto& seg : interaction.segments)
        {
            int localID1 = seg.cellID_1;
            int localID2 = seg.cellID_2;

            size_t globalID1 = indexer.offset[idx1] + localID1;
            size_t globalID2 = indexer.offset[idx2] + localID2;

            if (localID1 >= f1.fracCells.size() || localID2 >= f2.fracCells.size()) {
                connectionIdx++; continue;
            }

            Vector c1 = f1.fracCells[localID1].centroid;
            Vector c2 = f2.fracCells[localID2].centroid;

            double d1 = _pointToLineDistance(c1, seg.start, seg.end, seg.length);
            double d2 = _pointToLineDistance(c2, seg.start, seg.end, seg.length);

            if (d1 < 1e-6) d1 = 1e-6;
            if (d2 < 1e-6) d2 = 1e-6;

            double area = seg.length * w_harm;

            // --- Flow (半传导率调和平均) ---
            double T1_flow = KtF[globalID1] * area / (d1 + eps);
            double T2_flow = KtF[globalID2] * area / (d2 + eps);
            double trans_flow = (T1_flow * T2_flow) / (T1_flow + T2_flow + eps);

            // --- Heat (半传导率调和平均) ---
            double lam_eff_1 = PhiF[globalID1] * lam_fluid_ref + (1.0 - PhiF[globalID1]) * LamF[globalID1];
            double lam_eff_2 = PhiF[globalID2] * lam_fluid_ref + (1.0 - PhiF[globalID2]) * LamF[globalID2];

            double T1_heat = lam_eff_1 * area / (d1 + eps);
            double T2_heat = lam_eff_2 * area / (d2 + eps);
            double trans_heat = (T1_heat * T2_heat) / (T1_heat + T2_heat + eps);

            if (connectionIdx < T_ff_flow.size()) {
                T_ff_flow[connectionIdx] = trans_flow;
                T_ff_heat[connectionIdx] = trans_heat;
            }
            connectionIdx++;
        }
       
    }
    std::cout << "[Transmissibility] Calculated FF (F-F) for " << connectionIdx << " connections." << std::endl;
}