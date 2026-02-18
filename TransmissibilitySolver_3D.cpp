#include "TransmissibilitySolver_3D.h"
#include "SolverContrlStrName_op.h" // 包含 PhysicalProperties_string 定义

#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_map>

using namespace PhysicalProperties_string_op;

// =========================================================
// 内部辅助函数
// =========================================================
namespace Geometry_3D {

    // 手动实现点积 (避免依赖 Vector 类是否重载了 *)
    inline double Dot(const Vector& a, const Vector& b) {
        return a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z;
    }

    // 手动实现模长平方 (解决 "Variable3D 没有成员 MagSqr" 报错)
    inline double MagSqr(const Vector& v) {
        return v.m_x * v.m_x + v.m_y * v.m_y + v.m_z * v.m_z;
    }

    // 手动实现模长
    inline double Mag(const Vector& v) {
        return std::sqrt(MagSqr(v));
    }

    /**
     * @brief 计算点 P 到线段 AB 的最短物理距离 (Robust)
     * @details 自动处理投影落在线段外的情况，避免传导率奇异
     */
    inline double PointToSegmentDistance(const Vector& P, const Vector& A, const Vector& B) {
        Vector AB = B - A;
        Vector AP = P - A;

        double lenSq = MagSqr(AB); // [Fix] 使用本地函数替代 AB.MagSqr()

        // 退化情况处理：线段长度极小，退化为点
        if (lenSq < 1e-12) return Mag(AP);

        // 计算投影系数 t = (AP . AB) / |AB|^2
        double t = Dot(AP, AB) / lenSq;

        // 钳位操作 (Clamping)：强制 t 在 [0, 1] 范围内
        if (t < 0.0) t = 0.0;
        else if (t > 1.0) t = 1.0;

        // 计算最近点坐标
        Vector ClosestPoint = A + AB * t;

        // 返回点 P 到最近点的距离
        return Mag(P - ClosestPoint);
    }
}

// =========================================================================
// 静态传导率计算：NNC (Matrix - Fracture)
// =========================================================================
void TransmissibilitySolver_3D::Calculate_Transmissibility_NNC(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
{
    Rock rock_str;
    Fracture_string frac_str;
    Water waterStr;

    // --- 获取场数据 (Direct Pointers) ---
    auto p_Kxx = fieldMgr.getMatrixScalar(rock_str.k_xx_tag);
    auto p_Kyy = fieldMgr.getMatrixScalar(rock_str.k_yy_tag);
    auto p_Kzz = fieldMgr.getMatrixScalar(rock_str.k_zz_tag);
    auto p_Lam_m = fieldMgr.getMatrixScalar(rock_str.lambda_tag);

    if (!p_Kxx) { std::cerr << "[Error] Matrix Permeability missing!" << std::endl; return; }

    const std::vector<double>& Kxx = p_Kxx->data;
    const std::vector<double>& Kyy = (p_Kyy) ? p_Kyy->data : Kxx;
    const std::vector<double>& Kzz = (p_Kzz) ? p_Kzz->data : Kxx;
    const std::vector<double>* Lam_m = (p_Lam_m) ? &(p_Lam_m->data) : nullptr;

    auto p_Kf = fieldMgr.getFractureScalar(frac_str.k_n_tag);
    if (!p_Kf) p_Kf = fieldMgr.getFractureScalar(frac_str.k_t_tag);
    auto p_Wf = fieldMgr.getFractureScalar(frac_str.aperture_tag);
    auto p_Lam_f = fieldMgr.getFractureScalar(frac_str.lambda_tag); // Solid thermal cond
    auto p_Phi_f = fieldMgr.getFractureScalar(frac_str.phi_tag);    // [New] Porosity

    // [New] 获取流体热导率 lambda_w (Assume initialized)
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Kf || !p_Wf) { std::cerr << "[Error] Frac properties missing!" << std::endl; return; }

    const std::vector<double>& Kf = p_Kf->data;
    const std::vector<double>& Wf = p_Wf->data;
    const std::vector<double>* Lam_f = (p_Lam_f) ? &(p_Lam_f->data) : nullptr;
    const std::vector<double>* Phi_f = (p_Phi_f) ? &(p_Phi_f->data) : nullptr;    // [New]
    const std::vector<double>* LamFluid = (p_LamFluid) ? &(p_LamFluid->data) : nullptr; // [New]

    // --- 输出场 ---
    const std::string tag_T_NNC_Flow = "T_NNC_Flow";
    const std::string tag_T_NNC_Heat = "T_NNC_Heat";

    auto& T_Flow = fieldMgr.createNNCScalar(tag_T_NNC_Flow, 0.0)->data;
    auto& T_Heat = fieldMgr.createNNCScalar(tag_T_NNC_Heat, 0.0)->data;

    const auto& pairs = meshMgr.getInteractionPairs();
    const size_t numPairs = pairs.size();
    if (T_Flow.size() != numPairs) { T_Flow.resize(numPairs); T_Heat.resize(numPairs); }

    const double MIN_VAL = 1e-20;

#pragma omp parallel for schedule(static)
    for (long long i = 0; i < numPairs; ++i)
    {
        const auto& pair = pairs[i];
        int mIdx = pair.matrixSolverIndex;
        int fIdx = pair.fracCellSolverIndex;

        // 1. Matrix Directional Permeability
        double nx = pair.polygonNormal.m_x;
        double ny = pair.polygonNormal.m_y;
        double nz = pair.polygonNormal.m_z;
        double k_m_dir = (nx * nx * Kxx[mIdx]) + (ny * ny * Kyy[mIdx]) + (nz * nz * Kzz[mIdx]);

        // 2. Flow Transmissibility
        double dist = std::max(pair.distMatrixToFracPlane, 1e-6);
        double term_m = dist / (k_m_dir + MIN_VAL);
        double term_f = Wf[fIdx] / (2.0 * Kf[fIdx] + MIN_VAL);

        T_Flow[i] = pair.intersectionArea / (term_m + term_f + MIN_VAL);

        // 3. Heat Transmissibility [Modified]
        // 使用有效热导率：lam_eff = phi * lam_fluid + (1-phi) * lam_solid
        if (Lam_m && Lam_f && Phi_f && LamFluid) {
            double lam_m_val = (*Lam_m)[mIdx];

            double phi_val = (*Phi_f)[fIdx];
            double lam_fluid_val = (*LamFluid)[fIdx];
            double lam_solid_val = (*Lam_f)[fIdx];

            // 计算裂缝侧有效热导率
            double lam_f_eff = phi_val * lam_fluid_val + (1.0 - phi_val) * lam_solid_val;

            double term_m_h = dist / (lam_m_val + MIN_VAL);
            double term_f_h = Wf[fIdx] / (2.0 * lam_f_eff + MIN_VAL);

            T_Heat[i] = pair.intersectionArea / (term_m_h + term_f_h + MIN_VAL);
        }
        else if (Lam_m && Lam_f) // Fallback for pure solid conduction (Compatibility)
        {
            double lam_m_val = (*Lam_m)[mIdx];
            double lam_f_val = (*Lam_f)[fIdx];

            double term_m_h = dist / (lam_m_val + MIN_VAL);
            double term_f_h = Wf[fIdx] / (2.0 * lam_f_val + MIN_VAL);

            T_Heat[i] = pair.intersectionArea / (term_m_h + term_f_h + MIN_VAL);
        }

    }
    std::cout << "[Solver] NNC Done (" << numPairs << " pairs)." << std::endl;
}


// =========================================================================
// 静态传导率计算：FF (Fracture - Fracture)
// =========================================================================
void TransmissibilitySolver_3D::Calculate_Transmissibility_FF(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
{
    std::cout << "\n[Solver] Calculating FF Transmissibility (Series Model)..." << std::endl;

    Fracture_string fracStr;
    Water waterStr;
    const double eps = 1e-15;
    const double lam_fluid_ref = 0.6;

    const auto& frNet = meshMgr.fracture_network();
    const auto& ffIntersections = frNet.ffIntersections;
    const auto& fracturesVec = frNet.getFractures();

    // --- 获取场数据 ---
    auto& Kf = fieldMgr.getFractureScalar(fracStr.k_t_tag)->data;
    auto& LamF = fieldMgr.getFractureScalar(fracStr.lambda_tag)->data;
    auto& PhiF = fieldMgr.getFractureScalar(fracStr.phi_tag)->data;
    auto& Wf = fieldMgr.getFractureScalar(fracStr.aperture_tag)->data;

    // 2. [New] 获取流体属性 (从 Water 结构体中获取 "lambda_w")
    auto& LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag)->data;

    const std::string tag_T_NNC_Flow = "T_FF_Flow";                                // 渗流传导率标签
    const std::string tag_T_NNC_Heat = "T_FF_Heat";                                // 传热传导率标签

    auto& T_Flow = fieldMgr.createNNCScalar(tag_T_NNC_Flow, 0.0)->data;
    auto& T_Heat = fieldMgr.createNNCScalar(tag_T_NNC_Heat, 0.0)->data;

    // 构建 ID 映射
    std::unordered_map<int, size_t> idToIndexMap;
    idToIndexMap.reserve(fracturesVec.size());
    for (size_t k = 0; k < fracturesVec.size(); ++k) idToIndexMap[fracturesVec[k].id] = k;

    // 预计算 offset 以便并行
    std::vector<size_t> offsets(ffIntersections.size() + 1, 0);
    for (size_t i = 0; i < ffIntersections.size(); ++i) offsets[i + 1] = offsets[i] + ffIntersections[i].segments.size();

    if (T_Flow.size() != offsets.back()) {
        T_Flow.resize(offsets.back()); T_Heat.resize(offsets.back());
    }

#pragma omp parallel for schedule(dynamic)
    for (long long i = 0; i < ffIntersections.size(); ++i)
    {
        const auto& interaction = ffIntersections[i];
        size_t base_idx = offsets[i];

        auto it1 = idToIndexMap.find(interaction.fracID_1);
        auto it2 = idToIndexMap.find(interaction.fracID_2);
        if (it1 == idToIndexMap.end() || it2 == idToIndexMap.end()) continue;

        const auto& f1 = fracturesVec[it1->second];
        const auto& f2 = fracturesVec[it2->second];

        for (size_t j = 0; j < interaction.segments.size(); ++j)
        {
            const auto& seg = interaction.segments[j];
            int s1 = seg.solverIndex_1;
            int s2 = seg.solverIndex_2;

            if (s1 < 0 || s2 < 0 || s1 >= static_cast<int>(Kf.size()) || s2 >= static_cast<int>(Kf.size())) continue;

            // Geometry
            int lid1 = seg.cellID_1;
            int lid2 = seg.cellID_2;
            if (lid1 >= static_cast<int>(f1.fracCells.size()) || lid2 >= static_cast<int>(f2.fracCells.size())) continue;

            // 使用 Geometry_3D 辅助函数计算点到线段距离
            double d1 = std::max(Geometry_3D::PointToSegmentDistance(f1.fracCells[lid1].centroid, seg.start, seg.end), 1e-6);
            double d2 = std::max(Geometry_3D::PointToSegmentDistance(f2.fracCells[lid2].centroid, seg.start, seg.end), 1e-6);

            // --- Flow Transmissibility (Series Resistance) ---
            // T = L / ( d1/(w1*K1) + d2/(w2*K2) )
            // Resistance 1 per length: R1 = d1 / (Wf[s1] * Kf[s1])
            double cond1 = Wf[s1] * Kf[s1];
            double cond2 = Wf[s2] * Kf[s2];

            double res_term = (d1 / (cond1 + eps)) + (d2 / (cond2 + eps));
            T_Flow[base_idx + j] = seg.length / (res_term + eps);

            // --- Heat Transmissibility ---
            double lam_eff_1 = PhiF[s1] * LamFluid[s1] + (1.0 - PhiF[s1]) * LamF[s1];
            double lam_eff_2 = PhiF[s2] * LamFluid[s2] + (1.0 - PhiF[s2]) * LamF[s2];

            double cond1_h = Wf[s1] * lam_eff_1;
            double cond2_h = Wf[s2] * lam_eff_2;

            double res_term_h = (d1 / (cond1_h + eps)) + (d2 / (cond2_h + eps));
            T_Heat[base_idx + j] = seg.length / (res_term_h + eps);
        }
    }
    std::cout << "[Solver] FF Done (" << offsets.back() << " connections)." << std::endl;
}


////    获取场数据 (使用 const 版本的 meshMgr 接口)
   ////auto& Kxx = fieldMgr.getOrCreateMatrixScalar(rockStr.k_xx_tag, 1e-15)->data;    // 基岩渗透率场
   ////auto& Kyy = fieldMgr.getOrCreateMatrixScalar(rockStr.k_yy_tag, 1e-15)->data;    // 基岩渗透率场
   ////auto& Kzz = fieldMgr.getOrCreateMatrixScalar(rockStr.k_zz_tag, 1e-15)->data;    // 基岩渗透率场
   ////auto& LamM = fieldMgr.getOrCreateMatrixScalar(rockStr.lambda_tag, 2.0)->data;   // 基岩导热系数场
   ////auto& PhiM = fieldMgr.getOrCreateMatrixScalar(rockStr.phi_tag, 0.1)->data;      // 基岩孔隙度场
   ////auto& KnF = fieldMgr.getOrCreateFractureScalar(fracStr.k_n, 1e-12)->data;       // 裂缝法向渗透率场
   ////auto& LamF = fieldMgr.getOrCreateFractureScalar(fracStr.lambda_tag, 2.0)->data;	// 裂缝导热系数场
   ////auto& PhiF = fieldMgr.getOrCreateFractureScalar(fracStr.phi_tag, 1.0)->data;    // 裂缝孔隙度场
   ////const std::string tag_T_NNC_Flow = "T_NNC_Flow";                                // 渗流传导率标签
   ////const std::string tag_T_NNC_Heat = "T_NNC_Heat";                                // 传热传导率标签
////   auto& T_flow = fieldMgr.createNNCScalar(tag_T_NNC_Flow, 0.0)->data;
////   auto& T_heat = fieldMgr.createNNCScalar(tag_T_NNC_Heat, 0.0)->data;
////   const auto& pairs = meshMgr.getInteractionPairs();
////   const auto& frNet = meshMgr.fracture_network(); 
////   auto fracMap = frNet.buildFractureIDMap();
////   std::cout << "[Transmissibility] Calculating NNC (M-F) for " << pairs.size() << " pairs..." << std::endl;
////   for (size_t i = 0; i < pairs.size(); ++i)
////   {
////       const auto& pair = pairs[i];
   ////	int mGlobalID = pair.matrixCellGlobalID;                                    // 基岩单元全局ID
   ////	int mIdx = meshMgr.mesh().getCellIndex(mGlobalID);                          // 基岩单元索引
   ////	int fElemID = pair.fracElementGlobalID;     			                    // 裂缝单元全局ID      
////       double A = pair.intersectionArea;                                           // 该交互多边形的面积          
   ////	double d = pair.distMatrixToFracPlane;									    // 基岩单元中心到裂缝平面的距离
   ////	Vector n = pair.polygonNormal;              			                    // 多边形法向量（已归一化）
////       double w = 0.0;
////       auto it = fracMap.find(pair.fracMacroID);
   ////	if (it != fracMap.end()) w = it->second->aperture;                          // 裂缝开度
////       else continue;
////        --- Flow ---
////       double Km_n = Kxx[mIdx] * n.m_x * n.m_x +
////           Kyy[mIdx] * n.m_y * n.m_y +
////           Kzz[mIdx] * n.m_z * n.m_z;
////       double Kf_n = KnF[fElemID];
////        保持 4.0 系数 经典文献（如 Karimi-Fard et al. 2004, Moinfar 2013）推导出的形状因子 (Shape Factor) 对应的分母项确实是 $4 K_f$
////       double term_flow_m = d / (Km_n + eps);
////       double term_flow_f = w / (4.0 * Kf_n + eps);
////       T_flow[i] = A / (term_flow_m + term_flow_f + eps);
////        --- Heat ---
////       double lam_eff_m = PhiM[mIdx] * lam_fluid_ref + (1.0 - PhiM[mIdx]) * LamM[mIdx];
////       double lam_eff_f = PhiF[fElemID] * lam_fluid_ref + (1.0 - PhiF[fElemID]) * LamF[fElemID];
////       double term_heat_m = d / (lam_eff_m + eps);
////       double term_heat_f = w / (4.0 * lam_eff_f + eps);
////       T_heat[i] = A / (term_heat_m + term_heat_f + eps);
  // }

//void TransmissibilitySolver_3D::Calculate_Transmissibility_FF(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
//{
//    Fracture_string fracStr; // [修改]
//    const double eps = 1e-30;
//	const double lam_fluid_ref = 0.6; // [Static Ref]参考流体导热系数 (W/mK)
//
//    const auto& frNet = meshMgr.fracture_network();
//    const auto& ffIntersections = frNet.ffIntersections;
//    const auto& indexer = frNet.fracElemIndex;
//    const auto& fracturesVec = frNet.getFractures();
//
//    auto& KtF = fieldMgr.getOrCreateFractureScalar(fracStr.k_t, 1e-12)->data;
//    auto& LamF = fieldMgr.getOrCreateFractureScalar(fracStr.lambda_tag, 2.0)->data;
//    auto& PhiF = fieldMgr.getOrCreateFractureScalar(fracStr.phi_tag, 1.0)->data;
//
//    const std::string tag_T_FF_Flow = "T_FF_Flow";
//    const std::string tag_T_FF_Heat = "T_FF_Heat";
//
//    auto& T_ff_flow = fieldMgr.createFFScalar(tag_T_FF_Flow, 0.0)->data;
//    auto& T_ff_heat = fieldMgr.createFFScalar(tag_T_FF_Heat, 0.0)->data;
//
//    std::cout << "[Transmissibility] Calculating FF (F-F) connections..." << std::endl;
//
//	// 构建裂缝ID到索引的映射，便于快速查找
//    std::unordered_map<int, size_t> idToIndexMap;
//    idToIndexMap.reserve(fracturesVec.size());
//    for (size_t k = 0; k < fracturesVec.size(); ++k) {
//        idToIndexMap[fracturesVec[k].id] = k;
//    }
//
//    size_t connectionIdx = 0;
//
//    for (const auto& interaction : ffIntersections)
//    {
//        int id1 = interaction.fracID_1;
//        int id2 = interaction.fracID_2;
//
//		// 查找裂缝索引
//        auto it1 = idToIndexMap.find(id1);
//        auto it2 = idToIndexMap.find(id2);
//
//        if (it1 == idToIndexMap.end() || it2 == idToIndexMap.end()) {
//            connectionIdx += interaction.segments.size();
//            continue;
//        }
//
//        size_t idx1 = it1->second;
//        size_t idx2 = it2->second;
//        const auto& f1 = fracturesVec[idx1];
//        const auto& f2 = fracturesVec[idx2];
//
//        // 调和平均开度 (Harmonic Aperture)
//        double w_harm = 2.0 * f1.aperture * f2.aperture / (f1.aperture + f2.aperture + eps);
//
//        for (const auto& seg : interaction.segments)
//        {
//            int localID1 = seg.cellID_1;
//            int localID2 = seg.cellID_2;
//
//            size_t globalID1 = indexer.offset[idx1] + localID1;
//            size_t globalID2 = indexer.offset[idx2] + localID2;
//
//            if (localID1 >= f1.fracCells.size() || localID2 >= f2.fracCells.size()) {
//                connectionIdx++; continue;
//            }
//
//            Vector c1 = f1.fracCells[localID1].centroid;
//            Vector c2 = f2.fracCells[localID2].centroid;
//
//            double d1 = _pointToLineDistance(c1, seg.start, seg.end, seg.length);
//            double d2 = _pointToLineDistance(c2, seg.start, seg.end, seg.length);
//
//            if (d1 < 1e-6) d1 = 1e-6;
//            if (d2 < 1e-6) d2 = 1e-6;
//
//            double area = seg.length * w_harm;
//
//            // --- Flow (半传导率调和平均) ---
//            double T1_flow = KtF[globalID1] * area / (d1 + eps);
//            double T2_flow = KtF[globalID2] * area / (d2 + eps);
//            double trans_flow = (T1_flow * T2_flow) / (T1_flow + T2_flow + eps);
//
//            // --- Heat (半传导率调和平均) ---
//            double lam_eff_1 = PhiF[globalID1] * lam_fluid_ref + (1.0 - PhiF[globalID1]) * LamF[globalID1];
//            double lam_eff_2 = PhiF[globalID2] * lam_fluid_ref + (1.0 - PhiF[globalID2]) * LamF[globalID2];
//
//            double T1_heat = lam_eff_1 * area / (d1 + eps);
//            double T2_heat = lam_eff_2 * area / (d2 + eps);
//            double trans_heat = (T1_heat * T2_heat) / (T1_heat + T2_heat + eps);
//
//            if (connectionIdx < T_ff_flow.size()) {
//                T_ff_flow[connectionIdx] = trans_flow;
//                T_ff_heat[connectionIdx] = trans_heat;
//            }
//            connectionIdx++;
//        }
//       
//    }
//    std::cout << "[Transmissibility] Calculated FF (F-F) for " << connectionIdx << " connections." << std::endl;
//}