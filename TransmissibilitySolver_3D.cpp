#include "TransmissibilitySolver_3D.h"
#include "SolverContrlStrName_op.h" // 包含 PhysicalProperties_string 定义
#include "FVM_Ops.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <tuple>
#include <cstdint>
#include <string>


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

// ==================== 新增部分开始：3D Matrix 传导率计算 ====================
void TransmissibilitySolver_3D::Calculate_Transmissibility_Matrix(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
{
    // 1. 获取 3D 基岩网格数据
    const Mesh& mesh = meshMgr.mesh(); // [已修复] 调用正确的 mesh() 接口
    const auto& faces = mesh.getFaces();
    const auto& cells = mesh.getCells();
    const auto& cellId2Idx = mesh.getCellId2Index();

    // 2. 获取 3D 基岩物性场 (消除硬编码)
    PhysicalProperties_string_op::Rock rockStr; // 实例化名称结构体
    auto Kxx = fieldMgr.getMatrixScalar(rockStr.k_xx_tag);
    auto Kyy = fieldMgr.getMatrixScalar(rockStr.k_yy_tag);
    auto Kzz = fieldMgr.getMatrixScalar(rockStr.k_zz_tag);
    auto Lam_m = fieldMgr.getMatrixScalar(rockStr.lambda_tag);

    if (!Kxx || !Kyy || !Kzz) {
        std::cerr << "[Warning] Matrix Permeability (Kxx, Kyy, Kzz) not found! Skipping 3D Flow Transmissibility." << std::endl;
        return;
    }

    // 3. 分配面心场容器
    auto T_Flow = fieldMgr.getOrCreateMatrixFaceScalar("T_Matrix_Flow", 0.0);
    auto T_Heat = fieldMgr.getOrCreateMatrixFaceScalar("T_Matrix_Heat", 0.0);

    // 4. 面计算具有完美的独立性，安全开启 OpenMP 加速以应对百万级面计算
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(faces.size()); ++i)
    {
        const Face& face = faces[i];

        // 边界跳过机制
        if (face.isBoundary()) {
            continue;
        }

        size_t idxO = cellId2Idx.at(face.ownerCell);
        size_t idxN = cellId2Idx.at(face.neighborCell);

        double nx = face.normal.m_x;
        double ny = face.normal.m_y;
        double nz = face.normal.m_z;

        // 3D 渗透率对角张量基于面法向的标量投影
        double KO = (nx * nx * (*Kxx)[idxO]) + (ny * ny * (*Kyy)[idxO]) + (nz * nz * (*Kzz)[idxO]);
        double KN = (nx * nx * (*Kxx)[idxN]) + (ny * ny * (*Kyy)[idxN]) + (nz * nz * (*Kzz)[idxN]);

        Vector dVecO = face.midpoint - cells[idxO].center;
        Vector dVecN = cells[idxN].center - face.midpoint;

        // 3D 绝对法向投影距离计算
        double dO = std::max(std::abs(dVecO.m_x * nx + dVecO.m_y * ny + dVecO.m_z * nz), 1e-6);
        double dN = std::max(std::abs(dVecN.m_x * nx + dVecN.m_y * ny + dVecN.m_z * nz), 1e-6);

        // 3D 场景下，正交分解得到的 |E| 模长，本身即为真实的有效正交表面积
        double area = face.vectorE.Mag();

        // 写入结果 (数组索引 i 直接严格对齐 Face 的全局数组索引)
        (*T_Flow)[i] = FVM_Ops::Op_Math_Transmissibility(dO, KO, dN, KN, area);

        if (Lam_m) {
            double LO = (*Lam_m)[idxO];
            double LN = (*Lam_m)[idxN];
            (*T_Heat)[i] = FVM_Ops::Op_Math_Transmissibility(dO, LO, dN, LN, area);
        }
    }
}

// =========================================================================
// 静态传导率计算：FI (Fracture Internal in 3D)
// =========================================================================
void TransmissibilitySolver_3D::Calculate_Transmissibility_FractureInternal(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
{
    std::cout << "\n[Solver 3D] Calculating Fracture Internal (FI) Transmissibility..." << std::endl;

    PhysicalProperties_string_op::Fracture_string fracStr;
    PhysicalProperties_string_op::Water waterStr;

    const auto& frNet = meshMgr.fracture_network();
    const auto& globalEdges = frNet.getGlobalEdges();
    const auto& fracElements = frNet.getOrderedFractureElements();

    // [Fix-1] 使用真实 solver offset，而不是假定 nMat
    const int solverOffset = frNet.getSolverIndexOffset();

    // --- 获取场数据 ---
    auto p_Kt = fieldMgr.getFractureScalar(fracStr.k_t_tag);
    auto p_Wf = fieldMgr.getFractureScalar(fracStr.aperture_tag);
    auto p_LamF = fieldMgr.getFractureScalar(fracStr.lambda_tag);
    auto p_PhiF = fieldMgr.getFractureScalar(fracStr.phi_tag);
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Kt || !p_Wf) {
        std::cerr << "[Error] Critical fracture properties (Kt/Wf) missing for FI!" << std::endl;
        return;
    }

    const auto& Kt = p_Kt->data;
    const auto& Wf = p_Wf->data;
    const std::vector<double>* LamF = p_LamF ? &p_LamF->data : nullptr;
    const std::vector<double>* PhiF = p_PhiF ? &p_PhiF->data : nullptr;
    const std::vector<double>* LamFluid = p_LamFluid ? &p_LamFluid->data : nullptr;

    const size_t totalEdges = globalEdges.size();

    auto p_T_Flow = fieldMgr.getOrCreateFractureEdgeScalar("T_FI_Flow", 0.0);
    auto p_T_Heat = fieldMgr.getOrCreateFractureEdgeScalar("T_FI_Heat", 0.0);
    if (!p_T_Flow || !p_T_Heat) {
        std::cerr << "[Error] Failed to create/get FI output fields." << std::endl;
        return;
    }

    auto& T_FI_Flow = p_T_Flow->data;
    auto& T_FI_Heat = p_T_Heat->data;

    // [Fix-3] 每次调用先清零，避免 size 不变时残留旧值
    T_FI_Flow.assign(totalEdges, 0.0);
    T_FI_Heat.assign(totalEdges, 0.0);

    size_t validConnCount = 0;

    for (size_t i = 0; i < totalEdges; ++i)
    {
        const auto& edge = globalEdges[i];
        const int s1 = edge.ownerCell_solverIndex;
        const int s2 = edge.neighborCell_solverIndex;

        // [Fix-5] owner/neighbor 任一无效都跳过
        if (s1 < 0 || s2 < 0) continue;

        const int fLoc1 = s1 - solverOffset;
        const int fLoc2 = s2 - solverOffset;

        // [Fix-2] 增加 Wf.size() 与 fracElements.size() 的完整边界检查
        if (fLoc1 < 0 || fLoc2 < 0 ||
            fLoc1 >= static_cast<int>(Kt.size()) || fLoc2 >= static_cast<int>(Kt.size()) ||
            fLoc1 >= static_cast<int>(Wf.size()) || fLoc2 >= static_cast<int>(Wf.size()) ||
            fLoc1 >= static_cast<int>(fracElements.size()) || fLoc2 >= static_cast<int>(fracElements.size())) {
            continue;
        }

        const auto* pElem1 = fracElements[fLoc1];
        const auto* pElem2 = fracElements[fLoc2];
        if (!pElem1 || !pElem2) continue;

        // [Fix-4] 插值系数限幅 + 距离/面积最小值保护
        const double d_on = std::max(edge.ownerToNeighbor.Mag(), 1e-6);
        const double w = std::max(0.0, std::min(edge.f_linearInterpolationCoef, 1.0));
        const double d1 = std::max(d_on * w, 1e-6);
        const double d2 = std::max(d_on * (1.0 - w), 1e-6);
        const double area = std::max(edge.length, 1e-12);

        const double cond1 = Kt[fLoc1] * Wf[fLoc1];
        const double cond2 = Kt[fLoc2] * Wf[fLoc2];

        T_FI_Flow[i] = FVM_Ops::Op_Math_Transmissibility(d1, cond1, d2, cond2, area);

        if (LamF && PhiF && LamFluid &&
            fLoc1 < static_cast<int>(LamF->size()) && fLoc2 < static_cast<int>(LamF->size()) &&
            fLoc1 < static_cast<int>(PhiF->size()) && fLoc2 < static_cast<int>(PhiF->size()) &&
            fLoc1 < static_cast<int>(LamFluid->size()) && fLoc2 < static_cast<int>(LamFluid->size())) {

            const double lam_eff_1 = (*PhiF)[fLoc1] * (*LamFluid)[fLoc1] + (1.0 - (*PhiF)[fLoc1]) * (*LamF)[fLoc1];
            const double lam_eff_2 = (*PhiF)[fLoc2] * (*LamFluid)[fLoc2] + (1.0 - (*PhiF)[fLoc2]) * (*LamF)[fLoc2];
            const double h_cond1 = lam_eff_1 * Wf[fLoc1];
            const double h_cond2 = lam_eff_2 * Wf[fLoc2];

            T_FI_Heat[i] = FVM_Ops::Op_Math_Transmissibility(d1, h_cond1, d2, h_cond2, area);
        }

        validConnCount++;
    }

    std::cout << "[Solver 3D] FI Done (" << validConnCount
        << " internal connections on " << totalEdges << " total edges)." << std::endl;
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

    // 获取流体热导率 lambda_w (Assume initialized)
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
    const int nMat = static_cast<int>(meshMgr.mesh().getCells().size());

    const double MIN_VAL = 1e-20;

#pragma omp parallel for schedule(static)
    for (long long i = 0; i < numPairs; ++i)
    {
        const auto& pair = pairs[i];
        int mIdx = pair.matrixSolverIndex;
        int fIdx = pair.fracCellSolverIndex;

        int fLocIdx = fIdx - nMat;

        // 1. Matrix Directional Permeability
        double nx = pair.polygonNormal.m_x;
        double ny = pair.polygonNormal.m_y;
        double nz = pair.polygonNormal.m_z;
        double k_m_dir = (nx * nx * Kxx[mIdx]) + (ny * ny * Kyy[mIdx]) + (nz * nz * Kzz[mIdx]);

        // 2. Flow Transmissibility
        double dist = std::max(pair.distMatrixToFracPlane, 1e-6);
        // 调用底层算子，裂缝侧物理距离为 Wf/2.0
        T_Flow[i] = FVM_Ops::Op_Math_Transmissibility(dist, k_m_dir, Wf[fLocIdx] / 2.0, Kf[fLocIdx], pair.intersectionArea);

        // 3. Heat Transmissibility [Modified]
        // 使用有效热导率：lam_eff = phi * lam_fluid + (1-phi) * lam_solid
        if (Lam_m && Lam_f && Phi_f && LamFluid) {
            double lam_m_val = (*Lam_m)[mIdx];
            double phi_val = (*Phi_f)[fLocIdx];
            double lam_fluid_val = (*LamFluid)[fLocIdx];
            double lam_solid_val = (*Lam_f)[fLocIdx];

            // 计算裂缝侧有效热导率
            double lam_f_eff = phi_val * lam_fluid_val + (1.0 - phi_val) * lam_solid_val;
            T_Heat[i] = FVM_Ops::Op_Math_Transmissibility(dist, lam_m_val, Wf[fLocIdx] / 2.0, lam_f_eff, pair.intersectionArea);
        }
        else if (Lam_m && Lam_f) // Fallback for pure solid conduction (Compatibility)
        {
            double lam_m_val = (*Lam_m)[mIdx];
            double lam_f_val = (*Lam_f)[fLocIdx];

            T_Heat[i] = FVM_Ops::Op_Math_Transmissibility(dist, lam_m_val, Wf[fLocIdx] / 2.0, lam_f_val, pair.intersectionArea);
        }

    }
    std::cout << "[Solver] NNC Done (" << numPairs << " pairs)." << std::endl;
}


// =========================================================================
// 静态传导率计算：FF (Fracture - Fracture 3D Star-Delta Model)
// =========================================================================
/**
 * @brief 计算 3D 裂缝-裂缝 (FF) 星角变换传导率 (Industrial-Grade)
 * @details
 * 1. 端点量化无向线段键聚类 (Quantized Undirected Segment Key)，确保交线簇识别的绝对鲁棒。
 * 2. 几何要素同步锁定，保留交线簇中最大线段的严格起止坐标。
 * 3. 物理场解耦过滤：Flow 和 Heat 各自独立判定有效性，Heat 缺失自动回退置 0，绝不影响 Flow 拓扑。
 * 4. 强确定性排序 (Deterministic Ordering)：消除哈希迭代的随机性，保证 CSV 报表和矩阵装配绝对可复现。
 * @param meshMgr 3D 网格管理器
 * @param fieldMgr 3D 场管理器
 */
void TransmissibilitySolver_3D::Calculate_Transmissibility_FF(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
{
    std::cout << "\n[Solver 3D] Calculating FF Transmissibility (Deterministic Star-Delta)..." << std::endl;

    PhysicalProperties_string_op::Fracture_string fracStr;
    PhysicalProperties_string_op::Water waterStr;

    const auto& frNet = meshMgr.fracture_network();
    const auto& ffIntersections = frNet.ffIntersections;
    const auto& fracElements = frNet.getOrderedFractureElements();

    // 获取底层统一设定的索引偏移量
    const int offset = frNet.getSolverIndexOffset();

    // --- 获取并校验场数据 ---
    auto p_Kt = fieldMgr.getFractureScalar(fracStr.k_t_tag);
    auto p_Wf = fieldMgr.getFractureScalar(fracStr.aperture_tag);
    auto p_LamF = fieldMgr.getFractureScalar(fracStr.lambda_tag);
    auto p_PhiF = fieldMgr.getFractureScalar(fracStr.phi_tag);
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Kt || !p_Wf) {
        std::cerr << "[Error] Critical fracture properties (Kt/Wf) missing for FF!" << std::endl;
        return;
    }

    const auto& Kt = p_Kt->data;
    const auto& Wf = p_Wf->data;
    const auto& vLamF = (p_LamF) ? p_LamF->data : std::vector<double>();
    const auto& vPhiF = (p_PhiF) ? p_PhiF->data : std::vector<double>();
    const auto& vLamFluid = (p_LamFluid) ? p_LamFluid->data : std::vector<double>();

    bool hasGlobalHeat = (!vLamF.empty() && !vPhiF.empty() && !vLamFluid.empty());

    // =====================================================================
    // 步骤 1：基于端点量化的无向线段键聚类
    // =====================================================================
    struct JunctionCluster {
        Vector start;
        Vector end;
        double length = -1.0;
        std::vector<int> solverIndices;
    };
    std::unordered_map<std::string, JunctionCluster> clusterMap;
    const double TOL = 1e-4;

    auto quantize = [TOL](const Vector& v) -> std::tuple<long long, long long, long long> {
        return std::make_tuple(
            static_cast<long long>(std::floor(v.m_x / TOL + 0.5)),
            static_cast<long long>(std::floor(v.m_y / TOL + 0.5)),
            static_cast<long long>(std::floor(v.m_z / TOL + 0.5))
        );
        };

    auto getUndirectedKey = [&](const Vector& p1, const Vector& p2) -> std::string {
        auto q1 = quantize(p1), q2 = quantize(p2);
        if (q1 > q2) std::swap(q1, q2);
        return std::to_string(std::get<0>(q1)) + "_" + std::to_string(std::get<1>(q1)) + "_" + std::to_string(std::get<2>(q1)) + "-" +
            std::to_string(std::get<0>(q2)) + "_" + std::to_string(std::get<1>(q2)) + "_" + std::to_string(std::get<2>(q2));
        };

    for (const auto& inter : ffIntersections) {
        for (const auto& seg : inter.segments) {
            std::string key = getUndirectedKey(seg.start, seg.end);
            auto& cluster = clusterMap[key];

            // [修复 1] 几何一致性：只在捕捉到更长的代表性线段时，同步更新长度与起止坐标
            if (seg.length > cluster.length) {
                cluster.start = seg.start;
                cluster.end = seg.end;
                cluster.length = seg.length;
            }

            if (seg.solverIndex_1 >= 0) {
                if (std::find(cluster.solverIndices.begin(), cluster.solverIndices.end(), seg.solverIndex_1) == cluster.solverIndices.end())
                    cluster.solverIndices.push_back(seg.solverIndex_1);
            }
            if (seg.solverIndex_2 >= 0) {
                if (std::find(cluster.solverIndices.begin(), cluster.solverIndices.end(), seg.solverIndex_2) == cluster.solverIndices.end())
                    cluster.solverIndices.push_back(seg.solverIndex_2);
            }
        }
    }

    // =====================================================================
    // 步骤 2：确定性排序与有效性过滤 (Deterministic Ordering & Validation)
    // =====================================================================
    // [修复 3] 提取所有 Key 并字典序排序，保证后续计算和配对的绝对确定性可复现
    std::vector<std::string> sortedKeys;
    sortedKeys.reserve(clusterMap.size());
    for (const auto& kv : clusterMap) {
        sortedKeys.push_back(kv.first);
    }
    std::sort(sortedKeys.begin(), sortedKeys.end());

    std::vector<JunctionCluster> validClusters;
    size_t totalFFPairs = 0;

    for (const auto& key : sortedKeys) {
        auto& cluster = clusterMap[key];
        std::vector<int> validIndices;

        for (int sIdx : cluster.solverIndices) {
            int fLoc = sIdx - offset;

            // 基础渗流越界检查 (Flow 属性是拓扑存活的唯一判据)
            if (fLoc < 0 || fLoc >= (int)Kt.size() || fLoc >= (int)Wf.size() || fLoc >= (int)fracElements.size()) continue;

            const auto* pElem = fracElements[fLoc];
            if (!pElem) continue;

            // [修复 2] 将 Heat 的越界检查从外部循环剥离，仅作为局部的可用性 flag 参与运算
            validIndices.push_back(sIdx);
        }
        std::sort(validIndices.begin(), validIndices.end()); // <-- 新增
        cluster.solverIndices = validIndices;

        size_t n = validIndices.size();
        if (n >= 2) {
            totalFFPairs += n * (n - 1) / 2;
            validClusters.push_back(cluster); // 压入的顺序自然就是确定性的
        }
    }

    // --- 分配内存 ---
    auto p_T_Flow = fieldMgr.createFFScalar("T_FF_Flow", 0.0);
    auto p_T_Heat = fieldMgr.createFFScalar("T_FF_Heat", 0.0);
    auto& T_FF_Flow = p_T_Flow->data;
    auto& T_FF_Heat = p_T_Heat->data;

    if (T_FF_Flow.size() != totalFFPairs) {
        T_FF_Flow.assign(totalFFPairs, 0.0);
        T_FF_Heat.assign(totalFFPairs, 0.0);
    }

    // =====================================================================
    // 步骤 3：遍历合法枢纽，执行 Star-Delta 星角变换
    // =====================================================================
    size_t ffIdx = 0;
    for (const auto& cluster : validClusters) {
        size_t nElems = cluster.solverIndices.size();

        std::vector<double> half_T_Flow(nElems, 0.0);
        std::vector<double> half_T_Heat(nElems, 0.0);
        double sum_T_Flow = 0.0;
        double sum_T_Heat = 0.0;

        for (size_t i = 0; i < nElems; ++i) {
            int fLoc = cluster.solverIndices[i] - offset;
            const auto* pElem = fracElements[fLoc];

            double d = std::max(Geometry_3D::PointToSegmentDistance(pElem->centroid, cluster.start, cluster.end), 1e-6);

            // Flow
            double cond_flow = Kt[fLoc] * Wf[fLoc];
            double t_f = (cond_flow * cluster.length) / d;
            half_T_Flow[i] = t_f;
            sum_T_Flow += t_f;

            // [修复 2 延续] 独立的 Heat 判定：当前支路是否同时拥有完备的热学数据？
            bool branchHasHeat = hasGlobalHeat && (fLoc < (int)vLamF.size()) && (fLoc < (int)vPhiF.size()) && (fLoc < (int)vLamFluid.size());
            if (branchHasHeat) {
                double lam_eff = vPhiF[fLoc] * vLamFluid[fLoc] + (1.0 - vPhiF[fLoc]) * vLamF[fLoc];
                double cond_heat = lam_eff * Wf[fLoc];
                double t_h = (cond_heat * cluster.length) / d;
                half_T_Heat[i] = t_h;
                sum_T_Heat += t_h;
            }
            else {
                half_T_Heat[i] = 0.0; // 缺失热学数据，仅流体连通，热连通阻断
            }
        }

        // Star-Delta 展开
        fieldMgr.ff_topology.clear();
        fieldMgr.ff_topology.reserve(totalFFPairs);
        for (size_t i = 0; i < nElems; ++i) {
            for (size_t j = i + 1; j < nElems; ++j) {
                T_FF_Flow[ffIdx] = (sum_T_Flow > 1e-25) ? ((half_T_Flow[i] * half_T_Flow[j]) / sum_T_Flow) : 0.0;
                T_FF_Heat[ffIdx] = (sum_T_Heat > 1e-25) ? ((half_T_Heat[i] * half_T_Heat[j]) / sum_T_Heat) : 0.0;

                fieldMgr.ff_topology.emplace_back(cluster.solverIndices[i], cluster.solverIndices[j]);

                ffIdx++;
            }
        }
    }

    std::cout << "[Solver 3D] FF Done (" << totalFFPairs << " Deterministic Star-Delta pairs over "
        << validClusters.size() << " valid junctions)." << std::endl;
}