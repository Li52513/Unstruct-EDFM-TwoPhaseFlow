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
#include <stdexcept>


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
     * @brief 计算�?P 到线�?AB 的最短物理距�?(Robust)
     * @details 自动处理投影落在线段外的情况，避免传导率奇异
     */
    inline double PointToSegmentDistance(const Vector& P, const Vector& A, const Vector& B) {
        Vector AB = B - A;
        Vector AP = P - A;

        double lenSq = MagSqr(AB); // [Fix] 使用本地函数替代 AB.MagSqr()

        // 退化情况处理：线段长度极小，退化为�?
        if (lenSq < 1e-12) return Mag(AP);

        // 计算投影系数 t = (AP . AB) / |AB|^2
        double t = Dot(AP, AB) / lenSq;

        // 钳位操作 (Clamping)：强�?t �?[0, 1] 范围�?
        if (t < 0.0) t = 0.0;
        else if (t > 1.0) t = 1.0;

        // 计算最近点坐标
        Vector ClosestPoint = A + AB * t;

        // 返回�?P 到最近点的距�?
        return Mag(P - ClosestPoint);
    }
}

// ==================== 新增部分开始：3D Matrix 传导率计�?====================
void TransmissibilitySolver_3D::Calculate_Transmissibility_Matrix(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
{
    // 1. 获取 3D 基岩网格数据
    const Mesh& mesh = meshMgr.mesh(); // [已修复] 调用正确�?mesh() 接口
    const auto& faces = mesh.getFaces();
    const auto& cells = mesh.getCells();
    const auto& cellId2Idx = mesh.getCellId2Index();

    // 2. 获取 3D 基岩物性场 (消除硬编�?
    PhysicalProperties_string_op::Rock rockStr; // 实例化名称结构体
    PhysicalProperties_string_op::TransmissibilityFields tTags;
    auto Kxx = fieldMgr.getMatrixScalar(rockStr.k_xx_tag);
    auto Kyy = fieldMgr.getMatrixScalar(rockStr.k_yy_tag);
    auto Kzz = fieldMgr.getMatrixScalar(rockStr.k_zz_tag);
    auto Lam_m = fieldMgr.getMatrixScalar(rockStr.lambda_tag);

    if (!Kxx || !Kyy || !Kzz) {
        std::cerr << "[Warning] Matrix Permeability (Kxx, Kyy, Kzz) not found! Skipping 3D Flow Transmissibility." << std::endl;
        return;
    }

    // 3. 分配面心场容�?
    auto T_Flow = fieldMgr.getOrCreateMatrixFaceScalar(tTags.matrix_flow, 0.0);
    auto T_Heat = fieldMgr.getOrCreateMatrixFaceScalar(tTags.matrix_heat, 0.0);

    // 4. 面计算具有完美的独立性，安全开�?OpenMP 加速以应对百万级面计算
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

        // 3D 场景下，正交分解得到�?|E| 模长，本身即为真实的有效正交表面�?
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
    PhysicalProperties_string_op::TransmissibilityFields tTags;

    const auto& frNet = meshMgr.fracture_network();
    const auto& globalEdges = frNet.getGlobalEdges();
    const auto& fracElements = frNet.getOrderedFractureElements();

    // [Fix-1] 使用真实 solver offset，而不是假�?nMat
    const int solverOffset = frNet.getSolverIndexOffset();

    // --- 获取场数�?---
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

    auto p_T_Flow = fieldMgr.getOrCreateFractureEdgeScalar(tTags.fi_flow, 0.0);
    auto p_T_Heat = fieldMgr.getOrCreateFractureEdgeScalar(tTags.fi_heat, 0.0);
    if (!p_T_Flow || !p_T_Heat) {
        std::cerr << "[Error] Failed to create/get FI output fields." << std::endl;
        return;
    }

    auto& T_FI_Flow = p_T_Flow->data;
    auto& T_FI_Heat = p_T_Heat->data;

    // [Fix-3] 每次调用先清零，避免 size 不变时残留旧�?
    T_FI_Flow.assign(totalEdges, 0.0);
    T_FI_Heat.assign(totalEdges, 0.0);

    size_t validConnCount = 0;

    for (size_t i = 0; i < totalEdges; ++i)
    {
        const auto& edge = globalEdges[i];
        const int s1 = edge.ownerCell_solverIndex;
        const int s2 = edge.neighborCell_solverIndex;

        // [Fix-5] owner/neighbor 任一无效都跳�?
        if (s1 < 0 || s2 < 0) continue;

        const int fLoc1 = s1 - solverOffset;
        const int fLoc2 = s2 - solverOffset;

        // [Fix-2] 增加 Wf.size() �?fracElements.size() 的完整边界检�?
        if (fLoc1 < 0 || fLoc2 < 0 ||
            fLoc1 >= static_cast<int>(Kt.size()) || fLoc2 >= static_cast<int>(Kt.size()) ||
            fLoc1 >= static_cast<int>(Wf.size()) || fLoc2 >= static_cast<int>(Wf.size()) ||
            fLoc1 >= static_cast<int>(fracElements.size()) || fLoc2 >= static_cast<int>(fracElements.size())) {
            continue;
        }

        const auto* pElem1 = fracElements[fLoc1];
        const auto* pElem2 = fracElements[fLoc2];
        if (!pElem1 || !pElem2) continue;

        // [Fix-4] 插值系数限�?+ 距离/面积最小值保�?
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
    TransmissibilityFields tTags;

    auto p_Kxx = fieldMgr.getMatrixScalar(rock_str.k_xx_tag);
    auto p_Kyy = fieldMgr.getMatrixScalar(rock_str.k_yy_tag);
    auto p_Kzz = fieldMgr.getMatrixScalar(rock_str.k_zz_tag);
    auto p_Lam_m = fieldMgr.getMatrixScalar(rock_str.lambda_tag);

    if (!p_Kxx) {
        throw std::runtime_error("[Solver 3D] Matrix permeability field k_xx is required for NNC.");
    }

    const std::vector<double>& Kxx = p_Kxx->data;
    const std::vector<double>& Kyy = (p_Kyy) ? p_Kyy->data : Kxx;
    const std::vector<double>& Kzz = (p_Kzz) ? p_Kzz->data : Kxx;
    const std::vector<double>* Lam_m = (p_Lam_m) ? &(p_Lam_m->data) : nullptr;

    auto p_Kf = fieldMgr.getFractureScalar(frac_str.k_n_tag);
    if (!p_Kf) {
        throw std::runtime_error("[Solver 3D] Fracture normal permeability k_n is required for NNC.");
    }
    auto p_Wf = fieldMgr.getFractureScalar(frac_str.aperture_tag);
    auto p_Lam_f = fieldMgr.getFractureScalar(frac_str.lambda_tag);
    auto p_Phi_f = fieldMgr.getFractureScalar(frac_str.phi_tag);
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Wf) {
        throw std::runtime_error("[Solver 3D] Fracture aperture is required for NNC.");
    }

    const std::vector<double>& Kf = p_Kf->data;
    const std::vector<double>& Wf = p_Wf->data;
    const std::vector<double>* Lam_f = (p_Lam_f) ? &(p_Lam_f->data) : nullptr;
    const std::vector<double>* Phi_f = (p_Phi_f) ? &(p_Phi_f->data) : nullptr;
    const std::vector<double>* LamFluid = (p_LamFluid) ? &(p_LamFluid->data) : nullptr;

    if (Kf.size() != Wf.size()) {
        throw std::runtime_error("[Solver 3D] NNC requires same fracture field size for k_n and aperture.");
    }
    if (Lam_m && Lam_m->size() != Kxx.size()) {
        throw std::runtime_error("[Solver 3D] Matrix lambda size mismatch with matrix permeability size.");
    }
    if (Lam_f && Lam_f->size() != Wf.size()) {
        throw std::runtime_error("[Solver 3D] Fracture lambda size mismatch with aperture size.");
    }
    if (Phi_f && Phi_f->size() != Wf.size()) {
        throw std::runtime_error("[Solver 3D] Fracture porosity size mismatch with aperture size.");
    }
    if (LamFluid && LamFluid->size() != Wf.size()) {
        throw std::runtime_error("[Solver 3D] Fluid lambda size mismatch with aperture size.");
    }

    auto& T_Flow = fieldMgr.createNNCScalar(tTags.nnc_flow, 0.0)->data;
    auto& T_Heat = fieldMgr.createNNCScalar(tTags.nnc_heat, 0.0)->data;

    const auto& pairs = meshMgr.getInteractionPairs();
    const size_t numPairs = pairs.size();
    if (T_Flow.size() != numPairs) {
        T_Flow.resize(numPairs, 0.0);
        T_Heat.resize(numPairs, 0.0);
    }

    const int nMat = static_cast<int>(meshMgr.mesh().getCells().size());

    struct NNCGeom {
        int mIdx = -1;
        int fLocIdx = -1;
        double nx = 0.0;
        double ny = 0.0;
        double nz = 1.0;
        double dist = 1e-6;
        double area = 1e-12;
    };

    std::vector<NNCGeom> geom(numPairs);

    for (size_t i = 0; i < numPairs; ++i) {
        const auto& pair = pairs[i];
        const int mIdx = pair.matrixSolverIndex;
        const int fIdx = pair.fracCellSolverIndex;
        const int fLocIdx = fIdx - nMat;

        if (mIdx < 0 || mIdx >= static_cast<int>(Kxx.size())) {
            throw std::runtime_error("[Solver 3D] NNC matrix index out of range.");
        }
        if (fLocIdx < 0 || fLocIdx >= static_cast<int>(Wf.size())) {
            throw std::runtime_error("[Solver 3D] NNC fracture local index out of range.");
        }

        const double nxRaw = pair.polygonNormal.m_x;
        const double nyRaw = pair.polygonNormal.m_y;
        const double nzRaw = pair.polygonNormal.m_z;
        const double nNorm = std::sqrt(nxRaw * nxRaw + nyRaw * nyRaw + nzRaw * nzRaw);
        if (nNorm <= 1e-14) {
            throw std::runtime_error("[Solver 3D] NNC polygon normal magnitude is zero.");
        }

        NNCGeom g;
        g.mIdx = mIdx;
        g.fLocIdx = fLocIdx;
        g.nx = nxRaw / nNorm;
        g.ny = nyRaw / nNorm;
        g.nz = nzRaw / nNorm;
        g.dist = std::max(pair.distMatrixToFracPlane, 1e-6);
        g.area = std::max(pair.intersectionArea, 1e-12);
        geom[i] = g;
    }

#pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(numPairs); ++i)
    {
        const auto& g = geom[i];

        const double k_m_dir =
            (g.nx * g.nx * Kxx[g.mIdx]) +
            (g.ny * g.ny * Kyy[g.mIdx]) +
            (g.nz * g.nz * Kzz[g.mIdx]);

        T_Flow[i] = FVM_Ops::Op_Math_Transmissibility(g.dist, k_m_dir, Wf[g.fLocIdx] * 0.5, Kf[g.fLocIdx], g.area);

        if (Lam_m && Lam_f && Phi_f && LamFluid) {
            const double lam_m_val = (*Lam_m)[g.mIdx];
            const double phi_val = (*Phi_f)[g.fLocIdx];
            const double lam_fluid_val = (*LamFluid)[g.fLocIdx];
            const double lam_solid_val = (*Lam_f)[g.fLocIdx];
            const double lam_f_eff = phi_val * lam_fluid_val + (1.0 - phi_val) * lam_solid_val;

            T_Heat[i] = FVM_Ops::Op_Math_Transmissibility(g.dist, lam_m_val, Wf[g.fLocIdx] * 0.5, lam_f_eff, g.area);
        }
        else if (Lam_m && Lam_f) {
            const double lam_m_val = (*Lam_m)[g.mIdx];
            const double lam_f_val = (*Lam_f)[g.fLocIdx];
            T_Heat[i] = FVM_Ops::Op_Math_Transmissibility(g.dist, lam_m_val, Wf[g.fLocIdx] * 0.5, lam_f_val, g.area);
        }
    }

    std::cout << "[Solver 3D] NNC Done (" << numPairs << " pairs)." << std::endl;
}


// =========================================================================
// 静态传导率计算：FF (Fracture - Fracture 3D Star-Delta Model)
// =========================================================================
/**
 * @brief 计算 3D 裂缝-裂缝 (FF) 星角变换传导�?(Industrial-Grade)
 * @details
 * 1. 端点量化无向线段键聚�?(Quantized Undirected Segment Key)，确保交线簇识别的绝对鲁棒�?
 * 2. 几何要素同步锁定，保留交线簇中最大线段的严格起止坐标�?
 * 3. 物理场解耦过滤：Flow �?Heat 各自独立判定有效性，Heat 缺失自动回退�?0，绝不影�?Flow 拓扑�?
 * 4. 强确定性排�?(Deterministic Ordering)：消除哈希迭代的随机性，保证 CSV 报表和矩阵装配绝对可复现�?
 * @param meshMgr 3D 网格管理�?
 * @param fieldMgr 3D 场管理器
 */
void TransmissibilitySolver_3D::Calculate_Transmissibility_FF(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr)
{
    std::cout << "\n[Solver 3D] Calculating FF Transmissibility (Deterministic Star-Delta)..." << std::endl;

    PhysicalProperties_string_op::Fracture_string fracStr;
    PhysicalProperties_string_op::Water waterStr;
    PhysicalProperties_string_op::TransmissibilityFields tTags;

    const auto& frNet = meshMgr.fracture_network();
    const auto& ffIntersections = frNet.ffIntersections;
    const auto& fracElements = frNet.getOrderedFractureElements();
    const auto& globalEdges = frNet.getGlobalEdges();
    const int offset = frNet.getSolverIndexOffset();

    auto p_Kt = fieldMgr.getFractureScalar(fracStr.k_t_tag);
    auto p_Wf = fieldMgr.getFractureScalar(fracStr.aperture_tag);
    auto p_LamF = fieldMgr.getFractureScalar(fracStr.lambda_tag);
    auto p_PhiF = fieldMgr.getFractureScalar(fracStr.phi_tag);
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Kt || !p_Wf) {
        throw std::runtime_error("[Solver 3D] Critical fracture properties k_t/aperture missing for FF.");
    }

    const auto& Kt = p_Kt->data;
    const auto& Wf = p_Wf->data;
    const auto& vLamF = (p_LamF) ? p_LamF->data : std::vector<double>();
    const auto& vPhiF = (p_PhiF) ? p_PhiF->data : std::vector<double>();
    const auto& vLamFluid = (p_LamFluid) ? p_LamFluid->data : std::vector<double>();
    const bool hasGlobalHeat = (!vLamF.empty() && !vPhiF.empty() && !vLamFluid.empty());

    if (Kt.size() != Wf.size()) {
        throw std::runtime_error("[Solver 3D] FF requires same fracture field size for k_t and aperture.");
    }
    if (!vLamF.empty() && vLamF.size() != Wf.size()) {
        throw std::runtime_error("[Solver 3D] FF fracture lambda size mismatch.");
    }
    if (!vPhiF.empty() && vPhiF.size() != Wf.size()) {
        throw std::runtime_error("[Solver 3D] FF fracture porosity size mismatch.");
    }
    if (!vLamFluid.empty() && vLamFluid.size() != Wf.size()) {
        throw std::runtime_error("[Solver 3D] FF fluid lambda size mismatch.");
    }

    auto makePairKey = [](int a, int b) -> std::uint64_t {
        const std::uint32_t i = static_cast<std::uint32_t>(std::min(a, b));
        const std::uint32_t j = static_cast<std::uint32_t>(std::max(a, b));
        return (static_cast<std::uint64_t>(i) << 32) | static_cast<std::uint64_t>(j);
    };

    std::unordered_set<std::uint64_t> fiPairs;
    fiPairs.reserve(globalEdges.size() * 2 + 1);
    for (const auto& edge : globalEdges) {
        const int s1 = edge.ownerCell_solverIndex;
        const int s2 = edge.neighborCell_solverIndex;
        if (s1 < 0 || s2 < 0) {
            continue;
        }
        fiPairs.insert(makePairKey(s1, s2));
    }

    struct JunctionCluster {
        Vector start;
        Vector end;
        double length = -1.0;
        std::vector<int> solverIndices;
    };

    std::unordered_map<std::string, JunctionCluster> clusterMap;
    const double tol = 1e-4;

    auto quantize = [tol](const Vector& v) -> std::tuple<long long, long long, long long> {
        return std::make_tuple(
            static_cast<long long>(std::floor(v.m_x / tol + 0.5)),
            static_cast<long long>(std::floor(v.m_y / tol + 0.5)),
            static_cast<long long>(std::floor(v.m_z / tol + 0.5))
        );
    };

    auto getUndirectedKey = [&](const Vector& p1, const Vector& p2) -> std::string {
        auto q1 = quantize(p1);
        auto q2 = quantize(p2);
        if (q1 > q2) {
            std::swap(q1, q2);
        }
        return std::to_string(std::get<0>(q1)) + "_" + std::to_string(std::get<1>(q1)) + "_" + std::to_string(std::get<2>(q1)) + "-" +
            std::to_string(std::get<0>(q2)) + "_" + std::to_string(std::get<1>(q2)) + "_" + std::to_string(std::get<2>(q2));
    };

    for (const auto& inter : ffIntersections) {
        for (const auto& seg : inter.segments) {
            std::string key = getUndirectedKey(seg.start, seg.end);
            auto& cluster = clusterMap[key];

            if (seg.length > cluster.length) {
                cluster.start = seg.start;
                cluster.end = seg.end;
                cluster.length = seg.length;
            }

            if (seg.solverIndex_1 >= 0 &&
                std::find(cluster.solverIndices.begin(), cluster.solverIndices.end(), seg.solverIndex_1) == cluster.solverIndices.end()) {
                cluster.solverIndices.push_back(seg.solverIndex_1);
            }
            if (seg.solverIndex_2 >= 0 &&
                std::find(cluster.solverIndices.begin(), cluster.solverIndices.end(), seg.solverIndex_2) == cluster.solverIndices.end()) {
                cluster.solverIndices.push_back(seg.solverIndex_2);
            }
        }
    }

    std::vector<std::string> sortedKeys;
    sortedKeys.reserve(clusterMap.size());
    for (const auto& kv : clusterMap) {
        sortedKeys.push_back(kv.first);
    }
    std::sort(sortedKeys.begin(), sortedKeys.end());

    std::vector<JunctionCluster> validClusters;
    validClusters.reserve(sortedKeys.size());

    size_t totalFFPairs = 0;
    for (const auto& key : sortedKeys) {
        auto cluster = clusterMap[key];
        std::vector<int> validIndices;
        validIndices.reserve(cluster.solverIndices.size());

        for (int sIdx : cluster.solverIndices) {
            const int fLoc = sIdx - offset;
            if (fLoc < 0 || fLoc >= static_cast<int>(Kt.size()) || fLoc >= static_cast<int>(Wf.size()) || fLoc >= static_cast<int>(fracElements.size())) {
                continue;
            }
            if (!fracElements[fLoc]) {
                continue;
            }
            validIndices.push_back(sIdx);
        }

        std::sort(validIndices.begin(), validIndices.end());
        validIndices.erase(std::unique(validIndices.begin(), validIndices.end()), validIndices.end());

        if (validIndices.size() < 2) {
            continue;
        }

        cluster.solverIndices = std::move(validIndices);
        if (cluster.length <= 1e-14) {
            cluster.length = 1e-6;
        }

        for (size_t i = 0; i < cluster.solverIndices.size(); ++i) {
            for (size_t j = i + 1; j < cluster.solverIndices.size(); ++j) {
                if (fiPairs.find(makePairKey(cluster.solverIndices[i], cluster.solverIndices[j])) != fiPairs.end()) {
                    continue;
                }
                ++totalFFPairs;
            }
        }

        validClusters.push_back(std::move(cluster));
    }

    auto p_T_Flow = fieldMgr.createFFScalar(tTags.ff_flow, 0.0);
    auto p_T_Heat = fieldMgr.createFFScalar(tTags.ff_heat, 0.0);
    auto& T_FF_Flow = p_T_Flow->data;
    auto& T_FF_Heat = p_T_Heat->data;

    if (T_FF_Flow.size() != totalFFPairs) {
        T_FF_Flow.assign(totalFFPairs, 0.0);
        T_FF_Heat.assign(totalFFPairs, 0.0);
    }

    size_t ffIdx = 0;
    fieldMgr.ff_topology.clear();
    fieldMgr.ff_topology.reserve(totalFFPairs);

    for (const auto& cluster : validClusters) {
        const size_t nElems = cluster.solverIndices.size();
        const double segLength = std::max(cluster.length, 1e-6);

        std::vector<double> half_T_Flow(nElems, 0.0);
        std::vector<double> half_T_Heat(nElems, 0.0);
        double sum_T_Flow = 0.0;
        double sum_T_Heat = 0.0;

        for (size_t i = 0; i < nElems; ++i) {
            const int fLoc = cluster.solverIndices[i] - offset;
            const auto* pElem = fracElements[fLoc];

            const double d = std::max(Geometry_3D::PointToSegmentDistance(pElem->centroid, cluster.start, cluster.end), 1e-6);

            const double cond_flow = Kt[fLoc] * Wf[fLoc];
            const double t_flow = (cond_flow * segLength) / d;
            half_T_Flow[i] = t_flow;
            sum_T_Flow += t_flow;

            const bool branchHasHeat = hasGlobalHeat &&
                (fLoc < static_cast<int>(vLamF.size())) &&
                (fLoc < static_cast<int>(vPhiF.size())) &&
                (fLoc < static_cast<int>(vLamFluid.size()));

            if (branchHasHeat) {
                const double lam_eff = vPhiF[fLoc] * vLamFluid[fLoc] + (1.0 - vPhiF[fLoc]) * vLamF[fLoc];
                const double cond_heat = lam_eff * Wf[fLoc];
                const double t_heat = (cond_heat * segLength) / d;
                half_T_Heat[i] = t_heat;
                sum_T_Heat += t_heat;
            }
        }

        for (size_t i = 0; i < nElems; ++i) {
            for (size_t j = i + 1; j < nElems; ++j) {
                const int sI = cluster.solverIndices[i];
                const int sJ = cluster.solverIndices[j];
                if (fiPairs.find(makePairKey(sI, sJ)) != fiPairs.end()) {
                    continue;
                }

                T_FF_Flow[ffIdx] = (sum_T_Flow > 1e-25) ? ((half_T_Flow[i] * half_T_Flow[j]) / sum_T_Flow) : 0.0;
                T_FF_Heat[ffIdx] = (sum_T_Heat > 1e-25) ? ((half_T_Heat[i] * half_T_Heat[j]) / sum_T_Heat) : 0.0;

                fieldMgr.ff_topology.emplace_back(sI, sJ);
                ++ffIdx;
            }
        }
    }

    if (ffIdx != totalFFPairs) {
        throw std::runtime_error("[Solver 3D] FF topology size mismatch with filtered totalFFPairs.");
    }

    std::cout << "[Solver 3D] FF Done (" << totalFFPairs << " Deterministic Star-Delta pairs over "
              << validClusters.size() << " valid junctions)." << std::endl;
}
