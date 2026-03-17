#include "TransmissibilitySolver_2D.h"
#include "SolverContrlStrName_op.h" // 包含 PhysicalProperties_string 定义
#include "FVM_Ops.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <stdexcept>

using namespace PhysicalProperties_string_op;

// =========================================================
// 内部辅助函数
// =========================================================

// =========================================================================
// [Local] 本地定义几何辅助工具 (Geometry Helper for 2D)
// =========================================================================
namespace Geometry_2D {

    // 手动实现点积
    inline double Dot(const Vector& a, const Vector& b) {
        return a.m_x * b.m_x + a.m_y * b.m_y + a.m_z * b.m_z;
    }

    // 手动实现模长平方
    inline double MagSqr(const Vector& v) {
        return v.m_x * v.m_x + v.m_y * v.m_y + v.m_z * v.m_z;
    }

    // 手动实现模长
    inline double Mag(const Vector& v) {
        return std::sqrt(MagSqr(v));
    }

    /**
     * @brief 计算点 P 到线段 AB 的最短物理距离 (2D)
     * @details 用于 NNC 计算 (基岩中心到裂缝线段)
     */
    inline double PointToSegmentDistance(const Vector& P, const Vector& A, const Vector& B) {
        Vector AB = B - A;
        Vector AP = P - A;

        double lenSq = MagSqr(AB);
        if (lenSq < 1e-12) return Mag(AP); // 退化点

        double t = Dot(AP, AB) / lenSq;
        if (t < 0.0) t = 0.0;
        else if (t > 1.0) t = 1.0;

        Vector ClosestPoint = A + AB * t;
        return Mag(P - ClosestPoint);
    }
}
// ==================== 新增部分开始：2D Matrix 传导率计算 ====================
void TransmissibilitySolver_2D::Calculate_Transmissibility_Matrix(const MeshManager& meshMgr, FieldManager_2D& fieldMgr)
{
    // 1. 获取基岩网格数据与映射表
    const Mesh& mesh = meshMgr.mesh(); 
    const auto& faces = mesh.getFaces();
    const auto& cells = mesh.getCells();
    const auto& cellId2Idx = mesh.getCellId2Index();

    // 2. 获取基岩体心物性场 (消除硬编码)
    PhysicalProperties_string_op::Rock rockStr; // 实例化名称结构体
    auto Kxx = fieldMgr.getMatrixScalar(rockStr.k_xx_tag);
    auto Kyy = fieldMgr.getMatrixScalar(rockStr.k_yy_tag);
    auto Lam_m = fieldMgr.getMatrixScalar(rockStr.lambda_tag);

    // 严苛的安全检查
    if (!Kxx || !Kyy) {
        std::cerr << "[Warning] Matrix Permeability (Kxx, Kyy) not found! Skipping 2D Flow Transmissibility." << std::endl;
        return;
    }

    // 3. 获取或分配面心传导率场 (默认值 0.0)
    auto T_Flow = fieldMgr.getOrCreateMatrixFaceScalar("T_Matrix_Flow", 0.0);
    auto T_Heat = fieldMgr.getOrCreateMatrixFaceScalar("T_Matrix_Heat", 0.0);

    // 2D 模型的默认单位厚度
    const double thickness = 1.0;

    // 4. 遍历所有网格面，仅对内部面进行静态传导率计算
    for (size_t i = 0; i < faces.size(); ++i)
    {
        const Face& face = faces[i];

        // 边界面的流量受边界条件控制，此处静态传导率不作计算，维持为 0
        if (face.isBoundary()) {
            continue;
        }

        // 提取 owner 和 neighbor 在 vector 容器中的局部索引
        size_t idxO = cellId2Idx.at(face.ownerCell);
        size_t idxN = cellId2Idx.at(face.neighborCell);

        // 获取该面的单位法向量
        double nx = face.normal.m_x;
        double ny = face.normal.m_y;

        // 【核心】依据法向量计算张量渗透率在法向上的投影 (Directional Permeability)
        double KO = (nx * nx * (*Kxx)[idxO]) + (ny * ny * (*Kyy)[idxO]);
        double KN = (nx * nx * (*Kxx)[idxN]) + (ny * ny * (*Kyy)[idxN]);

        // 计算相邻网格体心到面的法向投影距离 (投影逻辑与插值权重解耦，确保绝对物理距离)
        Vector dVecO = face.midpoint - cells[idxO].center;
        Vector dVecN = cells[idxN].center - face.midpoint;

        // 限幅保护，避免极端畸变导致距离为0
        double dO = std::max(std::abs(dVecO.m_x * nx + dVecO.m_y * ny), 1e-6);
        double dN = std::max(std::abs(dVecN.m_x * nx + dVecN.m_y * ny), 1e-6);

        // 提取已算好的正交分解面积（2D中为有效正交长度）
        double normE = face.vectorE.Mag();
        double area = normE * thickness; // 转化为二维物理横截面积

        // 调用刚刚落地的 FVM 底层数学算子计算流体传导率
        (*T_Flow)[i] = FVM_Ops::Op_Math_Transmissibility(dO, KO, dN, KN, area);

        // 计算热传导率 (若存在)
        if (Lam_m) {
            double LO = (*Lam_m)[idxO];
            double LN = (*Lam_m)[idxN];
            (*T_Heat)[i] = FVM_Ops::Op_Math_Transmissibility(dO, LO, dN, LN, area);
        }
    }
}

/**
 * @brief 计算宏观裂缝内部相邻单元之间的静态传导率 (FI: Fracture Internal)
 * @details 对应于 FieldManager 中的 FractureFace 场。处理同一根宏观裂缝剖分后的线段单元间的连通性。
 * @param meshMgr 2D 网格管理器
 * @param fieldMgr 2D 场管理器
 */
void TransmissibilitySolver_2D::Calculate_Transmissibility_FractureInternal(const MeshManager& meshMgr, FieldManager_2D& fieldMgr)
{
    std::cout << "\n[Solver 2D] Calculating Fracture Internal (FI) Transmissibility..." << std::endl;

    // 1. 获取物性标签
    PhysicalProperties_string_op::Fracture_string fracStr;
    PhysicalProperties_string_op::Water waterStr;
    const double thickness = 1.0; // 2D 默认厚度

    const auto& frNet = meshMgr.fracture_network();
    const auto& macroFractures = frNet.fractures;
    const int nMat = static_cast<int>(meshMgr.mesh().getCells().size());

    // 2. 提取物理场 (注意：使用 getFractureScalar 获取单元属性)
    auto p_Kt = fieldMgr.getFractureScalar(fracStr.k_t_tag);
    auto p_Wf = fieldMgr.getFractureScalar(fracStr.aperture_tag);
    auto p_LamF = fieldMgr.getFractureScalar(fracStr.lambda_tag);
    auto p_PhiF = fieldMgr.getFractureScalar(fracStr.phi_tag);
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Kt || !p_Wf) {
        std::cerr << "[Error] Critical fracture properties (Kt or Wf) missing!" << std::endl;
        return;
    }

    const auto& Kt = p_Kt->data;
    const auto& Wf = p_Wf->data;
    const auto& vLamF = p_LamF->data;
    const auto& vPhiF = p_PhiF->data;
    const auto& vLamFluid = p_LamFluid->data;

    // 3. 创建输出场 (使用您确认的 createFractureFaceScalar 函数)
    // 注意：这里创建的是 faceScalarField，通常用于存储面上的传导率
    auto p_T_Flow = fieldMgr.createFractureFaceScalar("T_FI_Flow", 0.0);
    auto p_T_Heat = fieldMgr.createFractureFaceScalar("T_FI_Heat", 0.0);

    if (!p_T_Flow || !p_T_Heat) {
        std::cerr << "[Error] Failed to create FractureFace fields!" << std::endl;
        return;
    }

    auto& T_FI_Flow = p_T_Flow->data;
    auto& T_FI_Heat = p_T_Heat->data;

    // 4. 遍历宏观裂缝并计算
    // 在 2D 中，FI 连接的顺序严格遵循宏观裂缝 elements 的排列顺序
    size_t currentConnIdx = 0;
    for (size_t fId = 0; fId < macroFractures.size(); ++fId)
    {
        const auto& frac = macroFractures[fId];
        const auto& elements = frac.elements;
        if (elements.size() < 2) continue;

        for (size_t i = 0; i < elements.size() - 1; ++i)
        {
            const auto& e1 = elements[i];
            const auto& e2 = elements[i + 1];

            // 检查索引安全性 (此处使用 currentConnIdx 必须小于 T_FI_Flow.size())
            if (currentConnIdx >= T_FI_Flow.size()) {
                std::cerr << "[Warning] FI index out of range at Frac " << fId << std::endl;
                break;
            }

            // 全局索引转局部索引 (用于访问 Kf, Wf 等单元场)
            int s1 = e1.solverIndex - nMat;
            int s2 = e2.solverIndex - nMat;

            if (s1 < 0 || s2 < 0 || s1 >= (int)Kt.size() || s2 >= (int)Kt.size()) {
                currentConnIdx++;
                continue;
            }

            // 几何计算：1D FVM 单元中心到端点的距离
            double d1 = e1.length * 0.5;
            double d2 = e2.length * 0.5;

            // --- Flow Transmissibility ---
            // 传导能力 cond = K_tangential * W_f
            double cond1 = Kt[s1] * Wf[s1];
            double cond2 = Kt[s2] * Wf[s2];

            T_FI_Flow[currentConnIdx] = FVM_Ops::Op_Math_Transmissibility(d1, cond1, d2, cond2, thickness);

            // --- Heat Transmissibility ---
            if (!vLamF.empty() && !vPhiF.empty() && !vLamFluid.empty()) {
                double lam_eff_1 = vPhiF[s1] * vLamFluid[s1] + (1.0 - vPhiF[s1]) * vLamF[s1];
                double lam_eff_2 = vPhiF[s2] * vLamFluid[s2] + (1.0 - vPhiF[s2]) * vLamF[s2];

                double h_cond1 = lam_eff_1 * Wf[s1];
                double h_cond2 = lam_eff_2 * Wf[s2];

                T_FI_Heat[currentConnIdx] = FVM_Ops::Op_Math_Transmissibility(d1, h_cond1, d2, h_cond2, thickness);
            }

            currentConnIdx++;
        }
    }

    std::cout << "[Solver 2D] FI Done (" << currentConnIdx << " connections)." << std::endl;
}
// =========================================================================
// 静态传导率计算：NNC (Matrix - Fracture)
// =========================================================================
/**
 * @brief 计算基岩-裂缝 (NNC) 的静态传导率
 * @details 极致性能优化版，杜绝重复几何计算，直接调用预处理阶段积分生成的 avgDistance。
 * @param meshMgr 2D 网格管理器
 * @param fieldMgr 2D 场管理器
 */
void TransmissibilitySolver_2D::Calculate_Transmissibility_NNC(const MeshManager& meshMgr, FieldManager_2D& fieldMgr)
{
    Rock rock_str;
    Fracture_string frac_str;
    Water waterStr;

    const double thickness = 1.0;
    const double minDist = 1e-6;
    const double minArea = 1e-12;

    auto p_Kxx = fieldMgr.getMatrixScalar(rock_str.k_xx_tag);
    auto p_Kyy = fieldMgr.getMatrixScalar(rock_str.k_yy_tag);
    auto p_Lam_m = fieldMgr.getMatrixScalar(rock_str.lambda_tag);
    if (!p_Kxx) {
        throw std::runtime_error("[Solver 2D] Matrix permeability field k_xx is required for NNC.");
    }

    const std::vector<double>& Kxx = p_Kxx->data;
    const std::vector<double>& Kyy = (p_Kyy) ? p_Kyy->data : Kxx;
    const std::vector<double>* Lam_m = (p_Lam_m) ? &(p_Lam_m->data) : nullptr;

    auto p_Kf = fieldMgr.getFractureScalar(frac_str.k_n_tag);
    if (!p_Kf) {
        throw std::runtime_error("[Solver 2D] Fracture normal permeability k_n is required for NNC.");
    }
    auto p_Wf = fieldMgr.getFractureScalar(frac_str.aperture_tag);
    auto p_Lam_f = fieldMgr.getFractureScalar(frac_str.lambda_tag);
    auto p_Phi_f = fieldMgr.getFractureScalar(frac_str.phi_tag);
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Wf) {
        throw std::runtime_error("[Solver 2D] Fracture aperture is required for NNC.");
    }

    const std::vector<double>& Kf = p_Kf->data;
    const std::vector<double>& Wf = p_Wf->data;
    const std::vector<double>* Lam_f = (p_Lam_f) ? &(p_Lam_f->data) : nullptr;
    const std::vector<double>* Phi_f = (p_Phi_f) ? &(p_Phi_f->data) : nullptr;
    const std::vector<double>* LamFluid = (p_LamFluid) ? &(p_LamFluid->data) : nullptr;

    if (Kf.size() != Wf.size()) {
        throw std::runtime_error("[Solver 2D] NNC requires same fracture field size for k_n and aperture.");
    }
    if (Lam_m && Lam_m->size() != Kxx.size()) {
        throw std::runtime_error("[Solver 2D] Matrix lambda size mismatch with matrix permeability size.");
    }
    if (Lam_f && Lam_f->size() != Wf.size()) {
        throw std::runtime_error("[Solver 2D] Fracture lambda size mismatch with aperture size.");
    }
    if (Phi_f && Phi_f->size() != Wf.size()) {
        throw std::runtime_error("[Solver 2D] Fracture porosity size mismatch with aperture size.");
    }
    if (LamFluid && LamFluid->size() != Wf.size()) {
        throw std::runtime_error("[Solver 2D] Fluid lambda size mismatch with aperture size.");
    }

    auto& T_Flow = fieldMgr.createNNCScalar("T_NNC_Flow", 0.0)->data;
    auto& T_Heat = fieldMgr.createNNCScalar("T_NNC_Heat", 0.0)->data;

    struct NNCJob {
        int mIdx;
        int fSolverIdx;
        size_t nncIdx;
    };

    const auto& nncMap = meshMgr.getNNCTopologyMap();
    std::vector<NNCJob> jobs;

    size_t exactNNCCount = 0;
    for (const auto& kv : nncMap) {
        exactNNCCount += kv.second.size();
    }
    jobs.reserve(exactNNCCount);

    size_t nncIndex = 0;
    for (const auto& kv : nncMap) {
        const int mIdx = kv.first;
        for (int fSolverIdx : kv.second) {
            jobs.push_back({ mIdx, fSolverIdx, nncIndex++ });
        }
    }

    if (T_Flow.size() != jobs.size()) {
        T_Flow.resize(jobs.size(), 0.0);
        T_Heat.resize(jobs.size(), 0.0);
    }

    const auto& frNet = meshMgr.fracture_network();
    const auto& macroFractures = frNet.fractures;
    const auto& matrixCells = meshMgr.mesh().getCells();
    const int nMat = static_cast<int>(matrixCells.size());

    std::vector<int> mIdxBuf(jobs.size(), -1);
    std::vector<int> fLocIdxBuf(jobs.size(), -1);
    std::vector<double> nxBuf(jobs.size(), 0.0);
    std::vector<double> nyBuf(jobs.size(), 0.0);
    std::vector<double> dBuf(jobs.size(), minDist);
    std::vector<double> areaBuf(jobs.size(), minArea);

    for (size_t i = 0; i < jobs.size(); ++i) {
        const auto& job = jobs[i];
        const int mIdx = job.mIdx;
        const int fIdx = job.fSolverIdx;
        const int fLocIdx = fIdx - nMat;

        if (mIdx < 0 || mIdx >= static_cast<int>(Kxx.size())) {
            throw std::runtime_error("[Solver 2D] NNC matrix index out of range.");
        }
        if (fLocIdx < 0 || fLocIdx >= static_cast<int>(Wf.size())) {
            throw std::runtime_error("[Solver 2D] NNC fracture local index out of range.");
        }

        const FractureElement* pElem = meshMgr.getFractureElementBySolverIndex(fIdx);
        if (!pElem) {
            throw std::runtime_error("[Solver 2D] NNC cannot map solver index to fracture element.");
        }
        if (pElem->parentFractureID < 0 || pElem->parentFractureID >= static_cast<int>(macroFractures.size())) {
            throw std::runtime_error("[Solver 2D] NNC fracture element has invalid parentFractureID.");
        }

        const auto& parentFrac = macroFractures[pElem->parentFractureID];
        const Vector fracStart = parentFrac.start;
        const Vector fracVec = parentFrac.end - parentFrac.start;
        const Vector p1 = fracStart + fracVec * pElem->param0;
        const Vector p2 = fracStart + fracVec * pElem->param1;

        Vector tangent = p2 - p1;
        Vector normal(-tangent.m_y, tangent.m_x, 0.0);
        const double norm = Geometry_2D::Mag(normal);
        if (norm <= 1e-14) {
            throw std::runtime_error("[Solver 2D] NNC degenerate fracture segment yields zero normal.");
        }
        normal = normal * (1.0 / norm);

        mIdxBuf[i] = mIdx;
        fLocIdxBuf[i] = fLocIdx;
        nxBuf[i] = normal.m_x;
        nyBuf[i] = normal.m_y;
        dBuf[i] = std::max(pElem->avgDistance, minDist);
        areaBuf[i] = std::max(pElem->length * thickness, minArea);
    }

#pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(jobs.size()); ++i)
    {
        const int mIdx = mIdxBuf[i];
        const int fLocIdx = fLocIdxBuf[i];
        const double nx = nxBuf[i];
        const double ny = nyBuf[i];

        const double k_m_dir = (nx * nx * Kxx[mIdx]) + (ny * ny * Kyy[mIdx]);
        const double d_m = dBuf[i];
        const double area = areaBuf[i];

        T_Flow[jobs[i].nncIdx] = FVM_Ops::Op_Math_Transmissibility(d_m, k_m_dir, Wf[fLocIdx] * 0.5, Kf[fLocIdx], area);

        if (Lam_m && Lam_f && Phi_f && LamFluid) {
            const double lam_m_val = (*Lam_m)[mIdx];
            const double phi_val = (*Phi_f)[fLocIdx];
            const double lam_fluid_val = (*LamFluid)[fLocIdx];
            const double lam_solid_val = (*Lam_f)[fLocIdx];
            const double lam_f_eff = phi_val * lam_fluid_val + (1.0 - phi_val) * lam_solid_val;

            T_Heat[jobs[i].nncIdx] = FVM_Ops::Op_Math_Transmissibility(d_m, lam_m_val, Wf[fLocIdx] * 0.5, lam_f_eff, area);
        }
    }

    std::cout << "[Solver 2D] NNC Done (" << jobs.size() << " pairs)." << std::endl;
}
// =========================================================================
// 静态传导率计算：FF (Fracture - Fracture)
// =========================================================================
/**
 * @brief 路线 B：基于 Junction 的 FF 传导率重构 (Star-Delta 扁平化)
 * @details
 * 1. 从 MeshManager/FractureNetwork 中提取以交点 ID 为核心的 Junction 关联裂缝段集合。
 * 2. 对每个 Junction 计算每条分支到交点的星形导纳 (Flow/Heat)。
 * 3. 通过 Star-Delta 变换生成扁平化的两两 FF 连接传导率：
 *    T_ij = (G_i * G_j) / sum_k(G_k)
 * 4. 输出场 T_FF_Flow/T_FF_Heat 的长度为所有 Junction 的 pair 数总和。
 * @param meshMgr 网格管理器
 * @param fieldMgr 场管理器
 */
 // =========================================================================
 // 静态传导率计算：FF (Fracture - Fracture Star-Delta Transformation)
 // =========================================================================
 /**
  * @brief 基于星角变换 (Star-Delta) 计算微观裂缝单元间的交叉传导率
  * @details 将传统宏观的两两交叉升级为微观节点枢纽 (Junction) 的多路流量分配模型。
  * 消除虚拟节点压力，直接生成任意相交微段间的等效传导率配对。
  */
void TransmissibilitySolver_2D::Calculate_Transmissibility_FF(const MeshManager& meshMgr, FieldManager_2D& fieldMgr)
{
    std::cout << "\n[Solver 2D] Calculating FF Transmissibility (Star-Delta)..." << std::endl;

    PhysicalProperties_string_op::Fracture_string fracStr;
    PhysicalProperties_string_op::Water waterStr;
    const double thickness = 1.0;

    const auto& frNet = meshMgr.fracture_network();
    const auto& macroFractures = frNet.fractures;
    const int nMat = static_cast<int>(meshMgr.mesh().getCells().size());

    auto p_Kt = fieldMgr.getFractureScalar(fracStr.k_t_tag);
    auto p_Wf = fieldMgr.getFractureScalar(fracStr.aperture_tag);
    auto p_LamF = fieldMgr.getFractureScalar(fracStr.lambda_tag);
    auto p_PhiF = fieldMgr.getFractureScalar(fracStr.phi_tag);
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Kt || !p_Wf) {
        throw std::runtime_error("[Solver 2D] Critical fracture properties k_t/aperture missing for FF.");
    }

    const auto& KtF = p_Kt->data;
    const auto& Wf = p_Wf->data;
    const auto& LamF = p_LamF ? p_LamF->data : std::vector<double>();
    const auto& PhiF = p_PhiF ? p_PhiF->data : std::vector<double>();
    const auto& LamFluid = p_LamFluid ? p_LamFluid->data : std::vector<double>();

    if (KtF.size() != Wf.size()) {
        throw std::runtime_error("[Solver 2D] FF requires same fracture field size for k_t and aperture.");
    }
    if (!LamF.empty() && LamF.size() != Wf.size()) {
        throw std::runtime_error("[Solver 2D] FF fracture lambda size mismatch.");
    }
    if (!PhiF.empty() && PhiF.size() != Wf.size()) {
        throw std::runtime_error("[Solver 2D] FF fracture porosity size mismatch.");
    }
    if (!LamFluid.empty() && LamFluid.size() != Wf.size()) {
        throw std::runtime_error("[Solver 2D] FF fluid lambda size mismatch.");
    }

    auto makePairKey = [](int a, int b) -> std::uint64_t {
        const std::uint32_t i = static_cast<std::uint32_t>(std::min(a, b));
        const std::uint32_t j = static_cast<std::uint32_t>(std::max(a, b));
        return (static_cast<std::uint64_t>(i) << 32) | static_cast<std::uint64_t>(j);
    };

    std::unordered_set<std::uint64_t> fiPairs;
    for (const auto& frac : macroFractures) {
        if (frac.elements.size() < 2) continue;
        for (size_t i = 0; i + 1 < frac.elements.size(); ++i) {
            const int a = frac.elements[i].solverIndex;
            const int b = frac.elements[i + 1].solverIndex;
            if (a < 0 || b < 0) continue;
            fiPairs.insert(makePairKey(a, b));
        }
    }

    std::unordered_map<int, std::vector<int>> junctionMap;
    for (const auto& frac : macroFractures) {
        for (const auto& elem : frac.elements) {
            if (elem.isFFatStart && elem.gIDstart >= 0) {
                junctionMap[elem.gIDstart].push_back(elem.solverIndex);
            }
            if (elem.isFFatEnd && elem.gIDend >= 0) {
                junctionMap[elem.gIDend].push_back(elem.solverIndex);
            }
        }
    }

    struct JunctionEntry {
        int id = -1;
        std::vector<int> solverIndices;
    };

    std::vector<int> junctionIDs;
    junctionIDs.reserve(junctionMap.size());
    for (const auto& kv : junctionMap) {
        junctionIDs.push_back(kv.first);
    }
    std::sort(junctionIDs.begin(), junctionIDs.end());

    std::vector<JunctionEntry> validJunctions;
    validJunctions.reserve(junctionIDs.size());

    for (int jID : junctionIDs) {
        auto indices = junctionMap[jID];
        std::sort(indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());

        std::vector<int> valid;
        valid.reserve(indices.size());

        for (int sIdx : indices) {
            const int fLocIdx = sIdx - nMat;
            if (fLocIdx < 0 || fLocIdx >= static_cast<int>(Wf.size())) {
                continue;
            }

            const FractureElement* pElem = meshMgr.getFractureElementBySolverIndex(sIdx);
            if (!pElem) {
                continue;
            }
            if (pElem->length <= 1e-14) {
                continue;
            }

            valid.push_back(sIdx);
        }

        if (valid.size() >= 2) {
            validJunctions.push_back({ jID, std::move(valid) });
        }
    }

    size_t totalFFPairs = 0;
    for (const auto& junction : validJunctions) {
        const auto& elems = junction.solverIndices;
        for (size_t i = 0; i < elems.size(); ++i) {
            for (size_t j = i + 1; j < elems.size(); ++j) {
                if (fiPairs.find(makePairKey(elems[i], elems[j])) != fiPairs.end()) {
                    continue;
                }
                ++totalFFPairs;
            }
        }
    }

    auto& T_FF_Flow = fieldMgr.createFFScalar("T_FF_Flow", 0.0)->data;
    auto& T_FF_Heat = fieldMgr.createFFScalar("T_FF_Heat", 0.0)->data;
    if (T_FF_Flow.size() != totalFFPairs) {
        T_FF_Flow.assign(totalFFPairs, 0.0);
        T_FF_Heat.assign(totalFFPairs, 0.0);
    }

    size_t ffIdx = 0;
    fieldMgr.ff_topology.clear();
    fieldMgr.ff_topology.reserve(totalFFPairs);

    for (const auto& junction : validJunctions) {
        const auto& elemSolverIndices = junction.solverIndices;
        const size_t nElems = elemSolverIndices.size();

        std::vector<double> half_T_Flow(nElems, 0.0);
        std::vector<double> half_T_Heat(nElems, 0.0);
        double sum_T_Flow = 0.0;
        double sum_T_Heat = 0.0;

        for (size_t i = 0; i < nElems; ++i) {
            const int sIdx = elemSolverIndices[i];
            const int fLocIdx = sIdx - nMat;
            const FractureElement* pElem = meshMgr.getFractureElementBySolverIndex(sIdx);
            if (!pElem) {
                continue;
            }

            const double d = std::max(pElem->length * 0.5, 1e-12);
            const double condFlow = KtF[fLocIdx] * Wf[fLocIdx];
            const double tFlow = (condFlow * thickness) / d;
            half_T_Flow[i] = tFlow;
            sum_T_Flow += tFlow;

            if (!LamF.empty() && !PhiF.empty() && !LamFluid.empty()) {
                const double lamEff = PhiF[fLocIdx] * LamFluid[fLocIdx] + (1.0 - PhiF[fLocIdx]) * LamF[fLocIdx];
                const double condHeat = lamEff * Wf[fLocIdx];
                const double tHeat = (condHeat * thickness) / d;
                half_T_Heat[i] = tHeat;
                sum_T_Heat += tHeat;
            }
        }

        for (size_t i = 0; i < nElems; ++i) {
            for (size_t j = i + 1; j < nElems; ++j) {
                const int sI = elemSolverIndices[i];
                const int sJ = elemSolverIndices[j];
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
        throw std::runtime_error("[Solver 2D] FF generated pair count mismatch with pre-count.");
    }

    std::cout << "[Solver 2D] FF Done (" << totalFFPairs << " Star-Delta pairs over "
              << validJunctions.size() << " junctions)." << std::endl;
}