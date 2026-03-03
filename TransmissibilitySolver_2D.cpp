#include "TransmissibilitySolver_2D.h"
#include "SolverContrlStrName_op.h" // 包含 PhysicalProperties_string 定义
#include "FVM_Ops.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>

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

    const double MIN_VAL = 1e-15;
    const double thickness = 1.0;     // [Assumption] 2D 模拟单位厚度

    // --- 获取基岩场数据 ---
    auto p_Kxx = fieldMgr.getMatrixScalar(rock_str.k_xx_tag);
    auto p_Kyy = fieldMgr.getMatrixScalar(rock_str.k_yy_tag);
    auto p_Lam_m = fieldMgr.getMatrixScalar(rock_str.lambda_tag);

    if (!p_Kxx) { std::cerr << "[Error] Matrix Permeability missing!" << std::endl; return; }

    const std::vector<double>& Kxx = p_Kxx->data;
    const std::vector<double>& Kyy = (p_Kyy) ? p_Kyy->data : Kxx;
    const std::vector<double>* Lam_m = (p_Lam_m) ? &(p_Lam_m->data) : nullptr;

    // --- 获取裂缝场数据 ---
    auto p_Kf = fieldMgr.getFractureScalar(frac_str.k_n_tag);
    if (!p_Kf) p_Kf = fieldMgr.getFractureScalar(frac_str.k_t_tag); // Fallback
    auto p_Wf = fieldMgr.getFractureScalar(frac_str.aperture_tag);
    auto p_Lam_f = fieldMgr.getFractureScalar(frac_str.lambda_tag); // 裂缝固相热导率
    auto p_Phi_f = fieldMgr.getFractureScalar(frac_str.phi_tag);    // 裂缝孔隙度

    // 获取流体热导率 lambda_w
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Kf || !p_Wf) { std::cerr << "[Error] Frac properties missing!" << std::endl; return; }

    const std::vector<double>& Kf = p_Kf->data;
    const std::vector<double>& Wf = p_Wf->data;
    const std::vector<double>* Lam_f = (p_Lam_f) ? &(p_Lam_f->data) : nullptr;
    const std::vector<double>* Phi_f = (p_Phi_f) ? &(p_Phi_f->data) : nullptr;
    const std::vector<double>* LamFluid = (p_LamFluid) ? &(p_LamFluid->data) : nullptr;

    // --- 输出场 ---
    const std::string tag_T_NNC_Flow = "T_NNC_Flow";
    const std::string tag_T_NNC_Heat = "T_NNC_Heat";

    auto& T_Flow = fieldMgr.createNNCScalar(tag_T_NNC_Flow, 0.0)->data;
    auto& T_Heat = fieldMgr.createNNCScalar(tag_T_NNC_Heat, 0.0)->data;

    // --- 构建任务列表 (O(N_NNC) 极致性能优化版) ---
    struct NNCJob {
        int mIdx;         // Matrix Local Index
        int fSolverIdx;   // Frac Solver Index
        size_t nncIdx;    // Index in T_NNC array
    };

    const auto& nncMap = meshMgr.getNNCTopologyMap();
    std::vector<NNCJob> jobs;

    // 预分配准确内存，彻底杜绝 std::vector 动态扩容开销
    size_t exactNNCCount = 0;
    for (const auto& kv : nncMap) { exactNNCCount += kv.second.size(); }
    jobs.reserve(exactNNCCount);

    size_t nncIndex = 0;
    // 【核心提速】直接遍历含有裂缝的网格，跳过百万级无裂缝基岩
    for (const auto& kv : nncMap) {
        int mIdx = kv.first; // kv.first 严格对应基岩网格的局部索引 (Local Index)
        for (int fSolverIdx : kv.second) {
            jobs.push_back({ mIdx, fSolverIdx, nncIndex++ });
        }
    }

    // 动态校准 FieldManager 数组大小以匹配实际有效 NNC 数量
    if (T_Flow.size() != jobs.size()) {
        T_Flow.resize(jobs.size());
        T_Heat.resize(jobs.size());
    }

    // 获取宏观裂缝列表与基岩单元列表，以便在后续并行循环中 O(1) 访问几何
    const auto& macroFractures = meshMgr.fracture_network().fractures;
    const auto& matrixCells = meshMgr.mesh().getCells();

    // --- 计算循环 ---
#pragma omp parallel for schedule(static)
    for (long long i = 0; i < (long long)jobs.size(); ++i)
    {
        const auto& job = jobs[i];
        int mIdx = job.mIdx;
        int fIdx = job.fSolverIdx;

        // 【核心修复】将全局 solverIndex 转换为 FieldManager 中的局部裂缝索引
        int nMat = matrixCells.size();
        int fLocIdx = (fIdx >= nMat) ? (fIdx - nMat) : fIdx;

        // 越界安全保护
        if (mIdx < 0 || mIdx >= (int)Kxx.size()) continue;
        if (fLocIdx < 0 || fLocIdx >= (int)Wf.size()) continue;

        // 获取裂缝单元 (底层仍需使用全局 solverIndex 查询拓扑)
        const FractureElement* pElem = meshMgr.getFractureElementBySolverIndex(fIdx);
        if (!pElem) continue;

        // 2D FractureElement 没有直接存储坐标，需要通过 parentFractureID 和 params 还原
        if (pElem->parentFractureID < 0 || pElem->parentFractureID >= (int)macroFractures.size()) continue;

        const auto& parentFrac = macroFractures[pElem->parentFractureID];
        Vector fracStart = parentFrac.start;
        Vector fracVec = parentFrac.end - parentFrac.start;

        // 还原线段端点 p1, p2
        Vector p1 = fracStart + fracVec * pElem->param0;
        Vector p2 = fracStart + fracVec * pElem->param1;

        // 计算法向 (2D) 以用于获取渗透率张量的等效主值
        Vector tangent = p2 - p1;
        Vector normal(-tangent.m_y, tangent.m_x, 0.0);
        double len = Geometry_2D::Mag(normal);
        if (len > 1e-12) normal = normal * (1.0 / len);
        else normal = Vector(0, 0, 0);

        // 1. Matrix Directional Permeability
        double nx = normal.m_x;
        double ny = normal.m_y;
        double k_m_dir = (nx * nx * Kxx[mIdx]) + (ny * ny * Kyy[mIdx]);

        double d_m = std::max(pElem->avgDistance, 1e-6);

        double area = pElem->length * thickness;

        // 物理场必须使用局部索引 fLocIdx 访问！
        T_Flow[job.nncIdx] = FVM_Ops::Op_Math_Transmissibility(d_m, k_m_dir, Wf[fLocIdx] / 2.0, Kf[fLocIdx], area);

        // 3. Heat Transmissibility
        if (Lam_m && Lam_f && Phi_f && LamFluid) {
            double lam_m_val = (*Lam_m)[mIdx];
            double phi_val = (*Phi_f)[fLocIdx];
            double lam_fluid_val = (*LamFluid)[fLocIdx];
            double lam_solid_val = (*Lam_f)[fLocIdx];

            // 裂缝侧有效热导率
            double lam_f_eff = phi_val * lam_fluid_val + (1.0 - phi_val) * lam_solid_val;

            T_Heat[job.nncIdx] = FVM_Ops::Op_Math_Transmissibility(d_m, lam_m_val, Wf[fLocIdx] / 2.0, lam_f_eff, area);
        }
    }
    std::cout << "[Solver 2D] NNC Done (" << jobs.size() << " pairs)." << std::endl;
}// =========================================================================
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

    // --- 获取场数据 ---
    auto p_Kt = fieldMgr.getFractureScalar(fracStr.k_t_tag);
    auto p_Wf = fieldMgr.getFractureScalar(fracStr.aperture_tag);
    auto p_LamF = fieldMgr.getFractureScalar(fracStr.lambda_tag);
    auto p_PhiF = fieldMgr.getFractureScalar(fracStr.phi_tag);
    auto p_LamFluid = fieldMgr.getFractureScalar(waterStr.k_tag);

    if (!p_Kt || !p_Wf) {
        std::cerr << "[Error] Critical fracture properties missing for FF!" << std::endl;
        return;
    }

    const auto& KtF = p_Kt->data;
    const auto& Wf = p_Wf->data;
    const auto& LamF = p_LamF ? p_LamF->data : std::vector<double>();
    const auto& PhiF = p_PhiF ? p_PhiF->data : std::vector<double>();
    const auto& LamFluid = p_LamFluid ? p_LamFluid->data : std::vector<double>();

    // =====================================================================
    // 步骤 1：构建 Junction (枢纽) 到 FractureElement 的精确拓扑映射
    // Key: GlobalFFPoint ID
    // Value: vector of 全局 Solver Index (所有连接到该点的裂缝微段)
    // =====================================================================
    std::unordered_map<int, std::vector<int>> junctionMap;
    for (const auto& frac : macroFractures) {
        for (const auto& elem : frac.elements) {
            // 利用网格底层的打标，直接完成绝对精确的聚类
            if (elem.isFFatStart && elem.gIDstart >= 0) {
                junctionMap[elem.gIDstart].push_back(elem.solverIndex);
            }
            if (elem.isFFatEnd && elem.gIDend >= 0) {
                junctionMap[elem.gIDend].push_back(elem.solverIndex);
            }
        }
    }

    // =====================================================================
    // 步骤 2：统计所需的 Star-Delta 展开配对总数，并安全分配内存
    // 公式: sum( n * (n - 1) / 2 )
    // =====================================================================
    size_t totalFFPairs = 0;
    for (const auto& kv : junctionMap) {
        size_t n = kv.second.size();
        if (n > 1) totalFFPairs += n * (n - 1) / 2;
    }

    const std::string tag_T_FF_Flow = "T_FF_Flow";
    const std::string tag_T_FF_Heat = "T_FF_Heat";

    auto& T_FF_Flow = fieldMgr.createFFScalar(tag_T_FF_Flow, 0.0)->data;
    auto& T_FF_Heat = fieldMgr.createFFScalar(tag_T_FF_Heat, 0.0)->data;

    if (T_FF_Flow.size() != totalFFPairs) {
        T_FF_Flow.assign(totalFFPairs, 0.0);
        T_FF_Heat.assign(totalFFPairs, 0.0);
    }

    // =====================================================================
    // 步骤 3：遍历 Junction 节点，执行严格的星角变换数学计算
    // 为了保证科研求解在跨平台时的绝对一致性，必须对 unordered_map 的 Key 进行排序
    // =====================================================================
    std::vector<int> junctionIDs;
    junctionIDs.reserve(junctionMap.size());
    for (const auto& kv : junctionMap) junctionIDs.push_back(kv.first);
    std::sort(junctionIDs.begin(), junctionIDs.end());

    size_t ffIdx = 0;
    for (int jID : junctionIDs) {
        const auto& elemSolverIndices = junctionMap[jID];
        size_t nElems = elemSolverIndices.size();
        if (nElems < 2) continue; // 只有一条裂缝接触到此点 (无效交点) 不发生 FF 交换

        // 预计算该 Junction 下各个分支的“半传导率” (Half-Transmissibility, T_i)
        std::vector<double> half_T_Flow(nElems, 0.0);
        std::vector<double> half_T_Heat(nElems, 0.0);
        double sum_T_Flow = 0.0;
        double sum_T_Heat = 0.0;

        for (size_t i = 0; i < nElems; ++i) {
            int sIdx = elemSolverIndices[i];
            int fLocIdx = sIdx - nMat;

            if (fLocIdx < 0 || fLocIdx >= (int)KtF.size()) continue; // 安全保护

            const FractureElement* pElem = meshMgr.getFractureElementBySolverIndex(sIdx);
            if (!pElem) continue;

            // 物理距离：网格在交点处被打断，故交点正好在段的末端。中心到交点的距离严谨为 L/2
            double d = std::max(pElem->length * 0.5, 1e-12);

            // 支路渗流传导能力 T_i = (K_t * W_f * thickness) / d
            double cond = KtF[fLocIdx] * Wf[fLocIdx];
            double t_f = (cond * thickness) / d;
            half_T_Flow[i] = t_f;
            sum_T_Flow += t_f;

            // 支路热传导能力
            if (!LamF.empty() && !PhiF.empty() && !LamFluid.empty()) {
                double lam_eff = PhiF[fLocIdx] * LamFluid[fLocIdx] + (1.0 - PhiF[fLocIdx]) * LamF[fLocIdx];
                double cond_h = lam_eff * Wf[fLocIdx];
                double t_h = (cond_h * thickness) / d;
                half_T_Heat[i] = t_h;
                sum_T_Heat += t_h;
            }
        }
        fieldMgr.ff_topology.clear();
        fieldMgr.ff_topology.reserve(totalFFPairs);

        // 星角变换展开：T_ij = (T_i * T_j) / Sum(T_k)
        for (size_t i = 0; i < nElems; ++i) {
            for (size_t j = i + 1; j < nElems; ++j) {
                // Flow
                if (sum_T_Flow > 1e-25) {
                    T_FF_Flow[ffIdx] = (half_T_Flow[i] * half_T_Flow[j]) / sum_T_Flow;
                }
                else {
                    T_FF_Flow[ffIdx] = 0.0;
                }

                // Heat
                if (sum_T_Heat > 1e-25) {
                    T_FF_Heat[ffIdx] = (half_T_Heat[i] * half_T_Heat[j]) / sum_T_Heat;
                }
                else {
                    T_FF_Heat[ffIdx] = 0.0;
                }

                fieldMgr.ff_topology.emplace_back(elemSolverIndices[i], elemSolverIndices[j]);

                ffIdx++;
            }
        }
    }

    std::cout << "[Solver 2D] FF Done (" << totalFFPairs << " Star-Delta pairs over "
        << junctionIDs.size() << " junctions)." << std::endl;
}