#include "TransmissibilitySolver_2D.h"
#include "SolverContrlStrName_op.h" // 包含 PhysicalProperties_string 定义

#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

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

// =========================================================================
// 静态传导率计算：NNC (Matrix - Fracture)
// =========================================================================
void TransmissibilitySolver_2D::Calculate_Transmissibility_NNC(const MeshManager& meshMgr, FieldManager_2D& fieldMgr)
{
    Rock rock_str;
    Fracture_string frac_str;
    Water waterStr;

    const double eps = 1e-30;
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

    // --- 构建任务列表 ---
    const auto& matrixCells = meshMgr.mesh().getCells();
    struct NNCJob {
        int mIdx;         // Matrix Local Index
        int fSolverIdx;   // Frac Solver Index
        size_t nncIdx;    // Index in T_NNC array
    };

    std::vector<NNCJob> jobs;
    jobs.reserve(fieldMgr.numNNCPairs);

    size_t nncIndex = 0;
    for (const auto& cell : matrixCells) {
        int mGlobalID = cell.id;
        int mIdx = meshMgr.mesh().getCellIndex(mGlobalID);

        // [Fix] 调用 meshMgr.mesh() 来访问 getFracElemSolverIndexfromCellGlobalId
        std::vector<int> fracIndices = meshMgr.mesh().getFracElemSolverIndexfromCellGlobalId(mGlobalID);

        for (int fSolverIdx : fracIndices) {
            jobs.push_back({ mIdx, fSolverIdx, nncIndex++ });
        }
    }

    if (T_Flow.size() != jobs.size()) {
        T_Flow.resize(jobs.size());
        T_Heat.resize(jobs.size());
    }

    // 获取宏观裂缝列表以便并行中访问几何
    const auto& macroFractures = meshMgr.fracture_network().fractures;

    // --- 计算循环 ---
#pragma omp parallel for schedule(static)
    for (long long i = 0; i < jobs.size(); ++i)
    {
        const auto& job = jobs[i];
        int mIdx = job.mIdx;
        int fIdx = job.fSolverIdx;

        // 获取裂缝单元
        const FractureElement* pElem = meshMgr.getFractureElementBySolverIndex(fIdx);
        if (!pElem) continue;

        // [Fix] 2D FractureElement 没有直接存储坐标，需要通过 parentFractureID 和 params 还原
        // 确保 parentFractureID 有效
        if (pElem->parentFractureID < 0 || pElem->parentFractureID >= macroFractures.size()) continue;

        const auto& parentFrac = macroFractures[pElem->parentFractureID];
        Vector fracStart = parentFrac.start;
        Vector fracVec = parentFrac.end - parentFrac.start;

        // 还原线段端点 p1, p2
        Vector p1 = fracStart + fracVec * pElem->param0;
        Vector p2 = fracStart + fracVec * pElem->param1;

        // 计算法向 (2D)
        Vector tangent = p2 - p1;
        Vector normal(-tangent.m_y, tangent.m_x, 0.0);
        double len = Geometry_2D::Mag(normal);
        if (len > 1e-12) normal = normal * (1.0 / len);
        else normal = Vector(0, 0, 0);

        // 1. Matrix Directional Permeability
        double nx = normal.m_x;
        double ny = normal.m_y;
        double k_m_dir = (nx * nx * Kxx[mIdx]) + (ny * ny * Kyy[mIdx]);

        // 2. Flow Transmissibility
        double d_m = Geometry_2D::PointToSegmentDistance(matrixCells[mIdx].center, p1, p2);
        d_m = std::max(d_m, 1e-6);

        double term_m = d_m / (k_m_dir + MIN_VAL);
        double term_f = Wf[fIdx] / (2.0 * Kf[fIdx] + MIN_VAL);

        double area = pElem->length * thickness;
        T_Flow[job.nncIdx] = area / (term_m + term_f + MIN_VAL);

        // 3. Heat Transmissibility
        if (Lam_m && Lam_f && Phi_f && LamFluid) {
            double lam_m_val = (*Lam_m)[mIdx];

            double phi_val = (*Phi_f)[fIdx];
            double lam_fluid_val = (*LamFluid)[fIdx];
            double lam_solid_val = (*Lam_f)[fIdx];

            double lam_f_eff = phi_val * lam_fluid_val + (1.0 - phi_val) * lam_solid_val;

            double term_m_h = d_m / (lam_m_val + MIN_VAL);
            double term_f_h = Wf[fIdx] / (2.0 * lam_f_eff + MIN_VAL);
            T_Heat[job.nncIdx] = area / (term_m_h + term_f_h + MIN_VAL);
        }
    }
    std::cout << "[Solver 2D] NNC Done (" << jobs.size() << " pairs)." << std::endl;
}

// =========================================================================
// 静态传导率计算：FF (Fracture - Fracture)
// =========================================================================
void TransmissibilitySolver_2D::Calculate_Transmissibility_FF(const MeshManager& meshMgr, FieldManager_2D& fieldMgr)
{
    std::cout << "\n[Solver 2D] Calculating FF Transmissibility..." << std::endl;

    Fracture_string fracStr;
    Water waterStr;
    const double eps = 1e-15;
    const double thickness = 1.0;

    // [Fix] 调用 const 版本
    const auto& frNet = meshMgr.fracture_network();
    const auto& globalFFPts = frNet.globalFFPts;

    // --- 获取场数据 ---
    auto& KtF = fieldMgr.getOrCreateFractureScalar(fracStr.k_t_tag, 1e-12)->data;
    auto& LamF = fieldMgr.getOrCreateFractureScalar(fracStr.lambda_tag, 2.0)->data;
    auto& PhiF = fieldMgr.getOrCreateFractureScalar(fracStr.phi_tag, 1.0)->data;
    auto& Wf = fieldMgr.getOrCreateFractureScalar(fracStr.aperture_tag, 1e-3)->data;

    // 获取流体热导率
    auto& LamFluid = fieldMgr.getOrCreateFractureScalar(waterStr.k_tag, 0.6)->data;

    const std::string tag_T_FF_Flow = "T_FF_Flow";
    const std::string tag_T_FF_Heat = "T_FF_Heat";

    auto& T_FF_Flow = fieldMgr.createFFScalar(tag_T_FF_Flow, 0.0)->data;
    auto& T_FF_Heat = fieldMgr.createFFScalar(tag_T_FF_Heat, 0.0)->data;

    if (T_FF_Flow.size() != globalFFPts.size()) {
        T_FF_Flow.resize(globalFFPts.size());
        T_FF_Heat.resize(globalFFPts.size());
    }

#pragma omp parallel for schedule(dynamic)
    for (long long i = 0; i < globalFFPts.size(); ++i)
    {
        const auto& ffPt = globalFFPts[i];

        if (ffPt.fracA >= frNet.fractures.size() || ffPt.fracB >= frNet.fractures.size()) continue;

        // 这里需要原始宏观裂缝对象来计算几何
        const auto& f1 = frNet.fractures[ffPt.fracA];
        const auto& f2 = frNet.fractures[ffPt.fracB];

        int segIdx1 = f1.locateSegment(ffPt.paramA);
        int segIdx2 = f2.locateSegment(ffPt.paramB);

        if (segIdx1 < 0 || segIdx1 >= f1.elements.size()) continue;
        if (segIdx2 < 0 || segIdx2 >= f2.elements.size()) continue;

        const auto& elem1 = f1.elements[segIdx1];
        const auto& elem2 = f2.elements[segIdx2];

        int s1 = elem1.solverIndex;
        int s2 = elem2.solverIndex;

        if (s1 < 0 || s2 < 0 || s1 >= static_cast<int>(KtF.size()) || s2 >= static_cast<int>(KtF.size())) continue;

        // [Fix] 计算 elem1 和 elem2 的几何中心 (center)
        // 使用 elem1.param0, param1 和 f1 的起终点
        Vector v1 = f1.end - f1.start;
        Vector c1 = f1.start + v1 * ((elem1.param0 + elem1.param1) * 0.5);

        Vector v2 = f2.end - f2.start;
        Vector c2 = f2.start + v2 * ((elem2.param0 + elem2.param1) * 0.5);

        // 计算几何距离
        double d1 = Geometry_2D::Mag(c1 - ffPt.point);
        double d2 = Geometry_2D::Mag(c2 - ffPt.point);

        d1 = std::max(d1, 1e-6);
        d2 = std::max(d2, 1e-6);

        // --- Flow Transmissibility ---
        double cond1 = Wf[s1] * KtF[s1];
        double cond2 = Wf[s2] * KtF[s2];

        double res_term = (d1 / (cond1 + eps)) + (d2 / (cond2 + eps));
        T_FF_Flow[i] = thickness / (res_term + eps);

        // --- Heat Transmissibility ---
        double lam_eff_1 = PhiF[s1] * LamFluid[s1] + (1.0 - PhiF[s1]) * LamF[s1];
        double lam_eff_2 = PhiF[s2] * LamFluid[s2] + (1.0 - PhiF[s2]) * LamF[s2];

        double cond1_h = Wf[s1] * lam_eff_1;
        double cond2_h = Wf[s2] * lam_eff_2;

        double res_term_h = (d1 / (cond1_h + eps)) + (d2 / (cond2_h + eps));
        T_FF_Heat[i] = thickness / (res_term_h + eps);
    }
    std::cout << "[Solver 2D] FF Done (" << globalFFPts.size() << " connections)." << std::endl;
}