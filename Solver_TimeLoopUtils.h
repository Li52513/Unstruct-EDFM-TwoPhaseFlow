#pragma once
#include <iostream>
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "PhysicalProperties_CO2.h"
#include "WaterPropertyTable.h"
#include <cassert>
#include "Solver_AssemblerCOO.h" 
#include "Solver_Tools.h"


//*********补充工具函数**********//
// 1) 构建时间线性化点场（用于牛顿迭代的雅可比矩阵组装）
inline bool bulidTimeLinearizationPoint
(
    MeshManager& mgr, FieldRegistry& reg,
    const std::string& p_old_name, const std::string& T_old_name,   // n 层
    const std::string& p_iter_name = "", const std::string& T_iter_name = "", // 可选：迭代/预测层
    const std::string& p_lin_out = "p_time_lin",
    const std::string& T_lin_out = "T_time_lin"
)
{
    auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    if (cells.empty()) return false;

    //取出旧时层
    auto p0 = reg.get<volScalarField>(p_old_name);
    auto T0 = reg.get<volScalarField>(T_old_name);
    if (!p0 || !T0) { std::cerr << "[TimeLin] missing p_old/T_old\n"; return false; }

    //取出迭代层
    std::shared_ptr<volScalarField> pk, Tk;
    if (!p_iter_name.empty()) pk = reg.get<volScalarField>(p_iter_name);
    if (!T_iter_name.empty()) Tk = reg.get<volScalarField>(T_iter_name);

    // 生成线性化点场
    auto p_lin = reg.getOrCreate<volScalarField>(p_lin_out, cells.size(), 0.0);
    auto T_lin = reg.getOrCreate<volScalarField>(T_lin_out, cells.size(), 0.0);

    for (const auto& c : cells)
    {
        if (c.id < 0)continue;
        const size_t i = id2idx.at(c.id);
        double pn = (*p0)[i], Tn = (*T0)[i];
        double pl = pn, Tl = Tn;

        if (pk) {
            const double pik = (*pk)[i];
            const double Tik = Tk ? (*Tk)[i] : Tn; // 若无迭代温度场则用旧时层温度
			pl = 0.5 * (pn + pik);
			Tl = 0.5 * (Tn + Tik);
		}
        //Initializer::clampPT(pl, Tl);
        (*p_lin)[i] = pl;
        (*T_lin)[i] = Tl;
    
    }
    return true;
}

// 构造 p_eval, T_eval = (1-θ)*old + θ*iter
inline void buildEvalFields(
    FieldRegistry& reg,
    const std::string& p_old, const std::string& p_iter,
    const std::string& T_old, const std::string& T_iter,
    const std::string& p_eval_name, const std::string& T_eval_name,
    double theta_p = 1.0, double theta_T = 1.0
) {
    auto pOld = reg.get<volScalarField>(p_old);
    auto pIt = reg.get<volScalarField>(p_iter);
    auto TOld = reg.get<volScalarField>(T_old);
    auto TIt = reg.get<volScalarField>(T_iter);
    if (!pOld || !pIt || !TOld || !TIt) { throw std::runtime_error("[buildEvalFields] missing fields"); }

    theta_p = std::min(1.0, std::max(0.0, theta_p));
    theta_T = std::min(1.0, std::max(0.0, theta_T));

    auto pEval = reg.getOrCreate<volScalarField>(p_eval_name, pIt->data.size(), 0.0);
    auto TEval = reg.getOrCreate<volScalarField>(T_eval_name, TIt->data.size(), 0.0);

    for (size_t i = 0; i < pIt->data.size(); ++i) {
        (*pEval)[i] = (1.0 - theta_p) * (*pOld)[i] + theta_p * (*pIt)[i];
        (*TEval)[i] = (1.0 - theta_T) * (*TOld)[i] + theta_T * (*TIt)[i];
    }
}
