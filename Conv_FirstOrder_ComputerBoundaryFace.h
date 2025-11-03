#pragma once
#include <algorithm>
#include <cmath>
#include <string>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"
#include "TemperatureBCAdapter.h"

// —— 小工具：硬阈值/平滑Heaviside —— //
inline double heaviside(double m, double eps = 1e-30) {
    return (m > eps) ? 1.0 : 0.0; // 简洁稳健的版本；需要平滑可改成tanh
}

/***************计算单相流动条件下的能量守恒方程采用一阶迎风格式的离散系数和源项，边界面************/

// —— 边界面：用TPFA给的 m_f^B = a_f*pP + s_f 判断入/出流，入流取 T_in（来自温度BC）；只对边界面赋值 —— //

inline void Convective_FirstOrder_SinglePhase_Temperature_BoundaryFace
(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const TemperatureBCAdapter& Tbc,     // 用于取入流温度：Dirichlet时 c 即 T_in；Robin时近似 T_in=c/a
    const std::string& cp_field = "cp_w",
    const std::string& p_field = "p_w",
    const std::string& T_field = "T",
    const std::string& a_name = "a_f",  //注意这里名称与TPFA保持一致,输入扩散项离散系数面表，用于计算质量通量
    const std::string& s_name = "s_f",  //注意这里名称与TPFA保持一致，输入扩散项离散源项面表，用于计算质量通量
    const std::string& af_PP_name = "aPP_conv", //输出：对角线
    const std::string& af_PN_name = "aPN_conv", //输出：非对角线
    const std::string& bP_name = "bP_conv"  //输出：右端项
)
{
    double epsH = 1e-30;
	Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
	auto faces = const_cast<std::vector<Face>&>(mesh.getFaces());
    const auto& id2idx = mesh.getCellId2Index();

	//输出字段（面场）
    auto a_f = freg.get<faceScalarField>(a_name);
	auto s_f = freg.get<faceScalarField>(s_name);

    if (!a_f || !s_f) {
        std::cerr << "[Conv-Boundary] missing TPFA face fields '" << a_name << "' or '" << s_name << "'.\n";
        return;
    }
	auto m_f = freg.getOrCreate<faceScalarField>("m_f", faces.size(), 0.0); // 面质量流量
	auto theta_f = freg.getOrCreate<faceScalarField>("theta_f", faces.size(), 0.0); // 迎风选择
	auto cp_up = freg.getOrCreate<faceScalarField>("cp_up", faces.size(), 0.0); // 迎风比热
	auto aPP_conv = freg.getOrCreate<faceScalarField>(af_PP_name, faces.size(), 0.0); // 对角线
    auto aPN_conv = freg.getOrCreate<faceScalarField>(af_PN_name, faces.size(), 0.0); // 相邻项
	auto bP_conv = freg.getOrCreate<faceScalarField>(bP_name, faces.size(), 0.0); // 源项


    for (const auto& F : faces)
    {
		const int idx = F.id - 1;
        if (!F.isBoundary())  continue;

        const int P = F.ownerCell;
        const size_t iP = id2idx.at(P);

        const double pP = cellScalar(reg, mesh, p_field.c_str(), P, 0.0);
        const double mf = (*a_f)[idx] * pP - (*s_f)[idx]; // 边界质量通量（owner→外）
        (*m_f)[idx] = mf;

        const double cpP = cellScalar(reg, mesh, cp_field.c_str(), P, 4200.0);
        // 出/入流判别（owner→外 为正）
        const double theta = heaviside(mf, epsH); // 出流=1, 入流=0
        (*theta_f)[idx] = theta;
        if (theta >= 0.5) {
            // 出流：面对 P 方程，相当于把携能带走 → 对角加 mf*cpP
            (*cp_up)[idx] = cpP;
            (*aPP_conv)[idx] += mf * cpP;
            // aPN_conv无，bP_conv无
        }
        else {
            // 入流：需入流温度 T_in（来自温度BC）；cp_in 默认用 cpP（也可扩展单独给定）
            double a = 0, b = 0, c = 0;
            double T_in = cellScalar(reg, mesh, T_field.c_str(), P, 0.0); // 回退：若无BC，取本地T
            if (Tbc.getABC(F.id, a, b, c)) {
                if (std::abs(a) > 0.0) T_in = c / a; // Dirichlet: a=1,b=0 → c 即 T_b
                // Neumann: a=0,b=1 → 纯梯度，无法直接提供入口温度，这里仍用回退 T_in
            }
            const double cp_in = cpP; // 简洁一致：cp按入流同相同流体，先用owner的cp；若将来需要，可扩展 BC 里给定

            (*cp_up)[idx] = cp_in;
            (*bP_conv)[idx] -= mf * cp_in * T_in; // 注意 mf<=0：符号自然正确  为什么不是(*bP_conv)[idx] -= mf * cp_in * T_in
            // aPP_conv、aPN_conv对入流不加（上风在域外）
        }

    }
}