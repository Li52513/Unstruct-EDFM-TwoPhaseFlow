#pragma once
#include <algorithm>
#include <cmath>
#include <string>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"


// —— 小工具：硬阈值/平滑Heaviside —— //
inline double heaviside_inner(double m, double eps = 1e-30) {
    return (m > eps) ? 1.0 : 0.0; // 简洁稳健的版本；需要平滑可改成tanh
}

/***************计算单相流动条件下的能量守恒方程采用一阶迎风格式的离散系数和源项，内部面************/

// —— 内部面：生成一阶迎风格式下计算得到的 m_f （只对内部面赋值）—— //

inline void Convective_FirstOrder_SinglePhase_Temperature_InnerFace
(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
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
	double epsH = 1e-30; // Heaviside 阈值
    Mesh& mesh = const_cast<Mesh&>(mgr.mesh());
    auto faces = const_cast<std::vector<Face>&>(mesh.getFaces());

	//输入面表（来自之前的TPFA计算）
    auto a_f = freg.get<faceScalarField>(a_name);
    auto s_f = freg.get<faceScalarField>(s_name);
    if (!a_f || !s_f) {
        std::cerr << "[Conv-Inner] missing TPFA face fields '" << a_name << "' or '" << s_name << "'.\n";
        return;
    }

    // 输出（面账本）
	auto m_f = freg.getOrCreate<faceScalarField>("m_f", faces.size(), 0.0);      // 质量通量
	auto theta_f = freg.getOrCreate<faceScalarField>("theta_f", faces.size(), 0.0); // 迎风选择
	auto cp_up_f = freg.getOrCreate<faceScalarField>("cp_up", faces.size(), 0.0); // 迎风比热容
	auto aPP_conv = freg.getOrCreate<faceScalarField>(af_PP_name, faces.size(), 0.0); // 对角线
	auto aPN_conv = freg.getOrCreate<faceScalarField>(af_PN_name, faces.size(), 0.0); // 非对角线
	auto bP_conv = freg.getOrCreate<faceScalarField>(bP_name, faces.size(), 0.0); // 右端项

    // 清零（防叠加）
    std::fill(aPP_conv->data.begin(), aPP_conv->data.end(), 0.0);
    std::fill(aPN_conv->data.begin(), aPN_conv->data.end(), 0.0);
    std::fill(bP_conv->data.begin(), bP_conv->data.end(), 0.0);

    const auto& id2idx = mesh.getCellId2Index();

    for (const auto& F : faces) 
    {
        const int idx = F.id - 1;
        if (F.isBoundary()) { (*m_f)[idx] = 0.0; (*theta_f)[idx] = 0.0; (*cp_up_f)[idx] = 0.0; continue; }


        const int P = F.ownerCell;
        const int N = F.neighborCell;
        const size_t iP = id2idx.at(P);
        const size_t iN = id2idx.at(N);

        const double pP = cellScalar(reg, mesh, p_field.c_str(), P, 0.0);
        const double pN = cellScalar(reg, mesh, p_field.c_str(), N, 0.0);
        // —— 质量通量：完全复用 TPFA —— //
        const double mf = (*a_f)[idx] * (pP - pN) - (*s_f)[idx];
        //const double mf = (*a_f)[idx] * (pP - pN);
        (*m_f)[idx] = mf;

        // —— 迎风选择 —— //
        const double theta = heaviside_inner(mf, epsH);
        (*theta_f)[idx] = theta;

        const double cpP = cellScalar(reg, mesh, cp_field.c_str(), P, 4200.0);
        const double cpN = cellScalar(reg, mesh, cp_field.c_str(), N, 4200.0);
        const double cp_up = theta * cpP + (1.0 - theta) * cpN;
        (*cp_up_f)[idx] = cp_up;

        // —— 面对流子算子对 P 方程的系数 —— //
        if (theta >= 0.5) { // 出流：上风=P
            (*aPP_conv)[idx] = mf * cpP;
            (*aPN_conv)[idx] = 0.0;
        }
        else {            // 入流：上风=N
            (*aPP_conv)[idx] = 0.0;
            (*aPN_conv)[idx] = mf * cpN; // 注意 mf<=0，符号自然正确
        }
        (*bP_conv)[idx] = 0.0;
    }
}

