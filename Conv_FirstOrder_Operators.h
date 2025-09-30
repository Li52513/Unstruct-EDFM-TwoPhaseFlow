#pragma once
#include <algorithm>
#include <cmath>
#include <string>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FieldAcessForDiscre.h"
#include "TemperatureBCAdapter.h"
#include"Conv_FirstOrder_ComputerBoundaryFace.h"
#include"Conv_FirstOrder_ComputerInnerFace.h"




// ==================== 水单相能量对流项采用一阶迎风格式内部面和边界面离散系数和源项组装 ====================
inline void Convective_FirstOrder_SinglePhase_Temperature
(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const TemperatureBCAdapter& Tbc,
    const std::string& cp_field = "cp_w",
    const std::string& p_field = "p_w",
    const std::string& T_field = "T",
    const std::string& a_name = "a_f_Diff_p_w",  //注意这里名称与TPFA保持一致,输入扩散项离散系数面表，用于计算质量通量
    const std::string& s_name = "s_f_Diff_p_w",  //注意这里名称与TPFA保持一致，输入扩散项离散源项面表，用于计算质量通量
    const std::string& af_PP_name = "aPP_conv", //输出：对角线
    const std::string& af_PN_name = "aPN_conv", //输出：非对角线
    const std::string& bP_name = "bP_conv"  //输出：右端项
    
) 
{
    Convective_FirstOrder_SinglePhase_Temperature_InnerFace(mgr, reg, freg, cp_field, p_field, T_field, a_name, s_name, af_PP_name, af_PN_name, bP_name);
	Convective_FirstOrder_SinglePhase_Temperature_BoundaryFace(mgr, reg, freg, Tbc, cp_field, p_field, T_field, a_name, s_name, af_PP_name, af_PN_name, bP_name);
	

}