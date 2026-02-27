#include "3D_CO2_Properties.h"
#include <iostream>
#include <algorithm>

using namespace PhysicalProperties_string_op;

// =========================================================
// 构造函数
// =========================================================
CO2Properties_3D::CO2Properties_3D(FieldManager_3D& fieldMgr,
    const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig)
    : fieldMgr_(fieldMgr), pConfig_(pConfig), tConfig_(tConfig)
{
}

// =========================================================
// 1. Matrix Domain Implementation
// =========================================================
void CO2Properties_3D::UpdateMatrix_Constant(const CO2Properties& params) 
{
    CO2 cTags;
    auto rho = fieldMgr_.getOrCreateMatrixScalar(cTags.rho_tag);
    auto mu = fieldMgr_.getOrCreateMatrixScalar(cTags.mu_tag);
    auto cp = fieldMgr_.getOrCreateMatrixScalar(cTags.cp_tag);
    auto cv = fieldMgr_.getOrCreateMatrixScalar(cTags.cv_tag);
    auto h = fieldMgr_.getOrCreateMatrixScalar(cTags.h_tag);
    auto k = fieldMgr_.getOrCreateMatrixScalar(cTags.k_tag);

    int n = static_cast<int>(rho->data.size());
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        rho->data[i] = params.rho;  
        mu->data[i] = params.mu;   
        cp->data[i] = params.cp;  
        cv->data[i] = params.cv;   
        h->data[i] = params.h;   
        k->data[i] = params.k;    
    }
}

void CO2Properties_3D::UpdateMatrix_SpanWagner()
{
    // 1. 获取输入状态变量 (基于注入的配置)
    // 使用 getOrCreate 防止由空指针引发的崩溃，给予合理的默认值
    auto& P_vec = fieldMgr_.getOrCreateMatrixScalar(pConfig_.pressure_field, 101325.0)->data;
    auto& T_vec = fieldMgr_.getOrCreateMatrixScalar(tConfig_.temperatue_field, 293.15)->data;

    // 2. 获取输出物性场
    CO2 co2Str;
    auto& rho = fieldMgr_.getOrCreateMatrixScalar(co2Str.rho_tag)->data;
    auto& mu = fieldMgr_.getOrCreateMatrixScalar(co2Str.mu_tag)->data;
    auto& cp = fieldMgr_.getOrCreateMatrixScalar(co2Str.cp_tag)->data;
    auto& lam = fieldMgr_.getOrCreateMatrixScalar(co2Str.k_tag)->data;
    auto& h = fieldMgr_.getOrCreateMatrixScalar(co2Str.h_tag)->data;

    // 3. 获取插值表单例
    CO2PropertyTable& table = CO2PropertyTable::instance();

    size_t n = P_vec.size();

    // 4. 遍历更新
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p_val = P_vec[i];
        double t_val = T_vec[i];

        try {
            // 调用 Span-Wagner 插值接口
            auto props = table.getProperties(p_val, t_val);

            rho[i] = props.rho;
            mu[i] = props.mu;
            cp[i] = props.cp;
            lam[i] = props.k; // 注意：CO2Properties 结构体中成员名为 k
            h[i] = props.h;
        }
        catch (...) {
            // 异常兜底策略：保持原值或设为安全默认值
            // 此处保持静默，生产环境建议统计错误计数
        }
    }
}

void CO2Properties_3D::CalculateCompressibilityMatrix()
{
    CO2 co2Str;
    // 输入场
    auto& P_vec = fieldMgr_.getMatrixScalar(pConfig_.pressure_field)->data;
    auto& T_vec = fieldMgr_.getMatrixScalar(tConfig_.temperatue_field)->data;

    // 输出场
    auto& drho = fieldMgr_.getOrCreateMatrixScalar(co2Str.drho_dp_tag)->data;

    CO2PropertyTable& table = CO2PropertyTable::instance();
    size_t n = P_vec.size();
    double dp = deltaP_;

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = P_vec[i];
        double t = T_vec[i];

        // 中心差分: (rho(p+dp) - rho(p-dp)) / 2dp
        try {
            double r_plus = table.getProperties(p + dp, t).rho;
            double r_minus = table.getProperties(p - dp, t).rho;
            drho[i] = (r_plus - r_minus) / (2.0 * dp);
        }
        catch (...) {
            drho[i] = 0.0; // 边界异常设为 0 (不可压)
        }
    }
}

void CO2Properties_3D::UpdateEffectiveThermalPropertiesMatrix()
{
    // 标签定义
    CO2 co2Str;
    Rock rockStr;
    EffectiveProps effStr;

    // 获取输入场 (Matrix)
    // 必须确保 Step 1 (Rock) 已经执行，否则这些场可能为空或零
    auto& phi = fieldMgr_.getMatrixScalar(rockStr.phi_tag)->data;
    auto& rho_r = fieldMgr_.getMatrixScalar(rockStr.rho_tag)->data;
    auto& cp_r = fieldMgr_.getMatrixScalar(rockStr.cp_tag)->data;
    auto& lam_r = fieldMgr_.getMatrixScalar(rockStr.lambda_tag)->data;

    auto& rho_g = fieldMgr_.getMatrixScalar(co2Str.rho_tag)->data;
    auto& cp_g = fieldMgr_.getMatrixScalar(co2Str.cp_tag)->data;
    auto& lam_g = fieldMgr_.getMatrixScalar(co2Str.k_tag)->data;

    // 获取输出场
    auto& C_eff = fieldMgr_.getOrCreateMatrixScalar(effStr.C_eff_tag)->data;
    auto& L_eff = fieldMgr_.getOrCreateMatrixScalar(effStr.lambda_eff_tag)->data;

    size_t n = phi.size();

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = phi[i];

        // 体积热容混合: (1-phi)*rho_r*cp_r + phi*rho_g*cp_g
        double C_rock = (1.0 - p) * rho_r[i] * cp_r[i];
        double C_fluid = p * rho_g[i] * cp_g[i];
        C_eff[i] = C_rock + C_fluid;

        // 导热系数混合 (简单体积平均)
        L_eff[i] = (1.0 - p) * lam_r[i] + p * lam_g[i];
    }
}

// =========================================================
// 2. Fracture Domain Implementation
// =========================================================

void CO2Properties_3D::UpdateFracture_Constant(const CO2Properties& params) 
{
    CO2 cTags;
    auto rho = fieldMgr_.getOrCreateFractureScalar(cTags.rho_tag);
    auto mu = fieldMgr_.getOrCreateFractureScalar(cTags.mu_tag);
    auto cp = fieldMgr_.getOrCreateFractureScalar(cTags.cp_tag);
    auto cv = fieldMgr_.getOrCreateFractureScalar(cTags.cv_tag);
    auto h = fieldMgr_.getOrCreateFractureScalar(cTags.h_tag);
    auto k = fieldMgr_.getOrCreateFractureScalar(cTags.k_tag);

    int n = static_cast<int>(rho->data.size());
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        rho->data[i] = params.rho;  
        mu->data[i] = params.mu;   
        cp->data[i] = params.cp;   
        cv->data[i] = params.cv;   
        h->data[i] = params.h;    
        k->data[i] = params.k;    
    }
}

void CO2Properties_3D::UpdateFracture_SpanWagner()
{
    // 获取 Fracture 域的状态变量
    // 注意：Fracture 域的 P/T 场名称通常与 Matrix 域相同，但在 fractureFields 注册表中
    auto& P_vec = fieldMgr_.getOrCreateFractureScalar(pConfig_.pressure_field, 101325.0)->data;
    auto& T_vec = fieldMgr_.getOrCreateFractureScalar(tConfig_.temperatue_field, 293.15)->data;

    CO2 co2Str;
    auto& rho = fieldMgr_.getOrCreateFractureScalar(co2Str.rho_tag)->data;
    auto& mu = fieldMgr_.getOrCreateFractureScalar(co2Str.mu_tag)->data;
    auto& cp = fieldMgr_.getOrCreateFractureScalar(co2Str.cp_tag)->data;
    auto& lam = fieldMgr_.getOrCreateFractureScalar(co2Str.k_tag)->data;
    auto& h = fieldMgr_.getOrCreateFractureScalar(co2Str.h_tag)->data;

    CO2PropertyTable& table = CO2PropertyTable::instance();
    size_t n = P_vec.size();

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        try {
            auto props = table.getProperties(P_vec[i], T_vec[i]);
            rho[i] = props.rho;
            mu[i] = props.mu;
            cp[i] = props.cp;
            lam[i] = props.k;
            h[i] = props.h;
        }
        catch (...) {}
    }
}

void CO2Properties_3D::CalculateCompressibilityFracture()
{
    CO2 co2Str;
    auto& P_vec = fieldMgr_.getFractureScalar(pConfig_.pressure_field)->data;
    auto& T_vec = fieldMgr_.getFractureScalar(tConfig_.temperatue_field)->data;
    auto& drho = fieldMgr_.getOrCreateFractureScalar(co2Str.drho_dp_tag)->data;

    CO2PropertyTable& table = CO2PropertyTable::instance();
    size_t n = P_vec.size();
    double dp = deltaP_;

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        try {
            double r_plus = table.getProperties(P_vec[i] + dp, T_vec[i]).rho;
            double r_minus = table.getProperties(P_vec[i] - dp, T_vec[i]).rho;
            drho[i] = (r_plus - r_minus) / (2.0 * dp);
        }
        catch (...) { drho[i] = 0.0; }
    }
}

void CO2Properties_3D::UpdateEffectiveThermalPropertiesFracture()
{
    CO2 co2Str;
    Fracture_string fracStr; 
    EffectiveProps effStr;

    // 输入场
    auto& phi = fieldMgr_.getFractureScalar(fracStr.phi_tag)->data;
    auto& rho_r = fieldMgr_.getFractureScalar(fracStr.rho_tag)->data;
    auto& cp_r = fieldMgr_.getFractureScalar(fracStr.cp_tag)->data;
    auto& lam_r = fieldMgr_.getFractureScalar(fracStr.lambda_tag)->data;

    auto& rho_g = fieldMgr_.getFractureScalar(co2Str.rho_tag)->data;
    auto& cp_g = fieldMgr_.getFractureScalar(co2Str.cp_tag)->data;
    auto& lam_g = fieldMgr_.getFractureScalar(co2Str.k_tag)->data;

    // 输出场
    auto& C_eff = fieldMgr_.getOrCreateFractureScalar(effStr.C_eff_tag)->data;
    auto& L_eff = fieldMgr_.getOrCreateFractureScalar(effStr.lambda_eff_tag)->data;

    size_t n = phi.size();

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = phi[i];
        double C_rock = (1.0 - p) * rho_r[i] * cp_r[i];
        double C_fluid = p * rho_g[i] * cp_g[i];
        C_eff[i] = C_rock + C_fluid;
        L_eff[i] = (1.0 - p) * lam_r[i] + p * lam_g[i];
    }
}
