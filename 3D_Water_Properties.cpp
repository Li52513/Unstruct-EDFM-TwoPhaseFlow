#include "3D_Water_Properties.h"
#include <iostream>
#include <algorithm>

using namespace PhysicalProperties_string_op;

// =========================================================
// 构造函数
// =========================================================
WaterProperties_3D::WaterProperties_3D(FieldManager_3D& fieldMgr,
    const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig)
    : fieldMgr_(fieldMgr), pConfig_(pConfig), tConfig_(tConfig)
{
}

// =========================================================
// 1. Matrix Domain Implementation
// =========================================================

void WaterProperties_3D::UpdateMatrix_Constant(const WaterProperties& params) // [修正此处]
{
    Water wTags;
    auto rho = fieldMgr_.getOrCreateMatrixScalar(wTags.rho_tag);
    auto mu = fieldMgr_.getOrCreateMatrixScalar(wTags.mu_tag);
    auto cp = fieldMgr_.getOrCreateMatrixScalar(wTags.cp_tag);
    auto cv = fieldMgr_.getOrCreateMatrixScalar(wTags.cv_tag);
    auto h = fieldMgr_.getOrCreateMatrixScalar(wTags.h_tag);
    auto k = fieldMgr_.getOrCreateMatrixScalar(wTags.k_tag);

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

void WaterProperties_3D::UpdateMatrix_IAPWS()
{
    // 1. 获取输入状态变量 (基于注入的配置)
    auto& P_vec = fieldMgr_.getOrCreateMatrixScalar(pConfig_.pressure_field, 101325.0)->data;
    auto& T_vec = fieldMgr_.getOrCreateMatrixScalar(tConfig_.temperatue_field, 293.15)->data;

    // 2. 获取输出场
    Water waterStr;
    auto& rho = fieldMgr_.getOrCreateMatrixScalar(waterStr.rho_tag)->data;
    auto& mu = fieldMgr_.getOrCreateMatrixScalar(waterStr.mu_tag)->data;
    auto& cp = fieldMgr_.getOrCreateMatrixScalar(waterStr.cp_tag)->data;
    auto& lam = fieldMgr_.getOrCreateMatrixScalar(waterStr.k_tag)->data;
    auto& h = fieldMgr_.getOrCreateMatrixScalar(waterStr.h_tag)->data;

    // 3. 获取插值表单例
    WaterPropertyTable& table = WaterPropertyTable::instance();

    int n = static_cast<int>(P_vec.size());

    // 4. 遍历更新
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p_val = P_vec[i];
        double t_val = T_vec[i];

        try {
            auto props = table.getProperties(p_val, t_val);

            rho[i] = props.rho;
            mu[i] = props.mu;
            cp[i] = props.cp;
            lam[i] = props.k; // WaterProperties::k 是导热系数
            h[i] = props.h;
        }
        catch (...) {
            // 异常兜底 (保持原值或设为常数，这里不做破坏性操作)
        }
    }
}

void WaterProperties_3D::CalculateCompressibilityMatrix()
{
    Water waterStr;
    // 输入场
    auto& P_vec = fieldMgr_.getMatrixScalar(pConfig_.pressure_field)->data;
    auto& T_vec = fieldMgr_.getMatrixScalar(tConfig_.temperatue_field)->data;

    // 输出场
    auto& drho = fieldMgr_.getOrCreateMatrixScalar(waterStr.drho_dp_tag)->data;

    WaterPropertyTable& table = WaterPropertyTable::instance();
    int n = static_cast<int>(P_vec.size());
    double dp = deltaP_;

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = P_vec[i];
        double t = T_vec[i];

        try {
            double r_plus = table.getProperties(p + dp, t).rho;
            double r_minus = table.getProperties(p - dp, t).rho;
            drho[i] = (r_plus - r_minus) / (2.0 * dp);
        }
        catch (...) {
            drho[i] = 0.0;
        }
    }
}

void WaterProperties_3D::UpdateEffectiveThermalPropertiesMatrix()
{
    // 标签定义
    Water waterStr;
    Rock rockStr;
    EffectiveProps effStr;

    // 获取输入场 (Matrix)
    auto& phi = fieldMgr_.getMatrixScalar(rockStr.phi_tag)->data;
    auto& rho_r = fieldMgr_.getMatrixScalar(rockStr.rho_tag)->data;
    auto& cp_r = fieldMgr_.getMatrixScalar(rockStr.cp_tag)->data;
    auto& lam_r = fieldMgr_.getMatrixScalar(rockStr.lambda_tag)->data;

    auto& rho_w = fieldMgr_.getMatrixScalar(waterStr.rho_tag)->data;
    auto& cp_w = fieldMgr_.getMatrixScalar(waterStr.cp_tag)->data;
    auto& lam_w = fieldMgr_.getMatrixScalar(waterStr.k_tag)->data;

    // 获取输出场
    auto& C_eff = fieldMgr_.getOrCreateMatrixScalar(effStr.C_eff_tag)->data;
    auto& L_eff = fieldMgr_.getOrCreateMatrixScalar(effStr.lambda_eff_tag)->data;

    int n = static_cast<int>(phi.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = phi[i];

        // 体积热容混合 (假设孔隙充满水，或作为单相参考)
        double C_rock = (1.0 - p) * rho_r[i] * cp_r[i];
        double C_fluid = p * rho_w[i] * cp_w[i];
        C_eff[i] = C_rock + C_fluid;

        // 导热系数混合
        L_eff[i] = (1.0 - p) * lam_r[i] + p * lam_w[i];
    }
}

// =========================================================
// 2. Fracture Domain Implementation
// =========================================================

void WaterProperties_3D::UpdateFracture_Constant(const WaterProperties& params) // [修正此处]
{
    Water wTags;
    auto rho = fieldMgr_.getOrCreateFractureScalar(wTags.rho_tag);
    auto mu = fieldMgr_.getOrCreateFractureScalar(wTags.mu_tag);
    auto cp = fieldMgr_.getOrCreateFractureScalar(wTags.cp_tag);
    auto cv = fieldMgr_.getOrCreateFractureScalar(wTags.cv_tag);
    auto h = fieldMgr_.getOrCreateFractureScalar(wTags.h_tag);
    auto k = fieldMgr_.getOrCreateFractureScalar(wTags.k_tag);

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

void WaterProperties_3D::UpdateFracture_IAPWS()
{
    auto& P_vec = fieldMgr_.getOrCreateFractureScalar(pConfig_.pressure_field, 101325.0)->data;
    auto& T_vec = fieldMgr_.getOrCreateFractureScalar(tConfig_.temperatue_field, 293.15)->data;

    Water waterStr;
    auto& rho = fieldMgr_.getOrCreateFractureScalar(waterStr.rho_tag)->data;
    auto& mu = fieldMgr_.getOrCreateFractureScalar(waterStr.mu_tag)->data;
    auto& cp = fieldMgr_.getOrCreateFractureScalar(waterStr.cp_tag)->data;
    auto& lam = fieldMgr_.getOrCreateFractureScalar(waterStr.k_tag)->data;
    auto& h = fieldMgr_.getOrCreateFractureScalar(waterStr.h_tag)->data;

    WaterPropertyTable& table = WaterPropertyTable::instance();
    int n = static_cast<int>(P_vec.size());

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

void WaterProperties_3D::CalculateCompressibilityFracture()
{
    Water waterStr;
    auto& P_vec = fieldMgr_.getFractureScalar(pConfig_.pressure_field)->data;
    auto& T_vec = fieldMgr_.getFractureScalar(tConfig_.temperatue_field)->data;
    auto& drho = fieldMgr_.getOrCreateFractureScalar(waterStr.drho_dp_tag)->data;

    WaterPropertyTable& table = WaterPropertyTable::instance();
    int n = static_cast<int>(P_vec.size());
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

void WaterProperties_3D::UpdateEffectiveThermalPropertiesFracture()
{
    Water waterStr;
    Fracture_string fracStr; // [Fix] 使用 Fracture_string
    EffectiveProps effStr;

    // 输入场
    auto& phi = fieldMgr_.getFractureScalar(fracStr.phi_tag)->data;
    auto& rho_r = fieldMgr_.getFractureScalar(fracStr.rho_tag)->data;
    auto& cp_r = fieldMgr_.getFractureScalar(fracStr.cp_tag)->data;
    auto& lam_r = fieldMgr_.getFractureScalar(fracStr.lambda_tag)->data;

    auto& rho_w = fieldMgr_.getFractureScalar(waterStr.rho_tag)->data;
    auto& cp_w = fieldMgr_.getFractureScalar(waterStr.cp_tag)->data;
    auto& lam_w = fieldMgr_.getFractureScalar(waterStr.k_tag)->data;

    // 输出场
    auto& C_eff = fieldMgr_.getOrCreateFractureScalar(effStr.C_eff_tag)->data;
    auto& L_eff = fieldMgr_.getOrCreateFractureScalar(effStr.lambda_eff_tag)->data;

    int n = static_cast<int>(phi.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = phi[i];
        double C_rock = (1.0 - p) * rho_r[i] * cp_r[i];
        double C_fluid = p * rho_w[i] * cp_w[i];
        C_eff[i] = C_rock + C_fluid;
        L_eff[i] = (1.0 - p) * lam_r[i] + p * lam_w[i];
    }
}