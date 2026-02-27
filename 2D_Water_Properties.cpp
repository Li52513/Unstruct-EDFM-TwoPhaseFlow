/**
 * @file 2D_Water_Properties.cpp
 * @brief 2D 水相流体物理属性管理子模块实现
 */

#include "2D_Water_Properties.h"

using namespace PhysicalProperties_string_op;

WaterProperties_2D::WaterProperties_2D(FieldManager_2D& fieldMgr,
    const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig)
    : fieldMgr_(fieldMgr), pConfig_(pConfig), tConfig_(tConfig)
{
}

void WaterProperties_2D::UpdateMatrix_Constant(const WaterProperties& params)
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

void WaterProperties_2D::UpdateFracture_Constant(const WaterProperties& params)
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

void WaterProperties_2D::UpdateMatrix_IAPWS()
{
    Water wTags;
    auto P = fieldMgr_.getMatrixScalar(pConfig_.pressure_field);
    auto T = fieldMgr_.getMatrixScalar(tConfig_.temperatue_field);
    if (!P || !T) return;

    auto rho = fieldMgr_.getOrCreateMatrixScalar(wTags.rho_tag);
    auto mu = fieldMgr_.getOrCreateMatrixScalar(wTags.mu_tag);
    auto cp = fieldMgr_.getOrCreateMatrixScalar(wTags.cp_tag);
    auto cv = fieldMgr_.getOrCreateMatrixScalar(wTags.cv_tag);
    auto h = fieldMgr_.getOrCreateMatrixScalar(wTags.h_tag);
    auto k = fieldMgr_.getOrCreateMatrixScalar(wTags.k_tag);

    int n = static_cast<int>(P->data.size());
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        auto props = WaterPropertyTable::instance().getProperties(P->data[i], T->data[i]);
        rho->data[i] = props.rho;
        mu->data[i] = props.mu;
        cp->data[i] = props.cp;
        cv->data[i] = props.cv;
        h->data[i] = props.h;
        k->data[i] = props.k;
    }
}

void WaterProperties_2D::UpdateFracture_IAPWS()
{
    Water wTags;
    auto P = fieldMgr_.getFractureScalar(pConfig_.pressure_field);
    auto T = fieldMgr_.getFractureScalar(tConfig_.temperatue_field);
    if (!P || !T) return;

    auto rho = fieldMgr_.getOrCreateFractureScalar(wTags.rho_tag);
    auto mu = fieldMgr_.getOrCreateFractureScalar(wTags.mu_tag);
    auto cp = fieldMgr_.getOrCreateFractureScalar(wTags.cp_tag);
    auto cv = fieldMgr_.getOrCreateFractureScalar(wTags.cv_tag);
    auto h = fieldMgr_.getOrCreateFractureScalar(wTags.h_tag);
    auto k = fieldMgr_.getOrCreateFractureScalar(wTags.k_tag);

    int n = static_cast<int>(P->data.size());
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        auto props = WaterPropertyTable::instance().getProperties(P->data[i], T->data[i]);
        rho->data[i] = props.rho;
        mu->data[i] = props.mu;
        cp->data[i] = props.cp;
        cv->data[i] = props.cv;
        h->data[i] = props.h;
        k->data[i] = props.k;
    }
}

void WaterProperties_2D::CalculateCompressibilityMatrix()
{
    Water wTags;
    auto P = fieldMgr_.getMatrixScalar(pConfig_.pressure_field);
    auto T = fieldMgr_.getMatrixScalar(tConfig_.temperatue_field);
    auto rho = fieldMgr_.getMatrixScalar(wTags.rho_tag);
    if (!P || !T || !rho) return;

    auto c_w = fieldMgr_.getOrCreateMatrixScalar(wTags.c_w_tag);
    int n = static_cast<int>(P->data.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p_val = P->data[i];
        double t_val = T->data[i];
        double rho_val = rho->data[i];

        if (rho_val < 1e-3) {
            c_w->data[i] = 4.5e-10; // 物理兜底：标准水压缩系数
            continue;
        }

        double dP = std::max(10.0, std::min(std::abs(p_val) * 1e-6, 1000.0));
        try {
            double r_plus = WaterPropertyTable::instance().getProperties(p_val + dP, t_val).rho;
            double r_minus = WaterPropertyTable::instance().getProperties(p_val - dP, t_val).rho;
            c_w->data[i] = ((r_plus - r_minus) / (2.0 * dP)) / rho_val;
        }
        catch (...) {
            c_w->data[i] = 4.5e-10; // 越界兜底
        }
    }
}

void WaterProperties_2D::CalculateCompressibilityFracture()
{
    Water wTags;
    auto P = fieldMgr_.getFractureScalar(pConfig_.pressure_field);
    auto T = fieldMgr_.getFractureScalar(tConfig_.temperatue_field);
    auto rho = fieldMgr_.getFractureScalar(wTags.rho_tag);
    if (!P || !T || !rho) return;

    auto c_w = fieldMgr_.getOrCreateFractureScalar(wTags.c_w_tag);
    int n = static_cast<int>(P->data.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p_val = P->data[i];
        double t_val = T->data[i];
        double rho_val = rho->data[i];

        if (rho_val < 1e-3) {
            c_w->data[i] = 4.5e-10;
            continue;
        }

        double dP = std::max(10.0, std::min(std::abs(p_val) * 1e-6, 1000.0));
        try {
            double r_plus = WaterPropertyTable::instance().getProperties(p_val + dP, t_val).rho;
            double r_minus = WaterPropertyTable::instance().getProperties(p_val - dP, t_val).rho;
            c_w->data[i] = ((r_plus - r_minus) / (2.0 * dP)) / rho_val;
        }
        catch (...) {
            c_w->data[i] = 4.5e-10;
        }
    }
}

void WaterProperties_2D::UpdateEffectiveThermalPropertiesMatrix()
{
    EffectiveProps effTags;
    Rock rTags;
    Water wTags;

    auto phi = fieldMgr_.getMatrixScalar(rTags.phi_tag);
    auto rho_r = fieldMgr_.getMatrixScalar(rTags.rho_tag);
    auto cp_r = fieldMgr_.getMatrixScalar(rTags.cp_tag);
    auto k_r = fieldMgr_.getMatrixScalar(rTags.lambda_tag);
    auto rho_w = fieldMgr_.getMatrixScalar(wTags.rho_tag);
    auto cp_w = fieldMgr_.getMatrixScalar(wTags.cp_tag);
    auto k_w = fieldMgr_.getMatrixScalar(wTags.k_tag);

    if (!phi || !rho_w || !rho_r || !cp_r || !cp_w || !k_r || !k_w) return;

    auto C_eff = fieldMgr_.getOrCreateMatrixScalar(effTags.C_eff_tag);
    auto L_eff = fieldMgr_.getOrCreateMatrixScalar(effTags.lambda_eff_tag);
    int n = static_cast<int>(phi->data.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = phi->data[i];
        C_eff->data[i] = (1.0 - p) * rho_r->data[i] * cp_r->data[i] + p * rho_w->data[i] * cp_w->data[i];
        L_eff->data[i] = (1.0 - p) * k_r->data[i] + p * k_w->data[i];
    }
}

void WaterProperties_2D::UpdateEffectiveThermalPropertiesFracture()
{
    EffectiveProps effTags;
    Fracture_string fTags;
    Water wTags;

    auto phi = fieldMgr_.getFractureScalar(fTags.phi_tag);
    auto rho_r = fieldMgr_.getFractureScalar(fTags.rho_tag);
    auto cp_r = fieldMgr_.getFractureScalar(fTags.cp_tag);
    auto k_r = fieldMgr_.getFractureScalar(fTags.lambda_tag);
    auto rho_w = fieldMgr_.getFractureScalar(wTags.rho_tag);
    auto cp_w = fieldMgr_.getFractureScalar(wTags.cp_tag);
    auto k_w = fieldMgr_.getFractureScalar(wTags.k_tag);

    if (!phi || !rho_w || !rho_r || !cp_r || !cp_w || !k_r || !k_w) return;

    auto C_eff = fieldMgr_.getOrCreateFractureScalar(effTags.C_eff_tag);
    auto L_eff = fieldMgr_.getOrCreateFractureScalar(effTags.lambda_eff_tag);
    int n = static_cast<int>(phi->data.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = phi->data[i];
        C_eff->data[i] = (1.0 - p) * rho_r->data[i] * cp_r->data[i] + p * rho_w->data[i] * cp_w->data[i];
        L_eff->data[i] = (1.0 - p) * k_r->data[i] + p * k_w->data[i];
    }
}