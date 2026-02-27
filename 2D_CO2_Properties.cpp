/**
 * @file 2D_CO2_Properties.cpp
 * @brief 2D CO2流体物理属性管理子模块实现
 */

#include "2D_CO2_Properties.h"

using namespace PhysicalProperties_string_op;

CO2Properties_2D::CO2Properties_2D(FieldManager_2D& fieldMgr,
    const PressureEquation_String& pConfig,
    const TemperatureEquation_String& tConfig)
    : fieldMgr_(fieldMgr), tConfig_(tConfig)
{
    // [严谨修正] CO2 的物性压力必须使用气相压力 P_g，而非全局的 P_w
    TwoPhaseState_String tpState;
    pg_name_ = tpState.p_g_field;
}

void CO2Properties_2D::UpdateMatrix_Constant(const CO2Properties& params)
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

void CO2Properties_2D::UpdateFracture_Constant(const CO2Properties& params)
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

void CO2Properties_2D::UpdateMatrix_SpanWagner()
{
    CO2 cTags;
    auto P = fieldMgr_.getMatrixScalar(pg_name_);
    auto T = fieldMgr_.getMatrixScalar(tConfig_.temperatue_field);
    if (!P || !T) return;

    auto rho = fieldMgr_.getOrCreateMatrixScalar(cTags.rho_tag);
    auto mu = fieldMgr_.getOrCreateMatrixScalar(cTags.mu_tag);
    auto cp = fieldMgr_.getOrCreateMatrixScalar(cTags.cp_tag);
    auto cv = fieldMgr_.getOrCreateMatrixScalar(cTags.cv_tag);
    auto h = fieldMgr_.getOrCreateMatrixScalar(cTags.h_tag);
    auto k = fieldMgr_.getOrCreateMatrixScalar(cTags.k_tag);

    int n = static_cast<int>(P->data.size());
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        auto props = CO2PropertyTable::instance().getProperties(P->data[i], T->data[i]);
        rho->data[i] = props.rho;
        mu->data[i] = props.mu;
        cp->data[i] = props.cp;
        cv->data[i] = props.cv;
        h->data[i] = props.h;
        k->data[i] = props.k;
    }
}

void CO2Properties_2D::UpdateFracture_SpanWagner()
{
    CO2 cTags;
    auto P = fieldMgr_.getFractureScalar(pg_name_);
    auto T = fieldMgr_.getFractureScalar(tConfig_.temperatue_field);
    if (!P || !T) return;

    auto rho = fieldMgr_.getOrCreateFractureScalar(cTags.rho_tag);
    auto mu = fieldMgr_.getOrCreateFractureScalar(cTags.mu_tag);
    auto cp = fieldMgr_.getOrCreateFractureScalar(cTags.cp_tag);
    auto cv = fieldMgr_.getOrCreateFractureScalar(cTags.cv_tag);
    auto h = fieldMgr_.getOrCreateFractureScalar(cTags.h_tag);
    auto k = fieldMgr_.getOrCreateFractureScalar(cTags.k_tag);

    int n = static_cast<int>(P->data.size());
#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        auto props = CO2PropertyTable::instance().getProperties(P->data[i], T->data[i]);
        rho->data[i] = props.rho;
        mu->data[i] = props.mu;
        cp->data[i] = props.cp;
        cv->data[i] = props.cv;
        h->data[i] = props.h;
        k->data[i] = props.k;
    }
}

void CO2Properties_2D::CalculateCompressibilityMatrix()
{
    CO2 cTags;
    auto P = fieldMgr_.getMatrixScalar(pg_name_);
    auto T = fieldMgr_.getMatrixScalar(tConfig_.temperatue_field);
    auto rho = fieldMgr_.getMatrixScalar(cTags.rho_tag);
    if (!P || !T || !rho) return;

    auto c_g = fieldMgr_.getOrCreateMatrixScalar(cTags.c_g_tag);
    int n = static_cast<int>(P->data.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p_val = P->data[i];
        double t_val = T->data[i];
        double rho_val = rho->data[i];

        if (rho_val < 1e-3) {
            c_g->data[i] = 1e-7; // 物理兜底
            continue;
        }

        double dP = std::max(10.0, std::min(std::abs(p_val) * 1e-6, 1000.0));
        try {
            double r_plus = CO2PropertyTable::instance().getProperties(p_val + dP, t_val).rho;
            double r_minus = CO2PropertyTable::instance().getProperties(p_val - dP, t_val).rho;
            c_g->data[i] = ((r_plus - r_minus) / (2.0 * dP)) / rho_val;
        }
        catch (...) {
            c_g->data[i] = 1e-7;
        }
    }
}

void CO2Properties_2D::CalculateCompressibilityFracture()
{
    CO2 cTags;
    auto P = fieldMgr_.getFractureScalar(pg_name_);
    auto T = fieldMgr_.getFractureScalar(tConfig_.temperatue_field);
    auto rho = fieldMgr_.getFractureScalar(cTags.rho_tag);
    if (!P || !T || !rho) return;

    auto c_g = fieldMgr_.getOrCreateFractureScalar(cTags.c_g_tag);
    int n = static_cast<int>(P->data.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p_val = P->data[i];
        double t_val = T->data[i];
        double rho_val = rho->data[i];

        if (rho_val < 1e-3) {
            c_g->data[i] = 1e-7;
            continue;
        }

        double dP = std::max(10.0, std::min(std::abs(p_val) * 1e-6, 1000.0));
        try {
            double r_plus = CO2PropertyTable::instance().getProperties(p_val + dP, t_val).rho;
            double r_minus = CO2PropertyTable::instance().getProperties(p_val - dP, t_val).rho;
            c_g->data[i] = ((r_plus - r_minus) / (2.0 * dP)) / rho_val;
        }
        catch (...) {
            c_g->data[i] = 1e-7;
        }
    }
}

void CO2Properties_2D::UpdateEffectiveThermalPropertiesMatrix()
{
    EffectiveProps effTags;
    Rock rTags;
    CO2 cTags;

    auto phi = fieldMgr_.getMatrixScalar(rTags.phi_tag);
    auto rho_r = fieldMgr_.getMatrixScalar(rTags.rho_tag);
    auto cp_r = fieldMgr_.getMatrixScalar(rTags.cp_tag);
    auto k_r = fieldMgr_.getMatrixScalar(rTags.lambda_tag);
    auto rho_g = fieldMgr_.getMatrixScalar(cTags.rho_tag);
    auto cp_g = fieldMgr_.getMatrixScalar(cTags.cp_tag);
    auto k_g = fieldMgr_.getMatrixScalar(cTags.k_tag);

    if (!phi || !rho_g || !rho_r || !cp_r || !cp_g || !k_r || !k_g) return;

    auto C_eff = fieldMgr_.getOrCreateMatrixScalar(effTags.C_eff_tag);
    auto L_eff = fieldMgr_.getOrCreateMatrixScalar(effTags.lambda_eff_tag);
    int n = static_cast<int>(phi->data.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = phi->data[i];
        C_eff->data[i] = (1.0 - p) * rho_r->data[i] * cp_r->data[i] + p * rho_g->data[i] * cp_g->data[i];
        L_eff->data[i] = (1.0 - p) * k_r->data[i] + p * k_g->data[i];
    }
}

void CO2Properties_2D::UpdateEffectiveThermalPropertiesFracture()
{
    EffectiveProps effTags;
    Fracture_string fTags;
    CO2 cTags;

    auto phi = fieldMgr_.getFractureScalar(fTags.phi_tag);
    auto rho_r = fieldMgr_.getFractureScalar(fTags.rho_tag);
    auto cp_r = fieldMgr_.getFractureScalar(fTags.cp_tag);
    auto k_r = fieldMgr_.getFractureScalar(fTags.lambda_tag);
    auto rho_g = fieldMgr_.getFractureScalar(cTags.rho_tag);
    auto cp_g = fieldMgr_.getFractureScalar(cTags.cp_tag);
    auto k_g = fieldMgr_.getFractureScalar(cTags.k_tag);

    if (!phi || !rho_g || !rho_r || !cp_r || !cp_g || !k_r || !k_g) return;

    auto C_eff = fieldMgr_.getOrCreateFractureScalar(effTags.C_eff_tag);
    auto L_eff = fieldMgr_.getOrCreateFractureScalar(effTags.lambda_eff_tag);
    int n = static_cast<int>(phi->data.size());

#pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        double p = phi->data[i];
        C_eff->data[i] = (1.0 - p) * rho_r->data[i] * cp_r->data[i] + p * rho_g->data[i] * cp_g->data[i];
        L_eff->data[i] = (1.0 - p) * k_r->data[i] + p * k_g->data[i];
    }
}