#include "CaseCommon_Catalog.h"

#include "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_InjProd.h"
#include "Test_2D_EDFM_H_CO2_ConstPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.h"
#include "Test_2D_EDFM_H_CO2_ConstPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_VaryPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_VaryPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_VaryPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_InjProd.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_VaryPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_VaryPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_VaryPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_InjProd.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_VaryPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_VaryPP_ComplexFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_VaryPP_SingleFrac_NoWell.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <unordered_set>

namespace CaseCommon {
namespace {

bool IsAllowedToken(const std::string& value,
                    std::initializer_list<const char*> allowed) {
    for (const char* token : allowed) {
        if (value == token) {
            return true;
        }
    }
    return false;
}

std::string NormalizeKey(const std::string& text) {
    std::string out = text;
    std::transform(out.begin(), out.end(), out.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return out;
}

int RunA1(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_CO2_ConstPP::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_CO2_ConstPP::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_CO2_ConstPP::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_CO2_ConstPP::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for A1.");
    }
}

int RunA2(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_CO2_VaryPP::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_CO2_VaryPP::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_CO2_VaryPP::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_CO2_VaryPP::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for A2.");
    }
}

int RunA3(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_CO2_ConstPP_SingleFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_CO2_ConstPP_SingleFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_CO2_ConstPP_SingleFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_CO2_ConstPP_SingleFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for A3.");
    }
}

int RunA4(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_CO2_VaryPP_SingleFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_CO2_VaryPP_SingleFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_CO2_VaryPP_SingleFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_CO2_VaryPP_SingleFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for A4.");
    }
}

int RunA5(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_CO2_ConstPP_ComplexFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_CO2_ConstPP_ComplexFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_CO2_ConstPP_ComplexFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_CO2_ConstPP_ComplexFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for A5.");
    }
}

int RunA6(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_CO2_VaryPP_ComplexFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_CO2_VaryPP_ComplexFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_CO2_VaryPP_ComplexFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_CO2_VaryPP_ComplexFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for A6.");
    }
}

int RunA7(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_CO2_ConstPP_NoFrac_InjProd::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_CO2_ConstPP_NoFrac_InjProd::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_CO2_ConstPP_NoFrac_InjProd::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_CO2_ConstPP_NoFrac_InjProd::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for A7.");
    }
}

int RunB1(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_T_CO2_ConstPP_NoFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_T_CO2_ConstPP_NoFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_T_CO2_ConstPP_NoFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_T_CO2_ConstPP_NoFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for B1.");
    }
}

int RunB2(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_T_CO2_VaryPP_NoFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_T_CO2_VaryPP_NoFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_T_CO2_VaryPP_NoFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_T_CO2_VaryPP_NoFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for B2.");
    }
}

int RunB3(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_T_CO2_ConstPP_SingleFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_T_CO2_ConstPP_SingleFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_T_CO2_ConstPP_SingleFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_T_CO2_ConstPP_SingleFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for B3.");
    }
}

int RunB4(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_T_CO2_VaryPP_SingleFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_T_CO2_VaryPP_SingleFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_T_CO2_VaryPP_SingleFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_T_CO2_VaryPP_SingleFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for B4.");
    }
}

int RunB5(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_T_CO2_ConstPP_ComplexFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_T_CO2_ConstPP_ComplexFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_T_CO2_ConstPP_ComplexFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_T_CO2_ConstPP_ComplexFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for B5.");
    }
}

int RunB6(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_T_CO2_VaryPP_ComplexFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_T_CO2_VaryPP_ComplexFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_T_CO2_VaryPP_ComplexFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_T_CO2_VaryPP_ComplexFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for B6.");
    }
}

int RunB7(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_T_CO2_ConstPP_NoFrac_InjProd::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_T_CO2_ConstPP_NoFrac_InjProd::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_T_CO2_ConstPP_NoFrac_InjProd::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_T_CO2_ConstPP_NoFrac_InjProd::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for B7.");
    }
}

int RunC1(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_TP_CO2H2O_ConstPP_NoFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_TP_CO2H2O_ConstPP_NoFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_TP_CO2H2O_ConstPP_NoFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_TP_CO2H2O_ConstPP_NoFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for C1.");
    }
}

int RunC2(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_TP_CO2H2O_VaryPP_NoFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_TP_CO2H2O_VaryPP_NoFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_TP_CO2H2O_VaryPP_NoFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_TP_CO2H2O_VaryPP_NoFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for C2.");
    }
}

int RunC3(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_TP_CO2H2O_ConstPP_SingleFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_TP_CO2H2O_ConstPP_SingleFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_TP_CO2H2O_ConstPP_SingleFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_TP_CO2H2O_ConstPP_SingleFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for C3.");
    }
}

int RunC4(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_TP_CO2H2O_VaryPP_SingleFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_TP_CO2H2O_VaryPP_SingleFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_TP_CO2H2O_VaryPP_SingleFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_TP_CO2H2O_VaryPP_SingleFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for C4.");
    }
}

int RunC5(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_TP_CO2H2O_ConstPP_ComplexFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_TP_CO2H2O_ConstPP_ComplexFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_TP_CO2H2O_ConstPP_ComplexFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_TP_CO2H2O_ConstPP_ComplexFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for C5.");
    }
}

int RunC6(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_TP_CO2H2O_VaryPP_ComplexFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_TP_CO2H2O_VaryPP_ComplexFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_TP_CO2H2O_VaryPP_ComplexFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_TP_CO2H2O_VaryPP_ComplexFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for C6.");
    }
}

int RunC7(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_TP_CO2H2O_ConstPP_NoFrac_InjProd::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_TP_CO2H2O_ConstPP_NoFrac_InjProd::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_TP_CO2H2O_ConstPP_NoFrac_InjProd::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_TP_CO2H2O_ConstPP_NoFrac_InjProd::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for C7.");
    }
}

CaseDimension CaseDimensionFromFamily(char family) {
    return (family == 'A' || family == 'B' || family == 'C')
        ? CaseDimension::Dim2D
        : CaseDimension::Dim3D;
}

CaseEquationMode CaseEquationFromFamily(char family) {
    if (family == 'A' || family == 'D') return CaseEquationMode::N1;
    if (family == 'B' || family == 'E') return CaseEquationMode::N2;
    return CaseEquationMode::N3;
}

CasePropertyMode CasePropertyFromIndex(int number) {
    const int local = ((number - 1) % 6) + 1;
    return (local % 2 == 1) ? CasePropertyMode::ConstPP : CasePropertyMode::VaryPP;
}

CaseFractureMode CaseFractureFromIndex(int number) {
    const int local = ((number - 1) % 6) + 1;
    if (local <= 2) return CaseFractureMode::NoFrac;
    if (local <= 4) return CaseFractureMode::SingleFrac;
    return CaseFractureMode::ComplexFrac;
}

CaseWellMode CaseWellFromIndex(int number) {
    return (number <= 6) ? CaseWellMode::NoWell : CaseWellMode::InjProd;
}

std::string FluidToken(CaseEquationMode mode) {
    return (mode == CaseEquationMode::N3) ? "co2h2o" : "co2";
}

std::string BuildDispatcherKey(CaseDimension dim,
                               CaseEquationMode eq,
                               CasePropertyMode prop,
                               CaseFractureMode frac,
                               CaseWellMode well) {
    return std::string(ToString(dim)) + "_" + std::string(ToString(eq)) + "_" +
        FluidToken(eq) + "_" + ToString(prop) + "_" + ToString(frac) + "_" + ToString(well);
}

std::string BuildDescription(const CaseMetadata& meta) {
    std::string description = std::string(ToString(meta.dimension)) + " ";
    description += std::string(ToString(meta.equation_mode)) + " ";
    description += (meta.equation_mode == CaseEquationMode::N3)
        ? "CO2/H2O"
        : "CO2";
    description += " ";
    description += ToString(meta.property_mode);
    description += " ";
    description += ToString(meta.fracture_mode);
    description += " ";
    description += ToString(meta.well_mode);
    return description;
}

std::string BuildReferenceMode(const std::string& caseCode,
                               CaseEquationMode eq,
                               CaseWellMode well) {
    if (caseCode == "A1") return "analytical";
    if (caseCode == "B1") return "hybrid_analytical_comsol";
    if (eq == CaseEquationMode::N1 && well == CaseWellMode::NoWell) return "comsol";
    return "comsol";
}

std::string BuildOutputRoot() {
    const std::filesystem::path source_root = std::filesystem::path(__FILE__).parent_path();
    std::error_code ec;
    const std::filesystem::path normalized_root =
        std::filesystem::weakly_canonical(source_root, ec);
    const std::filesystem::path output_root =
        (ec ? source_root : normalized_root) / "Test" / "Transient" / "A1_F12_TemplateSystem";
    return output_root.lexically_normal().string();
}

std::string BuildWellControlPolicy(CaseWellMode well) {
    return (well == CaseWellMode::NoWell)
        ? "none"
        : "injector_rate_producer_bhp";
}

std::string BuildInjectionFluid(CaseWellMode well) {
    return (well == CaseWellMode::NoWell)
        ? "none"
        : "co2";
}

std::string BuildThermalInjectionPolicy(CaseEquationMode eq, CaseWellMode well) {
    if (well == CaseWellMode::NoWell) {
        return "none";
    }
    return (eq == CaseEquationMode::N1)
        ? "none"
        : "cold_injection";
}

FamilyAcceptancePolicy BuildPressureDiffusionAcceptancePolicy() {
    FamilyAcceptancePolicy policy;
    policy.family_group_key = "A/D";
    policy.required_profile_variables = {"pressure"};
    policy.required_monitor_variables = {"pressure"};
    policy.required_well_series_when_present = {"well_bhp", "well_rate"};
    return policy;
}

FamilyAcceptancePolicy BuildSinglePhaseThermalAcceptancePolicy() {
    FamilyAcceptancePolicy policy;
    policy.family_group_key = "B/E";
    policy.required_profile_variables = {"pressure", "temperature"};
    policy.required_monitor_variables = {"pressure", "temperature"};
    policy.required_well_series_when_present = {
        "well_bhp",
        "well_rate",
        "production_temperature"
    };
    return policy;
}

FamilyAcceptancePolicy BuildTwoPhaseThermalAcceptancePolicy() {
    FamilyAcceptancePolicy policy;
    policy.family_group_key = "C/F";
    policy.required_profile_variables = {"pressure", "temperature", "co2_saturation"};
    policy.required_monitor_variables = {"pressure", "temperature", "co2_saturation"};
    policy.required_well_series_when_present = {
        "well_bhp",
        "well_rate",
        "production_temperature",
        "phase_fraction",
        "water_cut",
        "co2_production_rate"
    };
    return policy;
}

struct BindingInfo {
    CaseStageRunner runner = nullptr;
    const char* implementation_status = "planned";
};

void ValidateMetadataContractOrThrow(const CaseCatalogEntry& entry) {
    const CaseMetadata& meta = entry.metadata;
    if (meta.case_code.empty()) {
        throw std::runtime_error("[CaseCatalog] empty case_code.");
    }
    if (meta.dispatcher_key.empty()) {
        throw std::runtime_error("[CaseCatalog] empty dispatcher_key for " + meta.case_code);
    }
    if (meta.case_slug.empty()) {
        throw std::runtime_error("[CaseCatalog] empty case_slug for " + meta.case_code);
    }
    if (meta.description.empty()) {
        throw std::runtime_error("[CaseCatalog] empty description for " + meta.case_code);
    }
    if (meta.output_root.empty()) {
        throw std::runtime_error("[CaseCatalog] empty output_root for " + meta.case_code);
    }
    if (!IsAllowedToken(meta.reference_mode, {"analytical", "hybrid_analytical_comsol", "comsol"})) {
        throw std::runtime_error("[CaseCatalog] invalid reference_mode for " + meta.case_code +
                                 ": " + meta.reference_mode);
    }
    if (!IsAllowedToken(meta.implementation_status, {"planned", "skeleton", "implemented"})) {
        throw std::runtime_error("[CaseCatalog] invalid implementation_status for " + meta.case_code +
                                 ": " + meta.implementation_status);
    }
    if (!IsAllowedToken(meta.well_control_policy, {"none", "injector_rate_producer_bhp"})) {
        throw std::runtime_error("[CaseCatalog] invalid well_control_policy for " + meta.case_code +
                                 ": " + meta.well_control_policy);
    }
    if (!IsAllowedToken(meta.injection_fluid, {"none", "co2"})) {
        throw std::runtime_error("[CaseCatalog] invalid injection_fluid for " + meta.case_code +
                                 ": " + meta.injection_fluid);
    }
    if (!IsAllowedToken(meta.thermal_injection_policy, {"none", "cold_injection"})) {
        throw std::runtime_error("[CaseCatalog] invalid thermal_injection_policy for " + meta.case_code +
                                 ": " + meta.thermal_injection_policy);
    }
    if (meta.case_slug != meta.dispatcher_key) {
        throw std::runtime_error("[CaseCatalog] case_slug must match dispatcher_key for " + meta.case_code);
    }
    if (meta.well_mode == CaseWellMode::NoWell) {
        if (meta.well_control_policy != "none" ||
            meta.injection_fluid != "none" ||
            meta.thermal_injection_policy != "none") {
            throw std::runtime_error("[CaseCatalog] no-well case must keep none-valued well policies: " +
                                     meta.case_code);
        }
    }
    else {
        if (meta.well_control_policy != "injector_rate_producer_bhp") {
            throw std::runtime_error("[CaseCatalog] inj-prod case must use injector_rate_producer_bhp: " +
                                     meta.case_code);
        }
        if (meta.injection_fluid != "co2") {
            throw std::runtime_error("[CaseCatalog] inj-prod case must inject co2: " + meta.case_code);
        }
        const std::string expected_thermal = (meta.equation_mode == CaseEquationMode::N1)
            ? "none"
            : "cold_injection";
        if (meta.thermal_injection_policy != expected_thermal) {
            throw std::runtime_error("[CaseCatalog] inconsistent thermal_injection_policy for " +
                                     meta.case_code + ": expected " + expected_thermal +
                                     ", got " + meta.thermal_injection_policy);
        }
    }
    if ((entry.run_stage == nullptr) != (meta.implementation_status == "planned")) {
        throw std::runtime_error("[CaseCatalog] implementation_status/run_stage mismatch for " +
                                 meta.case_code);
    }
}

BindingInfo GetBinding(const std::string& caseCode) {
    if (caseCode == "A1") return BindingInfo{&RunA1, "implemented"};
    if (caseCode == "A2") return BindingInfo{&RunA2, "implemented"};
    if (caseCode == "A3") return BindingInfo{&RunA3, "implemented"};
    if (caseCode == "A4") return BindingInfo{&RunA4, "implemented"};
    if (caseCode == "A5") return BindingInfo{&RunA5, "implemented"};
    if (caseCode == "A6") return BindingInfo{&RunA6, "implemented"};
    if (caseCode == "A7") return BindingInfo{&RunA7, "implemented"};
    if (caseCode == "B1") return BindingInfo{&RunB1, "implemented"};
    if (caseCode == "B2") return BindingInfo{&RunB2, "implemented"};
    if (caseCode == "B3") return BindingInfo{&RunB3, "implemented"};
    if (caseCode == "B4") return BindingInfo{&RunB4, "implemented"};
    if (caseCode == "B5") return BindingInfo{&RunB5, "implemented"};
    if (caseCode == "B6") return BindingInfo{&RunB6, "implemented"};
    if (caseCode == "B7") return BindingInfo{&RunB7, "implemented"};
    if (caseCode == "C1") return BindingInfo{&RunC1, "implemented"};
    if (caseCode == "C2") return BindingInfo{&RunC2, "implemented"};
    if (caseCode == "C3") return BindingInfo{&RunC3, "implemented"};
    if (caseCode == "C4") return BindingInfo{&RunC4, "implemented"};
    if (caseCode == "C5") return BindingInfo{&RunC5, "implemented"};
    if (caseCode == "C6") return BindingInfo{&RunC6, "implemented"};
    if (caseCode == "C7") return BindingInfo{&RunC7, "implemented"};
    return BindingInfo{};
}

std::vector<CaseCatalogEntry> BuildCatalog() {
    std::vector<CaseCatalogEntry> entries;
    entries.reserve(72);

    for (char family = 'A'; family <= 'F'; ++family) {
        for (int number = 1; number <= 12; ++number) {
            CaseCatalogEntry entry;
            entry.metadata.case_code = std::string(1, family) + std::to_string(number);
            entry.metadata.dimension = CaseDimensionFromFamily(family);
            entry.metadata.equation_mode = CaseEquationFromFamily(family);
            entry.metadata.property_mode = CasePropertyFromIndex(number);
            entry.metadata.fracture_mode = CaseFractureFromIndex(number);
            entry.metadata.well_mode = CaseWellFromIndex(number);
            entry.metadata.dispatcher_key = BuildDispatcherKey(
                entry.metadata.dimension,
                entry.metadata.equation_mode,
                entry.metadata.property_mode,
                entry.metadata.fracture_mode,
                entry.metadata.well_mode);
            entry.metadata.case_slug = entry.metadata.dispatcher_key;
            entry.metadata.description = BuildDescription(entry.metadata);
            entry.metadata.reference_mode = BuildReferenceMode(
                entry.metadata.case_code,
                entry.metadata.equation_mode,
                entry.metadata.well_mode);
            entry.metadata.output_root = BuildOutputRoot();
            entry.metadata.well_control_policy = BuildWellControlPolicy(entry.metadata.well_mode);
            entry.metadata.injection_fluid = BuildInjectionFluid(entry.metadata.well_mode);
            entry.metadata.thermal_injection_policy = BuildThermalInjectionPolicy(
                entry.metadata.equation_mode,
                entry.metadata.well_mode);

            const BindingInfo binding = GetBinding(entry.metadata.case_code);
            entry.metadata.implementation_status = binding.implementation_status;
            entry.run_stage = binding.runner;
            entries.push_back(entry);
        }
    }

    return entries;
}

} // namespace

const std::vector<CaseCatalogEntry>& GetCaseCatalog() {
    static const std::vector<CaseCatalogEntry> catalog = BuildCatalog();
    return catalog;
}

const CaseCatalogEntry* FindCaseCatalogEntry(const std::string& selector) {
    const std::string key = NormalizeKey(selector);
    const auto& catalog = GetCaseCatalog();
    for (const auto& entry : catalog) {
        if (key == NormalizeKey(entry.metadata.case_code) ||
            key == NormalizeKey(entry.metadata.dispatcher_key) ||
            key == NormalizeKey(entry.metadata.case_slug)) {
            return &entry;
        }
    }
    return nullptr;
}

int RunCatalogCase(const CaseCatalogEntry& entry, CaseStage stage) {
    if (!entry.run_stage) {
        throw std::runtime_error(
            "[CaseCatalog] case " + entry.metadata.case_code +
            " is registered in the A1-F12 catalog but not implemented yet.");
    }
    return entry.run_stage(stage);
}

void ValidateCaseCatalogOrThrow() {
    const auto& catalog = GetCaseCatalog();
    if (catalog.size() != 72u) {
        throw std::runtime_error("[CaseCatalog] expected 72 catalog entries.");
    }

    std::unordered_set<std::string> codes;
    std::unordered_set<std::string> keys;
    for (const auto& entry : catalog) {
        const std::string code = NormalizeKey(entry.metadata.case_code);
        const std::string key = NormalizeKey(entry.metadata.dispatcher_key);
        if (!codes.insert(code).second) {
            throw std::runtime_error("[CaseCatalog] duplicate case code: " + entry.metadata.case_code);
        }
        if (!keys.insert(key).second) {
            throw std::runtime_error("[CaseCatalog] duplicate dispatcher key: " + entry.metadata.dispatcher_key);
        }
        ValidateMetadataContractOrThrow(entry);
    }
}

const char* ToString(CaseDimension value) {
    switch (value) {
    case CaseDimension::Dim2D: return "2d";
    case CaseDimension::Dim3D: return "3d";
    default: return "unknown_dim";
    }
}

const char* ToString(CaseEquationMode value) {
    switch (value) {
    case CaseEquationMode::N1: return "n1";
    case CaseEquationMode::N2: return "n2";
    case CaseEquationMode::N3: return "n3";
    default: return "unknown_n";
    }
}

const char* ToString(CasePropertyMode value) {
    switch (value) {
    case CasePropertyMode::ConstPP: return "constpp";
    case CasePropertyMode::VaryPP: return "varypp";
    default: return "unknown_property";
    }
}

const char* ToString(CaseFractureMode value) {
    switch (value) {
    case CaseFractureMode::NoFrac: return "nofrac";
    case CaseFractureMode::SingleFrac: return "singlefrac";
    case CaseFractureMode::ComplexFrac: return "complexfrac";
    default: return "unknown_frac";
    }
}

const char* ToString(CaseWellMode value) {
    switch (value) {
    case CaseWellMode::NoWell: return "nowell";
    case CaseWellMode::InjProd: return "injprod";
    default: return "unknown_well";
    }
}

const FamilyAcceptancePolicy& GetFamilyAcceptancePolicy(CaseEquationMode value) {
    static const FamilyAcceptancePolicy pressure_diffusion = BuildPressureDiffusionAcceptancePolicy();
    static const FamilyAcceptancePolicy single_phase_thermal = BuildSinglePhaseThermalAcceptancePolicy();
    static const FamilyAcceptancePolicy two_phase_thermal = BuildTwoPhaseThermalAcceptancePolicy();

    switch (value) {
    case CaseEquationMode::N1: return pressure_diffusion;
    case CaseEquationMode::N2: return single_phase_thermal;
    case CaseEquationMode::N3: return two_phase_thermal;
    default: return pressure_diffusion;
    }
}

} // namespace CaseCommon
