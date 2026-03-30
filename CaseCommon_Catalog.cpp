#include "CaseCommon_Catalog.h"

#include "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_InjProd.h"
#include "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_NoWell_.h"
#include "Test_2D_EDFM_H_CO2_ConstPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_VaryPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_CO2_VaryPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_InjProd.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_NoWell.h"
#include "Test_2D_EDFM_H_T_CO2_ConstPP_SingleFrac_NoWell.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_InjProd.h"
#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_NoWell.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <unordered_set>

namespace CaseCommon {
namespace {

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

int RunB3(CaseStage stage) {
    switch (stage) {
    case CaseStage::SolveOnly: Test_H_T_CO2_ConstPP_SingleFrac::RunSolveOnly(); return 0;
    case CaseStage::PrepareReference: Test_H_T_CO2_ConstPP_SingleFrac::RunPrepareReference(); return 0;
    case CaseStage::ValidateOnly: Test_H_T_CO2_ConstPP_SingleFrac::RunValidateOnly(); return 0;
    case CaseStage::FullWorkflow: Test_H_T_CO2_ConstPP_SingleFrac::RunFullWorkflow(); return 0;
    default: throw std::runtime_error("[CaseCatalog] unsupported stage for B3.");
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

struct BindingInfo {
    CaseStageRunner runner = nullptr;
    const char* implementation_status = "planned";
};

BindingInfo GetBinding(const std::string& caseCode) {
    if (caseCode == "A1") return BindingInfo{&RunA1, "implemented"};
    if (caseCode == "A2") return BindingInfo{&RunA2, "implemented"};
    if (caseCode == "A3") return BindingInfo{&RunA3, "implemented"};
    if (caseCode == "A4") return BindingInfo{&RunA4, "implemented"};
    if (caseCode == "A7") return BindingInfo{&RunA7, "skeleton"};
    if (caseCode == "B1") return BindingInfo{&RunB1, "implemented"};
    if (caseCode == "B3") return BindingInfo{&RunB3, "implemented"};
    if (caseCode == "B7") return BindingInfo{&RunB7, "skeleton"};
    if (caseCode == "C1") return BindingInfo{&RunC1, "skeleton"};
    if (caseCode == "C7") return BindingInfo{&RunC7, "skeleton"};
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

} // namespace CaseCommon
