#pragma once

#include "CaseCommon_Artifacts.h"

#include <string>
#include <vector>

namespace CaseCommon {

enum class CaseDimension {
    Dim2D,
    Dim3D
};

enum class CaseEquationMode {
    N1,
    N2,
    N3
};

enum class CasePropertyMode {
    ConstPP,
    VaryPP
};

enum class CaseFractureMode {
    NoFrac,
    SingleFrac,
    ComplexFrac
};

enum class CaseWellMode {
    NoWell,
    InjProd
};

struct CaseMetadata {
    std::string case_code;
    std::string dispatcher_key;
    std::string case_slug;
    std::string description;
    std::string reference_mode;
    std::string implementation_status;
    std::string output_root;
    std::string well_control_policy;
    std::string injection_fluid;
    std::string thermal_injection_policy;
    CaseDimension dimension = CaseDimension::Dim2D;
    CaseEquationMode equation_mode = CaseEquationMode::N1;
    CasePropertyMode property_mode = CasePropertyMode::ConstPP;
    CaseFractureMode fracture_mode = CaseFractureMode::NoFrac;
    CaseWellMode well_mode = CaseWellMode::NoWell;
};

using CaseStageRunner = int (*)(CaseStage stage);

struct CaseCatalogEntry {
    CaseMetadata metadata;
    CaseStageRunner run_stage = nullptr;
};

const std::vector<CaseCatalogEntry>& GetCaseCatalog();
const CaseCatalogEntry* FindCaseCatalogEntry(const std::string& selector);
int RunCatalogCase(const CaseCatalogEntry& entry, CaseStage stage);
void ValidateCaseCatalogOrThrow();

const char* ToString(CaseDimension value);
const char* ToString(CaseEquationMode value);
const char* ToString(CasePropertyMode value);
const char* ToString(CaseFractureMode value);
const char* ToString(CaseWellMode value);
const FamilyAcceptancePolicy& GetFamilyAcceptancePolicy(CaseEquationMode value);

} // namespace CaseCommon
