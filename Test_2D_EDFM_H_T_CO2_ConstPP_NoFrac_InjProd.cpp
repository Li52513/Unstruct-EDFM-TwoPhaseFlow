#include "Test_2D_EDFM_H_T_CO2_ConstPP_NoFrac_InjProd.h"

#include "CaseCommon_Skeleton.h"

namespace Test_H_T_CO2_ConstPP_NoFrac_InjProd {
namespace {

CaseCommon::SkeletonTemplateContext BuildContext() {
    CaseCommon::SkeletonTemplateContext ctx;
    ctx.case_code = "B7";
    ctx.case_slug = "2d_n2_co2_constpp_nofrac_injprod";
    ctx.title = "2D N=2 CO2 const-property no-fracture injector-producer";
    ctx.notes.push_back("Default thermal well policy: cold CO2 injection under rate control, producer under BHP control.");
    ctx.notes.push_back("Well-level validation must include BHP, actual rate and produced temperature once the template is completed.");
    ctx.notes.push_back("This file reserves the thin-template boundary without embedding COMSOL execution.");
    return ctx;
}

} // namespace

void RunTestCase() {
    RunFullWorkflow();
}

void RunSolveOnly() {
    const CaseCommon::SkeletonTemplateContext ctx = BuildContext();
    CaseCommon::PrepareSkeletonArtifacts(ctx, CaseCommon::CaseStage::SolveOnly);
    CaseCommon::ThrowStageNotImplemented(
        ctx,
        CaseCommon::CaseStage::SolveOnly,
        "B7 well-enabled P-T assembly has not been split into the new staged template yet.");
}

void RunPrepareReference() {
    CaseCommon::PrepareSkeletonArtifacts(BuildContext(), CaseCommon::CaseStage::PrepareReference);
}

void RunValidateOnly() {
    const CaseCommon::SkeletonTemplateContext ctx = BuildContext();
    CaseCommon::PrepareSkeletonArtifacts(ctx, CaseCommon::CaseStage::ValidateOnly);
    CaseCommon::ThrowStageNotImplemented(
        ctx,
        CaseCommon::CaseStage::ValidateOnly,
        "B7 validation requires the future well-aware engineering outputs and COMSOL reference import path.");
}

void RunFullWorkflow() {
    const CaseCommon::SkeletonTemplateContext ctx = BuildContext();
    CaseCommon::PrepareSkeletonArtifacts(ctx, CaseCommon::CaseStage::FullWorkflow);
    CaseCommon::ThrowStageNotImplemented(
        ctx,
        CaseCommon::CaseStage::FullWorkflow,
        "B7 is scaffolded only; the production solver and validation chain are pending.");
}

} // namespace Test_H_T_CO2_ConstPP_NoFrac_InjProd
