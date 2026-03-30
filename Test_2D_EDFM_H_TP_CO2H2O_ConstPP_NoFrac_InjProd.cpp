#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_InjProd.h"

#include "CaseCommon_Skeleton.h"

namespace Test_H_TP_CO2H2O_ConstPP_NoFrac_InjProd {
namespace {

CaseCommon::SkeletonTemplateContext BuildContext() {
    CaseCommon::SkeletonTemplateContext ctx;
    ctx.case_code = "C7";
    ctx.case_slug = "2d_n3_co2h2o_constpp_nofrac_injprod";
    ctx.title = "2D N=3 CO2/H2O const-property no-fracture injector-producer";
    ctx.notes.push_back("Default phase policy: inject CO2, produce mixed CO2/H2O; matrix completion only.");
    ctx.notes.push_back("Future validation must add producer phase fraction / water cut / CO2 rate curves.");
    ctx.notes.push_back("This skeleton keeps COMSOL integration external and only reserves the artifact contract.");
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
        "C7 requires a two-phase well-enabled thin template that has not been ported yet.");
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
        "C7 validation-only mode is blocked until the two-phase well engineering outputs are standardized.");
}

void RunFullWorkflow() {
    const CaseCommon::SkeletonTemplateContext ctx = BuildContext();
    CaseCommon::PrepareSkeletonArtifacts(ctx, CaseCommon::CaseStage::FullWorkflow);
    CaseCommon::ThrowStageNotImplemented(
        ctx,
        CaseCommon::CaseStage::FullWorkflow,
        "C7 currently exists as a staged skeleton only.");
}

} // namespace Test_H_TP_CO2H2O_ConstPP_NoFrac_InjProd
