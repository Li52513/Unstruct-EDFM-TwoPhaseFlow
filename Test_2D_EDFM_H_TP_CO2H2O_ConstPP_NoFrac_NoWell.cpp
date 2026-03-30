#include "Test_2D_EDFM_H_TP_CO2H2O_ConstPP_NoFrac_NoWell.h"

#include "CaseCommon_Skeleton.h"

namespace Test_H_TP_CO2H2O_ConstPP_NoFrac {
namespace {

CaseCommon::SkeletonTemplateContext BuildContext() {
    CaseCommon::SkeletonTemplateContext ctx;
    ctx.case_code = "C1";
    ctx.case_slug = "2d_n3_co2h2o_constpp_nofrac_nowell";
    ctx.title = "2D N=3 CO2/H2O const-property no-fracture no-well";
    ctx.notes.push_back("Reserved as the first two-phase thin-template baseline in the A1-F12 system.");
    ctx.notes.push_back("Reference solving remains weakly coupled: external Java/PS1 COMSOL automation consumes reference specs from engineering/reference.");
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
        "Two-phase 2D no-well thin-template assembly has not been ported into the new staged contract yet.");
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
        "Validation-only replay requires the C1 engineering/reference readers, which are still pending.");
}

void RunFullWorkflow() {
    const CaseCommon::SkeletonTemplateContext ctx = BuildContext();
    CaseCommon::PrepareSkeletonArtifacts(ctx, CaseCommon::CaseStage::FullWorkflow);
    CaseCommon::ThrowStageNotImplemented(
        ctx,
        CaseCommon::CaseStage::FullWorkflow,
        "C1 is scaffolded but not yet connected to a production N=3 solver template.");
}

} // namespace Test_H_TP_CO2H2O_ConstPP_NoFrac
