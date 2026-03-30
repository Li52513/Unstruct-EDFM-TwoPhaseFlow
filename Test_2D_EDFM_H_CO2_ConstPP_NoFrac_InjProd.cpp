#include "Test_2D_EDFM_H_CO2_ConstPP_NoFrac_InjProd.h"

#include "CaseCommon_Skeleton.h"

namespace Test_H_CO2_ConstPP_NoFrac_InjProd {
namespace {

CaseCommon::SkeletonTemplateContext BuildContext() {
    CaseCommon::SkeletonTemplateContext ctx;
    ctx.case_code = "A7";
    ctx.case_slug = "2d_n1_co2_constpp_nofrac_injprod";
    ctx.title = "2D N=1 CO2 const-property no-fracture injector-producer";
    ctx.notes.push_back("Default well policy: injector rate control + producer BHP control.");
    ctx.notes.push_back("Default well layout: 1/4L injector and 3/4L producer on the domain midline; matrix completion only.");
    ctx.notes.push_back("Current blocker: the dedicated N=1 AD route still rejects wells inside FIM_TransientEngine/RunGeneric_impl.hpp.");
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
        "A7 requires true N=1 well support; the current pressure-only AD route still throws when wells are present.");
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
        "A7 validation depends on well-enabled engineering outputs that are not available until N=1 well assembly is completed.");
}

void RunFullWorkflow() {
    const CaseCommon::SkeletonTemplateContext ctx = BuildContext();
    CaseCommon::PrepareSkeletonArtifacts(ctx, CaseCommon::CaseStage::FullWorkflow);
    CaseCommon::ThrowStageNotImplemented(
        ctx,
        CaseCommon::CaseStage::FullWorkflow,
        "A7 has been scaffolded with directory and well-policy contracts, but solve/validate execution is still pending.");
}

} // namespace Test_H_CO2_ConstPP_NoFrac_InjProd
