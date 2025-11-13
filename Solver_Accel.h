#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <vector>
#include <unordered_map>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "PhysicalPropertiesManager.h"
#include "BCAdapter.h"
#include "TemperatureBCAdapter.h"
#include "WellConfig.h"
#include "FVM_WellDOF.h"
#include "FVM_SourceTerm_WellHeat.h"

#include "DiffusionCentral.h"
#include "ConvectionUpwind_Flux.h"
#include "ConvectionUpwind.h"
#include "Timeterm_BDF.h"

#include "Solver_AssemblerCOO.h"
#include "Solver_TimeLoopUtils.h"
#include "Solver_PostChecks.h"
#include "LinearSolver_Eigen.h"
#include "FaceSignMask.hpp"

namespace Solver {
namespace Accel {

inline long long makeRCKey(int r, int c)
{
    return (static_cast<long long>(r) << 32) ^ static_cast<unsigned long long>(c);
}

struct ReuseLinearSolver
{
    using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor>;
    using Trip = Eigen::Triplet<double>;

    SpMat Aeq;
    Eigen::VectorXd Dinvsqrt;
    bool useScaling = false;
    bool patternBuilt = false;
    bool haveFactor = false;
    int n = 0;
    size_t nnzUnique = 0;

    std::vector<std::pair<int, int>> entries;
    std::vector<double> valuesBuffer;
    std::unordered_map<long long, size_t> entryLookup;

    Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<double>> solver;
    double lastFill = -1.0;
    double lastDrop = -1.0;

    void reset()
    {
        Aeq.resize(0, 0);
        Dinvsqrt.resize(0);
        useScaling = false;
        patternBuilt = false;
        haveFactor = false;
        n = 0;
        nnzUnique = 0;
        entries.clear();
        valuesBuffer.clear();
        entryLookup.clear();
        lastFill = -1.0;
        lastDrop = -1.0;
    }

    bool buildPattern(const SparseSystemCOO& sys, const LinearSolverOptions& opt, bool report, const char* tag)
    {
        reset();
        if (sys.n <= 0) return false;

        n = sys.n;
        entryLookup.reserve(sys.A.size() * 2 + 1);
        std::vector<double> diag(n, 1.0);
        if (opt.equil) std::fill(diag.begin(), diag.end(), 0.0);

        for (const auto& t : sys.A) {
            long long key = makeRCKey(t.r, t.c);
            auto emplaceRes = entryLookup.emplace(key, entries.size());
            auto it = emplaceRes.first;
            bool inserted = emplaceRes.second;
            if (inserted) {
                entries.emplace_back(t.r, t.c);
                valuesBuffer.push_back(t.v);
            }
            else {
                valuesBuffer[it->second] += t.v;
            }
            if (opt.equil && t.r >= 0 && t.r < n && t.r == t.c) {
                diag[t.r] = std::max(diag[t.r], std::abs(t.v));
            }
        }

        Dinvsqrt = Eigen::VectorXd::Ones(n);
        if (opt.equil) {
            for (int i = 0; i < n; ++i) {
                double d = diag[i];
                if (d <= 1e-30) d = 1.0;
                Dinvsqrt[i] = 1.0 / std::sqrt(d);
            }
            useScaling = true;
        }
        else {
            useScaling = false;
        }

        std::vector<Trip> trips;
        trips.reserve(entries.size());
        for (size_t i = 0; i < entries.size(); ++i) {
            const auto& rc = entries[i];
            double v = valuesBuffer[i];
            if (useScaling) v *= Dinvsqrt[rc.first] * Dinvsqrt[rc.second];
            trips.emplace_back(rc.first, rc.second, v);
        }
        Aeq.resize(n, n);
        Aeq.setFromTriplets(trips.begin(), trips.end());
        Aeq.makeCompressed();

        nnzUnique = entries.size();
        valuesBuffer.assign(entries.size(), 0.0);
        patternBuilt = true;
        haveFactor = false;

        if (report) {
            std::cout << "[ReuseSolver] build pattern (" << tag << ") n=" << n << " nnz=" << nnzUnique << "\n";
        }
        return true;
    }

    bool updateValues(const SparseSystemCOO& sys)
    {
        if (!patternBuilt || sys.n != n) return false;
        if (valuesBuffer.size() != entries.size()) valuesBuffer.assign(entries.size(), 0.0);
        else std::fill(valuesBuffer.begin(), valuesBuffer.end(), 0.0);

        for (const auto& t : sys.A) {
            auto it = entryLookup.find(makeRCKey(t.r, t.c));
            if (it == entryLookup.end()) return false;
            valuesBuffer[it->second] += t.v;
        }

        for (size_t i = 0; i < entries.size(); ++i) {
            const auto& rc = entries[i];
            double v = valuesBuffer[i];
            if (useScaling) v *= Dinvsqrt[rc.first] * Dinvsqrt[rc.second];
            Aeq.coeffRef(rc.first, rc.second) = v;
        }
        Aeq.makeCompressed();
        return true;
    }

    bool ensurePattern(const SparseSystemCOO& sys, const LinearSolverOptions& opt, bool report, const char* tag, bool& rebuilt)
    {
        rebuilt = false;
        if (!patternBuilt || sys.n != n) {
            rebuilt = true;
            return buildPattern(sys, opt, report, tag);
        }
        if (!updateValues(sys)) {
            rebuilt = true;
            return buildPattern(sys, opt, report, tag);
        }
        return true;
    }

    bool ensurePreconditioner(const LinearSolverOptions& opt, bool rebuild, bool report, const char* tag)
    {
        solver.preconditioner().setDroptol(opt.iluDrop);
        solver.preconditioner().setFillfactor(opt.iluFill);
        solver.setMaxIterations(opt.maxIters);
        solver.setTolerance(opt.tol);

        bool mustRebuild = rebuild || !haveFactor
            || lastFill != opt.iluFill
            || lastDrop != opt.iluDrop;

        if (mustRebuild) {
            solver.compute(Aeq);
            haveFactor = (solver.info() == Eigen::Success);
            lastFill = opt.iluFill;
            lastDrop = opt.iluDrop;
            if (report) {
                std::cout << "[ReuseSolver] rebuild preconditioner (" << tag << ")\n";
            }
        }
        return haveFactor;
    }

    bool solveSystem(
        const SparseSystemCOO& sys,
        std::vector<double>& x,
        const LinearSolverOptions& opt,
        bool rebuildPrecond,
        bool report,
        const char* tag,
        double* errOut = nullptr)
    {
        if (!patternBuilt || sys.n != n) return false;
        if (!ensurePreconditioner(opt, rebuildPrecond, report, tag)) return false;

        Eigen::VectorXd beq = Eigen::VectorXd::Zero(n);
        for (int i = 0; i < sys.n && i < static_cast<int>(sys.b.size()); ++i) beq[i] = sys.b[i];
        if (useScaling) {
            for (int i = 0; i < n; ++i) beq[i] *= Dinvsqrt[i];
        }

        Eigen::VectorXd x0eq = Eigen::VectorXd::Zero(n);
        if (x.size() == static_cast<size_t>(n)) {
            for (int i = 0; i < n; ++i) {
                x0eq[i] = x[i];
                if (useScaling) x0eq[i] *= Dinvsqrt[i];
            }
        }

        Eigen::VectorXd sol_eq = solver.solveWithGuess(beq, x0eq);
        bool ok = (solver.info() == Eigen::Success);
        if (!ok) {
            solver.compute(Aeq);
            haveFactor = (solver.info() == Eigen::Success);
            if (haveFactor) {
                sol_eq = solver.solveWithGuess(beq, x0eq);
                ok = (solver.info() == Eigen::Success);
            }
        }
        if (!ok) return false;

        if (x.size() != static_cast<size_t>(n)) x.resize(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double val = sol_eq[i];
            if (useScaling) val *= Dinvsqrt[i];
            x[i] = val;
        }

        if (errOut) {
            Eigen::VectorXd resid = Aeq * sol_eq - beq;
            *errOut = resid.lpNorm<Eigen::Infinity>();
        }
        return true;
    }
};

struct OuterIterRuntime
{
        struct FieldHistory
        {
            std::vector<double> prev;      // x^{k}
            std::vector<double> prevPrev;  // x^{k-1}
            double omega = 1.0;            // 上一次的 Aitken 系数/占位符
            bool hasPrev = false;
            bool hasPrevPrev = false;
            int rejectCount = 0;
            int suspendIters = 0;
        };

    FieldHistory pHist;
    FieldHistory THist;

    FaceSignMask faceSignMask;         // 面质量通量符号掩码
    std::vector<int> flippedFaces;     // 最近一次符号翻转的面索引

    std::vector<double> scratchP;      // 通用工作缓冲
    std::vector<double> scratchT;

    ReuseLinearSolver reuseP;
    ReuseLinearSolver reuseT;

    int pressurePatternStamp = -1;
    int temperaturePatternStamp = -1;
    int precondCountdownP = 0;
    int precondCountdownT = 0;

    double lastDpInner = std::numeric_limits<double>::infinity();
    double lastCFL_T = 0.0;
    double lastDt = -1.0;
};

inline void ensureOuterRuntimeCapacity(const Mesh& mesh, OuterIterRuntime& rt)
{
    const size_t nCells = mesh.getCells().size();
    const size_t nFaces = mesh.getFaces().size();

    auto ensureHistory = [nCells](OuterIterRuntime::FieldHistory& h)
    {
        if (h.prev.size() != nCells) {
            h.prev.assign(nCells, 0.0);
            h.hasPrev = false;
        }
        if (h.prevPrev.size() != nCells) {
            h.prevPrev.assign(nCells, 0.0);
            h.hasPrevPrev = false;
        }
    };

    ensureHistory(rt.pHist);
    ensureHistory(rt.THist);

    if (rt.scratchP.size() != nCells) rt.scratchP.assign(nCells, 0.0);
    if (rt.scratchT.size() != nCells) rt.scratchT.assign(nCells, 0.0);

    rt.faceSignMask.ensureSize(nFaces);
    if (!rt.flippedFaces.empty()) rt.flippedFaces.clear();
}

inline double clampScalar(double value, double minVal, double maxVal)
{
    if (value < minVal) return minVal;
    if (value > maxVal) return maxVal;
    return value;
}

inline bool shouldUseReuseSolver(const LinearSolverOptions& opt, bool enabled)
{
    return enabled && opt.type == LinearSolverOptions::Type::BiCGSTAB;
}

inline void updateHistoryNoAitken(
    const std::vector<double>& values,
    OuterIterRuntime::FieldHistory& hist)
{
    const size_t n = values.size();
    if (hist.prev.size() != n) {
        hist.prev.assign(n, 0.0);
    }
    if (hist.prevPrev.size() != n) {
        hist.prevPrev.assign(n, 0.0);
    }
    if (hist.hasPrev) {
        hist.prevPrev = hist.prev;
        hist.hasPrevPrev = true;
    }
    hist.prev = values;
    hist.hasPrev = true;
    hist.omega = 1.0;
}

inline bool aitken_update_inplace(
    std::vector<double>& latest,
    OuterIterRuntime::FieldHistory& hist,
    double omegaMin,
    double omegaMax)
{
    const size_t n = latest.size();
    if (n == 0) return false;

    if (hist.prev.size() != n) {
        hist.prev.assign(n, 0.0);
        hist.hasPrev = false;
    }
    if (hist.prevPrev.size() != n) {
        hist.prevPrev.assign(n, 0.0);
        hist.hasPrevPrev = false;
    }

    if (!hist.hasPrev) {
        hist.prev = latest;
        hist.hasPrev = true;
        hist.hasPrevPrev = false;
        hist.omega = 1.0;
        return false;
    }
    if (!hist.hasPrevPrev) {
        hist.prevPrev = hist.prev;
        hist.hasPrevPrev = true;
        hist.prev = latest;
        return false;
    }

    double numer = 0.0;
    double denom = 0.0;

    for (size_t i = 0; i < n; ++i) {
        const double prev = hist.prev[i];
        const double prevPrev = hist.prevPrev[i];
        const double dk = latest[i] - prev;
        const double dkm1 = prev - prevPrev;
        const double diff = dk - dkm1;
        numer += dk * diff;
        denom += diff * diff;
    }

    if (denom <= 0.0 || !std::isfinite(numer) || !std::isfinite(denom)) {
        hist.prevPrev = hist.prev;
        hist.prev = latest;
        return false;
    }

    double omega = -numer / denom;
    if (!std::isfinite(omega)) omega = 1.0;
    omega = clampScalar(omega, omegaMin, omegaMax);

    for (size_t i = 0; i < n; ++i) {
        const double prev = hist.prev[i];
        const double dk = latest[i] - prev;
        latest[i] = prev + omega * dk;
    }

    hist.prevPrev = hist.prev;
    hist.prev = latest;
    hist.hasPrev = true;
    hist.hasPrevPrev = true;
    hist.omega = omega;
    return true;
}

inline void applyAitkenToField(
    FieldRegistry& reg,
    const std::string& fieldName,
    bool enable,
    double omegaMin,
    double omegaMax,
    OuterIterRuntime::FieldHistory& hist,
    std::vector<double>& scratch)
{
    auto field = reg.get<volScalarField>(fieldName);
    if (!field) return;
    scratch.assign(field->data.begin(), field->data.end());

    auto commitState = [&](const std::vector<double>& vec, double omegaVal) {
        field->data = vec;
        if (hist.hasPrev) {
            hist.prevPrev = hist.prev;
            hist.hasPrevPrev = hist.hasPrev;
        }
        hist.prev = vec;
        hist.hasPrev = true;
        hist.omega = omegaVal;
    };

    if (!enable) {
        commitState(scratch, 1.0);
        hist.rejectCount = 0;
        hist.suspendIters = 0;
        return;
    }

    if (hist.suspendIters > 0) {
        --hist.suspendIters;
        commitState(scratch, 1.0);
        return;
    }

    if (!hist.hasPrev) {
        commitState(scratch, 1.0);
        hist.hasPrevPrev = false;
        hist.rejectCount = 0;
        hist.suspendIters = 0;
        return;
    }
    if (!hist.hasPrevPrev) {
        hist.prevPrev = hist.prev;
        hist.hasPrevPrev = true;
        commitState(scratch, 1.0);
        hist.rejectCount = 0;
        hist.suspendIters = 0;
        return;
    }

    const size_t n = scratch.size();
    if (hist.prev.size() != n || hist.prevPrev.size() != n) {
        commitState(scratch, 1.0);
        return;
    }

    std::vector<double> dk(n, 0.0);
    double baseMetric = 0.0;
    for (size_t i = 0; i < n; ++i) {
        dk[i] = scratch[i] - hist.prev[i];
        baseMetric = std::max(baseMetric, std::abs(dk[i]));
    }
    if (baseMetric <= 0.0) {
        commitState(scratch, 1.0);
        hist.rejectCount = 0;
        hist.suspendIters = 0;
        return;
    }

    double numer = 0.0;
    double denom = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const double prev = hist.prev[i];
        const double prevPrev = hist.prevPrev[i];
        const double dkCurr = dk[i];
        const double dkPrev = prev - prevPrev;
        const double diff = dkCurr - dkPrev;
        numer += dkCurr * diff;
        denom += diff * diff;
    }

    if (denom <= 0.0 || !std::isfinite(numer) || !std::isfinite(denom)) {
        commitState(scratch, 1.0);
        hist.rejectCount = 0;
        hist.suspendIters = 0;
        return;
    }

    double omega = clampScalar(-numer / denom, omegaMin, omegaMax);
    const double minOmega = omegaMin;
    const int maxAttempts = 3;
    const double acceptanceTol = 1e-12;
    std::vector<double> candidate(n, 0.0);
    bool accepted = false;
    double acceptedOmega = omega;

    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        for (size_t i = 0; i < n; ++i) {
            candidate[i] = hist.prev[i] + omega * dk[i];
        }
        double metric = 0.0;
        for (size_t i = 0; i < n; ++i) {
            metric = std::max(metric, std::abs(candidate[i] - hist.prev[i]));
        }
        if (metric <= baseMetric + acceptanceTol) {
            accepted = true;
            acceptedOmega = omega;
            break;
        }
        omega *= 0.5;
        if (omega < minOmega * 0.25) break;
    }

    if (accepted) {
        scratch = candidate;
        commitState(scratch, acceptedOmega);
        hist.rejectCount = 0;
        hist.suspendIters = 0;
    }
    else {
        commitState(scratch, 1.0);
        hist.rejectCount++;
        if (hist.rejectCount >= 3) {
            hist.suspendIters = 2;
            hist.rejectCount = 0;
        }
    }
}

//------------------------------------------------------------------------------
// Pressure assemble+solve: 统一处理有井/无井并返回 dp_inf
//------------------------------------------------------------------------------
        inline bool assembleSolvePressureInner(
            MeshManager& mgr,
            FieldRegistry& reg,
            FaceFieldRegistry& freg,
            const PressureBCAdapter& Pbc,
            const Vector& g,
            double dt,
            const SolverControls& ctrl,
            const OperatorFieldNames& nmP,
            const std::vector<WellConfig>& wellsCfg_in,
            std::vector<WellConfig>& wellsCfg_work,
            std::vector<WellDOF>& wellsWork,
            double& dp_inf_out,
            SparseSystemCOO* lastSys = nullptr,
            OuterIterRuntime* runtime = nullptr)
{
    Mesh& mesh = mgr.mesh();
    const bool hasWells = !wellsCfg_in.empty();

    if (!FVM::Diffusion::build_FaceCoeffs_Central(
        mgr, reg, freg,
        nmP.a_f_diff, nmP.s_f_diff,
        "p_g",
        { "kxx:kxx","kyy:kyy","kzz:kzz","/mu_g","rho:rho_g" },
        "rho_g",
        FVM::Diffusion::RhoFaceMethod::Linear,
        g, Pbc,
        /*massForm=*/false, /*alpha_anisotropy=*/0))
        return false;

    if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Flow(
        mgr, reg, dt,
        "c_phi", "phi",
        "p_g_old", "rho_g",
        "p_g_prev", "rho_g",
        "Drho_Dp_g",
        nmP.a_time, nmP.b_time))
        return false;

    SparseSystemCOO sysp;
    if (!assemble_COO(mgr, reg, freg, "ddt+diffusion", nmP, &sysp)) return false;
    if (lastSys) *lastSys = sysp;

    int Np = 0; auto lid_p = buildUnknownMap(mesh, Np);
    auto pvec = gatherFieldToVec(reg, mesh, "p_g", lid_p, Np);

    if (hasWells) {
        wellsCfg_work = wellsCfg_in;
        build_masks_and_PI_for_all(mgr, reg, wellsCfg_work);

        wellsWork.clear();
        const int Nc = (int)mesh.getCells().size();
        const int Ntot = register_well_dofs_for_all(Nc, wellsCfg_work, wellsWork);
        extend_linear_system_size(sysp, Ntot);

        std::vector<int> lid_cell(Nc);
        for (int i = 0; i < Nc; ++i) lid_cell[i] = i;

        for (const auto& w : wellsWork) {
            add_peaceman_coupling_cell_rows(
                sysp, mesh, reg, w.PI_field, w.mask_field,
                lid_cell, w.lid, /*scaleByCellMeasure=*/false);
            add_well_row(
                sysp, mesh, reg, w.PI_field, w.mask_field,
                lid_cell, w.lid, w.mode, w.target,
                /*scaleByCellMeasure=*/false);
        }

        pvec.resize(Ntot, 0.0);
        for (const auto& w : wellsWork) {
            pvec[w.lid] = (w.mode == WellDOF::Mode::Pressure) ? w.target : pvec[0];
        }
    }

    auto optP = ctrl.lin_p;
    if (optP.tol <= 0.0) optP.tol = ctrl.tol_p_abs;

    double resP = 0.0; int itP = 0;
    bool solvedWithReuse = false;
    const bool canReuse = runtime && shouldUseReuseSolver(optP, ctrl.reuse_linear_pattern);
    if (canReuse) {
        bool rebuilt = false;
        if (runtime->reuseP.ensurePattern(sysp, optP, ctrl.reportPerIter, "pressure", rebuilt)) {
            if (rebuilt && ctrl.rebuild_precond_every > 0) runtime->precondCountdownP = 0;
            bool rebuildPrecond = (ctrl.rebuild_precond_every <= 0) ? true : (runtime->precondCountdownP <= 0);
            solvedWithReuse = runtime->reuseP.solveSystem(
                sysp, pvec, optP, rebuildPrecond, ctrl.reportPerIter, "pressure", &resP);
            if (solvedWithReuse && ctrl.rebuild_precond_every > 0) {
                runtime->precondCountdownP = rebuildPrecond ? ctrl.rebuild_precond_every : std::max(0, runtime->precondCountdownP - 1);
            }
        }
        else {
            runtime->reuseP.reset();
        }
    }

    if (!solvedWithReuse) {
        if (!solveCOO_Eigen(sysp, pvec, optP, &itP, &resP)) return false;
        if (runtime) {
            runtime->reuseP.reset();
            runtime->precondCountdownP = ctrl.rebuild_precond_every;
        }
    }

            if (hasWells) {
                std::vector<double> p_cells(Np);
                std::copy(pvec.begin(), pvec.begin() + Np, p_cells.begin());
                scatterVecToField(reg, mesh, "p_g", lid_p, p_cells);
                writeback_pw_fields_for_all(reg, wellsWork, pvec);
            }
            else {
                scatterVecToField(reg, mesh, "p_g", lid_p, pvec);
            }

            if (runtime) {
                applyAitkenToField(
                    reg,
                    "p_g",
                    ctrl.enable_Aitken_p,
                    ctrl.aitken_omega_min,
                    ctrl.aitken_omega_max,
                    runtime->pHist,
                    runtime->scratchP);
            }

            underRelaxInPlace(reg, "p_g", "p_g_prev", ctrl.urf_p);
            dp_inf_out = maxAbsDiff(reg, "p_g", "p_g_prev");
            updatePrevIterates(reg, { {"p_g","p_g_prev"} });
            return true;
        }

struct PressureStageResult
{
    double dpInf = std::numeric_limits<double>::infinity();
    int    sweepsPerformed = 0;
    bool   metInnerTol = false;
};

inline bool runPressureSweeps(
            MeshManager& mgr,
            FieldRegistry& reg,
            FaceFieldRegistry& freg,
            const PressureBCAdapter& Pbc,
            const Vector& g,
            double dt,
            const SolverControls& ctrl,
            const OperatorFieldNames& nmP,
            const std::vector<WellConfig>& wellsCfg_in,
            std::vector<WellConfig>& wellsCfg_work,
            std::vector<WellDOF>& wellsWork,
            PressureStageResult& outRes,
            SparseSystemCOO* lastSys = nullptr,
            OuterIterRuntime* runtime = nullptr)
{
    outRes = {};
    const int maxSweeps = std::max(1, ctrl.NsweepP_max);
    for (int sweep = 0; sweep < maxSweeps; ++sweep)
    {
        double dpInner = 0.0;
                if (!assembleSolvePressureInner(
                        mgr, reg, freg, Pbc, g, dt,
                        ctrl, nmP, wellsCfg_in,
                        wellsCfg_work, wellsWork,
                        dpInner, lastSys, runtime))
                    return false;

        outRes.dpInf = dpInner;
        outRes.sweepsPerformed = sweep + 1;
        if (dpInner < ctrl.tol_p_inner) 
        {
            outRes.metInnerTol = true;
            break;
        }
    }
    return true;
}

inline bool rebuildDarcyFluxes(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const OperatorFieldNames& nmP,
    const PressureBCAdapter& Pbc)
{
    if (!FVM::Convection::buildFlux_Darcy_Mass(
            mgr, reg, freg,
            nmP.a_f_diff, nmP.s_f_diff,
            "p_g", "rho_g",
            "mf_g", "Qf_g", "ufn_g",
            &Pbc, true))
        return false;

    debugCheckMassFlux(mgr, freg, "mf_g", 1e-20);
    return true;
}

inline FaceSignUpdateInfo refreshFaceSignMaskRuntime(
    FaceFieldRegistry& freg,
    const std::string& fluxField,
    double epsilon,
    OuterIterRuntime& runtime)
{
    auto flux = freg.get<faceScalarField>(fluxField.c_str());
    if (!flux) {
        std::cerr << "[Accel] Missing flux field '" << fluxField << "' for face-sign mask.\n";
        runtime.faceSignMask.reset();
        runtime.flippedFaces.clear();
        return {};
    }
    return updateFaceSignMask_fromFlux(*flux, epsilon, runtime.faceSignMask, runtime.flippedFaces);
}

inline double computeMaxCFL_T(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    double dt,
    double* out_sumInf = nullptr)
{
    Mesh& mesh = mgr.mesh();
    const auto& faces = mesh.getFaces();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();

    auto Qf = freg.get<faceScalarField>("Qf_g");
    auto mf = freg.get<faceScalarField>("mf_g");
    auto rho = reg.get<volScalarField>("rho_g");

    std::vector<double> sumFlux(cells.size(), 0.0);

    const bool hasQf = static_cast<bool>(Qf);
    const bool hasMf = static_cast<bool>(mf);

    for (const auto& F : faces) {
        const int iF = F.id - 1;
        double q = 0.0;
        if (hasQf) {
            q = std::abs((*Qf)[iF]);
        }
        else if (hasMf && rho) {
            const int P = F.ownerCell;
            const int N = F.neighborCell;
            const size_t iP = (P >= 0 ? id2idx.at(P) : 0);
            const size_t iN = (N >= 0 ? id2idx.at(N) : iP);
            const double rhoF = std::max(1e-12, 0.5 * ((*rho)[iP] + (*rho)[iN]));
            q = std::abs((*mf)[iF]) / rhoF;
        }
        else {
            continue; // no flux information
        }

        if (F.ownerCell >= 0) {
            const size_t iP = id2idx.at(F.ownerCell);
            sumFlux[iP] += q;
        }
        if (!F.isBoundary() && F.neighborCell >= 0) {
            const size_t iN = id2idx.at(F.neighborCell);
            sumFlux[iN] += q;
        }
    }

    double maxCFL = 0.0;
    double maxSum = 0.0;
    for (size_t i = 0; i < cells.size(); ++i) {
        const double V = std::max(0.0, cells[i].volume);
        if (V <= 0.0) continue;
        const double cfl = dt * sumFlux[i] / V;
        if (cfl > maxCFL) maxCFL = cfl;
        if (sumFlux[i] > maxSum) maxSum = sumFlux[i];
    }
    if (out_sumInf) *out_sumInf = maxSum;
    return maxCFL;
}

inline bool solveTemperatureStage(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const TemperatureBCAdapter& Tbc,
    const Vector& g,
    double dt,
    const SolverControls& ctrl,
    const OperatorFieldNames& nmT,
    double& dT_inf_out,
    SparseSystemCOO* lastSys,
    OuterIterRuntime* runtime,
    const std::vector<WellDOF>* wellsForHeat = nullptr)
{
    if (!FVM::Timeterm::TimeTerm_FullyImplicit_SinglePhase_Temperature(
            mgr, reg, dt,
            "C_eff", "T_old",
            nmT.a_time, nmT.b_time))
        return false;

    if (!FVM::Diffusion::build_FaceCoeffs_Central(
            mgr, reg, freg,
            nmT.a_f_diff, nmT.s_f_diff,
            "T",
            { "iso:lambda_eff" },
            "",
            FVM::Diffusion::RhoFaceMethod::Linear,
            g, Tbc,
            /*massForm=*/false, /*alpha_anisotropy=*/0))
        return false;

    debugCheckMassFlux(mgr, freg, "mf_g", 1e-20);

    if (!FVM::Convection::build_FaceCoeffs_Upwind(
            mgr, reg, freg,
            "T", "mf_g",
            { "cp_g" },
            nmT, Tbc))
        return false;

    bool hasWellHeat = false;
    if (wellsForHeat && !wellsForHeat->empty()) {
        bool first = true;
        for (const auto& w : *wellsForHeat) {
            const double Tin = (w.role == WellDOF::Role::Injector) ? w.Tin : 0.0;
            if (!FVM::SourceTerm::add_temperature_source_from_qm_single(
                    mgr, reg,
                    w.role,
                    w.PI_field.c_str(), w.mask_field.c_str(), w.name.c_str(),
                    Tin,
                    "p_g", "cp_g",
                    1.0,
                    !first,
                    false))
                return false;
            first = false;
            hasWellHeat = true;
        }
    }

    SparseSystemCOO sysT;
    const char* opName = hasWellHeat ? "ddt+diffusion+convection+source" : "ddt+diffusion+convection";
    if (!assemble_COO(mgr, reg, freg, opName, nmT, &sysT))
        return false;
    if (lastSys) *lastSys = sysT;

    Mesh& mesh = mgr.mesh();
    int Nt = 0; auto lid_t = buildUnknownMap(mesh, Nt);
    auto Tvec = gatherFieldToVec(reg, mesh, "T", lid_t, Nt);

    auto optT = ctrl.lin_T;
    if (optT.tol <= 0.0) optT.tol = ctrl.tol_T_abs;

    double resT = 0.0; int itT = 0;
    bool solvedWithReuse = false;
    const bool canReuse = runtime && shouldUseReuseSolver(optT, ctrl.reuse_linear_pattern);
    if (canReuse) {
        bool rebuilt = false;
        if (runtime->reuseT.ensurePattern(sysT, optT, ctrl.reportPerIter, "temperature", rebuilt)) {
            if (rebuilt && ctrl.rebuild_precond_every > 0) runtime->precondCountdownT = 0;
            bool rebuildPrecond = (ctrl.rebuild_precond_every <= 0) ? true : (runtime->precondCountdownT <= 0);
            solvedWithReuse = runtime->reuseT.solveSystem(
                sysT, Tvec, optT, rebuildPrecond, ctrl.reportPerIter, "temperature", &resT);
            if (solvedWithReuse && ctrl.rebuild_precond_every > 0) {
                runtime->precondCountdownT = rebuildPrecond ? ctrl.rebuild_precond_every : std::max(0, runtime->precondCountdownT - 1);
            }
        }
        else {
            runtime->reuseT.reset();
        }
    }

    if (!solvedWithReuse) {
        if (!solveCOO_Eigen(sysT, Tvec, optT, &itT, &resT)) return false;
        if (runtime) {
            runtime->reuseT.reset();
            runtime->precondCountdownT = ctrl.rebuild_precond_every;
        }
    }

    scatterVecToField(reg, mesh, "T", lid_t, Tvec);

    if (runtime) {
        applyAitkenToField(
            reg,
            "T",
            ctrl.enable_Aitken_T,
            ctrl.aitken_omega_min,
            ctrl.aitken_omega_max,
            runtime->THist,
            runtime->scratchT);
    }
    underRelaxInPlace(reg, "T", "T_prev", ctrl.urf_T);
    dT_inf_out = maxAbsDiff(reg, "T", "T_prev");
    updatePrevIterates(reg, { {"T","T_prev"} });
    return true;
}

inline bool solveTemperatureWithCFLGuard(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const TemperatureBCAdapter& Tbc,
    const Vector& g,
    double dt,
    const SolverControls& ctrl,
    const OperatorFieldNames& nmT,
    double& dT_inf_out,
    SparseSystemCOO* lastSys = nullptr,
    OuterIterRuntime* runtime = nullptr,
    const std::vector<WellDOF>* wellsForHeat = nullptr)
{
    if (!ctrl.enable_CFL_guard || ctrl.CFL_T_threshold <= 0.0) {
        if (runtime) runtime->lastCFL_T = 0.0;
        return solveTemperatureStage(mgr, reg, freg, Tbc, g, dt, ctrl, nmT, dT_inf_out, lastSys, runtime, wellsForHeat);
    }

    double sumInf = 0.0;
    const double maxCFL = computeMaxCFL_T(mgr, reg, freg, dt, &sumInf);
    if (runtime) runtime->lastCFL_T = maxCFL;

    int Nsub = 1;
    if (maxCFL > ctrl.CFL_T_threshold) {
        Nsub = static_cast<int>(std::ceil(maxCFL / std::max(1e-12, ctrl.CFL_T_threshold)));
    }
    const int NsubCap = 1000;
    if (Nsub > NsubCap) Nsub = NsubCap;

    double dt_sub = dt / std::max(1, Nsub);
    const double dt_min = dt * std::max(0.0, ctrl.CFL_dt_scale_min);
    if (dt_sub < dt_min) {
        dt_sub = dt_min;
        Nsub = static_cast<int>(std::ceil(dt / dt_sub));
    }

    if (ctrl.reportPerIter) {
        std::cout << "[CFLGuard] maxCFL_T=" << maxCFL
                  << "  Nsub=" << Nsub
                  << "  dt_sub=" << dt_sub << "\n";
    }

    bool ok = true;
    double dT_last = 0.0;
    for (int k = 0; k < Nsub; ++k) {
        SparseSystemCOO* outSys = (k == Nsub - 1) ? lastSys : nullptr;
        ok = solveTemperatureStage(mgr, reg, freg, Tbc, g, dt_sub, ctrl, nmT, dT_last, outSys, runtime, wellsForHeat);
        if (!ok) break;
    }
    dT_inf_out = dT_last;
    return ok;
}

inline bool rebuild_convection_face_one(
    Mesh& mesh,
    const std::vector<Face>& faces,
    const std::unordered_map<int, int>& id2idx,
    const std::vector<std::shared_ptr<volScalarField>>& chi_ptrs,
    faceScalarField& flux,
    faceScalarField& aPP,
    faceScalarField& aPN,
    faceScalarField& bP,
    const TemperatureBCAdapter& Tbc,
    int idx)
{
    if (idx < 0 || static_cast<size_t>(idx) >= faces.size()) return true;

    const Face& F = faces[idx];
    const double mf = flux[idx];

    double& aPPv = aPP[idx];
    double& aPNv = aPN[idx];
    double& bPv = bP[idx];
    aPPv = 0.0;
    aPNv = 0.0;
    bPv = 0.0;

    const double epsF = 1e-18;
    if (std::abs(mf) <= epsF) return true;

    const int ownerId = F.ownerCell;
    if (ownerId < 0) return true;
    const size_t iP = static_cast<size_t>(id2idx.at(ownerId));

    auto product_at_cell = [&](size_t iCell)->double {
        double w = 1.0;
        for (const auto& chi : chi_ptrs) w *= (*chi)[iCell];
        return w;
    };

    if (!F.isBoundary()) {
        const int neighborId = F.neighborCell;
        if (neighborId < 0) return true;
        const size_t iN = static_cast<size_t>(id2idx.at(neighborId));
        const size_t iUP = (mf >= 0.0) ? iP : iN;
        const double Wup = product_at_cell(iUP);
        const double beta = mf * Wup;
        if (beta >= 0.0) {
            aPPv = beta;
            aPNv = 0.0;
        }
        else {
            aPPv = 0.0;
            aPNv = beta;
        }
    }
    else {
        if (mf > 0.0) {
            const double WP = product_at_cell(iP);
            aPPv += mf * WP;
        }
        else {
            double Wb = 1.0;
            for (size_t k = 0; k < chi_ptrs.size(); ++k) {
                Wb *= (*chi_ptrs[k])[iP];
            }
            double phi_b = 0.0;
            if (getDirichletFromABC(Tbc, F.id, phi_b)) {
                bPv += (-mf) * Wb * phi_b;
            }
        }
    }
    return true;
}

inline bool rebuildConvectionFaceSubset(
    MeshManager& mgr,
    const FieldRegistry& reg,
    FaceFieldRegistry& freg,
    const TemperatureBCAdapter& Tbc,
    const OperatorFieldNames& nm,
    const std::string& fluxFieldName,
    const std::vector<std::string>& upwindFields,
    const std::vector<int>& faceIndices)
{
    if (faceIndices.empty()) return true;

    auto fluxF = freg.get<faceScalarField>(fluxFieldName.c_str());
    if (!fluxF) {
        std::cerr << "[Accel] Missing face flux '" << fluxFieldName << "' for partial convection rebuild.\n";
        return false;
    }

    auto aPP = freg.get<faceScalarField>(nm.aPP_conv.c_str());
    auto aPN = freg.get<faceScalarField>(nm.aPN_conv.c_str());
    auto bP = freg.get<faceScalarField>(nm.bP_conv.c_str());
    if (!aPP || !aPN || !bP) {
        std::cerr << "[Accel] Missing convection coefficient fields for partial rebuild.\n";
        return false;
    }

    std::vector<std::shared_ptr<volScalarField>> chi_ptrs;
    chi_ptrs.reserve(upwindFields.size());
    for (const auto& field : upwindFields) {
        auto chi = reg.get<volScalarField>(field.c_str());
        if (!chi) {
            std::cerr << "[Accel] Missing upwind scalar field '" << field << "' for partial rebuild.\n";
            return false;
        }
        chi_ptrs.emplace_back(chi);
    }

    Mesh& mesh = mgr.mesh();
    auto& faces = const_cast<std::vector<Face>&>(mesh.getFaces());
    const auto& id2idx = mesh.getCellId2Index();
    for (int idx : faceIndices) {
        if (!rebuild_convection_face_one(
                mesh, faces, id2idx, chi_ptrs,
                *fluxF, *aPP, *aPN, *bP, Tbc, idx))
            return false;
    }
    return true;
}

inline bool outerIter_constProperties_singlePhase_CO2_T_H_withWell_accel(
    MeshManager& mgr,
    FieldRegistry& reg,
    FaceFieldRegistry& freg,
    PhysicalPropertiesManager& ppm,
    const TemperatureBCAdapter& Tbc,
    const PressureBCAdapter& Pbc,
    const Vector& g,
    const std::vector<WellConfig>& wellsCfg_in,
    double dt,
    const SolverControls& ctrl,
    OuterIterRuntime& runtime)
{
    Mesh& mesh = mgr.mesh();
    ensureOuterRuntimeCapacity(mesh, runtime);

    double prev_dp_g = std::numeric_limits<double>::infinity();
    double prev_dT = std::numeric_limits<double>::infinity();

    std::vector<WellConfig> wellsCfg_work;
    std::vector<WellDOF> wellsWork;

    for (int it = 0; it < ctrl.maxOuter; ++it) {
        if (!startOuterIteration_scatter(reg, "p_g", "p_g_prev")) return false;
        if (!startOuterIteration_T(reg, "T", "T_prev")) return false;

        ppm.CO2Properties_test_constProperties_singlePhase_CO2(mgr, reg);
        ppm.RockProperties_test_constProperties_singlePhase_CO2(mgr, reg);
        ppm.ComputeEffectiveThermalProperties_constProperties_singlePhase_CO2_T_H(mgr, reg);

        const OperatorFieldNames nmP = makeNames("p_g");
        const OperatorFieldNames nmT = makeNames("T");

        SparseSystemCOO lastP;
        SparseSystemCOO lastT;

        const double prevDpInner = runtime.lastDpInner;
        PressureStageResult presRes;
        if (!runPressureSweeps(
                mgr, reg, freg, Pbc, g, dt, ctrl,
                nmP, wellsCfg_in, wellsCfg_work, wellsWork,
                presRes,
                ctrl.dumpMMOnLastIter ? &lastP : nullptr,
                &runtime))
            return false;
        runtime.lastDpInner = presRes.dpInf;

        if (!rebuildDarcyFluxes(mgr, reg, freg, nmP, Pbc)) return false;

        const auto maskInfo = refreshFaceSignMaskRuntime(
            freg, "mf_g", ctrl.flux_sign_epsilon, runtime);
        const bool hasPrevMask = maskInfo.hadPreviousMask;
        double dpRatio = std::numeric_limits<double>::infinity();
        if (std::isfinite(prevDpInner) && prevDpInner > 0.0) {
            dpRatio = presRes.dpInf / std::max(prevDpInner, 1e-30);
        }

        bool convectionOk = true;
        const std::vector<std::string> upwindChi = { "cp_g" };
        const bool incrementalCandidate =
            ctrl.enable_incremental_convection &&
            hasPrevMask &&
            (dpRatio < ctrl.kappa_p_for_flux_update);

        if (incrementalCandidate) {
            const double flipRatio = maskInfo.flipRatio();
            if (maskInfo.flippedCount == 0 && flipRatio <= ctrl.incremental_flip_ratio) {
                // reuse previous convection operator
                if (ctrl.reportPerIter) {
                    std::cout << "[Accel] reuse convection coefficients (dpRatio="
                              << dpRatio << ", flipRatio=" << flipRatio << ")\n";
                }
            }
            else if (flipRatio <= ctrl.incremental_flip_ratio) {
                if (ctrl.reportPerIter) {
                    std::cout << "[Accel] partial convection rebuild for "
                              << maskInfo.flippedCount << " faces (flipRatio="
                              << flipRatio << ")\n";
                }
                convectionOk = rebuildConvectionFaceSubset(
                    mgr, reg, freg, Tbc, nmT, "mf_g", upwindChi, runtime.flippedFaces);
            }
            else {
                if (ctrl.reportPerIter) {
                    std::cout << "[Accel] fallback to full convection rebuild (flipRatio="
                              << flipRatio << ")\n";
                }
                convectionOk = FVM::Convection::build_FaceCoeffs_Upwind(
                    mgr, reg, freg,
                    "T", "mf_g",
                    upwindChi,
                    nmT, Tbc);
            }
        }
        else {
            convectionOk = FVM::Convection::build_FaceCoeffs_Upwind(
                mgr, reg, freg,
                "T", "mf_g",
                upwindChi,
                nmT, Tbc);
        }

        if (!convectionOk) return false;

        const std::vector<WellDOF>* wellsForHeat = wellsWork.empty() ? nullptr : &wellsWork;
        double dT = 0.0;
        if (!solveTemperatureWithCFLGuard(
                mgr, reg, freg, Tbc, g, dt,
                ctrl, nmT, dT,
                ctrl.dumpMMOnLastIter ? &lastT : nullptr,
                &runtime,
                wellsForHeat))
            return false;

        std::cout << "[Outer " << it << "]  sweeps=" << presRes.sweepsPerformed
                  << "  |Δp|_inf=" << presRes.dpInf
                  << "  |ΔT|_inf=" << dT << "\n";

        if (it > 0) {
            double rT = dT / std::max(prev_dT, 1e-30);
            double uT = (rT < 0.7) ? +0.05 : (rT > 0.95 ? -0.05 : 0.0);
            auto& ctrl_mut = const_cast<SolverControls&>(ctrl);
            ctrl_mut.urf_T = std::min(0.7, std::max(0.15, ctrl_mut.urf_T + uT));

            double rp = presRes.dpInf / std::max(prev_dp_g, 1e-30);
            double up = (rp < 0.7) ? +0.05 : (rp > 0.95 ? -0.05 : 0.0);
            ctrl_mut.urf_p = std::min(0.7, std::max(0.15, ctrl_mut.urf_p + up));
        }

        prev_dp_g = presRes.dpInf;
        prev_dT = dT;

        auto maxAbsField = [&](const std::string& fld) -> double {
            double m = 0.0;
            auto f = reg.get<volScalarField>(fld);
            if (!f) return 1.0;
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();
            for (const auto& c : cells) {
                if (c.id < 0) continue;
                size_t i = id2idx.at(c.id);
                m = std::max(m, std::abs((*f)[i]));
            }
            return std::max(1.0, m);
        };

        const double PScale = maxAbsField("p_g");
        const double TScale = maxAbsField("T");

        const bool convP = (presRes.dpInf < std::max(ctrl.tol_p_abs, ctrl.tol_p_rel * PScale));
        const bool convT = (dT < std::max(ctrl.tol_T_abs, ctrl.tol_T_rel * TScale));

        if (convP && convT) {
            std::cout << "Converged at outer iter " << it << "\n";

            if (ctrl.dumpMMOnLastIter) {
    #if __cplusplus >= 201703L
                try { std::filesystem::create_directories("mm"); }
                catch (...) {}
    #endif
                PostChecks::dumpCOO_to_matrix_market(
                    lastP,
                    "mm/A_P_CO2_pTH_accel.mtx",
                    "mm/b_P_CO2_pTH_accel.txt",
                    /*sym=*/false);
                PostChecks::dumpCOO_to_matrix_market(
                    lastT,
                    "mm/A_T_CO2_pTH_accel.mtx",
                    "mm/b_T_CO2_pTH_accel.txt",
                    /*sym=*/false);
            }
            return true;
        }

        if (it == ctrl.maxOuter - 1) {
            std::cout << "Reached maxOuter without meeting P/T tolerances.\n";
        }
    }

    return true;
}

} // namespace Accel
} // namespace Solver

