# PLAN_7D_CHECKLIST

## 0. Rules (Must Follow)

- This file is the daily execution tracker: target -> strict acceptance -> PASS sign-off.
- A day can be checked only when all implementation items and all acceptance gates are PASS.
- End of day must include evidence: commit id, log path, result summary.
- If a change touches `--case`, regression commands, or pass criteria, update in the same change:
  - `REGRESSION_CASES.md`
  - `PROJECT_CONTEXT.md` (at least `Last Updated` and `Git Commit`)

---

## D0 Baseline (Current State)

- [x] `main.cpp` has dispatcher mode (`--case / --list / --help`)
- [x] `PROJECT_CONTEXT.md` is available for new-session handover
- [x] `REGRESSION_CASES.md` defines minimum suite `R1-R5`

Evidence:
- Commit: `1b74f9c (dirty working tree baseline)`
- Date: `2026-03-04`

---

## Day 1 - Architecture Freeze and Connection Consistency

Status: `[ ]`

Implementation checklist:
- [ ] Freeze DOF mapping conventions for single-phase and two-phase paths
- [ ] Unify MM/FI/NNC/FF connection key and deterministic ordering
- [ ] Freeze connection stats output format (CI-friendly)

Strict acceptance gates:
- [ ] R4 PASS (2D transmissibility + FIM conservation)
- [ ] R5 PASS (3D transmissibility + FIM conservation)
- [ ] Connection stats are reproducible for same input

Evidence (fill):
- Commit:
- Logs:
- Verdict: PASS / FAIL

---

## Day 2 - Potential, Upwind, and Phase-Flux Kernel

Status: `[ ]`

Implementation checklist:
- [ ] Finalize potential-difference interface (including gravity)
- [ ] Finalize upwind-by-potential logic with correct AD derivative direction
- [ ] Integrate two-phase phase-flux assembly on shared kernel

Strict acceptance gates:
- [ ] R1 PASS (FVM-AD operator tests)
- [ ] `--case=grad_all` PASS
- [ ] No non-physical net flow in hydrostatic/gravity equilibrium case

Evidence (fill):
- Commit:
- Logs:
- Verdict: PASS / FAIL

---

## Day 3 - Boundary Conditions and Leakoff Discretization

Status: `[ ]`

Implementation checklist:
- [ ] Integrate Dirichlet/Neumann into residual and Jacobian
- [ ] Add parameterized leakoff source/sink with switch
- [ ] Ensure boundary derivatives enter correct blocks

Strict acceptance gates:
- [ ] `--case=boundary_export` PASS
- [ ] R1 PASS
- [ ] Linear pressure-drop and adiabatic patch tests PASS

Evidence (fill):
- Commit:
- Logs:
- Verdict: PASS / FAIL

---

## Day 4 - Well Model and Injection/Production Control (WAG Base)

Status: `[ ]`

Implementation checklist:
- [ ] Integrate well terms (BHP/Rate) consistently into grid equations
- [ ] Verify coupled response with leakoff on/off
- [ ] Build switchable schedule control skeleton for WAG extension

Strict acceptance gates:
- [ ] Single-well pressure response is physically consistent
- [ ] Leakoff on/off trend is correct and explainable
- [ ] At least one of R4/R5 PASS (depending on impact scope)

Evidence (fill):
- Commit:
- Logs:
- Verdict: PASS / FAIL

---

## Day 5 - Global Assembly and Jacobian Verification

Status: `[ ]`

Implementation checklist:
- [ ] Finish key block assembly for single-phase and two-phase equations
- [ ] Ensure accumulation and flux derivatives are assembled consistently
- [ ] Provide FD vs AD Jacobian comparison case/script

Strict acceptance gates:
- [ ] FD vs AD max relative error `< 1e-6`
- [ ] R1 PASS
- [ ] At least one of R4/R5 PASS for impacted dimension

Evidence (fill):
- Commit:
- Logs:
- Verdict: PASS / FAIL

---

## Day 6 - Transient Nonlinear Solver Robustness

Status: `[ ]`

Implementation checklist:
- [ ] Stable Newton/FIM or IMPES transient stepping path
- [ ] Limiter/rollback/adaptive-dt mechanisms are traceable and triggered when needed
- [ ] Enforce non-physical saturation protection

Strict acceptance gates:
- [ ] At least 50 transient steps are stable
- [ ] Convergence logs are complete (iters, residuals, dt history)
- [ ] R1 PASS

Evidence (fill):
- Commit:
- Logs:
- Verdict: PASS / FAIL

---

## Day 7 - Two-Phase Thermal Closed Loop and WAG Regression

Status: `[ ]`

Implementation checklist:
- [ ] Run a full non-isothermal two-phase closed-loop case with EDFM exchange
- [ ] Script one reproducible WAG scenario
- [ ] Export core paper metrics (pressure/temperature/saturation/net heat extraction)

Strict acceptance gates:
- [ ] R1-R5 all PASS
- [ ] At least one full WAG case output is generated
- [ ] Metric evolution is physically explainable (no obvious numerical artifacts)

Evidence (fill):
- Commit:
- Logs:
- Verdict: PASS / FAIL

---

## Daily Sign-off Log (Required)

| Date | Day | Owner | Commit | Required Checks | Result | Notes |
|---|---|---|---|---|---|---|
| YYYY-MM-DD | Day N |  |  | R? / case? | PASS/FAIL |  |

