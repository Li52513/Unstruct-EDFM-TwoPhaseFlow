# COMSOL Transient Solver Mapping

| COMSOL Setting | Code Equivalent | Notes |
| --- | --- | --- |
| Time stepping = Free/Adaptive | `TimeStepAdapt` (`ts`) configuration in `main.cpp` | Step-doubling error estimator (`ts.useErrorEstimator`) with `ts.safety`, `ts.facMin`, `ts.facMax` reproduces COMSOL's local time-error control. |
| Initial / Minimum / Maximum step size | `ts.dt`, `ts.dt_min`, `ts.dt_max` | Values can be staged to mimic COMSOL's "initial step" and "automatic upper bound". |
| Absolute tolerance per variable | `SolverControls::tol_p_abs`, `SolverControls::tol_T_abs` | Multiply both by any global "tolerance factor" when tuning. |
| Relative tolerance per variable | `SolverControls::tol_p_rel`, `SolverControls::tol_T_rel` | Used together with absolute tolerance when forming the normalised error. |
| Nonlinear solver damping | `SolverControls::urf_p`, `SolverControls::urf_T` (auto-adjusted each outer iteration) | Matches COMSOL's damping factor with adaptive updates driven by residual ratios. |
| Time integration method (BDF order) | `SolverControls::timeIntegrator` (`TimeIntegratorMethod::BDF1` or `BDF2`) | BDF2 uses variable-step coefficients with full state history (`*_old`, `*_old2`). |
| Error safety factor / step-size limits | `ts.safety`, `ts.facMin`, `ts.facMax`, `ts.errorAccept` | Equivalent to COMSOL's safety factor and min/max scaling applied to accepted steps. |
| Error estimator cadence | `ts.estimatorEvery` | Run step-doubling every *N* accepted steps (default every step). |
| Event / target time handling | `runTransient_singlePhase_adaptive` clamps the next step so that `t + dt` never overshoots `T_end` | Additional event times can be enforced by passing a smaller `dt` before calling the driver. |
| Output control | `writeEvery`, `onWrite`, and `export*` flags | Mirrors COMSOL output steps (per-step TXT/CSV, Tecplot, Matrix Market dumps). |

These mappings assume the driver is called through `runTransient_singlePhase_adaptive`, which now maintains the additional BDF2 history (`*_old2`) and evaluates the step-doubling error indicator before accepting each step.
