/**
 * @file Test_ValidationSuite.cpp
 * @brief Systematic validation test suite for the 2D FIM transient EDFM solver.
 *
 * 8 physics scenarios × 6 variants (Const/Var × NoFrac/SingleFrac/CrossFrac) = 48 functions.
 * Const_NoFrac variants compare against analytical solutions.
 * Const_SingleFrac/CrossFrac and Var_* variants verify numerical stability.
 *
 * All scenarios use RunGenericFIMTransient<N=2 or 3>.
 * Final state is read from FieldManager_2D after the run.
 */

#include "Test_ValidationSuite.h"

#include "FIM_TransientCaseKit.hpp"
#include "2D_PostProcess.h"
#include "BoundaryConditionManager.h"
#include "MeshDefinitions.h"
#include "MeshManager.h"
#include "SolverContrlStrName_op.h"
#include "Well_WellControlTypes.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define VS_MKDIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define VS_MKDIR(p) mkdir(p, 0777)
#endif

namespace ValidationSuite {

namespace {

// ============================================================
// Constants
// ============================================================
constexpr double kPi = 3.14159265358979323846;

// Fluid (water)
constexpr double kMu_w   = 1.0e-3;    // Pa·s
constexpr double kRho_w  = 1000.0;    // kg/m³
constexpr double kCp_w   = 4200.0;    // J/(kg·K)
constexpr double kCt_w   = 4.5e-10;   // Pa⁻¹ (fluid compressibility)
constexpr double kLam_w  = 0.6;       // W/(m·K)

// Rock
constexpr double kRho_r  = 2650.0;    // kg/m³
constexpr double kCp_r   = 920.0;     // J/(kg·K)
constexpr double kLam_r  = 3.0;       // W/(m·K)

// Two-phase (CO₂/water Corey params for analytical BL)
constexpr double kKrw_max = 0.5;
constexpr double kKro_max = 0.8;
constexpr double kNw      = 2.0;
constexpr double kNo      = 2.0;
constexpr double kSwi     = 0.2;
constexpr double kSor     = 0.2;
constexpr double kMu_o    = 5.0e-5;   // CO₂ viscosity [Pa·s]

// ============================================================
// Directory helper
// ============================================================
static void EnsureDir(const std::string& raw) {
    if (raw.empty()) return;
    std::string p = raw;
    for (char& c : p) if (c == '\\') c = '/';
    std::stringstream ss(p); std::string tok, cur;
    while (std::getline(ss, tok, '/')) {
        if (tok.empty() || tok == ".") continue;
        if (!cur.empty()) cur += "/";
        cur += tok;
        VS_MKDIR(cur.c_str());
    }
}

static bool FileNonEmpty(const std::string& p) {
    std::ifstream f(p, std::ios::binary | std::ios::ate);
    return f.good() && f.tellg() > 0;
}

// ============================================================
// Analytical solutions
// ============================================================

/**
 * @brief 1D diffusion: Fourier series solution P(x,t) for [0,L] with
 *   Dirichlet BC at x=0 (P_left), x=L (P_right), IC P(x,0)=P0.
 *   D = hydraulic or thermal diffusivity.
 */
static double Analytic1D(double x, double t, double L,
                          double P_left, double P_right, double P0,
                          double D, int nterms = 100)
{
    const double dP     = P_right - P_left;
    const double P_ss   = P_left + dP * (x / L);          // steady-state
    const double phi0   = P0 - P_ss;                       // IC deviation from steady-state

    // phi0 = P0 - P_left - (P_right-P_left)*x/L
    // bn = (2/L) * integral_0^L phi0 * sin(n*pi*x/L) dx
    //    = 2*(P0-P_left)/(n*pi)*(1-(-1)^n) - 2*(P_right-P_left)*(-1)^n/(n*pi)
    // simplified: phi0 = a - b*x/L  with a=P0-P_left, b=P_right-P_left
    const double a = P0 - P_left;
    const double b = dP;
    double sum = 0.0;
    for (int n = 1; n <= nterms; ++n) {
        const double np = static_cast<double>(n) * kPi;
        const double sign_n = (n % 2 == 0) ? 1.0 : -1.0;
        // bn = (2/np)*(a*(1-sign_n) - b*sign_n)
        const double bn  = (2.0 / np) * (a * (1.0 - sign_n) - b * sign_n);
        const double lam = np / L;
        sum += bn * std::sin(lam * x) * std::exp(-D * lam * lam * t);
    }
    return P_ss + sum;
}

/**
 * @brief Theis well function W(u) = -Ei(-u), via series expansion.
 * W(u) = -gamma - ln(u) + u - u^2/8 + ...  for u < 1
 * For u >= 10: use exponential integral upper bound approximation
 */
static double TheisW(double u) {
    if (u <= 0.0) return 1.0e30;
    if (u >= 60.0) return 0.0;

    // Use series expansion for all u (converges well for u<10)
    // W(u) = -Ei(-u) = ∫_u^∞ exp(-t)/t dt
    // For u<1: power series is accurate
    // For u>=1: use asymptotic series (careful)
    if (u < 10.0) {
        // Series: W(u) = -γ - ln(u) + u - u²/(2·2!) + u³/(3·3!) - ...
        const double gamma = 0.5772156649015328606;
        double sum = 0.0;
        double term = 1.0;
        double sign = 1.0;
        for (int k = 1; k <= 50; ++k) {
            term *= -u / static_cast<double>(k);
            sum += sign * term / static_cast<double>(k);
            sign = -sign;
        }
        // Actually simpler formula: W(u) = -γ - ln(u) - Σ_{k=1}^∞ (-u)^k / (k·k!)
        // Re-implement cleanly:
        double s2 = 0.0;
        double uk = u;
        for (int k = 1; k <= 80; ++k) {
            double fac = 1.0;
            for (int j = 1; j <= k; ++j) fac *= j;
            double contrib = (k % 2 == 0 ? 1.0 : -1.0) * uk / (k * fac);
            s2 += contrib;
            uk *= u;
            if (std::abs(contrib) < 1e-15 * std::abs(s2)) break;
        }
        return -gamma - std::log(u) - s2;
    }
    else {
        // Asymptotic series: W(u) ≈ exp(-u)/u * (1 - 1/u + 2/u² - 6/u³ + ...)
        // Only first few terms are useful for moderate u
        double emu = std::exp(-u);
        if (emu < 1e-30) return 0.0;
        double s = 1.0;
        double term2 = 1.0;
        for (int k = 1; k <= 8; ++k) {
            term2 *= -static_cast<double>(k) / u;
            if (std::abs(term2) > std::abs(s)) break;
            s += term2;
        }
        return emu / u * s;
    }
}

// Corey relative permeability (for analytical BL)
static double Corey_krw(double Sw) {
    double Se = std::max(0.0, std::min(1.0, (Sw - kSwi) / (1.0 - kSwi - kSor)));
    return kKrw_max * std::pow(Se, kNw);
}
static double Corey_kro(double Sw) {
    double Se = std::max(0.0, std::min(1.0, (1.0 - Sw - kSor) / (1.0 - kSwi - kSor)));
    return kKro_max * std::pow(Se, kNo);
}
static double BL_fw(double Sw) {
    double lam_w = Corey_krw(Sw) / kMu_w;
    double lam_o = Corey_kro(Sw) / kMu_o;
    double denom = lam_w + lam_o;
    return (denom > 1e-30) ? lam_w / denom : 0.0;
}
static double BL_dfw_dSw(double Sw) {
    const double eps = 1e-5;
    return (BL_fw(Sw + eps) - BL_fw(Sw - eps)) / (2.0 * eps);
}

/**
 * @brief Welge tangent method: find frontal saturation Sw_f.
 * Returns Sw_f such that tangent from (Swi, fw(Swi)) touches fw curve.
 */
static double BL_FrontalSw() {
    const double fw_swi = BL_fw(kSwi);
    const int N = 2000;
    double Sw_f = kSwi + 0.01;
    double min_err = 1e30;
    for (int i = 1; i <= N; ++i) {
        double Sw = kSwi + (1.0 - kSwi - kSor - 1e-5) * i / N;
        double fw_sw = BL_fw(Sw);
        double chord = (Sw > kSwi + 1e-10) ? (fw_sw - fw_swi) / (Sw - kSwi) : 0.0;
        double slope = BL_dfw_dSw(Sw);
        double err = std::abs(slope - chord);
        if (err < min_err) { min_err = err; Sw_f = Sw; }
    }
    return Sw_f;
}

// ============================================================
// Error computation helpers
// ============================================================
struct ErrorMetrics {
    double L_inf = 0.0;
    double L2    = 0.0;
};

static ErrorMetrics ComputeError(const std::vector<double>& num,
                                  const std::vector<double>& ana,
                                  double scale) {
    if (num.size() != ana.size() || num.empty()) return {};
    ErrorMetrics em;
    double sum2 = 0.0;
    for (size_t i = 0; i < num.size(); ++i) {
        double e = std::abs(num[i] - ana[i]) / std::max(std::abs(scale), 1e-30);
        em.L_inf = std::max(em.L_inf, e);
        sum2 += e * e;
    }
    em.L2 = std::sqrt(sum2 / static_cast<double>(num.size()));
    return em;
}

static void PrintRelError(const std::string& tag,
                           const ErrorMetrics& em,
                           double tol_warn)
{
    const bool over = (em.L_inf > tol_warn);
    std::cout << "[VAL] " << tag
              << "  L_inf=" << std::scientific << std::setprecision(3) << em.L_inf
              << "  L2="    << std::scientific << std::setprecision(3) << em.L2
              << "  tol="   << tol_warn;
    if (over) std::cout << "  *** WARNING: L_inf > tol ***";
    std::cout << "\n";
}

static void PrintStability(const std::string& tag, bool ok) {
    std::cout << "[VAL] " << tag << (ok ? "  PASS (converged)" : "  *** WARNING: run failed ***") << "\n";
}

// ============================================================
// Field reading helpers
// ============================================================
static std::vector<double> ReadField(FieldManager_2D& fm,
                                      const std::string& name,
                                      int nCells) {
    auto f = fm.getOrCreateMatrixScalar(name, 0.0);
    if (!f || static_cast<int>(f->data.size()) < nCells) return {};
    return std::vector<double>(f->data.begin(), f->data.begin() + nCells);
}

// ============================================================
// Mesh / topology builders
// ============================================================
static void BuildMatrixOnly(MeshManager& mgr) {
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
}

static void AddSingleFrac(MeshManager& mgr, double Lx, double Ly) {
    // Diagonal fracture from (0.2Lx, 0.1Ly) to (0.8Lx, 0.9Ly)
    mgr.addFracture(Vector(0.2 * Lx, 0.1 * Ly, 0.0),
                    Vector(0.8 * Lx, 0.9 * Ly, 0.0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
}

static void AddCrossFrac(MeshManager& mgr, double Lx, double Ly) {
    // Horizontal frac: (0.1Lx, 0.5Ly) -> (0.9Lx, 0.5Ly)
    mgr.addFracture(Vector(0.1 * Lx, 0.5 * Ly, 0.0),
                    Vector(0.9 * Lx, 0.5 * Ly, 0.0));
    // Vertical frac: (0.5Lx, 0.1Ly) -> (0.5Lx, 0.9Ly)
    mgr.addFracture(Vector(0.5 * Lx, 0.1 * Ly, 0.0),
                    Vector(0.5 * Lx, 0.9 * Ly, 0.0));
    mgr.DetectAndSubdivideFractures(IntersectionSearchStrategy_2D::GridIndexing_BasedOn8DOP_DDA);
    mgr.BuildGlobalSystemIndexing();
    mgr.BuildFracturetoFractureTopology();
}

// ============================================================
// Preset builders
// ============================================================
enum class FracConfig { None, Single, Cross };
enum class FluidConfig { PressureOnly, ConstantWater, VarWater };

static FIM_CaseKit::PropertyPreset2D BuildPreset_SP(
    double phi, double k, double c_r,
    double rho_r, double cp_r, double k_r,
    FluidConfig fc)
{
    auto p = FIM_CaseKit::MakeDefaultPropertyPreset2D();
    p.enable_rock_region = false;
    p.rock_bg.phi_r  = phi;
    p.rock_bg.kxx    = k;
    p.rock_bg.kyy    = k;
    p.rock_bg.kzz    = k;
    p.rock_bg.compressibility = c_r;
    p.rock_bg.rho_r  = rho_r;
    p.rock_bg.cp_r   = cp_r;
    p.rock_bg.k_r    = k_r;
    p.frac.permeability   = 1.0e-11;
    p.frac.aperture       = 1.0e-3;
    p.frac.phi_f          = 0.6;
    p.frac.compressibility = c_r;
    p.frac.rho_f = rho_r;
    p.frac.cp_f  = cp_r;
    p.frac.k_f   = k_r;
    if (fc == FluidConfig::PressureOnly)
        p.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::ConstantWaterNoConvection;
    else if (fc == FluidConfig::ConstantWater)
        p.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::ConstantWater;
    else
        p.single_phase_fluid = FIM_Engine::SinglePhaseFluidModel::Water;
    return p;
}

static FIM_CaseKit::PropertyPreset2D BuildPreset_TP(
    double phi, double k, double c_r,
    double rho_r, double cp_r, double k_r)
{
    // Two-phase: always Water+CO2 EOS for N=3
    auto p = FIM_CaseKit::MakeDefaultPropertyPreset2D();
    p.enable_rock_region = false;
    p.rock_bg.phi_r  = phi;
    p.rock_bg.kxx    = k;
    p.rock_bg.kyy    = k;
    p.rock_bg.kzz    = k;
    p.rock_bg.compressibility = c_r;
    p.rock_bg.rho_r  = rho_r;
    p.rock_bg.cp_r   = cp_r;
    p.rock_bg.k_r    = k_r;
    p.frac.permeability   = 1.0e-11;
    p.frac.aperture       = 1.0e-3;
    p.frac.phi_f          = 0.6;
    p.frac.compressibility = c_r;
    p.frac.rho_f = rho_r;
    p.frac.cp_f  = cp_r;
    p.frac.k_f   = k_r;
    // vG params: Swr=0.2, Sgr=0.2 matching Swi/Sor
    p.vg.Swr    = kSwi;
    p.vg.Sgr    = kSor;
    p.vg.n      = 2.0;
    p.vg.alpha  = 5.0e-6;
    p.vg.Pc_max = 1.0e6;
    p.vg.Se_eps = 1.0e-3;
    p.rp.L      = 0.5;
    return p;
}

// ============================================================
// Solver parameters builder for validation
// ============================================================
static FIM_Engine::TransientSolverParams BuildValParams(
    bool twoPhase, double dt_init, double dt_max,
    double t_final, int max_steps)
{
    auto params = FIM_CaseKit::BuildSolverParams(twoPhase, max_steps, dt_init);
    params.gravity_vector = Vector(0.0, 0.0, 0.0);
    params.enable_non_orthogonal_correction = true;
    params.target_end_time_s = t_final;
    params.dt_max    = dt_max;
    params.dt_min    = dt_init * 0.01;
    params.diag_level = FIM_Engine::DiagLevel::Off;
    params.enable_ls_trace = false;
    params.lin_solver = twoPhase
        ? FIM_Engine::LinearSolverType::AMGCL_CPR
        : FIM_Engine::LinearSolverType::AMGCL;
    if (twoPhase) {
        params.max_dSw = 0.08;
        params.max_newton_iter = 20;
    }
    else {
        params.max_newton_iter = 14;
    }
    return params;
}

// ============================================================
// T1: No-well single-phase pressure diffusion
//     N=2, 40×4 grid, 400m×40m
//     ConstantWater model, T frozen at T_init
// ============================================================
struct T1Config {
    std::string tag;
    FracConfig  frac;
    FluidConfig fluid;
    bool        compare_analytic;
};

static void RunT1(const T1Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 400.0, Ly = 40.0;
    const int nx = 40, ny = 4;
    const double phi  = 0.10;
    const double k    = 1.0e-13;   // 100 mD
    const double c_r  = kCt_w;    // rock compressibility = fluid compressibility (ConstantWater)
    const double P0   = 10.0e6;
    const double P_L  = 12.0e6;
    const double P_R  = 8.0e6;
    const double T0   = 380.0;
    const double t_final = 2.0e5;
    const double D_h  = k / (kMu_w * phi * kCt_w);

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_SP(phi, k, c_r, kRho_r, kCp_r, kLam_r, cfg.fluid);

    BoundarySetting::BoundaryConditionManager bcP, bcT;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P_L);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P_R);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP,    0.0);

    // Freeze T: all sides Dirichlet at T0
    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = 1.0;

    auto params = BuildValParams(false, 100.0, 2.0e4, t_final, 50000);
    // ConstantWaterNoConvection: P equation is linear (const rho, const mu).
    // PTC with row_floor=1.0 adds a fictitious 1 kg/(Pa·s) diagonal (real Jacobian ~4e-7),
    // making dx[P] = -b/ptc ≈ 1.5 Pa (negligible) instead of the correct ~7.5 MPa update.
    // Disable PTC so SparseLU solves the unmodified linear system → 1-step convergence.
    params.enable_ptc = false;
    params.lin_solver = FIM_Engine::LinearSolverType::SparseLU;

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<2>(
            cfg.tag, mgr, fm, ic, {}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const auto& cells = mgr.mesh().getCells();
        const int nMat = mgr.getMatrixDOFCount();
        auto p_num = ReadField(fm, "P", nMat);
        if (p_num.empty()) {
            std::cout << "[VAL] " << cfg.tag << " WARNING: could not read P field.\n";
            return;
        }

        const double DP = P_L - P_R;
        std::vector<double> p_ana(nMat);
        for (int i = 0; i < nMat; ++i) {
            p_ana[i] = Analytic1D(cells[i].center.m_x, t_final,
                                   Lx, P_L, P_R, P0, D_h, 100);
        }
        auto em = ComputeError(p_num, p_ana, DP);
        PrintRelError(cfg.tag + " P_error/ΔP", em, 0.03);

        // Write CSV
        std::ofstream csv(caseDir + "/analytical_compare.csv");
        csv << "cell_id,x,y,p_num,p_ana,rel_err\n";
        for (int i = 0; i < nMat; ++i) {
            double e = std::abs(p_num[i] - p_ana[i]) / std::max(std::abs(DP), 1.0);
            csv << i << "," << cells[i].center.m_x << "," << cells[i].center.m_y
                << "," << p_num[i] << "," << p_ana[i] << "," << e << "\n";
        }
    }
}


// ============================================================
// T2: No-well single-phase thermal diffusion
//     N=2, k=1e-20 (no-flow), 40×4 grid, 400m×40m
//     ConstantWater model
// ============================================================
struct T2Config {
    std::string tag;
    FracConfig  frac;
    FluidConfig fluid;
    bool        compare_analytic;
};

static void RunT2(const T2Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 400.0, Ly = 40.0;
    const int nx = 40, ny = 4;
    const double phi  = 0.10;
    const double k    = 1.0e-13;   // still need some k to avoid solver issues
    const double c_r  = kCt_w;
    const double P0   = 10.0e6;
    const double T0   = 380.0;
    const double T_L  = 400.0;
    const double T_R  = 360.0;
    const double t_final = 1.0e7;  // ~115 days

    // Effective thermal diffusivity
    const double lam_eff = phi * kLam_w + (1.0 - phi) * kLam_r;
    const double C_eff   = phi * kRho_w * kCp_w + (1.0 - phi) * kRho_r * kCp_r;
    const double D_T     = lam_eff / C_eff;

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_SP(phi, k, c_r, kRho_r, kCp_r, kLam_r, cfg.fluid);

    BoundarySetting::BoundaryConditionManager bcP, bcT;
    // Same pressure left/right = no Darcy flow
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP,    0.0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T_L);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T_R);
    bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcT.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = 1.0;

    auto params = BuildValParams(false, 1000.0, 5.0e5, t_final, 100000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<2>(
            cfg.tag, mgr, fm, ic, {}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const auto& cells = mgr.mesh().getCells();
        const int nMat = mgr.getMatrixDOFCount();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        auto t_num = ReadField(fm, tCfg.temperatue_field, nMat);
        if (t_num.empty()) {
            std::cout << "[VAL] " << cfg.tag << " WARNING: could not read T field.\n";
            return;
        }

        const double DT = T_L - T_R;
        std::vector<double> t_ana(nMat);
        for (int i = 0; i < nMat; ++i) {
            t_ana[i] = Analytic1D(cells[i].center.m_x, t_final,
                                   Lx, T_L, T_R, T0, D_T, 100);
        }
        auto em = ComputeError(t_num, t_ana, DT);
        PrintRelError(cfg.tag + " T_error/ΔT", em, 0.03);

        std::ofstream csv(caseDir + "/analytical_compare.csv");
        csv << "cell_id,x,y,t_num,t_ana,rel_err\n";
        for (int i = 0; i < nMat; ++i) {
            double e = std::abs(t_num[i] - t_ana[i]) / std::max(std::abs(DT), 1.0);
            csv << i << "," << cells[i].center.m_x << "," << cells[i].center.m_y
                << "," << t_num[i] << "," << t_ana[i] << "," << e << "\n";
        }
    }
}


// ============================================================
// T3: No-well two-phase Buckley-Leverett
//     N=3, 40×4 grid, 400m×40m
// ============================================================
struct T3Config {
    std::string tag;
    FracConfig  frac;
    bool        compare_analytic;
};

static void RunT3(const T3Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 400.0, Ly = 40.0;
    const int nx = 40, ny = 4;
    const double phi = 0.20;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P_L = 13.0e6, P_R = 11.0e6, P0 = 12.0e6;
    const double T0  = 380.0;
    const double Sw0 = kSwi;
    const double t_final = 5.0e6;

    // Compute BL frontal saturation and position
    const double Sw_f = BL_FrontalSw();
    const double fw_f = BL_fw(Sw_f);
    const double dfw_f = BL_dfw_dSw(Sw_f);
    // Total flux approx: q_darcy = k * A * ΔP / (mu_w * L)
    const double A_cross = Ly * 1.0;  // unit depth
    const double q_vol   = k * A_cross * (P_L - P_R) / (kMu_w * Lx);  // rough average
    const double x_f_ana = dfw_f * q_vol * t_final / (phi * A_cross);

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(3);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_TP(phi, k, c_r, kRho_r, kCp_r, kLam_r);

    BoundarySetting::BoundaryConditionManager bcP, bcT, bcS;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P_L);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P_R);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP,    0.0);

    // Freeze T
    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    bcS.Clear();
    bcS.SetDirichletBC(MeshTags::LEFT,  1.0);   // inject water Sw=1
    bcS.SetDirichletBC(MeshTags::RIGHT, Sw0);
    bcS.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcS.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcS);

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = Sw0;

    auto params = BuildValParams(true, 1000.0, 1.0e5, t_final, 200000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<3>(
            cfg.tag, mgr, fm, ic, {}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const auto& cells = mgr.mesh().getCells();
        const int nMat = mgr.getMatrixDOFCount();
        auto sw_num = ReadField(fm, "S_w", nMat);
        if (sw_num.empty()) {
            std::cout << "[VAL] " << cfg.tag << " WARNING: could not read S_w field.\n";
            return;
        }

        // Find numerical frontal position: where Sw transitions from high to low
        // Use x where Sw = (1.0 + Sw0) / 2
        double Sw_mid = 0.5 * (1.0 + Sw0);
        double x_f_num = -1.0;
        for (int i = 0; i < nMat; ++i) {
            if (sw_num[i] < Sw_mid) {
                x_f_num = cells[i].center.m_x;
                break;
            }
        }

        if (x_f_num > 0.0) {
            double err = std::abs(x_f_num - x_f_ana) / Lx;
            std::cout << "[VAL] " << cfg.tag
                      << "  Sw_f=" << Sw_f << "  x_f_ana=" << x_f_ana
                      << "  x_f_num=" << x_f_num
                      << "  rel_err/L=" << err;
            if (err > 0.05) std::cout << "  *** WARNING: front position error > 5% ***";
            std::cout << "\n";
        }
        else {
            std::cout << "[VAL] " << cfg.tag << "  WARNING: front not detected in domain.\n";
        }
    }
}


// ============================================================
// T4: No-well two-phase + heat (BL + thermal front coupling)
//     N=3, 40×4, 400m×40m
// ============================================================
struct T4Config {
    std::string tag;
    FracConfig  frac;
    bool        compare_analytic;
};

static void RunT4(const T4Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 400.0, Ly = 40.0;
    const int nx = 40, ny = 4;
    const double phi = 0.20;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P_L = 13.0e6, P_R = 11.0e6, P0 = 12.0e6;
    const double T0  = 380.0, T_inj = 350.0;
    const double Sw0 = kSwi;
    const double t_final = 5.0e6;

    // Adiabatic thermal front velocity (piston approximation)
    // v_T = q / (phi * A) * (rho_w * cp_w) / (rho_w * cp_w + (1-phi)/phi * rho_r * cp_r)
    // Total flux q = k * A * ΔP / (mu * L)  [approximate with water mobility]
    const double A_cross = Ly * 1.0;
    const double q_vol   = k * A_cross * (P_L - P_R) / (kMu_w * Lx);
    const double v_q     = q_vol / (phi * A_cross);
    const double Rch     = (kRho_w * kCp_w) /
                           (kRho_w * kCp_w + (1.0 - phi) / phi * kRho_r * kCp_r);
    const double x_T_ana = v_q * Rch * t_final;

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(3);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_TP(phi, k, c_r, kRho_r, kCp_r, kLam_r);

    BoundarySetting::BoundaryConditionManager bcP, bcT, bcS;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P_L);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P_R);
    bcP.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcP.SetNeumannBC(MeshTags::TOP,    0.0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T_inj);   // cold injection
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcT.SetNeumannBC(MeshTags::TOP,    0.0);

    bcS.Clear();
    bcS.SetDirichletBC(MeshTags::LEFT,  1.0);
    bcS.SetDirichletBC(MeshTags::RIGHT, Sw0);
    bcS.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcS.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcS);

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = Sw0;

    auto params = BuildValParams(true, 1000.0, 1.0e5, t_final, 200000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<3>(
            cfg.tag, mgr, fm, ic, {}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const auto& cells = mgr.mesh().getCells();
        const int nMat = mgr.getMatrixDOFCount();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        auto t_num = ReadField(fm, tCfg.temperatue_field, nMat);
        if (t_num.empty()) { std::cout << "[VAL] " << cfg.tag << " WARNING: T field missing.\n"; return; }

        // Find numerical thermal front: where T = (T_inj + T0)/2
        double T_mid = 0.5 * (T_inj + T0);
        double x_T_num = -1.0;
        for (int i = 0; i < nMat; ++i) {
            if (t_num[i] < T_mid) { x_T_num = cells[i].center.m_x; break; }
        }

        if (x_T_num > 0.0) {
            double err = std::abs(x_T_num - x_T_ana) / Lx;
            std::cout << "[VAL] " << cfg.tag
                      << "  x_T_ana=" << x_T_ana << "  x_T_num=" << x_T_num
                      << "  rel_err/L=" << err;
            if (err > 0.10) std::cout << "  *** WARNING: thermal front error > 10% ***";
            std::cout << "\n";
        }
        else {
            std::cout << "[VAL] " << cfg.tag << "  INFO: thermal front not reached domain (ok if x_T_ana > L).\n";
        }
    }
}


// ============================================================
// T5: With-well single-phase Theis (pressure drawdown)
//     N=2, 50×50 grid, 2000m×2000m
//     Rate-controlled production well at center
// ============================================================
struct T5Config {
    std::string tag;
    FracConfig  frac;
    FluidConfig fluid;
    bool        compare_analytic;
};

static void RunT5(const T5Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 2000.0, Ly = 2000.0;
    const int nx = 50, ny = 50;
    const double phi = 0.10;
    const double k   = 1.0e-14;   // 10 mD
    const double c_r = kCt_w;
    const double h   = 1.0;       // unit thickness (2D, Lz=1)
    const double P0  = 30.0e6;
    const double T0  = 380.0;
    const double t_final = 3.0e5;
    const double Q_mass  = 0.1;   // kg/s production rate

    // Theis parameters
    // T_darcy = k * h / mu  [m²/(Pa·s)]
    const double T_darcy = k * h / kMu_w;
    // S = phi * ct * h  (storativity, dimensionless)
    const double S = phi * kCt_w * h;
    const double Q_vol  = Q_mass / kRho_w;  // m³/s

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_SP(phi, k, c_r, kRho_r, kCp_r, kLam_r, cfg.fluid);

    // All boundaries: Dirichlet P=P0 (constant pressure outer boundary)
    BoundarySetting::BoundaryConditionManager bcP, bcT;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetDirichletBC(MeshTags::BOTTOM, P0);
    bcP.SetDirichletBC(MeshTags::TOP,    P0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    // Find center cell for production well
    const auto& cells = mgr.mesh().getCells();
    int well_cell = 0;
    double best_d2 = 1e30;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        double dx = cells[i].center.m_x - 0.5 * Lx;
        double dy = cells[i].center.m_y - 0.5 * Ly;
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; well_cell = i; }
    }

    WellScheduleStep prod;
    prod.well_name     = "PROD_T5";
    prod.domain        = WellTargetDomain::Matrix;
    prod.control_mode  = WellControlMode::Rate;
    prod.target_value  = -Q_mass;   // negative = production
    prod.completion_id = well_cell;
    prod.component_mode = WellComponentMode::Water;
    prod.frac_w = 1.0;
    prod.frac_g = 0.0;

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = 1.0;

    auto params = BuildValParams(false, 100.0, 3.0e4, t_final, 50000);
    // Same PTC+AMGCL issues as T1: disable PTC and use SparseLU for exact Newton direction.
    params.enable_ptc = false;
    params.lin_solver = FIM_Engine::LinearSolverType::SparseLU;

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<2>(
            cfg.tag, mgr, fm, ic, {prod}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const int nMat = mgr.getMatrixDOFCount();
        auto p_num = ReadField(fm, "P", nMat);
        if (p_num.empty()) { std::cout << "[VAL] " << cfg.tag << " WARNING: P field missing.\n"; return; }

        const double wx = cells[well_cell].center.m_x;
        const double wy = cells[well_cell].center.m_y;

        // Theis comparison at several radii
        std::array<double,4> r_vals = { 56.0, 112.0, 200.0, 400.0 };
        double max_rel_err = 0.0;
        std::cout << "[VAL] " << cfg.tag << " Theis comparison at t=" << t_final << "s:\n";
        for (double r : r_vals) {
            // Find nearest cell at distance r from well
            int best_i = 0; double best_dr = 1e30;
            for (int i = 0; i < nMat; ++i) {
                double dist = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
                if (std::abs(dist - r) < best_dr) { best_dr = std::abs(dist - r); best_i = i; }
            }
            double r_actual = std::hypot(cells[best_i].center.m_x - wx, cells[best_i].center.m_y - wy);
            double u  = r_actual * r_actual * S / (4.0 * T_darcy * t_final);
            // Drawdown: s(r,t) = Q_vol/(4*pi*T_darcy) * W(u)
            // Note: Q_vol in [m³/s], T_darcy in [m²/(Pa·s)], so:
            // s = Q_vol/(4*pi*T_darcy) * W(u) has units [m³/s * Pa·s/m²] = Pa
            double s_ana = Q_vol / (4.0 * kPi * T_darcy) * TheisW(u);
            double P_ana = P0 - s_ana;
            double dP_max = Q_vol / (4.0 * kPi * T_darcy) * TheisW(r_vals[0] * r_vals[0] * S / (4.0 * T_darcy * t_final));
            double rel_err = (dP_max > 1.0) ? std::abs(p_num[best_i] - P_ana) / dP_max : 0.0;
            max_rel_err = std::max(max_rel_err, rel_err);
            std::cout << "    r=" << r_actual << "  u=" << u << "  W(u)=" << TheisW(u)
                      << "  P_ana=" << P_ana << "  P_num=" << p_num[best_i]
                      << "  rel_err=" << rel_err << "\n";
        }
        if (max_rel_err > 0.05) std::cout << "  *** WARNING: max Theis relative error > 5% ***\n";
    }
}


// ============================================================
// T6: With-well single-phase thermal front
//     N=2, 50×50, 2000m×2000m, injection well at center
// ============================================================
struct T6Config {
    std::string tag;
    FracConfig  frac;
    FluidConfig fluid;
    bool        compare_analytic;
};

static void RunT6(const T6Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 2000.0, Ly = 2000.0;
    const int nx = 50, ny = 50;
    const double phi = 0.10;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P0  = 30.0e6;
    const double T0  = 380.0, T_inj = 340.0;
    const double t_final = 5.0e6;
    const double Q_mass  = 0.1;   // kg/s injection

    const double Q_vol  = Q_mass / kRho_w;
    const double Rch    = kRho_w * kCp_w /
                          (phi * kRho_w * kCp_w + (1.0 - phi) * kRho_r * kCp_r);
    // Thermal front radius: r_T = sqrt(Q_vol * Rch * t / pi)
    const double r_T_ana = std::sqrt(Q_vol * Rch * t_final / kPi);

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(2);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_SP(phi, k, c_r, kRho_r, kCp_r, kLam_r, cfg.fluid);

    BoundarySetting::BoundaryConditionManager bcP, bcT;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetDirichletBC(MeshTags::BOTTOM, P0);
    bcP.SetDirichletBC(MeshTags::TOP,    P0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, nullptr);

    const auto& cells = mgr.mesh().getCells();
    int well_cell = 0; double best_d2 = 1e30;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        double dx = cells[i].center.m_x - 0.5 * Lx;
        double dy = cells[i].center.m_y - 0.5 * Ly;
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; well_cell = i; }
    }

    WellScheduleStep inj;
    inj.well_name          = "INJ_T6";
    inj.domain             = WellTargetDomain::Matrix;
    inj.control_mode       = WellControlMode::Rate;
    inj.target_value       = Q_mass;
    inj.completion_id      = well_cell;
    inj.component_mode     = WellComponentMode::Water;
    inj.frac_w             = 1.0;
    inj.frac_g             = 0.0;
    inj.injection_temperature = T_inj;

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = 1.0;

    auto params = BuildValParams(false, 100.0, 1.0e5, t_final, 100000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<2>(
            cfg.tag, mgr, fm, ic, {inj}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const int nMat = mgr.getMatrixDOFCount();
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        auto t_num = ReadField(fm, tCfg.temperatue_field, nMat);
        if (t_num.empty()) { std::cout << "[VAL] " << cfg.tag << " WARNING: T field missing.\n"; return; }

        const double wx = cells[well_cell].center.m_x;
        const double wy = cells[well_cell].center.m_y;
        double T_mid = 0.5 * (T_inj + T0);

        // Find radius where T = T_mid
        double r_T_num = -1.0;
        double max_r_seen = 0.0;
        for (int i = 0; i < nMat; ++i) {
            double r = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
            if (r > max_r_seen) max_r_seen = r;
            if (t_num[i] < T_mid && r > r_T_num) r_T_num = r;
        }

        std::cout << "[VAL] " << cfg.tag
                  << "  r_T_ana=" << r_T_ana
                  << "  r_T_num=" << r_T_num;
        if (r_T_num > 0.0) {
            double err = std::abs(r_T_num - r_T_ana) / r_T_ana;
            std::cout << "  rel_err=" << err;
            if (err > 0.10) std::cout << "  *** WARNING: thermal radius error > 10% ***";
        }
        std::cout << "\n";
    }
}


// ============================================================
// T7: With-well two-phase radial BL
//     N=3, 50×50, 2000m×2000m, water injection at center
// ============================================================
struct T7Config {
    std::string tag;
    FracConfig  frac;
    bool        compare_analytic;
};

static void RunT7(const T7Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 2000.0, Ly = 2000.0;
    const int nx = 50, ny = 50;
    const double phi = 0.20;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P0  = 25.0e6;
    const double T0  = 380.0;
    const double Sw0 = kSwi;
    const double t_final = 3.0e6;
    const double Q_mass  = 0.1;

    const double Q_vol  = Q_mass / kRho_w;
    // Radial BL front: r_f² ≈ Q_vol * fw(Sw_f) * t / (pi * phi)
    const double Sw_f   = BL_FrontalSw();
    const double fw_f   = BL_fw(Sw_f);
    const double r_f_ana = std::sqrt(Q_vol * fw_f * t_final / (kPi * phi));

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(3);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_TP(phi, k, c_r, kRho_r, kCp_r, kLam_r);

    BoundarySetting::BoundaryConditionManager bcP, bcT, bcS;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetDirichletBC(MeshTags::BOTTOM, P0);
    bcP.SetDirichletBC(MeshTags::TOP,    P0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    bcS.Clear();
    bcS.SetNeumannBC(MeshTags::LEFT,   0.0);
    bcS.SetNeumannBC(MeshTags::RIGHT,  0.0);
    bcS.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcS.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcS);

    const auto& cells = mgr.mesh().getCells();
    int well_cell = 0; double best_d2 = 1e30;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        double dx = cells[i].center.m_x - 0.5 * Lx;
        double dy = cells[i].center.m_y - 0.5 * Ly;
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; well_cell = i; }
    }

    WellScheduleStep inj;
    inj.well_name          = "INJ_T7";
    inj.domain             = WellTargetDomain::Matrix;
    inj.control_mode       = WellControlMode::Rate;
    inj.target_value       = Q_mass;
    inj.completion_id      = well_cell;
    inj.component_mode     = WellComponentMode::Water;
    inj.frac_w             = 1.0;
    inj.frac_g             = 0.0;
    inj.injection_is_co2   = false;

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = Sw0;

    auto params = BuildValParams(true, 1000.0, 1.0e5, t_final, 200000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<3>(
            cfg.tag, mgr, fm, ic, {inj}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const int nMat = mgr.getMatrixDOFCount();
        auto sw_num = ReadField(fm, "S_w", nMat);
        if (sw_num.empty()) { std::cout << "[VAL] " << cfg.tag << " WARNING: S_w field missing.\n"; return; }

        const double wx = cells[well_cell].center.m_x;
        const double wy = cells[well_cell].center.m_y;
        double Sw_mid = 0.5 * (1.0 + Sw0);

        // Find numerical frontal radius
        double r_f_num = -1.0;
        for (int i = 0; i < nMat; ++i) {
            double r = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
            if (sw_num[i] > Sw_mid && r > r_f_num) r_f_num = r;
        }

        std::cout << "[VAL] " << cfg.tag
                  << "  Sw_f=" << Sw_f << "  fw_f=" << fw_f
                  << "  r_f_ana=" << r_f_ana << "  r_f_num=" << r_f_num;
        if (r_f_num > 0.0 && r_f_ana > 0.0) {
            double err = std::abs(r_f_num - r_f_ana) / r_f_ana;
            std::cout << "  rel_err=" << err;
            if (err > 0.10) std::cout << "  *** WARNING: radial BL front error > 10% ***";
        }
        std::cout << "\n";
    }
}


// ============================================================
// T8: With-well two-phase radial BL + heat
//     N=3, 50×50, 2000m×2000m, cold water injection
// ============================================================
struct T8Config {
    std::string tag;
    FracConfig  frac;
    bool        compare_analytic;
};

static void RunT8(const T8Config& cfg) {
    const std::string caseDir = "Test/Transient/ValidationSuite/" + cfg.tag;
    EnsureDir(caseDir);
    std::cout << "\n[ValidationSuite] === " << cfg.tag << " ===\n";

    const double Lx = 2000.0, Ly = 2000.0;
    const int nx = 50, ny = 50;
    const double phi = 0.20;
    const double k   = 1.0e-13;
    const double c_r = kCt_w;
    const double P0  = 25.0e6;
    const double T0  = 380.0, T_inj = 340.0;
    const double Sw0 = kSwi;
    const double t_final = 3.0e6;
    const double Q_mass  = 0.1;

    const double Q_vol  = Q_mass / kRho_w;
    const double Sw_f   = BL_FrontalSw();
    const double fw_f   = BL_fw(Sw_f);
    const double r_f_ana = std::sqrt(Q_vol * fw_f * t_final / (kPi * phi));
    // Thermal front: two-phase heat capacity
    const double C_tp   = phi * (Sw0 * kRho_w * kCp_w + (1.0 - Sw0) * kRho_r * kCp_r * 0.0)
                        + (1.0 - phi) * kRho_r * kCp_r;  // approx Sw=Swi region
    // Simpler: use single-fluid capacity for approximate thermal front
    const double Rch_tp = kRho_w * kCp_w /
                          (phi * kRho_w * kCp_w + (1.0 - phi) * kRho_r * kCp_r);
    const double r_T_ana = std::sqrt(Q_vol * Rch_tp * t_final / kPi);

    MeshManager mgr(Lx, Ly, 0.0, nx, ny, 0, true, false);
    mgr.BuildSolidMatrixGrid_2D();
    if      (cfg.frac == FracConfig::Single) AddSingleFrac(mgr, Lx, Ly);
    else if (cfg.frac == FracConfig::Cross)  AddCrossFrac(mgr, Lx, Ly);
    else                                      BuildMatrixOnly(mgr);
    mgr.setNumDOFs(3);

    FieldManager_2D fm;
    FIM_CaseKit::InitFieldManager(mgr, fm);

    auto preset = BuildPreset_TP(phi, k, c_r, kRho_r, kCp_r, kLam_r);

    BoundarySetting::BoundaryConditionManager bcP, bcT, bcS;
    bcP.Clear();
    bcP.SetDirichletBC(MeshTags::LEFT,   P0);
    bcP.SetDirichletBC(MeshTags::RIGHT,  P0);
    bcP.SetDirichletBC(MeshTags::BOTTOM, P0);
    bcP.SetDirichletBC(MeshTags::TOP,    P0);

    bcT.Clear();
    bcT.SetDirichletBC(MeshTags::LEFT,   T0);
    bcT.SetDirichletBC(MeshTags::RIGHT,  T0);
    bcT.SetDirichletBC(MeshTags::BOTTOM, T0);
    bcT.SetDirichletBC(MeshTags::TOP,    T0);

    bcS.Clear();
    bcS.SetNeumannBC(MeshTags::LEFT,   0.0);
    bcS.SetNeumannBC(MeshTags::RIGHT,  0.0);
    bcS.SetNeumannBC(MeshTags::BOTTOM, 0.0);
    bcS.SetNeumannBC(MeshTags::TOP,    0.0);

    auto modules = FIM_CaseKit::BuildModules2D(preset, &bcP, &bcT, &bcS);

    const auto& cells = mgr.mesh().getCells();
    int well_cell = 0; double best_d2 = 1e30;
    for (int i = 0; i < static_cast<int>(cells.size()); ++i) {
        double dx = cells[i].center.m_x - 0.5 * Lx;
        double dy = cells[i].center.m_y - 0.5 * Ly;
        double d2 = dx * dx + dy * dy;
        if (d2 < best_d2) { best_d2 = d2; well_cell = i; }
    }

    WellScheduleStep inj;
    inj.well_name             = "INJ_T8";
    inj.domain                = WellTargetDomain::Matrix;
    inj.control_mode          = WellControlMode::Rate;
    inj.target_value          = Q_mass;
    inj.completion_id         = well_cell;
    inj.component_mode        = WellComponentMode::Water;
    inj.frac_w                = 1.0;
    inj.frac_g                = 0.0;
    inj.injection_is_co2      = false;
    inj.injection_temperature = T_inj;

    FIM_Engine::InitialConditions ic;
    ic.P_init  = P0;
    ic.T_init  = T0;
    ic.Sw_init = Sw0;

    auto params = BuildValParams(true, 1000.0, 1.0e5, t_final, 200000);

    bool run_ok = true;
    try {
        FIM_Engine::RunGenericFIMTransient<3>(
            cfg.tag, mgr, fm, ic, {inj}, params,
            FIM_Engine::SolverRoute::FIM, modules);
    }
    catch (const std::exception& ex) {
        run_ok = false;
        std::cout << "[VAL] " << cfg.tag << " EXCEPTION: " << ex.what() << "\n";
    }

    if (!run_ok) { PrintStability(cfg.tag, false); return; }
    PrintStability(cfg.tag, true);

    if (cfg.compare_analytic) {
        const int nMat = mgr.getMatrixDOFCount();
        auto sw_num = ReadField(fm, "S_w", nMat);
        const auto tCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        auto t_num = ReadField(fm, tCfg.temperatue_field, nMat);

        const double wx = cells[well_cell].center.m_x;
        const double wy = cells[well_cell].center.m_y;

        // BL front
        if (!sw_num.empty()) {
            double Sw_mid = 0.5 * (1.0 + Sw0);
            double r_f_num = -1.0;
            for (int i = 0; i < nMat; ++i) {
                double r = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
                if (sw_num[i] > Sw_mid && r > r_f_num) r_f_num = r;
            }
            std::cout << "[VAL] " << cfg.tag
                      << "  BL: r_f_ana=" << r_f_ana << "  r_f_num=" << r_f_num;
            if (r_f_num > 0.0 && r_f_ana > 0.0) {
                double err = std::abs(r_f_num - r_f_ana) / r_f_ana;
                std::cout << "  rel_err=" << err;
                if (err > 0.10) std::cout << "  *** WARNING ***";
            }
            std::cout << "\n";
        }

        // Thermal front
        if (!t_num.empty()) {
            double T_mid = 0.5 * (T_inj + T0);
            double r_T_num = -1.0;
            for (int i = 0; i < nMat; ++i) {
                double r = std::hypot(cells[i].center.m_x - wx, cells[i].center.m_y - wy);
                if (t_num[i] < T_mid && r > r_T_num) r_T_num = r;
            }
            std::cout << "[VAL] " << cfg.tag
                      << "  TH: r_T_ana=" << r_T_ana << "  r_T_num=" << r_T_num;
            if (r_T_num > 0.0 && r_T_ana > 0.0) {
                double err = std::abs(r_T_num - r_T_ana) / r_T_ana;
                std::cout << "  rel_err=" << err;
                if (err > 0.15) std::cout << "  *** WARNING ***";
            }
            std::cout << "\n";
        }
    }
}

} // namespace (anonymous)

// ============================================================
// Public wrapper definitions (outside anonymous namespace,
// matching declarations in Test_ValidationSuite.h)
// ============================================================

// T1
void Val_T1_Const_NoFrac()    { RunT1({"Val_T1_Const_NoFrac",    FracConfig::None,   FluidConfig::PressureOnly, true});  }
void Val_T1_Const_SingleFrac(){ RunT1({"Val_T1_Const_SingleFrac",FracConfig::Single, FluidConfig::PressureOnly, false}); }
void Val_T1_Const_CrossFrac() { RunT1({"Val_T1_Const_CrossFrac", FracConfig::Cross,  FluidConfig::PressureOnly, false}); }
void Val_T1_Var_NoFrac()      { RunT1({"Val_T1_Var_NoFrac",      FracConfig::None,   FluidConfig::VarWater,      false}); }
void Val_T1_Var_SingleFrac()  { RunT1({"Val_T1_Var_SingleFrac",  FracConfig::Single, FluidConfig::VarWater,      false}); }
void Val_T1_Var_CrossFrac()   { RunT1({"Val_T1_Var_CrossFrac",   FracConfig::Cross,  FluidConfig::VarWater,      false}); }

// T2
void Val_T2_Const_NoFrac()    { RunT2({"Val_T2_Const_NoFrac",    FracConfig::None,   FluidConfig::ConstantWater, true});  }
void Val_T2_Const_SingleFrac(){ RunT2({"Val_T2_Const_SingleFrac",FracConfig::Single, FluidConfig::ConstantWater, false}); }
void Val_T2_Const_CrossFrac() { RunT2({"Val_T2_Const_CrossFrac", FracConfig::Cross,  FluidConfig::ConstantWater, false}); }
void Val_T2_Var_NoFrac()      { RunT2({"Val_T2_Var_NoFrac",      FracConfig::None,   FluidConfig::VarWater,      false}); }
void Val_T2_Var_SingleFrac()  { RunT2({"Val_T2_Var_SingleFrac",  FracConfig::Single, FluidConfig::VarWater,      false}); }
void Val_T2_Var_CrossFrac()   { RunT2({"Val_T2_Var_CrossFrac",   FracConfig::Cross,  FluidConfig::VarWater,      false}); }

// T3
void Val_T3_Const_NoFrac()    { RunT3({"Val_T3_Const_NoFrac",    FracConfig::None,   true});  }
void Val_T3_Const_SingleFrac(){ RunT3({"Val_T3_Const_SingleFrac",FracConfig::Single, false}); }
void Val_T3_Const_CrossFrac() { RunT3({"Val_T3_Const_CrossFrac", FracConfig::Cross,  false}); }
void Val_T3_Var_NoFrac()      { RunT3({"Val_T3_Var_NoFrac",      FracConfig::None,   false}); }
void Val_T3_Var_SingleFrac()  { RunT3({"Val_T3_Var_SingleFrac",  FracConfig::Single, false}); }
void Val_T3_Var_CrossFrac()   { RunT3({"Val_T3_Var_CrossFrac",   FracConfig::Cross,  false}); }

// T4
void Val_T4_Const_NoFrac()    { RunT4({"Val_T4_Const_NoFrac",    FracConfig::None,   true});  }
void Val_T4_Const_SingleFrac(){ RunT4({"Val_T4_Const_SingleFrac",FracConfig::Single, false}); }
void Val_T4_Const_CrossFrac() { RunT4({"Val_T4_Const_CrossFrac", FracConfig::Cross,  false}); }
void Val_T4_Var_NoFrac()      { RunT4({"Val_T4_Var_NoFrac",      FracConfig::None,   false}); }
void Val_T4_Var_SingleFrac()  { RunT4({"Val_T4_Var_SingleFrac",  FracConfig::Single, false}); }
void Val_T4_Var_CrossFrac()   { RunT4({"Val_T4_Var_CrossFrac",   FracConfig::Cross,  false}); }

// T5
void Val_T5_Const_NoFrac()    { RunT5({"Val_T5_Const_NoFrac",    FracConfig::None,   FluidConfig::PressureOnly, true});  }
void Val_T5_Const_SingleFrac(){ RunT5({"Val_T5_Const_SingleFrac",FracConfig::Single, FluidConfig::PressureOnly, false}); }
void Val_T5_Const_CrossFrac() { RunT5({"Val_T5_Const_CrossFrac", FracConfig::Cross,  FluidConfig::PressureOnly, false}); }
void Val_T5_Var_NoFrac()      { RunT5({"Val_T5_Var_NoFrac",      FracConfig::None,   FluidConfig::VarWater,      false}); }
void Val_T5_Var_SingleFrac()  { RunT5({"Val_T5_Var_SingleFrac",  FracConfig::Single, FluidConfig::VarWater,      false}); }
void Val_T5_Var_CrossFrac()   { RunT5({"Val_T5_Var_CrossFrac",   FracConfig::Cross,  FluidConfig::VarWater,      false}); }

// T6
void Val_T6_Const_NoFrac()    { RunT6({"Val_T6_Const_NoFrac",    FracConfig::None,   FluidConfig::ConstantWater, true});  }
void Val_T6_Const_SingleFrac(){ RunT6({"Val_T6_Const_SingleFrac",FracConfig::Single, FluidConfig::ConstantWater, false}); }
void Val_T6_Const_CrossFrac() { RunT6({"Val_T6_Const_CrossFrac", FracConfig::Cross,  FluidConfig::ConstantWater, false}); }
void Val_T6_Var_NoFrac()      { RunT6({"Val_T6_Var_NoFrac",      FracConfig::None,   FluidConfig::VarWater,      false}); }
void Val_T6_Var_SingleFrac()  { RunT6({"Val_T6_Var_SingleFrac",  FracConfig::Single, FluidConfig::VarWater,      false}); }
void Val_T6_Var_CrossFrac()   { RunT6({"Val_T6_Var_CrossFrac",   FracConfig::Cross,  FluidConfig::VarWater,      false}); }

// T7
void Val_T7_Const_NoFrac()    { RunT7({"Val_T7_Const_NoFrac",    FracConfig::None,   true});  }
void Val_T7_Const_SingleFrac(){ RunT7({"Val_T7_Const_SingleFrac",FracConfig::Single, false}); }
void Val_T7_Const_CrossFrac() { RunT7({"Val_T7_Const_CrossFrac", FracConfig::Cross,  false}); }
void Val_T7_Var_NoFrac()      { RunT7({"Val_T7_Var_NoFrac",      FracConfig::None,   false}); }
void Val_T7_Var_SingleFrac()  { RunT7({"Val_T7_Var_SingleFrac",  FracConfig::Single, false}); }
void Val_T7_Var_CrossFrac()   { RunT7({"Val_T7_Var_CrossFrac",   FracConfig::Cross,  false}); }

// T8
void Val_T8_Const_NoFrac()    { RunT8({"Val_T8_Const_NoFrac",    FracConfig::None,   true});  }
void Val_T8_Const_SingleFrac(){ RunT8({"Val_T8_Const_SingleFrac",FracConfig::Single, false}); }
void Val_T8_Const_CrossFrac() { RunT8({"Val_T8_Const_CrossFrac", FracConfig::Cross,  false}); }
void Val_T8_Var_NoFrac()      { RunT8({"Val_T8_Var_NoFrac",      FracConfig::None,   false}); }
void Val_T8_Var_SingleFrac()  { RunT8({"Val_T8_Var_SingleFrac",  FracConfig::Single, false}); }
void Val_T8_Var_CrossFrac()   { RunT8({"Val_T8_Var_CrossFrac",   FracConfig::Cross,  false}); }

// ============================================================
// Grouped runners
// ============================================================
void Run_All_T1() {
    std::cout << "\n========== ValidationSuite T1: Pressure Diffusion ==========\n";
    Val_T1_Const_NoFrac();    Val_T1_Const_SingleFrac(); Val_T1_Const_CrossFrac();
    Val_T1_Var_NoFrac();      Val_T1_Var_SingleFrac();   Val_T1_Var_CrossFrac();
}
void Run_All_T2() {
    std::cout << "\n========== ValidationSuite T2: Thermal Diffusion ==========\n";
    Val_T2_Const_NoFrac();    Val_T2_Const_SingleFrac(); Val_T2_Const_CrossFrac();
    Val_T2_Var_NoFrac();      Val_T2_Var_SingleFrac();   Val_T2_Var_CrossFrac();
}
void Run_All_T3() {
    std::cout << "\n========== ValidationSuite T3: Buckley-Leverett ==========\n";
    Val_T3_Const_NoFrac();    Val_T3_Const_SingleFrac(); Val_T3_Const_CrossFrac();
    Val_T3_Var_NoFrac();      Val_T3_Var_SingleFrac();   Val_T3_Var_CrossFrac();
}
void Run_All_T4() {
    std::cout << "\n========== ValidationSuite T4: BL + Heat ==========\n";
    Val_T4_Const_NoFrac();    Val_T4_Const_SingleFrac(); Val_T4_Const_CrossFrac();
    Val_T4_Var_NoFrac();      Val_T4_Var_SingleFrac();   Val_T4_Var_CrossFrac();
}
void Run_All_T5() {
    std::cout << "\n========== ValidationSuite T5: Theis Well ==========\n";
    Val_T5_Const_NoFrac();    Val_T5_Const_SingleFrac(); Val_T5_Const_CrossFrac();
    Val_T5_Var_NoFrac();      Val_T5_Var_SingleFrac();   Val_T5_Var_CrossFrac();
}
void Run_All_T6() {
    std::cout << "\n========== ValidationSuite T6: Thermal Front Well ==========\n";
    Val_T6_Const_NoFrac();    Val_T6_Const_SingleFrac(); Val_T6_Const_CrossFrac();
    Val_T6_Var_NoFrac();      Val_T6_Var_SingleFrac();   Val_T6_Var_CrossFrac();
}
void Run_All_T7() {
    std::cout << "\n========== ValidationSuite T7: Radial BL ==========\n";
    Val_T7_Const_NoFrac();    Val_T7_Const_SingleFrac(); Val_T7_Const_CrossFrac();
    Val_T7_Var_NoFrac();      Val_T7_Var_SingleFrac();   Val_T7_Var_CrossFrac();
}
void Run_All_T8() {
    std::cout << "\n========== ValidationSuite T8: Radial BL + Heat ==========\n";
    Val_T8_Const_NoFrac();    Val_T8_Const_SingleFrac(); Val_T8_Const_CrossFrac();
    Val_T8_Var_NoFrac();      Val_T8_Var_SingleFrac();   Val_T8_Var_CrossFrac();
}

void Run_ValidationSuite_All() {
    std::cout << "\n#################### ValidationSuite ALL 48 CASES ####################\n";
    Run_All_T1(); Run_All_T2(); Run_All_T3(); Run_All_T4();
    Run_All_T5(); Run_All_T6(); Run_All_T7(); Run_All_T8();
    std::cout << "\n#################### ValidationSuite COMPLETE ####################\n";
}

} // namespace ValidationSuite
