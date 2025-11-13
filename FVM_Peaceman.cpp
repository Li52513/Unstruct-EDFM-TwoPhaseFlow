#include "FVM_Peaceman.h"
#include "FieldAcessForDiscre.h"   // cellScalar(...)
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cassert>

namespace {
    constexpr double kPI = 3.14159265358979323846;
}

// ========= Mesh 访问小工具（适配你的 Mesh 接口） =========
static inline int cell_id_by_index(const Mesh& mesh, int cidx)
{
    const auto& cells = mesh.getCells();
    assert(cidx >= 0 && cidx < (int)cells.size());
    return cells[cidx].id; // 与 getCellId2Index 的 key 对齐
}
static inline double cell_measure_2D(const Mesh& mesh, int cidx) {
    const int cid = cell_id_by_index(mesh, cidx);
    return mesh.getCellArea(cid); // 2D：面积为测度
}
static inline Vector cell_center(const Mesh& mesh, int cidx)
{
    const auto& cells = mesh.getCells();
    const Cell& c = cells[cidx];  // 若你用 c.center()，在此替换
    return Vector(c.center.m_x, c.center.m_y, c.center.m_z);
}

// ========= 从字段估算 k_h =========
static inline double kh_from_fields(const FieldRegistry& reg, const Mesh& mesh,
    int cid, double fallbackKh)
{
    const double kxx = cellScalar(reg, mesh, "kxx", cid, std::numeric_limits<double>::quiet_NaN());
    const double kyy = cellScalar(reg, mesh, "kyy", cid, std::numeric_limits<double>::quiet_NaN());
    if (std::isfinite(kxx) && std::isfinite(kyy) && kxx > 0.0 && kyy > 0.0)
        return std::sqrt(kxx * kyy);

    const double k = cellScalar(reg, mesh, "k", cid, std::numeric_limits<double>::quiet_NaN());
    if (std::isfinite(k) && k > 0.0) return k;

    return fallbackKh;
}

// ========= 选择穿孔单元（半径 or 最近 N）=========
static std::vector<int> pick_perforated_cells(const Mesh& mesh, const WellSpec& spec)
{
    const int Nc = (int)mesh.getCells().size();
    std::vector<int> hit;

    if (spec.perfRadius > 0.0) {
        const double R2 = spec.perfRadius * spec.perfRadius;
        hit.reserve(32);
        for (int i = 0; i < Nc; ++i) {
            const Vector cc = cell_center(mesh, i);
            const double dx = cc.m_x - spec.pos.m_x;
            const double dy = cc.m_y - spec.pos.m_y;
            const double dz = cc.m_z - spec.pos.m_z;
            if (dx * dx + dy * dy + dz * dz <= R2) hit.push_back(i);
        }
        if (!hit.empty()) return hit;
        // 半径内未命中则降级为最近 N
    }

    struct Item { int idx; double d2; };
    std::vector<Item> all; all.reserve(Nc);
    for (int i = 0; i < Nc; ++i) {
        const Vector cc = cell_center(mesh, i);
        const double dx = cc.m_x - spec.pos.m_x;
        const double dy = cc.m_y - spec.pos.m_y;
        const double dz = cc.m_z - spec.pos.m_z;
        all.push_back({ i, dx * dx + dy * dy + dz * dz });
    }
    const int take = std::max(1, std::min(spec.maxHitCells, Nc));
    std::nth_element(all.begin(), all.begin() + (take - 1), all.end(),
        [](const Item& a, const Item& b) { return a.d2 < b.d2; });

    hit.resize(take);
    for (int j = 0; j < take; ++j) hit[j] = all[j].idx;
    return hit;
}

// ========= 通用核心：根据 spec 构建 mask + PI（字段名由 out 传入）=========
static void build_mask_and_PI_impl(
    MeshManager& mgr, FieldRegistry& reg,
    const WellSpec& spec, const PeacemanParams& prm,
    const WellFieldNames& out)
{
    Mesh& mesh = mgr.mesh();
    const int Nc = (int)mesh.getCells().size();

    // 1) 穿孔单元集合
    std::vector<int> perf = pick_perforated_cells(mesh, spec);
    if (perf.empty()) {
        std::cerr << "[Peaceman] No perforated cells for '" << spec.name
            << "'. Fallback to cell 0.\n";
        perf = { 0 };
    }

    // 2) 输出场（不存在则创建；存在则清零）
    auto maskFld = reg.get<volScalarField>(out.mask);
    if (!maskFld) {
        maskFld = reg.create<volScalarField>(out.mask, Nc, 0.0);
    }
    else {
        // VolField<double> 无 begin()/end()，用显式下标清零
        for (int i = 0; i < Nc; ++i) { (*maskFld)[i] = 0.0; }
    }

    auto PIFld = reg.get<volScalarField>(out.PI);
    if (!PIFld) {
        PIFld = reg.create<volScalarField>(out.PI, Nc, 0.0);
    }
    else {
        for (int i = 0; i < Nc; ++i) { (*PIFld)[i] = 0.0; }
    }

    // 3) 逐穿孔单元计算 WI_i 与 PI_i
    //    Peaceman（2D 厚度 H）：
    //      WI_i = 2π * k_h_i * H / ( ln(r_e_i / r_w) + skin )
    //      PI(vol)  = WI / mu
    //      PI(mass) = rho * WI / mu
    //    r_e_i = reFactor * sqrt(A_i/π)
    double WI_total = 0.0;
    for (int cidx : perf) {
        const int cid = cell_id_by_index(mesh, cidx);
        const double Ai = cell_measure_2D(mesh, cidx);
        const double kh = kh_from_fields(reg, mesh, cid, prm.fallbackKh);

        const double re = prm.reFactor * std::sqrt(std::max(Ai, 1e-30) / kPI);
        const double lnre = std::log(std::max(re / spec.rw, 1.0 + 1e-12));
        const double denom = std::max(lnre + spec.skin, 1e-12);

        const double WI_i = (2.0 * kPI) * kh * spec.H / denom;
        double PI_i = WI_i / std::max(prm.mu, 1e-30);
        if (prm.PI_is_mass) PI_i *= std::max(prm.rho, 1e-12);

        (*maskFld)[cidx] = 1.0;
        (*PIFld)[cidx] = PI_i;
        WI_total += WI_i;
    }

    // 4) 诊断
    std::cout << "[Peaceman][" << out.mask << "] perfs=" << perf.size()
        << ", sum(WI)=" << WI_total
        << ", PI_is_mass=" << (prm.PI_is_mass ? "true" : "false")
        << ", reFactor=" << prm.reFactor << "\n";
}

// ========= 对外接口：注入井 =========
void build_injection_mask_and_PI(
    MeshManager& mgr, FieldRegistry& reg,
    const WellSpec& spec, const PeacemanParams& prm)
{
    WellFieldNames out; out.mask = "inj_mask"; out.PI = "PI_inj";
    build_mask_and_PI_impl(mgr, reg, spec, prm, out);
}

// ========= 对外接口：采出井 =========
void build_production_mask_and_PI(
    MeshManager& mgr, FieldRegistry& reg,
    const WellSpec& spec, const PeacemanParams& prm)
{
    WellFieldNames out; out.mask = "prod_mask"; out.PI = "PI_prod";
    build_mask_and_PI_impl(mgr, reg, spec, prm, out);
}


void build_well_mask_and_PI
(
    MeshManager& mgr, FieldRegistry& reg,
    const WellSpec& spec, const PeacemanParams& prm,
    const std::string& mask_name, const std::string& PI_name)
{
    Mesh& mesh = mgr.mesh();
    const int Nc = (int)mesh.getCells().size();

    // 1) 穿孔单元集合
    std::vector<int> perf = pick_perforated_cells(mesh, spec);
    if (perf.empty()) {
        std::cerr << "[Peaceman] No perforated cells for '" << spec.name
            << "'. Fallback to cell 0.\n";
        perf = { 0 };
    }

    // 2) 输出场（不存在则创建；存在则清零）
    auto maskFld = reg.get<volScalarField>(mask_name);
    if (!maskFld) {
        maskFld = reg.create<volScalarField>(mask_name, Nc, 0.0);
    }
    else {
        // VolField<double> 无 begin()/end()，用显式下标清零
        for (int i = 0; i < Nc; ++i) { (*maskFld)[i] = 0.0; }
    }

    auto PIFld = reg.get<volScalarField>(PI_name);
    if (!PIFld) {
        PIFld = reg.create<volScalarField>(PI_name, Nc, 0.0);
    }
    else {
        for (int i = 0; i < Nc; ++i) { (*PIFld)[i] = 0.0; }
    }

    // 3) 逐穿孔单元计算 WI_i 与 PI_i
    //    Peaceman（2D 厚度 H）：
    //      WI_i = 2π * k_h_i * H / ( ln(r_e_i / r_w) + skin )
    //      PI(vol)  = WI / mu
    //      PI(mass) = rho * WI / mu
    //    r_e_i = reFactor * sqrt(A_i/π)
    double WI_total = 0.0;
    for (int cidx : perf) {
        const int cid = cell_id_by_index(mesh, cidx);
        const double Ai = cell_measure_2D(mesh, cidx);
        const double kh = kh_from_fields(reg, mesh, cid, prm.fallbackKh);

        const double re = prm.reFactor * std::sqrt(std::max(Ai, 1e-30) / kPI);
        const double lnre = std::log(std::max(re / spec.rw, 1.0 + 1e-12));
        const double denom = std::max(lnre + spec.skin, 1e-12);

        const double WI_i = (2.0 * kPI) * kh * spec.H / denom;
        double PI_i = WI_i / std::max(prm.mu, 1e-30);
        if (prm.PI_is_mass) PI_i *= std::max(prm.rho, 1e-12);

        (*maskFld)[cidx] = 1.0;
        (*PIFld)[cidx] = PI_i;
        WI_total += WI_i;
    }

    // 4) 诊断
    //std::cout << "[Peaceman][" << mask_name << "] perfs=" << perf.size()
    //    << ", sum(WI)=" << WI_total
    //    << ", PI_is_mass=" << (prm.PI_is_mass ? "true" : "false")
    //    << ", reFactor=" << prm.reFactor << "\n";
}