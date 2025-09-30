#pragma once
#include <iostream>
#include <algorithm>
#include "Mesh.h"
#include "FieldRegistry.h"


// ——确保瞬态计算所需的 *_old 和 *_prev 场存在且尺寸正确 —— //
inline bool ensureTransientFields
(
    Mesh& mesh, FieldRegistry& reg,
    const std::string& p_name = "p_w",
    const std::string& T_name = "T",
    const std::string& p_old_name = "p_w_old",
    const std::string& T_old_name = "T_old",
    const std::string& p_prev_name = "p_w_prev",
    const std::string& T_prev_name = "T_prev"
)

{
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    const size_t n = cells.size();

    // 主场必须已存在
    auto p = reg.get<volScalarField>(p_name);
    auto T = reg.get<volScalarField>(T_name);
   
    if (!p || !T) 
    {
        std::cerr << "[ensureTransientFields] missing primary fields: '"
            << p_name << "' or '" << T_name << "'.\n";
        return false;
    }

    if (p->data.size() != n || T->data.size() != n) 
    {
        std::cerr << "[ensureTransientFields] primary field size mismatch with mesh cells.\n";
        return false;
    }

    // 小工具：获取或创建并确保尺寸；返回是否“新建或刚刚resize”
    auto ensureSized = [&]( const std::string& name,
        std::shared_ptr<volScalarField>& out) -> bool
        {
            out = reg.get<volScalarField>(name);
           
            bool created_or_resized = false;
           
            if (!out) 
            {
                out = reg.getOrCreate<volScalarField>(name, n, 0.0);
                created_or_resized = true;
            }
            
            else if (out->data.size() != n)
            {
                out->data.resize(n, 0.0);
                created_or_resized = true;
            }
            
            return created_or_resized;
        };

    std::shared_ptr<volScalarField> p_old, T_old, p_prev, T_prev;
    const bool need_init_p_old = ensureSized(p_old_name, p_old);
    const bool need_init_T_old = ensureSized(T_old_name, T_old);
    const bool need_init_p_prev = ensureSized(p_prev_name, p_prev);
    const bool need_init_T_prev = ensureSized(T_prev_name, T_prev);

    // 若新建或尺寸变化，用当前 p/T 初始化 *_old 与 *_prev
    if (need_init_p_old || need_init_T_old || need_init_p_prev || need_init_T_prev) {
        for (const auto& c : cells) {
            if (c.id < 0) continue;          // 跳过 ghost
            const size_t i = id2idx.at(c.id);
            if (need_init_p_old)  (*p_old)[i] = (*p)[i];
            if (need_init_T_old)  (*T_old)[i] = (*T)[i];
            if (need_init_p_prev) (*p_prev)[i] = (*p)[i];
            if (need_init_T_prev) (*T_prev)[i] = (*T)[i];
        }
    }

    return true;
}

// 1) 复制标量场：dst <- src
inline bool copyField
(
	FieldRegistry& reg, // 需要确保 src 和 dst 都存在
	const std::string& src_name, //被复制的场
	const std::string& dst_name  //目标场
)
{
    auto src = reg.get<volScalarField>(src_name);
    auto dst = reg.get<volScalarField>(dst_name);
    if (!src || !dst) return false;
    if (dst->data.size() != src->data.size()) dst->data.resize(src->data.size(), 0.0);
	std::copy(src->data.begin(), src->data.end(), dst->data.begin()); // 复制
    return true;
}

// 2) 时间步开始：把当前(p,T)存到 *_old，并用 *_old 初始化 *_prev
inline bool startTimeStep
(
    Mesh& mesh, FieldRegistry& reg,
    const std::string& p_name = "p_w",
    const std::string& T_name = "T",
    const std::string& p_old_name = "p_w_old",
    const std::string& T_old_name = "T_old",
    const std::string& p_prev_name = "p_w_prev",
    const std::string& T_prev_name = "T_prev"
)
{
    if (!ensureTransientFields(mesh, reg, p_name, T_name, p_old_name, T_old_name, p_prev_name, T_prev_name)) return false;

    // p^n, T^n
	if (!copyField(reg, p_name, p_old_name)) return false; //将p_w 复制到_old
	if (!copyField(reg, T_name, T_old_name)) return false; //将T 复制到_old

    // prev ← old (给 k=0 的外迭代作为初值)
	if (!copyField(reg, p_old_name, p_prev_name)) return false; //将p_w_old 复制到 p_w_prev
	if (!copyField(reg, T_old_name, T_prev_name)) return false; //将T_old 复制到 T_prev

    return true;
}


// 3) 外迭代开始：把 prev (k层) 拷到当前工作场 (p,T)
inline bool startOuterIteration
(
    FieldRegistry& reg,
    const std::string& p_name = "p_w",
    const std::string& T_name = "T",
    const std::string& p_prev_name = "p_w_prev",
    const std::string& T_prev_name = "T_prev"
)
{
    bool ok = true;
    ok = ok && copyField(reg, p_prev_name, p_name);
    ok = ok && copyField(reg, T_prev_name, T_name);
    return ok;
}


// 4) 欠松弛：x <- x_prev + alpha * (x - x_prev), 在 x 场内就地更新
inline bool underRelaxInPlace
(
    FieldRegistry& reg,
    const std::string& x_name,
    const std::string& x_prev_name,
    double alpha)  // 0<alpha<=1
{
    auto x = reg.get<volScalarField>(x_name);
    auto xp = reg.get<volScalarField>(x_prev_name);
    if (!x || !xp) return false;
    if (x->data.size() != xp->data.size()) return false;
    alpha = std::max(0.0, std::min(1.0, alpha));
    for (size_t i = 0; i < x->data.size(); ++i) 
    {
        x->data[i] = xp->data[i] + alpha * (x->data[i] - xp->data[i]);
    }
    return true;
}

// 5) 外迭代推进：把当前解写回 prev，作为下一次外迭代的“上一层”
inline bool updatePrevIterates
(
    FieldRegistry& reg,
    const std::string& p_name = "p_w",
    const std::string& T_name = "T",
    const std::string& p_prev_name = "p_w_prev",
    const std::string& T_prev_name = "T_prev")
{
    bool ok = true;
    ok = ok && copyField(reg, p_name, p_prev_name);
    ok = ok && copyField(reg, T_name, T_prev_name);
    return ok;
}

// 6) 收敛判据：∞-范数(最大绝对差)
inline double maxAbsDiff
(
    const FieldRegistry& reg,
    const std::string& a_name,
    const std::string& b_name)
{
    auto a = reg.get<volScalarField>(a_name);
    auto b = reg.get<volScalarField>(b_name);
    if (!a || !b || a->data.size() != b->data.size()) return std::numeric_limits<double>::infinity();
    double m = 0.0;
    for (size_t i = 0; i < a->data.size(); ++i) m = std::max(m, std::abs(a->data[i] - b->data[i]));
    return m;
}

// 7) 收敛判据：L2-范数(均方根)
inline double rmsDiff
(
	const FieldRegistry& reg,
	const std::string& a_name,
	const std::string& b_name)
{
	auto a = reg.get<volScalarField>(a_name);
	auto b = reg.get<volScalarField>(b_name);
	if (!a || !b || a->data.size() != b->data.size() || a->data.empty()) return std::numeric_limits<double>::infinity();
	double sum2 = 0.0;
	for (size_t i = 0; i < a->data.size(); ++i) {
		double d = a->data[i] - b->data[i];
		sum2 += d * d;
	}
	return std::sqrt(sum2 / static_cast<double>(a->data.size()));
}

// 8) 收敛判据：相对误差的 L2-范数(均方根)
inline double maxRelChange(const FieldRegistry& reg,
    const std::string& x_name,
    const std::string& x_prev_name,
    double eps = 1e-20)
{
    auto x = reg.get<volScalarField>(x_name);
    auto xp = reg.get<volScalarField>(x_prev_name);
    if (!x || !xp || x->data.size() != xp->data.size()) return std::numeric_limits<double>::infinity();
    double m = 0.0;
    for (size_t i = 0; i < x->data.size(); ++i) {
        double denom = std::max(std::abs(xp->data[i]), eps);
        m = std::max(m, std::abs(x->data[i] - xp->data[i]) / denom);
    }
    return m;
}


// 9) 时间步结束：把新时层(p,T)提交为下一步的 old
inline bool endTimeStep(FieldRegistry& reg,
    const std::string& p_name = "p_w",
    const std::string& T_name = "T",
    const std::string& p_old_name = "p_w_old",
    const std::string& T_old_name = "T_old")
{
    bool ok = true;
    ok = ok && copyField(reg, p_name, p_old_name);
    ok = ok && copyField(reg, T_name, T_old_name);
    return ok;
}



//*********补充工具函数**********//
// 10) 构建时间线性化点场（用于牛顿迭代的雅可比矩阵组装）
inline bool bulidTimeLinearizationPoint
(
    MeshManager& mgr, FieldRegistry& reg,
    const std::string& p_old_name, const std::string& T_old_name,   // n 层
    const std::string& p_iter_name = "", const std::string& T_iter_name = "", // 可选：迭代/预测层
    const std::string& p_lin_out = "p_time_lin",
    const std::string& T_lin_out = "T_time_lin"
)
{
    auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    if (cells.empty()) return false;

    //取出旧时层
    auto p0 = reg.get<volScalarField>(p_old_name);
    auto T0 = reg.get<volScalarField>(T_old_name);
    if (!p0 || !T0) { std::cerr << "[TimeLin] missing p_old/T_old\n"; return false; }

    //取出迭代层
    std::shared_ptr<volScalarField> pk, Tk;
    if (!p_iter_name.empty()) pk = reg.get<volScalarField>(p_iter_name);
    if (!T_iter_name.empty()) Tk = reg.get<volScalarField>(T_iter_name);

    // 生成线性化点场
    auto p_lin = reg.getOrCreate<volScalarField>(p_lin_out, cells.size(), 0.0);
    auto T_lin = reg.getOrCreate<volScalarField>(T_lin_out, cells.size(), 0.0);

    for (const auto& c : cells)
    {
        if (c.id < 0)continue;
        const size_t i = id2idx.at(c.id);
        double pn = (*p0)[i], Tn = (*T0)[i];
        double pl = pn, Tl = Tn;

        if (pk) {
            const double pik = (*pk)[i];
            const double Tik = Tk ? (*Tk)[i] : Tn; // 若无迭代温度场则用旧时层温度
			pl = 0.5 * (pn + pik);
			Tl = 0.5 * (Tn + Tik);
		}
        Initializer::clampPT(pl, Tl);
        (*p_lin)[i] = pl;
        (*T_lin)[i] = Tl;
    
    }

    return true;
}

// 11)在 (p_lin, T_lin) 处评估 ρ 与 ∂ρ/∂p（数值微分）
inline bool computeRhoAndDrhoDpAt
(
    MeshManager& mgr, FieldRegistry& reg,
    const std::string& p_lin_name = "p_time_lin",
    const std::string& T_lin_name = "T_time_lin",
    const std::string& phase = "water",          // "water" | "co2"
    const std::string& rho_out = "rho_time",
    const std::string& drhodp_out = "drho_dp_time",
    double dp_rel = 1e-4, double dp_abs_min = 10.0
)
{
    auto& mesh = mgr.mesh();
    const auto& cells = mesh.getCells();
    const auto& id2idx = mesh.getCellId2Index();
    if (cells.empty()) return false;

    auto pL = reg.get<volScalarField>(p_lin_name);
    auto TL = reg.get<volScalarField>(T_lin_name);
    if (!pL || !TL) { std::cerr << "[TimeLin] missing p_lin/T_lin\n"; return false; }

    auto rho = reg.getOrCreate<volScalarField>(rho_out, cells.size(), 0.0);
    auto dr = reg.getOrCreate<volScalarField>(drhodp_out, cells.size(), 0.0);

    auto wt = WaterPropertyTable::instance();
    auto gt = CO2PropertyTable::instance();
    const bool isCO2 = (phase == "co2" || phase == "CO2");

    for (const auto& c : cells) {
        if (c.id < 0) continue;
        const size_t i = id2idx.at(c.id);
        double p = (*pL)[i], T = (*TL)[i];
        Initializer::clampPT(p, T);

        double rho_lin = 1000.0;
        try { rho_lin = isCO2 ? gt.getProperties(p, T).rho : wt.getProperties(p, T).rho; }
        catch (...) { /* OOR 兜底 */ }

        const double dpa = std::max(dp_abs_min, std::abs(p) * dp_rel);
        double rp = rho_lin, rm = rho_lin;
        try {
            double pp = p + dpa, pm = p - dpa;
            Initializer::clampPT(pp, T); Initializer::clampPT(pm, T);
            rp = isCO2 ? gt.getProperties(pp, T).rho : wt.getProperties(pp, T).rho;
            rm = isCO2 ? gt.getProperties(pm, T).rho : wt.getProperties(pm, T).rho;
        }
        catch (...) { /* OOR 兜底 */ }

        (*rho)[i] = rho_lin;
        (*dr)[i] = (rp - rm) / std::max(2.0 * dpa, 1e-12);
    }
    return true;
}

