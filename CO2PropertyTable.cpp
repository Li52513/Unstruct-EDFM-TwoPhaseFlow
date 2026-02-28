#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include "CO2PropertyTable.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <experimental/filesystem>

namespace {
    /**
     * @brief 边界钳位函数 (Boundary Clamping)
     * @details 保证查表在到达物性表边缘时，索引不会越界。
     */
    inline size_t clamp_idx(int idx, int max_idx) {
        if (idx < 0) return 0;
        if (idx > max_idx) return static_cast<size_t>(max_idx);
        return static_cast<size_t>(idx);
    }

    /**
     * @brief 在升序数组中定位局部区间及归一化坐标 frac
     */
    inline void bracket2(const std::vector<double>& axis, double x,
        size_t& i0, size_t& i1, double& frac)
    {
        auto it1 = std::upper_bound(axis.begin(), axis.end(), x);
        size_t idx1 = size_t(it1 - axis.begin());
        if (idx1 == 0) idx1 = 1;
        if (idx1 >= axis.size()) idx1 = axis.size() - 1;
        i1 = idx1;  i0 = idx1 - 1;

        double x0 = axis[i0], x1 = axis[i1];
        frac = (x1 > x0) ? ((x - x0) / (x1 - x0)) : 0.0;
        if (frac < 0.0) frac = 0.0; else if (frac > 1.0) frac = 1.0;
    }

    /**
     * @brief 1D Catmull-Rom 三次样条插值核心
     * @param p0, p1, p2, p3 连续四个网格点的值 (p1, p2 为目标区间端点)
     * @param t 归一化局部坐标 t in [0, 1]
     */
    inline double catmull_rom_1d(double p0, double p1, double p2, double p3, double t) {
        double t2 = t * t;
        double t3 = t2 * t;
        return 0.5 * ((2.0 * p1) +
            (-p0 + p2) * t +
            (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2 +
            (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3);
    }

    /**
     * @brief 1D Catmull-Rom 三次样条一阶导数 (对归一化坐标 t 求导)
     */
    inline double catmull_rom_1d_deriv_t(double p0, double p1, double p2, double p3, double t) {
        double t2 = t * t;
        return 0.5 * ((-p0 + p2) +
            2.0 * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t +
            3.0 * (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t2);
    }

    /**
     * @brief 2D Catmull-Rom 双三次样条插值 (物理空间)
     * @param v 4x4 数值矩阵模板
     * @param a 第一维(Pressure)局部坐标 [0,1]
     * @param b 第二维(Temperature)局部坐标 [0,1]
     */
    inline double bicubic_spline_linear(const double v[4][4], double a, double b) {
        double arr[4];
        for (int i = 0; i < 4; ++i) {
            arr[i] = catmull_rom_1d(v[i][0], v[i][1], v[i][2], v[i][3], b); // 先沿 T 方向插值
        }
        return catmull_rom_1d(arr[0], arr[1], arr[2], arr[3], a); // 再沿 P 方向插值
    }

    /**
     * @brief 2D Catmull-Rom 双三次样条插值 (对数空间防负值机制)
     * @details 在 ln 空间内执行高阶插值，百分百避免龙格现象导致密度、黏度等物理量出现非物理的负值。
     */
    inline double bicubic_spline_logpos(const double v[4][4], double a, double b) {
        double v_log[4][4];
        bool all_positive = true;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (v[i][j] <= 0.0) {
                    all_positive = false; break;
                }
                v_log[i][j] = std::log(v[i][j]);
            }
            if (!all_positive) break;
        }
        // 正常状态在对数空间内计算
        if (all_positive) {
            double L = bicubic_spline_linear(v_log, a, b);
            return std::exp(L);
        }
        // 退化机制：如果含有0或负数，回退到普通物理空间插值
        return bicubic_spline_linear(v, a, b);
    }

    /**
     * @brief 2D Catmull-Rom 对数空间关于第一维局部坐标(a)的一阶偏导数
     * @return 返回 d(ln_Value) / da
     */
    inline double bicubic_spline_logpos_deriv_a(const double v[4][4], double a, double b) {
        double v_log[4][4];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (v[i][j] <= 0.0) return 0.0; // 越界触发安全兜底，返回 0 梯度
                v_log[i][j] = std::log(v[i][j]);
            }
        }
        double arr[4];
        for (int i = 0; i < 4; ++i) {
            arr[i] = catmull_rom_1d(v_log[i][0], v_log[i][1], v_log[i][2], v_log[i][3], b);
        }
        return catmull_rom_1d_deriv_t(arr[0], arr[1], arr[2], arr[3], a);
    }
}
namespace fs = std::experimental::filesystem;

CO2PropertyTable& CO2PropertyTable::instance() 
{
    static CO2PropertyTable inst("D:/dataBase_MatrialProperties/CO2_SpanWagner.txt");
    return inst;
}


// 构造函数：自动加载
CO2PropertyTable::CO2PropertyTable(const std::string& filename) 
{
    load(filename);
}

// ——— load(): 支持文本解析 + .cache 二进制缓存 ———
void CO2PropertyTable::load(const std::string& filename) {
    const std::string cacheFile = filename + ".cache";
    bool useCache = false;

    if (fs::exists(cacheFile) && fs::exists(filename)) {
        if (fs::last_write_time(cacheFile) > fs::last_write_time(filename)) {
            useCache = true;
        }
    }
    if (useCache)
    {
        std::ifstream binIn(cacheFile, std::ios::binary);
        if (binIn) {
            size_t np, nt;
            binIn.read(reinterpret_cast<char*>(&np), sizeof(np));
            pressures_.resize(np);
            binIn.read(reinterpret_cast<char*>(pressures_.data()), np * sizeof(double));

            binIn.read(reinterpret_cast<char*>(&nt), sizeof(nt));
            temps_.resize(nt);
            binIn.read(reinterpret_cast<char*>(temps_.data()), nt * sizeof(double));

            data_.assign(np, std::vector<CO2Properties>(nt));
            for (size_t i = 0; i < np; ++i) {
                binIn.read(reinterpret_cast<char*>(data_[i].data()), nt * sizeof(CO2Properties));
            }
            if (binIn) return;
        }
    }
    // 2) 文本解析
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open CO2 property file: " + filename);

    std::string line;
    // 跳过表头行，直到含有 “Density” 和 “Viscosity”
    while (std::getline(in, line)) {
        if (line.find("Density") != std::string::npos &&
            line.find("Viscosity") != std::string::npos)
        {
            break;
        }
    }

    struct Row { double P{}, T{}; CO2Properties w{}; };
    std::vector<Row> rows;

    // 逐行读取：P  T  rho  mu  cp  cv  h  k
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        Row r{};  // P,T 自动归零
        if (iss >> r.P >> r.T
            >> r.w.rho  // 密度
            >> r.w.mu   // 粘度
            >> r.w.cp   // 定压比热
            >> r.w.cv   // 定容比热
            >> r.w.h    // 比焓
            >> r.w.k)   // 导热系数
        {
            r.w.c = 0.0;
            r.w.dRho_dP = 0.0;
            rows.push_back(r);
        }
    }

    // 构造升序唯一 pressures_
    {
        std::vector<double> tmp;
        tmp.reserve(rows.size());
        for (auto& r : rows) tmp.push_back(r.P);
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
        pressures_ = std::move(tmp);
    }
    // 构造升序唯一 temps_（基于第一个压力值 P0）
    {
        double P0 = pressures_.front();
        std::vector<double> tmp;
        tmp.reserve(rows.size());
        for (auto& r : rows)
            if (r.P == P0) tmp.push_back(r.T);
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
        temps_ = std::move(tmp);
    }

    // 填充 data_
    size_t np = pressures_.size(), nt = temps_.size();
    data_.assign(np, std::vector<CO2Properties>(nt));
    for (auto& r : rows) {
        size_t iP = std::lower_bound(pressures_.begin(), pressures_.end(), r.P) - pressures_.begin();
        size_t iT = std::lower_bound(temps_.begin(), temps_.end(), r.T) - temps_.begin();
        data_[iP][iT] = r.w;
    }

    // 3) 写入二进制缓存
    {
        std::ofstream binOut(cacheFile, std::ios::binary);
        binOut.write(reinterpret_cast<char*>(&np), sizeof(np));
        binOut.write(reinterpret_cast<char*>(pressures_.data()), np * sizeof(double));
        binOut.write(reinterpret_cast<char*>(&nt), sizeof(nt));
        binOut.write(reinterpret_cast<char*>(temps_.data()), nt * sizeof(double));
        for (size_t i = 0; i < np; ++i) {
            for (size_t j = 0; j < nt; ++j) {
                auto& w = data_[i][j];
                binOut.write(reinterpret_cast<char*>(&w), sizeof(w));
            }
        }
    }
}

CO2Properties CO2PropertyTable::getProperties(double P, double T) const
{
    if (P < pressures_.front() || P > pressures_.back()
        || T < temps_.front() || T > temps_.back())
        throw std::out_of_range("CO2PropertyTable: (P,T) outside bounds");

    size_t iP0, iP1, iT0, iT1; double a, b;
    bracket2(pressures_, P, iP0, iP1, a);
    bracket2(temps_, T, iT0, iT1, b);

    int nP = static_cast<int>(pressures_.size());
    int nT = static_cast<int>(temps_.size());
    int ip_base = static_cast<int>(iP0);
    int it_base = static_cast<int>(iT0);

    // 构建 4x4 局部网格模板数据 (Stencil)
    double v_rho[4][4], v_mu[4][4], v_k[4][4];
    double v_cp[4][4], v_cv[4][4], v_h[4][4];

    for (int i = -1; i <= 2; ++i) {
        size_t idx_P = clamp_idx(ip_base + i, nP - 1); // 边界钳位保护
        for (int j = -1; j <= 2; ++j) {
            size_t idx_T = clamp_idx(it_base + j, nT - 1);
            const auto& w = data_[idx_P][idx_T];
            v_rho[i + 1][j + 1] = w.rho;
            v_mu[i + 1][j + 1] = w.mu;
            v_k[i + 1][j + 1] = w.k;
            v_cp[i + 1][j + 1] = w.cp;
            v_cv[i + 1][j + 1] = w.cv;
            v_h[i + 1][j + 1] = w.h;
        }
    }

    CO2Properties R;
    // 严格为正物理量应用对数空间双三次样条插值
    R.rho = bicubic_spline_logpos(v_rho, a, b);
    R.mu = bicubic_spline_logpos(v_mu, a, b);
    R.k = bicubic_spline_logpos(v_k, a, b);

    // 热力学物性采用物理空间双三次样条插值
    R.cp = bicubic_spline_linear(v_cp, a, b);
    R.cv = bicubic_spline_linear(v_cv, a, b);
    R.h = bicubic_spline_linear(v_h, a, b);

    // 动态计算 dRho/dP 和压缩系数 c (严密的链式法则)
    double deltaP_grid = pressures_[iP1] - pressures_[iP0];
    if (deltaP_grid > 1e-9) {
        // L = ln(rho)，求偏导： dL/dP = (dL/da) * (da/dP) = (dL/da) / deltaP_grid
        double dL_da = bicubic_spline_logpos_deriv_a(v_rho, a, b);
        double dL_dP = dL_da / deltaP_grid;

        R.c = dL_dP; // 因 c = 1/rho * d(rho)/dP = d(ln_rho)/dP
        R.dRho_dP = R.rho * R.c;
    }
    else {
        R.c = 0.0;
        R.dRho_dP = 0.0;
    }

    return R;
}