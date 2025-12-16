#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include "WaterPropertyTable.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <experimental/filesystem>

namespace 
{
    // 在 axis（升序）里定位 x 所在的小区间 [i0,i1] 以及分数 frac
    inline void bracket2(const std::vector<double>& axis, double x,
        size_t& i0, size_t& i1, double& frac)
    {
        // x==最大端时，upper_bound 会返回 end()，要手动夹住
        auto it1 = std::upper_bound(axis.begin(), axis.end(), x);
        size_t idx1 = size_t(it1 - axis.begin());
        if (idx1 == 0) idx1 = 1;                             // x<=axis[0]
        if (idx1 >= axis.size()) idx1 = axis.size() - 1;     // x>=axis.back()
        i1 = idx1;  i0 = idx1 - 1;

        double x0 = axis[i0], x1 = axis[i1];
        frac = (x1 > x0) ? ((x - x0) / (x1 - x0)) : 0.0;    // 防 0 除
        if (frac < 0) frac = 0; else if (frac > 1) frac = 1;
    }

    inline double bilerp_linear(double v00, double v10, double v01, double v11,
        double a, double b)
    {
        // (x: a in [0,1], y: b in [0,1])
        return (1 - a) * (1 - b) * v00 + a * (1 - b) * v10 + (1 - a) * b * v01 + a * b * v11;
    }

    inline double bilerp_logpos(double v00, double v10, double v01, double v11,
        double a, double b)
    {
        // 对正量做 log 双线性。若有 <=0 值，自动退回线性
        const double eps = 1e-300;
        if (v00 > 0 && v10 > 0 && v01 > 0 && v11 > 0) {
            double L00 = std::log(v00), L10 = std::log(v10);
            double L01 = std::log(v01), L11 = std::log(v11);
            double L = bilerp_linear(L00, L10, L01, L11, a, b);
            return std::exp(L);
        }
        return bilerp_linear(v00, v10, v01, v11, a, b);
    }
}
namespace fs = std::experimental::filesystem;

// 单例的定义
WaterPropertyTable& WaterPropertyTable::instance()
{
    static WaterPropertyTable inst("D:/dataBase_MatrialProperties/Water_LAPWS-95.txt");
    return inst;
}

// 构造时自动加载
WaterPropertyTable::WaterPropertyTable(const std::string& filename)
{
    load(filename);
}

void WaterPropertyTable::load(const std::string& filename)
{
    const std::string cacheFile = filename + ".cache";

    // 1) 如果 .cache 存在，则二进制加载
    if (fs::exists(cacheFile)) {
        std::ifstream binIn(cacheFile, std::ios::binary);
        if (binIn) {
            size_t np, nt;
            binIn.read(reinterpret_cast<char*>(&np), sizeof(np));
            pressures_.resize(np);
            binIn.read(reinterpret_cast<char*>(pressures_.data()), np * sizeof(double));

            binIn.read(reinterpret_cast<char*>(&nt), sizeof(nt));
            temps_.resize(nt);
            binIn.read(reinterpret_cast<char*>(temps_.data()), nt * sizeof(double));

            data_.assign(np, std::vector<WaterProperties>(nt));
            for (size_t i = 0; i < np; ++i) {
                for (size_t j = 0; j < nt; ++j) {
                    WaterProperties w;
                    binIn.read(reinterpret_cast<char*>(&w), sizeof(w));
                    data_[i][j] = w;
                }
            }
            return;
        }
        // 打开失败则继续文本解析
    }

    // 2) 文本解析
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open water property file: " + filename);

    std::string line;
    // 跳过表头行，直到包含 “Density” 与 “Viscosity”
    while (std::getline(in, line)) {
        if (line.find("Density") != std::string::npos &&
            line.find("Viscosity") != std::string::npos) {
            break;
        }
    }

    struct Row { double P{}, T{}; WaterProperties w{}; };
    std::vector<Row> rows;

    // 逐行读取：P  T  rho  mu  cp  cv  h  k
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        Row r{};  // P, T 自动归零
        if (iss >> r.P >> r.T
            >> r.w.rho  // 密度
            >> r.w.mu   // 粘度
            >> r.w.cp   // 定压比热
            >> r.w.cv   // 定容比热
            >> r.w.h    // 比焓
            >> r.w.k)   // 导热系数
        {
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
    data_.assign(np, std::vector<WaterProperties>(nt));
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

double WaterPropertyTable::cubicHermite(double y0, double y1, double y2, double y3, double t) {
    double m1 = 0.5 * (y2 - y0);
    double m2 = 0.5 * (y3 - y1);
    double t2 = t * t, t3 = t2 * t;
    return (2 * y1 - 2 * y2 + m2 + m1) * t3
        + (-3 * y1 + 3 * y2 - 2 * m1 - m2) * t2
        + m1 * t
        + y1;
}

double WaterPropertyTable::bicubicInterpolate(
    const std::vector<std::vector<WaterProperties>>& data,
    const std::vector<double>& X,
    const std::vector<double>& Y,
    size_t i, size_t j,
    double xFrac, double yFrac,
    const std::function<double(const WaterProperties&)>& field)
{
    double arr[4];
    // 沿 Y 对 4 行做 Hermite
    for (int di = -1; di <= 2; ++di) {
        int ii = int(clamp(double(i + di), 0.0, double(X.size() - 1)));
        double v[4];
        for (int dj = -1; dj <= 2; ++dj) {
            int jj = int(clamp(double(j + dj), 0.0, double(Y.size() - 1)));
            v[dj + 1] = field(data[ii][jj]);
        }
        arr[di + 1] = cubicHermite(v[0], v[1], v[2], v[3], yFrac);
    }
    // 再沿 X 对 4 个结果做一次 Hermite
    return cubicHermite(arr[0], arr[1], arr[2], arr[3], xFrac);
}

WaterProperties WaterPropertyTable::getProperties(const double& P, const double& T) const
{
    //if (P < pressures_.front() || P > pressures_.back()
    //    || T < temps_.front() || T > temps_.back())
    //{
    //    throw std::out_of_range("WaterPropertyTable: (P,T) outside bounds");
    //}

    //// 找压力区间 iP0,iP1 与分数 xFrac
    //auto itP1 = std::upper_bound(pressures_.begin(), pressures_.end(), P);
    //size_t iP1 = itP1 - pressures_.begin(), iP0 = iP1 - 1;
    //double P0 = pressures_[iP0], P1 = pressures_[iP1];
    //double xFrac = (P - P0) / (P1 - P0);

    //// 找温度区间 iT0,iT1 与分数 yFrac
    //auto itT1 = std::upper_bound(temps_.begin(), temps_.end(), T);
    //size_t iT1 = itT1 - temps_.begin(), iT0 = iT1 - 1;
    //double T0 = temps_[iT0], T1 = temps_[iT1];
    //double yFrac = (T - T0) / (T1 - T0);

    //WaterProperties R;
    //// 对每个字段分别做 bicubicInterpolate
    //R.rho = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
    //    [&](const WaterProperties& w) { return w.rho; });
    //R.mu = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
    //    [&](const WaterProperties& w) { return w.mu;  });
    //R.cp = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
    //    [&](const WaterProperties& w) { return w.cp;  });
    //R.cv = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
    //    [&](const WaterProperties& w) { return w.cv;  });
    //R.h = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
    //    [&](const WaterProperties& w) { return w.h;   });
    //R.k = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
    //    [&](const WaterProperties& w) { return w.k;   });

    //return R;

    if (P < pressures_.front() || P > pressures_.back()
        || T < temps_.front() || T > temps_.back())
        throw std::out_of_range("WaterPropertyTable: (P,T) outside bounds");

    size_t iP0, iP1, iT0, iT1; double a, b;
    bracket2(pressures_, P, iP0, iP1, a);   // a: P 向的分数
    bracket2(temps_, T, iT0, iT1, b);   // b: T 向的分数

    const auto& V00 = data_[iP0][iT0];
    const auto& V10 = data_[iP1][iT0];
    const auto& V01 = data_[iP0][iT1];
    const auto& V11 = data_[iP1][iT1];

    WaterProperties R;
    // 水的 μ 变化也挺陡，你可以二选一：
    // 1) 完全线性（更“保守”、与老结果接近）：
    //R.mu  = bilerp_linear (V00.mu , V10.mu , V01.mu , V11.mu , a, b);
    // 2) 用 log-lin（更贴近真实数量级变化，推荐）：
    R.mu = bilerp_logpos(V00.mu, V10.mu, V01.mu, V11.mu, a, b);

    R.rho = bilerp_linear(V00.rho, V10.rho, V01.rho, V11.rho, a, b);
    R.cp = bilerp_linear(V00.cp, V10.cp, V01.cp, V11.cp, a, b);
    R.cv = bilerp_linear(V00.cv, V10.cv, V01.cv, V11.cv, a, b);
    R.h = bilerp_linear(V00.h, V10.h, V01.h, V11.h, a, b);
    R.k = bilerp_linear(V00.k, V10.k, V01.k, V11.k, a, b);
    return R;


}

