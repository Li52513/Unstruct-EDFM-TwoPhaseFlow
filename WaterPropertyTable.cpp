#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include "WaterPropertyTable.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <experimental/filesystem>
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
    if (P < pressures_.front() || P > pressures_.back()
        || T < temps_.front() || T > temps_.back())
    {
        throw std::out_of_range("WaterPropertyTable: (P,T) outside bounds");
    }

    // 找压力区间 iP0,iP1 与分数 xFrac
    auto itP1 = std::upper_bound(pressures_.begin(), pressures_.end(), P);
    size_t iP1 = itP1 - pressures_.begin(), iP0 = iP1 - 1;
    double P0 = pressures_[iP0], P1 = pressures_[iP1];
    double xFrac = (P - P0) / (P1 - P0);

    // 找温度区间 iT0,iT1 与分数 yFrac
    auto itT1 = std::upper_bound(temps_.begin(), temps_.end(), T);
    size_t iT1 = itT1 - temps_.begin(), iT0 = iT1 - 1;
    double T0 = temps_[iT0], T1 = temps_[iT1];
    double yFrac = (T - T0) / (T1 - T0);

    WaterProperties R;
    // 对每个字段分别做 bicubicInterpolate
    R.rho = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
        [&](const WaterProperties& w) { return w.rho; });
    R.mu = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
        [&](const WaterProperties& w) { return w.mu;  });
    R.cp = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
        [&](const WaterProperties& w) { return w.cp;  });
    R.cv = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
        [&](const WaterProperties& w) { return w.cv;  });
    R.h = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
        [&](const WaterProperties& w) { return w.h;   });
    R.k = bicubicInterpolate(data_, pressures_, temps_, iP0, iT0, xFrac, yFrac,
        [&](const WaterProperties& w) { return w.k;   });

    return R;
}

