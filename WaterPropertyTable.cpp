#include "WaterPropertyTable.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>


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
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open water property file: " + filename);

    std::string line;
    // 跳过到 header 行（包含 "T"、"P"、"rho"）
    static const char* keys[] = { "P","T","rho","cp","cv","k","mu","h" };
    while (std::getline(in, line))
    {
        int hit = 0;
        for (auto k : keys) if (line.find(k) != std::string::npos) ++hit;
        if (hit >= 5) break;   // 5 个以上关键字命中 → 基本可判定为表头
    }

    struct Row { double P{ 0.0 }, T{ 0.0 }; WaterProperties w{}; };
    std::vector<Row> rows;

    // 逐行读取：P  T  rho  cp  cv  k  mu  h
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        Row r;
        if (iss >> r.P >> r.T
            >> r.w.rho >> r.w.cp >> r.w.cv
            >> r.w.k >> r.w.mu >> r.w.h) {
            rows.push_back(r);
        }
    }

    // 构造升序唯一 pressures_
    {
        std::vector<double> tmp;
        for (auto& r : rows) tmp.push_back(r.P);
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
        pressures_ = std::move(tmp);
    }
    // 构造 temps_（对第一个压力值 P0）
    {
        double P0 = pressures_.front();
        std::vector<double> tmp;
        for (auto& r : rows) if (r.P == P0) tmp.push_back(r.T);
        std::sort(tmp.begin(), tmp.end());
        tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
        temps_ = std::move(tmp);
    }

    // 分配 data_
    data_.assign(pressures_.size(),
        std::vector<WaterProperties>(temps_.size()));

    // 填充 data_
    for (auto& r : rows) {
        size_t iP = std::lower_bound(pressures_.begin(), pressures_.end(), r.P)
            - pressures_.begin();
        size_t iT = std::lower_bound(temps_.begin(), temps_.end(), r.T)
            - temps_.begin();
        data_[iP][iT] = r.w;
    }
}

WaterProperties WaterPropertyTable::getProperties(const double& P, const double& T) const
{
    if (P < pressures_.front() || P > pressures_.back()
        || T < temps_.front() || T > temps_.back()) 
    {
        throw std::out_of_range("WaterPropertyTable: (P,T) outside bounds");
    }

    // 找 P0, P1 的索引
    auto itP1 = std::upper_bound(pressures_.begin(), pressures_.end(), P);
    size_t iP1 = itP1 - pressures_.begin(), iP0 = iP1 - 1;
    double P0 = pressures_[iP0], P1 = pressures_[iP1];

    // 找 T0, T1 的索引
    auto itT1 = std::upper_bound(temps_.begin(), temps_.end(), T);
    size_t iT1 = itT1 - temps_.begin(), iT0 = iT1 - 1;
    double T0 = temps_[iT0], T1 = temps_[iT1];

    // 四角物性
    const auto& W00 = data_[iP0][iT0]; // at (P0, T0)
    const auto& W01 = data_[iP0][iT1]; // at (P0, T1)
    const auto& W10 = data_[iP1][iT0]; // at (P1, T0)
    const auto& W11 = data_[iP1][iT1]; // at (P1, T1)

    // 先沿 P 方向做一次插值（得到关于 T 的两个中间状态）
    double pFrac = (P - P0) / (P1 - P0);
    WaterProperties R0, R1;
    R0.rho = W00.rho + (W10.rho - W00.rho) * pFrac;
    R0.cp = W00.cp + (W10.cp - W00.cp) * pFrac;
    R0.cv = W00.cv + (W10.cv - W00.cv) * pFrac;
    R0.k = W00.k + (W10.k - W00.k) * pFrac;
    R0.mu = W00.mu + (W10.mu - W00.mu) * pFrac;
    R0.h = W00.h + (W10.h - W00.h) * pFrac;

    R1.rho = W01.rho + (W11.rho - W01.rho) * pFrac;
    R1.cp = W01.cp + (W11.cp - W01.cp) * pFrac;
    R1.cv = W01.cv + (W11.cv - W01.cv) * pFrac;
    R1.k = W01.k + (W11.k - W01.k) * pFrac;
    R1.mu = W01.mu + (W11.mu - W01.mu) * pFrac;
    R1.h = W01.h + (W11.h - W01.h) * pFrac;

    // 再沿 T 方向插值
    double tFrac = (T - T0) / (T1 - T0);
    WaterProperties result;
    result.rho = R0.rho + (R1.rho - R0.rho) * tFrac;
    result.cp = R0.cp + (R1.cp - R0.cp) * tFrac;
    result.cv = R0.cv + (R1.cv - R0.cv) * tFrac;
    result.k = R0.k + (R1.k - R0.k) * tFrac;
    result.mu = R0.mu + (R1.mu - R0.mu) * tFrac;
    result.h = R0.h + (R1.h - R0.h) * tFrac;

    return result;
}

