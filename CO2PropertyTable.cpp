#include "CO2PropertyTable.h"
#include <fstream>
#include <sstream>
#include <algorithm>


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

void CO2PropertyTable::load(const std::string& filename) 
{
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open CO2 property file: " + filename);

    // 跳表头，直到含 “Density” 
    std::string line;
    while (std::getline(in, line)) {
        if (line.find("Density") != std::string::npos &&
            line.find("Viscosity") != std::string::npos) {
            break;
        }
    }

    struct Row { double P{ 0.0 }, T{ 0.0 }; CO2Properties w{}; };
    std::vector<Row> rows;

    // 逐行读取：P  T  rho  mu  cp  cv  h  k
    while (std::getline(in, line)) 
    {
        if (line.empty()) continue;
        std::istringstream iss(line);
        Row r;
        if (iss >> r.P >> r.T
            >> r.w.rho >> r.w.mu >> r.w.cp
            >> r.w.cv >> r.w.h >> r.w.k) 
        {
           
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
        std::vector<CO2Properties>(temps_.size()));

    // 填充 data_
    for (auto& r : rows) {
        size_t iP = std::lower_bound(pressures_.begin(), pressures_.end(), r.P) - pressures_.begin();
        size_t iT = std::lower_bound(temps_.begin(), temps_.end(), r.T) - temps_.begin();
        data_[iP][iT] = r.w;
    }
}

template<class F>
static inline double bilinear(double x, double y, double x0, double x1, double y0, double y1, F f00, F f10, F f01, F f11)
{
    const double tx = (x - x0) / (x1 - x0);
    const double ty = (y - y0) / (y1 - y0);
    return (1 - ty) * ((1 - tx) * f00 + tx * f10) +
        ty * ((1 - tx) * f01 + tx * f11);
}

CO2Properties CO2PropertyTable::getProperties(double P, double T) const
{
    /* ---------- 1. 边界检查 ---------- */
    if (P < pressures_.front() || P > pressures_.back()
        || T < temps_.front() || T > temps_.back())
        throw std::out_of_range("CO2PropertyTable: (P,T) outside bounds");

    /* ---------- 2. 找到包围区间 ---------- */
    auto itP1 = std::upper_bound(pressures_.begin(), pressures_.end(), P);
    size_t iP1 = itP1 - pressures_.begin(), iP0 = iP1 - 1;
    double P0 = pressures_[iP0], P1 = pressures_[iP1];

    auto itT1 = std::upper_bound(temps_.begin(), temps_.end(), T);
    size_t iT1 = itT1 - temps_.begin(), iT0 = iT1 - 1;
    double T0 = temps_[iT0], T1 = temps_[iT1];

    const CO2Properties& W00 = data_[iP0][iT0];   // (P0,T0)
    const CO2Properties& W10 = data_[iP1][iT0];   // (P1,T0)
    const CO2Properties& W01 = data_[iP0][iT1];   // (P0,T1)
    const CO2Properties& W11 = data_[iP1][iT1];   // (P1,T1)

    CO2Properties res;

    /* ---------- 3. 线性双线性插值的量 ---------- */
    res.cp = bilinear(P, T, P0, P1, T0, T1,
        W00.cp, W10.cp, W01.cp, W11.cp);
    res.cv = bilinear(P, T, P0, P1, T0, T1,
        W00.cv, W10.cv, W01.cv, W11.cv);
    res.h = bilinear(P, T, P0, P1, T0, T1,
        W00.h, W10.h, W01.h, W11.h);

    /* ---------- 4. 对数空间插值的量 (ρ, μ, k) ---------- */
    auto logInterp = [&](double v00, double v10, double v01, double v11)
        {
            return std::exp(
                bilinear(P, T, P0, P1, T0, T1,
                    std::log(v00), std::log(v10),
                    std::log(v01), std::log(v11)));
        };

    res.rho = logInterp(W00.rho, W10.rho, W01.rho, W11.rho);
    res.mu = logInterp(W00.mu, W10.mu, W01.mu, W11.mu);
    res.k = logInterp(W00.k, W10.k, W01.k, W11.k);

    return res;
}