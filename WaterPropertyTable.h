#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include "PropertiesSummary.h"
#include <functional>
using namespace std;

/**
 * @class WaterPropertyTable
 * @brief 负责：①读取水物性表格文件；②在给定 P,T 下插值返回水物性
 */

class WaterPropertyTable 
{
public:
    /// 构造时加载文件
        // 单例访问
    static WaterPropertyTable& instance();
    /// 在给定压力 P (Pa) 和 温度 T (K) 下插值返回全部水物性
    WaterProperties getProperties(const double& P, const double& T) const;

private:
    WaterPropertyTable(const std::string& filename);
    void load(const std::string& filename);

    vector<double> pressures_;    ///< 升序压力列表
    vector<double> temps_;        ///< 对第一个压力值升序温度列表

    /// data_[iP][iT] 存放 pressures_[iP], temps_[iT] 对应的 WaterProperties
    vector<vector<WaterProperties>> data_;

    // 单纯量 clamp (兼容任意 C++ 标准)
    static double clamp(double v, double lo, double hi)
    {
        return v < lo ? lo : (v > hi ? hi : v);
    }

    // Cubic Hermite 插值——声明为私有静态成员
    static double cubicHermite(double y0, double y1, double y2, double y3, double t);

    // Bicubic Catmull–Rom 插值，针对任意字段的通用接口
    static double bicubicInterpolate
    (
        const std::vector<std::vector<WaterProperties>>& data,
        const std::vector<double>& X,
        const std::vector<double>& Y,
        size_t i, size_t j,
        double xFrac, double yFrac,
        const std::function<double(const WaterProperties&)>& field
    );
};




