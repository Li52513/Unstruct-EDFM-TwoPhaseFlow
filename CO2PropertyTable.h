#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include "PropertiesSummary.h"
#include <functional>

/// CO2 物性表格：先插值压力，再插值温度
class CO2PropertyTable 
{
public:
 
   static CO2PropertyTable& instance();
   CO2Properties getProperties(double P, double T) const;

private:
  
   CO2PropertyTable(const std::string& filename);
   void load(const std::string& filename);

   std::vector<double> pressures_;       ///< 升序压力列表
   std::vector<double> temps_;           ///< 对 pressures_[0] 升序温度列表

   /// data_[iP][iT] 存放 pressures_[iP], temps_[iT] 对应的 CO2Properties
   std::vector<std::vector<CO2Properties>> data_;

   /// clamp: 防止下标越界
   static double clamp(double v, double lo, double hi) {
       return v < lo ? lo : (v > hi ? hi : v);
   }

   /// 1D Catmull–Rom 三次 Hermite 插值
   static double cubicHermite(double y0, double y1, double y2, double y3, double t);

   /// 2D 双三次 Catmull–Rom 插值，用于任意 CO2Properties 字段
   static double bicubicInterpolate(
       const std::vector<std::vector<CO2Properties>>& data,
       const std::vector<double>& X,
       const std::vector<double>& Y,
       size_t i, size_t j,
       double xFrac, double yFrac,
       const std::function<double(const CO2Properties&)>& field // Fixed type declaration
   );
};
