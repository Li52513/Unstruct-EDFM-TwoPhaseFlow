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

   // Table bounds helpers for robust AD evaluator
   double minPressure() const { return pressures_.empty() ? 0.0 : pressures_.front(); }
   double maxPressure() const { return pressures_.empty() ? 0.0 : pressures_.back(); }
   double minTemperature() const { return temps_.empty() ? 0.0 : temps_.front(); }
   double maxTemperature() const { return temps_.empty() ? 0.0 : temps_.back(); }
   double clampPressure(double P) const { return clamp(P, minPressure(), maxPressure()); }
   double clampTemperature(double T) const { return clamp(T, minTemperature(), maxTemperature()); }
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
};

