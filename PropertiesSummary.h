#pragma once

/**
 * @brief 物理属性：水相与 CO₂ 相共用结构
 * rho = 密度 [kg/m³]
 * mu  = 黏度 [Pa·s]
 * cp  = 定压比热 [J/(kg·K)]
 * cv  = 定容比热 [J/(kg·K)]
 * k   = 导热系数 [W/(m·K)]
 * h   = 比焓 [J/kg]
 */
struct WaterProperties 
{
    double rho = 0.0;
    double mu = 0.0;
    double cp = 0.0;
    double cv = 0.0;
    double k = 0.0;
    double h = 0.0;
};
using CO2Properties = WaterProperties;  // 字段完全一致

struct SolidProperties_RockMatrix 
{
    double porosity=0.0;
    double permeability = 0.0;
    double compressibility = 0.0;
    double rho_s = 0.0;
    double cp_s = 0.0;
    double k_s = 0.0;
};

struct SolidProperties_Frac 
{
    double porosity = 0.0;
    double permeability = 0.0;
    double compressibility = 0.0;
    double aperture = 0.0;
    double rho_s = 0.0;
    double cp_s = 0.0;
    double k_s = 0.0;
};
