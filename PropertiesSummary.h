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
    double dRho_dP = 0.0;
};

struct CO2Properties
{
    double rho = 0.0;
    double mu = 0.0;
    double cp = 0.0;
    double cv = 0.0;
    double k = 0.0;
    double h = 0.0;
    double dRho_dP = 0.0;
};

struct SolidProperties_RockMatrix 
{
    double phi_r=0.0;
    double kxx = 0.0;
    double kyy = 0.0;
    double kzz = 0.0;
    double compressibility = 0.0;
    double rho_r = 0.0;
    double cp_r = 0.0;
    double k_r = 0.0;
};

struct SolidProperties_Frac 
{
    double phi_f= 0.0;
    double permeability = 0.0;
    double compressibility = 0.0;
    double aperture = 0.0;
    double rho_f = 0.0;
    double cp_f = 0.0;
    double k_f = 0.0;
};
