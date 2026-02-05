#pragma once

/**
 * @brief 物理属性：
 * rho = 密度 [kg/m³]
 * mu  = 黏度 [Pa·s]
 * cp  = 定压比热 [J/(kg·K)]
 * cv  = 定容比热 [J/(kg·K)]
 * k   = 导热系数 [W/(m·K)]
 * c   = 可压缩系数 [1/Pa]
 * h   = 比焓 [J/kg]
 */
struct WaterProperties 
{
	double rho = 0.0; // 密度
	double mu = 0.0;  // 黏度
	double cp = 0.0;  // 定压比热
	double cv = 0.0;  // 定容比热
	double k = 0.0; // 导热系数
	double h = 0.0; // 比焓
	double c = 0.0; // compressibility, 1/Pa
	double dRho_dP = 0.0; // 密度对压力的导数, kg/(m³·Pa)
};

struct CO2Properties
{
    double rho = 0.0;
    double mu = 0.0;
    double cp = 0.0;
    double cv = 0.0;
    double k = 0.0;
    double h = 0.0;
	double c = 0.0;
    double dRho_dP = 0.0;
};

struct SolidProperties_RockMatrix 
{
	double phi_r = 0.0; // 孔隙度
	double kxx = 0.0;   // 渗透率，m²
	double kyy = 0.0;   // 渗透率，m²
	double kzz = 0.0;   // 渗透率，m²
	double compressibility = 0.0; // 可压缩系数，1/Pa
	double rho_r = 0.0; // 密度
	double cp_r = 0.0;  // 比热容
	double k_r = 0.0;   // 导热系数
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
