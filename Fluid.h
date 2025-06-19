#pragma once

class Fluid {
public:
    double fluid_rho;     // 流体的密度 (kg/m^3)
    double fluid_viscosity; // 流体的粘度 (Pa.s)

    // 构造函数，初始化流体物性参数
    Fluid(double rho = 1000.0, double viscosity = 0.001); // 默认值为水的常规物性参数

    // 设置流体物性参数
    void setProperties(double rho, double viscosity);
};