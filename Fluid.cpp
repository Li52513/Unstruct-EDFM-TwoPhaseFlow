#include "Fluid.h"

// 构造函数，初始化流体物性参数
Fluid::Fluid(double rho, double viscosity)
    : fluid_rho(rho), fluid_viscosity(viscosity) {
}

// 设置流体物性参数
void Fluid::setProperties(double rho, double viscosity)
{
    fluid_rho = rho;
    fluid_viscosity = viscosity;
}
