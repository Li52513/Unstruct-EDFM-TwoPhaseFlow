#include "RockSolidProperties.h"

// 具体常数或经验公式
static constexpr double BASE_RHO = 2500.0;
static constexpr double BASE_CP = 1000.0;
static constexpr double BASE_K = 2.0;
static constexpr double BASE_COMP = 1e-9;

SolidProperties_RockMatrix rock::computeSolidProperties(Cell::RegionType region, double P, double T) 
{
    SolidProperties_RockMatrix s;
    s.rho_s = BASE_RHO;
    s.cp_s = BASE_CP;
    s.k_s = BASE_K;
    s.compressibility = BASE_COMP;

    //现在均是常数，后期可以优化放入不同的方程
    switch (region) 
    {
    case Cell::RegionType::Low:
        s.porosity = 0.05;  s.permeability = 1e-14;  break;
    case Cell::RegionType::Medium:
        s.porosity = 0.15;  s.permeability = 1e-12;  break;
    case Cell::RegionType::High:
        s.porosity = 0.30;  s.permeability = 1e-10;  break;
    }
    return s;
}