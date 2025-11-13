#pragma once
#include "PropertiesSummary.h" // 里面定义 SolidProperties_Frac 等

enum class FractureElementType { Blocking, Conductive };

namespace fracture 
{
    
    // 具体常数
    static constexpr double BASE_RHO_F = 2500.0;
    static constexpr double BASE_CP_F = 1000.0;
    static constexpr double BASE_K_F = 2.0;
    static constexpr double BASE_COMP_F = 1e-10;

    /// 根据 elem.type 和 (P,T) 返回裂缝固相物性
    inline SolidProperties_Frac computeSolidProperties(FractureElementType type, double P, double T)
    {
        SolidProperties_Frac f;
        f.rho_f = BASE_RHO_F;
        f.cp_f = BASE_CP_F;
        f.k_f = BASE_K_F;
        f.compressibility = BASE_COMP_F;

        switch (type)
        {
        case FractureElementType::Blocking:
            f.phi_f = 0.05;  f.permeability = 1e-15; f.aperture = 1e-5;  break;
        case FractureElementType::Conductive:
            f.phi_f = 0.20;  f.permeability = 1e-10; f.aperture = 1e-3;  break;
        }
        return f;
    }

}
