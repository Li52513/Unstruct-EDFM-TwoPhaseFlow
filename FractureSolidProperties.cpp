#include "FractureSolidProperties.h"

static constexpr double BASE_RHO_F = 2500.0;
static constexpr double BASE_CP_F = 1000.0;
static constexpr double BASE_K_F = 2.0;
static constexpr double BASE_COMP_F = 1e-10;

SolidProperties_Frac
fracture::computeSolidProperties(FractureElementType type, double P, double T)
{
    SolidProperties_Frac f;
    f.rho_s = BASE_RHO_F;
    f.cp_s = BASE_CP_F;
    f.k_s = BASE_K_F;
    f.compressibility = BASE_COMP_F;

    switch (type) 
    {
    case FractureElementType::Blocking:
        f.porosity = 0.05;  f.permeability = 1e-15; f.aperture = 1e-5;  break;
    case FractureElementType::Conductive:
        f.porosity = 0.20;  f.permeability = 1e-10; f.aperture = 1e-3;  break;
    }
    return f;
}