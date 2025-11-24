#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "BCAdapter.h"

namespace IMPES_Iteration
{
    struct PressureAssemblyConfig
    {
        std::string operator_tag = "p_w_IMPES"; // pressure operator tag for nm
        std::string pressure_field = "p_w";    // current eval pressure field
        std::string pressure_old_field = "p_w_old";
        std::string pressure_prev_field = "p_w_prev";
        std::string pressure_g = "p_g";   //current CO2 pressure
        std::string Pc_field = "Pc";
        Vector gravity = { 0.0, 0.0, 0.0 };
        bool enable_buoyancy = false;
        int gradient_smoothing = 0;
    };

    //¹¤¾ßº¯Êý
    namespace detailed
    {

    }

}