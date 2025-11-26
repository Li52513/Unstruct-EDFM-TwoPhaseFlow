#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"




namespace IMPES_Iteration
{
    struct SaturationTransportConfig
    {
        std::string saturation = "s_w";
        std::string saturation_old = "s_w_old";
        std::string saturation_prev = "s_w_prev";
        VGParams vg_params;
        RelPermParams rp_params;
    };
}