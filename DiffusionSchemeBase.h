#pragma once
#include <string>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "UserDefineVarType.h"

namespace FVM {
    namespace Diffusion {

        struct FaceCoeffNames 
        {
            std::string a_vol;  // a_f（体积形，不含ρ）
            std::string s_vol;  // s_f（体积形，交叉项+浮力等记账）
        };

    }
}