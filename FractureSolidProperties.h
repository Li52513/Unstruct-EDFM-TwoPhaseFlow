#pragma once
#include "FractureTypes.h"    // 这里拿到枚举
#include "PropertiesSummary.h" // 里面定义 SolidProperties_Frac 等

namespace fracture 
{

    /// 根据 elem.type 和 (P,T) 返回裂缝固相物性
    SolidProperties_Frac computeSolidProperties(FractureElementType type, double P, double T);

}
