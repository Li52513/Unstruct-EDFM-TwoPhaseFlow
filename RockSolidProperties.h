#pragma once
#include "Cell.h"
#include "PhysicalPropertiesManager.h"  // SolidProperties_RockMatrix

namespace rock
{

    /// 根据 cell.region 和 (P,T) 返回基岩固相物性
    SolidProperties_RockMatrix computeSolidProperties ( Cell::RegionType region,double P, double T);


}