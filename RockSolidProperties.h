#pragma once
#include "Cell.h"
#include"PropertiesSummary.h"

namespace rock
{
    /// 根据 cell.region 和 (P,T) 返回基岩固相物性
    SolidProperties_RockMatrix computeSolidProperties(Cell::RegionType region, double P, double T);

    // 具体常数或经验公式
    static constexpr double BASE_RHO = 2650.0;  // 密度
    static constexpr double BASE_CP = 1000.0;   // 比热容
    static constexpr double BASE_K = 2.5;       // 导热系数
    static constexpr double BASE_COMP = 1e-9;   // 可压缩系数

    // 函数定义
    inline SolidProperties_RockMatrix computeSolidProperties(Cell::RegionType region, double P, double T)
    {
        SolidProperties_RockMatrix s;
        s.rho_r = BASE_RHO;
        s.cp_r = BASE_CP;
        s.k_r = BASE_K;
        s.compressibility = BASE_COMP;

        // 现在均是常数，后期可以优化放入不同的方程
        switch (region)
        {
        case Cell::RegionType::Low:
            s.phi_r = 0.3;  s.kxx = 1e-15; s.kyy = 1e-15; s.kzz = 1e-15;  break;
        case Cell::RegionType::Medium:
            s.phi_r = 0.3;  s.kxx = 1e-13; s.kyy = 1e-13; s.kzz = 1e-13; break;
        case Cell::RegionType::High:
            s.phi_r = 0.3;  s.kxx = 1e-11; s.kyy = 1e-11; s.kzz = 1e-11; break;
        }
        return s;
    }
}