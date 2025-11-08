#pragma once
#include <string>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "UserDefineVarType.h"

// ----------------- 井几何参数（注/采共用） -----------------
struct WellSpec 
{
    std::string name = "WELL"; // 井名（日志/导出用）
    Vector pos;                // 井在 2D/3D 中的投影坐标（本实现用于 2D 平面）
    double rw = 0.1;         // 井筒半径 [m]
    double skin = 0.0;         // skin 系数 [-]
    double H = 1.0;         // 有效厚度（2D 模型厚度）[m]
    double perfRadius = 0.0;   // 穿孔半径（>0 用半径选单元；=0 用最近 N 个）
    int    maxHitCells = 1;    // perfRadius 为 0 时，选最近的 N 个单元
};

// ----------------- 物性/口径参数（注/采共用） -----------------
struct PeacemanParams
{
    // 黏度/密度：用于从体积导纳→质量导纳
    double mu = 1.0e-4;  // [Pa·s]
    double rho = 700.0;   // [kg/m3]
    // 等效半径 r_e = reFactor * sqrt(A/pi)
    double reFactor = 0.28;
    // k_h 来源：优先 kxx/kyy，否则 k，否则 fallbackKh
    double fallbackKh = 1e-14; // [m^2]
    // 输出导纳口径：true=质量导纳(kg/(Pa·s))；false=体积导纳(m3/(Pa·s))
    bool PI_is_mass = true;
};


// ----------------- 输出字段名（与步骤4完全匹配） -----------------
struct WellFieldNames 
{
    std::string mask; // "inj_mask" 或 "prod_mask"
    std::string PI;   // "PI_inj"   或 "PI_prod"
};

// 注入井：构建 inj_mask + PI_inj
void build_injection_mask_and_PI
(
    MeshManager& mgr, FieldRegistry& reg,
    const WellSpec& spec, const PeacemanParams& prm
);

// 采出井：构建 prod_mask + PI_prod
void build_production_mask_and_PI
(
    MeshManager& mgr, FieldRegistry& reg,
    const WellSpec& spec, const PeacemanParams& prm
);



void build_well_mask_and_PI
(
    MeshManager& mgr, FieldRegistry& reg,
    const WellSpec& spec, const PeacemanParams& prm,
    const std::string& mask_name, const std::string& PI_name
);
