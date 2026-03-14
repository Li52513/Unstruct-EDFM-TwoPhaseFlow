/**
 * @file Well_WellControlTypes.h
 * @brief 井模型控制类型与调度步骤定义 (Production-Ready)
 */
#ifndef WELL_WELLCONTROLTYPES_H
#define WELL_WELLCONTROLTYPES_H

#include <string>

enum class WellControlMode {
    BHP,    ///< 定井底流压 (Bottom Hole Pressure)
    Rate    ///< 定流量 (Mass Rate or Volume Rate)
};

enum class WellTargetDomain {
    Matrix,     ///< 基岩网格
    Fracture    ///< 裂缝网格
};

enum class WellComponentMode {
    Water,  ///< 仅水相
    Gas,    ///< 仅气相 (CO2)
    Total   ///< 总相
};

enum class WellAxis {
    X,
    Y,
    Z,
    None    ///< 适用于 2D 或纯用户指定覆盖 (Override) 的情况
};

struct WellScheduleStep {
    double t_start = 0.0;
    double t_end = 0.0;
    std::string well_name = "";
    WellTargetDomain domain = WellTargetDomain::Matrix;
    WellControlMode control_mode = WellControlMode::BHP;
    double target_value = 0.0;
    WellComponentMode component_mode = WellComponentMode::Total;
    double rw = 0.1;
    double skin = 0.0;
    WellAxis well_axis = WellAxis::None;

    // Matrix: cell local index | Fracture: solverIndex or fracture-local index
    int completion_id = -1;

    double wi_override = -1.0;
    double L_override = -1.0;
    // For WellComponentMode::Total:
    // - If frac_w + frac_g > 0, assembler normalizes and uses them.
    // - If both are 0, assembler falls back to mobility-weighted auto split.
    double frac_w = 0.0;
    double frac_g = 0.0;

    // Optional injection stream temperature (K).
    // Effective for injection branch in energy equation when > 0.
    double injection_temperature = -1.0;
    bool injection_is_co2 = false;
};

#endif // WELL_WELLCONTROLTYPES_H
