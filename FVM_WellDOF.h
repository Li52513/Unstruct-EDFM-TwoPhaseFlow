#pragma once
#include <string>
#include <vector>

struct WellDOF {
    enum class Mode { Pressure, Rate };
    enum class Role { Injector, Producer };   // 新增：井角色

    std::string name;                         // 井名，用于生成 p_w_<name>
    Mode   mode = Mode::Pressure;
    Role   role = Role::Injector;             // 新增：默认当作注井
    double target = 0.0;                      // BHP[Pa] 或 质量流量[kg/s]
    double Tin = 0.0;                         // 仅注井有效：注入温度[K]（多井时各不相同）

    // 新增：每口井自己的掩码/PI 字段名（如 "mask_INJ1" / "PI_INJ1"）
    std::string mask_field;
    std::string PI_field;

    int    lid = -1;                          // 扩展未知里的位置
};

// 保持原注册函数不变
inline int register_well_unknowns(int Nc, std::vector<WellDOF>& wells)
{
    int lid = Nc;
    for (auto& w : wells) w.lid = lid++;
    return lid;
}
