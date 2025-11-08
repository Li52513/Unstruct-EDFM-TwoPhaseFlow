#pragma once
#include <string>
#include <vector>
#include "FVM_WellDOF.h"     // 约定你已按我们前面建议扩展了 Role/mask_field/PI_field/Tin 等
#include "FVM_Peaceman.h"    // WellSpec / PeacemanParams

// 统一的井配置：几何 + Peaceman + 运行模式/目标 + 字段名 + lid
struct WellConfig {
    std::string name;                 // 井名（用于拼字段名）
    WellDOF::Role role;               // Injector / Producer
    WellSpec       geom;              // 几何（pos/rw/skin/H/perfRadius/maxHitCells）
    PeacemanParams pm;                // Peaceman 参数（每井可独立，或全局同一份拷贝）

    // 运行控制（与 DoF 对齐）
    WellDOF::Mode mode;               // Pressure / Rate
    double target;                    // BHP[Pa] 或 质量流量[kg/s]
    double Tin;                       // 注入温度[K]（仅注井使用；采井可设为0）

    // 字段名（留空则根据 name 自动生成）
    std::string mask_name;            // 建议：mask_<name>
    std::string PI_name;              // 建议：PI_<name>
    std::string pw_name;              // 建议：p_w_<name>（单值场）

    int lid;                          // 注册后的未知量索引（求解后写回）

    WellConfig() : role(WellDOF::Role::Injector),
        mode(WellDOF::Mode::Pressure),
        target(0.0), Tin(0.0), lid(-1) {
    }

    // 自动补齐字段名
    void derive_names_if_empty() {
        if (mask_name.empty()) mask_name = std::string("mask_") + name;
        if (PI_name.empty())   PI_name = std::string("PI_") + name;
        if (pw_name.empty())   pw_name = std::string("p_w_") + name;
    }

    // 转成 WellDOF（供注册与耦合函数使用）
    WellDOF toWellDOF() const {
        WellDOF w;
        w.name = pw_name;       // 这里用“井压场名”承载（与后续写回一致）
        w.role = role;
        w.mode = mode;
        w.target = target;
        w.Tin = Tin;
        w.mask_field = mask_name;
        w.PI_field = PI_name;
        w.lid = -1;
        return w;
    }
};

// ---------- 批处理辅助：一次性处理一组 WellConfig ----------

// 1) 为所有井构建 mask + PI（Peaceman）
//    你已有 build_well_mask_and_PI(...)，这里做一层循环封装
inline void build_masks_and_PI_for_all(MeshManager& mgr, FieldRegistry& reg,
    std::vector<WellConfig>& wells)
{
    for (size_t i = 0; i < wells.size(); ++i) {
        wells[i].derive_names_if_empty();
        build_well_mask_and_PI(mgr, reg, wells[i].geom, wells[i].pm,
            wells[i].mask_name, wells[i].PI_name);
    }
}

// 2) 从配置生成 WellDOF 并注册未知量，回填 lid
inline int register_well_dofs_for_all(int Nc,
    const std::vector<WellConfig>& wellsCfg,
    std::vector<WellDOF>& wells)
{
    wells.clear();
    wells.reserve(wellsCfg.size());
    for (size_t i = 0; i < wellsCfg.size(); ++i) {
        wells.push_back(wellsCfg[i].toWellDOF());
    }
    const int Ntot = register_well_unknowns(Nc, wells);
    return Ntot;
}

// 3) 求解后，将井底压写回单值场 p_w_<name>
inline void writeback_pw_fields_for_all(FieldRegistry& reg,
    const std::vector<WellDOF>& wells,
    const std::vector<double>& x)
{
    for (size_t i = 0; i < wells.size(); ++i) {
        const std::string& pw_name = wells[i].name; // 这里我们把 name 设成了 pw_name
        auto pw = reg.get<volScalarField>(pw_name.c_str());
        if (!pw) pw = reg.create<volScalarField>(pw_name.c_str(), 1, 0.0);
        (*pw)[0] = x[wells[i].lid];
    }
}
