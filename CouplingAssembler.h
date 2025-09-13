#pragma once

#include <cstddef>

// 前置声明，避免重编译风暴
class MeshManager;
struct FieldRegistry;
struct VGParams;
struct RelPermParams;

// 创建/确保裂缝耦合场（CIw/CIg）
void ensureFracCouplingFields(FieldRegistry& reg_fr, std::size_t nSeg);

// 计算/更新：基岩-裂缝两相耦合系数 CIw / CIg
// upwind: 是否启用迎风（基于相位势）；include_gravity: 是否计入重力项；g: 重力加速度
void updateMatrixFractureCI(MeshManager& mgr,
    const FieldRegistry& Rm, FieldRegistry& Rf,
    const VGParams& vg,
    bool upwind,
    bool include_gravity, double g);

// 计算/更新：裂缝-裂缝两相耦合系数（交点 Star–Delta）TIw / TIg
void updateFractureFractureTI(MeshManager& mgr,
     FieldRegistry& Rf,
    const VGParams& vg, const RelPermParams& rp,
    bool include_gravity, double g);

// —— 只读诊断/导出（不会重新计算） —— 
void printCI_Diagnostics(
    MeshManager& mgr,
    const FieldRegistry& Rm,
    const FieldRegistry& Rf,
    std::size_t maxFracs = (std::size_t)-1,
    std::size_t maxSegsPerFrac = (std::size_t)-1
);

//void exportCI_CSV(
//    MeshManager& mgr,
//    const FieldRegistry& Rf,
//    const std::string& outdir,
//    int step,
//    double time
//);

void printTI_Diagnostics(
    MeshManager& mgr,
    const FieldRegistry& Rf,
    std::size_t maxFracs = (std::size_t)-1
);