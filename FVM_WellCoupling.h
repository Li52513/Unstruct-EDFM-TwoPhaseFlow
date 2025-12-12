#pragma once
#include <string>
#include <vector>
#include "FVM_WellDOF.h"
#include "Solver_AssemblerCOO.h"
#include "MeshManager.h"
#include "FieldRegistry.h"

// 关键：使用你已有的按 cellId 访问函数
#include "FieldAcessForDiscre.h"  // cellScalar(reg, mesh, name, cellId, fallback)

// 由“单元索引 cidx（0..Nc-1）”取得该单元的全局 id（与 getCellId2Index 的 key 一致）
inline int cell_id_by_index(const Mesh& mesh, int cidx) {
    const auto& cells = mesh.getCells();
    assert(cidx >= 0 && cidx < (int)cells.size());
    return cells[cidx].id;  // Cell 内部保存的全局 id
}

// 由“单元索引 cidx”取得单元测度：当前项目为 2D => 使用 getCellArea(cellID)
// 若未来扩展 3D，请在 Mesh 中提供 getCellVolume 并在此改为调用体积
inline double cell_measure(const Mesh& mesh, int cidx) {
    const int cid = cell_id_by_index(mesh, cidx);
    return mesh.getCellArea(cid);  // 2D: 面积即测度
}

/**
 * Peaceman 列耦合（单元行）：
 * 对掩码命中的 cell i：
 *   A[i,i]       +=  PI_i * w_i
 *   A[i, well_lid]+= -PI_i * w_i
 * 其中 w_i = (scaleByCellMeasure ? |cell| : 1.0)
 * 注意：调用前应已 extend_linear_system_size(sys, Ntot)
 */
inline void add_peaceman_coupling_cell_rows(
    SparseSystemCOO& sys,
    const Mesh& mesh, const FieldRegistry& reg,
    const std::string& PI_name,
    const std::string& mask_name,
    const std::vector<int>& lid_cell,  // cell index -> unknown id
    int well_lid,
    bool scaleByCellMeasure = false
) {
    const int Nc = (int)lid_cell.size();
    for (int cidx = 0; cidx < Nc; ++cidx) {
        const int cid = cell_id_by_index(mesh, cidx);                      // cellId
        const double m = cellScalar(reg, mesh, mask_name.c_str(), cid, 0); // 掩码
        if (m <= 0.5) continue;

        const double PIi = cellScalar(reg, mesh, PI_name.c_str(), cid, 0.0);
        if (PIi == 0.0) continue;

        const double wi = scaleByCellMeasure ? cell_measure(mesh, cidx) : 1.0;
        const int    li = lid_cell[cidx];

        sys.addA(li, li, +PIi * wi);
        sys.addA(li, well_lid, -PIi * wi);
    }
}

/**
 * 追加井行：
 * - Pressure:  A[w,w]=1;                b[w]+= target_p
 * - Rate:      A[w,w]+=Σ(PI_i w_i);     A[w,i]+=-PI_i w_i;  b[w]+= target_Q
 * 注意：调用前应已 extend_linear_system_size(sys, Ntot)
 */
inline void add_well_row(
    SparseSystemCOO& sys,
    const Mesh& mesh, const FieldRegistry& reg,
    const std::string& PI_name, const std::string& mask_name,
    const std::vector<int>& lid_cell, int well_lid,
    WellDOF::Mode mode, double target,
    bool scaleByCellMeasure = false
) {
    if (mode == WellDOF::Mode::Pressure) {
        sys.addA(well_lid, well_lid, 1.0);
        sys.addb(well_lid, target);
        return;
    }

    // Rate 模式
    const int Nc = (int)lid_cell.size();
    double diag = 0.0;

    for (int cidx = 0; cidx < Nc; ++cidx) {
        const int cid = cell_id_by_index(mesh, cidx);
        const double m = cellScalar(reg, mesh, mask_name.c_str(), cid, 0);
        if (m <= 0.5) continue;

        const double PIi = cellScalar(reg, mesh, PI_name.c_str(), cid, 0.0);
        if (PIi == 0.0) continue;

        const double wi = scaleByCellMeasure ? cell_measure(mesh, cidx) : 1.0;
        const int    li = lid_cell[cidx];

        diag += PIi * wi;
        sys.addA(well_lid, li, -PIi * wi);
    }

    if (diag != 0.0) sys.addA(well_lid, well_lid, diag);
    sys.addb(well_lid, target);
}
