#pragma once
#include <vector>
#include <map>
#include "Cell.h"
#include "Fracture.h"

/// 初始化所有基岩单元和裂缝段的压力
/// @param cells_in     输入的基岩单元数组
/// @param cellId2idx   Cell.id -> cells_in 的索引映射
/// @param fracture_network    所有裂缝
inline void initializePressure
(
    vector<Cell>& cells,
    const map<int, int>& cellId2idx,
    vector<Fracture>& fractures,
    double                 p_init = 1e5   // 默认初始压强：1e5 Pa
)
{
    // 1) 给所有基岩单元一个统一初始压力
    for (auto& C : cells)
    {
        C.pressure = p_init;
    }

    // 2) 给每条裂缝的每个段一个初始压力
    //    这里我们简单取它所在宿主单元的压力作为段的初始 p_fr
    for (auto& F : fractures) {
        for (auto& E : F.elements) {
            int hostID = E.cellID;                       // 段所属的 Cell.id
            int idx = cellId2idx.at(hostID);          // 找到在 cells[] 中的下标
            E.p_fr = cells[idx].pressure;            // 初始裂缝段压
        }
    }
};

