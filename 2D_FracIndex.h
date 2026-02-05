#pragma once

#include <vector>
#include <cstddef> // for size_t

// 前向声明
class FractureNetwork_2D;

/**
 * @brief 2D 裂缝网络的单元全局索引辅助结构
 * @details 用于将 (FractureID, LocalCellID) 映射为全局唯一的 Matrix Index
 */

struct FracElemIndex_2D
{
    std::vector<std::size_t> offset;  ///< 偏移量数组，size = fractures.size() + 1
    std::size_t total = 0;            ///< 总微观单元数量
};

/**
 * @brief 构建 2D 裂缝网络的全局索引
 * @param frNet 2D 裂缝网络对象
 * @return 填充好的索引结构
 */
FracElemIndex_2D buildFracElemIndex_2D(const FractureNetwork_2D& frNet);