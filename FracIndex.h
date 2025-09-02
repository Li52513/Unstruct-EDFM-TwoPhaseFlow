#pragma once
#include <vector>
#include <cstddef>

// 前置声明，不需要完整定义
class FractureNetwork;

struct FracElemIndex {
    std::vector<std::size_t> offset;  // size = fractures.size() + 1
    std::size_t total = 0;
};


FracElemIndex buildFracElemIndex(const FractureNetwork& frNet);