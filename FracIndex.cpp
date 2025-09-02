#include "FracIndex.h"
#include "FractureNetwork.h"  // 这里再包含完整定义
#include <cstddef>

FracElemIndex buildFracElemIndex(const FractureNetwork& frNet)
{
    FracElemIndex idx;
    idx.offset.resize(frNet.fractures.size() + 1, 0);
    std::size_t acc = 0;
    for (std::size_t f = 0; f < frNet.fractures.size(); ++f) {
        idx.offset[f] = acc;
        acc += frNet.fractures[f].elements.size();
    }
    idx.offset.back() = acc;
    idx.total = acc;
    return idx;
}