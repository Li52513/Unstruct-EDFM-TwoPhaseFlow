#include "2D_FracIndex.h"
#include "2D_FractureNetwork.h" // 需要访问 fractures 容器

FracElemIndex_2D buildFracElemIndex_2D(const FractureNetwork_2D& frNet)
{
    FracElemIndex_2D idx;
    const auto& fracs = frNet.getFractures(); // 需在 Network 中提供此访问器

    idx.offset.resize(fracs.size() + 1);
    idx.offset[0] = 0;

    for (size_t i = 0; i < fracs.size(); ++i) {
        // 累加每个裂缝包含的微观单元数量 (fracCells)
        idx.offset[i + 1] = idx.offset[i] + fracs[i].fracCells.size();
    }

    idx.total = idx.offset.back();
    return idx;
}