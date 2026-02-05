#pragma once
#include <functional>
struct MatrixEdge
{
    int id;         // 临时编号
    int n1, n2;     // 两个端点 ID (确保 n1 < n2 以便去重)

    // 构造函数：自动排序节点 ID
    MatrixEdge(int node1, int node2, int edgeID = -1)
    {
        if (node1 < node2) { n1 = node1; n2 = node2; }
        else { n1 = node2; n2 = node1; }
        id = edgeID;
    }

    // 重载 == 和 Hash，用于 unordered_set 去重
    bool operator==(const MatrixEdge& other) const {
        return n1 == other.n1 && n2 == other.n2;
    }
};

// 哈希函数，用于将 MatrixEdge 放入 unordered_map/set
struct MatrixEdgeHash {
    size_t operator()(const MatrixEdge& e) const {
        // 简单的 Cantor pairing 或位移组合
        return std::hash<int>()(e.n1) ^ (std::hash<int>()(e.n2) << 1);
    }
};