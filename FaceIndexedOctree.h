#pragma once
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>


#include "Face.h" 
#include "AABB.h"

// =========================================================
// 八叉树参数配置
// =========================================================
const int OCTREE_MAX_ITEMS_PER_LEAF = 40; // 阈值：叶子节点超过多少个面时分裂
const int OCTREE_MAX_DEPTH = 12;          // 阈值：最大深度，防止死循环或过细

// =========================================================
// 八叉树节点类
// =========================================================
class OctreeNode
{
public:
    AABB box;                       // 节点包围盒
    std::vector<int> faceIndices;   // 存储的面 ID 列表 (仅叶子节点有效)
    OctreeNode* children[8];        // 子节点指针
    bool isLeaf;                    // 是否为叶子

    OctreeNode(const AABB& bbox) : box(bbox), isLeaf(true) {
        for (int i = 0; i < 8; ++i) children[i] = nullptr;
    }

    ~OctreeNode() {
        for (int i = 0; i < 8; ++i) {
            if (children[i]) delete children[i];
        }
    }
};

// =========================================================
// 面索引八叉树类 (Face Indexed Octree)
// =========================================================
class FaceIndexedOctree
{
public:
    FaceIndexedOctree() : root_(nullptr) {}

    // 带参数构造函数：直接利用给定的全局包围盒和面列表构建树
    FaceIndexedOctree(const AABB& globalBox, const std::vector<Face>& faces) : root_(nullptr)
    {
        if (faces.empty()) return;

        // 1. 初始化根节点
        root_ = new OctreeNode(globalBox);

        // 2. 插入所有面
        // 直接调用内部递归插入函数
        for (size_t i = 0; i < faces.size(); ++i) {
            insertRecursive(root_, static_cast<int>(i), faces[i].boundingBox, faces, 0);
        }

        std::cout << "[Octree] Constructed with " << faces.size() << " faces." << std::endl;
    }

    ~FaceIndexedOctree() { if (root_) delete root_; }

    // 禁止拷贝 (防止深拷贝带来的性能问题)
    FaceIndexedOctree(const FaceIndexedOctree&) = delete;
    FaceIndexedOctree& operator=(const FaceIndexedOctree&) = delete;

    // -----------------------------------------------------
    // 构建树
    // faces: 传入网格的所有面，用于计算包围盒及分裂时重分配
    // -----------------------------------------------------
    void build(const std::vector<Face>& faces)
    {
        if (faces.empty()) return;
        if (root_) { delete root_; root_ = nullptr; }
        // 1. 计算全局包围盒 (Global Bounding Box)
        Vector globalMin(1e30, 1e30, 1e30);
        Vector globalMax(-1e30, -1e30, -1e30);

        for (const auto& f : faces) {
            updateMinMax(globalMin, globalMax, f.boundingBox.min);
            updateMinMax(globalMin, globalMax, f.boundingBox.max);
        }

        // 添加微小余量，避免边界点判定误差
        Vector margin(1e-4, 1e-4, 1e-4);
        AABB globalBox(globalMin - margin, globalMax + margin);

        root_ = new OctreeNode(globalBox);

        // 2. 插入所有面
        // 注意：这里需要传入 faces 引用，因为节点分裂时需要回溯旧数据的包围盒
        for (size_t i = 0; i < faces.size(); ++i) {
            insertRecursive(root_, static_cast<int>(i), faces[i].boundingBox, faces, 0);
        }

        std::cout << "[Octree] Built successfully. Root Box: "
            << "(" << globalMin.m_x << "," << globalMin.m_y << "," << globalMin.m_z << ") -> "
            << "(" << globalMax.m_x << "," << globalMax.m_y << "," << globalMax.m_z << ")" << std::endl;
    }

    // -----------------------------------------------------
    // 空间查询
    // queryBox: 搜索区域 (如裂缝面的包围盒)
    // results:  输出相交的候选 Face ID 列表
    // -----------------------------------------------------
    void query(const AABB& queryBox, std::vector<int>& results) const
    {
        results.clear();
        if (!root_) return;
        queryRecursive(root_, queryBox, results);

        // 结果去重 (因为一个面可能跨越多个八叉树节点，会被多次添加)
        if (!results.empty()) {
            std::sort(results.begin(), results.end());
            auto last = std::unique(results.begin(), results.end());
            results.erase(last, results.end());
        }
    }

private:
    OctreeNode* root_;

    // 辅助：更新极值
    void updateMinMax(Vector& minP, Vector& maxP, const Vector& p) {
        minP.m_x = std::min(minP.m_x, p.m_x); minP.m_y = std::min(minP.m_y, p.m_y); minP.m_z = std::min(minP.m_z, p.m_z);
        maxP.m_x = std::max(maxP.m_x, p.m_x); maxP.m_y = std::max(maxP.m_y, p.m_y); maxP.m_z = std::max(maxP.m_z, p.m_z);
    }

    // 辅助：AABB 相交检测
    bool checkOverlap(const AABB& a, const AABB& b) const {
        if (a.max.m_x < b.min.m_x || a.min.m_x > b.max.m_x) return false;
        if (a.max.m_y < b.min.m_y || a.min.m_y > b.max.m_y) return false;
        if (a.max.m_z < b.min.m_z || a.min.m_z > b.max.m_z) return false;
        return true;
    }

    // -----------------------------------------------------
    // 递归插入逻辑 (核心)
    // -----------------------------------------------------
    void insertRecursive(OctreeNode* node, int faceIdx, const AABB& faceBox, const std::vector<Face>& allFaces, int depth)
    {
        // 1. 如果面与当前节点不相交，直接返回
        if (!checkOverlap(node->box, faceBox)) return;

        // 2. 如果是叶子节点
        if (node->isLeaf) {
            // 情况 A: 容量未满 或 深度已达极限 -> 直接存储
            if (node->faceIndices.size() < OCTREE_MAX_ITEMS_PER_LEAF || depth >= OCTREE_MAX_DEPTH) {
                node->faceIndices.push_back(faceIdx);
            }
            // 情况 B: 容量已满 -> 分裂 (Split)
            else {
                node->isLeaf = false; // 标记为内部节点

                // 初始化 8 个子节点
                Vector c = (node->box.min + node->box.max) * 0.5; // 中心点
                for (int i = 0; i < 8; ++i) {
                    Vector subMin, subMax;
                    // 位运算确定子象限边界: 0=min, 1=max (相对于center)
                    subMin.m_x = (i & 1) ? c.m_x : node->box.min.m_x;
                    subMax.m_x = (i & 1) ? node->box.max.m_x : c.m_x;
                    subMin.m_y = (i & 2) ? c.m_y : node->box.min.m_y;
                    subMax.m_y = (i & 2) ? node->box.max.m_y : c.m_y;
                    subMin.m_z = (i & 4) ? c.m_z : node->box.min.m_z;
                    subMax.m_z = (i & 4) ? node->box.max.m_z : c.m_z;

                    node->children[i] = new OctreeNode(AABB(subMin, subMax));
                }

                // [关键步骤] 数据重分发 (Re-distribute)
                // 将当前节点积压的旧面 ID，重新插入到子节点中
                for (int oldIdx : node->faceIndices) {
                    // 需要从 allFaces 获取旧面的 AABB
                    const AABB& oldBox = allFaces[oldIdx].boundingBox;
                    for (int k = 0; k < 8; ++k) {
                        insertRecursive(node->children[k], oldIdx, oldBox, allFaces, depth + 1);
                    }
                }
                node->faceIndices.clear(); // 清空当前节点数据 (已下沉)

                // 最后，插入本次调用的新面
                for (int k = 0; k < 8; ++k) {
                    insertRecursive(node->children[k], faceIdx, faceBox, allFaces, depth + 1);
                }
            }
        }
        // 3. 如果是内部节点 -> 递归下发
        else {
            for (int i = 0; i < 8; ++i) {
                insertRecursive(node->children[i], faceIdx, faceBox, allFaces, depth + 1);
            }
        }
    }

    // -----------------------------------------------------
    // 递归查询逻辑
    // -----------------------------------------------------
    void queryRecursive(OctreeNode* node, const AABB& queryBox, std::vector<int>& results) const
    {
        // 剪枝：如果不相交，这整棵子树都不需要看
        if (!checkOverlap(node->box, queryBox)) return;

        // 如果是叶子节点，收集数据
        if (node->isLeaf) {
            results.insert(results.end(), node->faceIndices.begin(), node->faceIndices.end());
        }
        // 否则递归子节点
        else {
            for (int i = 0; i < 8; ++i) {
                queryRecursive(node->children[i], queryBox, results);
            }
        }
    }

    
};