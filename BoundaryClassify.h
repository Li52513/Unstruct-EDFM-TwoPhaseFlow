
#pragma once
#include <vector>
#include <cmath>
#include "Mesh.h"

namespace BoundaryClassify 
{

    /// 按几何位置把所有边界面分成 bottom/top/left/right 四组
    struct FaceGroups 
    {
        std::vector<int> bottom;
        std::vector<int> top;
        std::vector<int> left;
        std::vector<int> right;
    };


    /// 提取所有真边界面的 ID
    inline std::vector<int> getBoundaryFaceIDs(const Mesh& mesh)
    {
        std::vector<int> ids;
        for (auto const& f : mesh.getFaces())
        {
            if (f.isBoundary())
                ids.push_back(f.id);
        }
        return ids;
    }


    /**
     * @param mesh      已构建好的 Mesh 对象
     * @param lengthX   区域在 X 方向的长度（用于判定左/右）
     * @param lengthY   区域在 Y 方向的长度（用于判定上/下）
     * @param tol       坐标比较容差，默认 1e-8
     */
    inline FaceGroups ClassifySolidMatrixMeshBoundaryFaces(
        const Mesh& mesh,
        double lengthX,
        double lengthY,
        double tol = 1e-8)
    {
        FaceGroups groups;
        auto bfaceIDs = getBoundaryFaceIDs(mesh);
        for (int fid : bfaceIDs) {
            const auto& f = mesh.getFaces()[fid - 1]; // face.id 从 1 开始
            double x = f.midpoint.m_x;
            double y = f.midpoint.m_y;

            if (std::abs(y - 0.0) < tol) groups.bottom.push_back(fid);
            else if (std::abs(y - lengthY) < tol) groups.top.push_back(fid);
            else if (std::abs(x - 0.0) < tol) groups.left.push_back(fid);
            else if (std::abs(x - lengthX) < tol) groups.right.push_back(fid);
        }
        return groups;
    }

} 
