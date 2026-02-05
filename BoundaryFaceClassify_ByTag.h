#pragma once
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


#include "Mesh.h" 
#include "MeshDefinitions.h"

namespace BoundaryFaceClassify_byTag
{
    // ---------------------------- 结果类型 ----------------------------
// 统一 2D/3D 的六面分组（2D 时 z0/zL 为空）
    struct FaceGroups
    {
        std::vector<int> x0, xL;  // x=0 / x=Lx
        std::vector<int> y0, yL;  // y=0 / y=Ly
        std::vector<int> z0, zL;  // z=0 / z=Lz（Lz<=0 时为空）

        // 用于存储标准盒子以外的边界（如内部障碍物）
        // Key: Tag ID, Value: Face ID List
        std::map<int, std::vector<int>> custom;
    };

    // =========================================================
    // 基于物理组 Tag 进行 2D 边界分类
    // =========================================================
    // 功能：根据 Face 的 physicalGroupId 将其分派到 FaceGroups 的 x0/xL/y0/yL/z0/zL 中
    // 依赖：MeshDefinitions.h 中的 Tag 定义以及 Mesh::GenerateMesh3D 中的标记逻辑
    inline FaceGroups ClassifyBoundaryFacesByTag_2D(const Mesh& mesh)
    {
        FaceGroups groups;
        const auto& faces = mesh.getFaces();

        for (const auto& f : faces)
        {
            // 只有标记了物理组的才是边界 (BuildMesh2D 中非边界默认 -1)
            if (f.physicalGroupId == -1) continue;

            int pid = f.physicalGroupId;

            // 将 Tag 映射回标准的 x0, xL 结构 (为了兼容旧代码)
            switch (pid)
            {
            case MeshTags::LEFT: groups.x0.push_back(f.id); break;
            case MeshTags::RIGHT:  groups.xL.push_back(f.id); break;
            case MeshTags::BOTTOM: groups.y0.push_back(f.id); break;
            case MeshTags::TOP:    groups.yL.push_back(f.id); break;
            default:
                // 其他自定义边界（如障碍物）存入 map
                groups.custom[pid].push_back(f.id);
                break;
            }
        }
        return groups;
    }

    // =========================================================
    // 基于物理组 Tag 进行 3D 边界分类
    // =========================================================
    // 功能：根据 Face 的 physicalGroupId 将其分派到 FaceGroups 的 x0/xL/y0/yL/z0/zL 中
    // 依赖：MeshDefinitions.h 中的 Tag 定义以及 Mesh::GenerateMesh3D 中的标记逻辑
    inline FaceGroups ClassifyBoundaryFacesByTag_3D(const Mesh& mesh)
    {
        FaceGroups groups;
        const auto& faces = mesh.getFaces();

        for (const auto& f : faces)
        {
            // 1. 过滤：只有被 GenerateMesh3D 标记过的面才处理
            //    (内部面默认为 -1，未标记的边界也会是 -1)
            if (f.physicalGroupId == -1) continue;

            int pid = f.physicalGroupId;

            // 2. 映射：将 Tag ID 分派到对应的边界列表
            switch (pid)
            {
                // --- 侧面 (Sides) ---
            case MeshTags::LEFT:    groups.x0.push_back(f.id); break; // X=0
            case MeshTags::RIGHT:   groups.xL.push_back(f.id); break; // X=Lx
            case MeshTags::BOTTOM:  groups.y0.push_back(f.id); break; // Y=0
            case MeshTags::TOP:     groups.yL.push_back(f.id); break; // Y=Ly

                // --- 3D 特有底面与顶面 (Front/Back) ---
                // 对应 GenerateMesh3D 中: 
                // gmsh::model::geo::addPhysicalGroup(2, { baseSurface }, MeshTags::TAG_FRONT, "Z0_Face");
            case MeshTags::TAG_FRONT: groups.z0.push_back(f.id); break; // Z=0

                // 对应 GenerateMesh3D 中:
                // gmsh::model::geo::addPhysicalGroup(2, { topTag }, MeshTags::TAG_BACK, "ZL_Face");
            case MeshTags::TAG_BACK:  groups.zL.push_back(f.id); break; // Z=Lz
                // --- 其他自定义边界 (如内部障碍物) ---
            default:
                groups.custom[pid].push_back(f.id);
                break;
            }
        }
        return groups;
    }
}