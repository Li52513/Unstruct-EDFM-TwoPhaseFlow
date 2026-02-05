#pragma once
// 定义物理组 ID，确保 Mesh 和 BC 设置使用相同的标准
namespace MeshTags {
    const int LEFT = 1; // 对应 x0
    const int RIGHT = 2; // 对应 xL
    const int BOTTOM = 3; // 对应 y0
    const int TOP = 4; // 对应 yL
    const int TAG_FRONT = 5; // 3D 底面
    const int TAG_BACK = 6; // 3D 顶面
    const int OBSTACLE = 7; // 示例：内部障碍物
    const int FLUID = 10;
}