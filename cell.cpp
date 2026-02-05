#include "Cell.h"
#include "GeometryCalculate.h"
#include <cmath>
#include <limits>
#include <algorithm>

//构造函数
Cell::Cell(int id, const std::vector<int>& nodeIDs)
	: id(id), CellNodeIDs(nodeIDs), volume(0.0), center(0.0, 0.0, 0.0) /*pressure(6.531e7), pressureGradient(0.0, 0.0, 0.0), temperature(597.65), saturation_water(0.8), saturation_CO2(0.2),*/ /*error(1e-6)*/ // 初始化物性参数
{
    // 物性参数通过 materialProps 成员进行初始化
}

// =========================================================
// 【2D】找到组成2D网格cell的面
// =========================================================
vector<vector<int>> Cell::getLocalFaces_2D(const std::unordered_map<int, Node>& allNodes) const
{
    const auto& cn = CellNodeIDs;
    std::vector<std::vector<int>> faces;

    auto is2DCell = [&](double eps = 1e-12)->bool
    {
        double zmin = 1e300, zmax = -1e300;
        for (int nid : cn) {
            const Vector& p = allNodes.at(nid).coord;
            zmin = std::min(zmin, p.m_z);
            zmax = std::max(zmax, p.m_z);
        }
        return (zmax - zmin) < eps; 
    };
    // 2D：边（按环连接）
    if ((cn.size() == 3 || cn.size() == 4) && is2DCell()) {
        const int m = static_cast<int>(cn.size());
        faces.reserve(m);
        for (int i = 0; i < m; ++i) {
            int j = (i + 1) % m;
            faces.push_back({ cn[i], cn[j] });
        }
    }
    return faces;
}

// =========================================================
// 【3D】找到组成3D网格cell的面，四面体4nodes / 棱柱Prism 6nodes/ 金字塔Pyramid5 nodes / 六面体Hexahedron 8nodes
// =========================================================
vector<vector<int>> Cell::getLocalFaces_3D() const
{
    const auto& cn = CellNodeIDs;
    std::vector<std::vector<int>> faces;

    // 3D：典型单元
    if (cn.size() == 4)
    { // tetra
        faces = {
            {cn[0], cn[1], cn[2]},
            {cn[0], cn[1], cn[3]},
            {cn[1], cn[2], cn[3]},
            {cn[0], cn[2], cn[3]}
        };
    }
    else if (cn.size() == 5)
    { // pyramid
        faces = {
            {cn[0], cn[1], cn[2], cn[3]},
            {cn[0], cn[1], cn[4]},
            {cn[1], cn[2], cn[4]},
            {cn[2], cn[3], cn[4]},
            {cn[3], cn[0], cn[4]}
        };
    }
    else if (cn.size() == 6)
    { // prism
        faces = {
            {cn[0], cn[1], cn[2]},
            {cn[3], cn[4], cn[5]},
            {cn[0], cn[1], cn[4], cn[3]},
            {cn[1], cn[2], cn[5], cn[4]},
            {cn[2], cn[0], cn[3], cn[5]}
        };
    }
    else if (cn.size() == 8)
    { // hexa
        faces =
        {
            {cn[0], cn[1], cn[2], cn[3]},
            {cn[4], cn[5], cn[6], cn[7]},
            {cn[0], cn[1], cn[5], cn[4]},
            {cn[1], cn[2], cn[6], cn[5]},
            {cn[2], cn[3], cn[7], cn[6]},
            {cn[3], cn[0], cn[4], cn[7]}
        };
    }
    return faces;
}

// =========================================================
// 【2D】计算平面单元的面积与质心
// =========================================================
void Cell::computeCenterAndVolume_2D(const std::unordered_map<int, Node>& allNodes)
{
    // 定义 Z 轴容差
    const double Z_FLAT_TOL = 1e-12;
    const auto& cn = CellNodeIDs;

    // 内部 Lambda: 快速检查是否共面 (Flat check)
    // 只有 Z 值几乎一致才视为 2D 单元处理
    auto is2D = [&]() -> bool {
        if (cn.empty()) return false;
        double zmin = 1e300, zmax = -1e300;
        for (int nid : cn) {
            const double z = allNodes.at(nid).coord.m_z;
            zmin = std::min(zmin, z);
            zmax = std::max(zmax, z);
        }
        return (zmax - zmin) < Z_FLAT_TOL;
        };

    // 逻辑：必须是三角形(3)或四边形(4)，且 Z 轴高度差极小
    if ((cn.size() == 3 || cn.size() == 4) && is2D())
    {
        double A = 0.0;
        Vector C(0, 0, 0);

        // 调用您的辅助函数 (内部调用 GeomCalculate::polygonCentroid)
        calculatePolygonProps(cn, allNodes, A, C);

        this->center = C;
        this->volume = A; // 2D 单元将面积存储在 volume 字段中
    }
    else
    {
        // 如果调用了 2D 函数但单元实际上是 3D 的，或者节点数不对，可以在此报错或置零
        std::cerr << "Warning: computeCenterAndVolume (2D) called on non-flat or invalid cell.\n";
        this->volume = 0.0;
        this->center = Vector(0, 0, 0);
    }
}

// =========================================================
// 【3D】辅助函数：计算多边形的面积与质心 (基于 2D 投影)
// =========================================================
void Cell::calculatePolygonProps(const std::vector<int>& ids,
    const std::unordered_map<int, Node>& allNodes, double& outArea, Vector& outCentroid)
{
    // 1. 数据转换：从 ID 提取坐标，同时手动计算 Z 总和
    std::vector<Vector> points;
    points.reserve(ids.size());

    double zSum = 0.0;
    for (int id : ids) {
        const Vector& v = allNodes.at(id).coord;
        points.push_back(v);
        zSum += v.m_z;
    }

    // 2. 调用库函数计算 2D 投影面积和 XY 质心
    // 注意：这里的 cent.z 会是 0.0
    Vector cent = GeomCalculate::polygonCentroid(points, outArea);

    // 3. 补全 Z 轴逻辑
    // 将库函数计算的 XY 与我们手动计算的 Z 平均值结合
    if (!points.empty()) { cent.m_z = zSum / points.size(); }

    // 4. 输出结果
    outCentroid = cent;
}
// =========================================================
// 【3D】计算体单元的体积与质心 (基于挤出假设)
// =========================================================
void Cell::computeCenterAndVolume_3D(const std::unordered_map<int, Node>& allNodes)
{
    const auto& cn = CellNodeIDs;
    const size_t n = cn.size();

    // ---------- Case 1: 四面体 (Tetrahedron) ----------
    // 四面体不使用投影面积，而是使用混合积解析解
    if (n == 4) {
        const Vector& p0 = allNodes.at(cn[0]).coord;
        const Vector& p1 = allNodes.at(cn[1]).coord;
        const Vector& p2 = allNodes.at(cn[2]).coord;
        const Vector& p3 = allNodes.at(cn[3]).coord;

        this->center = (p0 + p1 + p2 + p3) / 4.0;

        Vector v1 = p0 - p3;
        Vector v2 = p1 - p3;
        Vector v3 = p2 - p3;
        // 注意：保留原有的运算符重载逻辑 (&为叉乘, *为点乘)
        this->volume = std::fabs((v1 & v2) * v3) / 6.0;
        return;
    }

    // ---------- Case 2: 金字塔 (Pyramid - 5 nodes) ----------
    else if (n == 5) {
        // 拓扑假设: 底面(0,1,2,3) + 顶点(4)
        std::vector<int> base = { cn[0], cn[1], cn[2], cn[3] };
        const Vector& apex = allNodes.at(cn[4]).coord;

        double A_base = 0.0;
        Vector C_base(0, 0, 0);

        // 使用辅助函数计算底面投影性质
        calculatePolygonProps(base, allNodes, A_base, C_base);

        // 计算底面平均 Z
        double zbase = 0.0;
        for (int id : base) zbase += allNodes.at(id).coord.m_z;
        zbase /= 4.0;

        // 计算高度 (挤出高度)
        double h = apex.m_z - zbase;

        this->volume = (A_base * std::fabs(h)) / 3.0;

        // 质心 Z 轴修正: 底面 Z + 1/4 高度 (重心性质)
        // 质心 XY: 使用底面形心坐标近似
        this->center = Vector(C_base.m_x, C_base.m_y, zbase + h * 0.25);
        return;
    }

    // ---------- Case 3: 棱柱 (Prism - 6 nodes) ----------
    else if (n == 6) {
        // 拓扑假设: 底三角(0,1,2) + 顶三角(3,4,5)
        std::vector<int> bot = { cn[0], cn[1], cn[2] };

        double A_bot = 0.0;
        Vector C_bot(0, 0, 0);

        // 替换：计算底面
        calculatePolygonProps(bot, allNodes, A_bot, C_bot);

        // 计算平均高度 (对应节点 Z 差值的平均)
        double h = 0.0;
        h += allNodes.at(cn[3]).coord.m_z - allNodes.at(cn[0]).coord.m_z;
        h += allNodes.at(cn[4]).coord.m_z - allNodes.at(cn[1]).coord.m_z;
        h += allNodes.at(cn[5]).coord.m_z - allNodes.at(cn[2]).coord.m_z;
        h /= 3.0;

        this->volume = A_bot * std::fabs(h);
        // 质心: 底面形心 + 半高
        this->center = Vector(C_bot.m_x, C_bot.m_y, C_bot.m_z + h * 0.5);
        return;
    }

    // ---------- Case 4: 六面体 (Hexahedron - 8 nodes) ----------
    else if (n == 8) {
        // 拓扑假设: 底四边形(0..3) + 顶四边形(4..7)
        std::vector<int> bot = { cn[0], cn[1], cn[2], cn[3] };

        double A_bot = 0.0;
        Vector C_bot(0, 0, 0);

        // 替换：计算底面
        calculatePolygonProps(bot, allNodes, A_bot, C_bot);

        // 计算平均高度
        double h = 0.0;
        h += allNodes.at(cn[4]).coord.m_z - allNodes.at(cn[0]).coord.m_z;
        h += allNodes.at(cn[5]).coord.m_z - allNodes.at(cn[1]).coord.m_z;
        h += allNodes.at(cn[6]).coord.m_z - allNodes.at(cn[2]).coord.m_z;
        h += allNodes.at(cn[7]).coord.m_z - allNodes.at(cn[3]).coord.m_z;
        h /= 4.0;

        this->volume = A_bot * std::fabs(h);
        this->center = Vector(C_bot.m_x, C_bot.m_y, C_bot.m_z + h * 0.5);
        return;
    }

    // ---------- Case 5: 不支持的类型 ----------
    std::cerr << "不支持的 3D 单元类型，无法计算几何中心和体积。 (nodes=" << n << ")\n";
}

// =========================================================
// 【3D】计算体单元的AABB 后期由于构建 3D 均匀背景网格索引
// =========================================================
void Cell::computeAABB(const std::unordered_map<int, Node>& allNodes)
{
    if (CellNodeIDs.empty()) return;

    Vector minP(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
    Vector maxP(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

    for (int nodeID : CellNodeIDs)
    {
        // 安全检查：确保节点存在
        auto it = allNodes.find(nodeID);
        if (it != allNodes.end())
        {
            const Vector& p = it->second.coord;
            if (p.m_x < minP.m_x) minP.m_x = p.m_x;
            if (p.m_y < minP.m_y) minP.m_y = p.m_y;
            if (p.m_z < minP.m_z) minP.m_z = p.m_z;

            if (p.m_x > maxP.m_x) maxP.m_x = p.m_x;
            if (p.m_y > maxP.m_y) maxP.m_y = p.m_y;
            if (p.m_z > maxP.m_z) maxP.m_z = p.m_z;
        }
    }

    // 稍微外扩一点点，避免浮点数精度的边界接触问题
    double tolerance = 1e-6;
    minP = minP - Vector(tolerance, tolerance, tolerance);
    maxP = maxP + Vector(tolerance, tolerance, tolerance);

    this->boundingBox = AABB(minP, maxP);
}

