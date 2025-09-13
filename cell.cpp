#include "Cell.h"
#include <cmath>

Cell::Cell(int id, const std::vector<int>& nodeIDs)
	: id(id), CellNodeIDs(nodeIDs), volume(0.0), center(0.0, 0.0, 0.0), sourceTerm(0.0), faceDiscreCoef(0.0), /*pressure(6.531e7), pressureGradient(0.0, 0.0, 0.0), temperature(597.65), saturation_water(0.8), saturation_CO2(0.2),*/ error(1e-6), SolidMaterialProps(), WaterMaterialProps(), CO2MaterialProps()  // 初始化物性参数
{
    // 物性参数通过 materialProps 成员进行初始化
}

vector<vector<int>> Cell::getLocalFaces(const std::unordered_map<int, Node>& allNodes) const
{
    const auto& cn = CellNodeIDs;
    std::vector<std::vector<int>> faces;

    auto is2DCell = [&](double eps = 1e-12)->bool {
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
        return faces;
    }

    // 3D：典型单元
    if (cn.size() == 4) { // tetra
        faces = {
            {cn[0], cn[1], cn[2]},
            {cn[0], cn[1], cn[3]},
            {cn[1], cn[2], cn[3]},
            {cn[0], cn[2], cn[3]}
        };
    }
    else if (cn.size() == 5) { // pyramid
        faces = {
            {cn[0], cn[1], cn[2], cn[3]},
            {cn[0], cn[1], cn[4]},
            {cn[1], cn[2], cn[4]},
            {cn[2], cn[3], cn[4]},
            {cn[3], cn[0], cn[4]}
        };
    }
    else if (cn.size() == 6) { // prism
        faces = {
            {cn[0], cn[1], cn[2]},
            {cn[3], cn[4], cn[5]},
            {cn[0], cn[1], cn[4], cn[3]},
            {cn[1], cn[2], cn[5], cn[4]},
            {cn[2], cn[0], cn[3], cn[5]}
        };
    }
    else if (cn.size() == 8) { // hexa
        faces = {
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

void Cell::computeCenterAndVolume(const std::unordered_map<int, Node>& allNodes)
{
    const auto& cn = CellNodeIDs;
    auto polygonAreaCentroid2D = [&](const std::vector<int>& ids, double& A, Vector& C)
        {
            // 2D 多边形（按顺/逆时针）面积与质心，忽略 z
            const size_t n = ids.size();
            double area2 = 0.0, Cx = 0.0, Cy = 0.0;
            double zavg = 0.0;
            for (size_t i = 0; i < n; ++i) {
                const Vector& Pi = allNodes.at(ids[i]).coord;
                const Vector& Pj = allNodes.at(ids[(i + 1) % n]).coord;
                const double cross = Pi.m_x * Pj.m_y - Pj.m_x * Pi.m_y;
                area2 += cross;
                Cx += (Pi.m_x + Pj.m_x) * cross;
                Cy += (Pi.m_y + Pj.m_y) * cross;
                zavg += Pi.m_z;
            }
            double A_signed = 0.5 * area2;     // 可正可负
            A = std::fabs(A_signed);
            if (std::fabs(A_signed) < 1e-30) { C = Vector(0, 0, 0); return; }
            const double inv = 1.0 / (6.0 * A_signed);
            C = Vector(Cx * inv, Cy * inv, zavg / n); // z 取顶点平均
        };

    auto is2D = [&]() 
    {
        double zmin = 1e300, zmax = -1e300;
        for (int nid : cn) {
            const double z = allNodes.at(nid).coord.m_z;
            zmin = std::min(zmin, z);
            zmax = std::max(zmax, z);
        }
        return (zmax - zmin) < 1e-12;
    };

    // ---------- 2D：三角/四边形（面积与质心） ----------
    if (is2D() && (cn.size() == 3 || cn.size() == 4))
    {
        double A = 0.0; Vector C(0, 0, 0);
        polygonAreaCentroid2D(cn, A, C);
        center = C;
        volume = A; // 2D 用面积作为 "volume"
        return;
    }
    // ---------- 3D：典型单元 ----------
    if (cn.size() == 4) {
        // 四面体：混合积
        const Vector& p0 = allNodes.at(cn[0]).coord;
        const Vector& p1 = allNodes.at(cn[1]).coord;
        const Vector& p2 = allNodes.at(cn[2]).coord;
        const Vector& p3 = allNodes.at(cn[3]).coord;
        center = (p0 + p1 + p2 + p3) / 4.0;
        Vector v1 = p0 - p3, v2 = p1 - p3, v3 = p2 - p3;
        volume = std::fabs((v1 & v2) * v3) / 6.0; // & 叉乘，* 点乘
        return;
    }
    else if (cn.size() == 5) {
        // 金字塔：底面四边形 + 顶点（默认 cn[4] 为 apex）
        std::vector<int> base = { cn[0], cn[1], cn[2], cn[3] };
        const Vector& apex = allNodes.at(cn[4]).coord;

        double A = 0.0; Vector Cb(0, 0, 0);
        polygonAreaCentroid2D(base, A, Cb);

        double zbase = 0.0; for (int id : base) zbase += allNodes.at(id).coord.m_z; zbase /= 4.0;
        double h = apex.m_z - zbase; // 与 z 轴一致的挤出假设

        volume = (A * std::fabs(h)) / 3.0;
        center = Vector(Cb.m_x, Cb.m_y, (zbase + apex.m_z * 1.0) / 2.0); // 近似：介于1/4~1/2处；若需精确：z = zbase + h*1/4
        center.m_z = zbase + h * 0.25; // 精确：金字塔质心距底面 1/4 高
        return;
    }
    else if (cn.size() == 6) {
        // 棱柱：底三角(0,1,2) + 顶三角(3,4,5)，来自直线挤出
        std::vector<int> bot = { cn[0], cn[1], cn[2] };
        std::vector<int> top = { cn[3], cn[4], cn[5] };

        double Ab = 0.0; Vector Cb(0, 0, 0);
        polygonAreaCentroid2D(bot, Ab, Cb);

        double h = 0.0;
        h += allNodes.at(cn[3]).coord.m_z - allNodes.at(cn[0]).coord.m_z;
        h += allNodes.at(cn[4]).coord.m_z - allNodes.at(cn[1]).coord.m_z;
        h += allNodes.at(cn[5]).coord.m_z - allNodes.at(cn[2]).coord.m_z;
        h /= 3.0;

        volume = Ab * std::fabs(h);
        center = Vector(Cb.m_x, Cb.m_y, Cb.m_z + h * 0.5);
        return;
    }
    else if (cn.size() == 8) {
        // 六面体：底四边形(0..3) + 顶(4..7)，来自直线挤出
        std::vector<int> bot = { cn[0], cn[1], cn[2], cn[3] };

        double Ab = 0.0; Vector Cb(0, 0, 0);
        polygonAreaCentroid2D(bot, Ab, Cb);

        double h = 0.0;
        h += allNodes.at(cn[4]).coord.m_z - allNodes.at(cn[0]).coord.m_z;
        h += allNodes.at(cn[5]).coord.m_z - allNodes.at(cn[1]).coord.m_z;
        h += allNodes.at(cn[6]).coord.m_z - allNodes.at(cn[2]).coord.m_z;
        h += allNodes.at(cn[7]).coord.m_z - allNodes.at(cn[3]).coord.m_z;
        h /= 4.0;

        volume = Ab * std::fabs(h);
        center = Vector(Cb.m_x, Cb.m_y, Cb.m_z + h * 0.5);
        return;
    }

    // 若有其它多面体类型，建议用“通用多面体体积/质心”算法（基于三角面拼合），此处略。
    std::cerr << "不支持的单元类型，无法计算几何中心和体积/面积。 (nodes=" << cn.size() << ")\n";
        
}


//vector<vector<int>> Cell::getLocalFaces()const
//{
//	const auto& cn = CellNodeIDs; // 获取单元的节点编号
//    vector<vector<int>> faces;
//
//    if (cn.size() == 3) {
//        // 2D三角形
//        faces = {
//            {cn[0], cn[1]},
//            {cn[1], cn[2]},
//            {cn[2], cn[0]}
//        };
//    }
//    else if (cn.size() == 4) {
//        // 3D四面体
//        faces = {
//            {cn[0], cn[1], cn[2]},
//            {cn[0], cn[1], cn[3]},
//            {cn[1], cn[2], cn[3]},
//            {cn[0], cn[2], cn[3]}
//        };
//    }
//    else if (cn.size() == 5) {
//        // 金字塔
//        faces = {
//            {cn[0], cn[1], cn[2], cn[3]}, // 底面
//            {cn[0], cn[1], cn[4]},
//            {cn[1], cn[2], cn[4]},
//            {cn[2], cn[3], cn[4]},
//            {cn[3], cn[0], cn[4]}
//        };
//    }
//    else if (cn.size() == 6) {
//        // 棱柱
//        faces = {
//            {cn[0], cn[1], cn[2]},
//            {cn[3], cn[4], cn[5]},
//            {cn[0], cn[1], cn[4], cn[3]},
//            {cn[1], cn[2], cn[5], cn[4]},
//            {cn[2], cn[0], cn[3], cn[5]}
//        };
//    }
//    else if (cn.size() == 8) {
//        // 六面体
//        faces = {
//            {cn[0], cn[1], cn[2], cn[3]},
//            {cn[4], cn[5], cn[6], cn[7]},
//            {cn[0], cn[1], cn[5], cn[4]},
//            {cn[1], cn[2], cn[6], cn[5]},
//            {cn[2], cn[3], cn[7], cn[6]},
//            {cn[3], cn[0], cn[4], cn[7]}
//        };
//    }
//    // 可继续扩展其它多面体类型
//    return faces;
//}

//void Cell::computeCenterAndVolume(const unordered_map<int, Node>& allNodes)
//{
//    if (CellNodeIDs.size() == 3) 
//    {
//        Vector p0 = allNodes.at(CellNodeIDs[0]).coord;
//        Vector p1 = allNodes.at(CellNodeIDs[1]).coord;
//        Vector p2 = allNodes.at(CellNodeIDs[2]).coord;
//
//        center = (p0 + p1 + p2) / 3.0;
//		//cout << "Cell " << id << " center: [" << center.m_x <<"," << center.m_y << "," << center.m_z << "]" << endl;
//        Vector v1 = p1 - p0;
//        Vector v2 = p2 - p0;
//		volume = 0.5 * fabs(v1.m_x * v2.m_y - v1.m_y * v2.m_x); //fabs表示绝对值
//		//cout << "Cell " << id << " area: " << volume << endl; // 输出面积
//    }
//    else if (CellNodeIDs.size() == 4) 
//    {
//        // 3D 四面体单元
//        Vector p0 = allNodes.at(CellNodeIDs[0]).coord;
//        Vector p1 = allNodes.at(CellNodeIDs[1]).coord;
//        Vector p2 = allNodes.at(CellNodeIDs[2]).coord;
//        Vector p3 = allNodes.at(CellNodeIDs[3]).coord;
//
//        center = (p0 + p1 + p2 + p3) / 4.0;
//
//        // 使用混合积法计算体积：V = |(a−d)·((b−d)×(c−d))| / 6
//        Vector a = p0, b = p1, c = p2, d = p3;
//        Vector v1 = a - d;
//        Vector v2 = b - d;
//        Vector v3 = c - d;
//        volume = fabs((v1 & v2) * v3) / 6.0; //  &是叉乘，* 是点乘
//	}
//    else
//    {
//        cerr << "不支持的单元类型，无法计算几何中心和面积。" << endl;
//        
//    }
//
//}



