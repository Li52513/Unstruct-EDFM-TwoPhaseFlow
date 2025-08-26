#include "Cell.h"
#include <cmath>

Cell::Cell(int id, const std::vector<int>& nodeIDs)
	: id(id), CellNodeIDs(nodeIDs), volume(0.0), center(0.0, 0.0, 0.0), sourceTerm(0.0), faceDiscreCoef(0.0), /*pressure(6.531e7), pressureGradient(0.0, 0.0, 0.0), temperature(597.65), saturation_water(0.8), saturation_CO2(0.2),*/ error(1e-6), SolidMaterialProps(), WaterMaterialProps(), CO2MaterialProps()  // 初始化物性参数
{
    // 物性参数通过 materialProps 成员进行初始化
}

vector<vector<int>> Cell::getLocalFaces()const
{
	const auto& cn = CellNodeIDs; // 获取单元的节点编号
    vector<vector<int>> faces;
    if (cn.size() == 3) {
        // 2D三角形
        faces = {
            {cn[0], cn[1]},
            {cn[1], cn[2]},
            {cn[2], cn[0]}
        };
    }
    else if (cn.size() == 4) {
        // 3D四面体
        faces = {
            {cn[0], cn[1], cn[2]},
            {cn[0], cn[1], cn[3]},
            {cn[1], cn[2], cn[3]},
            {cn[0], cn[2], cn[3]}
        };
    }
    else if (cn.size() == 5) {
        // 金字塔
        faces = {
            {cn[0], cn[1], cn[2], cn[3]}, // 底面
            {cn[0], cn[1], cn[4]},
            {cn[1], cn[2], cn[4]},
            {cn[2], cn[3], cn[4]},
            {cn[3], cn[0], cn[4]}
        };
    }
    else if (cn.size() == 6) {
        // 棱柱
        faces = {
            {cn[0], cn[1], cn[2]},
            {cn[3], cn[4], cn[5]},
            {cn[0], cn[1], cn[4], cn[3]},
            {cn[1], cn[2], cn[5], cn[4]},
            {cn[2], cn[0], cn[3], cn[5]}
        };
    }
    else if (cn.size() == 8) {
        // 六面体
        faces = {
            {cn[0], cn[1], cn[2], cn[3]},
            {cn[4], cn[5], cn[6], cn[7]},
            {cn[0], cn[1], cn[5], cn[4]},
            {cn[1], cn[2], cn[6], cn[5]},
            {cn[2], cn[3], cn[7], cn[6]},
            {cn[3], cn[0], cn[4], cn[7]}
        };
    }
    // 可继续扩展其它多面体类型
    return faces;
}

void Cell::computeCenterAndVolume(const unordered_map<int, Node>& allNodes)
{
    if (CellNodeIDs.size() == 3) 
    {
        Vector p0 = allNodes.at(CellNodeIDs[0]).coord;
        Vector p1 = allNodes.at(CellNodeIDs[1]).coord;
        Vector p2 = allNodes.at(CellNodeIDs[2]).coord;

        center = (p0 + p1 + p2) / 3.0;
		//cout << "Cell " << id << " center: [" << center.m_x <<"," << center.m_y << "," << center.m_z << "]" << endl;
        Vector v1 = p1 - p0;
        Vector v2 = p2 - p0;
		volume = 0.5 * fabs(v1.m_x * v2.m_y - v1.m_y * v2.m_x); //fabs表示绝对值
		//cout << "Cell " << id << " area: " << volume << endl; // 输出面积
    }
    else if (CellNodeIDs.size() == 4) 
    {
        // 3D 四面体单元
        Vector p0 = allNodes.at(CellNodeIDs[0]).coord;
        Vector p1 = allNodes.at(CellNodeIDs[1]).coord;
        Vector p2 = allNodes.at(CellNodeIDs[2]).coord;
        Vector p3 = allNodes.at(CellNodeIDs[3]).coord;

        center = (p0 + p1 + p2 + p3) / 4.0;

        // 使用混合积法计算体积：V = |(a−d)·((b−d)×(c−d))| / 6
        Vector a = p0, b = p1, c = p2, d = p3;
        Vector v1 = a - d;
        Vector v2 = b - d;
        Vector v3 = c - d;
        volume = fabs((v1 & v2) * v3) / 6.0; //  &是叉乘，* 是点乘
	}
    else
    {
        cerr << "不支持的单元类型，无法计算几何中心和面积。" << endl;
        
    }

}



