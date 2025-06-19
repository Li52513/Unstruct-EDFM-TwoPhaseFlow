#include "Cell.h"
#include <cmath>

Cell::Cell(int id, const std::vector<int>& nodeIDs)
	: id(id), nodeIDs(nodeIDs), volume(0.0), center(0.0, 0.0, 0.0), sourceTerm(0.0), faceDiscreCoef(0.0), pressure(6.5e6), pressureGradient(0.0, 0.0, 0.0), temperature(483.15), saturation_water(0.8), saturation_CO2(0.2), error(1e-6), SolidMaterialProps(), WaterMaterialProps(), CO2MaterialProps()  // 初始化物性参数
{
    // 物性参数通过 materialProps 成员进行初始化
}


void Cell::computeCenterAndVolume(const map<int, Node>& allNodes) 
{
    if (nodeIDs.size() == 3) 
    {
        Vector p0 = allNodes.at(nodeIDs[0]).coord;
        Vector p1 = allNodes.at(nodeIDs[1]).coord;
        Vector p2 = allNodes.at(nodeIDs[2]).coord;

        center = (p0 + p1 + p2) / 3.0;
        Vector v1 = p1 - p0;
        Vector v2 = p2 - p0;
		volume = 0.5 * fabs(v1.m_x * v2.m_y - v1.m_y * v2.m_x); //fabs表示绝对值
    }
    else if (nodeIDs.size() == 4) 
    {
        // 3D 四面体单元
        Vector p0 = allNodes.at(nodeIDs[0]).coord;
        Vector p1 = allNodes.at(nodeIDs[1]).coord;
        Vector p2 = allNodes.at(nodeIDs[2]).coord;
        Vector p3 = allNodes.at(nodeIDs[3]).coord;

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



