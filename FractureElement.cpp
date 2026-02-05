#include "FractureElement.h"

//有参构造函数实现
FractureElement::FractureElement(int id_, int cellID_, double len, double avgDist)
	: id(id_)
	, cellID(cellID_)
	, length(len)
	, avgDistance(avgDist)
{
    // 初始化关键拓扑索引为无效值，防止脏数据
    solverIndex = -1;
    parentFractureID = -1;
    gIDstart = -1;
    gIDend = -1;

    // 初始化其他物理参数
    aperture = 1.0;
    geomCI = 0.0;
    geomAlpha = 0.0;
    alpha_fr = 0.0;
}