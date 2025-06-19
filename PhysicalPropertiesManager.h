#pragma once
#include <vector>
#include <unordered_map>
#include <utility>
#include "FractureTypes.h" // FractureElement::Type
#include "Cell.h"  // for Cell::RegionType
#include "PropertiesSummary.h"
#include "WaterPropertyTable.h"
#include "CO2PropertyTable.h"
#include "UserDefineVarType.h"  // for Vector


using namespace std;
class MeshManager;  // 前向声明 用于后面调用


class PhysicalPropertiesManager 
{
public:
    
    struct RegionGeometry
    {
        Cell::RegionType      type;     ///< 要贴的区域标签
		vector<Vector>   vertices; ///< 凸多边形顶点列表（至少3个点）
    };
    void classifyRockRegionsByGeometry(MeshManager& mgr, const vector<RegionGeometry>& regionGeomes, Cell::RegionType defaultRegion );

	void classifyFractureElementsByGeometry(MeshManager& mgr, int fracID,  const Vector& regionStart, const Vector& regionEnd, FractureElementType insideType, FractureElementType outsideType);

    
    
    // 1)【初始化】—在开始计算前根据初始的 P 和 T 初始化物性
    void InitializeRockMatrixProperties(MeshManager& mgr);
    void InitializeFractureElementsProperties(MeshManager& mgr);

	// 2) 【更新】—每一时间步迭代结束后，根据最新的 P 和 T 更新物性
	void UpdateMatrixProperties(MeshManager& mgr);
	void UpdateFractureSolidProperties(MeshManager& mgr);
    //【更新】更新基岩单元的分类（低/中/高渗透-孔隙）
    void classifyRockRegions(MeshManager& mgr, double poroLow, double poroHigh, double permLow, double permHigh);

    //【更新】更新裂缝单元的分类（低/中/高渗透-孔隙）
    void classifyFractureElements(MeshManager& mgr, double permThreshold);

    // 3） 【调试】 打印所有基岩单元和所有裂缝段的物性参数
    void debugPrintProperties( MeshManager& mgr) const;

};