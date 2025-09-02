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
#include "Initializer.h"
#include "FracIndex.h"   // 新增


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
    void classifyRockRegionsByGeometry(MeshManager& mgr, const vector<RegionGeometry>& regionGeomes, Cell::RegionType defaultRegion);

    void classifyFractureElementsByGeometry(MeshManager& mgr, int fracID, const Vector& regionStart, const Vector& regionEnd, FractureElementType insideType, FractureElementType outsideType);


	//************固相物性参数初始化与更新************//
    // 1)【初始化】—在开始计算前根据初始的 P 和 T 初始化固相物性
    void InitializeRockMatrixProperties(MeshManager& mgr, FieldRegistry& reg_r);
    void InitializeFractureElementsProperties(MeshManager& mgr, FieldRegistry& reg_fr);

    // 2) 【更新】—每一时间步迭代结束后，根据最新的 P 和 T 更新固相物性
    void UpdateMatrixProperties(MeshManager& mgr, FieldRegistry& reg_r);
    void UpdateFractureSolidProperties(MeshManager& mgr, FieldRegistry& reg_fr);


    //【更新】更新基岩单元的分类（低/中/高渗透-孔隙）
    void classifyRockRegions(MeshManager& mgr, double poroLow, double poroHigh, double permLow, double permHigh);

    //【更新】更新裂缝单元的分类（低/中/高渗透-孔隙）
    void classifyFractureElements(MeshManager& mgr, double permThreshold);

    // 3） 【调试】 打印所有基岩单元和所有裂缝段的物性参数（待更新）
    void debugPrintProperties(MeshManager& mgr,
        const FieldRegistry& reg,
        const FieldRegistry& reg_fr,
        std::size_t maxPrint = 20) const;


	//**************流体物性参数初始化与更新**************//

    void InitializeMatrixFluidProperties(MeshManager& mgr,FieldRegistry& reg,const VGParams& vg);
	void InitializeFractureFluidProperties(MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr, const VGParams& vg);


    void UpdateMatrixFluidProperties(MeshManager& mgr,FieldRegistry& reg,const VGParams& vg);
    void UpdateFractureFluidProperties(MeshManager& mgr, FieldRegistry& reg, FieldRegistry& reg_fr, const VGParams& vg);






};