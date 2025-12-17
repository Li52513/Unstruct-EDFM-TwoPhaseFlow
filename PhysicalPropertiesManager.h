#pragma once
#include <vector>
#include <unordered_map>
#include <utility>
#include "Cell.h"  // for Cell::RegionType
#include "PropertiesSummary.h"
#include  "FractureSolidProperties.h"
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

    // 3） 【调试】 打印所有基岩单元和所有裂缝段的物性参数（待更新）
    void debugPrintProperties(MeshManager& mgr,
        const FieldRegistry& reg,
        const FieldRegistry& reg_fr,
        std::size_t maxPrint = 20) const;

    // ===== 统一的“按字段名更新”方法（推荐在外迭代/时间步里用） =====
    //基岩-固相参数更新
    // void UpdateMatrixRockAt(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field);
	 bool UpdateRockProperties(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field); //调用之前需要对基岩区域进行划分，不冉默认是Medium区域

     //基岩-水相物性参数更新
	 bool UpdateWaterProperties(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field);

	 //基岩-二氧化碳物性参数更新
	 bool UpdateCO2Properties(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field);

	 //单相有效热物性参数更新
     //水相
     bool ComputeEffectiveThermalProperties_Water(MeshManager& mgr, FieldRegistry& reg);
	 //二氧化碳相
	 bool ComputeEffectiveThermalProperties_CO2(MeshManager& mgr, FieldRegistry& reg);

	 //单相密度对压力导数更新
	 //水相
	 bool ComputeDrho_dp_Water(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field);
	 //二氧化碳相
	 bool ComputeDrho_dp_CO2(MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field);


	//裂缝-固相
	void UpdateFractureRockAt(MeshManager& mgr, FieldRegistry& reg_fr, FieldRegistry& reg, const std::string& pf_field, const std::string& Tf_field);

	//裂缝-流体相：phase ∈ {"water","co2","both"}
	void UpdateFractureFluidAt(MeshManager& mgr, FieldRegistry& reg_fr, FieldRegistry& reg, const std::string& pf_field, const std::string& Tf_field, const std::string& phase);
 





};