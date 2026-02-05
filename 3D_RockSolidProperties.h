#pragma once

#include <vector>
#include <string>
#include <memory>
#include <algorithm>

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "UserDefineVarType.h" 

/**
 * @struct RockPropertyParams
 * @brief 单一岩石类型的物性参数集合
 * @details 包含孔隙度、渗透率张量、热物性等。
 * 输入渗透率单位建议为 mD，计算时会自动转换为 m^2。
 */

struct RockPropertyParams
{
    double k_xx = 0.0;      ///< X方向渗透率 [mD]
    double k_yy = 0.0;      ///< Y方向渗透率 [mD]
    double k_zz = 0.0;      ///< Z方向渗透率 [mD]
    double phi = 0.0;       ///< 孔隙度 [-]
    double rho = 0.0;       ///< 岩石密度 [kg/m^3]
    double cp = 0.0;        ///< 岩石比热容 [J/(kg*K)]
    double lambda = 0.0;    ///< 岩石导热系数 [W/(m*K)]
    double c_r = 0.0;       ///< 岩石压缩系数 [1/Pa]

    // 默认构造
    RockPropertyParams() = default;

    // 全参构造
    RockPropertyParams(double kx, double ky, double kz, double p, double r, double c, double l, double cr)
        : k_xx(kx), k_yy(ky), k_zz(kz), phi(p), rho(r), cp(c), lambda(l), c_r(cr) {
    }
};

/**
 * @struct BoundingBox3D
 * @brief 3D 轴对齐包围盒 (AABB) 区域定义
 */
struct BoundingBox3D
{
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;

    bool contains(const Vector& p) const
    {
        return (p.m_x >= x_min && p.m_x <= x_max &&
            p.m_y >= y_min && p.m_y <= y_max &&
            p.m_z >= z_min && p.m_z <= z_max);
    }
};

/**
 * @struct RockRegion3D
 * @brief 描述一个具有特定物性的空间区域
 */
struct RockRegion3D
{
    std::string name;           ///< 区域名称 (e.g. "HighPermLayer")
    BoundingBox3D bounds;       ///< 空间范围
    RockPropertyParams props;   ///< 该区域的物性参数
};

/**
 * @class RockSolidProperties_3D
 * @brief 3D 基岩静态物性管理器
 * @details 负责管理非均质岩石属性，并将参数分发到 FieldManager 的场数据中。
 * 支持 "背景值 + 多个局部区域" 的赋值策略。
 */
class RockSolidProperties_3D
{
public:
    /**
     * @brief 构造函数
     * @param meshMgr 3D网格管理器 (用于获取单元中心坐标)
     * @param fieldMgr 3D场管理器 (用于写入物性场)
     */
    RockSolidProperties_3D(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr);

    /**
     * @brief 设置全场默认背景物性
     * @param params 默认岩石参数
     */
    void setBackgroundProperties(const RockPropertyParams& params);

    /**
     * @brief 添加一个局部物性区域 (后添加的区域优先级更高，即覆盖之前的)
     * @param name 区域名称
     * @param box 空间范围
     * @param params 区域岩石参数
     */
    void addRegion(const std::string& name, const BoundingBox3D& box, const RockPropertyParams& params);

    /**
     * @brief 执行物性初始化
     * @details 遍历所有基岩单元，根据其中心坐标确定所属区域，并将物性写入场。
     * 渗透率在此处由 mD 转换为 m^2。
     */
    void InitializeRockProperties();

private:
    const MeshManager_3D& meshMgr_;
    FieldManager_3D& fieldMgr_;

    RockPropertyParams backgroundProps_;        ///< 背景物性
    std::vector<RockRegion3D> regions_;         ///< 局部区域列表

    // 单位转换常数: 1 mD = 9.869233e-16 m^2
    // 工程上常取 1e-15 或精确值，这里采用精确值以保证严谨性
    static constexpr double MD_TO_M2 = 9.869233e-16;
};