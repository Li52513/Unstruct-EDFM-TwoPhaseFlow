/**
 * @file BoundaryConditionManager.h
 * @brief 边界条件统一管理器
 * @details
 * 它基于 Gmsh 的 Physical Tag (Int) 来管理边界条件，并提供统一的 a-b-c 系数供 FVM 离散化模块调用。
 * * 统一方程形式 (The a-b-c Framework):
 * a * Phi_b + b * Flux_b = c
 * * - Dirichlet (定值): a=1, b=0, c=Value
 * - Neumann (定流):   a=0, b=1, c=Flux_Value
 */

#pragma once

#include <map> 
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <functional>

 // 引入 Mesh 定义，确保 Tag ID 类型一致
#include "MeshDefinitions.h"
#include "UserDefineVarType.h"

namespace BoundarySetting {

    /**
     * @enum BoundaryType
     * @brief 边界条件类型枚举
     */

    enum class BoundaryType {
        Dirichlet,  ///< 第一类边界条件 (定值，如定压、定温)
        Neumann     ///< 第二类边界条件 (定流，如定流量、绝热)
    };

    /**
     * @struct BCCoefficients
     * @brief a-b-c 系数结构体
     * @details 用于 FVM 离散化组装时的统一接口
     */
    struct BCCoefficients {
        double a;               ///< 变量系数
        double b;               ///< 通量系数
		double c;               ///< 右端项值 (RHS)
		BoundaryType type;      ///< 记录原始类型用于调试或特殊处理

        BCCoefficients() : a(0), b(0), c(0), type(BoundaryType::Dirichlet) {}
        BCCoefficients(double _a, double _b, double _c, BoundaryType _t)
            : a(_a), b(_b), c(_c), type(_t) {
        }
    };

    /**
     * @brief 边界值计算函数签名
     * @param coords 面中心坐标
     * @return double 物理量的值 (Value or Flux)
     */
    using BoundaryValueFunc = std::function<double(const Vector&)>;

    /**
     * @class BoundaryConditionManager
     * @brief 边界条件管理器 (Generic Manager)
     */

    class BoundaryConditionManager {
    public:

        /**
         * @brief 构造函数
         */
        BoundaryConditionManager();

        /**
         * @brief 析构函数
         */
        ~BoundaryConditionManager();

        // =========================================================
        // 设置接口 (供初始化模块调用)
        // =========================================================

        /**
         * @brief 设置 Dirichlet 边界条件 (定值)
         * @param tagID Gmsh Physical Tag ID (定义在 MeshDefinitions.h)
         * @param value 边界值 (如压力 Pa, 温度 K)
         */
        void SetDirichletBC(int tagID, double value);

        /**
         * @brief [New] 设置线性分布 Dirichlet BC
         * @details Val = RefVal + Gradient * (Coord[axis] - RefCoord)
         * 例如：P(z) = 30MPa + 1e4 * (z - 0)
         * @param tagID 边界 Tag
         * @param refValue 参考点的值 (RefVal)
         * @param refCoord 参考点的坐标值 (RefCoord)
         * @param gradient 梯度 (Slope, e.g., Pa/m)
         * @param axis 变化轴: 0=X, 1=Y, 2=Z
         */
        void SetLinearDirichletBC(int tagID, double refValue, double refCoord, double gradient, int axis);

        /**
         * @brief 设置 Neumann 边界条件 (定流)
         * @param tagID Gmsh Physical Tag ID
         * @param fluxValue 通量值 (如流速 m/s 或热流密度 W/m2，视物理方程定义而定)
         * @details
         * 对于 Neumann 边界: 0 * Phi + 1 * Flux = c
         * 此时 a=0, b=1, c=fluxValue
         */
        void SetNeumannBC(int tagID, double fluxValue);

        // =========================================================
        // 查询接口 (供 Discretization 模块调用)
        // =========================================================

/**
         * @brief [Updated] 获取指定位置的 a-b-c 系数
         * @param tagID 边界 Tag
         * @param faceCenter 面中心坐标 (用于计算线性分布的具体值)
         */
        BCCoefficients GetBCCoefficients(int tagID, const Vector& faceCenter) const;
        /**
         * @brief 检查某个 Tag 是否已定义边界条件
         * @param tagID Gmsh Physical Tag ID
         * @return true 已定义, false 未定义
         */
        bool HasBC(int tagID) const;

        /**
         * @brief 清除所有边界条件
         */
        void Clear();

        /**
         * @brief 打印当前所有边界条件配置 (调试用)
         * @param name 物理场名称 (如 "Pressure")
         */
        void PrintSummary(const std::string& name) const;

    private:
        // 内部存储结构升级
        struct BCDefinition {
            double a;
            double b;
            BoundaryValueFunc valFunc; // [New] 存储计算逻辑，而非死值
            BoundaryType type;
            std::string desc; // 用于打印信息
        };

        /// @brief 存储 Tag 到 BC 系数的映射
        std::map<int, BCDefinition> bcRegistry_;
    };
}// namespace BoundarySetting