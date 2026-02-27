/**
 * @file 2D_RockSolidProperties.h
 * @brief 2D 基岩静态物理属性管理子模块
 * @details 负责 2D-EDFM 框架下基岩网格 (Matrix) 物理属性的初始化与内存分配。
 * 采用状态机模式 (Stateful Builder Pattern)，支持背景参数设置与 AABB 局部区域非均质覆盖。
 */

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <omp.h>

#include "MeshManager.h"
#include "2D_FieldManager.h"
#include "PropertiesSummary.h"
#include "SolverContrlStrName_op.h"
#include "AABB.h"

 /**
  * @class RockSolidProperties_2D
  * @brief 2D 基岩固体属性管理器
  */
class RockSolidProperties_2D
{
public:
    /**
     * @brief 构造函数
     * @param meshMgr 2D网格管理器常量引用
     * @param fieldMgr 2D场管理器引用
     */
    RockSolidProperties_2D(const MeshManager& meshMgr, FieldManager_2D& fieldMgr);

    ~RockSolidProperties_2D() = default;

    /**
     * @brief 设置基岩背景 (全局默认) 物理属性
     * @param bg 包含渗透率、孔隙度、密度等参数的结构体
     */
    void setBackgroundProperties(const SolidProperties_RockMatrix& bg);

    /**
     * @brief 添加局部非均质区域
     * @details 后添加的区域具有更高的覆盖优先级。
     * @param name 区域自定义名称
     * @param box 定义该区域空间范围的 2D AABB 包围盒
     * @param props 该区域内的专属物理属性
     */
    void addRegion(const std::string& name, const AABB& box, const SolidProperties_RockMatrix& props);

    /**
     * @brief 执行基岩物理属性初始化
     * @details 遍历所有基岩网格，根据几何位置与 AABB 的包含关系，将对应属性写入 FieldManager。
     */
    void InitializeRockProperties();

private:
    const MeshManager& meshMgr_;   ///< 拓扑网格管理器引用
    FieldManager_2D& fieldMgr_;    ///< 物理场数据管理器引用

    SolidProperties_RockMatrix bg_; ///< 背景物理属性
    std::vector<std::pair<std::string, std::pair<AABB, SolidProperties_RockMatrix>>> regions_; ///< 局部区域集合
};