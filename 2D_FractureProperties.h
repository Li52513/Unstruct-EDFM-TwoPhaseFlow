/**
 * @file 2D_FractureProperties.h
 * @brief 2D 裂缝静态物理属性管理子模块
 * @details 负责 2D-EDFM 框架下离散裂缝网格 (Fracture) 物理属性的初始化与内存分配。
 * 提供对裂缝孔隙度、渗透率、孔径及岩石热物性的统一安全管理。
 */

#pragma once

#include <vector>
#include <iostream>
#include <omp.h>

#include "MeshManager.h"
#include "2D_FieldManager.h"
#include "PropertiesSummary.h"
#include "SolverContrlStrName_op.h"

 /**
  * @class FractureProperties_2D
  * @brief 2D 裂缝固体属性管理器
  */
class FractureProperties_2D
{
public:
    /**
     * @brief 构造函数
     * @param meshMgr 2D网格管理器常量引用
     * @param fieldMgr 2D场管理器引用
     */
    FractureProperties_2D(const MeshManager& meshMgr, FieldManager_2D& fieldMgr);

    ~FractureProperties_2D() = default;

    /**
     * @brief 设置裂缝全局基础物理属性
     * @param params 包含渗透率、孔隙度、孔径等参数的结构体
     */
    void setGlobalProperties(const SolidProperties_Frac& params);

    /**
     * @brief 执行裂缝物理属性初始化
     * @details 遍历所有离散裂缝段，计算安全偏移量，并将属性写入 FieldManager。
     * 若几何网格对象自身带有孔径属性，则优先使用几何网格对象的孔径。
     */
    void InitializeFractureProperties();

private:
    const MeshManager& meshMgr_; ///< 拓扑网格管理器引用
    FieldManager_2D& fieldMgr_;  ///< 物理场数据管理器引用

    SolidProperties_Frac params_; ///< 裂缝全局物理属性
};