#pragma once

#include "UserDefineVarType.h" // 获取 Vector 定义

/**
 * @struct LinearInitParams
 * @brief 线性分布初始化参数包 (纯数学结构，2D/3D通用)
 * @details 用于描述物理场在空间中的线性分布:
 * Val(x) = refVal + grad_x * (x - x_ref) + grad_y * (y - y_ref) + grad_z * (z - z_ref)
 */
struct LinearInitParams
{
    double refVal = 0.0;        ///< 参考点的值
    double x_ref = 0.0;         ///< 参考点 X 坐标
    double y_ref = 0.0;         ///< 参考点 Y 坐标
    double z_ref = 0.0;         ///< 参考点 Z 坐标

    double grad_x = 0.0;        ///< X 方向梯度
    double grad_y = 0.0;        ///< Y 方向梯度
    double grad_z = 0.0;        ///< Z 方向梯度 (2D 模拟中保持为 0 即可)

    LinearInitParams() = default;
    explicit LinearInitParams(double val) : refVal(val) {}
    LinearInitParams(double val, double zRef, double gradZ)
        : refVal(val), z_ref(zRef), grad_z(gradZ) {
    }
};

/**
 * @brief 根据空间坐标和参数包计算线性分布值
 * @param pos 空间坐标点 (Vector)
 * @param p 线性参数包
 * @return double 计算得到的物理量初值
 */
inline double ComputeLinearValue(const Vector& pos, const LinearInitParams& p)
{
    return p.refVal +
        p.grad_x * (pos.m_x - p.x_ref) +
        p.grad_y * (pos.m_y - p.y_ref) +
        p.grad_z * (pos.m_z - p.z_ref);
}