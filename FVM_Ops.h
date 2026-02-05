/**
 * @file FVM_Ops.h
 * @brief 有限体积法 (FVM) 基础通用算子库
 * @details 包含插值、梯度计算、迎风格式流向判定等核心离散算子。
 * 设计为 Header-only 模板库，以支持标量、向量及张量场。
 * @version 1.0.0
 */

#pragma once

#include <cmath>		// std::abs, std::sqrt
#include <algorithm>	// std::max, std::min
#include <limits>		// std::numeric_limits
#include <iostream>		// std::cerr

#include "mesh.h"           // 包含 Mesh, Face, Cell, Node
#include "VolField.h"       // 包含 VolField<T>
#include "UserDefineVarType.h"   // 包含 Vector, Tensor

namespace FVM_Ops
{
    // =========================================================
    // 常量定义
    // =========================================================
    static constexpr double kEpsilon = 1e-12; ///! 防止除零的小量

    // =========================================================
    // 第一组：插值类算子 (Interpolation Operators)
    // =========================================================

    /**
     * @brief 算术平均插值 (Arithmetic Interpolation)
     * @details 基于几何权重将单元中心值插值到面中心：
     * phi_f = w * phi_O + (1 - w) * phi_N
     * @tparam T 物理量类型 (double, Vector, etc.)
     * @param mesh 网格对象
     * @param field 待插值的体心场
     * @param faceIdx 面索引 在输入前需要确保索引正确
     * @return T 插值后的面心值
     */
    template<typename T>
    inline T Op_Interp_Arithmetic(const Mesh& mesh, const VolField<T>& field, int faceIdx)
    {
		const Face& face = mesh.getFaces()[faceIdx]; 
        int owner = face.ownerCell;
		int neighbor = face.neighborCell;
		size_t idxO = mesh.getCellId2Index().at(owner);
		const T& phiO = field[idxO];
		if (face.isBoundary())
		{
			return phiO; // 边界面直接取 owner 值
		}
		else
		{
			size_t idxN = mesh.getCellId2Index().at(neighbor);
			const T& phiN = field[idxN];
			double w = face.f_linearInterpolationCoef;
			return w * phiN + (1.0 - w) * phiO;
		}

    }

}
