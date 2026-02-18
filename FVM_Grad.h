/**
 * @file FVM_Grad.h
 * @brief 有限体积法梯度计算算子 (Gradient Operator) - Final Production Ready
 * @details
 * 提供基于非结构化网格的标量场梯度计算功能。
 * 核心特性：
 * 1. 鲁棒的 Least-Squares 算法：SVD 伪逆处理，兼容 3D/2D/1D 拓扑。
 * 2. 物理一致性：强制执行切空间投影 (Tangential Projection)。
 * 3. 工程健壮性：包含尺寸一致性检查、奇异矩阵防御及边界镇定机制。
 * 4. FieldManager 集成：支持场数据的自动注册与管理。
 */

#pragma once

#include <vector>
#include <string>
#include <memory>
#include <iostream>

 // 引入数学与网格库
#include <Eigen/Dense>
#include "UserDefineVarType.h"
#include "Mesh.h"
#include "VolField.h"
#include "BoundaryConditionManager.h"

// 前向声明
class FractureNetwork;      // 2D-EDFM
class FractureNetwork_2D;   // 3D-EDFM
class FieldManager_3D;
class FieldManager_2D;

class FVM_Grad
{
public:
    enum class Method
    {
        GreenGauss,     ///< 格林-高斯法 (基岩推荐)
        LeastSquares    ///< 最小二乘法 (裂缝必选，扭曲基岩推荐)
    };

    /**
     * @brief 构造函数
     * @param mesh 基岩网格
     * @param fracNet3D 3D-EDFM 网络 (可选)
     * @param fracNet2D 2D-EDFM 网络 (可选)
     * @param bcMgr 边界条件管理器 (可选)
     */
    FVM_Grad(const Mesh& mesh,
        const FractureNetwork_2D* fracNet3D = nullptr,
        const FractureNetwork* fracNet2D = nullptr,
        const BoundarySetting::BoundaryConditionManager* bcMgr = nullptr);

    ~FVM_Grad() = default;

    // =========================================================
    // 1. Matrix (基岩) 接口
    // =========================================================
    void precomputeLS();

    std::shared_ptr<volVectorField> compute(const volScalarField& scalarField, Method method = Method::LeastSquares);

    std::shared_ptr<volVectorField> compute(const std::string& fieldName, FieldManager_3D& fm, Method method = Method::LeastSquares);
    std::shared_ptr<volVectorField> compute(const std::string& fieldName, FieldManager_2D& fm, Method method = Method::LeastSquares);

    // =========================================================
    // 2. Fracture (裂缝) 接口
    // =========================================================
    /**
     * @brief 预计算裂缝 LS 矩阵
     * @details 包含物理投影、边界镇定及奇异值截断保护
     */
    void precomputeLS_Fracture();

    std::shared_ptr<volVectorField> compute(const volScalarField& fracField,
        const FractureNetwork_2D& frNet,
        Method method = Method::LeastSquares);

    std::shared_ptr<volVectorField> compute(const volScalarField& fracField,
        const FractureNetwork& frNet,
        Method method = Method::LeastSquares);

    std::shared_ptr<volVectorField> computeFractureGrad(const std::string& fieldName,
        FieldManager_3D& fm,
        Method method = Method::LeastSquares);

    std::shared_ptr<volVectorField> computeFractureGrad(const std::string& fieldName,
        FieldManager_2D& fm,
        Method method = Method::LeastSquares);

private:
    // =========================================================
    // 成员变量
    // =========================================================
    const Mesh& mesh_;
    const FractureNetwork_2D* frNet3D_;
    const FractureNetwork* frNet2D_;
    const BoundarySetting::BoundaryConditionManager* bcMgr_;

    // Matrix Cache
    std::vector<Eigen::Matrix3d> ls_pinv_matrices_;
    bool ls_precomputed_ = false;

    // Fracture Cache
    // Index: Local Fracture Index (Global Solver Index - Offset)
    std::vector<Eigen::Matrix3d> ls_frac_pinv_matrices_;
    bool ls_frac_precomputed_ = false;

    // =========================================================
    // 内部实现
    // =========================================================
    Eigen::Vector3d toEigen(const Vector& v) const;
    Vector toUser(const Eigen::Vector3d& v) const;

    Vector _computeGrad_GG_Cell(int cellIndex, const volScalarField& field) const;
    Vector _computeGrad_LS_Cell(int cellIndex, const volScalarField& field) const;
    double _getBoundaryFaceValue_GG(int faceIndex, double cellValue) const;

    Vector _computeGrad_LS_FractureCell(int localFracIdx, const volScalarField& field) const;
};