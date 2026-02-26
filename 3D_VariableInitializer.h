#pragma once

#include <vector>
#include <string>
#include <memory>
#include <omp.h>
#include <iostream>

#include "3D_MeshManager.h"
#include "3D_FieldManager.h"
#include "SolverContrlStrName_op.h"
#include "CapRelPerm_HD.h" // [新增] 引入相渗/毛管力模型

/**
 * @struct LinearInitParams
 * @brief 线性分布初始化参数包
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
    double grad_z = 0.0;        ///< Z 方向梯度 (如 rho * g)

    LinearInitParams() = default;
    explicit LinearInitParams(double val) : refVal(val) {}
    LinearInitParams(double val, double zRef, double gradZ)
        : refVal(val), z_ref(zRef), grad_z(gradZ) {
    }
};

// =========================================================
// 辅助函数: 线性分布计算
// =========================================================
inline double ComputeLinearValue(const Vector& pos, const LinearInitParams& p)
{
    return p.refVal +
        p.grad_x * (pos.m_x - p.x_ref) +
        p.grad_y * (pos.m_y - p.y_ref) +
        p.grad_z * (pos.m_z - p.z_ref);
}

/**
 * @class VariableInitializer_3D
 * @brief 3D 主变量场构建与初始化器
 */
class VariableInitializer_3D
{
public:
    VariableInitializer_3D(const MeshManager_3D& meshMgr, FieldManager_3D& fieldMgr);

    // =========================================================
    // 单相流初始化 (Single Phase)
    // =========================================================

    /**
     * @brief 初始化单相流主变量
     */
    void InitSinglePhaseState(const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig,
        const LinearInitParams& pInit,
        const LinearInitParams& tInit);

    // =========================================================
    // IMPES 两相流初始化 (Two Phase)
    // =========================================================

    /**
     * @brief 初始化 IMPES 两相流主变量及辅助变量 (Pc, Kr)
     * @param pConfig 压力方程配置
     * @param tConfig 温度方程配置
     * @param sConfig 饱和度方程配置
     * @param pInit 压力初始化参数
     * @param tInit 温度初始化参数
     * @param sInit 饱和度初始化参数
     * @param vgParams VG模型参数 (用于计算初始 Pc)
     * @param rpParams 相渗参数 (用于计算初始 Kr)
     */
    void InitIMPESState(const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig,
        const PhysicalProperties_string_op::SaturationEquation_String& sConfig,
        const LinearInitParams& pInit,
        const LinearInitParams& tInit,
        const LinearInitParams& sInit,
        const CapRelPerm::VGParams& vgParams,        
        const CapRelPerm::RelPermParams& rpParams);

    // =========================================================
    // FIM 全隐式两相流/单相流初始化 (自动梯度播种)
    // =========================================================
    /**
     * @brief 初始化 FIM 框架下的主变量并完成雅可比矩阵的梯度播种
     * @tparam N 独立自变量的数量 (单相N=2, 两相N=3)
     * @param pIdx 压力 P 在系统中的局部自由度索引 (如 0)
     * @param sIdx 饱和度 Sw 的自由度索引 (如 1, 单相流传入 -1 忽略播种)
     * @param tIdx 温度 T 的自由度索引 (如 2)
     */
    template<int N>
    void InitFIMState(const PhysicalProperties_string_op::PressureEquation_String& pConfig,
        const PhysicalProperties_string_op::TemperatureEquation_String& tConfig,
        const PhysicalProperties_string_op::SaturationEquation_String& sConfig,
        const LinearInitParams& pInit,
        const LinearInitParams& tInit,
        const LinearInitParams& sInit,
        int pIdx = 0, int sIdx = 1, int tIdx = 2)
    {
        std::cout << "[VarInit_3D] Initializing FIM State (N=" << N << ") with AD Seeding..." << std::endl;

        // 1. 初始化压力场 (带导数播种)
        CreateAndInitMatrixADScalar<N>(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit, pIdx);
        CreateAndInitFractureADScalar<N>(pConfig.pressure_field, pConfig.pressure_old_field, pConfig.pressure_prev_field, pInit, pIdx);

        // 2. 初始化温度场 (带导数播种)
        CreateAndInitMatrixADScalar<N>(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit, tIdx);
        CreateAndInitFractureADScalar<N>(tConfig.temperatue_field, tConfig.temperatue_old_field, tConfig.temperatue_prev_field, tInit, tIdx);

        // 3. 初始化饱和度场 (仅在 sIdx >= 0 时播种，用于两相流)
        if (sIdx >= 0) {
            CreateAndInitMatrixADScalar<N>(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit, sIdx);
            CreateAndInitFractureADScalar<N>(sConfig.saturation, sConfig.saturation_old, sConfig.saturation_prev, sInit, sIdx);
        }

        std::cout << "[VarInit_3D] FIM State Seeding Completed." << std::endl;
    }

private:
    // =========================================================
    // 内部通用实现
    // =========================================================

    /**
     * @brief 通用标量场初始化 (Matrix) - 强制创建
     */
    void CreateAndInitMatrixScalar(const std::string& tagName,
        const std::string& oldTagName,
        const std::string& prevTagName,
        const LinearInitParams& params);

    /**
     * @brief 通用标量场初始化 (Fracture) - 强制创建
     */
    void CreateAndInitFractureScalar(const std::string& tagName,
        const std::string& oldTagName,
        const std::string& prevTagName,
        const LinearInitParams& params);

    /**
     * @brief 计算并填充两相流衍生变量 (Pc, Kr)
     * @details 基于已初始化的 Sw 场，计算 Pc, krw, krg 并填入场管理器
     */
    void CalculateDerivedTwoPhaseFields(const PhysicalProperties_string_op::SaturationEquation_String& sConfig,
        const CapRelPerm::VGParams& vg,
        const CapRelPerm::RelPermParams& rp);

    // =========================================================
    // AD 场分配与播种底层通用实现
    // =========================================================
    template<int N>
    void CreateAndInitMatrixADScalar(const std::string& tagName, const std::string& oldTagName, const std::string& prevTagName, const LinearInitParams& params, int seedIdx)
    {
        auto& field = fieldMgr_.createMatrixADScalar<N>(tagName)->data;
        auto& fieldOld = fieldMgr_.createMatrixADScalar<N>(oldTagName)->data;
        auto& fieldPrev = fieldMgr_.createMatrixADScalar<N>(prevTagName)->data;

        const auto& cells = meshMgr_.mesh().getCells();
        int nCells = static_cast<int>(cells.size());

#pragma omp parallel for
        for (int i = 0; i < nCells; ++i)
        {
            const Vector& center = cells[i].center;
            double val = ComputeLinearValue(center, params);
            ADVar<N> adVal(val, seedIdx);
            field[i] = adVal;
            fieldOld[i] = adVal;
            fieldPrev[i] = adVal;
        }
    }

    template<int N>
    void CreateAndInitFractureADScalar(const std::string& tagName, const std::string& oldTagName, const std::string& prevTagName, const LinearInitParams& params, int seedIdx)
    {
        auto& field = fieldMgr_.createFractureADScalar<N>(tagName)->data;
        auto& fieldOld = fieldMgr_.createFractureADScalar<N>(oldTagName)->data;
        auto& fieldPrev = fieldMgr_.createFractureADScalar<N>(prevTagName)->data;

        const auto& frNet = meshMgr_.fracture_network();
        const auto& fractures = frNet.getFractures();
        const auto& indexer = frNet.fracElemIndex;

        if (indexer.offset.size() < fractures.size()) return;
        size_t totalFracCells = field.size();

        for (size_t i = 0; i < fractures.size(); ++i)
        {
            const auto& frac = fractures[i];
            size_t startOffset = indexer.offset[i];
            size_t numElems = frac.fracCells.size();

#pragma omp parallel for
            for (int localIdx = 0; localIdx < static_cast<int>(numElems); ++localIdx)
            {
                size_t globalIdx = startOffset + localIdx;
                const Vector& center = frac.fracCells[localIdx].centroid;
                double val = ComputeLinearValue(center, params);

                ADVar<N> adVal(val, seedIdx);
                if (globalIdx < totalFracCells) {
                    field[globalIdx] = adVal;
                    fieldOld[globalIdx] = adVal;
                    fieldPrev[globalIdx] = adVal;
                }
            }
        }
    }

    const MeshManager_3D& meshMgr_;
    FieldManager_3D& fieldMgr_;
};