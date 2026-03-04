/**
 * @file BoundaryAssembler.h
 * @brief 边界条件与漏失项装配门面 (Day 3 Production Entry)
 * @details 负责遍历基岩网格的物理边界，调用 AD 算子核计算通量/漏失残差，并将其累加至全局残差向量和 Jacobian 对角线对应的方程块中。
 */
#ifndef BOUNDARY_ASSEMBLER_H
#define BOUNDARY_ASSEMBLER_H

#include <vector>
#include <string>
#include "MeshManager.h"
#include "3D_MeshManager.h"
#include "BoundaryConditionManager.h"
#include "2D_FieldManager.h"
#include "3D_FieldManager.h"

struct BoundaryAssemblyStats {
    int matrixBCCount = 0;
    int fractureBCCount = 0;
    double sumResidual = 0.0;
    double sumJacobianDiag = 0.0;
};

class BoundaryAssembler {
public:
    /**
     * @brief 2D 边界与漏失主路径装配入口 (支持指定 DOF 偏移)
     * @param mgr 2D 网格管理器
     * @param bcMgr 边界条件管理器
     * @param dofOffset 当前场对应的方程行偏移 (如 Pressure=0, Temperature=1)
     * @param fm 场管理器
     * @param fieldName 需要读取的主变量场名
     * @param residual 全局残差向量
     * @param jacobianDiag 全局 Jacobian 主对角线向量
     */
    static BoundaryAssemblyStats Assemble_2D(
        MeshManager& mgr,
        const BoundarySetting::BoundaryConditionManager& bcMgr,
        int dofOffset,
        FieldManager_2D& fm,
        const std::string& fieldName,
        std::vector<double>& residual,
        std::vector<double>& jacobianDiag
    );

    /**
     * @brief 3D 边界与漏失主路径装配入口 (支持指定 DOF 偏移)
     */
    static BoundaryAssemblyStats Assemble_3D(
        MeshManager_3D& mgr,
        const BoundarySetting::BoundaryConditionManager& bcMgr,
        int dofOffset,
        FieldManager_3D& fm,
        const std::string& fieldName,
        std::vector<double>& residual,
        std::vector<double>& jacobianDiag
    );
};

#endif // BOUNDARY_ASSEMBLER_H