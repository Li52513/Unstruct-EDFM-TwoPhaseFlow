/**
 * @file BoundaryAssembler.h
 * @brief Boundary and well assembly facade
 */
#ifndef BOUNDARY_ASSEMBLER_H
#define BOUNDARY_ASSEMBLER_H

#include <vector>
#include <string>
#include <array>
#include <map>
#include "MeshManager.h"
#include "3D_MeshManager.h"
#include "BoundaryConditionManager.h"
#include "2D_FieldManager.h"
#include "3D_FieldManager.h"
#include "Well_WellControlTypes.h"
#include "CapRelPerm_HD.h"

 /**
  * @struct BoundaryTagStats
  * @brief 按边界 Tag 与类型聚合的详细统计
  */
struct BoundaryTagStats {
    int faces = 0;       ///< 识别到的物理面数量
    int applied = 0;     ///< 成功应用到 Jacobian 的方程行数
    int skipped = 0;     ///< 判定为无效或强置零跳过的行数
    double sumR = 0.0;   ///< 残差绝对值累加
    double sumDiag = 0.0;///< 对角线绝对值累加
};

struct BoundaryAssemblyStats {
    int matrixBCCount = 0;
    int fractureBCCount = 0;
    double sumResidual = 0.0;
    double sumJacobianDiag = 0.0;

    int visitedEqRows = 0;  ///< 访问的总方程行数
    int nonzeroEqRows = 0;  ///< 产生非零贡献的方程行数
    int zeroEqRows = 0;     ///< 因 Dirichlet 或无流边界导致贡献为零的行数
    int invalidEqRows = 0;  ///< 索引越界等非法行数
    std::map<std::string, BoundaryTagStats> perTagType; ///< 以 "Tag_Type" 为键的细粒度统计
};

class BoundaryAssembler {
public:
        static BoundaryAssemblyStats Assemble_2D(
        MeshManager& mgr,
        const BoundarySetting::BoundaryConditionManager& bcMgr,
        int dofOffset,
        FieldManager_2D& fm,
        const std::string& fieldName,
        std::vector<double>& residual,
        std::vector<double>& jacobianDiag
    );

    static BoundaryAssemblyStats Assemble_2D_FullJac(
        MeshManager& mgr,
        const BoundarySetting::BoundaryConditionManager& bcMgr,
        int dofOffset,
        FieldManager_2D& fm,
        const std::string& fieldName,
        std::vector<double>& residual,
        std::vector<std::array<double, 3>>& jacobianFull,
        const BoundarySetting::BoundaryConditionManager* coupledPressureBC = nullptr,
        bool single_phase_use_co2 = false,
        const CapRelPerm::VGParams& vg = CapRelPerm::VGParams(),
        const CapRelPerm::RelPermParams& rp = CapRelPerm::RelPermParams()
    );

    static BoundaryAssemblyStats Assemble_3D(
        MeshManager_3D& mgr,
        const BoundarySetting::BoundaryConditionManager& bcMgr,
        int dofOffset,
        FieldManager_3D& fm,
        const std::string& fieldName,
        std::vector<double>& residual,
        std::vector<double>& jacobianDiag
    );

    static BoundaryAssemblyStats Assemble_3D_FullJac(
        MeshManager_3D& mgr,
        const BoundarySetting::BoundaryConditionManager& bcMgr,
        int dofOffset,
        FieldManager_3D& fm,
        const std::string& fieldName,
        std::vector<double>& residual,
        std::vector<std::array<double, 3>>& jacobianFull,
        const BoundarySetting::BoundaryConditionManager* coupledPressureBC = nullptr,
        bool single_phase_use_co2 = false,
        const CapRelPerm::VGParams& vg = CapRelPerm::VGParams(),
        const CapRelPerm::RelPermParams& rp = CapRelPerm::RelPermParams()
    );

    static BoundaryAssemblyStats Assemble_Wells_2D(
        MeshManager& mgr,
        FieldManager_2D& fm,
        const std::vector<WellScheduleStep>& active_steps,
        int dofOffset_P,
        int dofOffset_W,
        int dofOffset_G,
        int dofOffset_E,
        std::vector<double>& residual,
        std::vector<double>& jacobianDiag
    );

    static BoundaryAssemblyStats Assemble_Wells_2D_FullJac(
        MeshManager& mgr,
        FieldManager_2D& fm,
        const std::vector<WellScheduleStep>& active_steps,
        int dofOffset_P,
        int dofOffset_W,
        int dofOffset_G,
        int dofOffset_E,
        std::vector<double>& residual,
        std::vector<std::array<double, 3>>& jacobianFull,
        bool single_phase_use_co2 = false,
        const CapRelPerm::VGParams& vg = CapRelPerm::VGParams(),
        const CapRelPerm::RelPermParams& rp = CapRelPerm::RelPermParams()
    );

    static BoundaryAssemblyStats Assemble_Wells_3D(
        MeshManager_3D& mgr,
        FieldManager_3D& fm,
        const std::vector<WellScheduleStep>& active_steps,
        int dofOffset_P,
        int dofOffset_W,
        int dofOffset_G,
        int dofOffset_E,
        std::vector<double>& residual,
        std::vector<double>& jacobianDiag
    );

    static BoundaryAssemblyStats Assemble_Wells_3D_FullJac(
        MeshManager_3D& mgr,
        FieldManager_3D& fm,
        const std::vector<WellScheduleStep>& active_steps,
        int dofOffset_P,
        int dofOffset_W,
        int dofOffset_G,
        int dofOffset_E,
        std::vector<double>& residual,
        std::vector<std::array<double, 3>>& jacobianFull,
        bool single_phase_use_co2 = false,
        const CapRelPerm::VGParams& vg = CapRelPerm::VGParams(),
        const CapRelPerm::RelPermParams& rp = CapRelPerm::RelPermParams()
    );
};

#endif // BOUNDARY_ASSEMBLER_H
