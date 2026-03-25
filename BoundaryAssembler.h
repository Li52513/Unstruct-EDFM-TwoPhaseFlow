/**
 * @file BoundaryAssembler.h
 * @brief Boundary and well assembly facade
 */
#ifndef BOUNDARY_ASSEMBLER_H
#define BOUNDARY_ASSEMBLER_H

#include <array>
#include <map>
#include <string>
#include <vector>

#include "MeshManager.h"
#include "3D_MeshManager.h"
#include "BoundaryConditionManager.h"
#include "2D_FieldManager.h"
#include "3D_FieldManager.h"
#include "Well_WellControlTypes.h"
#include "CapRelPerm_HD.h"

struct BoundaryTagStats {
    int faces = 0;
    int applied = 0;
    int skipped = 0;
    double sumR = 0.0;
    double sumDiag = 0.0;
};

struct BoundaryAssemblyStats {
    int matrixBCCount = 0;
    int fractureBCCount = 0;
    double sumResidual = 0.0;
    double sumJacobianDiag = 0.0;
    int inferredTagFail = 0;
    int satCompositionUsed = 0;
    int satDriveIgnored = 0;

    int visitedEqRows = 0;
    int nonzeroEqRows = 0;
    int zeroEqRows = 0;
    int invalidEqRows = 0;
    std::map<std::string, BoundaryTagStats> perTagType;
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
        const BoundarySetting::BoundaryConditionManager* coupledSaturationBC = nullptr,
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
        const BoundarySetting::BoundaryConditionManager* coupledSaturationBC = nullptr,
        bool single_phase_use_co2 = false,
        const CapRelPerm::VGParams& vg = CapRelPerm::VGParams(),
        const CapRelPerm::RelPermParams& rp = CapRelPerm::RelPermParams()
    );

    static BoundaryAssemblyStats Assemble_2D_CoupledN3_FullJac(
        MeshManager& mgr,
        FieldManager_2D& fm,
        const BoundarySetting::BoundaryConditionManager* pressureBC,
        const BoundarySetting::BoundaryConditionManager* saturationBC,
        const BoundarySetting::BoundaryConditionManager* temperatureBC,
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

    static BoundaryAssemblyStats Assemble_3D_CoupledN3_FullJac(
        MeshManager_3D& mgr,
        FieldManager_3D& fm,
        const BoundarySetting::BoundaryConditionManager* pressureBC,
        const BoundarySetting::BoundaryConditionManager* saturationBC,
        const BoundarySetting::BoundaryConditionManager* temperatureBC,
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
