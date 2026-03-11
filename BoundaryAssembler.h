/**
 * @file BoundaryAssembler.h
 * @brief Boundary and well assembly facade
 */
#ifndef BOUNDARY_ASSEMBLER_H
#define BOUNDARY_ASSEMBLER_H

#include <vector>
#include <string>
#include <array>
#include "MeshManager.h"
#include "3D_MeshManager.h"
#include "BoundaryConditionManager.h"
#include "2D_FieldManager.h"
#include "3D_FieldManager.h"
#include "Well_WellControlTypes.h"

struct BoundaryAssemblyStats {
    int matrixBCCount = 0;
    int fractureBCCount = 0;
    double sumResidual = 0.0;
    double sumJacobianDiag = 0.0;
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

    static BoundaryAssemblyStats Assemble_3D(
        MeshManager_3D& mgr,
        const BoundarySetting::BoundaryConditionManager& bcMgr,
        int dofOffset,
        FieldManager_3D& fm,
        const std::string& fieldName,
        std::vector<double>& residual,
        std::vector<double>& jacobianDiag
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
        std::vector<std::array<double, 3>>& jacobianFull
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
        std::vector<std::array<double, 3>>& jacobianFull
    );
};

#endif // BOUNDARY_ASSEMBLER_H
