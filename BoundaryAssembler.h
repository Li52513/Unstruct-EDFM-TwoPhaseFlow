/**
 * @file BoundaryAssembler.h
 * @brief Boundary and well assembly facade
 */
#ifndef BOUNDARY_ASSEMBLER_H
#define BOUNDARY_ASSEMBLER_H

#include <array>
#include <map>
#include <string>
#include <unordered_map>
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

struct FluidConstantProperties {
    double rho = 1000.0;
    double mu = 1.0e-3;
    double cp = 4200.0;
    double cv = 4182.0;
    double k = 0.6;
};

struct FluidPropertyEvalConfig {
    bool enable_single_phase_constant = false;
    bool single_phase_is_co2 = false;
    bool single_phase_no_convection = false;

    bool enable_two_phase_constant = false;

    FluidConstantProperties water = FluidConstantProperties{};
    FluidConstantProperties gas = FluidConstantProperties{ 700.0, 5.0e-5, 1200.0, 900.0, 0.08 };
};

struct WellEquationLinearization {
    bool valid = false;
    int eq_index = -1;         ///< global equation index in reservoir system
    int eq_dof = -1;           ///< equation DOF within cell block
    double q = 0.0;            ///< source contribution (outflow-positive convention)
    std::array<double, 3> dq_dcell{ {0.0, 0.0, 0.0} }; ///< d/d(P,Sw,T) local cell state
    double dq_dPbh = 0.0;      ///< derivative w.r.t well BHP
};

struct WellCompletionLinearization {
    size_t step_index = 0;     ///< index in active_steps passed to assembler
    std::string well_name;
    int completion_solver_index = -1;
    WellTargetDomain domain = WellTargetDomain::Matrix;
    double pbh_eval = 0.0;

    WellEquationLinearization water;
    WellEquationLinearization gas;
    WellEquationLinearization energy;

    double q_total = 0.0;
    std::array<double, 3> dqtotal_dcell{ {0.0, 0.0, 0.0} };
    double dqtotal_dPbh = 0.0;
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
        const FluidPropertyEvalConfig& fluid_cfg = FluidPropertyEvalConfig(),
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
        const FluidPropertyEvalConfig& fluid_cfg = FluidPropertyEvalConfig(),
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
        const FluidPropertyEvalConfig& fluid_cfg = FluidPropertyEvalConfig(),
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
        const FluidPropertyEvalConfig& fluid_cfg = FluidPropertyEvalConfig(),
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
        const FluidPropertyEvalConfig& fluid_cfg = FluidPropertyEvalConfig(),
        const CapRelPerm::VGParams& vg = CapRelPerm::VGParams(),
        const CapRelPerm::RelPermParams& rp = CapRelPerm::RelPermParams(),
        std::vector<WellCompletionLinearization>* completionLinearizations = nullptr,
        const std::unordered_map<std::string, double>* wellBhpByName = nullptr
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
        const FluidPropertyEvalConfig& fluid_cfg = FluidPropertyEvalConfig(),
        const CapRelPerm::VGParams& vg = CapRelPerm::VGParams(),
        const CapRelPerm::RelPermParams& rp = CapRelPerm::RelPermParams(),
        std::vector<WellCompletionLinearization>* completionLinearizations = nullptr,
        const std::unordered_map<std::string, double>* wellBhpByName = nullptr
    );
};

#endif // BOUNDARY_ASSEMBLER_H
