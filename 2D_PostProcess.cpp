/**
 * @file 2D_PostProcess.cpp
 * @brief 2D-EDFM 锟斤拷锟斤拷锟斤拷锟接伙拷模锟斤拷 (Tecplot Exporter) 实锟斤拷
 * @details 锟斤拷锟斤拷 2D 锟斤拷锟斤拷锟斤拷锟?1D 锟窖凤拷锟斤拷募锟斤拷锟斤拷锟斤拷思锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟轿?Tecplot BLOCK 锟斤拷式锟斤拷
 * 支锟街伙拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟?锟侥憋拷锟轿ｏ拷锟斤拷锟斤拷应锟斤拷洌拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟窖凤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟解。
 */

#include "2D_PostProcess.h"
#include "SolverContrlStrName_op.h"
#include "AD_FluidEvaluator.h"
#include "CapRelPerm_HD_AD.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <array>
#include <limits>
#include <set>
#include <unordered_map>
#include <vector>
#include <string>

 // ==============================================================================
 // 锟斤拷锟届函锟斤拷
 // ==============================================================================
PostProcess_2D::PostProcess_2D(const MeshManager& meshMgr,
                               const FieldManager_2D& fieldMgr,
                               const VTKBoundaryVisualizationContext* bcVizCtx)
    : meshMgr_(meshMgr), fieldMgr_(fieldMgr), bcVizCtx_(bcVizCtx)
{
}

namespace {
    inline double SafeFieldValue(const std::shared_ptr<volScalarField>& field, int idx, double fallback) {
        if (!field) return fallback;
        if (idx < 0 || idx >= static_cast<int>(field->data.size())) return fallback;
        const double v = field->data[idx];
        return std::isfinite(v) ? v : fallback;
    }

    struct VtkBCFieldCache2D {
        std::shared_ptr<volScalarField> kxx;
        std::shared_ptr<volScalarField> kyy;
        std::shared_ptr<volScalarField> kzz;
        std::shared_ptr<volScalarField> lambda_eff;
        std::shared_ptr<volScalarField> lambda_r;
        std::shared_ptr<volScalarField> rho_w;
        std::shared_ptr<volScalarField> mu_w;
        std::shared_ptr<volScalarField> lamw_mob;
        std::shared_ptr<volScalarField> rho_g;
        std::shared_ptr<volScalarField> mu_g;
        std::shared_ptr<volScalarField> lamg_mob;
    };

    VtkBCFieldCache2D BuildVtkBCFieldCache2D(const FieldManager_2D& fieldMgr) {
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::EffectiveProps eff;
        const PhysicalProperties_string_op::Water water;
        const PhysicalProperties_string_op::CO2 gas;

        VtkBCFieldCache2D cache;
        cache.kxx = fieldMgr.getMatrixScalar(rock.k_xx_tag);
        cache.kyy = fieldMgr.getMatrixScalar(rock.k_yy_tag);
        cache.kzz = fieldMgr.getMatrixScalar(rock.k_zz_tag);
        cache.lambda_eff = fieldMgr.getMatrixScalar(eff.lambda_eff_tag);
        cache.lambda_r = fieldMgr.getMatrixScalar(rock.lambda_tag);
        cache.rho_w = fieldMgr.getMatrixScalar(water.rho_tag);
        cache.mu_w = fieldMgr.getMatrixScalar(water.mu_tag);
        cache.lamw_mob = fieldMgr.getMatrixScalar(water.lambda_w_tag);
        cache.rho_g = fieldMgr.getMatrixScalar(gas.rho_tag);
        cache.mu_g = fieldMgr.getMatrixScalar(gas.mu_tag);
        cache.lamg_mob = fieldMgr.getMatrixScalar(gas.lambda_g_tag);
        return cache;
    }

    double ComputeKn2D(const VtkBCFieldCache2D& cache, int owner, const Vector& normal, double floor_v) {
        const double kxx = std::max(SafeFieldValue(cache.kxx, owner, 1.0e-13), 0.0);
        const double kyy = std::max(SafeFieldValue(cache.kyy, owner, 1.0e-13), 0.0);
        const double kzz = std::max(SafeFieldValue(cache.kzz, owner, 1.0e-13), 0.0);
        const double kn_raw =
            normal.m_x * normal.m_x * kxx +
            normal.m_y * normal.m_y * kyy +
            normal.m_z * normal.m_z * kzz;
        if (kn_raw > floor_v) return kn_raw;
        return std::max({ kxx, kyy, kzz, floor_v });
    }

    double ComputeWaterMassMobility2D(const VtkBCFieldCache2D& cache, int owner, double floor_v) {
        const double rho_w = std::max(SafeFieldValue(cache.rho_w, owner, 1000.0), floor_v);
        const double lamw_mob = SafeFieldValue(cache.lamw_mob, owner, std::numeric_limits<double>::quiet_NaN());
        if (std::isfinite(lamw_mob) && lamw_mob > 0.0) {
            return std::max(rho_w * lamw_mob, floor_v);
        }
        const double mu_w = std::max(SafeFieldValue(cache.mu_w, owner, 1.0e-3), floor_v);
        return std::max(rho_w / mu_w, floor_v);
    }

    double ComputeGasMassMobility2D(const VtkBCFieldCache2D& cache, int owner, double floor_v) {
        const double rho_g = std::max(SafeFieldValue(cache.rho_g, owner, 1.0), floor_v);
        const double lamg_mob = SafeFieldValue(cache.lamg_mob, owner, std::numeric_limits<double>::quiet_NaN());
        if (std::isfinite(lamg_mob) && lamg_mob > 0.0) {
            return std::max(rho_g * lamg_mob, floor_v);
        }
        const double mu_g = std::max(SafeFieldValue(cache.mu_g, owner, 1.0e-5), floor_v);
        return std::max(rho_g / mu_g, floor_v);
    }

    double ComputeTransportCoeff2D(const VTKBCVariableBinding& binding,
                                   const VtkBCFieldCache2D& cache,
                                   int owner,
                                   const Vector& normal,
                                   double floor_v) {
        const double kn = ComputeKn2D(cache, owner, normal, floor_v);
        if (binding.transport_kind == VTKBCTransportKind::Temperature) {
            const double lambda_eff = SafeFieldValue(cache.lambda_eff, owner, std::numeric_limits<double>::quiet_NaN());
            if (std::isfinite(lambda_eff) && lambda_eff > floor_v) return lambda_eff;
            return std::max(SafeFieldValue(cache.lambda_r, owner, 2.0), floor_v);
        }

        if (binding.transport_kind == VTKBCTransportKind::Saturation) {
            return std::max(kn * ComputeGasMassMobility2D(cache, owner, floor_v), floor_v);
        }

        const bool use_gas_pressure = (binding.field_name == "p_g");
        const double mass_mob = use_gas_pressure
            ? ComputeGasMassMobility2D(cache, owner, floor_v)
            : ComputeWaterMassMobility2D(cache, owner, floor_v);
        return std::max(kn * mass_mob, floor_v);
    }

    void ReconstructBCPointField2D(const MeshManager& meshMgr,
                                   const FieldManager_2D& fieldMgr,
                                   const std::shared_ptr<volScalarField>& matCellField,
                                   const VTKBoundaryVisualizationContext& ctx,
                                   const VTKBCVariableBinding& binding,
                                   const std::unordered_map<int, int>& nodeID2Index,
                                   const std::vector<double>& baseMatPointVals,
                                   std::vector<double>& bcMatPointVals,
                                   std::vector<int>& bcMask) {
        const auto& mesh = meshMgr.mesh();
        const auto& faces = mesh.getFaces();
        const auto& cells = mesh.getCells();
        const VtkBCFieldCache2D field_cache = BuildVtkBCFieldCache2D(fieldMgr);

        bcMatPointVals = baseMatPointVals;
        bcMask.assign(baseMatPointVals.size(), 0);
        if (!binding.bc || !matCellField) return;

        std::vector<double> weighted_sum(baseMatPointVals.size(), 0.0);
        std::vector<double> weighted_area(baseMatPointVals.size(), 0.0);
        int near_singular_count = 0;
        int invalid_value_count = 0;
        int area_fallback_count = 0;

        for (const auto& face : faces) {
            if (!face.isBoundary()) continue;
            if (!binding.bc->HasBC(face.physicalGroupId)) continue;

            const int owner = face.ownerCell_index;
            if (owner < 0 || owner >= static_cast<int>(cells.size())) continue;

            const auto bc = binding.bc->GetBCCoefficients(face.physicalGroupId, face.midpoint);
            const double phi_owner = SafeFieldValue(matCellField, owner, 0.0);

            const Vector d_owner = face.midpoint - cells[owner].center;
            const double dist = std::max(std::abs(d_owner * face.normal), ctx.geom_floor);
            double area = std::max(face.vectorE.Mag(), 0.0);
            if (area <= ctx.coeff_floor) {
                area = std::max(face.length, 0.0);
                ++area_fallback_count;
            }
            if (area <= ctx.geom_floor) continue;
            const double T_geom = area / dist;

            const double coeff = ComputeTransportCoeff2D(binding, field_cache, owner, face.normal, ctx.coeff_floor);
            const double K = std::max(coeff * T_geom, ctx.coeff_floor);
            const double denom = bc.a - bc.b * K;

            double phi_boundary = phi_owner;
            if (std::isfinite(denom) && std::abs(denom) > ctx.denom_eps) {
                phi_boundary = (bc.c - bc.b * K * phi_owner) / denom;
            }
            else {
                ++near_singular_count;
            }
            if (!std::isfinite(phi_boundary)) {
                phi_boundary = phi_owner;
                ++invalid_value_count;
            }

            for (int nodeId : face.FaceNodeIDs) {
                auto it = nodeID2Index.find(nodeId);
                if (it == nodeID2Index.end()) continue;
                const int pIdx = it->second;
                if (pIdx < 0 || pIdx >= static_cast<int>(bcMatPointVals.size())) continue;
                weighted_sum[pIdx] += area * phi_boundary;
                weighted_area[pIdx] += area;
                bcMask[pIdx] = 1;
            }
        }

        for (size_t i = 0; i < bcMatPointVals.size(); ++i) {
            if (weighted_area[i] > 0.0) {
                bcMatPointVals[i] = weighted_sum[i] / weighted_area[i];
            }
        }
        if (near_singular_count > 0 || invalid_value_count > 0 || area_fallback_count > 0) {
            std::cout << "[PostProcess_2D][BC-Reconstruct] field=" << binding.field_name
                << " singular_denom=" << near_singular_count
                << " invalid_phi=" << invalid_value_count
                << " area_fallback=" << area_fallback_count << "\n";
        }
    }

    void ProjectBoundaryOwnerField2D(const MeshManager& meshMgr,
                                     const std::shared_ptr<volScalarField>& matCellField,
                                     const BoundarySetting::BoundaryConditionManager* bcMgr,
                                     const std::unordered_map<int, int>& nodeID2Index,
                                     const std::vector<double>& baseMatPointVals,
                                     std::vector<double>& outMatPointVals,
                                     std::vector<int>& outMask,
                                     double geom_floor,
                                     double coeff_floor) {
        outMatPointVals = baseMatPointVals;
        outMask.assign(baseMatPointVals.size(), 0);
        if (!matCellField) return;

        const auto& mesh = meshMgr.mesh();
        const auto& faces = mesh.getFaces();
        const auto& cells = mesh.getCells();

        std::vector<double> weighted_sum(baseMatPointVals.size(), 0.0);
        std::vector<double> weighted_area(baseMatPointVals.size(), 0.0);

        for (const auto& face : faces) {
            if (!face.isBoundary()) continue;
            if (bcMgr && !bcMgr->HasBC(face.physicalGroupId)) continue;

            const int owner = face.ownerCell_index;
            if (owner < 0 || owner >= static_cast<int>(cells.size())) continue;
            const double phi_owner = SafeFieldValue(matCellField, owner, 0.0);

            double area = std::max(face.vectorE.Mag(), 0.0);
            if (area <= coeff_floor) area = std::max(face.length, 0.0);
            if (area <= geom_floor) continue;

            for (int nodeId : face.FaceNodeIDs) {
                auto it = nodeID2Index.find(nodeId);
                if (it == nodeID2Index.end()) continue;
                const int pIdx = it->second;
                if (pIdx < 0 || pIdx >= static_cast<int>(outMatPointVals.size())) continue;
                weighted_sum[pIdx] += area * phi_owner;
                weighted_area[pIdx] += area;
                outMask[pIdx] = 1;
            }
        }

        for (size_t i = 0; i < outMatPointVals.size(); ++i) {
            if (weighted_area[i] > 0.0) {
                outMatPointVals[i] = weighted_sum[i] / weighted_area[i];
            }
        }
    }

    inline bool IsSameField(const std::string& a, const std::string& b) {
        return a == b;
    }

    bool BuildDerivedFluidBCPointField2D(
        const std::string& fieldName,
        const VTKBoundaryVisualizationContext& ctx,
        const std::vector<double>& pBase,
        const std::vector<double>& tBase,
        const std::vector<double>& swBase,
        const std::vector<double>* pBC,
        const std::vector<double>* tBC,
        const std::vector<double>* swBC,
        const std::vector<int>& boundaryMask,
        std::vector<double>& outVals)
    {
        const PhysicalProperties_string_op::Water w;
        const PhysicalProperties_string_op::CO2 g;
        const PhysicalProperties_string_op::TwoPhaseState_String tp;

        auto hasPressure = [&]() { return pBC && pBC->size() == outVals.size(); };
        auto hasTemp = [&]() { return tBC && tBC->size() == outVals.size(); };
        auto hasSw = [&]() { return swBC && swBC->size() == outVals.size(); };

        const bool isWaterField =
            IsSameField(fieldName, w.rho_tag) || IsSameField(fieldName, w.mu_tag) ||
            IsSameField(fieldName, w.h_tag) || IsSameField(fieldName, w.cp_tag) ||
            IsSameField(fieldName, w.cv_tag) || IsSameField(fieldName, w.k_tag) ||
            IsSameField(fieldName, w.drho_dp_tag) || IsSameField(fieldName, w.c_w_tag) ||
            IsSameField(fieldName, w.k_rw_tag) || IsSameField(fieldName, w.lambda_w_tag);

        const bool isGasField =
            IsSameField(fieldName, g.rho_tag) || IsSameField(fieldName, g.mu_tag) ||
            IsSameField(fieldName, g.h_tag) || IsSameField(fieldName, g.cp_tag) ||
            IsSameField(fieldName, g.cv_tag) || IsSameField(fieldName, g.k_tag) ||
            IsSameField(fieldName, g.drho_dp_tag) || IsSameField(fieldName, g.c_g_tag) ||
            IsSameField(fieldName, g.k_rg_tag) || IsSameField(fieldName, g.lambda_g_tag);

        const bool isSatAlias =
            IsSameField(fieldName, tp.p_g_field) ||
            IsSameField(fieldName, "S_w") ||
            IsSameField(fieldName, "s_w");

        if (!isWaterField && !isGasField && !isSatAlias) return false;

        const CapRelPerm::VGParams vg;
        const CapRelPerm::RelPermParams rp;
        const double eps = std::max(1.0e-20, ctx.coeff_floor);

        for (size_t i = 0; i < outVals.size(); ++i) {
            if (i >= boundaryMask.size() || boundaryMask[i] == 0) continue;

            const double Pw = hasPressure() ? (*pBC)[i] : ((i < pBase.size()) ? pBase[i] : 0.0);
            const double T = hasTemp() ? (*tBC)[i] : ((i < tBase.size()) ? tBase[i] : 300.0);
            const double Sw_raw = hasSw() ? (*swBC)[i] : ((i < swBase.size()) ? swBase[i] : 1.0);
            const double Sw = std::max(0.0, std::min(1.0, Sw_raw));

            ADVar<1> Pw_ad(Pw); Pw_ad.grad(0) = 1.0;
            ADVar<1> T_ad(T);
            const auto wProps = AD_Fluid::Evaluator::evaluateWater<1>(Pw_ad, T_ad);

            ADVar<1> Sw_ad(Sw);
            ADVar<1> krw_ad, krg_ad;
            CapRelPerm::kr_Mualem_vG<1>(Sw_ad, vg, rp, krw_ad, krg_ad);
            const double krw = std::max(0.0, std::min(1.0, krw_ad.val));
            const double krg = std::max(0.0, std::min(1.0, krg_ad.val));

            const ADVar<1> Pc_ad = CapRelPerm::pc_vG<1>(Sw_ad, vg);
            ADVar<1> Pg_ad(Pw + Pc_ad.val); Pg_ad.grad(0) = 1.0;
            const auto gProps = AD_Fluid::Evaluator::evaluateCO2<1>(Pg_ad, T_ad);

            if (IsSameField(fieldName, "S_w") || IsSameField(fieldName, "s_w")) {
                outVals[i] = Sw;
            }
            else if (IsSameField(fieldName, tp.p_g_field)) {
                outVals[i] = Pg_ad.val;
            }
            else if (IsSameField(fieldName, w.rho_tag)) outVals[i] = wProps.rho.val;
            else if (IsSameField(fieldName, w.mu_tag)) outVals[i] = wProps.mu.val;
            else if (IsSameField(fieldName, w.h_tag)) outVals[i] = wProps.h.val;
            else if (IsSameField(fieldName, w.cp_tag)) outVals[i] = wProps.cp.val;
            else if (IsSameField(fieldName, w.cv_tag)) outVals[i] = wProps.cv.val;
            else if (IsSameField(fieldName, w.k_tag)) outVals[i] = wProps.k.val;
            else if (IsSameField(fieldName, w.drho_dp_tag)) outVals[i] = wProps.rho.grad[0];
            else if (IsSameField(fieldName, w.c_w_tag)) outVals[i] = wProps.rho.grad[0] / std::max(wProps.rho.val, eps);
            else if (IsSameField(fieldName, w.k_rw_tag)) outVals[i] = krw;
            else if (IsSameField(fieldName, w.lambda_w_tag)) outVals[i] = krw / std::max(wProps.mu.val, eps);
            else if (IsSameField(fieldName, g.rho_tag)) outVals[i] = gProps.rho.val;
            else if (IsSameField(fieldName, g.mu_tag)) outVals[i] = gProps.mu.val;
            else if (IsSameField(fieldName, g.h_tag)) outVals[i] = gProps.h.val;
            else if (IsSameField(fieldName, g.cp_tag)) outVals[i] = gProps.cp.val;
            else if (IsSameField(fieldName, g.cv_tag)) outVals[i] = gProps.cv.val;
            else if (IsSameField(fieldName, g.k_tag)) outVals[i] = gProps.k.val;
            else if (IsSameField(fieldName, g.drho_dp_tag)) outVals[i] = gProps.rho.grad[0];
            else if (IsSameField(fieldName, g.c_g_tag)) outVals[i] = gProps.rho.grad[0] / std::max(gProps.rho.val, eps);
            else if (IsSameField(fieldName, g.k_rg_tag)) outVals[i] = krg;
            else if (IsSameField(fieldName, g.lambda_g_tag)) outVals[i] = krg / std::max(gProps.mu.val, eps);
        }

        return true;
    }
}

// ==============================================================================
// 锟斤拷锟侥碉拷锟斤拷锟斤拷锟斤拷
// ==============================================================================
void PostProcess_2D::ExportTecplot(const std::string& filename, double time) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("[PostProcess_2D Error] Failed to open file for Tecplot export: " + filename);
    }

    std::cout << ">>> Exporting 2D Tecplot data to: " << filename << " (Time = " << time << ") ..." << std::endl;

    // 1. 锟斤拷取全锟斤拷唯一锟斤拷锟斤拷锟斤拷
    std::vector<std::string> varNames = GetAllUniqueFieldNames();

    // 2. 写锟斤拷全锟斤拷 Header
    file << "TITLE = \"2D EDFM Multiphase Simulation Result\"\n";
    file << "VARIABLES = \"X\", \"Y\"";
    for (const auto& name : varNames) {
        file << ", \"" << name << "\"";
    }
    file << "\n";

    // 锟斤拷锟矫匡拷学锟斤拷锟斤拷锟斤拷锟斤拷锟竭撅拷锟斤拷锟斤拷锟?
    file << std::scientific << std::setprecision(8);

    // =========================================================
    // Zone 1: Matrix Domain (锟斤拷锟斤拷锟斤拷 - 2D Cell)
    // =========================================================
    const auto& matrixNodes = meshMgr_.mesh().getNodes();
    const auto& matrixCells = meshMgr_.mesh().getCells();
    size_t numMatrixNodes = matrixNodes.size();
    size_t numMatrixCells = matrixCells.size();

    if (numMatrixCells > 0) {
        // 锟叫讹拷锟斤拷锟斤拷锟斤拷锟斤拷锟?Tecplot 锟斤拷元锟斤拷锟斤拷 (Tri / Quad)
        int maxNodesPerCell = 0;
        for (const auto& cell : matrixCells) {
            maxNodesPerCell = std::max(maxNodesPerCell, static_cast<int>(cell.CellNodeIDs.size()));
        }
        std::string matrixZoneType = GetTecplotElementType(maxNodesPerCell);

        // 写锟斤拷 Matrix Zone Header (BLOCK 锟斤拷式)
        file << "ZONE T=\"Matrix\", SOLUTIONTIME=" << time << ", DATAPACKING=BLOCK, "
            << "NODES=" << numMatrixNodes << ", ELEMENTS=" << numMatrixCells << ", "
            << "ZONETYPE=" << matrixZoneType;

        // 锟斤拷锟斤拷锟斤拷些锟斤拷锟斤拷锟斤拷锟节碉拷元锟斤拷锟斤拷 (X,Y 为 Node锟斤拷锟斤拷锟斤拷为 Cell-Centered)
        if (!varNames.empty()) {
            file << ", VARLOCATION=([3-" << 2 + varNames.size() << "]=CELLCENTERED)";
        }
        file << "\n";

        // (1.1) 锟斤拷锟斤拷锟斤拷医诘锟?X 锟斤拷锟斤拷
        for (const auto& node : matrixNodes) {
            file << node.coord.m_x << "\n";
        }
        // (1.2) 锟斤拷锟斤拷锟斤拷医诘锟?Y 锟斤拷锟斤拷
        for (const auto& node : matrixNodes) {
            file << node.coord.m_y << "\n";
        }

        // (1.3) 锟斤拷锟斤拷锟斤拷冶锟斤拷锟斤拷锟?(CELLCENTERED)
        for (const auto& name : varNames) {
            auto field = fieldMgr_.getMatrixScalar(name);
            if (field) {
                // 锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟绞碉拷锟斤拷锟?
                for (size_t i = 0; i < numMatrixCells; ++i) {
                    file << field->data[i] << "\n";
                }
            }
            else {
                // 锟斤拷锟斤拷锟斤拷锟斤拷冢锟斤拷锟斤拷锟街伙拷锟斤拷逊锟斤拷锟斤拷械锟斤拷锟斤拷裕锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷
                for (size_t i = 0; i < numMatrixCells; ++i) {
                    file << 0.0 << "\n";
                }
            }
        }

        // (1.4) 锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟酵拷锟?(Connectivity)
        // 锟斤拷锟斤拷 Gmsh Global Node ID 锟斤拷 锟街诧拷锟斤拷锟斤拷锟斤拷锟斤拷锟?1-based Index 映锟斤拷
        std::unordered_map<int, int> gmshIdToLocalIdx;
        for (size_t i = 0; i < numMatrixNodes; ++i) {
            gmshIdToLocalIdx[matrixNodes[i].id] = static_cast<int>(i) + 1; // Tecplot 锟斤拷 1-based
        }

        for (const auto& cell : matrixCells) {
            size_t nNodes = cell.CellNodeIDs.size();
            for (int i = 0; i < maxNodesPerCell; ++i) {
                // 锟斤拷锟斤拷锟角帮拷锟皆拷锟斤拷锟斤拷锟斤拷危锟?锟斤拷锟节点）锟斤拷锟斤拷锟斤拷锟斤拷 Zone 锟斤拷锟侥憋拷锟轿ｏ拷max=4锟斤拷锟斤拷锟斤拷锟截革拷锟斤拷锟揭伙拷锟斤拷诘锟?
                int nodeGmshId = (i < nNodes) ? cell.CellNodeIDs[i] : cell.CellNodeIDs.back();

                // 锟斤拷全锟斤拷取锟斤拷锟斤拷锟?1-based 锟斤拷锟斤拷
                auto it = gmshIdToLocalIdx.find(nodeGmshId);
                if (it != gmshIdToLocalIdx.end()) {
                    file << it->second << " ";
                }
                else {
                    file << "1 "; // 锟斤拷锟剿凤拷锟斤拷锟襟备凤拷锟斤拷
                }
            }
            file << "\n";
        }
    }

    // =========================================================
    // Zone 2: Fracture Domain (锟窖凤拷锟斤拷 - 1D Segments)
    // =========================================================
    const auto& fractures = meshMgr_.fracture_network().fractures;

    // 统锟斤拷锟窖凤拷微元锟斤拷
    size_t totalFracCells = 0;
    for (const auto& frac : fractures) {
        totalFracCells += frac.elements.size();
    }

    // 锟斤拷锟矫讹拷锟斤拷锟节碉拷锟斤拷裕锟矫匡拷锟轿⒃拷锟教拷锟斤拷锟?2 锟斤拷锟剿点，锟杰节碉拷锟斤拷为微元锟斤拷 * 2
    size_t totalFracNodes = totalFracCells * 2;

    if (totalFracCells > 0) {
        // 写锟斤拷 Fracture Zone Header
        file << "ZONE T=\"Fractures\", SOLUTIONTIME=" << time << ", DATAPACKING=BLOCK, "
            << "NODES=" << totalFracNodes << ", ELEMENTS=" << totalFracCells << ", "
            << "ZONETYPE=" << GetTecplotElementType(2);

        if (!varNames.empty()) {
            file << ", VARLOCATION=([3-" << 2 + varNames.size() << "]=CELLCENTERED)";
        }
        file << "\n";

        // 锟秸硷拷锟斤拷锟叫碉拷锟斤拷锟斤拷锟节碉拷锟斤拷锟斤拷 (锟斤拷锟斤拷锟斤拷循锟斤拷锟斤拷锟斤拷 I/O 锟斤拷位)
        std::vector<double> fracX;
        std::vector<double> fracY;
        fracX.reserve(totalFracNodes);
        fracY.reserve(totalFracNodes);

        for (const auto& frac : fractures) {
            // 锟斤拷锟斤拷锟斤拷锟窖凤拷锟?X, Y 锟斤拷锟?
            double dx = frac.end.m_x - frac.start.m_x;
            double dy = frac.end.m_y - frac.start.m_y;

            for (const auto& cell : frac.elements) {
                // 锟斤拷锟节癸拷一锟斤拷锟斤拷锟斤拷 [param0, param1] 锟斤拷锟皆诧拷值锟斤拷锟斤拷微元锟剿碉拷锟斤拷锟斤拷
                fracX.push_back(frac.start.m_x + dx * cell.param0); // 锟节碉拷1 X
                fracX.push_back(frac.start.m_x + dx * cell.param1); // 锟节碉拷2 X

                fracY.push_back(frac.start.m_y + dy * cell.param0); // 锟节碉拷1 Y
                fracY.push_back(frac.start.m_y + dy * cell.param1); // 锟节碉拷2 Y
            }
        }

        // (2.1) 锟斤拷锟斤拷逊锟节碉拷 X 锟斤拷锟斤拷
        for (double x : fracX) file << x << "\n";

        // (2.2) 锟斤拷锟斤拷逊锟节碉拷 Y 锟斤拷锟斤拷
        for (double y : fracY) file << y << "\n";

        // (2.3) 锟斤拷锟斤拷逊锟斤拷锟斤拷锟斤拷 (CELLCENTERED)
        for (const auto& name : varNames) {
            auto field = fieldMgr_.getFractureScalar(name);
            if (field) {
                for (size_t i = 0; i < totalFracCells; ++i) {
                    file << field->data[i] << "\n";
                }
            }
            else {
                for (size_t i = 0; i < totalFracCells; ++i) {
                    file << 0.0 << "\n";
                }
            }
        }

        // (2.4) 锟斤拷锟斤拷逊锟斤拷锟酵拷锟?(Connectivity)
        // 锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟剿ｏ拷锟斤拷 i 锟斤拷微元锟较革拷锟斤拷锟接节碉拷 (2*i + 1) 锟斤拷 (2*i + 2) [Tecplot要锟斤拷1-based]
        for (size_t i = 0; i < totalFracCells; ++i) {
            file << (2 * i + 1) << " " << (2 * i + 2) << "\n";
        }
    }

    file.close();
    std::cout << ">>> 2D Tecplot export completed successfully." << std::endl;
}

// ==============================================================================
// 锟斤拷锟斤拷锟斤拷 ParaView (VTK)
// ==============================================================================
void PostProcess_2D::ExportVTK(const std::string& filename, double time) const
{
    std::cout << "[PostProcess_2D] Exporting ParaView VTK file: " << filename << " (t=" << time << ")" << std::endl;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << std::endl;
        return;
    }

    out << std::scientific << std::setprecision(6);
    std::vector<std::string> varNames = GetAllUniqueFieldNames();

    const auto& matrixNodes = meshMgr_.mesh().getNodes();
    const auto& matrixCells = meshMgr_.mesh().getCells();
    size_t numMatrixNodes = matrixNodes.size();
    size_t numMatrixCells = matrixCells.size();

    // 锟斤拷锟斤拷锟斤拷锟揭节点到 0-based 锟斤拷锟斤拷锟斤拷映锟斤拷 (VTK 锟斤拷锟斤拷锟?0 锟斤拷始)
    std::unordered_map<int, int> nodeID2Index;
    for (size_t i = 0; i < numMatrixNodes; ++i) {
        nodeID2Index[matrixNodes[i].id] = static_cast<int>(i);
    }

    const auto& fractures = meshMgr_.fracture_network().fractures;
    size_t totalFracCells = 0;
    for (const auto& frac : fractures) {
        totalFracCells += frac.elements.size();
    }
    // 1D 锟窖凤拷锟斤拷 VTK 锟斤拷锟斤拷要锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷锟节碉拷锟斤拷锟斤拷染锟竭讹拷
    size_t totalFracNodes = totalFracCells * 2;

    size_t totalNodes = numMatrixNodes + totalFracNodes;
    size_t totalCells = numMatrixCells + totalFracCells;

    // 锟斤拷锟斤拷 CELLS 锟叫憋拷锟斤拷锟杰筹拷锟斤拷
    size_t cellListSize = 0;
    for (const auto& cell : matrixCells) {
        cellListSize += (1 + cell.CellNodeIDs.size());
    }
    cellListSize += totalFracCells * 3; // 锟窖凤拷锟竭讹拷: 1(锟节碉拷锟斤拷2) + 2(锟节碉拷锟斤拷锟斤拷) = 3

    // 1. VTK Header
    out << "# vtk DataFile Version 3.0\n";
    out << "2D-EDFM Simulation Result (Time: " << time << ")\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // 2. POINTS
    out << "POINTS " << totalNodes << " double\n";
    int col = 0;
    for (const auto& node : matrixNodes) {
        out << node.coord.m_x << " " << node.coord.m_y << " 0.0 "; // 2D 锟斤拷锟斤拷 Z=0
        if (++col % 3 == 0) out << "\n";
    }
    for (const auto& frac : fractures) {
        double dx = frac.end.m_x - frac.start.m_x;
        double dy = frac.end.m_y - frac.start.m_y;
        for (const auto& cell : frac.elements) {
            out << (frac.start.m_x + dx * cell.param0) << " " << (frac.start.m_y + dy * cell.param0) << " 0.0 ";
            if (++col % 3 == 0) out << "\n";
            out << (frac.start.m_x + dx * cell.param1) << " " << (frac.start.m_y + dy * cell.param1) << " 0.0 ";
            if (++col % 3 == 0) out << "\n";
        }
    }
    if (col % 3 != 0) out << "\n";

    // 3. CELLS
    out << "CELLS " << totalCells << " " << cellListSize << "\n";
    for (const auto& cell : matrixCells) {
        out << cell.CellNodeIDs.size() << " ";
        for (int id : cell.CellNodeIDs) {
            auto it = nodeID2Index.find(id);
            if (it != nodeID2Index.end()) {
                out << it->second << " ";
            }
            else {
                // [锟较革拷模式] 锟斤拷锟斤拷缺失锟节碉拷直锟接憋拷锟斤拷锟桔断ｏ拷锟斤拷锟斤拷锟节革拷锟斤拷锟斤拷锟斤拷锟剿达拷锟斤拷
                out.close();
                throw std::runtime_error("[PostProcess_2D Error] Node ID " + std::to_string(id) + " not found in Matrix Mesh!");
            }
        }
        out << "\n";
    }

    size_t fracNodeOffset = numMatrixNodes;
    for (size_t i = 0; i < totalFracCells; ++i) {
        out << "2 " << fracNodeOffset << " " << (fracNodeOffset + 1) << "\n";
        fracNodeOffset += 2;
    }

    // 4. CELL_TYPES
    out << "CELL_TYPES " << totalCells << "\n";
    for (const auto& cell : matrixCells) {
        size_t n = cell.CellNodeIDs.size();
        if (n == 3) {
            out << "5\n";      // VTK_TRIANGLE
        }
        else if (n == 4) {
            out << "9\n";      // VTK_QUAD
        }
        else {
            // [锟较革拷模式] 锟杰撅拷 2D Matrix 锟叫筹拷锟街凤拷 3/4 锟节碉拷锟斤拷嘶锟斤拷锟斤拷锟斤拷蔚锟皆拷锟斤拷锟街癸拷锟斤拷锟斤拷锟斤拷锟侥ｏ拷锟?
            out.close();
            throw std::runtime_error("[PostProcess_2D Error] Invalid matrix cell detected with " + std::to_string(n) + " nodes. Only Triangle (3) and Quad (4) are strictly supported.");
        }
    }
    for (size_t i = 0; i < totalFracCells; ++i) {
        out << "3\n";                  // VTK_LINE
    }

    // 5. CELL_DATA (cell-centered fields for accurate well/fracture visualization)
    out << "CELL_DATA " << totalCells << "\n";

    // Domain tag (0: matrix, 1: fracture)
    out << "SCALARS DomainID int 1\nLOOKUP_TABLE default\n";
    for (size_t i = 0; i < numMatrixCells; ++i) out << "0\n";
    for (size_t i = 0; i < totalFracCells; ++i) out << "1\n";

    // Cell-centered physical fields (no interpolation 鈥?exact solver values)
    for (const auto& name : varNames) {
        out << "SCALARS cell_" << name << " double 1\nLOOKUP_TABLE default\n";
        auto matField = fieldMgr_.getMatrixScalar(name);
        for (size_t c = 0; c < numMatrixCells; ++c) {
            double val = (matField && c < matField->data.size()) ? matField->data[c] : 0.0;
            if (!std::isfinite(val)) val = 0.0;
            out << val << "\n";
        }
        auto fracField = fieldMgr_.getFractureScalar(name);
        for (size_t c = 0; c < totalFracCells; ++c) {
            double val = (fracField && c < fracField->data.size()) ? fracField->data[c] : 0.0;
            if (!std::isfinite(val)) val = 0.0;
            out << val << "\n";
        }
    }

    // 6. POINT_DATA (node-averaged interpolation)
    out << "POINT_DATA " << totalNodes << "\n";

    // Build base point fields first so BC-aware fields can reuse the same baseline.
    std::unordered_map<std::string, std::vector<double>> pointFieldCache;
    pointFieldCache.reserve(varNames.size());

    for (const auto& name : varNames) {
        std::vector<double> pointVals(totalNodes, 0.0);
        std::vector<int> matCounts(numMatrixNodes, 0);
        auto matField = fieldMgr_.getMatrixScalar(name);

        if (matField) {
            for (size_t c = 0; c < numMatrixCells; ++c) {
                double val = (c < matField->data.size()) ? matField->data[c] : 0.0;
                if (!std::isfinite(val)) val = 0.0;
                for (int id : matrixCells[c].CellNodeIDs) {
                    auto it = nodeID2Index.find(id);
                    if (it != nodeID2Index.end()) {
                        const int pIdx = it->second;
                        if (pIdx >= 0 && pIdx < static_cast<int>(numMatrixNodes)) {
                            pointVals[pIdx] += val;
                            matCounts[pIdx]++;
                        }
                    }
                }
            }
            for (size_t p = 0; p < numMatrixNodes; ++p) {
                if (matCounts[p] > 0) pointVals[p] /= matCounts[p];
            }
        }

        auto fracField = fieldMgr_.getFractureScalar(name);
        size_t fracPointOffset = numMatrixNodes;
        if (fracField) {
            for (size_t c = 0; c < totalFracCells; ++c) {
                double val = (c < fracField->data.size()) ? fracField->data[c] : 0.0;
                if (!std::isfinite(val)) val = 0.0;
                if (fracPointOffset + 1 < pointVals.size()) {
                    pointVals[fracPointOffset] = val;
                    pointVals[fracPointOffset + 1] = val;
                }
                fracPointOffset += 2;
            }
        }

        pointFieldCache.emplace(name, std::move(pointVals));
    }

    // Original point fields (backward-compatible)
    for (const auto& name : varNames) {
        out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
        const auto it = pointFieldCache.find(name);
        if (it != pointFieldCache.end()) {
            for (double v : it->second) out << v << "\n";
        }
        else {
            for (size_t p = 0; p < totalNodes; ++p) out << "0.0\n";
        }
    }

    // BC-aware point fields for all exported variables: "<field>_bc_point" (+ "<field>_bc_mask")
    if (bcVizCtx_ && !bcVizCtx_->bindings.empty()) {
        const auto pEqCfg = PhysicalProperties_string_op::PressureEquation_String::FIM();
        const auto tEqCfg = PhysicalProperties_string_op::TemperatureEquation_String::FIM();
        const auto sEqCfg = PhysicalProperties_string_op::SaturationEquation_String::FIM();

        const auto* pressureBinding = bcVizCtx_->findFirstByKind(VTKBCTransportKind::Pressure);
        const auto* temperatureBinding = bcVizCtx_->findFirstByKind(VTKBCTransportKind::Temperature);
        const auto* saturationBinding = bcVizCtx_->findFirstByKind(VTKBCTransportKind::Saturation);

        auto loadBaseMatPoint = [&](const std::string& fieldName, std::vector<double>& outBase) -> bool {
            auto it = pointFieldCache.find(fieldName);
            if (it == pointFieldCache.end()) return false;
            if (it->second.size() < numMatrixNodes) return false;
            outBase.assign(numMatrixNodes, 0.0);
            for (size_t i = 0; i < numMatrixNodes; ++i) outBase[i] = it->second[i];
            return true;
            };

        auto reconstructPrimary = [&](const std::string& fieldName,
                                      const VTKBCVariableBinding* binding,
                                      std::vector<double>& baseVals,
                                      std::vector<double>& bcVals,
                                      std::vector<int>& bcMask) -> bool {
            if (!binding) return false;
            auto matField = fieldMgr_.getMatrixScalar(fieldName);
            if (!matField) return false;
            if (!loadBaseMatPoint(fieldName, baseVals)) return false;
            ReconstructBCPointField2D(meshMgr_, fieldMgr_, matField, *bcVizCtx_, *binding, nodeID2Index, baseVals, bcVals, bcMask);
            return true;
            };

        auto markBoundaryNodes = [&](const BoundarySetting::BoundaryConditionManager* bcMgr, std::vector<int>& mask) {
            const auto& faces = meshMgr_.mesh().getFaces();
            for (const auto& face : faces) {
                if (!face.isBoundary()) continue;
                if (bcMgr && !bcMgr->HasBC(face.physicalGroupId)) continue;
                for (int nodeId : face.FaceNodeIDs) {
                    auto it = nodeID2Index.find(nodeId);
                    if (it == nodeID2Index.end()) continue;
                    const int idx = it->second;
                    if (idx >= 0 && idx < static_cast<int>(mask.size())) mask[idx] = 1;
                }
            }
            };

        std::vector<double> pBaseMat, pBcMat, tBaseMat, tBcMat, swBaseMat, swBcMat;
        std::vector<int> pMaskMat, tMaskMat, swMaskMat;

        bool hasP = reconstructPrimary(pEqCfg.pressure_field, pressureBinding, pBaseMat, pBcMat, pMaskMat);
        if (!hasP) {
            hasP = reconstructPrimary("p", pressureBinding, pBaseMat, pBcMat, pMaskMat);
        }
        const bool hasT = reconstructPrimary(tEqCfg.temperatue_field, temperatureBinding, tBaseMat, tBcMat, tMaskMat);
        bool hasSw = reconstructPrimary(sEqCfg.saturation, saturationBinding, swBaseMat, swBcMat, swMaskMat);
        if (!hasSw) {
            hasSw = reconstructPrimary("S_w", saturationBinding, swBaseMat, swBcMat, swMaskMat);
        }

        std::vector<int> boundaryUnionMask(numMatrixNodes, 0);
        for (size_t i = 0; i < numMatrixNodes; ++i) {
            const int mp = (i < pMaskMat.size()) ? pMaskMat[i] : 0;
            const int mt = (i < tMaskMat.size()) ? tMaskMat[i] : 0;
            const int ms = (i < swMaskMat.size()) ? swMaskMat[i] : 0;
            boundaryUnionMask[i] = (mp || mt || ms) ? 1 : 0;
        }
        if (std::all_of(boundaryUnionMask.begin(), boundaryUnionMask.end(), [](int v) { return v == 0; })) {
            if (pressureBinding) markBoundaryNodes(pressureBinding->bc, boundaryUnionMask);
            if (temperatureBinding) markBoundaryNodes(temperatureBinding->bc, boundaryUnionMask);
            if (saturationBinding) markBoundaryNodes(saturationBinding->bc, boundaryUnionMask);
        }

        const BoundarySetting::BoundaryConditionManager* fallbackBcMgr = nullptr;
        if (pressureBinding) fallbackBcMgr = pressureBinding->bc;
        else if (temperatureBinding) fallbackBcMgr = temperatureBinding->bc;
        else if (saturationBinding) fallbackBcMgr = saturationBinding->bc;

        for (const auto& name : varNames) {
            auto itBase = pointFieldCache.find(name);
            if (itBase == pointFieldCache.end()) continue;

            std::vector<double> baseMatVals(numMatrixNodes, 0.0);
            for (size_t i = 0; i < numMatrixNodes; ++i) baseMatVals[i] = itBase->second[i];

            std::vector<double> bcMatVals = baseMatVals;
            std::vector<int> bcMaskMat(numMatrixNodes, 0);

            const auto* binding = bcVizCtx_->findBinding(name);
            auto matField = fieldMgr_.getMatrixScalar(name);

            if (binding && matField) {
                ReconstructBCPointField2D(meshMgr_, fieldMgr_, matField, *bcVizCtx_, *binding, nodeID2Index, baseMatVals, bcMatVals, bcMaskMat);
            }
            else if (matField) {
                bcMaskMat = boundaryUnionMask;
                const bool derivedDone = BuildDerivedFluidBCPointField2D(
                    name, *bcVizCtx_,
                    pBaseMat, tBaseMat, swBaseMat,
                    hasP ? &pBcMat : nullptr,
                    hasT ? &tBcMat : nullptr,
                    hasSw ? &swBcMat : nullptr,
                    boundaryUnionMask,
                    bcMatVals);

                if (!derivedDone) {
                    ProjectBoundaryOwnerField2D(
                        meshMgr_, matField, fallbackBcMgr, nodeID2Index,
                        baseMatVals, bcMatVals, bcMaskMat,
                        bcVizCtx_->geom_floor, bcVizCtx_->coeff_floor);
                }
            }

            std::vector<double> bcPointVals = itBase->second;
            for (size_t i = 0; i < std::min(numMatrixNodes, bcMatVals.size()); ++i) {
                bcPointVals[i] = bcMatVals[i];
            }

            out << "SCALARS " << name << "_bc_point double 1\nLOOKUP_TABLE default\n";
            for (double v : bcPointVals) out << v << "\n";

            if (bcVizCtx_->export_bc_mask) {
                out << "SCALARS " << name << "_bc_mask int 1\nLOOKUP_TABLE default\n";
                for (size_t i = 0; i < numMatrixNodes; ++i) {
                    const int m = (i < bcMaskMat.size()) ? bcMaskMat[i] : 0;
                    out << m << "\n";
                }
                for (size_t i = numMatrixNodes; i < totalNodes; ++i) out << "0\n";
            }
        }
    }
    out.close();
    std::cout << "[PostProcess_2D] VTK Export completed successfully with Point Interpolation.\n"; }

// ==============================================================================
// ExportVTU  (Step 4: ParaView XML Unstructured Grid for time-series animation)
// ==============================================================================
void PostProcess_2D::ExportVTU(const std::string& filename, double time) const
{
    std::cout << "[PostProcess_2D] Exporting VTU: " << filename << "  (t=" << time << ")" << std::endl;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << std::endl;
        return;
    }
    out << std::scientific << std::setprecision(8);

    const std::vector<std::string> varNames = GetAllUniqueFieldNames();
    const auto& matrixNodes = meshMgr_.mesh().getNodes();
    const auto& matrixCells = meshMgr_.mesh().getCells();
    const size_t numMatrixNodes = matrixNodes.size();
    const size_t numMatrixCells = matrixCells.size();

    // 0-based node ID 鈫?point-array index (VTK requires 0-based)
    std::unordered_map<int, int> nodeID2Index;
    nodeID2Index.reserve(numMatrixNodes * 2);
    for (size_t i = 0; i < numMatrixNodes; ++i)
        nodeID2Index[matrixNodes[i].id] = static_cast<int>(i);

    const auto& fractures = meshMgr_.fracture_network().fractures;
    size_t totalFracCells = 0;
    for (const auto& frac : fractures) totalFracCells += frac.elements.size();
    const size_t totalFracNodes = totalFracCells * 2;  // 2 endpoints per segment

    const size_t totalPoints = numMatrixNodes + totalFracNodes;
    const size_t totalCells  = numMatrixCells + totalFracCells;

    // 鈹€鈹€ XML header 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <!-- Time = " << time << " s -->\n";
    out << "    <Piece NumberOfPoints=\"" << totalPoints
        << "\" NumberOfCells=\"" << totalCells << "\">\n";

    // 鈹€鈹€ Points 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& nd : matrixNodes)
        out << "          " << nd.coord.m_x << " " << nd.coord.m_y << " 0.0\n";
    for (const auto& frac : fractures) {
        const double dx = frac.end.m_x - frac.start.m_x;
        const double dy = frac.end.m_y - frac.start.m_y;
        for (const auto& elem : frac.elements) {
            out << "          "
                << (frac.start.m_x + dx * elem.param0) << " "
                << (frac.start.m_y + dy * elem.param0) << " 0.0\n";
            out << "          "
                << (frac.start.m_x + dx * elem.param1) << " "
                << (frac.start.m_y + dy * elem.param1) << " 0.0\n";
        }
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // 鈹€鈹€ Cells 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    out << "      <Cells>\n";

    // connectivity
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& cell : matrixCells) {
        out << "          ";
        for (int id : cell.CellNodeIDs) {
            auto it = nodeID2Index.find(id);
            out << (it != nodeID2Index.end() ? it->second : 0) << " ";
        }
        out << "\n";
    }
    size_t fracPtBase = numMatrixNodes;
    for (size_t i = 0; i < totalFracCells; ++i) {
        out << "          " << (fracPtBase + 2*i) << " " << (fracPtBase + 2*i + 1) << "\n";
    }
    out << "        </DataArray>\n";

    // offsets (cumulative node count per cell)
    out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    out << "          ";
    int offset = 0;
    for (const auto& cell : matrixCells) {
        offset += static_cast<int>(cell.CellNodeIDs.size());
        out << offset << " ";
    }
    for (size_t i = 0; i < totalFracCells; ++i) {
        offset += 2;
        out << offset << " ";
    }
    out << "\n        </DataArray>\n";

    // types
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    out << "          ";
    for (const auto& cell : matrixCells) {
        const size_t n = cell.CellNodeIDs.size();
        if      (n == 3) out << "5 ";   // VTK_TRIANGLE
        else if (n == 4) out << "9 ";   // VTK_QUAD
        else             out << "7 ";   // VTK_POLYGON (fallback)
    }
    for (size_t i = 0; i < totalFracCells; ++i) out << "3 "; // VTK_LINE
    out << "\n        </DataArray>\n";

    out << "      </Cells>\n";

    // 鈹€鈹€ CellData 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    out << "      <CellData>\n";

    // DomainID (0=matrix, 1=fracture)
    out << "        <DataArray type=\"Int32\" Name=\"DomainID\" format=\"ascii\">\n";
    out << "          ";
    for (size_t i = 0; i < numMatrixCells; ++i) out << "0 ";
    for (size_t i = 0; i < totalFracCells;  ++i) out << "1 ";
    out << "\n        </DataArray>\n";

    // Physical fields
    for (const auto& name : varNames) {
        out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        out << "          ";

        auto matField  = fieldMgr_.getMatrixScalar(name);
        auto fracField = fieldMgr_.getFractureScalar(name);

        for (size_t i = 0; i < numMatrixCells; ++i) {
            double v = (matField && i < matField->data.size()) ? matField->data[i] : 0.0;
            if (!std::isfinite(v)) v = 0.0;
            out << v << " ";
        }
        size_t fi = 0;
        for (const auto& frac : fractures) {
            for (size_t ei = 0; ei < frac.elements.size(); ++ei, ++fi) {
                double v = (fracField && fi < fracField->data.size()) ? fracField->data[fi] : 0.0;
                if (!std::isfinite(v)) v = 0.0;
                out << v << " ";
            }
        }
        out << "\n        </DataArray>\n";
    }
    out << "      </CellData>\n";

    // 鈹€鈹€ Close 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";

    out.close();
    std::cout << "[PostProcess_2D] VTU export completed: " << filename << std::endl;
}

// ==============================================================================
// ExportPVD  (Step 4: ParaView time-series collection file)
// ==============================================================================
/*static*/
void PostProcess_2D::ExportPVD(const std::string& pvd_filename,
                                const std::vector<std::string>& vtu_filenames,
                                const std::vector<double>& times)
{
    if (vtu_filenames.size() != times.size()) {
        std::cerr << "[PostProcess_2D] ExportPVD: vtu_filenames and times size mismatch!\n";
        return;
    }
    std::ofstream out(pvd_filename);
    if (!out.is_open()) {
        std::cerr << "[PostProcess_2D] ExportPVD: cannot open " << pvd_filename << "\n";
        return;
    }
    out << std::scientific << std::setprecision(8);
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    out << "  <Collection>\n";
    for (size_t i = 0; i < vtu_filenames.size(); ++i) {
        out << "    <DataSet timestep=\"" << times[i]
            << "\" group=\"\" part=\"0\" file=\"" << vtu_filenames[i] << "\"/>\n";
    }
    out << "  </Collection>\n";
    out << "</VTKFile>\n";
    out.close();
    std::cout << "[PostProcess_2D] PVD collection written: " << pvd_filename
              << " (" << times.size() << " snapshots)\n";
}

// ==============================================================================
// 锟节诧拷锟斤拷锟斤拷锟斤拷锟斤拷
// ==============================================================================

std::vector<std::string> PostProcess_2D::GetAllUniqueFieldNames() const
{
    std::set<std::string> uniqueNames;

    // 锟斤拷取锟斤拷锟斤拷锟斤拷锟斤拷锟斤拷双锟斤拷锟饺憋拷锟斤拷锟斤拷
    for (const auto& pair : fieldMgr_.matrixFields.fields) {
        // 锟斤拷锟斤拷 dynamic_pointer_cast 锟斤拷全锟斤拷锟?volScalarField (锟脚筹拷 Vector/ADVar 锟斤拷)
        if (std::dynamic_pointer_cast<volScalarField>(pair.second)) {
            uniqueNames.insert(pair.first);
        }
    }

    // 锟斤拷取锟窖凤拷锟斤拷锟斤拷锟斤拷双锟斤拷锟饺憋拷锟斤拷锟斤拷
    for (const auto& pair : fieldMgr_.fractureFields.fields) {
        if (std::dynamic_pointer_cast<volScalarField>(pair.second)) {
            uniqueNames.insert(pair.first);
        }
    }

    return std::vector<std::string>(uniqueNames.begin(), uniqueNames.end());
}

std::string PostProcess_2D::GetTecplotElementType(int numNodes) const
{
    if (numNodes == 2) {
        return "FELINESEG";         // 1D 锟竭讹拷 (锟窖凤拷)
    }
    else if (numNodes == 3) {
        return "FETRIANGLE";        // 锟斤拷锟斤拷锟斤拷锟轿伙拷锟斤拷
    }
    else if (numNodes == 4) {
        return "FEQUADRILATERAL";   // 锟侥憋拷锟轿伙拷锟斤拷 锟斤拷 锟斤拷锟?锟斤拷锟斤拷锟斤拷+锟侥憋拷锟斤拷)
    }

    // 锟斤拷锟斤拷锟斤拷锟斤拷
    return "FEQUADRILATERAL";
}
