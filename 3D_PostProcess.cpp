#include "3D_PostProcess.h"
#include "SolverContrlStrName_op.h"
#include "AD_FluidEvaluator.h"
#include "CapRelPerm_HD_AD.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <set>
#include <unordered_map> 
#include <algorithm> 
#include <limits>

// =========================================================
// 辅助工具：获取 Tecplot 标准 ZoneType 及目标节点数
// =========================================================
struct TecplotElemInfo {
    std::string typeName;
    int targetNodeCount;
};

TecplotElemInfo GetTecplotInfo(const std::string& internalGridType) {
    if (internalGridType == "FETETRA")      return { "FETETRAHEDRON", 4 };
    if (internalGridType == "FEPYRAMID")    return { "FEBRICK", 8 };
    if (internalGridType == "FEPRISM")      return { "FEBRICK", 8 };
    if (internalGridType == "FEBRICK")      return { "FEBRICK", 8 };
    if (internalGridType == "FETRIANGLE")   return { "FETRIANGLE", 3 };
    if (internalGridType == "FEQUAD")       return { "FEQUADRILATERAL", 4 };
    return { "FETETRAHEDRON", 4 };
}

PostProcess_3D::PostProcess_3D(const MeshManager_3D& meshMgr,
                               const FieldManager_3D& fieldMgr,
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

    struct VtkBCFieldCache3D {
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

    VtkBCFieldCache3D BuildVtkBCFieldCache3D(const FieldManager_3D& fieldMgr) {
        const PhysicalProperties_string_op::Rock rock;
        const PhysicalProperties_string_op::EffectiveProps eff;
        const PhysicalProperties_string_op::Water water;
        const PhysicalProperties_string_op::CO2 gas;

        VtkBCFieldCache3D cache;
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

    double ComputeKn3D(const VtkBCFieldCache3D& cache, int owner, const Vector& normal, double floor_v) {
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

    double ComputeWaterMassMobility3D(const VtkBCFieldCache3D& cache, int owner, double floor_v) {
        const double rho_w = std::max(SafeFieldValue(cache.rho_w, owner, 1000.0), floor_v);
        const double lamw_mob = SafeFieldValue(cache.lamw_mob, owner, std::numeric_limits<double>::quiet_NaN());
        if (std::isfinite(lamw_mob) && lamw_mob > 0.0) {
            return std::max(rho_w * lamw_mob, floor_v);
        }
        const double mu_w = std::max(SafeFieldValue(cache.mu_w, owner, 1.0e-3), floor_v);
        return std::max(rho_w / mu_w, floor_v);
    }

    double ComputeGasMassMobility3D(const VtkBCFieldCache3D& cache, int owner, double floor_v) {
        const double rho_g = std::max(SafeFieldValue(cache.rho_g, owner, 1.0), floor_v);
        const double lamg_mob = SafeFieldValue(cache.lamg_mob, owner, std::numeric_limits<double>::quiet_NaN());
        if (std::isfinite(lamg_mob) && lamg_mob > 0.0) {
            return std::max(rho_g * lamg_mob, floor_v);
        }
        const double mu_g = std::max(SafeFieldValue(cache.mu_g, owner, 1.0e-5), floor_v);
        return std::max(rho_g / mu_g, floor_v);
    }

    double ComputeTransportCoeff3D(const VTKBCVariableBinding& binding,
                                   const VtkBCFieldCache3D& cache,
                                   int owner,
                                   const Vector& normal,
                                   double floor_v) {
        const double kn = ComputeKn3D(cache, owner, normal, floor_v);
        if (binding.transport_kind == VTKBCTransportKind::Temperature) {
            const double lambda_eff = SafeFieldValue(cache.lambda_eff, owner, std::numeric_limits<double>::quiet_NaN());
            if (std::isfinite(lambda_eff) && lambda_eff > floor_v) return lambda_eff;
            return std::max(SafeFieldValue(cache.lambda_r, owner, 2.0), floor_v);
        }
        if (binding.transport_kind == VTKBCTransportKind::Saturation) {
            return std::max(kn * ComputeGasMassMobility3D(cache, owner, floor_v), floor_v);
        }

        const bool use_gas_pressure = (binding.field_name == "p_g");
        const double mass_mob = use_gas_pressure
            ? ComputeGasMassMobility3D(cache, owner, floor_v)
            : ComputeWaterMassMobility3D(cache, owner, floor_v);
        return std::max(kn * mass_mob, floor_v);
    }

    void ReconstructBCPointField3D(const MeshManager_3D& meshMgr,
                                   const FieldManager_3D& fieldMgr,
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
        const VtkBCFieldCache3D field_cache = BuildVtkBCFieldCache3D(fieldMgr);

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

            const double coeff = ComputeTransportCoeff3D(binding, field_cache, owner, face.normal, ctx.coeff_floor);
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
            std::cout << "[PostProcess_3D][BC-Reconstruct] field=" << binding.field_name
                << " singular_denom=" << near_singular_count
                << " invalid_phi=" << invalid_value_count
                << " area_fallback=" << area_fallback_count << "\n";
        }
    }

    void ProjectBoundaryOwnerField3D(const MeshManager_3D& meshMgr,
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

    bool BuildDerivedFluidBCPointField3D(
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
        const bool usePrimaryForWFamily = ctx.usePrimaryForWaterFamily();
        const bool primaryIsCO2 = (ctx.primary_fluid_model == VTKBCPrimaryFluidModel::CO2);

        for (size_t i = 0; i < outVals.size(); ++i) {
            if (i >= boundaryMask.size() || boundaryMask[i] == 0) continue;

            const double Pw = hasPressure() ? (*pBC)[i] : ((i < pBase.size()) ? pBase[i] : 0.0);
            const double T = hasTemp() ? (*tBC)[i] : ((i < tBase.size()) ? tBase[i] : 300.0);
            const double Sw_raw = hasSw() ? (*swBC)[i] : ((i < swBase.size()) ? swBase[i] : 1.0);
            const double Sw = std::max(0.0, std::min(1.0, Sw_raw));

            ADVar<1> Pw_ad(Pw); Pw_ad.grad(0) = 1.0;
            ADVar<1> T_ad(T);
            // Keep legacy *_w_* output names for compatibility; physical source is policy-driven.
            const auto wProps = (usePrimaryForWFamily && primaryIsCO2)
                ? AD_Fluid::Evaluator::evaluateCO2<1>(Pw_ad, T_ad)
                : AD_Fluid::Evaluator::evaluateWater<1>(Pw_ad, T_ad);

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
std::vector<std::string> PostProcess_3D::GetAllUniqueFieldNames() const
{
    std::set<std::string> uniqueNames;
    for (const auto& pair : fieldMgr_.matrixFields.fields) {
        if (std::dynamic_pointer_cast<volScalarField>(pair.second)) uniqueNames.insert(pair.first);
    }
    for (const auto& pair : fieldMgr_.fractureFields.fields) {
        if (std::dynamic_pointer_cast<volScalarField>(pair.second)) uniqueNames.insert(pair.first);
    }
    return std::vector<std::string>(uniqueNames.begin(), uniqueNames.end());
}

std::string PostProcess_3D::GetTecplotElementType(int numNodes) const {
    switch (numNodes) {
    case 4: return "FETETRA";
    case 5: return "FEPYRAMID";
    case 6: return "FEPRISM";
    case 8: return "FEBRICK";
    default: return "FETETRA";
    }
}

// =========================================================
// 核心导出逻辑
// =========================================================
void PostProcess_3D::ExportTecplot(const std::string& filename, double time)
{
    std::cout << "[PostProcess] Exporting Ultra-Robust Tecplot file: " << filename << " (t=" << time << ")" << std::endl;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << std::endl;
        return;
    }

    out << std::scientific << std::setprecision(6);
    std::vector<std::string> allVarNames = GetAllUniqueFieldNames();

    const auto& nodes = meshMgr_.mesh().getNodes();
    std::unordered_map<int, int> nodeID2Index;
    for (size_t i = 0; i < nodes.size(); ++i) {
        nodeID2Index[nodes[i].id] = static_cast<int>(i) + 1;
    }

    // Header
    out << "TITLE = \"3D-EDFM Simulation Results\"\n";
    out << "VARIABLES = \"X\", \"Y\", \"Z\"";
    for (const auto& name : allVarNames) out << ", \"" << name << "\"";
    out << "\n";

    // ---------------------------------------------------------
    // Zone 1: Matrix
    // ---------------------------------------------------------
    const auto& cells = meshMgr_.mesh().getCells();
    size_t numCells = cells.size();

    if (numCells > 0) {
        std::string internalDesc = meshMgr_.mesh().getMatrixElementType();
        if (internalDesc == "Unknown" || internalDesc.empty()) {
            internalDesc = GetTecplotElementType(static_cast<int>(cells[0].CellNodeIDs.size()));
        }
        TecplotElemInfo tecInfo = GetTecplotInfo(internalDesc);

        // 强制写入 FEBRICK 类型
        out << "ZONE T=\"Matrix\"\n";
        out << " STRANDID=1, SOLUTIONTIME=" << time << "\n";
        out << " NODES=" << nodes.size() << ", ELEMENTS=" << numCells
            << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n";

        if (!allVarNames.empty()) {
            out << " VARLOCATION=([4-" << (3 + allVarNames.size()) << "]=CELLCENTERED)\n";
        }

        int col = 0;
        // Geometry Output 
        col = 0; for (const auto& node : nodes) { out << node.coord.m_x << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        col = 0; for (const auto& node : nodes) { out << node.coord.m_y << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        col = 0; for (const auto& node : nodes) { out << node.coord.m_z << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";

        // Fields Output (带 NaN 洗钱)
        for (const auto& name : allVarNames) {
            auto fieldPtr = fieldMgr_.getMatrixScalar(name);
            col = 0;
            if (fieldPtr) {
                for (double val : fieldPtr->data) {
                    out << (std::isfinite(val) ? val : 0.0) << " ";
                    if (++col % 10 == 0) out << "\n";
                }
            }
            else {
                for (size_t i = 0; i < numCells; ++i) {
                    out << "0.0 ";
                    if (++col % 10 == 0) out << "\n";
                }
            }
            if (col % 10 != 0) out << "\n";
        }

        // =========================================================
        // Connectivity (Mapping + Padding with Chirality Auto-Correction)
        // =========================================================
        for (const auto& cell : cells) {
            const auto& ids = cell.CellNodeIDs;
            int n = static_cast<int>(ids.size());

            std::vector<int> mappedIndices;
            mappedIndices.reserve(n);
            for (int id : ids) {
                if (nodeID2Index.find(id) != nodeID2Index.end()) mappedIndices.push_back(nodeID2Index[id]);
                else mappedIndices.push_back(1);
            }

            // --- 【核心引擎】：三维单元体积手性自动校验与修复 ---
            if (n >= 4) {
                int idx0 = mappedIndices[0] - 1;
                int idx1 = mappedIndices[1] - 1;
                int idx2 = mappedIndices[2] - 1;
                int idx3 = -1;

                if (n == 8) idx3 = mappedIndices[4] - 1;
                else if (n == 6) idx3 = mappedIndices[3] - 1;
                else if (n == 5) idx3 = mappedIndices[4] - 1;
                else if (n == 4) idx3 = mappedIndices[3] - 1;

                if (idx3 >= 0 && idx3 < nodes.size()) {
                    double ax = nodes[idx1].coord.m_x - nodes[idx0].coord.m_x;
                    double ay = nodes[idx1].coord.m_y - nodes[idx0].coord.m_y;
                    double az = nodes[idx1].coord.m_z - nodes[idx0].coord.m_z;

                    double bx = nodes[idx2].coord.m_x - nodes[idx0].coord.m_x;
                    double by = nodes[idx2].coord.m_y - nodes[idx0].coord.m_y;
                    double bz = nodes[idx2].coord.m_z - nodes[idx0].coord.m_z;

                    double cx = nodes[idx3].coord.m_x - nodes[idx0].coord.m_x;
                    double cy = nodes[idx3].coord.m_y - nodes[idx0].coord.m_y;
                    double cz = nodes[idx3].coord.m_z - nodes[idx0].coord.m_z;

                    double nx = ay * bz - az * by;
                    double ny = az * bx - ax * bz;
                    double nz = ax * by - ay * bx;

                    double vol = nx * cx + ny * cy + nz * cz;

                    // 如果体积小于0，强行翻转手性
                    if (vol < 0.0) {
                        if (n == 8) {
                            std::swap(mappedIndices[0], mappedIndices[4]);
                            std::swap(mappedIndices[1], mappedIndices[5]);
                            std::swap(mappedIndices[2], mappedIndices[6]);
                            std::swap(mappedIndices[3], mappedIndices[7]);
                        }
                        else if (n == 6) {
                            std::swap(mappedIndices[1], mappedIndices[2]);
                            std::swap(mappedIndices[4], mappedIndices[5]);
                        }
                        else if (n == 5) {
                            std::swap(mappedIndices[1], mappedIndices[3]);
                        }
                        else if (n == 4) {
                            std::swap(mappedIndices[1], mappedIndices[2]);
                        }
                    }
                }
            }

            // 执行退化映射
            if (n == 8) {
                for (int k = 0; k < 8; ++k) out << mappedIndices[k] << " ";
            }
            else if (n == 6) { // 三棱柱 6 -> 8
                out << mappedIndices[0] << " " << mappedIndices[1] << " " << mappedIndices[2] << " " << mappedIndices[2] << " "
                    << mappedIndices[3] << " " << mappedIndices[4] << " " << mappedIndices[5] << " " << mappedIndices[5] << " ";
            }
            else if (n == 5) { // 金字塔 5 -> 8
                out << mappedIndices[0] << " " << mappedIndices[1] << " " << mappedIndices[2] << " " << mappedIndices[3] << " "
                    << mappedIndices[4] << " " << mappedIndices[4] << " " << mappedIndices[4] << " " << mappedIndices[4] << " ";
            }
            else if (n == 4) { // 四面体 4 -> 8
                out << mappedIndices[0] << " " << mappedIndices[1] << " " << mappedIndices[2] << " " << mappedIndices[2] << " "
                    << mappedIndices[3] << " " << mappedIndices[3] << " " << mappedIndices[3] << " " << mappedIndices[3] << " ";
            }
            else {
                for (int k = 0; k < 8; ++k) out << mappedIndices[std::min(k, n - 1)] << " ";
            }
            out << "\n";
        }
    }

    // ---------------------------------------------------------
    // Zone 2: Fracture
    // ---------------------------------------------------------
    const auto& fractures = meshMgr_.fracture_network().getFractures();
    size_t totalFracElems = 0;
    size_t totalFracNodes = 0;

    for (const auto& frac : fractures) {
        totalFracNodes += frac.fracNodes.size();
        totalFracElems += frac.fracCells.size();
    }

    if (totalFracElems > 0) {
        // 强制写入 FEQUADRILATERAL 格式
        out << "ZONE T=\"Fracture\"\n";
        out << " STRANDID=2, SOLUTIONTIME=" << time << "\n";
        out << " NODES=" << totalFracNodes << ", ELEMENTS=" << totalFracElems
            << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n";

        if (!allVarNames.empty()) out << " VARLOCATION=([4-" << (3 + allVarNames.size()) << "]=CELLCENTERED)\n";

        int col = 0;
        // Fracture Geometry 
        col = 0; for (const auto& frac : fractures) for (const auto& n : frac.fracNodes) { out << n.coord.m_x << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        col = 0; for (const auto& frac : fractures) for (const auto& n : frac.fracNodes) { out << n.coord.m_y << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";
        col = 0; for (const auto& frac : fractures) for (const auto& n : frac.fracNodes) { out << n.coord.m_z << " "; if (++col % 10 == 0) out << "\n"; } if (col % 10 != 0) out << "\n";

        // Fracture Fields 
        for (const auto& name : allVarNames) {
            auto fieldPtr = fieldMgr_.getFractureScalar(name);
            col = 0;
            if (fieldPtr) {
                for (double val : fieldPtr->data) {
                    out << (std::isfinite(val) ? val : 0.0) << " ";
                    if (++col % 10 == 0) out << "\n";
                }
            }
            else {
                for (size_t i = 0; i < totalFracElems; ++i) {
                    out << "0.0 ";
                    if (++col % 10 == 0) out << "\n";
                }
            }
            if (col % 10 != 0) out << "\n";
        }

        // Connectivity
        size_t nodeOffset = 0;
        for (const auto& frac : fractures) {
            for (const auto& cell : frac.fracCells) {
                const auto& ids = cell.nodeIndices;
                int n = static_cast<int>(ids.size());

                for (int k = 0; k < 4; ++k) {
                    out << (nodeOffset + ids[std::min(k, n - 1)] + 1) << " ";
                }
                out << "\n";
            }
            nodeOffset += frac.fracNodes.size();
        }
    }

    out.close();
    std::cout << "[PostProcess] Export completed.\n";
}

// =========================================================
// 新一代渲染引擎：全量导出至 ParaView (VTK 格式)
// =========================================================
void PostProcess_3D::ExportVTK(const std::string& filename, double time)
{
    std::cout << "[PostProcess] Exporting ParaView VTK file: " << filename << " (t=" << time << ")" << std::endl;

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[Error] Cannot open file: " << filename << std::endl;
        return;
    }

    out << std::scientific << std::setprecision(6);
    std::vector<std::string> allVarNames = GetAllUniqueFieldNames();

    // 1. 提取基岩数据
    const auto& nodes = meshMgr_.mesh().getNodes();
    const auto& cells = meshMgr_.mesh().getCells();
    size_t numMatrixNodes = nodes.size();
    size_t numMatrixCells = cells.size();

    // 建立基岩节点到 0-based 索引的映射 (VTK 必须从 0 开始)
    std::unordered_map<int, int> nodeID2Index;
    for (size_t i = 0; i < numMatrixNodes; ++i) {
        nodeID2Index[nodes[i].id] = static_cast<int>(i);
    }

    // 2. 提取裂缝数据
    const auto& fractures = meshMgr_.fracture_network().getFractures();
    size_t numFracNodes = 0;
    size_t numFracCells = 0;
    for (const auto& frac : fractures) {
        numFracNodes += frac.fracNodes.size();
        numFracCells += frac.fracCells.size();
    }

    // 总统计
    size_t totalNodes = numMatrixNodes + numFracNodes;
    size_t totalCells = numMatrixCells + numFracCells;

    // 计算 CELLS 列表的总长度 (每个单元记录：节点数 n + n个索引)
    size_t cellListSize = 0;
    for (const auto& cell : cells) cellListSize += (1 + cell.CellNodeIDs.size());
    for (const auto& frac : fractures) {
        for (const auto& cell : frac.fracCells) {
            cellListSize += (1 + cell.nodeIndices.size());
        }
    }

    // ==========================================
    // VTK Header
    // ==========================================
    out << "# vtk DataFile Version 3.0\n";
    out << "3D-EDFM Polyhedron Simulation (Time: " << time << ")\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // ==========================================
    // POINTS (节点坐标)
    // ==========================================
    out << "POINTS " << totalNodes << " double\n";
    int col = 0;
    // 基岩节点
    for (const auto& n : nodes) {
        out << n.coord.m_x << " " << n.coord.m_y << " " << n.coord.m_z << " ";
        if (++col % 3 == 0) out << "\n";
    }
    // 裂缝节点
    for (const auto& frac : fractures) {
        for (const auto& n : frac.fracNodes) {
            out << n.coord.m_x << " " << n.coord.m_y << " " << n.coord.m_z << " ";
            if (++col % 3 == 0) out << "\n";
        }
    }
    if (col % 3 != 0) out << "\n";

    // ==========================================
    // CELLS (单元连接表)
    // ==========================================
    out << "CELLS " << totalCells << " " << cellListSize << "\n";

    // 基岩单元
    for (const auto& cell : cells) {
        out << cell.CellNodeIDs.size() << " ";
        for (int id : cell.CellNodeIDs) {
            out << nodeID2Index[id] << " ";
        }
        out << "\n";
    }

    // 裂缝单元 (节点索引需加上基岩节点总数作为偏移量)
    size_t fracNodeOffset = numMatrixNodes;
    for (const auto& frac : fractures) {
        for (const auto& cell : frac.fracCells) {
            out << cell.nodeIndices.size() << " ";
            for (int localIdx : cell.nodeIndices) {
                out << (fracNodeOffset + localIdx) << " ";
            }
            out << "\n";
        }
        fracNodeOffset += frac.fracNodes.size();
    }

    // ==========================================
    // CELL_TYPES (VTK 单元类型定义)
    // ==========================================
    // 10=Tetra, 12=Hex, 13=Wedge(Prism), 14=Pyramid, 5=Triangle, 9=Quad
    out << "CELL_TYPES " << totalCells << "\n";

    // 基岩类型
    for (const auto& cell : cells) {
        size_t n = cell.CellNodeIDs.size();
        if (n == 8) out << "12\n";
        else if (n == 6) out << "13\n"; // 原生支持三棱柱，完美！
        else if (n == 5) out << "14\n";
        else if (n == 4) out << "10\n";
        else out << "10\n";
    }

    // 裂缝类型
    for (const auto& frac : fractures) {
        for (const auto& cell : frac.fracCells) {
            if (cell.nodeIndices.size() == 4) out << "9\n"; // Quad
            else out << "5\n"; // Triangle
        }
    }

    // ==========================================
    // CELL_DATA (物理场输出)
    // ==========================================
    out << "CELL_DATA " << totalCells << "\n";

    out << "SCALARS DomainID int 1\n";
    out << "LOOKUP_TABLE default\n";
    // 基岩单元全部标记为 0
    for (size_t i = 0; i < numMatrixCells; ++i) out << "0\n";
    // 裂缝单元全部标记为 1
    for (size_t i = 0; i < numFracCells; ++i) out << "1\n";

    for (const auto& name : allVarNames) {
        out << "SCALARS " << name << " double 1\n";
        out << "LOOKUP_TABLE default\n";

        // 输出基岩场
        auto matField = fieldMgr_.getMatrixScalar(name);
        if (matField) {
            for (double val : matField->data) out << (std::isfinite(val) ? val : 0.0) << "\n";
        }
        else {
            for (size_t i = 0; i < numMatrixCells; ++i) out << "0.0\n";
        }

        // 输出裂缝场
        auto fracField = fieldMgr_.getFractureScalar(name);
        if (fracField) {
            for (double val : fracField->data) out << (std::isfinite(val) ? val : 0.0) << "\n";
        }
        else {
            for (size_t i = 0; i < numFracCells; ++i) out << "0.0\n";
        }
    }

    // ==========================================
    // POINT_DATA (node-averaged fields + BC-aware fields)
    // ==========================================
    out << "POINT_DATA " << totalNodes << "\n";

    std::unordered_map<std::string, std::vector<double>> pointFieldCache;
    pointFieldCache.reserve(allVarNames.size());

    for (const auto& name : allVarNames) {
        std::vector<double> pointVals(totalNodes, 0.0);

        // Matrix: cell-to-node arithmetic average
        std::vector<int> matCounts(numMatrixNodes, 0);
        auto matField = fieldMgr_.getMatrixScalar(name);
        if (matField) {
            for (size_t c = 0; c < numMatrixCells; ++c) {
                const double val = (c < matField->data.size() && std::isfinite(matField->data[c])) ? matField->data[c] : 0.0;
                for (int id : cells[c].CellNodeIDs) {
                    auto it = nodeID2Index.find(id);
                    if (it == nodeID2Index.end()) continue;
                    const int pIdx = it->second;
                    if (pIdx < 0 || pIdx >= static_cast<int>(numMatrixNodes)) continue;
                    pointVals[pIdx] += val;
                    matCounts[pIdx]++;
                }
            }
            for (size_t p = 0; p < numMatrixNodes; ++p) {
                if (matCounts[p] > 0) pointVals[p] /= matCounts[p];
            }
        }

        // Fracture: average over fracture cells sharing each fracture node
        auto fracField = fieldMgr_.getFractureScalar(name);
        if (fracField && numFracNodes > 0) {
            std::vector<double> fracSum(numFracNodes, 0.0);
            std::vector<int> fracCounts(numFracNodes, 0);
            size_t fracElemGlobal = 0;
            size_t fracNodeBase = 0;
            for (const auto& frac : fractures) {
                for (const auto& fcell : frac.fracCells) {
                    const double val = (fracElemGlobal < fracField->data.size() && std::isfinite(fracField->data[fracElemGlobal]))
                        ? fracField->data[fracElemGlobal]
                        : 0.0;
                    for (int localIdx : fcell.nodeIndices) {
                        const size_t gIdx = fracNodeBase + static_cast<size_t>(std::max(localIdx, 0));
                        if (gIdx < numFracNodes) {
                            fracSum[gIdx] += val;
                            fracCounts[gIdx]++;
                        }
                    }
                    ++fracElemGlobal;
                }
                fracNodeBase += frac.fracNodes.size();
            }
            for (size_t i = 0; i < numFracNodes; ++i) {
                const size_t pIdx = numMatrixNodes + i;
                if (fracCounts[i] > 0) pointVals[pIdx] = fracSum[i] / fracCounts[i];
            }
        }

        pointFieldCache.emplace(name, std::move(pointVals));
    }

    // Original point fields (backward-compatible)
    for (const auto& name : allVarNames) {
        out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
        const auto it = pointFieldCache.find(name);
        if (it != pointFieldCache.end()) {
            for (double v : it->second) out << v << "\n";
        }
        else {
            for (size_t p = 0; p < totalNodes; ++p) out << "0.0\n";
        }
    }

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
            ReconstructBCPointField3D(meshMgr_, fieldMgr_, matField, *bcVizCtx_, *binding, nodeID2Index, baseVals, bcVals, bcMask);
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

        for (const auto& name : allVarNames) {
            auto itBase = pointFieldCache.find(name);
            if (itBase == pointFieldCache.end()) continue;

            std::vector<double> baseMatVals(numMatrixNodes, 0.0);
            for (size_t i = 0; i < numMatrixNodes; ++i) baseMatVals[i] = itBase->second[i];

            std::vector<double> bcMatVals = baseMatVals;
            std::vector<int> bcMaskMat(numMatrixNodes, 0);

            const auto* binding = bcVizCtx_->findBinding(name);
            auto matField = fieldMgr_.getMatrixScalar(name);

            if (binding && matField) {
                ReconstructBCPointField3D(meshMgr_, fieldMgr_, matField, *bcVizCtx_, *binding, nodeID2Index, baseMatVals, bcMatVals, bcMaskMat);
            }
            else if (matField) {
                bcMaskMat = boundaryUnionMask;
                const bool derivedDone = BuildDerivedFluidBCPointField3D(
                    name, *bcVizCtx_,
                    pBaseMat, tBaseMat, swBaseMat,
                    hasP ? &pBcMat : nullptr,
                    hasT ? &tBcMat : nullptr,
                    hasSw ? &swBcMat : nullptr,
                    boundaryUnionMask,
                    bcMatVals);

                if (!derivedDone) {
                    ProjectBoundaryOwnerField3D(
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
    std::cout << "[PostProcess] VTK Export completed successfully.\n";
}
