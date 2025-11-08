#include "FVM_SourceTerm_WellHeat.h"
#include "FieldAcessForDiscre.h"
#include <iostream>
#include <algorithm>
#include <cassert>

namespace FVM {
    namespace SourceTerm {

        static inline int cell_id_by_index(const Mesh& mesh, int cidx)
        {
            const auto& cells = mesh.getCells();
            assert(cidx >= 0 && cidx < (int)cells.size());
            return cells[cidx].id;
        }

        bool add_temperature_source_from_qm_single(
            MeshManager& mgr,
            FieldRegistry& reg,
            WellDOF::Role role,
            const char* PI_name,
            const char* mask_name,
            const char* pw_field,
            double Tin,
            const char* p_cell,
            const char* cp_field,
            double thickness,
            bool accumulate,
            bool verbose,
            const InjectorDirichletPinOptions& pin)
        {
            Mesh& mesh = mgr.mesh();
            const int Nc = (int)mesh.getCells().size();

            auto p = reg.get<volScalarField>(p_cell);
            auto cp = reg.get<volScalarField>(cp_field);
            auto PI = reg.get<volScalarField>(PI_name);
            auto mk = reg.get<volScalarField>(mask_name);
            auto pwf = reg.get<volScalarField>(pw_field);
            if (!p || !cp || !PI || !mk || !pwf) {
                std::cerr << "[WellHeatSingle] missing field(s) for "
                    << mask_name << "/" << PI_name << "\n";
                return false;
            }

            auto aT = reg.get<volScalarField>("a_src_T");
            auto bT = reg.get<volScalarField>("b_src_T");
            if (!aT) aT = reg.create<volScalarField>("a_src_T", Nc, 0.0);
            if (!bT) bT = reg.create<volScalarField>("b_src_T", Nc, 0.0);
            if (!accumulate) {
                for (int i = 0; i < Nc; ++i) { (*aT)[i] = 0.0; (*bT)[i] = 0.0; }
            }

            const double pw = (*pwf)[0];

            double Qm_in = 0.0, Qm_out = 0.0;
            double Qh_in = 0.0;
            int nPerf = 0;

            for (int cidx = 0; cidx < Nc; ++cidx) {
                if ((*mk)[cidx] <= 0.5) continue;

                const double PIi = (*PI)[cidx];
                if (PIi <= 0.0) continue;

                const double p_i = (*p)[cidx];
                const double cp_i = (*cp)[cidx];
                const double qraw = PIi * (pw - p_i);
                const double q_in = (qraw > 0.0) ? qraw : 0.0;
                const double q_out = (qraw < 0.0) ? -qraw : 0.0;

                const double diagHeat = cp_i * (q_in + q_out) * thickness;
                (*aT)[cidx] += diagHeat;

                if (q_in > 0.0) {
                    const double rhs = cp_i * q_in * Tin * thickness;
                    (*bT)[cidx] += rhs;
                    Qh_in += rhs;
                    Qm_in += q_in;

                    // --- new: Ç¿Ô¼Êø -------------------------------------------------
                    if (pin.enable && role == WellDOF::Role::Injector) {
                        double weight = 0.0;
                        if (pin.relativeWeight > 0.0 && Tin != 0.0)
                            weight = pin.relativeWeight * rhs / Tin;
                        weight = std::max(weight, pin.minWeight);
                        if (weight > 0.0) {
                            (*aT)[cidx] += weight;
                            (*bT)[cidx] += weight * Tin;
                        }
                    }
                    // -----------------------------------------------------------------
                }

                if (q_out > 0.0) Qm_out += q_out;
                ++nPerf;
            }

            if (verbose) {
                std::cout << "[WellHeatSingle][" << mask_name << "] role="
                    << (role == WellDOF::Role::Injector ? "INJ" : "PROD")
                    << ", perfs=" << nPerf
                    << ", Qm_in=" << Qm_in << " kg/s"
                    << ", Qm_out=" << Qm_out << " kg/s"
                    << ", Qh_in=" << Qh_in << " W\n";
            }
            return true;
        }

    } // namespace SourceTerm
} // namespace FVM