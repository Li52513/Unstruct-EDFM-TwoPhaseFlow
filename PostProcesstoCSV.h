#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif
#if __cplusplus >= 201703L
#include <filesystem>
#endif
#include <limits>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "WellConfig.h"
#include "FVM_WellCoupling.h"

/**
 * @brief Output helpers for SinglePhase (HT) CSV snapshots.
 */
namespace SinglePhase
{
    namespace OutputCSV
    {
        inline std::string parent_path_of(const std::string& path)
        {
            const size_t pos = path.find_last_of("/\\");
            if (pos == std::string::npos) return {};
            return path.substr(0, pos);
        }

        inline void ensureDirectory(const std::string& dir)
        {
            if (dir.empty()) return;
#if __cplusplus >= 201703L
            try
            {
                std::error_code ec;
                std::filesystem::create_directories(dir, ec);
            }
            catch (...) {}
#else
#if defined(_WIN32)
            auto mkdir_impl = [](const std::string& d) { _mkdir(d.c_str()); };
#else
            auto mkdir_impl = [](const std::string& d) { mkdir(d.c_str(), 0755); };
#endif
            std::string partial;
            partial.reserve(dir.size());
            for (char ch : dir)
            {
                if (ch == '/' || ch == '\\')
                {
                    if (!partial.empty())
                    {
                        mkdir_impl(partial);
                    }
                }
                partial.push_back(ch);
            }
            if (!partial.empty())
            {
                if (partial.back() == '/' || partial.back() == '\\') partial.pop_back();
                if (!partial.empty()) mkdir_impl(partial);
            }
#endif
        }

        inline void ensureParentDirectory(const std::string& filePath)
        {
            ensureDirectory(parent_path_of(filePath));
        }

        /// @brief Build a step-stamped csv filename: <prefix>_<tag>_step_00010.csv
        inline std::string makeStepFilenameCSV(
            const std::string& prefix, const std::string& tag, int step)
        {
            std::ostringstream fn;
            fn << prefix << "_" << tag << "_step_"
                << std::setw(5) << std::setfill('0') << step << ".csv";
            return fn.str();
        }

        /**
         * @brief Write one snapshot CSV containing cell id/center and (p,T).
         *
         * CSV columns:
         *  - cell_id, x, y, z, p, T
         *
         * @param prefix   Path prefix (can include directory), e.g. "output/csv/state"
         * @param step     Time-step index
         * @param simTime  Physical time (written as a comment line)
         * @param mgr      MeshManager
         * @param reg      FieldRegistry
         * @param p_cell   Cell pressure field name (e.g. "p")
         * @param T_cell   Cell temperature field name (e.g. "T")
         */
        inline bool writeCellPTSnapshotCSV(
            const std::string& prefix,
            int step,
            double simTime,
            MeshManager& mgr,
            FieldRegistry& reg,
            const std::string& p_cell,
            const std::string& T_cell)
        {
            if (prefix.empty()) return true;

            auto p = reg.get<volScalarField>(p_cell.c_str());
            auto T = reg.get<volScalarField>(T_cell.c_str());
            if (!p || !T) {
                std::cerr << "[SinglePhase][CSV] missing fields for cell snapshot: "
                    << p_cell << " / " << T_cell << "\n";
                return false;
            }

            const std::string filename = makeStepFilenameCSV(prefix, "cells", step);
            ensureParentDirectory(filename);

            std::ofstream ofs(filename);
            if (!ofs) {
                std::cerr << "[SinglePhase][CSV] cannot write '" << filename << "'.\n";
                return false;
            }

            Mesh& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const int Nc = static_cast<int>(cells.size());

            ofs << "# time=" << std::setprecision(16) << simTime << "\n";
            ofs << "cell_id,x,y,z," << p_cell << "," << T_cell << "\n";

            ofs << std::setprecision(16);
            for (int cidx = 0; cidx < Nc; ++cidx)
            {
                const auto& cell = cells[cidx];
                if (cell.id < 0) continue;

                ofs << cell.id << ","
                    << cell.center.m_x << ","
                    << cell.center.m_y << ","
                    << cell.center.m_z << ","
                    << (*p)[cidx] << ","
                    << (*T)[cidx] << "\n";
            }
            return true;
        }
        /// @brief Strip leading "p_w_" if WellDOF::name is actually the pw-field name.
        inline std::string wellDisplayName(const WellDOF& w)
        {
            const std::string pref = "p_w_";
            if (w.name.rfind(pref, 0) == 0 && w.name.size() > pref.size())
                return w.name.substr(pref.size());
            return w.name;
        }

        inline const char* toStr(WellDOF::Role r)
        {
            return (r == WellDOF::Role::Injector) ? "INJ" : "PROD";
        }
        inline const char* toStr(WellDOF::Mode m)
        {
            return (m == WellDOF::Mode::Pressure) ? "BHP" : "RATE";
        }

        /**
         * @brief Write one snapshot CSV containing per-well pw/Tbh/Qm.
         *
         * We reuse the same sign convention as your well heat source:
         *   qraw = PI_i * (pw - p_i)
         *   Qin  = sum(max(qraw,0)), Qout = sum(max(-qraw,0)), Qnet = sum(qraw).
         *
         * CSV columns:
         *  - well, role, mode, pw, T_bh, Qm_net, Qm_in, Qm_out, perfs
         *
         * @param prefix   Path prefix, e.g. "output/csv/state"
         * @param step     Time-step index
         * @param simTime  Physical time
         * @param mgr      MeshManager
         * @param reg      FieldRegistry
         * @param wells    WellDOF list (with per-well mask_field/PI_field/pw field name)
         * @param p_cell   Cell pressure field name
         * @param T_cell   Cell temperature field name
         */
        inline bool writeWellBHP_T_Qm_SnapshotCSV(
            const std::string& prefix,
            int step,
            double simTime,
            MeshManager& mgr,
            FieldRegistry& reg,
            const std::vector<WellDOF>& wells,
            const std::string& p_cell,
            const std::string& T_cell)
        {
            if (prefix.empty()) return true;

            const std::string filename = makeStepFilenameCSV(prefix, "wells", step);
            ensureParentDirectory(filename);

            std::ofstream ofs(filename);
            if (!ofs) {
                std::cerr << "[SinglePhase][CSV] cannot write '" << filename << "'.\n";
                return false;
            }

            auto p = reg.get<volScalarField>(p_cell.c_str());
            auto T = reg.get<volScalarField>(T_cell.c_str());
            if (!p || !T) {
                std::cerr << "[SinglePhase][CSV] missing fields for well snapshot: "
                    << p_cell << " / " << T_cell << "\n";
                return false;
            }

            Mesh& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const int Nc = static_cast<int>(cells.size());

            ofs << "# time=" << std::setprecision(16) << simTime << "\n";
            ofs << "well,role,mode,pw,T_bh,Qm_net,Qm_in,Qm_out,perfs\n";
            ofs << std::setprecision(16);

            // Even if wells is empty, we still generate a valid header-only snapshot file.
            for (const auto& w : wells)
            {
                // pw is stored in a single-value volScalarField named by w.name (= pw_name).:contentReference[oaicite:6]{index=6}
                auto pwf = reg.get<volScalarField>(w.name.c_str());
                const double pw = (pwf && pwf->size > 0) ? (*pwf)[0] : w.target;

                auto mk = reg.get<volScalarField>(w.mask_field.c_str());
                auto PI = reg.get<volScalarField>(w.PI_field.c_str());

                double Qnet = 0.0, Qin = 0.0, Qout = 0.0;
                double Tprod_num = 0.0;
                double Tperf_sum = 0.0;
                int    nPerf = 0;

                if (mk && PI)
                {
                    for (int cidx = 0; cidx < Nc; ++cidx)
                    {
                        if ((*mk)[cidx] <= 0.5) continue;
                        const double PIi = (*PI)[cidx];
                        if (PIi <= 0.0) continue;

                        const double qraw = PIi * (pw - (*p)[cidx]);
                        const double q_in = (qraw > 0.0) ? qraw : 0.0;
                        const double q_out = (qraw < 0.0) ? -qraw : 0.0;

                        Qnet += qraw;
                        Qin += q_in;
                        Qout += q_out;

                        Tperf_sum += (*T)[cidx];
                        if (q_out > 0.0) Tprod_num += q_out * (*T)[cidx];
                        ++nPerf;
                    }
                }

                const double Tperf_avg = (nPerf > 0) ? (Tperf_sum / nPerf)
                    : std::numeric_limits<double>::quiet_NaN();

                double Tbh = Tperf_avg;
                if (w.role == WellDOF::Role::Injector) {
                    // injection: prefer Tin, fallback to perforation average
                    Tbh = (std::isfinite(w.Tin) && w.Tin > 0.0) ? w.Tin : Tperf_avg;
                }
                else {
                    // production: use outflow-weighted temp, fallback to perforation average
                    Tbh = (Qout > 0.0) ? (Tprod_num / Qout) : Tperf_avg;
                }

                ofs << wellDisplayName(w) << ","
                    << toStr(w.role) << ","
                    << toStr(w.mode) << ","
                    << pw << ","
                    << Tbh << ","
                    << Qnet << ","
                    << Qin << ","
                    << Qout << ","
                    << nPerf << "\n";
            }
            return true;
        }

    } // namespace OutputCSV
} // namespace SinglePhase
