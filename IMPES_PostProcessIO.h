#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
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

namespace IMPES
{
    namespace Output
    {
        namespace detail
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
        } // namespace detail

        struct MassFluxSummary
        {
            double minFlux = 0.0;
            double maxFlux = 0.0;
            double sumFlux = 0.0;
            size_t positives = 0;
            size_t negatives = 0;
            size_t zeros = 0;
            double boundaryInflow = 0.0;
            double boundaryOutflow = 0.0;
        };

        inline bool computeMassFluxSummary(
            MeshManager& mgr,
            FaceFieldRegistry& freg,
            const std::string& mf_name,
            MassFluxSummary& summary,
            double eps = 1e-18)
        {
            auto mfF = freg.get<faceScalarField>(mf_name.c_str());
            if (!mfF)
            {
                std::cerr << "[IMPES][Output] mass-flux field '" << mf_name << "' not found for summary.\n";
                return false;
            }

            const auto& faces = mgr.mesh().getFaces();
            summary.minFlux = std::numeric_limits<double>::infinity();
            summary.maxFlux = -std::numeric_limits<double>::infinity();
            summary.sumFlux = 0.0;
            summary.positives = summary.negatives = summary.zeros = 0;
            summary.boundaryInflow = summary.boundaryOutflow = 0.0;

            for (const auto& F : faces)
            {
                const double flux = (*mfF)[F.id - 1];
                summary.minFlux = std::min(summary.minFlux, flux);
                summary.maxFlux = std::max(summary.maxFlux, flux);
                summary.sumFlux += flux;

                if (flux > eps) ++summary.positives;
                else if (flux < -eps) ++summary.negatives;
                else ++summary.zeros;

                if (F.isBoundary())
                {
                    if (flux > eps) summary.boundaryOutflow += flux;
                    if (flux < -eps) summary.boundaryInflow += -flux;
                }
            }
            if (!std::isfinite(summary.minFlux)) summary.minFlux = 0.0;
            if (!std::isfinite(summary.maxFlux)) summary.maxFlux = 0.0;
            return true;
        }

        struct TimeSeriesRecord
        {
            int step = 0;
            double simTime = 0.0;
            int outerIterations = 0;
            double dp_inf = 0.0;
            double dT_inf = 0.0;
            double satCFL = 0.0;
            double satSuggestedDt = 0.0;
            int pressureLinIters = 0;
            int temperatureLinIters = 0;
            MassFluxSummary flux;
        };

        inline bool appendTimeSeriesRecord(
            const std::string& file,
            const TimeSeriesRecord& rec)
        {
            if (file.empty()) return true;

            detail::ensureParentDirectory(file);
            bool needHeader = true;
            {
                std::ifstream ifs(file, std::ios::binary);
                if (ifs)
                {
                    ifs.seekg(0, std::ios::end);
                    needHeader = (ifs.tellg() == 0);
                }
            }

            std::ofstream ofs(file, std::ios::app);
            if (!ofs)
            {
                std::cerr << "[IMPES][Output] cannot open time-series file '" << file << "'.\n";
                return false;
            }

            if (needHeader)
            {
                ofs << "step,time,outer_iters,dp_inf,dT_inf,sat_CFL,suggested_dt,"
                    << "p_lin_iters,T_lin_iters,"
                    << "flux_min,flux_max,flux_sum,flux_pos,flux_neg,flux_zero,"
                    << "bdry_in,bdry_out\n";
            }

            ofs << rec.step << ','
                << std::setprecision(12) << rec.simTime << ','
                << rec.outerIterations << ','
                << rec.dp_inf << ','
                << rec.dT_inf << ','
                << rec.satCFL << ','
                << rec.satSuggestedDt << ','
                << rec.pressureLinIters << ','
                << rec.temperatureLinIters << ','
                << rec.flux.minFlux << ','
                << rec.flux.maxFlux << ','
                << rec.flux.sumFlux << ','
                << rec.flux.positives << ','
                << rec.flux.negatives << ','
                << rec.flux.zeros << ','
                << rec.flux.boundaryInflow << ','
                << rec.flux.boundaryOutflow << '\n';

            return true;
        }

        inline bool writeFieldSnapshotCSV(
            const std::string& prefix,
            int step,
            double simTime,
            MeshManager& mgr,
            FieldRegistry& reg,
            const std::vector<std::string>& scalarFields)
        {
            if (prefix.empty() || scalarFields.empty()) return true;

            std::ostringstream fn;
            fn << prefix << "_step_" << std::setw(5) << std::setfill('0') << step << ".csv";
            const std::string filename = fn.str();
            detail::ensureParentDirectory(filename);

            std::ofstream ofs(filename);
            if (!ofs)
            {
                std::cerr << "[IMPES][Output] cannot write snapshot '" << filename << "'.\n";
                return false;
            }

            std::vector<std::shared_ptr<volScalarField>> fieldPtrs;
            fieldPtrs.reserve(scalarFields.size());
            for (const auto& name : scalarFields)
            {
                auto fld = reg.get<volScalarField>(name.c_str());
                if (!fld)
                {
                    std::cerr << "[IMPES][Output] missing scalar field '" << name
                              << "' for snapshot output.\n";
                    return false;
                }
                fieldPtrs.emplace_back(fld);
            }

            ofs << "# time=" << std::setprecision(12) << simTime << "\n";
            ofs << "cell_id,x,y,z";
            for (const auto& name : scalarFields) ofs << ',' << name;
            ofs << '\n';

            const auto& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();

            for (const auto& cell : cells)
            {
                if (cell.id < 0) continue;
                const size_t idx = id2idx.at(cell.id);
                ofs << cell.id << ','
                    << cell.center.m_x << ','
                    << cell.center.m_y << ','
                    << cell.center.m_z;
                for (const auto& fld : fieldPtrs)
                {
                    ofs << ',' << (*fld)[idx];
                }
                ofs << '\n';
            }
            return true;
        }
    } // namespace Output
} // namespace IMPES


namespace IMPES_reviesd
{
    namespace Output
    {
        namespace detail
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
        } // namespace detail

        struct MassFluxSummary
        {
            double minFlux = 0.0;
            double maxFlux = 0.0;
            double sumFlux = 0.0;
            size_t positives = 0;
            size_t negatives = 0;
            size_t zeros = 0;
            double boundaryInflow = 0.0;
            double boundaryOutflow = 0.0;
        };

        inline bool computeMassFluxSummary(
            MeshManager& mgr,
            FaceFieldRegistry& freg,
            const std::string& mf_name,
            MassFluxSummary& summary,
            double eps = 1e-18)
        {
            auto mfF = freg.get<faceScalarField>(mf_name.c_str());
            if (!mfF)
            {
                std::cerr << "[IMPES][Output] mass-flux field '" << mf_name << "' not found for summary.\n";
                return false;
            }

            const auto& faces = mgr.mesh().getFaces();
            summary.minFlux = std::numeric_limits<double>::infinity();
            summary.maxFlux = -std::numeric_limits<double>::infinity();
            summary.sumFlux = 0.0;
            summary.positives = summary.negatives = summary.zeros = 0;
            summary.boundaryInflow = summary.boundaryOutflow = 0.0;

            for (const auto& F : faces)
            {
                const double flux = (*mfF)[F.id - 1];
                summary.minFlux = std::min(summary.minFlux, flux);
                summary.maxFlux = std::max(summary.maxFlux, flux);
                summary.sumFlux += flux;

                if (flux > eps) ++summary.positives;
                else if (flux < -eps) ++summary.negatives;
                else ++summary.zeros;

                if (F.isBoundary())
                {
                    if (flux > eps) summary.boundaryOutflow += flux;
                    if (flux < -eps) summary.boundaryInflow += -flux;
                }
            }
            if (!std::isfinite(summary.minFlux)) summary.minFlux = 0.0;
            if (!std::isfinite(summary.maxFlux)) summary.maxFlux = 0.0;
            return true;
        }

        struct TimeSeriesRecord
        {
            int step = 0;
            double simTime = 0.0;
            int outerIterations = 0;
            double dp_inf = 0.0;
            double dT_inf = 0.0;
            double satCFL = 0.0;
            double satSuggestedDt = 0.0;
            int pressureLinIters = 0;
            int temperatureLinIters = 0;
            MassFluxSummary flux;
        };

        inline bool appendTimeSeriesRecord(
            const std::string& file,
            const TimeSeriesRecord& rec)
        {
            if (file.empty()) return true;

            detail::ensureParentDirectory(file);
            bool needHeader = true;
            {
                std::ifstream ifs(file, std::ios::binary);
                if (ifs)
                {
                    ifs.seekg(0, std::ios::end);
                    needHeader = (ifs.tellg() == 0);
                }
            }

            std::ofstream ofs(file, std::ios::app);
            if (!ofs)
            {
                std::cerr << "[IMPES][Output] cannot open time-series file '" << file << "'.\n";
                return false;
            }

            if (needHeader)
            {
                ofs << "step,time,outer_iters,dp_inf,dT_inf,sat_CFL,suggested_dt,"
                    << "p_lin_iters,T_lin_iters,"
                    << "flux_min,flux_max,flux_sum,flux_pos,flux_neg,flux_zero,"
                    << "bdry_in,bdry_out\n";
            }

            ofs << rec.step << ','
                << std::setprecision(12) << rec.simTime << ','
                << rec.outerIterations << ','
                << rec.dp_inf << ','
                << rec.dT_inf << ','
                << rec.satCFL << ','
                << rec.satSuggestedDt << ','
                << rec.pressureLinIters << ','
                << rec.temperatureLinIters << ','
                << rec.flux.minFlux << ','
                << rec.flux.maxFlux << ','
                << rec.flux.sumFlux << ','
                << rec.flux.positives << ','
                << rec.flux.negatives << ','
                << rec.flux.zeros << ','
                << rec.flux.boundaryInflow << ','
                << rec.flux.boundaryOutflow << '\n';

            return true;
        }

        inline bool writeFieldSnapshotCSV(
            const std::string& prefix,
            int step,
            double simTime,
            MeshManager& mgr,
            FieldRegistry& reg,
            const std::vector<std::string>& scalarFields)
        {
            if (prefix.empty() || scalarFields.empty()) return true;

            std::ostringstream fn;
            fn << prefix << "_step_" << std::setw(5) << std::setfill('0') << step << ".csv";
            const std::string filename = fn.str();
            detail::ensureParentDirectory(filename);

            std::ofstream ofs(filename);
            if (!ofs)
            {
                std::cerr << "[IMPES][Output] cannot write snapshot '" << filename << "'.\n";
                return false;
            }

            std::vector<std::shared_ptr<volScalarField>> fieldPtrs;
            fieldPtrs.reserve(scalarFields.size());
            for (const auto& name : scalarFields)
            {
                auto fld = reg.get<volScalarField>(name.c_str());
                if (!fld)
                {
                    std::cerr << "[IMPES][Output] missing scalar field '" << name
                        << "' for snapshot output.\n";
                    return false;
                }
                fieldPtrs.emplace_back(fld);
            }

            ofs << "# time=" << std::setprecision(12) << simTime << "\n";
            ofs << "cell_id,x,y,z";
            for (const auto& name : scalarFields) ofs << ',' << name;
            ofs << '\n';

            const auto& mesh = mgr.mesh();
            const auto& cells = mesh.getCells();
            const auto& id2idx = mesh.getCellId2Index();

            for (const auto& cell : cells)
            {
                if (cell.id < 0) continue;
                const size_t idx = id2idx.at(cell.id);
                ofs << cell.id << ','
                    << cell.center.m_x << ','
                    << cell.center.m_y << ','
                    << cell.center.m_z;
                for (const auto& fld : fieldPtrs)
                {
                    ofs << ',' << (*fld)[idx];
                }
                ofs << '\n';
            }
            return true;
        }
    } // namespace Output
} // namespace IMPES
