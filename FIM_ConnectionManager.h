#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <string>
#include <array>
#include "UserDefineVarType.h"  // for Vector (non-orthogonal correction vector)

enum class ConnectionType {
    Matrix_Matrix = 0,
    Fracture_Internal = 1,
    Matrix_Fracture = 2, 
    Fracture_Fracture = 3
};

struct Connection {
    int nodeI;
    int nodeJ;
    double T_Flow;
    double T_Heat;
    double aux_area;  
    double aux_dist;
    ConnectionType type;
    /// Non-orthogonal correction vector T = Aj - vectorE (only non-zero for MM connections).
    /// Always oriented from nodeI to nodeJ (negated if original i>j to maintain consistency).
    Vector vectorT;

    Connection(int i, int j, double tf, double th, double a, double d, ConnectionType t,
               Vector vT = Vector(0.0, 0.0, 0.0))
        : nodeI(std::min(i, j)), nodeJ(std::max(i, j)),
        T_Flow(tf), T_Heat(th), aux_area(a), aux_dist(d), type(t)
    {
        // When i and j are swapped (nodeI=min), the face direction reverses → negate vectorT
        vectorT = (i <= j) ? vT : Vector(-vT.m_x, -vT.m_y, -vT.m_z);
    }
};

struct ConnectionKey {
    ConnectionType type;
    int nodeI;
    int nodeJ;
    bool operator==(const ConnectionKey& o) const {
        return type == o.type && nodeI == o.nodeI && nodeJ == o.nodeJ;
    }
};

struct ConnectionKeyHash {
    std::size_t operator()(const ConnectionKey& k) const {
        return std::hash<int>()(static_cast<int>(k.type)) ^
            (std::hash<int>()(k.nodeI) << 1) ^
            (std::hash<int>()(k.nodeJ) << 2);
    }
};

class FIM_ConnectionManager {
public:
    FIM_ConnectionManager() = default;

    struct AggregationStats {
        size_t raw_count = 0;
        size_t unique_count = 0;
        size_t merged_count = 0;
        std::array<size_t, 4> merged_by_type{ 0, 0, 0, 0 };
        std::array<size_t, 4> physical_parallel_by_type{ 0, 0, 0, 0 };
        std::array<size_t, 4> geometric_duplicate_by_type{ 0, 0, 0, 0 };
    };

    const AggregationStats& GetAggregationStats() const { return stats_; }

    void PrintAggregationStats(std::ostream& os = std::cout) const {
        os << "[FIM Aggregation] raw=" << stats_.raw_count
            << ", unique=" << stats_.unique_count
            << ", merged=" << stats_.merged_count << "\n";
        for (int t = 0; t < 4; ++t) {
            os << "  - " << TypeName(static_cast<ConnectionType>(t))
                << ": merged=" << stats_.merged_by_type[t]
                << ", physical_parallel=" << stats_.physical_parallel_by_type[t]
                << ", geometric_duplicate=" << stats_.geometric_duplicate_by_type[t]
                << "\n";
        }
    }
    void ReserveRawConnections(size_t expected) {
        rawBuffer_.reserve(expected);
    }

    void PushConnection(int nodeI, int nodeJ, double t_flow, double t_heat, double aux_area, double aux_dist, ConnectionType type,
                        Vector vectorT = Vector(0.0, 0.0, 0.0)) {
        if (nodeI < 0 || nodeJ < 0 || nodeI == nodeJ) return;

        auto validateAndClip = [](double& val, const std::string& name) {
            if (std::isnan(val) || std::isinf(val))
                throw std::runtime_error("[FIM Error] " + name + " is NaN/Inf!");
            if (val < 0.0) {
                if (val > -1e-14) val = 0.0;
                else throw std::runtime_error("[FIM Error] " + name + " is strictly negative!");
            }
            };

        validateAndClip(t_flow, "T_Flow");
        validateAndClip(t_heat, "T_Heat");
        validateAndClip(aux_area, "Aux_Area");
        validateAndClip(aux_dist, "Aux_Dist");

        rawBuffer_.emplace_back(nodeI, nodeJ, t_flow, t_heat, aux_area, aux_dist, type, vectorT);

        if (std::isnan(vectorT.m_x) || std::isnan(vectorT.m_y) || std::isnan(vectorT.m_z)) 
        {
            throw std::runtime_error("[FIM Error] vectorT contains NaN!");
        }

    }

    void FinalizeAndAggregate() {
        stats_ = AggregationStats{};
        stats_.raw_count = rawBuffer_.size();

        std::unordered_map<ConnectionKey, Connection, ConnectionKeyHash> aggMap;

        struct Sig {
            double area;
            double dist;
            double tf;
            double th;
        };
        std::unordered_map<ConnectionKey, std::vector<Sig>, ConnectionKeyHash> signatures;

        // pair = (aux_area, aux_dist)

        for (const auto& raw : rawBuffer_) {
            ConnectionKey key{ raw.type, raw.nodeI, raw.nodeJ };
            auto it = aggMap.find(key);

            if (it == aggMap.end()) {
                aggMap.insert({ key, raw });
                signatures[key].push_back({ raw.aux_area, raw.aux_dist, raw.T_Flow, raw.T_Heat });
                continue;
            }

            const size_t tid = TypeToIdx(raw.type);

            bool sameGeom = false;
            bool sameGeomAndTrans = false;
            auto& sigs = signatures[key];

            for (const auto& s : sigs) {
                const bool g = NearlyEqual(s.area, raw.aux_area, 1e-6, 1e-10) &&
                    NearlyEqual(s.dist, raw.aux_dist, 1e-6, 1e-10);
                if (g) {
                    sameGeom = true;
                    const bool t = NearlyEqual(s.tf, raw.T_Flow, 1e-6, 1e-12) &&
                        NearlyEqual(s.th, raw.T_Heat, 1e-6, 1e-12);
                    if (t) {
                        sameGeomAndTrans = true;
                        break;
                    }
                }
            }

            if (sameGeomAndTrans) {
                stats_.geometric_duplicate_by_type[tid]++;
                continue;

            }
            else {
                stats_.physical_parallel_by_type[tid]++;
                sigs.push_back({ raw.aux_area, raw.aux_dist, raw.T_Flow, raw.T_Heat });
            }
            if (raw.type == ConnectionType::Matrix_Fracture || raw.type == ConnectionType::Fracture_Fracture) {
                it->second.T_Flow += raw.T_Flow;
                it->second.T_Heat += raw.T_Heat;
                double totalArea = it->second.aux_area + raw.aux_area;
                if (totalArea > 1e-14) {
                    it->second.aux_dist =
                        (it->second.aux_dist * it->second.aux_area + raw.aux_dist * raw.aux_area) / totalArea;
                }
                it->second.aux_area = totalArea;

                stats_.merged_count++;
                stats_.merged_by_type[tid]++;
            }
            else {
                throw std::logic_error(
                    "[FIM Topology Error] Duplicate Matrix/FI connection at (" +
                    std::to_string(raw.nodeI) + "," + std::to_string(raw.nodeJ) + ")"
                );
            }
        }

        globalConnections_.clear();
        globalConnections_.reserve(aggMap.size());
        for (const auto& kv : aggMap) globalConnections_.push_back(kv.second);
        rawBuffer_.clear();

        std::sort(globalConnections_.begin(), globalConnections_.end(),
            [](const Connection& a, const Connection& b) {
                if (a.type != b.type) return static_cast<int>(a.type) < static_cast<int>(b.type);
                if (a.nodeI != b.nodeI) return a.nodeI < b.nodeI;
                return a.nodeJ < b.nodeJ;
            });

        stats_.unique_count = globalConnections_.size();
    }

    const std::vector<Connection>& GetConnections() const { return globalConnections_; }

private:
    std::vector<Connection> rawBuffer_;
    std::vector<Connection> globalConnections_;

    AggregationStats stats_;

    static size_t TypeToIdx(ConnectionType t) {
        return static_cast<size_t>(t);
    }

    static const char* TypeName(ConnectionType t) {
        switch (t) {
        case ConnectionType::Matrix_Matrix: return "Matrix_Matrix";
        case ConnectionType::Fracture_Internal: return "Fracture_Internal";
        case ConnectionType::Matrix_Fracture: return "Matrix_Fracture";
        case ConnectionType::Fracture_Fracture: return "Fracture_Fracture";
        default: return "Unknown";
        }
    }

    static bool NearlyEqual(double a, double b, double rel = 1e-8, double abs = 1e-12) {
        return std::fabs(a - b) <= std::max(abs, rel * std::max(std::fabs(a), std::fabs(b)));
    }

};