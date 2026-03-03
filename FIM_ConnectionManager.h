/**
 * @file FIM_ConnectionManager.h
 * @brief 维度无关的纯代数连接聚合引擎 (支持 FF 多通道并联)
 */

#pragma once

#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <string>

enum class ConnectionType {
    Matrix_Matrix = 0,
    Fracture_Internal = 1,
    Matrix_Fracture = 2, ///< NNC (基岩-裂缝)
    Fracture_Fracture = 3  ///< FF (裂缝-裂缝)
};

struct Connection {
    int nodeI;
    int nodeJ;
    double T_Flow;
    double T_Heat;
    double aux_area;  ///< 辅助几何量：面积 (用作权重)
    double aux_dist;  ///< 辅助几何量：特征距离
    ConnectionType type;

    Connection(int i, int j, double tf, double th, double a, double d, ConnectionType t)
        : nodeI(std::min(i, j)), nodeJ(std::max(i, j)),
        T_Flow(tf), T_Heat(th), aux_area(a), aux_dist(d), type(t) {
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

    void PushConnection(int nodeI, int nodeJ, double t_flow, double t_heat, double aux_area, double aux_dist, ConnectionType type) {
        if (nodeI < 0 || nodeJ < 0 || nodeI == nodeJ) return;

        auto validateAndClip = [](double& val, const std::string& name) {
            if (std::isnan(val) || std::isinf(val))
                throw std::runtime_error("[FIM Error] " + name + " is NaN/Inf!");
            if (val < 0.0) {
                if (val > -1e-14) val = 0.0; // 数值噪声安全裁剪
                else throw std::runtime_error("[FIM Error] " + name + " is strictly negative!");
            }
            };

        validateAndClip(t_flow, "T_Flow");
        validateAndClip(t_heat, "T_Heat");
        validateAndClip(aux_area, "Aux_Area");
        validateAndClip(aux_dist, "Aux_Dist");

        rawBuffer_.emplace_back(nodeI, nodeJ, t_flow, t_heat, aux_area, aux_dist, type);
    }

    void FinalizeAndAggregate() {
        std::unordered_map<ConnectionKey, Connection, ConnectionKeyHash> aggMap;

        for (const auto& raw : rawBuffer_) {
            ConnectionKey key{ raw.type, raw.nodeI, raw.nodeJ };
            auto it = aggMap.find(key);

            if (it == aggMap.end()) {
                aggMap.insert({ key, raw });
            }
            else {
                // 【核心物理规则】：仅允许 NNC 与 FF 出现多通道并联，并进行物理传导率精准累加
                if (raw.type == ConnectionType::Matrix_Fracture || raw.type == ConnectionType::Fracture_Fracture) {
                    it->second.T_Flow += raw.T_Flow;
                    it->second.T_Heat += raw.T_Heat;
                    double totalArea = it->second.aux_area + raw.aux_area;
                    if (totalArea > 1e-14) {
                        it->second.aux_dist = (it->second.aux_dist * it->second.aux_area + raw.aux_dist * raw.aux_area) / totalArea;
                    }
                    it->second.aux_area = totalArea;
                }
                else {
                    // Matrix_Matrix 和 Fracture_Internal 在拓扑上绝对不应该出现重复连接
                    throw std::logic_error("[FIM Topology Error] Duplicate Matrix or FI connection detected at nodes (" +
                        std::to_string(raw.nodeI) + ", " + std::to_string(raw.nodeJ) + ")!");
                }
            }
        }

        globalConnections_.clear();
        globalConnections_.reserve(aggMap.size());
        for (const auto& kv : aggMap) globalConnections_.push_back(kv.second);
        rawBuffer_.clear();

        // 强确定性排序
        std::sort(globalConnections_.begin(), globalConnections_.end(),
            [](const Connection& a, const Connection& b) {
                if (a.type != b.type) return static_cast<int>(a.type) < static_cast<int>(b.type);
                if (a.nodeI != b.nodeI) return a.nodeI < b.nodeI;
                return a.nodeJ < b.nodeJ;
            });
    }

    const std::vector<Connection>& GetConnections() const { return globalConnections_; }

private:
    std::vector<Connection> rawBuffer_;
    std::vector<Connection> globalConnections_;
};