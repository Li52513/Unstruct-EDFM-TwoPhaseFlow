#pragma once
#include <unordered_map>
#include <tuple>
#include "TemperatureBC.h"

// 用法与 PressureBCAdapter 完全对称：从 faceId 取“原始”(未面积化) a,b,c
struct TemperatureBCAdapter {
    const TemperatureBC::Registry* reg = nullptr;
    const std::unordered_map<int, TemperatureBC::BoundaryCoefficient>* tableBC = nullptr;
    const std::unordered_map<int, std::tuple<double, double, double>>* tableTuple = nullptr;

    TemperatureBCAdapter() = default;
    explicit TemperatureBCAdapter(const TemperatureBC::Registry& r) : reg(&r) {}
    explicit TemperatureBCAdapter(const std::unordered_map<int, TemperatureBC::BoundaryCoefficient>& t) : tableBC(&t) {}
    explicit TemperatureBCAdapter(const std::unordered_map<int, std::tuple<double, double, double>>& t) : tableTuple(&t) {}

    bool getABC(int faceId, double& a, double& b, double& c) const {
        if (reg) {
            if (const auto* bc = reg->find(faceId)) { a = bc->a; b = bc->b; c = bc->c; return true; }
        }
        if (tableBC) {
            auto it = tableBC->find(faceId);
            if (it != tableBC->end()) { a = it->second.a; b = it->second.b; c = it->second.c; return true; }
        }
        if (tableTuple) {
            auto it = tableTuple->find(faceId);
            if (it != tableTuple->end()) { std::tie(a, b, c) = it->second; return true; }
        }
        return false;
    }
};