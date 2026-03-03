/**
 * @file Test_FIM_Topology.h
 * @brief FIM 全连接图守恒与纯净度测试
 */

#pragma once
#include <iostream>
#include <unordered_set>
#include <stdexcept>
#include <cmath>
#include "FIM_ConnectionManager.h"
#include "FIM_TopologyBuilder3D.h"
#include "FIM_TopologyBuilder2D.h"

inline void Benchmark_FIM_Topology_Pipeline(const MeshManager_3D& meshMgr, const FieldManager_3D& fieldMgr) {
    std::cout << "\n========== [Day 1] FIM Global Topology Ironclad Benchmark ==========" << std::endl;

    auto pNNC_Flow = fieldMgr.getNNCScalar("T_NNC_Flow");
    auto pNNC_Heat = fieldMgr.getNNCScalar("T_NNC_Heat");
    double raw_NNC_Flow_Sum = 0.0, raw_NNC_Heat_Sum = 0.0;

    if (pNNC_Flow && pNNC_Heat) {
        for (double f : pNNC_Flow->data) raw_NNC_Flow_Sum += f;
        for (double h : pNNC_Heat->data) raw_NNC_Heat_Sum += h;
    }

    FIM_ConnectionManager connMgr;
    FIM_TopologyBuilder3D::LoadAllConnections(connMgr, meshMgr, fieldMgr);
    connMgr.FinalizeAndAggregate();
    connMgr.PrintAggregationStats();
    const auto& conns = connMgr.GetConnections();

    std::unordered_set<ConnectionKey, ConnectionKeyHash> uniqueKeySet;
    double agg_NNC_Flow_Sum = 0.0, agg_NNC_Heat_Sum = 0.0;
    int countMM = 0, countFI = 0, countNNC = 0, countFF = 0;

    for (size_t i = 0; i < conns.size(); ++i) {
        const auto& c = conns[i];
        if (c.nodeI >= c.nodeJ) throw std::logic_error("[Fatal] Undirected edge not normalized.");

        ConnectionKey key{ c.type, c.nodeI, c.nodeJ };
        if (!uniqueKeySet.insert(key).second) {
            throw std::logic_error("[Fatal] Duplicate Mechanism detected in Final Output!");
        }

        if (i > 0) {
            const auto& prev = conns[i - 1];
            bool sorted = (prev.type < c.type) ||
                (prev.type == c.type && prev.nodeI < c.nodeI) ||
                (prev.type == c.type && prev.nodeI == c.nodeI && prev.nodeJ < c.nodeJ);
            if (!sorted) throw std::logic_error("[Fatal] Array is not deterministically sorted.");
        }

        switch (c.type) {
        case ConnectionType::Matrix_Matrix: countMM++; break;
        case ConnectionType::Fracture_Internal: countFI++; break;
        case ConnectionType::Fracture_Fracture: countFF++; break;
        case ConnectionType::Matrix_Fracture:
            countNNC++;
            agg_NNC_Flow_Sum += c.T_Flow;
            agg_NNC_Heat_Sum += c.T_Heat;
            break;
        }
    }

    if (std::abs(raw_NNC_Flow_Sum - agg_NNC_Flow_Sum) > 1e-8 || std::abs(raw_NNC_Heat_Sum - agg_NNC_Heat_Sum) > 1e-8) {
        throw std::logic_error("[Fatal] Transmissibility Conservation (Flow or Heat) Broken!");
    }

    std::cout << "[PASS] Multi-Mechanism Aggregation & Deduplication Successful." << std::endl;
    std::cout << "[PASS] Flow & Heat Transmissibility Conservation Validated." << std::endl;
    std::cout << "  |- Matrix-Matrix : " << countMM << " connections\n"
        << "  |- Frac-Internal : " << countFI << " connections\n"
        << "  |- Matrix-Frac(NNC): " << countNNC << " clustered connections\n"
        << "  |- Frac-Frac(FF) : " << countFF << " connections (Successfully injected)" << std::endl;
    std::cout << "====================================================================\n" << std::endl;
}

inline void Benchmark_FIM_Topology_Pipeline_2D(const MeshManager& meshMgr, const FieldManager_2D& fieldMgr) {
    std::cout << "\n========== [Day 1] FIM Global Topology Ironclad Benchmark (2D) ==========" << std::endl;

    auto pNNC_Flow = fieldMgr.getNNCScalar("T_NNC_Flow");
    auto pNNC_Heat = fieldMgr.getNNCScalar("T_NNC_Heat");
    double raw_NNC_Flow_Sum = 0.0, raw_NNC_Heat_Sum = 0.0;

    if (pNNC_Flow && pNNC_Heat) {
        for (double f : pNNC_Flow->data) raw_NNC_Flow_Sum += f;
        for (double h : pNNC_Heat->data) raw_NNC_Heat_Sum += h;
    }

    // 调用 2D 装载器
    FIM_ConnectionManager connMgr;
    FIM_TopologyBuilder2D::LoadAllConnections(connMgr, meshMgr, fieldMgr);
    connMgr.FinalizeAndAggregate();
    connMgr.PrintAggregationStats();
    const auto& conns = connMgr.GetConnections();

    std::unordered_set<ConnectionKey, ConnectionKeyHash> uniqueKeySet;
    double agg_NNC_Flow_Sum = 0.0, agg_NNC_Heat_Sum = 0.0;
    int countMM = 0, countNNC = 0, countFF = 0;

    for (size_t i = 0; i < conns.size(); ++i) {
        const auto& c = conns[i];
        if (c.nodeI >= c.nodeJ) throw std::logic_error("[Fatal] Undirected edge not normalized.");

        ConnectionKey key{ c.type, c.nodeI, c.nodeJ };
        if (!uniqueKeySet.insert(key).second) {
            throw std::logic_error("[Fatal] Duplicate Mechanism detected in Final Output!");
        }

        if (i > 0) {
            const auto& prev = conns[i - 1];
            bool sorted = (prev.type < c.type) ||
                (prev.type == c.type && prev.nodeI < c.nodeI) ||
                (prev.type == c.type && prev.nodeI == c.nodeI && prev.nodeJ < c.nodeJ);
            if (!sorted) throw std::logic_error("[Fatal] Array is not deterministically sorted.");
        }

        switch (c.type) {
        case ConnectionType::Matrix_Matrix: countMM++; break;
        case ConnectionType::Fracture_Fracture: countFF++; break;
        case ConnectionType::Matrix_Fracture:
            countNNC++;
            agg_NNC_Flow_Sum += c.T_Flow;
            agg_NNC_Heat_Sum += c.T_Heat;
            break;
        default: break;
        }
    }

    if (std::abs(raw_NNC_Flow_Sum - agg_NNC_Flow_Sum) > 1e-8 || std::abs(raw_NNC_Heat_Sum - agg_NNC_Heat_Sum) > 1e-8) {
        throw std::logic_error("[Fatal] 2D Transmissibility Conservation (Flow or Heat) Broken!");
    }

    std::cout << "[PASS] 2D Multi-Mechanism Aggregation & Deduplication Successful." << std::endl;
    std::cout << "[PASS] 2D Flow & Heat Transmissibility Conservation Validated." << std::endl;
    std::cout << "  |- Matrix-Matrix : " << countMM << " connections\n"
        << "  |- Matrix-Frac(NNC): " << countNNC << " clustered connections\n"
        << "  |- Frac-Frac(FF) : " << countFF << " connections" << std::endl;
    std::cout << "=========================================================================\n" << std::endl;
}