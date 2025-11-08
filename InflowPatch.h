// InflowProviderUtils_Told_Conditional.h
#pragma once
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "BoundaryFaceClassify.h"
#include "BCAdapter.h"  // TemperatureBCAdapter / PressureBCAdapter
#include "TemperatureBCAdapter.h"

namespace BoundaryUtils {

    // —— 复用你已有的 PatchFaceLookup（若已定义，可删掉此重复声明） —— //
    struct PatchFaceLookup {
        std::unordered_map<std::string, std::unordered_set<int>> byName;
        inline bool isIn(const std::string& name, int fid) const {
            auto it = byName.find(name);
            return it != byName.end() && it->second.count(fid);
        }
        inline bool isInAny(const std::vector<std::string>& names, int fid) const {
            for (const auto& n : names) if (isIn(n, fid)) return true;
            return false;
        }
    };
    inline PatchFaceLookup makePatchFaceLookup(const BoundaryFaceClassify::FaceGroups& G) {
        auto mk = [](const std::vector<int>& v) {
            return std::unordered_set<int>(v.begin(), v.end());
            };
        PatchFaceLookup L;
        L.byName["x0"] = mk(G.x0);
        L.byName["xL"] = mk(G.xL);
        L.byName["y0"] = mk(G.y0);
        L.byName["yL"] = mk(G.yL);
        // 3D 如有：
        // L.byName["z0"] = mk(G.z0);
        // L.byName["zL"] = mk(G.zL);
        return L;
    }

    // —— 通用：判断“零通量（绝热/零法向导数）”的 ABC —— //
    template<class ABC>
    inline bool isZeroFluxGradBC(const ABC& bc, int faceId, double eps = 1e-30) {
        double a = 0, b = 0, c = 0;
        if (!bc.getABC(faceId, a, b, c)) return false;
        return (std::abs(a) <= eps && std::abs(b) > eps && std::abs(c) <= eps);
    }

    // 核心：仅当 (1) 面在 patches 内；(2) 温度为绝热；(3) 压力非零通量；(4) 此面对“温度对流”为入流(flux<0)
    // 时，提供 Tb = T_prev(owner)。否则返回 false。
    inline std::function<bool(int, double&)>
        makeConditionalPrevTimeValueProviderForPatches
        (
            MeshManager& mgr,
            const FieldRegistry& reg,
            FaceFieldRegistry& freg,
            const PatchFaceLookup& L,
            const std::vector<std::string>& patches,   // 例如 {"xL"}
            const std::string& fieldName_prev,         // 一般为 "T_prev"
            const TemperatureBCAdapter& Tbc,
            const PressureBCAdapter& Pbc,
            const std::string& flux_face_name,         // 用于温度对流的面通量，例："mf_g" 或 "Qf_g"
            double epsFlux = 1e-18
        ) 
    {
        Mesh& mesh = mgr.mesh();
        const auto& faces = mesh.getFaces();
        const auto& id2idx = mesh.getCellId2Index();

        auto Tprev = reg.get<volScalarField>(fieldName_prev.c_str());
        if (!Tprev) {
            std::cerr << "[Provider] Missing vol field '" << fieldName_prev
                << "'. Conditional T_old provider disabled.\n";
            return [](int, double&)->bool { return false; };
        }
        auto fluxF = freg.get<faceScalarField>(flux_face_name.c_str());
        if (!fluxF) {
            std::cerr << "[Provider] Missing face flux '" << flux_face_name
                << "'. Conditional T_old provider disabled.\n";
            return [](int, double&)->bool { return false; };
        }

        const auto patchesCopy = patches; // 捕获
        return [&faces, &id2idx, Tprev, fluxF, L, patchesCopy, &Tbc, &Pbc, epsFlux](int fid, double& Tb)->bool {
            // 1) patch 命中
            if (!L.isInAny(patchesCopy, fid)) return false;

            // 2) 温度为绝热（零法向导数）
            if (!isZeroFluxGradBC(Tbc, fid)) return false;

            // 3) 压力不是零通量（即该面可能有达西流动）
            if (isZeroFluxGradBC(Pbc, fid)) return false;

            // 4) 对温度对流是入流（flux<0）
            const int iF = fid - 1;  // 若 face id 非连续需改用 id->index
            const double flux = (*fluxF)[iF];
            if (!(flux < -epsFlux)) return false;

            // 5) 入边温度 = 上一时层(或外迭代)的 owner 单元温度
            const auto& F = faces.at(iF);
            const int    P = F.ownerCell;
            const size_t iP = id2idx.at(P);

            Tb = (*Tprev)[iP];
            return true;
            };
    }

} // namespace BoundaryUtils
