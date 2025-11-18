```C++
#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "FaceSignMask.hpp"

namespace IMPES
{
    template<typename T>
    inline T clampValue(T v, T lo, T hi)
    {
        if (v < lo) return lo;
        if (v > hi) return hi;
        return v;
    }

    struct FluxSplitConfig
    {
        std::string total_mass_flux = "mf_total";    
        std::string water_mass_flux = "mf_w";        
        std::string gas_mass_flux = "mf_g";          
        std::string fractional_flow_face = "fw_face";
        std::string lambda_gas = "lambda_g";         
        std::string saturation = "s_w";              
        double flux_sign_epsilon = 1e-15;            
    
    struct FluxSplitResult
    {
        std::shared_ptr<faceScalarField> mf_w;
        std::shared_ptr<faceScalarField> mf_g;
        std::shared_ptr<faceScalarField> fw_face;
        FaceSignUpdateInfo signInfo;
        std::vector<int> flippedFaces;
    };
    
        
    inline bool splitTwoPhaseMassFlux(
        MeshManager& mgr,
        FieldRegistry& reg,
        FaceFieldRegistry& freg,
        const FluxSplitConfig& cfg,
        FaceSignMask* mask = nullptr,
        FluxSplitResult* result = nullptr)
    {
        auto mf_total = freg.get<faceScalarField>(cfg.total_mass_flux.c_str());
        if (!mf_total)
        {
            std::cerr << "[IMPES][Flux] missing total mass flux field '" << cfg.total_mass_flux << "'.\n";
            return false;
        }
    
        auto lambda_w = reg.get<volScalarField>(cfg.lambda_water.c_str());
        auto lambda_g = reg.get<volScalarField>(cfg.lambda_gas.c_str());
        if (!lambda_w || !lambda_g)
        {
            std::cerr << "[IMPES][Flux] missing lambda fields '" << cfg.lambda_water
                      << "' or '" << cfg.lambda_gas << "'.\n";
            return false;
        }
    
        auto mf_w = freg.getOrCreate<faceScalarField>(cfg.water_mass_flux.c_str(), mf_total->data.size(), 0.0);
        auto mf_g = freg.getOrCreate<faceScalarField>(cfg.gas_mass_flux.c_str(), mf_total->data.size(), 0.0);
        std::shared_ptr<faceScalarField> fw_face;
        if (!cfg.fractional_flow_face.empty())
        {
            fw_face = freg.getOrCreate<faceScalarField>(cfg.fractional_flow_face.c_str(), mf_total->data.size(), 0.0);
        }
    
        Mesh& mesh = mgr.mesh();
        const auto& faces = mesh.getFaces();
        const auto& id2idx = mesh.getCellId2Index();
    
        for (const auto& F : faces)
        {
            const int iF = F.id - 1;
            const double flux = (*mf_total)[iF];
            double fw_value = 0.0;
            double mfW = 0.0;
            double mfG = 0.0;
    
            if (std::abs(flux) > cfg.flux_sign_epsilon)
            {
                const int P = F.ownerCell;
                const int N = F.neighborCell;
                int upCell = -1;
                if (flux >= 0.0 || N < 0)
                {
                    upCell = (P >= 0) ? static_cast<int>(id2idx.at(P)) : -1;
                }
                else
                {
                    upCell = (N >= 0) ? static_cast<int>(id2idx.at(N)) : (P >= 0 ? static_cast<int>(id2idx.at(P)) : -1);
                }
    
                if (upCell >= 0)
                {
                    const double lamW = std::max((*lambda_w)[upCell], 0.0);
                    const double lamG = std::max((*lambda_g)[upCell], 0.0);
                    const double denom = std::max(lamW + lamG, cfg.min_lambda);
                    fw_value = clampValue(lamW / denom, 0.0, 1.0);
                }
                else
                {
                    fw_value = 0.0;
                }
    
                mfW = fw_value * flux;
                mfG = flux - mfW;
            }
    
            (*mf_w)[iF] = mfW;
            (*mf_g)[iF] = mfG;
            if (fw_face)
            {
                (*fw_face)[iF] = fw_value;
            }
        }
    
        FaceSignUpdateInfo info;
        if (mask)
        {
            std::vector<int> flippedTmp;
            auto& flipped = result ? result->flippedFaces : flippedTmp;
            info = updateFaceSignMask_fromFlux(*mf_total, cfg.flux_sign_epsilon, *mask, flipped);
        }
    
        if (result)
        {
            result->mf_w = mf_w;
            result->mf_g = mf_g;
            result->fw_face = fw_face;
            result->signInfo = info;
            if (!mask)
            {
                result->flippedFaces.clear();
            }
        }
    
        return true;
    }

} 
```

