```c++
#pragma once

#include <string>
#include <iostream>
#include <algorithm>

#include "MeshManager.h"
#include "FieldRegistry.h"

namespace IMPES
{
inline bool TimeTerm_IMPES_Pressure(
        MeshManager& mgr,
        FieldRegistry& reg,
        double dt,
        const std::string& phi_name,
        const std::string& p_old_name,
        const std::string& p_eval_name,
        const std::string& rho_old_name,
        const std::string& rho_eval_name,
        const std::string& drho_dp_name,
        double rock_compressibility,
        const std::string& aC_name,
        const std::string& bC_name,
        const std::vector<char>* strong_mask = nullptr)
    {
        if (dt <= 0.0)
        {
            std::cerr << "[IMPES][TimeTerm] invalid dt.\n";
            return false;
        }

        auto& mesh = mgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();
    
        auto phi = reg.get<volScalarField>(phi_name);
        auto p_old = reg.get<volScalarField>(p_old_name);
        auto p_eval = reg.get<volScalarField>(p_eval_name);
        auto rho_old = reg.get<volScalarField>(rho_old_name);
        auto rho_eval = reg.get<volScalarField>(rho_eval_name);
        auto drho_dp = reg.get<volScalarField>(drho_dp_name);
    
        if (!phi || !p_old || !p_eval || !rho_old || !rho_eval || !drho_dp)
        {
            std::cerr << "[IMPES][TimeTerm] missing fields for pressure accumulation.\n";
            return false;
        }
    
        auto aC = reg.getOrCreate<volScalarField>(aC_name.c_str(), cells.size(), 0.0);
        auto bC = reg.getOrCreate<volScalarField>(bC_name.c_str(), cells.size(), 0.0);
        std::fill(aC->data.begin(), aC->data.end(), 0.0);
        std::fill(bC->data.begin(), bC->data.end(), 0.0);
    
        const double inv_dt = 1.0 / dt;
    
        for (const auto& c : cells)
        {
            if (c.id < 0) continue;
            const size_t i = id2idx.at(c.id);
            if (strong_mask && (*strong_mask)[i]) continue;
    
            const double V = std::max(0.0, c.volume);
            const double phi_i = std::max(0.0, std::min(1.0, (*phi)[i]));
            const double rho_n = std::max(0.0, (*rho_old)[i]);
            const double rho_star = std::max(0.0, (*rho_eval)[i]);
            const double drdp = std::max(0.0, (*drho_dp)[i]);
            const double p_n = (*p_old)[i];
            const double p_star = (*p_eval)[i];
    
            const double a = V * inv_dt * phi_i * (drdp + rho_star * rock_compressibility);
            const double b = V * inv_dt * (phi_i * rho_n
                - phi_i * rho_star
                + phi_i * drdp * p_star
                + phi_i * rho_star * rock_compressibility * p_n);
    
            (*aC)[i] = a;
            (*bC)[i] = b;
        }
    
        return true;
    }

} 

```

