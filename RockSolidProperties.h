#pragma once
#include "Cell.h"
#include "MeshManager.h"
#include "FieldRegistry.h"
#include"PropertiesSummary.h"
#include "SolverContrlStrName.h"

namespace rock
{
    struct Rock_heterogeneous_Parameters
    {

        double phi_r_high = 0.4;         // 孔隙度
        double kxx_high = 1e-13;           // 渗透率，m²
        double kyy_high = 1e-13;
        double kzz_high = 1e-13;

        double phi_r_medium = 0.3;         // 孔隙度
        double kxx_medium = 1e-14;           // 渗透率，m²
        double kyy_medium = 1e-14;
        double kzz_medium = 1e-14;

        double phi_r_low = 0.2;         // 孔隙度
        double kxx_low = 1e-15;           // 渗透率，m²
        double kyy_low = 1e-15;
        double kzz_low = 1e-15;
    };


    /// 根据 cell.region 和 (P,T) 返回基岩固相物性
    SolidProperties_RockMatrix computeSolidProperties(Cell::RegionType region, double P, double T);

    // 具体常数或经验公式
    static constexpr double BASE_RHO = 2650.0;  // 密度
    static constexpr double BASE_CP = 1000.0;   // 比热容
    static constexpr double BASE_K = 2.5;       // 导热系数
    static constexpr double BASE_COMP = 1e-8;   // 可压缩系数

    inline bool ensure_RockProp_Fields(FieldRegistry& reg, std::size_t n)
    {
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Rock().rho_tag, n, BASE_RHO);        // rho，kg/m³
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Rock().cp_tag, n, BASE_CP);          // Cp，J/(kg·K)
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Rock().lambda_tag, n, BASE_K);       // k，W/(m·K)
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Rock().c_r_tag, n, BASE_COMP);       // 可压缩系数，1/Pa
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Rock().phi_tag, n, 0.2);            // 孔隙度
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Rock().k_xx_tag, n, 1e-14);          // 渗透率，m²
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Rock().k_yy_tag, n, 1e-14);
		reg.getOrCreate<volScalarField>(PhysicalProperties_string::Rock().k_zz_tag, n, 1e-14);
		return true;
    }

	// 按区域划分后，计算每个cell的基岩物性参数
    inline SolidProperties_RockMatrix computeSolidProperties(Cell::RegionType region, double P, double T)
    {
        Rock_heterogeneous_Parameters r;
        SolidProperties_RockMatrix s;
        s.rho_r = BASE_RHO;
        s.cp_r = BASE_CP;
        s.k_r = BASE_K;
        s.compressibility = BASE_COMP;

        // 现在均是常数，后期可以优化放入不同的方程
        switch (region)
        {
        case Cell::RegionType::Low:
            s.phi_r = r.phi_r_low;  s.kxx = r.kxx_low; s.kyy = r.kyy_low; s.kzz = r.kzz_low;  break;
        case Cell::RegionType::Medium:
            s.phi_r = r.phi_r_medium;  s.kxx = r.kxx_medium; s.kyy = r.kyy_medium; s.kzz = r.kzz_medium; break;
        case Cell::RegionType::High:
            s.phi_r = r.phi_r_high;  s.kxx = r.kxx_high; s.kyy = r.kyy_high; s.kzz = r.kzz_high; break;
        }
        return s;
    }

	inline bool computer_rock_properties (MeshManager& mgr, FieldRegistry& reg, const std::string& p_field, const std::string& T_field)
	{

		auto& mesh = mgr.mesh();
		const auto& cells = mesh.getCells();
		const size_t n = cells.size();
		auto TF = reg.get<volScalarField>(T_field);
		auto pF = reg.get<volScalarField>(p_field);

		ensure_RockProp_Fields(reg, n);

		auto rho_rF = reg.get<volScalarField>(PhysicalProperties_string::Rock().rho_tag);
		auto cp_rF = reg.get<volScalarField>(PhysicalProperties_string::Rock().cp_tag);
		auto k_rF = reg.get<volScalarField>(PhysicalProperties_string::Rock().lambda_tag);
		auto c_rF = reg.get<volScalarField>(PhysicalProperties_string::Rock().c_r_tag);
		auto phi_rF = reg.get<volScalarField>(PhysicalProperties_string::Rock().phi_tag);
		auto kxx_rF = reg.get<volScalarField>(PhysicalProperties_string::Rock().k_xx_tag);
		auto kyy_rF = reg.get<volScalarField>(PhysicalProperties_string::Rock().k_yy_tag);
		auto kzz_rF = reg.get<volScalarField>(PhysicalProperties_string::Rock().k_zz_tag);
		for (size_t ic = 0; ic < cells.size(); ++ic)
		{
			const auto& c = cells[ic];
			const size_t i = mesh.getCellId2Index().at(c.id);
			double T = (*TF)[i];
			double P = (*pF)[i];
			const auto sp = computeSolidProperties(c.region, P, T);
			(*rho_rF)[i] = sp.rho_r;
			(*cp_rF)[i] = sp.cp_r;
			(*k_rF)[i] = sp.k_r;
			(*c_rF)[i] = sp.compressibility;
			(*phi_rF)[i] = sp.phi_r;
			(*kxx_rF)[i] = sp.kxx;
			(*kyy_rF)[i] = sp.kyy;
			(*kzz_rF)[i] = sp.kzz;
		}
		return true;
	}
}