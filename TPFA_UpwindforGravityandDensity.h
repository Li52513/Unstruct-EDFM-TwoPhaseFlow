#pragma once
#include<vector>
#include<string>
#include<algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FieldAcessForDiscre.h"

// 密度
//   double rhoUp(const Mesh&, const FieldRegistry&, int P, int N,
//                double pP, double pN, const Vector& CP, const Vector& CN,
//                const GravUpwind&) const;
//   double rhoBar(const Mesh&, const FieldRegistry&, int P, int N) const;

//===========小工具==================//
//开关：重力是否用迎风格式，计算重力势能
struct GravUpwind {
	Vector g;
	bool use_potential = true;
	inline double dot_pos(const Vector& p)const { return g.m_x * p.m_x + g.m_y * p.m_y + g.m_z * p.m_z; }
};



struct UpwindDensityByPotential_water {
	std::string rho_name = "rho_w";
	bool use_potential = true;
	explicit UpwindDensityByPotential_water(bool usePhi = true) : use_potential(usePhi) {}

	double rhoBar(const Mesh& mesh, const FieldRegistry& reg, int P, int N, double gamma) const
	{
		const double rP = cellScalar(reg, mesh, rho_name.c_str(), P, 0.0);
		const double rN = cellScalar(reg, mesh, rho_name.c_str(), N, 0.0);
		return (1.0 - gamma) * rP + gamma * rN;
	}
	double rhoUp(const Mesh& mesh, const FieldRegistry& reg, int P, int N,
		double pP, double pN, const Vector& CP, const Vector& CN,
		const GravUpwind& gu) const
	{
		const double rP = cellScalar(reg, mesh, rho_name.c_str(), P, 0.0);
		const double rN = cellScalar(reg, mesh, rho_name.c_str(), N, 0.0);
		if (!use_potential) return (pP - pN > 0.0) ? rP : rN;
		const double rBar = 0.5 * (rP + rN);
		const double phiP = pP - rBar * gu.dot_pos(CP);
		const double phiN = pN - rBar * gu.dot_pos(CN);
		return (phiP - phiN > 0.0) ? rP : rN;
	}
};

struct UpwindDensityByPotential_CO2 {
	std::string rho_name = "rho_g";
	bool use_potential = true;
	explicit UpwindDensityByPotential_CO2(bool usePhi = true) : use_potential(usePhi) {}
	double rhoBar(const Mesh& mesh, const FieldRegistry& reg, int P, int N,double gamma) const 
	{
		const double rP = cellScalar(reg, mesh, rho_name.c_str(), P, 0.0);
		const double rN = cellScalar(reg, mesh, rho_name.c_str(), N, 0.0);
		return (1.0 - gamma) * rP + gamma * rN;
	}
	double rhoUp(const Mesh& mesh, const FieldRegistry& reg, int P, int N,
		double pP, double pN, const Vector& CP, const Vector& CN,
		const GravUpwind& gu) const
	{
		const double rP = cellScalar(reg, mesh, rho_name.c_str(), P, 0.0);
		const double rN = cellScalar(reg, mesh, rho_name.c_str(), N, 0.0);
		if (!use_potential) return (pP - pN > 0.0) ? rP : rN;
		const double rBar = 0.5 * (rP + rN);
		const double phiP = pP - rBar * gu.dot_pos(CP);
		const double phiN = pN - rBar * gu.dot_pos(CN);
		return (phiP - phiN > 0.0) ? rP : rN;
	}
};

// 无密度（纯扩散）策略  取rhoBar=0，rhoUp=1,即可消除重力项和迎风密度
struct NoDensity {
	double rhoBar(const Mesh&, const FieldRegistry&, int, int,double) const { return 0.0; }
	double rhoUp(const Mesh&, const FieldRegistry&, int, int,
		double, double, const Vector&, const Vector&,
		const GravUpwind&) const {
		return 1.0;
	}
};