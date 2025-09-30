#pragma once
#include<vector>
#include<string>
#include<algorithm>
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FieldAcessForDiscre.h"
#include "Diff_TPFA_PermeabilityOperation.h"

// 移动率/导热率提供器
//   double mobilityAlong(const Mesh&, const FieldRegistry&, int cellId,
//                        const Vector& ehat) const;
// 返回 λ_e（例如 k_e/μ、k_r*k/μ、或导热率 λ）

struct DarcyWaterMobility_singlePhase
{
	double k_default;
	explicit DarcyWaterMobility_singlePhase(double k_iso_default) : k_default(k_iso_default) {}
	double mobilityAlong(const Mesh& mesh, const FieldRegistry& reg, int cellId,
		const Vector& ehat) const
	{
		const double mu = std::max(cellScalar(reg, mesh, "mu_w", cellId, 1.0), 1e-30);
		double kxx, kyy, kzz;
		getKdiag(reg, mesh, cellId, k_default, kxx, kyy, kzz);
		const double ke = std::max(kEffAlong(ehat, kxx, kyy, kzz), 1e-30);
		return ke / mu;
	}
};

struct DarcyCO2Mobility_singlePhase {
	double k_default;
	explicit DarcyCO2Mobility_singlePhase(double k_iso_default) : k_default(k_iso_default) {}
	double mobilityAlong(const Mesh& mesh, const FieldRegistry& reg, int cellId,
		const Vector& ehat) const
	{
		const double mu = std::max(cellScalar(reg, mesh, "mu_g", cellId, 1.0), 1e-30);
		double kxx, kyy, kzz;
		getKdiag(reg, mesh, cellId, k_default, kxx, kyy, kzz);
		const double ke = std::max(kEffAlong(ehat, kxx, kyy, kzz), 1e-30);
		return ke / mu;
	}
};

//由于两相流中各相渗透率需要考虑迎风性，故不再提供各向异性两相流流度计算示例

// 纯扩散（例如导热率）示例：直接从 cell 标量场 "lambda_eff" 读取
struct IsotropicLambdaFromField
{
	std::string name; // e.g. "lambda_eff"
	explicit IsotropicLambdaFromField(std::string n) : name(std::move(n)) {}
	double mobilityAlong(const Mesh& mesh, const FieldRegistry& reg, int cellId,
		const Vector& /*ehat*/) const
	{
		return std::max(cellScalar(reg, mesh, name.c_str(), cellId, 0.0), 0.0);
	}
};