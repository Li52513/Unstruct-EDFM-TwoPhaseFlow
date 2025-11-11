#pragma once
#include <string> // std::string
#include <memory> // std::shared_ptr
#include <tuple> // std::tie
#include <utility> // std::forward
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"

// ===== 工具：从单元场取标量 =====
inline double cellScalar(const FieldRegistry& reg, const Mesh& mesh,
	const char* name, int cellId, double fallback = 0.0)
{
	auto fld = reg.get<volScalarField>(name);
	if (!fld) return fallback; // field not found
	auto it = mesh.getCellId2Index().find(cellId);
	if (it == mesh.getCellId2Index().end()) return fallback; // cellId not found
	return (*fld)[it->second];
}

// ===== 工具：从单元场取向量 =====
inline Vector cellVector(const FieldRegistry& reg, const Mesh& mesh,
	const char* name, int cellId, const Vector& fallback = Vector(0.0, 0.0, 0.0))
{
	auto fld = reg.get<volVectorField>(name);
	if (!fld) return fallback; // field not found
	auto it = mesh.getCellId2Index().find(cellId);
	if (it == mesh.getCellId2Index().end()) return fallback; // cellId not found
	return (*fld)[it->second];
}


