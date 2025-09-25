#pragma once
#include <string> // std::string
#include <memory> // std::shared_ptr
#include <tuple> // std::tie
#include <utility> // std::forward
#include "MeshManager.h"
#include "FieldRegistry.h"
#include "FieldAcessForDiscre.h"

// ===== 工具：K 的对角三主值 + 方向等效 =====
inline void getKdiag(const FieldRegistry& reg, const Mesh& mesh, int cellId, double k_default,
    double& kxx, double& kyy, double& kzz)
{
    kxx = cellScalar(reg, mesh, "kxx", cellId, k_default);
    kyy = cellScalar(reg, mesh, "kyy", cellId, k_default);
    kzz = cellScalar(reg, mesh, "kzz", cellId, k_default);
}

// 将各相渗透率等效到e方向上 k_e = e^T K e（K=diag）
inline double kEffAlong(const Vector& e, double kxx, double kyy, double kzz)
{
    return kxx * e.m_x * e.m_x + kyy * e.m_y * e.m_y + kzz * e.m_z * e.m_z;
}