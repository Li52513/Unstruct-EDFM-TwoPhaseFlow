/**
 * @file 3D_BoundaryConditionManager.cpp
 * @brief 3D 边界条件统一管理器实现
 */
#include <iomanip>
#include "3D_BoundaryConditionManager.h"

namespace BoundarySetting {

    // =========================================================
    // 构造与析构
    // =========================================================
    BoundaryConditionManager_3D::BoundaryConditionManager_3D() {
        // 初始化为空
    }

    BoundaryConditionManager_3D::~BoundaryConditionManager_3D() {
        bcRegistry_.clear();
    }

    // =========================================================
    // 设置接口实现
    // =========================================================

    void BoundaryConditionManager_3D::SetDirichletBC(int tagID, double constantValue) {
        // 使用 Lambda 封装常数返回
        BCDefinition def;
        def.a = 1.0;
        def.b = 0.0;
        def.type = BoundaryType::Dirichlet;
        def.valFunc = [constantValue](const Vector&) { return constantValue; };
        def.desc = "Dirichlet (Constant=" + std::to_string(constantValue) + ")";

        bcRegistry_[tagID] = def;
    }

    void BoundaryConditionManager_3D::SetLinearDirichletBC(int tagID, double refValue, double refCoord, double gradient, int axis) {
        // 使用 Lambda 封装线性方程
        // Val = Ref + Grad * (Pos[axis] - Ref)
        BCDefinition def;
        def.a = 1.0;
        def.b = 0.0;
        def.type = BoundaryType::Dirichlet;

        // 捕获所有参数
        def.valFunc = [=](const Vector& p) -> double {
            double pos = 0.0;
            if (axis == 0) pos = p.m_x;
            else if (axis == 1) pos = p.m_y;
            else pos = p.m_z;

            return refValue + gradient * (pos - refCoord);
            };

        std::string axisStr = (axis == 0 ? "X" : (axis == 1 ? "Y" : "Z"));
        def.desc = "Dirichlet (Linear along " + axisStr + ", Grad=" + std::to_string(gradient) + ")";

        bcRegistry_[tagID] = def;
    }



    void BoundaryConditionManager_3D::SetNeumannBC(int tagID, double constantFlux) {
        BCDefinition def;
        def.a = 0.0;
        def.b = 1.0;
        def.type = BoundaryType::Neumann;
        def.valFunc = [constantFlux](const Vector&) { return constantFlux; };
        def.desc = "Neumann (Constant=" + std::to_string(constantFlux) + ")";

        bcRegistry_[tagID] = def;
    }
    // =========================================================
    // 查询接口实现
    // =========================================================
    BCCoefficients BoundaryConditionManager_3D::GetBCCoefficients(int tagID, const Vector& faceCenter) const {
        auto it = bcRegistry_.find(tagID);
        if (it != bcRegistry_.end()) {
            const auto& def = it->second;
            // [Core] 调用函数对象，计算当前坐标下的 c 值
            double c_val = def.valFunc(faceCenter);
            return BCCoefficients(def.a, def.b, c_val, def.type);
        }
        else {
            throw std::runtime_error("[Error] No BC found for Tag ID " + std::to_string(tagID));
        }
    }

    bool BoundaryConditionManager_3D::HasBC(int tagID) const {
        return bcRegistry_.find(tagID) != bcRegistry_.end();
    }

    void BoundaryConditionManager_3D::Clear() { bcRegistry_.clear(); }

    void BoundaryConditionManager_3D::PrintSummary(const std::string& name) const {
        std::cout << "--- Boundary Conditions Summary for [" << name << "] ---" << std::endl;
        if (bcRegistry_.empty()) {
            std::cout << "  (No BCs set)" << std::endl;
            return;
        }
        for (const auto& pair : bcRegistry_) {
            std::cout << "  Tag " << std::setw(3) << pair.first << " : " << pair.second.desc << std::endl;
        }
        std::cout << "----------------------------------------------------" << std::endl;
    }

}