#pragma once

#include <vector>
#include <string>

// 引入核心基础类型
#include "UserDefineVarType.h" 
// [关键修正] 引入 VolField 定义，因为函数参数现在直接依赖它
#include "VolField.h" 

// =========================================================
// 前向声明 (Forward Declarations)
// =========================================================
class Mesh;
class MeshManager_3D;
class FieldManager_3D;

namespace FVM
{
    /**
     * @class EDFM_GradientsOperation_3D
     * @brief 3D-EDFM 统一梯度计算器 (Type-Safe Version)
     * @details
     * 修正了类型转换错误，直接操作 volScalarField。
     */
    class EDFM_GradientsOperation_3D
    {
    public:
        /**
         * @brief 计算基岩单元的梯度 (Matrix Gradient)
         */
        static std::vector<Vector> ComputeMatrixGradients(
            const ::MeshManager_3D& meshMgr,
            const ::FieldManager_3D& fieldMgr,
            const std::string& scalarName,
            int smoothIters = 0);

        /**
         * @brief 计算裂缝微元的梯度 (Fracture Gradient)
         */
        static std::vector<Vector> ComputeFractureGradients(
            const ::MeshManager_3D& meshMgr,
            const ::FieldManager_3D& fieldMgr,
            const std::string& scalarName);

    private:
        struct Mat3;

        // [关键修正] 参数改为 const volScalarField&
        static Vector _computeGG_Matrix(
            const ::Mesh& mesh,
            const volScalarField& phi,
            int cellID);

        static void _performSmoothing(
            const ::Mesh& mesh,
            const std::vector<Vector>& inputGrad,
            std::vector<Vector>& outputGrad);
    };
}