#pragma once

#include <string>
#include <memory>
#include <vector>

#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "VolField.h"

/**
 * @class FieldManager_2D
 * @brief 2D-EDFM 专用场数据管理器
 * @details 统一管理 2D 基岩域、1D 裂缝域及交互域（NNC/FF）的物理场注册表。
 * 适配 MeshManager (2D-EDFM) 的拓扑结构。
 */
class FieldManager_2D
{
public:
    // =========================================================
    // 成员变量：三个域的独立注册表
    // =========================================================

    /// @brief 基岩体心场 (Size = N_matrix_cells)
    FieldRegistry matrixFields;

    /// @brief 裂缝体心场 (Size = N_total_frac_elements, 即 1D 线段总数)
    FieldRegistry fractureFields;

    /// @brief 交互域场 (用于 NNC 传导率 Size = N_pairs, 或 FF 传导率 Size = N_ff_connections)
    FieldRegistry nncFields;

    /// @brief 基岩网格面（边）数据 (Size = N_matrix_faces / edges)
    FaceFieldRegistry matrixFaceFields;

    /// @brief 裂缝内部连接面（节点）数据 (Size = N_frac_internal_nodes)
    /// @note 对于 1D 裂缝，"Face" 实际上是裂缝段之间的连接节点 (Vertex)
    FaceFieldRegistry fractureFaceFields;

    // =========================================================
    // 尺寸记录 (用于创建新场时的默认大小)
    // =========================================================
    size_t numMatrixCells = 0;      ///< 基岩单元数量
    size_t numFracCells = 0;        ///< 裂缝微元（线段）总数量
    size_t numNNCPairs = 0;         ///< NNC (基岩-裂缝) 交互对数量
    size_t numFFConnections = 0;    ///< F-F (裂缝-裂缝) 连接数量

    size_t numMatrixFaces = 0;      ///< 基岩面（2D中的边）数量
    size_t numFracFaces = 0;        ///< 裂缝面（1D中的点）数量

    // =========================================================
    // 初始化与配置
    // =========================================================
    /**
     * @brief 初始化各域的网格数量
     * @details 在 MeshManager 完成网格生成和拓扑构建后调用此函数，
     * 确保后续创建场时能分配正确的内存大小。
     * * @param nMatrix 基岩单元数量 (mesh_.getCells().size())
     * @param nFrac 裂缝微元总数量 (MeshManager::getTotalDOFCount() - nMatrix)
     * @param nNNC NNC 交互对数量 (即 CellLocalIndexToFracElemSolverIndexMap_ 的元素总和)
     * @param nFF F-F 连接数量 (FractureNetwork 中记录的有效连接数)
     * @param nMatFaces 基岩面数量 (mesh_.getFaces().size())
     * @param nFracFaces 裂缝内部连接面数量 (即裂缝段之间的通量接口总数)
     */
    void InitSizes(size_t nMatrix, size_t nFrac, size_t nNNC, size_t nFF, size_t nMatFaces, size_t nFracFaces)
    {
        numMatrixCells = nMatrix;
        numFracCells = nFrac;
        numNNCPairs = nNNC;
        numFFConnections = nFF;

        numMatrixFaces = nMatFaces;
        numFracFaces = nFracFaces;
    }

    // =========================================================
    // 辅助封装接口 (Helpers)
    // =========================================================

    // ---------------------------------------------------------
    // 1. Matrix Domain (基岩)
    // ---------------------------------------------------------

    /**
     * @brief 创建基岩标量场 (double)
     * @param name 场名称 (如 "p", "T")
     * @param initVal 初始值 (默认 0.0)
     * @return 指向创建场的智能指针
     */
    std::shared_ptr<volScalarField> createMatrixScalar(const std::string& name, double initVal = 0.0)
    {
        return matrixFields.create<volScalarField>(name, numMatrixCells, initVal);
    }

    /**
     * @brief 获取基岩标量场
     * @param name 场名称
     * @return 智能指针 (如果不存在返回 nullptr)
     */
    std::shared_ptr<volScalarField> getMatrixScalar(const std::string& name) const
    {
        return matrixFields.get<volScalarField>(name);
    }

    /**
     * @brief 获取或创建基岩标量场
     * @details 若场已存在则返回，否则按默认值创建
     */
    std::shared_ptr<volScalarField> getOrCreateMatrixScalar(const std::string& name, double initVal = 0.0)
    {
        return matrixFields.getOrCreate<volScalarField>(name, numMatrixCells, initVal);
    }

    /**
     * @brief 检查基岩场是否存在
     */
    bool hasMatrixField(const std::string& name) const
    {
        return matrixFields.has(name);
    }

    // --- 基岩面数据 (Flux / Transmissibility) ---

    /**
     * @brief 创建基岩面标量场 (如面通量)
     */
    std::shared_ptr<faceScalarField> createMatrixFaceScalar(const std::string& name, double initVal = 0.0)
    {
        return matrixFaceFields.create<faceScalarField>(name, numMatrixFaces, initVal);
    }

    std::shared_ptr<faceScalarField> getMatrixFaceScalar(const std::string& name) const
    {
        return matrixFaceFields.get<faceScalarField>(name);
    }

    std::shared_ptr<faceScalarField> getOrCreateMatrixFaceScalar(const std::string& name, double initVal = 0.0)
    {
        return matrixFaceFields.getOrCreate<faceScalarField>(name, numMatrixFaces, initVal);
    }

    bool hasMatrixFaceField(const std::string& name) const
    {
        return matrixFaceFields.has(name);
    }

    // =========================================================
    // [新增] 基岩域 (Matrix) AD 泛化接口
    // =========================================================

    /**
     * @brief 创建基岩自动微分标量场 (ADVar)
     * @tparam N 独立自变量的数量
     * @param name 场名称 (如 "p_AD", "T_AD")
     * @param initVal 初始自动微分值 (默认物理值为0，梯度为零向量)
     * @return 指向创建的 AD 场的智能指针
     */
    template<int N>
    std::shared_ptr<volADField<N>> createMatrixADScalar(const std::string& name, const ADVar<N>& initVal = ADVar<N>(0.0))
    {
        return matrixFields.create<volADField<N>>(name, numMatrixCells, initVal);
    }

    /**
     * @brief 获取基岩自动微分标量场
     * @tparam N 独立自变量的数量
     * @param name 场名称
     * @return 智能指针 (如果不存在返回 nullptr)
     */
    template<int N>
    std::shared_ptr<volADField<N>> getMatrixADScalar(const std::string& name) const
    {
        return matrixFields.get<volADField<N>>(name);
    }

    /**
     * @brief 获取或创建基岩自动微分标量场
     * @tparam N 独立自变量的数量
     */
    template<int N>
    std::shared_ptr<volADField<N>> getOrCreateMatrixADScalar(const std::string& name, const ADVar<N>& initVal = ADVar<N>(0.0))
    {
        return matrixFields.getOrCreate<volADField<N>>(name, numMatrixCells, initVal);
    }

    /**
     * @brief 创建基岩面 (边) 的自动微分标量场
     * @tparam N 独立自变量的数量
     */
    template<int N>
    std::shared_ptr<faceADField<N>> createMatrixFaceADScalar(const std::string& name, const ADVar<N>& initVal = ADVar<N>(0.0))
    {
        return matrixFaceFields.create<faceADField<N>>(name, numMatrixFaces, initVal);
    }

    /**
     * @brief 获取基岩面 (边) 的自动微分标量场
     * @tparam N 独立自变量的数量
     */
    template<int N>
    std::shared_ptr<faceADField<N>> getMatrixFaceADScalar(const std::string& name) const
    {
        return matrixFaceFields.get<faceADField<N>>(name);
    }

    // ---------------------------------------------------------
    // 2. Fracture Domain (裂缝 - 1D Segments)
    // ---------------------------------------------------------

    /**
     * @brief 创建裂缝标量场
     * @details 对应 1D 裂缝线段上的物理量
     */
    std::shared_ptr<volScalarField> createFractureScalar(const std::string& name, double initVal = 0.0)
    {
        return fractureFields.create<volScalarField>(name, numFracCells, initVal);
    }

    std::shared_ptr<volScalarField> getFractureScalar(const std::string& name) const
    {
        return fractureFields.get<volScalarField>(name);
    }

    std::shared_ptr<volScalarField> getOrCreateFractureScalar(const std::string& name, double initVal = 0.0)
    {
        return fractureFields.getOrCreate<volScalarField>(name, numFracCells, initVal);
    }

    bool hasFractureField(const std::string& name) const
    {
        return fractureFields.has(name);
    }

    // --- 裂缝连接面数据 (Internal Nodes) ---

    /**
     * @brief 创建裂缝内部连接面（节点）标量场
     * @details 对应裂缝段之间连接点处的通量或传导率
     */
    std::shared_ptr<faceScalarField> createFractureFaceScalar(const std::string& name, double initVal = 0.0)
    {
        return fractureFaceFields.create<faceScalarField>(name, numFracFaces, initVal);
    }

    std::shared_ptr<faceScalarField> getFractureFaceScalar(const std::string& name) const
    {
        return fractureFaceFields.get<faceScalarField>(name);
    }

    std::shared_ptr<faceScalarField> getOrCreateFractureFaceScalar(const std::string& name, double initVal = 0.0)
    {
        return fractureFaceFields.getOrCreate<faceScalarField>(name, numFracFaces, initVal);
    }

    bool hasFractureFaceField(const std::string& name) const
    {
        return fractureFaceFields.has(name);
    }

    // =========================================================
    // [新增] 裂缝域 (Fracture) AD 泛化接口
    // =========================================================

    template<int N>
    std::shared_ptr<volADField<N>> createFractureADScalar(const std::string& name, const ADVar<N>& initVal = ADVar<N>(0.0))
    {
        return fractureFields.create<volADField<N>>(name, numFracCells, initVal);
    }

    template<int N>
    std::shared_ptr<volADField<N>> getFractureADScalar(const std::string& name) const
    {
        return fractureFields.get<volADField<N>>(name);
    }

    template<int N>
    std::shared_ptr<volADField<N>> getOrCreateFractureADScalar(const std::string& name, const ADVar<N>& initVal = ADVar<N>(0.0))
    {
        return fractureFields.getOrCreate<volADField<N>>(name, numFracCells, initVal);
    }

    template<int N>
    std::shared_ptr<faceADField<N>> createFractureFaceADScalar(const std::string& name, const ADVar<N>& initVal = ADVar<N>(0.0))
    {
        return fractureFaceFields.create<faceADField<N>>(name, numFracFaces, initVal);
    }

    template<int N>
    std::shared_ptr<faceADField<N>> getFractureFaceADScalar(const std::string& name) const
    {
        return fractureFaceFields.get<faceADField<N>>(name);
    }

    // ---------------------------------------------------------
    // 3. Interaction Domain (NNC & FF)
    // ---------------------------------------------------------

    // --- NNC (Matrix-Fracture) ---
    /**
     * @brief 创建 NNC 交互标量场 (如 T_NNC)
     */
    std::shared_ptr<volScalarField> createNNCScalar(const std::string& name, double initVal = 0.0)
    {
        return nncFields.create<volScalarField>(name, numNNCPairs, initVal);
    }

    std::shared_ptr<volScalarField> getNNCScalar(const std::string& name) const
    {
        return nncFields.get<volScalarField>(name);
    }


    // --- FF (Fracture-Fracture) ---
    /**
     * @brief 创建 F-F 交互标量场
     * @details 共享 nncFields 注册表容器，但使用 numFFConnections 作为大小。
     * 只要命名不冲突 (如 "T_FF" vs "T_NNC")，可以安全共存。
     */
    std::shared_ptr<volScalarField> createFFScalar(const std::string& name, double initVal = 0.0)
    {
        return nncFields.create<volScalarField>(name, numFFConnections, initVal);
    }

    std::shared_ptr<volScalarField> getFFScalar(const std::string& name) const
    {
        return nncFields.get<volScalarField>(name);
    }

    // =========================================================
    // [新增] 交互域 (NNC & FF) AD 泛化接口
    // =========================================================

    /**
     * @brief 创建 NNC (基岩-裂缝) 交互自动微分标量场
     * @tparam N 独立自变量的数量
     */
    template<int N>
    std::shared_ptr<volADField<N>> createNNCADScalar(const std::string& name, const ADVar<N>& initVal = ADVar<N>(0.0))
    {
        return nncFields.create<volADField<N>>(name, numNNCPairs, initVal);
    }

    template<int N>
    std::shared_ptr<volADField<N>> getNNCADScalar(const std::string& name) const
    {
        return nncFields.get<volADField<N>>(name);
    }

    /**
     * @brief 创建 FF (裂缝-裂缝) 交互自动微分标量场
     * @tparam N 独立自变量的数量
     */
    template<int N>
    std::shared_ptr<volADField<N>> createFFADScalar(const std::string& name, const ADVar<N>& initVal = ADVar<N>(0.0))
    {
        return nncFields.create<volADField<N>>(name, numFFConnections, initVal);
    }

    template<int N>
    std::shared_ptr<volADField<N>> getFFADScalar(const std::string& name) const
    {
        return nncFields.get<volADField<N>>(name);
    }
};