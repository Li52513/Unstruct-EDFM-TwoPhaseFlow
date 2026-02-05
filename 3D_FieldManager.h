#pragma once

#include <string>
#include <memory>
#include <vector>

#include "FieldRegistry.h"
#include "FaceFieldRegistry.h"
#include "VolField.h"

/**
 * @class FieldManager_3D
 * @brief 3D-EDFM 专用场数据管理器
 * @details 管理三个独立的域 (Matrix, Fracture, Interaction) 的场注册表。
 */
class FieldManager_3D
{
public:
    // =========================================================
    // 成员变量：三个域的独立注册表
    // =========================================================
    FieldRegistry matrixFields;   ///< 基岩域场 (Size = N_matrix_cells)
    FieldRegistry fractureFields; ///< 裂缝域场 (Size = N_total_frac_micro_elements)
    FieldRegistry nncFields;      ///< 交互域场 (NNC: Size = N_pairs, FF: Size = N_ff_segs)

    FaceFieldRegistry matrixFaceFields;    ///< 基岩面数据
    FaceFieldRegistry fractureEdgeFields;  ///< 裂缝边数据

    // =========================================================
    // 尺寸记录 (用于创建新场时的默认大小)
    // =========================================================
    size_t numMatrixCells = 0;
    size_t numFracCells = 0;
    size_t numNNCPairs = 0;
    size_t numFFConnections = 0;

    size_t numMatrixFaces = 0;
    size_t numFracEdges = 0;

    // =========================================================
    // 初始化与配置
    // =========================================================
    /**
     * @brief 初始化各域的网格数量
     * @param nMatrix 基岩单元数量
     * @param nFrac 裂缝微元总数量
     * @param nNNC NNC 交互对数量
     * @param nFF F-F 连接数量
     */
    void InitSizes(size_t nMatrix, size_t nFrac, size_t nNNC, size_t nFF, size_t nMatFaces, size_t nFracEdges)
    {
        numMatrixCells = nMatrix;
        numFracCells = nFrac;
        numNNCPairs = nNNC;
        numFFConnections = nFF;

        numMatrixFaces = nMatFaces;
        numFracEdges = nFracEdges;

    }

    // =========================================================
    // 辅助封装接口 (Helpers)
    // 严格匹配 FieldRegistry.h 的 create/get/has 接口
    // =========================================================

    // ---------------------------------------------------------
    // 1. Matrix Domain (基岩)
    // ---------------------------------------------------------

    /**
     * @brief 创建基岩标量场 (double)
     * @param name 场名称
     * @param initVal 初始值 (默认 0.0)
     * @return 指向创建场的智能指针
     */
    std::shared_ptr<volScalarField> createMatrixScalar(const std::string& name, double initVal = 0.0)
    {
        // FieldRegistry::create 内部调用 make_shared<FieldType>(name, args...)
        // VolField 构造函数: VolField(name, size, initVal)
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
     * @details 如果场已存在直接返回，否则按默认值创建
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

    /**
     * @brief 创建基岩网格面的标量场 (double)
     * @param name 场名称
     * @param initVal 初始值 (默认 0.0)
     * @return 指向创建场的智能指针
     */
    std::shared_ptr<faceScalarField> createMatrixFaceScalar(const std::string& name, double initVal = 0.0) {
        return matrixFaceFields.create<faceScalarField>(name, numMatrixFaces, initVal);
    }

    /**
    * @brief 获取基岩网格面的标量场
    * @param name 场名称
    * @return 智能指针 (如果不存在返回 nullptr)
    */
    std::shared_ptr<faceScalarField> getMatrixFaceScalar(const std::string& name) const {
        return matrixFaceFields.get<faceScalarField>(name);
    }

    /**
     * @brief 获取或创建基岩网格面的标量场
     * @details 如果场已存在直接返回，否则按默认值创建
     */
    std::shared_ptr<faceScalarField> getOrCreateMatrixFaceScalar(const std::string& name, double initVal = 0.0)
    {
        return matrixFaceFields.getOrCreate<faceScalarField>(name, numMatrixFaces, initVal);
    }

    /**
     * @brief 检查基岩网格面场是否存在
     */
    bool hasMatrixFaceField(const std::string& name) const
    {
        return matrixFaceFields.has(name);
    }

    // ---------------------------------------------------------
    // 2. Fracture Domain (裂缝)
    // ---------------------------------------------------------
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

    std::shared_ptr<faceScalarField> createFractureEdgeScalar(const std::string& name, double initVal = 0.0) {
        return fractureEdgeFields.create<faceScalarField>(name, numFracEdges, initVal);
    }

    std::shared_ptr<faceScalarField> getFractureEdgeScalar(const std::string& name) const {
        return fractureEdgeFields.get<faceScalarField>(name);
    }

    std::shared_ptr<faceScalarField> getOrCreateFractureEdgeScalar(const std::string& name, double initVal = 0.0) {
        return fractureEdgeFields.getOrCreate<faceScalarField>(name, numFracEdges, initVal);
    }

    bool hasFractureEdgeField(const std::string& name) const
    {
        return fractureEdgeFields.has(name);
    }

    // ---------------------------------------------------------
    // 3. Interaction Domain (NNC & FF)
    // ---------------------------------------------------------
    // 注：NNC 和 FF 数据通常使用一维标量场存储传导率

    // --- NNC (Matrix-Fracture) ---
    std::shared_ptr<volScalarField> createNNCScalar(const std::string& name, double initVal = 0.0)
    {
        return nncFields.create<volScalarField>(name, numNNCPairs, initVal);
    }

    std::shared_ptr<volScalarField> getNNCScalar(const std::string& name) const
    {
        return nncFields.get<volScalarField>(name);
    }


    // --- FF (Fracture-Fracture) ---
    // 共享 nncFields 注册表，但使用 numFFConnections 尺寸
    std::shared_ptr<volScalarField> createFFScalar(const std::string& name, double initVal = 0.0)
    {
        // 注意：这里使用的是同一个 Registry (nncFields)，但创建的场大小是 F-F 连接数
        // 只要场名称不冲突，完全可以共存
        return nncFields.create<volScalarField>(name, numFFConnections, initVal);
    }

    std::shared_ptr<volScalarField> getFFScalar(const std::string& name) const
    {
        return nncFields.get<volScalarField>(name);
    }
};