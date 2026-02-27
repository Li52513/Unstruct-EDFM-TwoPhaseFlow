#pragma once

#include <string>
#include <vector>
#include <memory>
#include <set>
#include <algorithm>

// 引入核心依赖 (适配 2D 环境)
#include "MeshManager.h" 
#include "2D_FieldManager.h"

/**
 * @class PostProcess_2D
 * @brief 2D-EDFM 后处理可视化模块 (Tecplot Exporter)
 * @details 将 2D 基岩(Triangle/Quad) 和 1D 裂缝(LineSeg) 及其物理场导出为 Tecplot 格式。
 */
class PostProcess_2D
{
public:
    /**
     * @brief 构造函数
     * @param meshMgr 2D网格管理器
     * @param fieldMgr 2D场管理器
     */
    PostProcess_2D(const MeshManager& meshMgr, const FieldManager_2D& fieldMgr);

    /**
     * @brief 导出全场数据到 Tecplot
     * @param filename 输出文件名 (如 "Result_2D.dat")
     * @param time 当前模拟时间
     */
    void ExportTecplot(const std::string& filename, double time = 0.0) const;

    /**
     * @brief [FIM 数据降维解耦] 将 2D 自动微分场 (ADVar) 同步为纯双精度标量场 (double)
     * @tparam N 独立自变量的数量 (如单相为2，两相为3)
     * @param fieldMgr 2D 场管理器引用 (非 const，因需创建/修改新场)
     * @param adFieldName 原始带梯度的自动微分场名称 (例如 "Pressure_AD")
     * @param outScalarName 降维后用于导出的纯标量场名称 (例如 "Pressure")
     * @details
     * 自动遍历基岩域与裂缝域，安全剥离 ADVar.val 物理值并拷贝至 volScalarField。
     * 后处理模块 ExportTecplot 仅需检索 outScalarName 即可实现零越界风险导出。
     */
    template<int N>
    static void SyncADFieldToScalar(FieldManager_2D& fieldMgr, const std::string& adFieldName, const std::string& outScalarName)
    {
        // 1. 同步基岩域 (Matrix)
        auto adMat = fieldMgr.getMatrixADScalar<N>(adFieldName);
        if (adMat) {
            auto scMat = fieldMgr.getOrCreateMatrixScalar(outScalarName);
            const size_t nCells = adMat->data.size();
            for (size_t i = 0; i < nCells; ++i) {
                scMat->data[i] = adMat->data[i].val; // 精确访问 ADVar 成员变量 val
            }
        }

        // 2. 同步裂缝域 (Fracture)
        auto adFrac = fieldMgr.getFractureADScalar<N>(adFieldName);
        if (adFrac) {
            auto scFrac = fieldMgr.getOrCreateFractureScalar(outScalarName);
            const size_t nFracCells = adFrac->data.size();
            for (size_t i = 0; i < nFracCells; ++i) {
                scFrac->data[i] = adFrac->data[i].val; // 精确访问 ADVar 成员变量 val
            }
        }
    }

private:
    const MeshManager& meshMgr_;
    const FieldManager_2D& fieldMgr_;

    // =========================================================
    // 内部辅助函数
    // =========================================================

    /**
     * @brief 获取所有域 (2D Matrix + 1D Fracture) 中注册的唯一标量场名称
     */
    std::vector<std::string> GetAllUniqueFieldNames() const;

    /**
     * @brief 根据节点数判定 Tecplot 单元类型
     * @param numNodes 单元顶点数 (2: FELINESEG, 3: FETRIANGLE, 4: FEQUADRILATERAL)
     */
    std::string GetTecplotElementType(int numNodes) const;
};