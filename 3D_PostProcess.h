#pragma once

#include <string>
#include <vector>
#include <memory>
#include <set>
#include <algorithm>

// 引入核心依赖
#include "3D_MeshManager.h"
#include "3D_FieldManager.h"

/**
 * @class PostProcess_3D
 * @brief 3D 后处理可视化模块 (Tecplot Exporter)
 * @details 负责将计算结果导出为 Tecplot ASCII (.dat) 格式。
 * * 核心特性：
 * 1. 自适应拓扑：自动识别混合网格类型 (Tetra/Prism/Hex)。
 * 2. 全场遍历 (Assumption-Free)：自动获取基岩与裂缝所有物理场的并集，
 * 即使两者的场定义不完全一致，也能通过补零安全输出。
 * 3. 嵌入式裂缝可视化：将裂缝面作为独立的 Zone 输出。
 */
class PostProcess_3D
{
public:
    /**
     * @brief 构造函数
     * @param meshMgr 网格管理器 (提供几何拓扑)
     * @param fieldMgr 场管理器 (提供物理量数据)
     */
    PostProcess_3D(const MeshManager_3D& meshMgr, const FieldManager_3D& fieldMgr);

    /**
     * @brief 导出全场数据到 Tecplot
     * @param filename 输出文件名 (建议以 .dat 结尾)
     * @param time 当前模拟时间 (用于 Zone Header 的 STRANDID 和 SOLUTIONTIME)
     */
    void ExportTecplot(const std::string& filename, double time = 0.0);

    // 在 PostProcess_3D 类中，与 ExportTecplot 并列
    void ExportVTK(const std::string& filename, double time = 0.0);

    /**
     * @brief [FIM 数据降维解耦] 将 3D 自动微分场 (ADVar) 同步为纯双精度标量场 (double)
     * @tparam N 独立自变量的数量 (如单相为2，两相为3)
     * @param fieldMgr 3D 场管理器引用 (非 const，因需创建/修改新场)
     * @param adFieldName 原始带梯度的自动微分场名称 (例如 "Temperature_AD")
     * @param outScalarName 降维后用于导出的纯标量场名称 (例如 "Temperature")
     * @details
     * 自动遍历基岩域与裂缝域，安全剥离 ADVar.val 物理值并拷贝至 volScalarField。
     * 配合 3D_PostProcess::GetAllUniqueFieldNames 中的 dynamic_pointer_cast 过滤器，
     * 彻底解耦 FIM 数值求解层与 Tecplot 文件 I/O 层的耦合。
     */
    template<int N>
    static void SyncADFieldToScalar(FieldManager_3D& fieldMgr, const std::string& adFieldName, const std::string& outScalarName)
    {
        // 1. 同步基岩域 (Matrix)
        auto adMat = fieldMgr.getMatrixADScalar<N>(adFieldName);
        if (adMat) {
            auto scMat = fieldMgr.getOrCreateMatrixScalar(outScalarName);
            const size_t nCells = adMat->data.size();
            for (size_t i = 0; i < nCells; ++i) {
                scMat->data[i] = adMat->data[i].val; // 根据 ADVar.hpp 解析，物理值为 val
            }
        }

        // 2. 同步裂缝域 (Fracture)
        auto adFrac = fieldMgr.getFractureADScalar<N>(adFieldName);
        if (adFrac) {
            auto scFrac = fieldMgr.getOrCreateFractureScalar(outScalarName);
            const size_t nFracCells = adFrac->data.size();
            for (size_t i = 0; i < nFracCells; ++i) {
                scFrac->data[i] = adFrac->data[i].val; // 根据 ADVar.hpp 解析，物理值为 val
            }
        }
    }

private:
    const MeshManager_3D& meshMgr_;
    const FieldManager_3D& fieldMgr_;

    // =========================================================
    // 内部辅助函数
    // =========================================================

    /**
     * @brief 获取所有域 (Matrix + Fracture) 中注册的唯一标量场名称
     * @details 取并集并排序，确保 Header 的变量列表涵盖所有可能的数据。
     * @return 已排序的变量名列表
     */
    std::vector<std::string> GetAllUniqueFieldNames() const;

    /**
     * @brief 根据节点数判定 Tecplot 单元类型
     * @param numNodes 单元顶点数
     * @return Tecplot 关键字 (e.g., "FETETRA", "FEPRISM")
     */
    std::string GetTecplotElementType(int numNodes) const;
};