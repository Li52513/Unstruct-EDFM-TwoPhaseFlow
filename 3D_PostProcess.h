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