#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <unordered_map>

// 引入基础类型
#include "UserDefineVarType.h" // 包含 Vector 定义及运算符重载
#include "Node.h"              // 包含 Node 定义
#include "AABB.h"              // 包含 AABB 定义

// 引入下一层级（Level 2）的定义
#include "2D_FractureMeshElement.h"
#include "2D_FractureEdge.h"
#include "Face.h"

// 前向声明，避免循环引用 (Level 3 & Mesh Manager)
class Mesh;
class FaceIndexedOctree;

// =========================================================
// 宏观裂缝对象 (Macro Fracture Surface)
// 对应层级: Level 1
// 物理含义: 定义裂缝的宏观几何，并管理生成其内部的独立网格 (Level 2)
// =========================================================
class Fracture_2D
{
public:
	// =========================================================
	// 1. 宏观几何属性 (Macro Geometry)
	// =========================================================
	int id;                               ///< 裂缝全局唯一编号 1-based
	std::vector<Vector> boundaryVertices; ///< 定义大裂缝边界的顶点 (通常4个，建议逆时针顺序)
	

	// --- 物性参数 --
	double conductivity;                  ///< 裂缝导流能力 (Fcd) 或 渗透率
	double aperture;                      ///< 平均开度 (Aperture)

	// =========================================================
	// 2. 独立网格数据 (Independent Mesh Data)
	// =========================================================
	std::vector<Node> fracNodes;                ///< 裂缝自身的节点列表 (Level 2 的顶点)
	std::vector<FractureElement_2D> fracCells;  ///< 裂缝自身的单元列表 (Level 2 的网格单元)
	std::vector<FractureEdge_2D> fracEdges;
	// =========================================================
	// 3. 构造与初始化
	// =========================================================
	
	/**
	 * @brief 构造函数
	 * @param id 裂缝编号
	 * @param corners 裂缝边界顶点 (建议 4 个点: P0, P1, P2, P3)
	 * @param cond 导流能力 (默认 1.0)
	 * @param width 开度 (默认 0.001)
	 */

	Fracture_2D(int id, const std::vector<Vector>& corners, double cond = 1.0, double width = 0.001);

	// =========================================================
	// 4. 核心行为 (Core Behaviors)
	// =========================================================

/**
	 * @brief 独立网格生成器
	 * @param nU 第一方向划分数
	 * @param nV 第二方向划分数
	 * @param method 非正交修正方法，默认 OrthogonalCorrection
	 */
	void MeshFractureSurface(int nU, int nV, NormalVectorCorrectionMethod method = NormalVectorCorrectionMethod::OrthogonalCorrection);

	// =========================================================
	// 5. 调试与 I/O (Debug & I/O)
	// =========================================================

	/**
	 * @brief 打印裂缝基本统计信息至控制台
	 */
	void printInfo() const;

	/**
	 * @brief 导出裂缝独立网格数据至 TXT 文件
	 * @details 用于 MATLAB 后处理可视化。将生成两个文件：
	 * 1. {prefix}_Frac{id}_Nodes.txt : 局部索引 x y z
	 * 2. {prefix}_Frac{id}_Elements.txt : 单元ID 节点数 n1 n2... 面积 cx cy cz
	 * @param filePrefix 文件名前缀 (例如 "Check")
	 */

	void exportToTxt(const std::string& filePrefix) const;
	void exportToTxt_improved(const std::string& filePrefix) const;

	/**
	 * @brief 导出内部边及其 FVM 向量信息至 CSV
	 * @details 输出格式: EdgeID, Owner, Neighbor, Length, Nx, Ny, Nz, Ex, Ey, Ez, Tx, Ty, Tz
	 * @param filePrefix 文件前缀
	 */
	void exportEdgesDetailToCSV(const std::string& filePrefix) const;

	/**
	 * @brief  获取与给定 AABB 相交的所有微观单元 ID
	 * @details 用于 F-F 求交时的微观映射加速
	 * @param box 目标包围盒 (通常是交线段的包围盒)
	 * @return 相交的单元局部 ID 列表
	 */
	std::vector<int> findCellsIntersectingBox(const AABB& box) const;

	// =========================================================
	// 6. 索引管理 (Index Management)
	// =========================================================
	/**
	 * @brief 通过单元 ID 获取其局部索引 (Safe Lookup)
	 * @param elemID 单元全局 ID (1-based)
	 * @return 局部索引 (0-based) in fracCells; 若未找到返回 -1
	 */
	int getElemIndex(int elemID) const;

private:
	// =========================================================
	// 私有辅助函数
	// ========================================================

	/// 单元 ID 到 局部索引 的快速查找表
	/// Key: FractureElement_2D::id (GID)
	/// Value: index in fracCells (LID)
	std::unordered_map<int, int> elemId2Index_;

	/**
	* @brief 构建裂缝内部边的拓扑与几何
	* @param method 传递给 computeFVMVectors 的修正方法
	*/
	void _buildInternalEdges(int nU, int nV, NormalVectorCorrectionMethod method);

	/**
	 * @brief 构建裂缝的物理边界边
	 * @details 处理 U 方向和 V 方向的首尾边界
	 * @param nU U方向网格数
	 * @param nV V方向网格数
	 * @param method 法向修正方法（用于 FVM 向量计算）
	 */
	void _buildBoundaryEdges(int nU, int nV, NormalVectorCorrectionMethod method);

	void _reorderCellEdges();
};
