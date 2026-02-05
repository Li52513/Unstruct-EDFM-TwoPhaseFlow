```C++
#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

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
	int id;                               ///< 裂缝全局唯一编号
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

	/**
	 * @brief 3D 求交核心接口 (预留，对应 Level 3 的生成)
	 * @details 内部调度 Type 1, Type 2, Type 3 的检测算法，生成交互多边形
	 * @param matrixMesh 基岩网格对象
	 * @param octree 基岩网格面的空间索引
	 */
	void SolveIntersection3D(const Mesh& matrixMesh, const FaceIndexedOctree& octree);

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

private:
	// =========================================================
	// 私有辅助函数
	// ========================================================

	/**
	* @brief 构建裂缝内部边的拓扑与几何
	* @param method 传递给 computeFVMVectors 的修正方法
	*/
	void _buildInternalEdges(int nU, int nV, NormalVectorCorrectionMethod method);

};

```

```C++
#include "2D_Fracture.h"
#include "Mesh.h"
#include "FaceIndexedOctree.h"

#include <iostream>
#include <numeric>
#include <fstream>
#include <iomanip>

// =========================================================
// 构造函数实现
// =========================================================

Fracture_2D::Fracture_2D(int id, const std::vector<Vector>& corners, double cond, double width)
    : id(id), boundaryVertices(corners), conductivity(cond), aperture(width)
{
}




// =========================================================
// 核心行为实现
// =========================================================
void Fracture_2D::MeshFractureSurface(int nU, int nV, NormalVectorCorrectionMethod method)
{
    // 1. 校验：仅支持 4 顶点裂缝的结构化划分
    if (boundaryVertices.size() != 4) {
        std::cerr << "[Error] Fracture_2D ID " << id
            << ": Structured meshing currently requires exactly 4 boundary vertices." << std::endl;
        return;
    }

    // 2. 清理旧数据
    fracNodes.clear();
    fracCells.clear();
    fracEdges.clear();

    // 预分配内存以提升性能
    fracNodes.reserve((nU + 1) * (nV + 1));
    fracCells.reserve(nU * nV);

    // 获取四个角点 (假设输入顺序为逆时针环绕)
    // P0(0,0) --u--> P1(1,0)
    //   |              |
    //   v              v
    // P3(0,1) --u--> P2(1,1)
    // 注意：这里的 u,v 是参数坐标，不是空间坐标
    const Vector& p0 = boundaryVertices[0]; // Bottom-Left
    const Vector& p1 = boundaryVertices[1]; // Bottom-Right
    const Vector& p2 = boundaryVertices[2]; // Top-Right
    const Vector& p3 = boundaryVertices[3]; // Top-Left

    // 3. 生成节点 (基于参数坐标 u,v 的双线性插值)
    // 外部循环 v (0 -> 1), 内部循环 u (0 -> 1)
    int nodeIDCounter = 0;
    for (int j = 0; j <= nV; ++j)
    {
        double v = (nV > 0) ? static_cast<double>(j) / nV : 0.0;
        // pLeft:  从 P0 到 P3 的插值
        // pRight: 从 P1 到 P2 的插值
        Vector pLeft = p0 * (1.0 - v) + p3 * v;
        Vector pRight = p1 * (1.0 - v) + p2 * v;
        for (int i = 0; i <= nU; ++i)
        {
            double u = (nU > 0) ? static_cast<double>(i) / nU : 0.0;

            // 在 pLeft 和 pRight 之间进行线性插值得到内部网格点
            Vector p = pLeft * (1.0 - u) + pRight * u;

            // 创建节点 (ID 从 0 开始的局部索引)
            fracNodes.emplace_back(nodeIDCounter++, p);
        }

    }

    // 4. 生成四边形单元 (构建拓扑)
    // 节点排列是行优先 (Row-major): index = j * (nU+1) + i
    int nNodeU = nU + 1; // 每一行的节点数
    int cellIDCounter = 1; // 单元 ID 通常从 1 开始

    for (int j = 0; j < nV; ++j)
    {
        for (int i = 0; i < nU; ++i)
        {
            // 计算当前网格四个角点的局部索引 (逆时针顺序)
            // n_TL (j+1, i) --- n_TR (j+1, i+1)
            //       |               |
            // n_BL (j, i)   --- n_BR (j, i+1)

            int n_BL = j * nNodeU + i;               // Bottom-Left
            int n_BR = n_BL + 1;                     // Bottom-Right
            int n_TR = (j + 1) * nNodeU + (i + 1);   // Top-Right
            int n_TL = (j + 1) * nNodeU + i;         // Top-Left

            // 创建 FractureElement_2D (Level 2 对象)
            FractureElement_2D face(cellIDCounter++, { n_BL, n_BR, n_TR, n_TL });

            // 立即计算该微观单元的几何属性 (面积、法向、包围盒)
            // 传入刚刚生成的 fracNodes 列表供其查询坐标
            face.computeGeometry(fracNodes);
            fracCells.push_back(face);  
        }
    }
    std::cout << "[Fracture_2D ID: " << id << "] Generating elements done. Now building edges..." << std::endl;

    // 5. 构建内部边 (Edges)
    _buildInternalEdges(nU, nV,method);

    std::cout << "[Fracture_2D ID: " << id << "] Meshing Complete: "
        << nU << "x" << nV << " grid, "
        << fracCells.size() << " elements, "
        << fracEdges.size() << " internal edges." << std::endl;
}

void Fracture_2D::_buildInternalEdges(int nU, int nV, NormalVectorCorrectionMethod method)
{
    fracEdges.clear();
    int edgeIDCounter = 0;
    size_t nNodeU = static_cast<size_t>(nU + 1);

    // Lambda: 安全添加边
    auto addEdge = [&](size_t c1_idx, size_t c2_idx, size_t n1_idx, size_t n2_idx)
        {
            // --- 核心修复：越界安全检查 ---
            if (c1_idx >= fracCells.size() || c2_idx >= fracCells.size()) {
                std::cerr << "[Error] Fracture_2D::_buildInternalEdges - Cell Index Out of Range! "
                    << "c1=" << c1_idx << ", c2=" << c2_idx << ", Size=" << fracCells.size() << std::endl;
                return;
            }
            if (n1_idx >= fracNodes.size() || n2_idx >= fracNodes.size()) {
                std::cerr << "[Error] Fracture_2D::_buildInternalEdges - Node Index Out of Range! " << std::endl;
                return;
            }

            // 使用 0-based 索引构造边，内部存储为 int
            FractureEdge_2D edge(edgeIDCounter++, static_cast<int>(n1_idx), static_cast<int>(n2_idx));

            // Owner/Neighbor ID = Index + 1 (FVM 惯例)
            edge.ownerCellID = static_cast<int>(c1_idx + 1);
            edge.neighborCellID = static_cast<int>(c2_idx + 1);

            // 1. 计算几何
            edge.computeGeometry(fracNodes, fracCells[c1_idx], fracCells[c2_idx]);

            // 2. 计算 FVM 向量 (使用传入的 method)
            const Vector& Cp = fracCells[c1_idx].centroid;
            const Vector& Cn = fracCells[c2_idx].centroid;
            edge.computeFVMVectors(Cp, Cn, method);

            fracEdges.push_back(edge);
        };

    // 1. 竖直方向内部边 (Vertical)
    // 强制使用 size_t 避免 int * int 溢出
    for (int j = 0; j < nV; ++j) {
        for (int i = 0; i < nU - 1; ++i) {
            // 索引计算：确保每一步都是 size_t 运算
            size_t idx_j = static_cast<size_t>(j);
            size_t idx_i = static_cast<size_t>(i);
            size_t idx_nU = static_cast<size_t>(nU);

            size_t cLeftIdx = idx_j * idx_nU + idx_i;
            size_t cRightIdx = idx_j * idx_nU + (idx_i + 1);

            size_t nBottom = idx_j * nNodeU + (idx_i + 1);
            size_t nTop = (idx_j + 1) * nNodeU + (idx_i + 1);

            addEdge(cLeftIdx, cRightIdx, nBottom, nTop);
        }
    }

    // 2. 水平方向内部边 (Horizontal)
    for (int j = 0; j < nV - 1; ++j) {
        for (int i = 0; i < nU; ++i) {
            size_t idx_j = static_cast<size_t>(j);
            size_t idx_i = static_cast<size_t>(i);
            size_t idx_nU = static_cast<size_t>(nU);

            size_t cBottomIdx = idx_j * idx_nU + idx_i;
            size_t cTopIdx = (idx_j + 1) * idx_nU + idx_i;

            size_t nLeft = (idx_j + 1) * nNodeU + idx_i;
            size_t nRight = (idx_j + 1) * nNodeU + (idx_i + 1);

            addEdge(cBottomIdx, cTopIdx, nLeft, nRight);
        }
    }
}

// 3D 求交核心接口 (框架实现)
void Fracture_2D::SolveIntersection3D(const Mesh& matrixMesh, const FaceIndexedOctree& octree)
{
    std::cout << "[Fracture_2D ID: " << id << "] Calculating 3D Intersection (Method: Octree Accelerated)..." << std::endl;

    // -------------------------------------------------------------
    // 3D-EDFM 几何求交伪代码框架 (Level 3 生成逻辑)
    // 后续将在 MeshManager 中完善具体实现
    // -------------------------------------------------------------
    /*
    // 1. 遍历裂缝自身的每一个微观单元 (Level 2)
    for (auto& fracFace : fracCells)
    {
        // 2. 利用八叉树进行 AABB 粗筛，找到可能相交的基岩面
        std::vector<int> candidateFaceIDs;
        // octree.query(fracFace.boundingBox, candidateFaceIDs);

        if (candidateFaceIDs.empty()) continue;

        // 3. 进一步获取相关的基岩单元 (Matrix Cells)
        std::set<int> candidateCellIDs;
        // ... (从 FaceID 反查 CellID)

        // 4. 对每一对 (FracElement, MatrixCell) 执行精确几何求交
        for (int mCellID : candidateCellIDs) {
             // 4.1 检测 Type 1, 2, 3 交点 (Level 5)
             // 4.2 构建交互多边形 (Level 4)
             // 4.3 计算 NNC 传导率
        }
    }
    */

    std::cout << "  -> Intersection calculation logic to be implemented in MeshManager." << std::endl;
}

// =========================================================
// 调试与 I/O 实现
// =========================================================

void Fracture_2D::printInfo() const
{
    std::cout << "----------- Fracture_2D Info (ID: " << id << ") -----------" << std::endl;

    std::cout << "  Geometry Type : Planar Polygon (Boundary Vertices: " << boundaryVertices.size() << ")" << std::endl;
    for (int i = 0; i < boundaryVertices.size(); i++)
    {
        std::cout << "  Nodes  coords : (" << boundaryVertices[i].m_x << "," << boundaryVertices[i].m_y << "," << boundaryVertices[i].m_z << ",)" << endl;
    }
    
    std::cout << "  Properties    : Cond = " << conductivity << ", Aperture = " << aperture << std::endl;
    std::cout << "  Mesh Status   : " << std::endl;
    std::cout << "    - Nodes     : " << fracNodes.size() << std::endl;
    std::cout << "    - Elements  : " << fracCells.size() << std::endl;

    // 简单的几何检查
    double totalArea = 0.0;
    for (const auto& face : fracCells) {
        totalArea += face.area;
    }
    std::cout << "    - Total Area: " << totalArea << " m^2" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
}

void Fracture_2D::exportToTxt(const std::string& filePrefix) const
{
    if (fracNodes.empty() || fracCells.empty()) {
        std::cerr << "[Warning] Fracture " << id << " has no mesh data to export." << std::endl;
        return;
    }

    // 1. 导出节点 (Nodes)
    // 文件名格式: Prefix_FracID_Nodes.txt
    std::string nodeFileName = filePrefix + "_Frac" + std::to_string(id) + "_Nodes.txt";
    std::ofstream fn(nodeFileName);

    if (!fn.is_open())
    {
        std::cerr << "[Error] Failed to open file: " << nodeFileName << std::endl;
        return;
    }
    // 设置精度
    fn << std::fixed << std::setprecision(8);

    // 格式: LocalIndex x y z (LocalIndex 用于后续 Element 的索引)
    for (size_t i = 0; i < fracNodes.size(); i++)
    {
        const Vector& p = fracNodes[i].coord;
        fn << i << " " << p.m_x << " " << p.m_y << " " << p.m_z << "\n";
    }
    fn.close();

    // 2. 导出单元 (Elements)
    std::string elemFileName= filePrefix + "_Frac" + std::to_string(id) + "_Elements.txt";
    std::ofstream fe(elemFileName);

    if (!fe.is_open())
    {
        std::cerr << "[Error] Failed to open file: " << nodeFileName << std::endl;
        return;
    }

    fe << std::fixed << std::setprecision(8);
    // 格式: ID NumNodes n1 n2 ... Area cx cy cz
    // 采用变长格式以兼容三角形或任意多边形网格
    for (const auto& face : fracCells) {
        fe << face.id << " " << face.nodeIndices.size();
        // 输出节点索引 (对应上面的 LocalIndex)
        for (int nid : face.nodeIndices) {
            fe << " " << nid;
        }
        // 输出几何属性用于检查
        fe << " " << face.area
            << " " << face.centroid.m_x
            << " " << face.centroid.m_y
            << " " << face.centroid.m_z << "\n";
    }
    fe.close();

    std::cout << "[Export] Fracture " << id << " mesh exported to:\n"
        << "         " << nodeFileName << "\n"
        << "         " << elemFileName << std::endl;
}

// =========================================================
// 新增：导出 Edge FVM 详情
// =========================================================
void Fracture_2D::exportEdgesDetailToCSV(const std::string& filePrefix) const
{
    std::string fileName = filePrefix + "_Frac" + std::to_string(id) + "_EdgeVectors.csv";
    std::ofstream fs(fileName);

    if (!fs.is_open()) {
        std::cerr << "[Error] Cannot open file " << fileName << std::endl;
        return;
    }

    // 写入表头
    fs << "EdgeID,OwnerCell,NeighborCell,Length,"
        << "Nx,Ny,Nz," // 面内法向
        << "Ex,Ey,Ez," // 正交分量
        << "Tx,Ty,Tz\n"; // 非正交分量

    fs << std::fixed << std::setprecision(8);

    for (const auto& edge : fracEdges) {
        fs << edge.id << ","
            << edge.ownerCellID << "," << edge.neighborCellID << ","
            << edge.length << ","
            << edge.normal.m_x << "," << edge.normal.m_y << "," << edge.normal.m_z << ","
            << edge.vectorE.m_x << "," << edge.vectorE.m_y << "," << edge.vectorE.m_z << ","
            << edge.vectorT.m_x << "," << edge.vectorT.m_y << "," << edge.vectorT.m_z << "\n";
    }
    fs.close();
    std::cout << "[Export] Edge vectors exported to: " << fileName << std::endl;
}

// =========================================================
// 新增：AABB 粗筛辅助函数 (用于 F-F 微观映射)
// =========================================================
std::vector<int> Fracture_2D::findCellsIntersectingBox(const AABB& queryBox) const
{
    std::vector<int> result;
    // 遍历所有单元 (如果网格量巨大，建议后续用八叉树优化；目前线性遍历即可)
    for (const auto& cell : fracCells) {
        if (cell.boundingBox.overlaps(queryBox)) {
            result.push_back(cell.id); // 注意：这里返回的是全局/局部 ID 取决于 FractureElement_2D 定义，通常建议用局部索引
        }
    }
    return result;
}


```

