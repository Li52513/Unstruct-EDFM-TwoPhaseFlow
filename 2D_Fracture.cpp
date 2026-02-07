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
// [New] 索引查找实现
// =========================================================
int Fracture_2D::getElemIndex(int elemID) const
{
    auto it = elemId2Index_.find(elemID);
    if (it != elemId2Index_.end()) {
        return it->second;
    }
    return -1; // Not found
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
    elemId2Index_.clear();

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

            // [新增] 建立索引闭环
            // -------------------------------------------------------------
            // 1. 获取当前在 vector 中的物理下标 (0-based)
            int currentLocalIndex = static_cast<int>(fracCells.size());

            // 2. 赋予单元自知身份
            face.localIndex = currentLocalIndex;
            face.parentFractureID = this->id;

            // 3. 注册到映射表 (GID -> LID)
            elemId2Index_[face.id] = currentLocalIndex;
            // -------------------------------------------------------------

            // 立即计算该微观单元的几何属性 (面积、法向、包围盒)
            // 传入刚刚生成的 fracNodes 列表供其查询坐标
            face.computeGeometry(fracNodes);
            fracCells.push_back(face);  
        }
    }
    std::cout << "[Fracture_2D ID: " << id << "] Generating elements done. Now building edges..." << std::endl;

    // 5. 构建内部边 (Edges)
    _buildInternalEdges(nU, nV,method);

    // 6. 构建边界边
    _buildBoundaryEdges(nU, nV, method); // 追加边界边

    // 7. 对边界进行重新排序
    _reorderCellEdges();


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

            edge.parentFractureID = this->id;
            edge.ownerCellLocalIndex = static_cast<int>(c1_idx); // 0-based Local Index
            edge.neighborCellLocalIndex = static_cast<int>(c2_idx); // 0-based Local Index

            // 将此边注册给 Owner
            fracCells[c1_idx].connectedEdgeIndices.push_back(edge.id);

            // 将此边注册给 Neighbor (如果存在)
            if (c2_idx < fracCells.size()) { // 确保 neighbor 也是有效单元（虽然addEgde逻辑保证了这点）
                fracCells[c2_idx].connectedEdgeIndices.push_back(edge.id);
            }

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

void Fracture_2D::_buildBoundaryEdges(int nU, int nV, NormalVectorCorrectionMethod method)
{
    // 注意：不要 clear fracEdges，我们是在追加
    // 获取当前的 Edge ID 计数器起始值
    int edgeIDCounter = static_cast<int>(fracEdges.size());

    size_t nNodeU = static_cast<size_t>(nU + 1);

    // Lambda: 添加边界边通用逻辑
    auto addBoundaryEdge = [&](size_t c_idx, size_t n1_idx, size_t n2_idx)
        {
            // 越界检查
            if (c_idx >= fracCells.size() || n1_idx >= fracNodes.size() || n2_idx >= fracNodes.size()) {
                std::cerr << "[Error] Fracture_2D::_buildBoundaryEdges - Index Out of Range!" << std::endl;
                return;
            }

            FractureEdge_2D edge(edgeIDCounter++, static_cast<int>(n1_idx), static_cast<int>(n2_idx));

            // 拓扑连接
            edge.ownerCellID = static_cast<int>(c_idx + 1); // 1-based Global ID placeholder
            edge.neighborCellID = -1;                       // -1 明确标识为边界

            edge.parentFractureID = this->id;
            edge.ownerCellLocalIndex = static_cast<int>(c_idx);
            edge.neighborCellLocalIndex = -1;

            // 注册到 Owner 单元
            fracCells[c_idx].connectedEdgeIndices.push_back(edge.id);

            // 1. 计算几何 (使用专用边界几何函数)
            edge.computeBoundaryGeometry(fracNodes, fracCells[c_idx]);

            // 2. FVM 向量 (边界处理)
            // 边界边的 d 向量通常定义为从 Owner 中心到 Edge 中点
            const Vector& Cp = fracCells[c_idx].centroid;
            const Vector& Cf = edge.midpoint; // 边界“邻居”中心即为面心

            // 调用 computeFVMVectors，但在边界情况下，Neighbor Center 传入 Edge Midpoint
            // 这将计算出 Owner 到 边界 的几何参数
            edge.computeFVMVectors(Cp, Cf, method);

            // 特殊处理：边界边的 ownerToNeighbor 向量
            // 对于边界，它等于 midpoint - owner.centroid
            edge.ownerToNeighbor = edge.midpoint - Cp;

            fracEdges.push_back(edge);
        };

    // -------------------------------------------------
    // 按顺序生成四条边界 (Left, Right, Bottom, Top)
    // -------------------------------------------------

   // 1. 左边界 (Left, i=0)
    for (int j = 0; j < nV; ++j) {
        size_t cIdx = j * nU;
        size_t n1 = j * nNodeU;
        size_t n2 = (j + 1) * nNodeU;
        addBoundaryEdge(cIdx, n2, n1);
    }

    // 2. 右边界 (Right, i=nU-1)
    for (int j = 0; j < nV; ++j) {
        size_t cIdx = j * nU + (nU - 1);
        size_t n1 = j * nNodeU + nU;
        size_t n2 = (j + 1) * nNodeU + nU;
        addBoundaryEdge(cIdx, n1, n2);
    }

    // 3. 下边界 (Bottom, j=0)
    for (int i = 0; i < nU; ++i) {
        size_t cIdx = i;
        size_t n1 = i;
        size_t n2 = i + 1;
        addBoundaryEdge(cIdx, n1, n2);
    }

    // 4. 上边界 (Top, j=nV-1)
    for (int i = 0; i < nU; ++i) {
        size_t cIdx = (nV - 1) * nU + i;
        size_t n1 = nV * nNodeU + i;
        size_t n2 = nV * nNodeU + (i + 1);
        addBoundaryEdge(cIdx, n2, n1);
    }

    // 5. 再次对边界单元进行几何重排序 (确保 connectedEdgeIndices 有序)
    // 复用之前写好的逻辑，或者直接调用一次全局重排序
    // 这里简单起见，建议在 MeshFractureSurface 统一调用一次重排序
}

void Fracture_2D::_reorderCellEdges()
{
    std::cout << "  [Fracture] Reordering ALL edges (Internal + Boundary) to match node winding order..." << std::endl;

    for (auto& cell : fracCells)
    {
        if (cell.nodeIndices.empty()) continue;
        if (cell.connectedEdgeIndices.empty()) continue;

        std::vector<int> sortedEdges;
        sortedEdges.reserve(cell.nodeIndices.size());

        for (size_t i = 0; i < cell.nodeIndices.size(); ++i)
        {
            int n1 = cell.nodeIndices[i];
            int n2 = cell.nodeIndices[(i + 1) % cell.nodeIndices.size()];

            int foundEdgeID = -1;
            for (int edgeID : cell.connectedEdgeIndices)
            {
                if (edgeID >= 0 && edgeID < static_cast<int>(fracEdges.size()))
                {
                    const auto& edge = fracEdges[edgeID];
                    bool match = (edge.nodeIndices[0] == n1 && edge.nodeIndices[1] == n2) ||
                        (edge.nodeIndices[0] == n2 && edge.nodeIndices[1] == n1);
                    if (match) {
                        foundEdgeID = edgeID;
                        break;
                    }
                }
            }
            sortedEdges.push_back(foundEdgeID);
        }
        cell.connectedEdgeIndices = sortedEdges;
    }
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

void Fracture_2D::exportToTxt_improved(const std::string& filePrefix) const
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
    std::string elemFileName = filePrefix + "_Frac" + std::to_string(id) + "_Elements.txt";
    std::ofstream fe(elemFileName);
    // ... (打开文件检查) ...

    fe << std::fixed << std::setprecision(8);

    for (const auto& face : fracCells) {
        // 基础信息
        fe << face.localIndex << " "
            << face.id << " "
            << face.solverIndex << " "
            << face.nodeIndices.size();

        for (int nid : face.nodeIndices) fe << " " << nid;

        fe << " " << face.area
            << " " << face.centroid.m_x
            << " " << face.centroid.m_y
            << " " << face.centroid.m_z;

        // [New] 直接输出预存的边信息 (Zero Logic)
        const auto& edges = face.connectedEdgeIndices;

        fe << " " << edges.size(); // 输出边数
        for (int edgeIdx : edges) {
            fe << " " << edgeIdx;
        }
        fe << "\n";
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

