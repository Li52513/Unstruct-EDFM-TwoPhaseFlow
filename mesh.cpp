#include "Mesh.h"
#include "gmsh.h"
#include "Node.h"
#include "Cell.h"
#include "Face.h"
#include "Fracture.h"   // 访问 Fracture::elements
#include "FracIndex.h"  // 访问 offset
#include "MeshDefinitions.h" // 包含 MeshTags 定义
#include "FaceIndexedOctree.h"
#include "EDFM_Geometry_3D.h"

 
#include <iostream>
#include <algorithm>
#include <map>
#include <gmsh.h>
#include <set>
#include <cmath>
#include <unordered_set>
#include <chrono> // 用于计时对比
#include <cmath>     // std::floor, std::abs
#include <limits>    // std::numeric_limits

// =========================================================
// 构造函数
// =========================================================
Mesh::Mesh() : gridCount_(0) {}

// =========================================================
// 1. 【基岩网格】核心数据访问接口 (Accessors)
// =========================================================
// 【2D/3D通用】 实体访问
const std::vector<Node>& Mesh::getNodes() const { return nodes_; }
const std::vector<Face>& Mesh::getFaces() const { return faces_; }
const std::vector<Cell>& Mesh::getCells() const { return cells_; }
std::vector<Cell>& Mesh::getCells() { return cells_; }                //提供非 const 引用以便求解器修改状态

// 【适用3D】 获取指定单元的所有棱
std::vector<MatrixEdge> Mesh::getCellEdges(int cellIndex) const
{
    std::vector<MatrixEdge> edges;
    const Cell& cell = cells_[cellIndex];
    const auto& cNodes = cell.CellNodeIDs;
    size_t n = cNodes.size();

    // 根据单元类型生成局部棱 (逻辑与 getLocalFaces_3D 类似)
    // 4=Tet, 5=Hex, 6=Prism, 7=Pyramid
    if (n == 4) { // Tetra
        // 6 edges: (0,1), (1,2), (2,0), (0,3), (1,3), (2,3)
        edges = { {cNodes[0],cNodes[1]}, {cNodes[1],cNodes[2]}, {cNodes[2],cNodes[0]},
                  {cNodes[0],cNodes[3]}, {cNodes[1],cNodes[3]}, {cNodes[2],cNodes[3]} };
    }
    else if (n == 8) { // Hex
        // 12 edges: 底面4, 顶面4, 侧棱4
        edges = { {cNodes[0],cNodes[1]}, {cNodes[1],cNodes[2]}, {cNodes[2],cNodes[3]}, {cNodes[3],cNodes[0]}, // Bottom
                  {cNodes[4],cNodes[5]}, {cNodes[5],cNodes[6]}, {cNodes[6],cNodes[7]}, {cNodes[7],cNodes[4]}, // Top
                  {cNodes[0],cNodes[4]}, {cNodes[1],cNodes[5]}, {cNodes[2],cNodes[6]}, {cNodes[3],cNodes[7]} }; // Side
    }
    else if (n == 6) { // Prism
        // 9 edges: 底3, 顶3, 侧3
        edges = { {cNodes[0],cNodes[1]}, {cNodes[1],cNodes[2]}, {cNodes[2],cNodes[0]}, // Bottom
                  {cNodes[3],cNodes[4]}, {cNodes[4],cNodes[5]}, {cNodes[5],cNodes[3]}, // Top
                  {cNodes[0],cNodes[3]}, {cNodes[1],cNodes[4]}, {cNodes[2],cNodes[5]} }; // Side
    }
    else if (n == 5) { // Pyramid
        // 8 edges: 底4, 侧4
        edges = { {cNodes[0],cNodes[1]}, {cNodes[1],cNodes[2]}, {cNodes[2],cNodes[3]}, {cNodes[3],cNodes[0]}, // Base
                  {cNodes[0],cNodes[4]}, {cNodes[1],cNodes[4]}, {cNodes[2],cNodes[4]}, {cNodes[3],cNodes[4]} }; // Side
    }

    return edges;

}
// 【适用3D】 提取全局唯一的基岩棱 (Matrix Edges)
std::vector<MatrixEdge> Mesh::extractUniqueMatrixEdges() const
{
    std::cout << "[EDFM] Building unique matrix edges for intersection detection..." << std::endl;

    // 使用 HashSet 去重
    std::unordered_set<MatrixEdge, MatrixEdgeHash> uniqueSet;
    // 预估容量：大约是 Cell 数的 1.5 到 2 倍（非结构化网格经验值）
    uniqueSet.reserve(cells_.size() * 2);

    for (size_t i = 0; i < cells_.size(); ++i) {
        std::vector<MatrixEdge> localEdges = getCellEdges(static_cast<int>(i));
        for (const auto& e : localEdges) {
            uniqueSet.insert(e);
        }
    }

    // 转为 Vector 并赋予 ID
    std::vector<MatrixEdge> result;
    result.reserve(uniqueSet.size());
    int edgeID = 0;
    for (const auto& e : uniqueSet) {
        result.emplace_back(e.n1, e.n2, edgeID++);
    }

    std::cout << "[EDFM] Extracted " << result.size() << " unique matrix edges." << std::endl;
    return result;
}

// 【2D/3D通用】 获取私有成员节点、面、单元索引映射表的外部访问
const std::unordered_map<int, Node>& Mesh::getNodesMap() const { return nodesMap_; }          // 经_fetchNodesFromGmsh()构建的 Node ID 到 Node 对象的映射，key值与 Gmsh nodeTag 一致，1-based
const std::unordered_map<int, int>& Mesh::getCellId2Index() const { return cellId2index_; }   // 经_fetchCells2D()/3D()构建的 Cell Global ID（1based) 到 Local Vector Index(0-based) 的映射，key值与 Gmsh elementTag 一致，1-based
const std::unordered_map<int, int>& Mesh::getFaceId2Index()const { return faceId2Index_; }    // 经 buildFaceIndexMap构建的 Face Global ID（1based) 到 Local Vector Index(0-based) 的映射

// 【2D/3D通用】 元数据访问 
int Mesh::getGridCount() const { return gridCount_; } // 获取网格总数
std::string Mesh::getMatrixElementType() const { return matrixElementType_; }


// =========================================================
// 2. 【基岩网格】 ID 转换安全接口 (Safe ID Lookups)
// =========================================================

/**
 * @brief 【2D/3D通用】通过 Cell Global ID (Gmsh Tag) 获取 Local Vector Index
 * @return 成功返回 index [0, N-1], 失败返回 -1
 */
int Mesh::getCellIndex(int cellGlobalID) const
{
    auto it = cellId2index_.find(cellGlobalID);
    if (it != cellId2index_.end()) {
        return it->second;
    }
    return -1; // Not found
}

/**
 * @brief 【2D/3D通用】 通过 Face Global ID 获取 Local Vector Index
 * @return 成功返回 index [0, N-1], 失败返回 -1
 */
int Mesh::getFaceIndex(int globalFaceID) const
{
    auto it = faceId2Index_.find(globalFaceID);
    if (it != faceId2Index_.end()) {
        return it->second;
    }
    std::cerr << "[Error] Face ID " << globalFaceID << " not found in map!" << std::endl;
    return -1;
}

// =========================================================
// 3. 【基岩网格】 网格生成 (Generation)
// =========================================================
/**
 * @brief: 【适用2D EDFM】生成 2D 非结构化基岩网格
 * 流程: Gmsh生成 -> 读取节点 -> 构建2D单元 -> 构建边(Face)拓扑 -> 映射物理边界 -> 构建网格face的Bin索引 -> 区分内部/边界单元 -> 填充face的owner neighbor 单元索引
 */
void Mesh::GenerateMesh2D(double lengthX, double lengthY, int nx, int ny, bool useQuadBase, std::string jobName)
{
    std::cout << "=========================================================" << std::endl;
    std::cout << "[Mesh] Generating 2D Mesh: " << jobName << "..." << std::endl;
    std::cout << "       Domain: " << lengthX << " x " << lengthY << " x "<< std::endl;
    std::cout << "       Grid  : " << nx << " x " << ny << std::endl;

    // 1. 初始化 Gmsh
    gmsh::initialize();
    gmsh::model::add(jobName);

    double lc = std::min(lengthX / std::max(nx, 1), lengthY / std::max(ny, 1));

    // 2. 几何建模 (Z=0 平面)
    int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
    int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
    int p3 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
    int p4 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);

    int l1 = gmsh::model::geo::addLine(p1, p2); // Bottom (y=0)
    int l2 = gmsh::model::geo::addLine(p2, p3); // Right (x=L)
    int l3 = gmsh::model::geo::addLine(p3, p4); // Top (y=H)
    int l4 = gmsh::model::geo::addLine(p4, p1); // Left (x=0)

    int loop = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 });
    int surface = gmsh::model::geo::addPlaneSurface({ loop });

    // 3. 标记物理组
    // 1D 物理组用于边界标记
    gmsh::model::geo::addPhysicalGroup(1, { l4 }, MeshTags::LEFT, "Left");
    gmsh::model::geo::addPhysicalGroup(1, { l2 }, MeshTags::RIGHT, "Right");
    gmsh::model::geo::addPhysicalGroup(1, { l1 }, MeshTags::BOTTOM, "Bottom");
    gmsh::model::geo::addPhysicalGroup(1, { l3 }, MeshTags::TOP, "Top");
    // 2D 物理组用于流体域
    gmsh::model::geo::addPhysicalGroup(2, { surface }, MeshTags::FLUID, "FluidDomain");

    // 4. 网格参数设置
    if (useQuadBase)
    {
        gmsh::option::setNumber("Mesh.RecombineAll", 1);
        gmsh::option::setNumber("Mesh.Algorithm", 8); // Frontal-Delaunay for Quads
        gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2); // Blossom
    }
    else 
    {
        gmsh::option::setNumber("Mesh.RecombineAll", 0);
        gmsh::option::setNumber("Mesh.Algorithm", 6); // Frontal-Delaunay
    }

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);

    // 5. 执行流水线 (Pipeline)
    _fetchNodesFromGmsh(); // 1. 读取节点
    _fetchCells2D();       // 2. 读取单元
    _buildTopology2D();    // 3. 构建拓扑 & 映射边界 Tag
	buildFaceGlobalId2LocalIndexMap(); // 4. 构建 Face ID 到 Index 映射表

    // 6. 后处理
    buildFaceBins_2D();         // 2D 需要空间索引
    ClassifySolidMatrixCells(); // 区分内部/边界单元

    _calcuteOwnerandNeighborCell_indexofFace();

    gmsh::write(jobName + "_.msh");

    gmsh::finalize();
    std::cout << "[Mesh] 2D Generation Complete. Nodes: " << nodes_.size()
        << ", Cells: " << cells_.size() << ", Faces: " << faces_.size() << std::endl;
}

/**
 * @brief: 【适用3D EDFM】生成 3D 非结构化基岩网格
 * 流程: Gmsh生成 -> 映射物理边界-> 读取节点 -> 构建3D单元 -> 构建边(Face)拓扑 ->  区分内部/边界单元-> 构建网格face的Bin索引 -> 填充face的owner neighbor 单元索引 ->提取基岩棱 (Matrix Edges)-> 构建空间索引 (Octree)
 */
void Mesh::GenerateMesh3D(double lengthX, double lengthY, double lengthZ, int nx, int ny, int nz, bool usePrism, bool useQuadBase, std::string jobName)
{
    std::cout << "=========================================================" << std::endl;
    std::cout << "[Mesh] Starting 3D Mesh Generation: " << jobName << std::endl;
    std::cout << "       Domain: " << lengthX << " x " << lengthY << " x " << lengthZ << std::endl;
    std::cout << "       Grid  : " << nx << " x " << ny << " x " << nz << std::endl;

    // 1. [关键修改] 预先确定并记录网格类型
    // 这是解决 Tecplot 报错的第一步：明确我们期望什么类型的网格
    if (usePrism) {
        if (useQuadBase) {
            matrixElementType_ = "FEBRICK"; // Hex (8 nodes)
        }
        else {
            matrixElementType_ = "FEPRISM"; // Prism (6 nodes)
        }
    }
    else {
        matrixElementType_ = "FETETRA";     // Tet (4 nodes)
    }
    std::cout << " -> Target Element Type: " << matrixElementType_ << std::endl;


    // 1. Gmsh 初始化与几何建模 (Extrude 流程)
    gmsh::initialize();
    gmsh::model::add(jobName);

    double lc = std::min(lengthX / std::max(nx, 1), lengthY / std::max(ny, 1));

    // Base Points (z=0)
    int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
    int p2 = gmsh::model::geo::addPoint(lengthX, 0, 0, lc);
    int p3 = gmsh::model::geo::addPoint(lengthX, lengthY, 0, lc);
    int p4 = gmsh::model::geo::addPoint(0, lengthY, 0, lc);

    // Base Lines
    int l1 = gmsh::model::geo::addLine(p1, p2);
    int l2 = gmsh::model::geo::addLine(p2, p3);
    int l3 = gmsh::model::geo::addLine(p3, p4);
    int l4 = gmsh::model::geo::addLine(p4, p1);

    int loop = gmsh::model::geo::addCurveLoop({ l1, l2, l3, l4 });
    int baseSurface = gmsh::model::geo::addPlaneSurface({ loop });

    // 网格算法设置
    if (useQuadBase) {
        gmsh::option::setNumber("Mesh.RecombineAll", 1);
        gmsh::option::setNumber("Mesh.Algorithm", 8);    // Frontal-Delaunay for Quads
        gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2);
    }
    else {
        gmsh::option::setNumber("Mesh.RecombineAll", 0);
        gmsh::option::setNumber("Mesh.Algorithm", 6);    // Frontal-Delaunay
    }

    // 3D 扫掠 (Extrude)
    std::vector<std::pair<int, int>> extrudeIn = { {2, baseSurface} };
    std::vector<std::pair<int, int>> extrudeOut;
    std::vector<int> numLayers = { std::max(nz, 1) };
    std::vector<double> heights = { 1.0 };
    bool recombine = true; // 尝试合并以生成棱柱或六面体

    // 执行 Extrude
    gmsh::model::geo::extrude(extrudeIn, 0, 0, lengthZ, extrudeOut, numLayers, heights, recombine);

    // 同步几何实体
    gmsh::model::geo::synchronize();

    // =========================================================
    // 2. 物理组标记 (Physical Groups) - 完整六面标记
    // =========================================================
    // 确保 extrudeOut 大小符合预期 (1 Top + 1 Vol + 4 Sides = 6 elements)
    if (extrudeOut.size() >= 6)
    {
        // 1. 标记 Volume (体域)
        int volTag = extrudeOut[1].second;
        gmsh::model::addPhysicalGroup(3, { volTag }, MeshTags::FLUID, "FluidVolume");

        // 2. 标记 Front / Back (Z 方向)
        // Base Surface (Z=0) 是我们手动创建的 baseSurface
        gmsh::model::addPhysicalGroup(2, { baseSurface }, MeshTags::TAG_FRONT, "Z0_Face");
        // Top Surface (Z=L) 是 extrudeOut[0]
        int topTag = extrudeOut[0].second;
        gmsh::model::addPhysicalGroup(2, { topTag }, MeshTags::TAG_BACK, "ZL_Face");

        // 3. 标记 Sides (XY平面四周)
        // 顺序对应 loop = {l1, l2, l3, l4}
        int sideBottom = extrudeOut[2].second; // 对应 l1 (y=0)
        int sideRight = extrudeOut[3].second; // 对应 l2 (x=L)
        int sideTop = extrudeOut[4].second; // 对应 l3 (y=L)
        int sideLeft = extrudeOut[5].second; // 对应 l4 (x=0)

        gmsh::model::addPhysicalGroup(2, { sideBottom }, MeshTags::BOTTOM, "Y0_Face");
        gmsh::model::addPhysicalGroup(2, { sideRight }, MeshTags::RIGHT, "XL_Face");
        gmsh::model::addPhysicalGroup(2, { sideTop }, MeshTags::TOP, "YL_Face");
        gmsh::model::addPhysicalGroup(2, { sideLeft }, MeshTags::LEFT, "X0_Face");

        std::cout << "[Mesh] Physical Groups tagged: FLUID, FRONT(Z0), BACK(ZL), BOTTOM(Y0), RIGHT(XL), TOP(YL), LEFT(X0)" << std::endl;
    }
    else
    {
        std::cerr << "[Error] Extrude did not return expected number of entities. Physical groups may be incomplete." << std::endl;
    }

    // 生成网格
    gmsh::model::mesh::generate(3);

    // 2. 数据读取与拓扑构建 (流水线)
    std::cout << "[Mesh] Fetching Nodes..." << std::endl;
    _fetchNodesFromGmsh();

    std::cout << "[Mesh] Building 3D Cells..." << std::endl;
    _fetchCells3D();

    std::cout << "[Mesh] Building Face Topology..." << std::endl;
    _buildTopology3D();

    buildFaceGlobalId2LocalIndexMap();

    std::cout << "[Mesh] Classifying Matrix Cells..." << std::endl;
    ClassifySolidMatrixCells();

    _calcuteOwnerandNeighborCell_indexofFace();

    // =========================================================
    // 3. EDFM 预处理核心：提取基岩棱 (Matrix Edges)
    // =========================================================
    std::cout << "[Mesh] Extracting Global Unique Matrix Edges (for Type-3 Intersection)..." << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<MatrixEdge> globalEdges = extractUniqueMatrixEdges();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "       -> Extracted " << globalEdges.size() << " unique edges." << std::endl;
    std::cout << "       -> Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;

    // =========================================================
    // 4. 构建空间索引 (Octree) 与 效率对比测试桩
    // =========================================================
    std::cout << "[Mesh] Building Face Indexed Octree..." << std::endl;
    t1 = std::chrono::high_resolution_clock::now();

    FaceIndexedOctree octree;
    octree.build(faces_);

    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "       -> Octree Build Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " ms" << std::endl;

    // ---------------------------------------------------------
    // 5. 拓扑结构自检 (Topology Verification)
    // ---------------------------------------------------------
    std::cout << "\n---------------------------------------------------------" << std::endl;
    std::cout << "          3D Mesh Topology Verification Report           " << std::endl;
    std::cout << "---------------------------------------------------------" << std::endl;

    long nV = nodes_.size();
    long nE = globalEdges.size();
    long nF = faces_.size();
    long nC = cells_.size();

    // 欧拉示性数 (Euler Characteristic): V - E + F - C
    // 理论值：对于单连通的简单多面体区域，通常为 1。
    // 注意：如果有孔洞或非流形结构，值会变化。此指标主要用于监测异常巨大的偏差。
    long euler = nV - nE + nF - nC;

    std::cout << "1. Entity Counts:" << std::endl;
    std::cout << "   - Nodes (V) : " << nV << std::endl;
    std::cout << "   - Edges (E) : " << nE << std::endl;
    std::cout << "   - Faces (F) : " << nF << std::endl;
    std::cout << "   - Cells (C) : " << nC << std::endl;
    std::cout << "2. Euler (V-E+F-C) : " << euler << " (Expected approx. 1)" << std::endl;

    // 几何完整性检查
    double totalVol = 0.0;
    for (const auto& c : cells_) totalVol += c.volume;
    double expectedVol = lengthX * lengthY * lengthZ;
    double volErr = std::abs(totalVol - expectedVol);

    std::cout << "3. Volume Consistency:" << std::endl;
    std::cout << "   - Total Volume : " << totalVol << std::endl;
    std::cout << "   - Expected     : " << expectedVol << std::endl;
    std::cout << "   - Error        : " << volErr << (volErr < 1e-6 ? " [OK]" : " [WARNING]") << std::endl;

    std::cout << "---------------------------------------------------------" << std::endl;

    gmsh::write("UnstructuredMesh-EDFM.msh");

    // 6. 清理
    gmsh::finalize();
}

// =========================================================
// 4. 【基岩网格】 拓扑关系与网格面几何信息计算 (Topology)
// =========================================================

/**
 * @brief 【2D/3D通用】 构建 Face Global ID -> Local Index 的映射表
 * @details 必须GenerateMesh*D之后即在网格生成完毕后调用
 */
void Mesh::buildFaceGlobalId2LocalIndexMap()
{
    faceId2Index_.clear();
    for (size_t i = 0; i < faces_.size(); ++i) {
        faceId2Index_[faces_[i].id] = static_cast<int>(i);
    }
    std::cout << "[Mesh] Built face index map for " << faces_.size() << " faces." << std::endl;
}

/**
 * @brief 【2D/3D通用】对基岩单元进行分类 (Inner / Boundary)
 * @details 在GenerateMesh*D函数中调用，实现区分 Inner (内部) 和 Boundary (边界) 单元
 */
void Mesh::ClassifySolidMatrixCells()
{
    // 构造一个 faceId 到 Face 的映射

    for (auto& cell : cells_)
    {
        bool isBoundary = false;
        for (int faceId : cell.CellFaceIDs)
        {
            // 转换 ID 到索引 (假设 1-based ID 对应 0-based Index)
            int faceIndex = getFaceIndex(faceId);
            //int faceIndex = faceId - 1;

            if (faceIndex >= 0 && faceIndex < faces_.size()) {
                if (faces_[faceIndex].neighborCell == -1) {
                    isBoundary = true;
                    break;
                }
            }
        }
        cell.location = isBoundary ? Cell::LocationType::Boundary : Cell::LocationType::Inner;
    }
}

/**
 * @brief 【2D/3D通用】 计算基岩网格面的几何信息
 * @param method 法向向量修正方法 (最小修正/正交修正/超松弛)
 * @details 计算面法向、面积矢量、正交(E)与非正交(T)分量
 */
void Mesh::ComputeMatrixMeshFaceGeometricInfor(NormalVectorCorrectionMethod method)
{
    for (auto& face : faces_)
    {
        // [优化] 直接使用 index 访问，不再查表
        if (face.ownerCell_index == -1) continue;

        const Vector& Cp = cells_[face.ownerCell_index].center;

        if (face.neighborCell_index != -1)
        {
            const Vector& Cn = cells_[face.neighborCell_index].center;
            face.computeFaceVectors(Cp, Cn, method);
        }
        else
        {
            face.computeFaceVectorsBoundary(Cp, method);
        }
    }
}

// =========================================================
// 5.1 【2D-EDFM】 2D基岩网格与1D裂缝快速求交算法
// =========================================================

/**
 * @brief 【适用2D EDFM】 为基岩网格整体构建第一次AABB索引
 * @details 为利用背景网格索引快速获取候选面建立基础，即Fracture::DetectFracturetoMeshFaceIntersections中的方法3
 */
 // 构建均匀背景网格索引，建立数组/哈希表逻辑（直接映射） 且支持 DDA 快速遍历
void Mesh::buildFaceBins_2D()
{
    if (faces_.empty()) return;

    // 1. 计算全局包围盒 (一次性计算，避免查询时重复遍历)
    globalMinX_ = std::numeric_limits<double>::max();
    globalMinY_ = std::numeric_limits<double>::max();
    double maxX = -std::numeric_limits<double>::max();
    double maxY = -std::numeric_limits<double>::max();

    for (const auto& f : faces_) {
        globalMinX_ = std::min(globalMinX_, f.boundingBox.min.m_x);  //找到几何域的左下角x
        globalMinY_ = std::min(globalMinY_, f.boundingBox.min.m_y);  //找到几何域的左下角y
        maxX = std::max(maxX, f.boundingBox.max.m_x);                //找到几何域的右上角x
        maxY = std::max(maxY, f.boundingBox.max.m_y);                //找到几何域的右上角y
    }

    // 2. 自动计算 Bin 数量
    int targetFacesPerBin = TARGET_FACES_PER_BIN;
    int totalFaces = static_cast<int>(faces_.size());
    int totalBins = std::max(1, totalFaces / targetFacesPerBin);
    int sideBins = static_cast<int>(std::sqrt(totalBins));

    // 限制网格规模，防止过细 (Clamp 5 ~ 200)
    binCountX_ = std::min(std::max(sideBins, 5), 1000);
    binCountY_ = std::min(std::max(sideBins, 5), 1000);

    // 3. 计算步长 & 鲁棒性保护
    double width = maxX - globalMinX_;
    double height = maxY - globalMinY_;

    // 防止宽度为0 (例如 1D 网格)
    if (width < 1e-9) width = 1.0;
    if (height < 1e-9) height = 1.0;

    binSizeX_ = width / binCountX_;
    binSizeY_ = height / binCountY_;

    // 4. 清空并填充 Bins
    faceBins_.clear();

    ///为每个面片建立空间索引卡片
    for (const auto& f : faces_) {
        const AABB& box = f.boundingBox; //取出面的box

        int ix_min = static_cast<int>((box.min.m_x - globalMinX_) / binSizeX_);
        int iy_min = static_cast<int>((box.min.m_y - globalMinY_) / binSizeY_);
        int ix_max = static_cast<int>((box.max.m_x - globalMinX_) / binSizeX_);
        int iy_max = static_cast<int>((box.max.m_y - globalMinY_) / binSizeY_);

        // 边界钳制
        ix_min = std::max(0, std::min(binCountX_ - 1, ix_min));
        iy_min = std::max(0, std::min(binCountY_ - 1, iy_min));
        ix_max = std::max(0, std::min(binCountX_ - 1, ix_max));
        iy_max = std::max(0, std::min(binCountY_ - 1, iy_max));

        for (int ix = ix_min; ix <= ix_max; ++ix) {
            for (int iy = iy_min; iy <= iy_max; ++iy) {
                int binID = computeBinID(ix, iy);
                faceBins_[binID].push_back(f.id);
            }
        }
    }
    // ---------- 5. 日志 ----------
    std::size_t filledBins = 0;
    for (const auto& kv : faceBins_) if (!kv.second.empty()) ++filledBins;

    std::cout << "[buildFaceBins] 总面数 = " << faces_.size()
        << ", Bin 划分 = " << binCountX_ << "×" << binCountY_
        << ", 平均每 bin 面数 ≈ "
        << static_cast<double>(faces_.size()) / std::max<std::size_t>(1, filledBins)
        << std::endl;
}

/**
 * @brief  【适用2D EDFM】 通过2D基岩网格面（1D）face与待几何离散的裂缝的AABB包围盒，判断是否需要进行下一步（二次AABB求交或直接精准计算是否相交）
 * @param   待几何离散的裂缝的AABB包围盒
 * @details 由于在2D情况，裂缝网格是通过1D宏观裂缝与2D基岩网格face的交点、裂缝的起终点以及F-F交点进行网格划分的
 * 这里裂缝的AABB包围盒，经buildFaceBins_2D构建基岩网格面的AABB后这里判断两AABB是否重合，形成需要第二次粗筛的基岩网格面候选面 ID
 */
 //构建与输入的AABB需要第二次粗筛的网格面候选面 ID 列表 method 3 核心
std::vector<int> Mesh::getCandidateFacesFromBins(const AABB& box) const
{
    std::vector<int> candidates;
    if (faceBins_.empty()) return candidates;

    // 1. 计算查询输入 AABB 覆盖的 Bin 范围
    int ix_min = static_cast<int>((box.min.m_x - globalMinX_) / binSizeX_); //计算输入的AABB盒子的最左端位于第几列
    int iy_min = static_cast<int>((box.min.m_y - globalMinY_) / binSizeY_); //计算输入的AABB盒子的最下端位于第几行
    int ix_max = static_cast<int>((box.max.m_x - globalMinX_) / binSizeX_); //计算输入的AABB盒子的最右端位于第几列
    int iy_max = static_cast<int>((box.max.m_y - globalMinY_) / binSizeY_); //计算输入的AABB盒子的最上端位于第几行

    // 2. 快速剔除：如果完全在网格域外
    if (ix_max < 0 || iy_max < 0 || ix_min >= binCountX_ || iy_min >= binCountY_) {
        return candidates;
    }

    // 3. 边界钳制
    ix_min = std::max(0, std::min(binCountX_ - 1, ix_min));
    iy_min = std::max(0, std::min(binCountY_ - 1, iy_min));
    ix_max = std::max(0, std::min(binCountX_ - 1, ix_max));
    iy_max = std::max(0, std::min(binCountY_ - 1, iy_max));

    // 4. 收集候选 ID
    // 预估容量，减少内存重分配
    candidates.reserve(50);

    for (int ix = ix_min; ix <= ix_max; ++ix) {
        for (int iy = iy_min; iy <= iy_max; ++iy) {
            int binID = computeBinID(ix, iy);
            auto it = faceBins_.find(binID);
            if (it != faceBins_.end()) {
                // 将该 Bin 内的所有面加入
                candidates.insert(candidates.end(), it->second.begin(), it->second.end());
            }
        }
    }

    // 5. 排序与去重 (因为一个面可能跨越多个 Bin)
    if (!candidates.empty()) {
        std::sort(candidates.begin(), candidates.end());
        auto last = std::unique(candidates.begin(), candidates.end());
        candidates.erase(last, candidates.end());
    }
    return candidates;
}

/**
 * @brief 【适用2D EDFM】 基于 DDA (Digital Differential Analyzer) 算法的裂缝候选面筛选
 * @param p1 裂缝线段起点
 * @param p2 裂缝线段终点
 * @return 潜在相交的基岩网格面 ID 列表
 * @details
 * 相比于 getCandidateFacesFromBins 的矩形区域查询，本方法利用 DDA 算法
 * 仅遍历线段实际经过的网格 Bin (Rasterization)，将搜索复杂度从 O(N*M) 降低至 O(max(N,M))。
 * 这极大地减少了对角线裂缝带来的“虚空”候选面。
 */
 //  DDA 算法实现
std::vector<int> Mesh::getCandidateFacesFromBins_DDA(const Vector& p1, const Vector& p2) const
{
    std::vector<int> candidates;

    // 0. 基础检查：若索引未构建，直接返回空
    if (faceBins_.empty()) return candidates;

    // 1. 坐标归一化：映射到 Bin 的局部浮点坐标系
    // p = (world_pos - origin) / bin_size
    // 注意：需保护 binSize 为 0 的情况（虽然 buildFaceBins_2D 已有保护，但此处做防御性编程）
    double invSizeX = (binSizeX_ > 1e-12) ? (1.0 / binSizeX_) : 0.0;
    double invSizeY = (binSizeY_ > 1e-12) ? (1.0 / binSizeY_) : 0.0;

    double startX = (p1.m_x - globalMinX_) * invSizeX;
    double startY = (p1.m_y - globalMinY_) * invSizeY;
    double endX = (p2.m_x - globalMinX_) * invSizeX;
    double endY = (p2.m_y - globalMinY_) * invSizeY;

    // 2. 计算离散网格坐标 (Grid Integer Coordinates)
    int x = static_cast<int>(std::floor(startX));
    int y = static_cast<int>(std::floor(startY));
    int xEnd = static_cast<int>(std::floor(endX));
    int yEnd = static_cast<int>(std::floor(endY));

    // 3. DDA 步进参数初始化
    double dx = endX - startX;
    double dy = endY - startY;

    // 步进方向 (Step Direction)
    int stepX = (dx >= 0) ? 1 : -1;
    int stepY = (dy >= 0) ? 1 : -1;

    // 步进单位增量 (tDelta): 移动一个网格单位所需的归一化步长
    // 避免除零：若 dx 为 0，则 tDeltaX 设为无穷大
    double tDeltaX = (std::abs(dx) > 1e-12) ? std::abs(1.0 / dx) : std::numeric_limits<double>::max();
    double tDeltaY = (std::abs(dy) > 1e-12) ? std::abs(1.0 / dy) : std::numeric_limits<double>::max();

    // 初始步进距离 (tMax): 从起点到达下一个网格边界所需的归一化步长
    // 公式说明：
    // 若向正向移动：(floor(start) + 1 - start) * tDelta
    // 若向负向移动：(start - floor(start)) * tDelta
    double distX = (stepX > 0) ? (std::floor(startX) + 1.0 - startX) : (startX - std::floor(startX));
    double distY = (stepY > 0) ? (std::floor(startY) + 1.0 - startY) : (startY - std::floor(startY));

    double tMaxX = (std::abs(dx) > 1e-12) ? (distX * tDeltaX) : std::numeric_limits<double>::max();
    double tMaxY = (std::abs(dy) > 1e-12) ? (distY * tDeltaY) : std::numeric_limits<double>::max();

    // 4. 辅助 Lambda：边界检查与候选收集
    // 仅当 (ix, iy) 在有效网格范围内时才查询
    auto collectFacesInBin = [&](int ix, int iy)
        {
            if (ix >= 0 && ix < binCountX_ && iy >= 0 && iy < binCountY_)
            {
                int binID = computeBinID(ix, iy);
                auto it = faceBins_.find(binID);
                if (it != faceBins_.end())
                {
                    candidates.insert(candidates.end(), it->second.begin(), it->second.end());
                }
            }
        };

    // [优化] 内存预留策略
    // 理论容量 = 穿过的 Bin 数量 * 平均每个 Bin 的面数
    // 使用 TARGET_FACES_PER_BIN 替代硬编码的 10，保证策略一致性
    // 额外 +1 是为了安全冗余
    int estimatedBins = std::abs(xEnd - x) + std::abs(yEnd - y) + 1;
    candidates.reserve(estimatedBins * TARGET_FACES_PER_BIN);

    // 5. DDA 迭代主循环
    // 最大迭代次数保护：防止浮点误差导致的死循环
    // 理论最大步数 = binCountX + binCountY
    int maxSteps = binCountX_ + binCountY_ + 10;

    for (int i = 0; i < maxSteps; ++i)
    {
        // 5.1 收集当前 Bin 的面
        collectFacesInBin(x, y);

        // 5.2 终止条件：已到达终点所在的 Bin
        if (x == xEnd && y == yEnd) break;

        // 5.3 迈向下一个 Bin
        // 比较 tMaxX 和 tMaxY，谁小就沿谁的方向走一步
        if (tMaxX < tMaxY)
        {
            tMaxX += tDeltaX;
            x += stepX;
        }
        else
        {
            tMaxY += tDeltaY;
            y += stepY;
        }
    }

    // 6. 后处理：排序与去重
    // 必须步骤：因为同一个大面可能跨越多个相邻的 Bin，会被多次添加
    if (!candidates.empty())
    {
        std::sort(candidates.begin(), candidates.end());
        auto last = std::unique(candidates.begin(), candidates.end());
        candidates.erase(last, candidates.end());
    }

    return candidates;
}


// =========================================================
// 5.2 【3D-EDFM】 3D基岩网格与2D裂缝快速求交算法
// =========================================================

// @brief【适用3D EDFM】 触发所有 Cell 的 AABB 计算，计算所有Cell的AABB包围盒,为后续构建背景网格索引做准备
void Mesh::computeAllCellAABBs()
{
    std::cout << "[Mesh] Computing AABB for all " << cells_.size() << " cells..." << std::endl;
    for (auto& cell : cells_)
    {
        cell.computeAABB(nodesMap_);
    }
}

/**
 * @brief 【适用3D EDFM】构建 3D 均匀背景网格索引 (Binning)
 * @details 将所有基岩单元根据其 AABB 注册到覆盖的 Bin 中
 * @param resolutionX X 方向期望划分的份数
 * @param resolutionY Y 方向期望划分的份数
 * @param resolutionZ Z 方向期望划分的份数
 */
void Mesh::buildCellBins(int resolutionX, int resolutionY, int resolutionZ)
{
    std::cout << "[Mesh] Building 3D Uniform Background Grid ("
        << resolutionX << "x" << resolutionY << "x" << resolutionZ << ")..." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // 1. 重置旧数据
    cellBins_.clear();
    gridParams_.reset();

    if (cells_.empty()) {
        std::cerr << "[Warning] Empty mesh cells, skipping background grid." << std::endl;
        return;
    }

    // 2. 计算全局包围盒 (Global Mesh AABB)
    // 虽然 cellId2index_ 构建时可能算过，这里为了独立性重新快速计算
    // 既然 cells 都有 AABB 了，直接用 Cell AABB 的并集更快
    Vector globalMin(1e30, 1e30, 1e30);
    Vector globalMax(-1e30, -1e30, -1e30);

    for (const auto& cell : cells_) {
        const AABB& b = cell.boundingBox;
        // 简单合并 min/max
        if (b.min.m_x < globalMin.m_x) globalMin.m_x = b.min.m_x;
        if (b.min.m_y < globalMin.m_y) globalMin.m_y = b.min.m_y;
        if (b.min.m_z < globalMin.m_z) globalMin.m_z = b.min.m_z;

        if (b.max.m_x > globalMax.m_x) globalMax.m_x = b.max.m_x;
        if (b.max.m_y > globalMax.m_y) globalMax.m_y = b.max.m_y;
        if (b.max.m_z > globalMax.m_z) globalMax.m_z = b.max.m_z;
    }

    // 3. 稍微扩大边界 (Epsilon Padding)
    // 避免点恰好落在边界上导致的索引越界或精度问题
    Vector margin = (globalMax - globalMin) * 0.001;
    if (margin.Mag() < 1e-6) margin = Vector(0.1, 0.1, 0.1);

    gridParams_.minBound = globalMin - margin;
    gridParams_.maxBound = globalMax + margin;

    // 4. 计算 Bin 尺寸
    Vector extent = gridParams_.maxBound - gridParams_.minBound;
    gridParams_.nx = std::max(1, resolutionX);
    gridParams_.ny = std::max(1, resolutionY);
    gridParams_.nz = std::max(1, resolutionZ);

    gridParams_.binSize.m_x = extent.m_x / gridParams_.nx;
    gridParams_.binSize.m_y = extent.m_y / gridParams_.ny;
    gridParams_.binSize.m_z = extent.m_z / gridParams_.nz;

    // 5. 将所有 Cell 注册到 Bins 中
    // 策略：Cell AABB 覆盖法 (Broad Phase)
    // 只要 Cell 的 AABB 覆盖了某个 Bin，就将 Cell ID 放入该 Bin
    for (size_t idx = 0; idx < cells_.size(); ++idx)
    {
        const Cell& cell = cells_[idx];
        const AABB& box = cell.boundingBox;

        // 计算 Cell AABB 覆盖的索引范围 [min_ijk, max_ijk]
        int i_min, j_min, k_min;
        int i_max, j_max, k_max;

        gridParams_.getIJK(box.min, i_min, j_min, k_min);
        gridParams_.getIJK(box.max, i_max, j_max, k_max);

        // 遍历覆盖的所有 Bin 并注册
        for (int k = k_min; k <= k_max; ++k) {
            for (int j = j_min; j <= j_max; ++j) {
                for (int i = i_min; i <= i_max; ++i) {
                    int binID = gridParams_.getBinIndex(i, j, k);
                    cellBins_[binID].push_back(static_cast<int>(idx));
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms = end - start;

    std::cout << "       [Done] Grid Params: Origin("
        << gridParams_.minBound.m_x << "," << gridParams_.minBound.m_y << "," << gridParams_.minBound.m_z << ")"
        << " BinSize(" << gridParams_.binSize.m_x << "," << gridParams_.binSize.m_y << "," << gridParams_.binSize.m_z << ")"
        << std::endl;
    std::cout << "       [Time] " << ms.count() << " ms. Total Bins with Data: " << cellBins_.size() << std::endl;
}

/**
 * @brief 在buildCellBins之后，基于Rasterization构建候选cell列表，其他两种方式则是利用八叉树和暴力遍历
 */
void Mesh::getCandidateCells_Rasterization(
    const Vector& p1, const Vector& p2, const Vector& p3,
    std::vector<int>& outCellIDs) const
{
    outCellIDs.clear();

    // 0. 安全检查：如果背景网格未构建，直接返回
    if (gridParams_.nx == 0 || gridParams_.ny == 0 || gridParams_.nz == 0) {
        // 生产级容错：如果未构建，理论上应回退或报错，这里留空由上层处理
        return;
    }

    // 1. 计算三角形的 AABB (用于确定光栅化遍历范围)
    // 这一步能快速缩小需要检测的 Bin 范围 (Broad Phase inside Grid)
    double minX = std::min({ p1.m_x, p2.m_x, p3.m_x });
    double minY = std::min({ p1.m_y, p2.m_y, p3.m_y });
    double minZ = std::min({ p1.m_z, p2.m_z, p3.m_z });
    double maxX = std::max({ p1.m_x, p2.m_x, p3.m_x });
    double maxY = std::max({ p1.m_y, p2.m_y, p3.m_y });
    double maxZ = std::max({ p1.m_z, p2.m_z, p3.m_z });

    // 2. 获取覆盖的 Bin 索引范围 [min_ijk, max_ijk]
    int i_min, j_min, k_min;
    int i_max, j_max, k_max;

    // 使用 AABB min/max 获取范围
    gridParams_.getIJK(Vector(minX, minY, minZ), i_min, j_min, k_min);
    gridParams_.getIJK(Vector(maxX, maxY, maxZ), i_max, j_max, k_max);

    // 3. 准备 Tribox 算法所需的 Box HalfSize (所有 Bin 尺寸相同)
    double binEps = 1e-5;
    Vector binHalfSize = gridParams_.binSize * 0.5 + Vector(binEps, binEps, binEps);

    // 使用 HashSet 临时去重 (一个 Cell 可能跨越多个 Bin)
    // 也可以用 vector + sort + unique，视 Cell 密度而定
    std::vector<int> rawCandidates;

    // 4. 遍历范围内的所有 Bin
    for (int k = k_min; k <= k_max; ++k) {
        for (int j = j_min; j <= j_max; ++j) {
            for (int i = i_min; i <= i_max; ++i) {

                // [Core] 精确光栅化判断: 三角形是否真的切过这个 Bin?
                // 计算当前 Bin 的几何中心
                // Center = MinBound + (Index + 0.5) * Size
                double cx = gridParams_.minBound.m_x + (i + 0.5) * gridParams_.binSize.m_x;
                double cy = gridParams_.minBound.m_y + (j + 0.5) * gridParams_.binSize.m_y;
                double cz = gridParams_.minBound.m_z + (k + 0.5) * gridParams_.binSize.m_z;
                Vector binCenter(cx, cy, cz);

                // 调用 Task 1 实现的 Tribox 算法
                if (EDFM_Geometry_3D::TriBoxOverlap(binCenter, binHalfSize, p1, p2, p3))
                {
                    // 如果相交，获取该 Bin 内注册的所有 Cell ID
                    int binID = gridParams_.getBinIndex(i, j, k);
                    auto it = cellBins_.find(binID);
                    if (it != cellBins_.end()) {
                        const auto& cellsInBin = it->second;
                        rawCandidates.insert(rawCandidates.end(), cellsInBin.begin(), cellsInBin.end());
                    }
                }
            }
        }
    }

    // 5. 去重与输出
    if (!rawCandidates.empty()) {
        std::sort(rawCandidates.begin(), rawCandidates.end());
        auto last = std::unique(rawCandidates.begin(), rawCandidates.end());
        outCellIDs.assign(rawCandidates.begin(), last);
    }
}

// 获取背景网格参数 (供调试或外部算法使用)
const UniformGrid3DParams& Mesh::getGridParams() const { return gridParams_; }

// =========================================================
// 6. 【2D-EDFM】2D基岩网格与1D裂缝拓扑关系
// =========================================================

/**
 * @brief   【适用2D EDFM】 重建 "Cell ID -> 全局裂缝段 ID 列表" 的映射
 * @param   储存全部裂缝的裂缝网络容器
 * @param   输出全局裂缝段全局编号类FracElemIndex
 * @details  基于FracIndex.h中定义的全局裂缝段全局编号类FracElemIndex，构建与基岩网格单元cell Id之间的拓扑关系
 */
void Mesh::rebuildFractureMap(const std::vector<Fracture>& fractures)
{
    // 1. 清空旧映射
    CellLocalIndexToFracElemSolverIndexMap_.clear();
    CellLocalIndexToFracElemLocalIdMap_.clear();
    int mappedCount = 0;

    // 2. 遍历所有裂缝
    for (const auto& frac : fractures)
    {
        for (const auto& elem : frac.elements)
        {
            // 2.1 基础有效性检查
            if (elem.cellID < 0) continue;

            // 2.2 [安全检查] 确保 Solver Index 已分配
            if (elem.solverIndex == -1)
            {
                std::cerr << "[Error] Fracture Element (ID=" << elem.id
                    << ") has invalid solverIndex (-1). "
                    << "Please call BuildGlobalSystemIndexing() BEFORE rebuildFractureMap()."
                    << std::endl;
                continue; // 跳过无效单元
            }

            auto it = cellId2index_.find(elem.cellID);
            if (it != cellId2index_.end())
            {
                int matrixLocalIndex = it->second;

                CellLocalIndexToFracElemSolverIndexMap_[matrixLocalIndex].push_back(elem.solverIndex);
                CellLocalIndexToFracElemLocalIdMap_[matrixLocalIndex].push_back(elem.id);
                mappedCount++;
            }
            else
            {
                std::cerr << "[Warning] Fracture Element refers to unknown Matrix Cell ID: "
                    << elem.cellID << std::endl;
            }
        }
    }
    std::cout << "[Mesh] Fracture map rebuilt. Mapped " << CellLocalIndexToFracElemSolverIndexMap_.size()
        << " active matrix cells (Total embedded elements: " << mappedCount << ")." << std::endl;
}

/**
 * @brief   【适用2D EDFM】 获取特定基岩网格内存在的裂缝段
 * @param   特定网格单元的全局ID
 * @details  访问经rebuildFractureMap函数构建与基岩网格单元cell Id之间的拓扑关系之后，并获取CellLocalIndexToFracElemSolverIndexMap_的value值
 */
std::vector<int> Mesh::getFracElemSolverIndexfromCellGlobalId(int cellGlobalID) const
{
    // 1. 先将输入的 Global ID 转换为 Local Index
    auto itIdx = cellId2index_.find(cellGlobalID);
    if (itIdx == cellId2index_.end())
    {
        return {}; // 找不到该单元
    }
    int localIndex = itIdx->second;

    // 2. 使用 Local Index 查询映射表
    auto it = CellLocalIndexToFracElemSolverIndexMap_.find(localIndex);
    if (it != CellLocalIndexToFracElemSolverIndexMap_.end())
    {
        return it->second;
    }
    return {};
}

/**
 * @brief   【适用2D EDFM】 获取特定基岩网格内存在的裂缝段
 * @param   特定网格单元的全局ID
 * @details 访问经rebuildFractureMap函数构建与基岩网格单元cell Id之间的拓扑关系之后，并获取CellLocalIndexToFracElemLocalIdMap_的value值
 */
std::vector<int> Mesh::getFracElemLocalIdfromCellGlobalId(int cellGlobalID) const
{
    // 1. 先将输入的 Global ID 转换为 Local Index
    auto itIdx = cellId2index_.find(cellGlobalID);
    if (itIdx == cellId2index_.end())
    {
        return {}; // 找不到该单元
    }
    int localIndex = itIdx->second;

    // 2. 使用 Local Index 查询映射表
    auto it = CellLocalIndexToFracElemLocalIdMap_.find(localIndex);
    if (it != CellLocalIndexToFracElemLocalIdMap_.end())
    {
        return it->second;
    }
    return {};

}

// =========================================================
// 7. I/O 与 调试接口 (Input/Output & Debug)
// =========================================================

/**
 * @brief 【适用2D/3D EDFM】在终端打印网格信息
 */
void Mesh::printMeshInfo()
{
    std::cout << "================ Mesh Statistics ================" << std::endl;
    std::cout << "总节点数 (Nodes) : " << nodes_.size() << endl;
    std::cout << "总单元数 (Cells) : " << cells_.size() << endl;
    std::cout << "总面数   (Faces) : " << faces_.size() << endl;


    int nInner = 0, nBoundary = 0;
    for (const auto& cell : cells_) {
        if (cell.location == Cell::LocationType::Inner) ++nInner;
        else if (cell.location == Cell::LocationType::Boundary) ++nBoundary;
    }
    std::cout << "内部单元数: " << nInner << endl;
    std::cout << "边界单元数: " << nBoundary << endl;

    // 仅当网格规模较小时打印详细信息，防止刷屏
    if (cells_.size() < 100)
    {
        std::cout << "\n---------------- Node Info ----------------" << endl;
        for (const auto& node : nodes_)
        {
            std::cout << "ID: " << node.id << " \t("
                << node.coord.m_x << ", " << node.coord.m_y << ", " << node.coord.m_z << ")"
                << endl;
        }

        std::cout << "\n---------------- Cell Info ----------------" << endl;
        for (const auto& cell : cells_)
        {
            if (cell.id < 0) continue;
            std::cout << "ID: " << cell.id
                << " \tType: " << (cell.location == Cell::LocationType::Inner ? "Inner" : "Bndry")
                << " \tVol/Area: " << cell.volume
                << " \tCenter: (" << cell.center.m_x << ", " << cell.center.m_y << ", " << cell.center.m_z << ")"
                << endl;
            std::cout << "    Faces: [ ";
            for (int fid : cell.CellFaceIDs) std::cout << fid << " ";
            std::cout << "]" << endl;
        }

        std::cout << "\n---------------- Face Info ----------------" << endl;
        for (const auto& face : faces_)
        {
            std::cout << "ID: " << face.id
                << " \tOwner: " << face.ownerCell
                << " \tNeighbor: " << face.neighborCell
                << " \tMeasure: " << face.length // 2D是长度，3D是面积
                << endl;
            std::cout << "    Normal: (" << face.normal.m_x << ", " << face.normal.m_y << ", " << face.normal.m_z << ")" << endl;
            std::cout << "    Nodes: ";
            for (int nid : face.FaceNodeIDs) std::cout << nid << " ";
            std::cout << endl;
        }
        std::cout << "=================================================" << std::endl;
    }
}

/**
 * @brief 【适用2D/3D EDFM】导出网格数据到文本文件
 * @details 包含 _nodes.txt, _faces.txt (变长格式), _cells.txt
 */
void Mesh::exportToTxt(const std::string& prefix) const
{
    // 1. Nodes (格式: id x y z)
    std::ofstream fn(prefix + "_nodes.txt");
    // fn << "id x y z\n"; // 为了Matlab读取方便，建议不输出表头，或者Matlab端跳过一行
    for (const auto& n : nodes_)
        fn << n.id << " "
        << n.coord.m_x << " "
        << n.coord.m_y << " "
        << n.coord.m_z << "\n";
    fn.close();

    // 2. Faces (格式: id num_nodes n1 n2 ... [last_node] owner neighbor)
    // 这种变长格式可以兼容 2D(2节点) 和 3D(3-4节点)
    std::ofstream ff(prefix + "_faces.txt");
    for (const auto& f : faces_) {
        ff << faceId2Index_.at(f.id) << " " << f.FaceNodeIDs.size(); // 输出节点数
        for (int nid : f.FaceNodeIDs) {
            ff << " " << nid;
        }
        ff << " " << f.ownerCell_index << " " << f.neighborCell_index << "\n";
    }
    ff.close();

    // 3. Cells (格式: id volume cx cy cz type)
    // type: 0=Inner, 1=Boundary

    std::ofstream fc(prefix + "_cells.txt");
    for (const auto& c : cells_) {
        fc << cellId2index_.at(c.id) << " "
            << c.volume << " "
            << c.center.m_x << " "
            << c.center.m_y << " "
            << c.center.m_z << " "
            << (c.location == Cell::LocationType::Boundary ? 1 : 0) << "\n";
    }
    fc.close();

    std::cout << "[Export] Mesh data exported to " << prefix << "_*.txt (Nodes, Faces, Cells)" << std::endl;
}

/**
* @brief 【适用2D/3D EDFM】获取特定网格单元cell的面积
* @param 特定网格单元编号，读取经Cell中computeCenterAndVolume_*D行为计算得到的cell.volume
*/
double Mesh::getCellArea(int cellID) const
{
    // 利用 map 加速查找 (前提是 BuildMesh 中已构建 cellId2index_)
    auto it = cellId2index_.find(cellID);
    if (it != cellId2index_.end()) {
        return cells_[it->second].volume;
    }
    return 0.0;
}

// 私有辅助函数 (内部逻辑解耦)
// =========================================================
// 【2D/3D通用】私有辅助函数 ：1. 构建网格Node拓扑关系
// =========================================================
void Mesh::_fetchNodesFromGmsh()
{
    std::vector<size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);

    nodes_.clear(); nodes_.reserve(nodeTags.size());
    nodesMap_.clear(); nodesMap_.reserve(nodeTags.size());

    for (size_t i = 0; i < nodeTags.size(); i++) {
        int id = static_cast<int>(nodeTags[i]);
        Vector coord(nodeCoords[3 * i], nodeCoords[3 * i + 1], nodeCoords[3 * i + 2]);
        nodes_.emplace_back(id, coord);
        nodesMap_[id] = Node(id, coord);
    }
}

// =========================================================
// 【2D】 私有辅助函数：2.1 构建2D网格cell拓扑关系
// =========================================================
void Mesh::_fetchCells2D()
{
    std::vector<int> elementTypes;
    std::vector<std::vector<size_t>> elementTags, nodeTagsPerElement;
    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElement, 2);

    cells_.clear();
    cellId2index_.clear();

    for (size_t i = 0; i < elementTypes.size(); i++) {
        int type = elementTypes[i];
        if (type != 2 && type != 3) continue; // Tri(2), Quad(3)

        size_t numElems = elementTags[i].size();
        int nodesPerElem = (type == 2) ? 3 : 4;

        for (size_t j = 0; j < numElems; j++)
        {
            int cellId = static_cast<int>(elementTags[i][j]);
            std::vector<int> ids;
            ids.reserve(nodesPerElem);
            for (int k = 0; k < nodesPerElem; k++) {
                ids.push_back(static_cast<int>(nodeTagsPerElement[i][nodesPerElem * j + k]));
            }

            Cell cell(cellId, ids);
            // 2D 专用计算
            cell.computeCenterAndVolume_2D(nodesMap_);
            cells_.push_back(cell);
        }
    }
    for (size_t i = 0; i < cells_.size(); ++i) cellId2index_[cells_[i].id] = static_cast<int>(i);
    gridCount_ = static_cast<int>(cells_.size());
}

// =========================================================
// 【3D】 私有辅助函数：2.2 构建3D网格cell拓扑关系
// =========================================================
void Mesh::_fetchCells3D()
{
    cells_.clear();
    cellId2index_.clear();

    // 显式指定 dim=3，只获取体单元，过滤掉所有 2D 边界单元
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, 3);

    for (const auto& entity : entities) {
        int dim = entity.first;
        int tag = entity.second;

        std::vector<int> elemTypes;
        std::vector<std::vector<size_t>> elemTags;
        std::vector<std::vector<size_t>> elemNodeTags;

        gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

        for (size_t i = 0; i < elemTypes.size(); ++i) {
            int type = elemTypes[i];
            int nodesPerElem = 0;
            // 根据 Gmsh 文档映射节点数
            if (type == 4) nodesPerElem = 4;      // Tet
            else if (type == 5) nodesPerElem = 8; // Hex
            else if (type == 6) nodesPerElem = 6; // Prism
            else if (type == 7) nodesPerElem = 5; // Pyramid
            else continue; // 跳过非体单元

            size_t numElems = elemTags[i].size();
            for (size_t j = 0; j < numElems; j++) {
                int cellId = static_cast<int>(elemTags[i][j]);
                std::vector<int> ids;
                ids.reserve(nodesPerElem);

                for (int k = 0; k < nodesPerElem; k++) {
                    // 转换 nodeTag 到 cells_ 内部存储的 ID (假设与 nodesMap key 一致)
                    ids.push_back(static_cast<int>(elemNodeTags[i][nodesPerElem * j + k]));
                }

                Cell cell(cellId, ids);
                cell.computeCenterAndVolume_3D(nodesMap_);
                cells_.push_back(cell);
            }
        }
    }

    gridCount_ = static_cast<int>(cells_.size());
    for (int i = 0; i < gridCount_; ++i) {
        cellId2index_[cells_[i].id] = i;
    }
}

// =========================================================
// 【2D】 私有辅助函数：3.1 构建2D网格face的拓扑关系 & 实现边界Tag映射
// =========================================================
void Mesh::_buildTopology2D()
{
    struct FaceAcc { std::vector<int> ordered; std::vector<int> cells; };
    std::unordered_map<std::vector<int>, FaceAcc, VectorHash> faceMap;

    // 1. 遍历 Cell 收集 Edge
    for (const auto& cell : cells_) {
        auto localFaces = cell.getLocalFaces_2D(nodesMap_);
        for (const auto& fNodes : localFaces) {
            std::vector<int> key = fNodes;
            std::sort(key.begin(), key.end());
            auto& acc = faceMap[key];
            if (acc.ordered.empty()) acc.ordered = fNodes;
            acc.cells.push_back(cell.id);
        }
    }

    faces_.clear();
    // 映射表：Key(Sorted Nodes) -> Face Index in faces_
    std::map<std::vector<int>, int> nodeSetToFaceIndex;

    int faceIdCounter = 1;
    for (auto& kv : faceMap) {
        const auto& orderedNodes = kv.second.ordered;
        std::vector<Vector> coords;
        for (int nid : orderedNodes) coords.push_back(nodesMap_.at(nid).coord);

        Face face(faceIdCounter, orderedNodes, coords);
        face.ownerCell = kv.second.cells[0];
        face.neighborCell = (kv.second.cells.size() > 1) ? kv.second.cells[1] : -1;

        // 法向校正
        if (face.neighborCell != -1) {
            Vector co = cells_[cellId2index_[face.ownerCell]].center;
            Vector cn = cells_[cellId2index_[face.neighborCell]].center;
            if ((face.normal * (cn - co)) < 0)
            {
                // 翻转
                std::vector<int> revNodes = orderedNodes;
                std::reverse(revNodes.begin(), revNodes.end());
                std::vector<Vector> revCoords = coords;
                std::reverse(revCoords.begin(), revCoords.end());
                face = Face(faceIdCounter, revNodes, revCoords);
                face.ownerCell = kv.second.cells[0];
                face.neighborCell = kv.second.cells[1];
            }
        }

        // 记录边界面索引用于 Tag 匹配
        if (face.neighborCell == -1) {
            // 【新增】边界面：确保指向外部 (Owner -> Out)
            // 逻辑：法向 必须与 (Face中心 - Owner中心) 同向
            Vector co = cells_[cellId2index_[face.ownerCell]].center;
            Vector faceCenter = face.midpoint; // 假设 Face 有 midpoint 或 center 属性
			if ((face.normal * (faceCenter - co)) < 0)
			{
				// 翻转
				std::vector<int> revNodes = orderedNodes;
				std::reverse(revNodes.begin(), revNodes.end());
				std::vector<Vector> revCoords = coords;
				std::reverse(revCoords.begin(), revCoords.end());
				face = Face(faceIdCounter, revNodes, revCoords);
				face.ownerCell = kv.second.cells[0];
				face.neighborCell = -1;
			}
            nodeSetToFaceIndex[kv.first] = static_cast<int>(faces_.size());
        }

        faces_.push_back(face);
        // 回填 Cell
        for (int cid : kv.second.cells) cells_[cellId2index_[cid]].CellFaceIDs.push_back(faceIdCounter);
        faceIdCounter++;
    }

    // 2. 物理组映射 (Line -> Edge Tag)
    gmsh::vectorpair physGroups;
    gmsh::model::getPhysicalGroups(physGroups, 1); // 1D Entities

    for (auto& p : physGroups) {
        int tag = p.second;
        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(1, tag, entities);

        for (int ent : entities) {
            std::vector<int> types;
            std::vector<std::vector<size_t>> tags, nodes;
            gmsh::model::mesh::getElements(types, tags, nodes, 1, ent);

            for (size_t k = 0; k < types.size(); k++) {
                if (types[k] != 1) continue; // Only 2-node lines
                size_t nElems = tags[k].size();
                for (size_t j = 0; j < nElems; j++) {
                    std::vector<int> key = { static_cast<int>(nodes[k][2 * j]), static_cast<int>(nodes[k][2 * j + 1]) };
                    std::sort(key.begin(), key.end());

                    if (nodeSetToFaceIndex.count(key)) {
                        faces_[nodeSetToFaceIndex[key]].physicalGroupId = tag;
                    }
                }
            }
        }
    }

}

// =========================================================
// 【3D】 私有辅助函数：3.2 构建3D网格face的拓扑关系 & 实现边界Tag映射
// =========================================================
void Mesh::_buildTopology3D()
{
    struct FaceAcc { std::vector<int> ordered; std::vector<int> cells; };
    std::unordered_map<std::vector<int>, FaceAcc, VectorHash> faceMap;

    // 1. 遍历 Cell 收集 Face
    for (const auto& cell : cells_) {
        auto localFaces = cell.getLocalFaces_3D();
        for (const auto& fNodes : localFaces) {
            std::vector<int> key = fNodes;
            std::sort(key.begin(), key.end());
            auto& acc = faceMap[key];
            if (acc.ordered.empty()) acc.ordered = fNodes;
            acc.cells.push_back(cell.id);
        }
    }

    faces_.clear();
    std::map<std::vector<int>, int> nodeSetToFaceIndex;

    int faceIdCounter = 1;
    for (auto& kv : faceMap) {
        const auto& orderedNodes = kv.second.ordered;
        std::vector<Vector> coords;
        for (int nid : orderedNodes) coords.push_back(nodesMap_.at(nid).coord);

        Face face(faceIdCounter, orderedNodes, coords);
        face.ownerCell = kv.second.cells[0];
        face.neighborCell = (kv.second.cells.size() > 1) ? kv.second.cells[1] : -1;

        // 法向校正
        if (face.neighborCell != -1) {
            Vector co = cells_[cellId2index_[face.ownerCell]].center;
            Vector cn = cells_[cellId2index_[face.neighborCell]].center;
            if ((face.normal * (cn - co)) < 0) {
                // 翻转
                std::vector<int> revNodes = orderedNodes;
                std::reverse(revNodes.begin(), revNodes.end());
                std::vector<Vector> revCoords = coords;
                std::reverse(revCoords.begin(), revCoords.end());
                face = Face(faceIdCounter, revNodes, revCoords);
                face.ownerCell = kv.second.cells[0];
                face.neighborCell = kv.second.cells[1];
            }
        }

        if (face.neighborCell == -1) {
            // 边界面：确保指向外部 (Owner -> Out)
            Vector co = cells_[cellId2index_[face.ownerCell]].center;
            Vector faceCenter = face.midpoint; // 假设 Face 有 midpoint 或 center 属性
            if ((face.normal * (faceCenter - co)) < 0)
            {
                // 翻转节点顺序
                std::vector<int> revNodes = orderedNodes;
                std::reverse(revNodes.begin(), revNodes.end());
                std::vector<Vector> revCoords = coords;
                std::reverse(revCoords.begin(), revCoords.end());

                // 重新构建 Face (更新法向)
                face = Face(faceIdCounter, revNodes, revCoords);
                face.ownerCell = kv.second.cells[0];
                face.neighborCell = -1; // 保持为边界
            }
            nodeSetToFaceIndex[kv.first] = static_cast<int>(faces_.size());
        }

        faces_.push_back(face);
        for (int cid : kv.second.cells) cells_[cellId2index_[cid]].CellFaceIDs.push_back(faceIdCounter);
        faceIdCounter++;
    }

    // 2. 物理组映射 (Surface -> Polygon Face Tag)
    gmsh::vectorpair physGroups;
    gmsh::model::getPhysicalGroups(physGroups, 2); // 2D Entities (Surfaces)

    for (auto& p : physGroups) {
        int tag = p.second;
        std::vector<int> entities;
        gmsh::model::getEntitiesForPhysicalGroup(2, tag, entities);

        for (int ent : entities) {
            std::vector<int> types;
            std::vector<std::vector<size_t>> tags, nodes;
            gmsh::model::mesh::getElements(types, tags, nodes, 2, ent);

            for (size_t k = 0; k < types.size(); k++) {
                int npe = 0;
                if (types[k] == 2) npe = 3;      // Triangle
                else if (types[k] == 3) npe = 4; // Quad
                else continue;

                size_t nElems = tags[k].size();
                for (size_t j = 0; j < nElems; j++) {
                    std::vector<int> key;
                    key.reserve(npe);
                    for (int m = 0; m < npe; m++) {
                        key.push_back(static_cast<int>(nodes[k][npe * j + m]));
                    }
                    std::sort(key.begin(), key.end());

                    if (nodeSetToFaceIndex.count(key)) {
                        faces_[nodeSetToFaceIndex[key]].physicalGroupId = tag;
                    }
                }
            }
        }
    }
}

// =========================================================
// 【2D/3D通用】私有辅助函数:利用 cellId2index_填充 Face 的 owner/neighbor 索引缓存
// =========================================================
void Mesh::_calcuteOwnerandNeighborCell_indexofFace()
{
    int errCount = 0;

    for (auto& face : faces_)
    {
        // --- Owner ---
        auto itO = cellId2index_.find(face.ownerCell);
        if (itO != cellId2index_.end()) {
            face.ownerCell_index = itO->second;
        }
        else {
            // 严重错误：Face 的 Owner ID 不在 Cell 列表中
            face.ownerCell_index = -1;
            errCount++;
        }
        // --- Neighbor ---
        if (face.neighborCell != -1) {
            auto itN = cellId2index_.find(face.neighborCell);
            if (itN != cellId2index_.end()) {
                face.neighborCell_index = itN->second;
            }
            else {
                // 可能是 ID 错误，也可能是边界处理问题
                face.neighborCell_index = -1;
                // 如果确定不是边界(-1)，却找不到 Index，则记录错误
                errCount++;
            }
        }
        else {
            face.neighborCell_index = -1; // 明确标记为边界
        }
    }
    if (errCount > 0) {
        std::cerr << "[Warning] Found " << errCount << " faces with invalid Cell IDs during index population." << std::endl;
    }
}













