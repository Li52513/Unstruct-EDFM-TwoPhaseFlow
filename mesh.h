#pragma once

// =========================================================
// 标准库与第三方库
// =========================================================
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>

// =========================================================
// 项目头文件
// =========================================================
#include "Node.h"
#include "Face.h"
#include "Cell.h"
#include "MatrixEdge.h"
#include "VectorHash.h"

// ==========================================
// 前向声明 (避免循环引用 Fracture.h)
// ==========================================
class Fracture;
struct FracElemIndex;

// =========================================================
// 【3D-EDFM】 3D 均匀背景网格参数结构体
// =========================================================
struct UniformGrid3DParams
{
    Vector minBound;        // 背景网格全局最小点 (Grid Origin)
    Vector maxBound;        // 背景网格全局最大点
    Vector binSize;         // 每个 Bin 的物理尺寸 (dx, dy, dz)
    int nx = 0;             // X 方向 Bin 数量
    int ny = 0;             // Y 方向 Bin 数量
    int nz = 0;             // Z 方向 Bin 数量

    // 辅助：根据物理坐标计算 Bin 的三维索引 (i, j, k)
    // 包含边界安全检查 (Clamp)
    void getIJK(const Vector& p, int& i, int& j, int& k) const
    {
        // 相对坐标
        double rx = p.m_x - minBound.m_x;
        double ry = p.m_y - minBound.m_y;
        double rz = p.m_z - minBound.m_z;

        // 计算索引 (向下取整)
        i = static_cast<int>(std::floor(rx / binSize.m_x));
        j = static_cast<int>(std::floor(ry / binSize.m_y));
        k = static_cast<int>(std::floor(rz / binSize.m_z));

        // 边界钳制 (Clamp to [0, n-1])
        if (i < 0) i = 0; else if (i >= nx) i = nx - 1;
        if (j < 0) j = 0; else if (j >= ny) j = ny - 1;
        if (k < 0) k = 0; else if (k >= nz) k = nz - 1;
    }

    // 辅助：计算 Flattened 1D 索引
    int getBinIndex(int i, int j, int k) const {
        return i + j * nx + k * nx * ny;
    }

    // 重置
    void reset() { nx = ny = nz = 0; }
};

/**
 * @brief 网格管理类 (Mesh)
 * @details 负责基岩网格的生成(Gmsh)、拓扑构建、几何计算及数据导出。
 * 支持 2D (Tri/Quad) 和 3D (Prism/Hex/Tet/Pyramid) 混合网格。
 */
class Mesh 
{
public:
    // =========================================================
    // 构造与析构
    // =========================================================
    Mesh();
    ~Mesh() = default;

    // =========================================================
    // 1. 【基岩网格】核心数据访问接口 (Accessors)
    // =========================================================

    // 【2D/3D通用】 实体访问
    const std::vector<Node>& getNodes() const;
    const std::vector<Face>& getFaces() const;
    const std::vector<Cell>& getCells() const;
    std::vector<Cell>& getCells();              

    /**
    * @brief 【适用3D EDFM】获取指定单元的所有棱
    * @details 返回一个储存在MatrixEdge.h中声明的MatrixEdge类容器，适用// 4=Tet, 5=Hex, 6=Prism, 7=Pyramid
    */
    std::vector<MatrixEdge> getCellEdges(int cellIndex) const;

    /**
     * @brief 【适用3D EDFM】提取全局唯一的基岩棱 (Matrix Edges)
     * @return 全局去重后的棱列表，用于 EDFM Type-3 (裂缝面-基岩棱) 求交检测
     */
    std::vector<MatrixEdge> extractUniqueMatrixEdges() const;


	// 【2D/3D通用】 获取私有成员节点、面、单元索引映射表的外部访问
    const std::unordered_map<int, Node>& getNodesMap() const;      // 经_fetchNodesFromGmsh()构建的 Node ID 到 Node 对象的映射，key值与 Gmsh nodeTag 一致，1-based
    const std::unordered_map<int, int>& getCellId2Index() const;   // 经_fetchCells2D()/3D()构建的 Cell Global ID（1based) 到 Local Vector Index(0-based) 的映射，key值与 Gmsh elementTag 一致，1-based
    const std::unordered_map<int, int>& getFaceId2Index()const;    // 经 buildFaceIndexMap构建的 Face Global ID（1based) 到 Local Vector Index(0-based) 的映射
    
    // 【2D/3D通用】 元数据访问 
    int getGridCount() const; // 获取网格总数
    std::string getMatrixElementType() const;


    // =========================================================
    // 2. 【基岩网格】 ID 转换安全接口 (Safe ID Lookups)
    // =========================================================

    /**
     * @brief 【2D/3D通用】通过 Cell Global ID (Gmsh Tag) 获取 Local Vector Index
     * @return 成功返回 index [0, N-1], 失败返回 -1
     */
    int getCellIndex(int cellGlobalID) const;

    /**
     * @brief 【2D/3D通用】 通过 Face Global ID 获取 Local Vector Index
     * @return 成功返回 index [0, N-1], 失败返回 -1
     */
    int getFaceIndex(int faceGlobalID) const;

    // =========================================================
    // 3. 【基岩网格】 网格生成 (Generation)
    // =========================================================
    /**
     * @brief: 【适用2D EDFM】生成 2D 非结构化基岩网格
	 * 流程: Gmsh生成 -> 读取节点 -> 构建2D单元 -> 构建边(Face)拓扑 -> 映射物理边界 -> 构建网格face的Bin索引 -> 区分内部/边界单元 -> 填充face的owner neighbor 单元索引
     */
    void GenerateMesh2D(double lengthX, double lengthY, int nx, int ny, bool useQuadBase, std::string jobName = "EDFM_2D");

    /**
     * @brief: 【适用3D EDFM】生成 3D 非结构化基岩网格
     * 流程: Gmsh生成 -> 映射物理边界-> 读取节点 -> 构建3D单元 -> 构建边(Face)拓扑 ->  区分内部/边界单元-> 构建网格face的Bin索引 -> 填充face的owner neighbor 单元索引 ->提取基岩棱 (Matrix Edges)-> 构建空间索引 (Octree)
     */
    void  GenerateMesh3D(double lengthX, double lengthY, double lengthZ, int nx, int ny, int nz, bool usePrism, bool useQuadBase, std::string jobName = "EDFM_3D");



    // =========================================================
	// 4. 【基岩网格】 拓扑关系与网格面几何信息计算 (Topology)
    // =========================================================
    /**
     * @brief 【2D/3D通用】 构建 Face Global ID -> Local Index 的映射表
     * @details 必须GenerateMesh*D之后即在网格生成完毕后调用
     */
    void buildFaceGlobalId2LocalIndexMap();     
    
    /**
     * @brief 【2D/3D通用】对基岩单元进行分类 (Inner / Boundary)
     * @details 在GenerateMesh*D函数中调用，实现区分 Inner (内部) 和 Boundary (边界) 单元
     */
    void ClassifySolidMatrixCells();
    
    /**
     * @brief 【2D/3D通用】 计算基岩网格面的几何信息
     * @param method 法向向量修正方法 (最小修正/正交修正/超松弛)
     * @details 计算面法向、面积矢量、正交(E)与非正交(T)分量
     */
    void ComputeMatrixMeshFaceGeometricInfor(NormalVectorCorrectionMethod method);
    
    // =========================================================
	// 5.1 【2D-EDFM】 2D基岩网格与1D裂缝快速求交算法
    // =========================================================
    
    /**
     * @brief 【适用2D EDFM】 为基岩网格整体构建第一次AABB索引
     * @details 为利用背景网格索引快速获取候选面建立基础，即Fracture::DetectFracturetoMeshFaceIntersections中的方法3
     */
	void buildFaceBins_2D();  //通过调整TARGET_FACES_PER_BIN大小，可以影响背景网格的细化程度，从而影响后续裂缝面求交的效率

    /**
     * @brief  【适用2D EDFM】 通过2D基岩网格面（1D）face与待几何离散的裂缝的AABB包围盒，判断是否需要进行下一步（二次AABB求交或直接精准计算是否相交）
     * @param   待几何离散的裂缝的AABB包围盒
     * @details 由于在2D情况，裂缝网格是通过1D宏观裂缝与2D基岩网格face的交点、裂缝的起终点以及F-F交点进行网格划分的
     * 这里裂缝的AABB包围盒，经buildFaceBins_2D构建基岩网格面的AABB后这里判断两AABB是否重合，形成需要第二次粗筛的基岩网格面候选面 ID
     */
    std::vector<int> getCandidateFacesFromBins(const AABB& box) const;

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
    std::vector<int> getCandidateFacesFromBins_DDA(const Vector& p1, const Vector& p2) const;


    // =========================================================
    // 5.2 【3D-EDFM】 3D基岩网格与2D裂缝快速求交算法
    // =========================================================
     
	// @brief【适用3D EDFM】 触发所有 Cell 的 AABB 计算，计算所有Cell的AABB包围盒,为后续构建背景网格索引做准备
    void computeAllCellAABBs();

    /**
     * @brief 【适用3D EDFM】构建 3D 均匀背景网格索引 (Binning)
     * @details 将所有基岩单元根据其 AABB 注册到覆盖的 Bin 中
     * @param resolutionX X 方向期望划分的份数 
     * @param resolutionY Y 方向期望划分的份数
     * @param resolutionZ Z 方向期望划分的份数
     */
    void buildCellBins(int resolutionX, int resolutionY, int resolutionZ);
   
    /**
     * @brief 在buildCellBins之后，基于Rasterization构建候选cell列表，其他两种方式则是利用八叉树和暴力遍历
     */
    void getCandidateCells_Rasterization(
        const Vector& p1, const Vector& p2, const Vector& p3,
        std::vector<int>& outCellIDs) const;

    // 获取背景网格参数 (供调试或外部算法使用)
    const UniformGrid3DParams& getGridParams() const;

    // =========================================================
    // 6. 【2D-EDFM】2D基岩网格与1D裂缝拓扑关系
    // =========================================================
    /**
     * @brief   【适用2D EDFM】 重建 "Cell ID -> 全局裂缝段 ID 列表" 的映射
     * @param   储存全部裂缝的裂缝网络容器
     * @param   输出全局裂缝段全局编号类FracElemIndex
     * @details  基于FracIndex.h中定义的全局裂缝段全局编号类FracElemIndex，构建与基岩网格单元cell Id之间的拓扑关系
     */
    void rebuildFractureMap(const std::vector<Fracture>& fractures);

    /**
     * @brief   【适用2D EDFM】 获取特定基岩网格内存在的裂缝段
     * @param   特定网格单元的全局ID
     * @details  访问经rebuildFractureMap函数构建与基岩网格单元cell Id之间的拓扑关系之后，并获取CellLocalIndexToFracElemSolverIndexMap_的value值
     */
    std::vector<int> getFracElemSolverIndexfromCellGlobalId(int cellGlobalID) const;

    /**
     * @brief   【适用2D EDFM】 获取特定基岩网格内存在的裂缝段
     * @param   特定网格单元的全局ID
     * @details 访问经rebuildFractureMap函数构建与基岩网格单元cell Id之间的拓扑关系之后，并获取CellLocalIndexToFracElemLocalIdMap_的value值
     */
    std::vector<int> getFracElemLocalIdfromCellGlobalId(int cellGlobalID) const;

    // =========================================================
    // 7. I/O 与 调试接口 (Input/Output & Debug)
    // =========================================================
   
    /**
     * @brief 【适用2D/3D EDFM】在终端打印网格信息
     */
     void printMeshInfo();
 
    /**
     * @brief 【适用2D/3D EDFM】导出网格数据到文本文件
     * @details 包含 _nodes.txt, _faces.txt (变长格式), _cells.txt
     */
    void exportToTxt(const string& prefix) const;

    /**
    * @brief 【适用2D/3D EDFM】获取特定网格单元cell的面积
    * @param 特定网格单元编号，读取经Cell中computeCenterAndVolume_*D行为计算得到的cell.volume
    */
    // 计算网格单元 cell 的面积 (输入参数：int Cell编号id)
    double getCellArea(int cellID) const;
    



private:

    // 私有成员变量 (Private Members)

    // =========================================================
    // 【基岩网格】
    // =========================================================

    // 【2D/3D通用】核心拓扑数据 --- 
    vector<Node> nodes_;
	vector<Face> faces_;
	vector<Cell> cells_;

    //【2D/3D通用】
    unordered_map<int, Node> nodesMap_;
    unordered_map<int, int> cellId2index_;
    unordered_map<int, int> faceId2Index_;

    //【2D/3D通用】
    int gridCount_ = 0; // 网格总数
    
    //【3D】记录基岩网格的主类型
    std::string matrixElementType_ = "Unknown";

    // =========================================================
    // 【EDFM拓扑关系】
    // =========================================================

    //【2D-EDFM】裂缝段与所在cell之间的索引
    /// Key: CellID
    /// Value: 列表包含该 Cell 内所有裂缝段的 [全局 ID]
    std::unordered_map<int, std::vector<int>> CellLocalIndexToFracElemSolverIndexMap_;
    std::unordered_map<int, std::vector<int>> CellLocalIndexToFracElemLocalIdMap_;

    // 【2D-EDFM】空间索引数据 (2D-EDFM，裂缝快速划分方法 Method 3 核心变量) ---
    double globalMinX_ = 0.0;
    double globalMinY_ = 0.0;
	int binCountX_=0; // X 方向网格划分数
	int binCountY_ = 0; // Y 方向网格划分数
    double binSizeX_ = 0;   // X 方向实际 bin 宽度
    double binSizeY_ = 0;   // Y 方向实际 bin 高度

    // 【2D-EDFM】
	static constexpr int TARGET_FACES_PER_BIN = 10; // 目标每个 bin 包含的面数
    
    // 【2D-EDFM】使用 unordered_map 存储稀疏网格 (Key: binID, Value: FaceIDs)
	unordered_map<int, vector<int>> faceBins_; // 面空间划分的 bin

    // 【3D-EDFM】  光栅化数据成员
    UniformGrid3DParams gridParams_;

    // 【3D-EDFM】 空间哈希表
    // Key: Flattened Bin ID (i + j*nx + k*nx*ny)
    // Value: 该 Bin 内包含的 Matrix Cell 索引列表 (Local Index)
    // 使用 unordered_map 以支持稀疏存储 (大量空 Bin 不占内存)
    std::unordered_map<int, std::vector<int>> cellBins_;

    
    // 私有辅助函数 (内部逻辑解耦)
 
    // =========================================================
    // 【基岩网格】
    // =========================================================
   
    // 1. 【2D/3D-EDFM通用】基础数据读取 (通用)
    void _fetchNodesFromGmsh();

    // 2. 【2D/3D-EDFM通用】单元读取 (维度分离)
    void _fetchCells2D();
    void _fetchCells3D();

    // 3. 拓扑构建 (维度分离)
    //    2D: 识别 Edge, 处理 Physical Curves
    void _buildTopology2D();
    //    3D: 识别 Polygon Face, 处理 Physical Surfaces
    void _buildTopology3D();

    /**
     * @brief  【2D/3D-EDFM通用】 填充 Face 的 owner/neighbor 索引缓存
     * @details 遍历所有 Face，利用 cellId2index_ 将 Global ID 映射为 Local Index
     * 必须在 buildFaceGlobalId2LocalIndexMap() 或网格生成最后一步调用。
     */
    void _calcuteOwnerandNeighborCell_indexofFace();

    // =========================================================
    // 【EDFM拓扑关系】
    // =========================================================

    // 【2D-EDFM】计算 binID
    inline int computeBinID(int ix, int iy) const 
    {
        return ix * binCountY_ + iy;
    }

};