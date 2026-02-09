#pragma once

#include <vector>
#include <string>
#include <sstream>

#include "Mesh.h"
#include "BoundaryFaceClassify_ByTag.h"
#include "2D_FractureNetwork.h"
#include "2D_Fracture.h"




// ====================================================================================================================================================
// 3D EDFM (3D基岩-2D裂缝，3D基岩-2D裂缝独立网格划分，利用三种类型交点对2D裂缝网格进一步在基岩网格内构建交互多边形，形成拓扑关系，建立质量交换系数
// ====================================================================================================================================================


// =========================================================
// 结构：交互对 (Interaction Pair) - Level 3 Data Structure
// =========================================================
/**
 * @brief 描述一个有效的 基岩-裂缝 交互对
 * @details 包含了几何重构后的多边形信息，用于后续 CI 计算。
 * 解决了 "One Matrix - Multi Fracs" 和 "One Frac - Multi Matrices" 的多对多映射问题。
 */
struct InteractionPair
{
    // =========================================================
    // 1. 身份标识 (IDs for Debug/IO)
    // =========================================================
    int matrixCellGlobalID;             ///< 基岩单元全局 ID
    int fracElementGlobalID;            ///< 裂缝微元全局 ID (需确保 FractureElement_2D 有全局唯一ID)
    int fracMacroID;                    ///< 所属宏观裂缝 ID (用于分组 TI 计算)

    // =========================================================
    // 2. 求解器索引 (Indices for Matrix Assembly)
    // =========================================================
    int matrixSolverIndex;              ///< 基岩单元在全局矩阵 A 中的行号 (0-based)
    int fracCellSolverIndex;                ///< 裂缝单元在全局矩阵 A 中的行号 (0-based)


    std::vector<Vector> polygonPoints;  ///< 有序的多边形顶点 (CCW 逆时针排列)
    Vector polygonNormal;               ///< 多边形法向 (通常继承自裂缝面)
    Vector polygonCenter;               ///< 多边形几何中心 (用于 NNC 距离计算)
    double intersectionArea;            ///< 交互面积 (用于 NNC 传导率计算)

    double distMatrixToFracPlane;       ///< 基岩中心到裂缝平面的垂直距离 (d_NNC)

    // 构造函数
    InteractionPair(int mID, int fElemID, int fMacroID)
        : matrixCellGlobalID(mID), fracElementGlobalID(fElemID), fracMacroID(fMacroID),
        matrixSolverIndex(-1), fracCellSolverIndex(-1),
        intersectionArea(0.0), polygonNormal(0, 0, 0), polygonCenter(0, 0, 0), distMatrixToFracPlane(0.0)
    {
    }
};


/**
 * @class MeshManager_3D
 * @brief 统一管理基岩 & 裂缝网格及几何CI计算
 */

class MeshManager_3D
{
public:
    // =========================================================
    // 构造与初始化
    // =========================================================

     /**
     * @brief 构造函数
     * @param lx 区域长度 X
     * @param ly 区域长度 Y
     * @param lz 区域长度 Z
     * @param nx X 方向划分份数
     * @param ny Y 方向划分份数
     * @param nz Z 方向划分份数
     * @param usePrism 是否使用棱柱单元（true）或四面体（false）
     * @param useQuadBase 棱柱底面保留四边形（true）或分割为三角形（false）
     */
    MeshManager_3D(double lx, double ly, double lz, int nx, int ny, int nz, bool usePrism , bool useQuadBase); //默认有参构造函数


    // =========================================================
    // 定义求交策略枚举
    // =========================================================
    enum class IntersectionStrategy
    {
        BruteForce,             ///< 暴力遍历 (基准测试用)
        Octree_Optimized,       ///< 八叉树 AABB 预筛 (传统方法)
        Octree_With_14DOP,      ///< 八叉树 + 14-DOP 精筛
        Rasterization_14DOP     ///< [Task 3 新增] 3D 光栅化 + 14-DOP (极速方法)
    };

    // =========================================================
    // 基岩网格构建 & 预处理 (Grid Construction)
    // =========================================================

    //【基岩】
     /**
     * @brief 网格生成-利用gmsh构建3D基岩非结构化网格并计算3D基岩网格面的非正交修正，利用 Tag 分类并识别边界,边界会返回bcGroups_
     * @param NormalVectorCorrectionMethod：网格面非正交修正计算方法
     * 最小修正法 MinimumCorrection；
     *OrthogonalCorrection 正交修正法；
     *OverRelaxed 超松弛修正法
    */

    void BuildSolidMatrixGrid_3D(NormalVectorCorrectionMethod corr = NormalVectorCorrectionMethod::OrthogonalCorrection, std::string matrix_mesh_Name = "3D_Mesh_test");

    // =========================================================
    // [New] 全局索引管理接口
    // =========================================================
    /**
     * @brief 初始化全局 Solver Index 体系
     * @details
     * 1. 确认基岩网格数量 (Matrix DOF)。
     * 2. 调用 FractureNetwork_2D::distributeSolverIndices 分配裂缝 DOF。
     * @note 必须在 addFracture 和 meshAllFractures 之后，且在 SolveIntersection 之前调用。
     */
    void setupGlobalIndices();

    /**
     * @brief 交互数据清洗与去重 (Data Cleanup)
     * @details 执行两步清洗：
     * 1. 几何过滤：剔除“幽灵交互”，即距离超过基岩网格半对角线的错误关联。
     * 2. 智能去重：剔除 (MatrixID, FracID, Area) 完全一致的冗余计算。
     * * @note 必须在 SolveIntersection 之后，buildTopologyMaps 之前调用。
     */
    void removeDuplicateInteractions();

    /**
     * @brief 解决共面/边界导致的双重计算问题 (Step 3 of Cleanup)
     * @details 针对共面裂缝，如果总交互面积超过裂缝单元面积，
     * 通过几何中心竞争机制剔除冗余项，并进行最终归一化。
     */
    void resolveCoplanarInteractions();

    // =========================================================
    // [New] 拓扑映射构建接口
    // =========================================================
    /**
     * @brief 构建基岩与裂缝的双向交互拓扑映射
     * @details 建立 MatrixID -> Pairs 和 FracID -> Pairs 的快速索引。
     * 建议在 SolveIntersection 完成后调用。
     */
    void buildTopologyMaps();

    // =========================================================
    // [New] 拓扑查询接口 (供 Solver 或 Benchmark 使用)
    // =========================================================

    /**
     * @brief 获取指定基岩单元(SolverIndex)关联的所有交互对
     * @return 交互对指针列表 (若无交互返回空列表)
     */
    const std::vector<const InteractionPair*>& getInteractionsOfMatrix(int matrixSolverIdx) const;

    /**
     * @brief 获取指定裂缝单元(SolverIndex)关联的所有交互对
     * @return 交互对指针列表
     */
    const std::vector<const InteractionPair*>& getInteractionsOfFracture(int fracSolverIdx) const;

    const Fracture_2D* findFractureByID(const FractureNetwork_2D& net, int fracID);


    //【基岩】
    /**
     * @brief 接收BuildSolidMatrixGrid_*函数中返回的边界面分组bcGroups_提供读取接口
    */
    const BoundaryFaceClassify_byTag::FaceGroups& boundaryFaces_byTag() const { return bcGroups_byTag_; }

    //【基岩】
    /**
      * @brief 接收BuildSolidMatrixGrid_*函数中返回的边界面个数boundaryCount_提供读取接口
    */
    size_t countBoundaryFaces() const { return boundaryCount_; }

    //【裂缝】
    /**
     * @brief  添加2D裂缝至裂缝网络
     * @param  单个2D裂缝平面
    */
    void addFracturetoFractureNetwork(const Fracture_2D& frac);

    //【裂缝】
    /**
     * @brief 对裂缝网络内部所有2D宏观裂缝进行网格独立划分并重新构建全局裂缝网格段索引
     * @param u方向网格数
     * @param v方向网格数
     * @param method  非正交修正方法
    */
    void meshAllFracturesinNetwork(int nU, int nV, NormalVectorCorrectionMethod method = NormalVectorCorrectionMethod::OrthogonalCorrection);
     
    //【裂缝】
    /**
     * @brief 检测宏观裂缝网络内的裂缝-裂缝相交 (Element-wise)
     * @details 支持非共面裂缝，生成微观线段集合
     * @param strategy 策略开关 Octree_Optimized 八叉树优化  BruteForce 暴力遍历
     */
    void DetectFractureFractureIntersectionsInNetwork(FFIntersectionStrategy strategy = FFIntersectionStrategy::Octree_Optimized);

    /**
     * @brief 重建裂缝全局索引 (Wrapper)
     * @details 调用内部 FractureNetwork_2D 的 rebuildGlobalIndex，
     * 用于在划分网格后更新微元与宏观裂缝的映射关系 (offset 数组)。
     */
    void rebuildFractureGlobalIndex()
    {
        frNet_2D_.rebuildGlobalIndex();
    }
    // =========================================================
    //  3D-EDFM 核心求交接口 (F-M Intersection)-实施 三角化策略 可以有效解决扭转
    // =========================================================
    /**
     * @brief 执行 3D 基岩与 2D 裂缝的几何求交并重构交互多边形，仅适用于非扭曲裂缝网格，如果处理非扭曲裂缝，该函数即可
     * @details
     * 1. 策略筛选 (Octree/BruteForce) 确定候选对。
     * 2. 调用 EDFM_Geometry_3D 执行 Type 1/2/3 检测。
     * 3. 收集交点，去重，排序重构闭合多边形。
     * 4. 计算面积，剔除微小交互，存入 interactionPairs_。
     */
    void SolveIntersection3D_improved(IntersectionStrategy strategy = IntersectionStrategy::Octree_Optimized);

    //(Per-Triangle Exact Consistency)实施了 "Per-Sub-Triangle" 策略，将每一个裂缝网格再次划分为两个三角形，可以解决扭曲网格交互单元面积和裂缝网格面积不等的问题
    void SolveIntersection3D_improved_twist_accleration(IntersectionStrategy strategy = IntersectionStrategy::Octree_Optimized);

    // 获取交互对数据的接口 (供后续 NNC 计算调用)
    const std::vector<InteractionPair>& getInteractionPairs() const { return interactionPairs_; }

    // =========================================================
    // 输出 & 调试 (Export & Debug)
    // =========================================================
    //【基岩】调用Mesh类已有的exportToTxt函数，导出基岩网格的几何信息至txt，用于matlab脚本进行可视化后处理
    void exportMeshInfortoTxt(const std::string& prefix = " EDFM_test_") const;

    //【裂缝】调用 Fracture_2D 类已有的 exportToTxt 成员函数导出裂缝网络中所有宏观裂缝的几何信息另外导出交线信息至txt 用于matlab脚本进行可视化后处理
    void exportFracturesNetworkInfortoTxt(const std::string& prefix = " EDFM_test_") const;
    void exportFracturesNetworkInfortoTxt_improved(const std::string& prefix = " EDFM_test_") const;
    
    //【裂缝】调用FractureNetwork_2D类已有的inspectFractureEdges函数导出特定宏观裂缝裂缝微观裂缝的网格面非正交信息至.csv
    void inspectFractureEdges_non_orthogonalInfor(int fracID, const std::string& prefix = " EDFM_test_") const;

    //【裂缝】调用FractureNetwork_2D类已有的inspectIntersections函数检查并导出 F-F 交线及微观网格对应关系至 CSV
    void inspectIntersections_FracturetoFracture(const std::string& prefix = " EDFM_test_") const;

    //【EDFM】导出交互对信息用于验证
    void exportInteractionPairsToCSV(const std::string& filename) const;

    //【EDFM】导出基岩-裂缝交互多边形，供 MATLAB 可视化
    // 格式: NumVertices x1 y1 z1 x2 y2 z2 ...
    void exportInteractionPolygonsToTxt(const std::string& filename) const;
    void exportInteractionPolygonsToTxt_improved(const std::string& filename) const;

    //【EDFM】导出搜索空间验证数据
    void exportSearchSpaceToTxt(const std::string& filename, int targetFracID, IntersectionStrategy strategy);


   // =========================================================
   // 访问底层对象 (Accessors)
   // =========================================================
    Mesh& mesh() { return mesh_; }
    FractureNetwork_2D& fracture_network() { return frNet_2D_; }

    //  Const 重载版本 (供 TransmissibilitySolver 使用)
    const Mesh& mesh() const { return mesh_; }
    const FractureNetwork_2D& fracture_network() const { return frNet_2D_; }


private:
    // =========================================================
    // 私有成员变量
    // =========================================================
    Mesh mesh_; // 网格成员对象
    FractureNetwork_2D frNet_2D_;               // 裂缝网络成员对象

    //【基岩】
    double lx_, ly_, lz_;                       // 区域长度成员对象
    int nx_, ny_, nz_;                          // 网格数量成员对象
    bool usePrism_, useQuadBase_;               // 网格单元类型成员对象
    BoundaryFaceClassify_byTag::FaceGroups bcGroups_byTag_; //边界面组
    size_t boundaryCount_ = 0;                  //边界面数量
    
    //【3D-EDFM专用】
    // =========================================================
    // 存储所有有效的交互对
    // =========================================================
    std::vector<InteractionPair> interactionPairs_;

    // =========================================================
    // [New] 拓扑加速索引 (Topology Acceleration Indices)
    // =========================================================
    // Key: SolverIndex, Value: List of pointers to InteractionPair
    // 使用 vector<vector> 存储 Matrix 映射 (因为 Matrix Index 是稠密的 0~Nm)
    std::vector<std::vector<const InteractionPair*>> mat2InteractionMap_;

    // 使用 unordered_map 存储 Fracture 映射 (因为 Frac Index 是稀疏的或带偏移的)
    std::unordered_map<int, std::vector<const InteractionPair*>> frac2InteractionMap_;

    // 空列表缓存 (用于返回空引用)
    const std::vector<const InteractionPair*> emptyPairList_;

    // =========================================================
    // 内部辅助函数
    // =========================================================
    /**
     * @brief 为指定 Cell 构建唯一的棱 (Edge) 集合
     * @details 由于 Mesh 不存储 Edge，需从 Faces 动态构建并去重
     */
    void _buildLocalMatrixEdges(const Cell& cell, std::vector<MatrixEdge>& outEdges) const;

    /**
     * @brief 对散乱的交点进行排序和重构
     * @param rawPoints 输入的散乱点集
     * @param normal 参考法向 (裂缝法向)
     * @param outPair [输出] 填充其中的 points, area, center
     * @return true 如果构成了有效多边形 (Area > Tolerance)
     */
    bool _reconstructPolygon(const std::vector<Vector>& rawPoints, const Vector& normal, InteractionPair& outPair) const;
};