#pragma once

#include <vector>
#include <string>
#include "Mesh.h"
#include "FractureNetwork.h"
#include "BoundaryFaceClassify.h"
#include "BoundaryFaceClassify_byTag.h"
#include "FractureCommon.h"


// =================================================================================================================
// 2D EDFM (2D基岩-1D裂缝，通过1D裂缝与2D基岩网格face的交点、1D裂缝的起终点以及1D裂缝-裂缝的交点实现1D裂缝网格的划分
// =================================================================================================================

/**
 * @class MeshManager
 * @brief 统一管理基岩 & 裂缝网格及几何CI计算
 */
class MeshManager
{
public:

    // =========================================================
    // 构造与初始化
    // =========================================================

    /**
     * @brief 构造函数
     * @param lx 区域长度 X
     * @param ly 区域长度 Y
     * @param lz 区域长度 Z（0 表示纯二维）
     * @param nx X 方向划分份数
     * @param ny Y 方向划分份数
     * @param nz Z 方向划分份数
     * @param usePrism 是否使用棱柱单元（true）或四面体（false）
     * @param useQuadBase 棱柱底面保留四边形（true）或分割为三角形（false）
     */
    MeshManager(double lx, double ly, double lz, int nx, int ny, int nz, bool usePrism, bool useQuadBase); //默认有参构造函数

    // =========================================================
    // 基岩网格构建 & 预处理 (Grid Construction)
    // =========================================================
    
    //【基岩】 
    /**
     * @brief 网格生成-利用gmsh构建2D基岩非结构化网格，利用 Tag 分类并识别边界,边界会返回bcGroups_byTag_
     * @param NormalVectorCorrectionMethod：网格面非正交修正计算方法 
     * 最小修正法 MinimumCorrection；
     *OrthogonalCorrection 正交修正法；
     *OverRelaxed 超松弛修正法
    */
    void BuildSolidMatrixGrid_2D(NormalVectorCorrectionMethod corr = NormalVectorCorrectionMethod::OrthogonalCorrection);

    /**
     * @brief 接收BuildSolidMatrixGrid_*函数中返回的边界面分组bcGroups_byTag_提供读取接口,利用gmsh自生成group
    */
    const BoundaryFaceClassify_byTag::FaceGroups& boundaryFaces_byTag() const { return bcGroups_byTag_; }

    //【基岩】
    /**
     * @brief 接收BuildSolidMatrixGrid_*函数中返回的边界面个数boundaryCount_提供读取接口
    */
    size_t countBoundaryFaces() const { return boundaryCount_; }

    //【裂缝】
    /**
     * @brief 裂缝添加-输入裂缝线段的起点和终点，实现裂缝的添加
     * @param  Vector& start ： 裂缝起点
     * @param  Vector& end ：   裂缝终点
    */
    void addFracture(const Vector& start, const Vector& end); //手动插入裂缝用于调试

    //【裂缝】
    // --- DFN 裂缝模型生成 ---
    /**
     * @brief 裂缝添加-生成随机种子
     * @param  unsigned seed：设置随机种子
    */
    void setDFNRandomSeed(unsigned seed); 
  
    /**
     * @brief 裂缝添加-利用随机种子自动生成DFN裂缝网络模型
    */
    void generateDFN(int N,
        const Vector& minPoint,
        const Vector& maxPoint,
        double Lmin,
        double Lmax,
        double alpha,
        double kappa,
        bool avoidOverlap); //生成随机裂缝网络


    //【裂缝】
     /**
     * @brief 获取设置CI计算时，基岩网格单元到内部裂缝网格的距离d的度量方式AreaWeight CellCenter NodeAverage CrossAwareGauss
     * @param DistanceMetric 
     * CellCenter       使用单元中心点到线段的距离
     * NodeAverage      使用单元节点平均位置距离
     * AreaWeight       使用面积加权的数值积分平均
     * CrossAwareGauss  考虑穿越情形的高斯积分法
     */
    // 
    void setDistanceMetric(DistanceMetric m) { distanceMetric_ = m; }
    DistanceMetric distanceMetric() const { return distanceMetric_; }


    //【裂缝】
    /**
     * @brief 处理裂缝几何：交点 -> 排序 -> 细分 (Subdivide)
     * @param IntersectionSearchStrategy_2D 裂缝线段与网格面求交计算算法：
     * IntersectionSearchStrategy_2D::BruteForce,    // 方法1：暴力遍历 (无 AABB 预筛)
     * IntersectionSearchStrategy_2D::GlobalAABB,    // 方法2：全局遍历 + AABB 预筛
     * IntersectionSearchStrategy_2D::GridIndexing   // 方法3：背景网格索引 (Grid-based Indexing)
     */
    void DetectAndSubdivideFractures(IntersectionSearchStrategy_2D searchStrategy);
    
    //【裂缝】
    /**
    * @brief 计算基岩-裂缝质量交换几何系数CI_geom 和Alpha_alpha
    */
    void ComputeFractureGeometryCouplingCoefficient();


    // =========================================================
	// 核心功能：系统自由度索引管理 (System DOF Management) 按照先基岩后裂缝的顺序编号
    // =========================================================

    /**
     * @brief 构建用于构建系数方程组的索引 (Solver Indexing)
     * @details
     * 此函数负责统一规划基岩与裂缝在求解器中的行号 (Row Indices)，且裂缝编号从 n_matrix_dofs_ 开始连续递增
     * 编号策略 (Flattening Strategy):
     * 1. 基岩单元 (Matrix Cells): [0, N_m - 1]，对应cell local index
     * - 直接对应 mesh_.cells_ 的 vector 下标。
     * 2. 裂缝单元 (Fracture Elements): [N_m, N_m + N_f - 1]
     * - 遍历所有 FractureElement 并赋值其 solverIndex和parentFractureID。
     * * 前置条件:
     * - 基岩网格已生成 (mesh_.cells_ 非空)
     * - 裂缝网络已离散化 (fracture_network().elements 已生成)
     */
    void BuildGlobalSystemIndexing();

    // 构建裂缝-裂缝拓扑关系
    // 前置条件：必须先调用 BuildGlobalSystemIndexing() 确保 solverIndex 已分配
    void BuildFracturetoFractureTopology();

    /**
     * @brief 通过 SolverIndex 快速获取裂缝单元指针
     * @param solverIdx 全局求解器索引
     * @return const FractureElement* 指向裂缝单元的指针。
     * 如果 solverIdx 属于基岩范围或越界，返回 nullptr。
     */
    const FractureElement* getFractureElementBySolverIndex(int solverIdx) const;

    /**
     * @brief 获取基岩自由度总数 (Matrix DOFs)
     * @return int 基岩网格单元总数，也是裂缝索引的起始 Offset
     */
    int getMatrixDOFCount() const { return n_matrix_dofs_; }

    /**
     * @brief 获取全局总自由度数 (Total DOFs)
     * @return int Matrix + Fracture 总单元数
     */
    int getTotalDOFCount() const;

    // =========================================================
    // 输出 & 调试 (Export & Debug)
    // =========================================================
    void exportMesh(const std::string& prefix) const;
    void exportFractures(const std::string& prefix) const;
    void printFractureInfo() const;
  
    // =========================================================
    // 访问底层对象 (Accessors)
    // =========================================================
    Mesh& mesh() { return mesh_; }
    FractureNetwork& fracture_network() { return frNet_; }


    // 【旧接口】
    /**
     * @brief 网格生成-利用gmsh构建2D基岩非结构化网格，并利用几何信息进行边界识别及分类,仍然依赖BoundaryFaceClassify.h 已被 BuildSolidMatrixGrid_2D替代
     * @param NormalVectorCorrectionMethod：网格面非正交修正计算方法 最小修正法 MinimumCorrection；OrthogonalCorrection 正交修正法；OverRelaxed 超松弛修正法
     */
    void BuildSolidMatrixGrid(NormalVectorCorrectionMethod corr = NormalVectorCorrectionMethod::OrthogonalCorrection);

    //  【旧接口】
    /**
	  * @brief 接收BuildSolidMatrixGrid_*函数中返回的边界面分组bcGroups_提供读取接口，后期将被bcGroups_byTag_替代
      */
    const BoundaryFaceClassify::FaceGroups& boundaryFaces() const { return bcGroups_; }


private:
    // =========================================================
    // 私有成员变量
    // =========================================================
    Mesh mesh_; // 网格成员对象
    FractureNetwork frNet_; // 裂缝网络成员对象

    double lx_, ly_, lz_; // 区域长度成员对象
    int nx_, ny_, nz_; // 网格数量成员对象
    bool usePrism_, useQuadBase_; // 网格单元类型成员对象

    DistanceMetric distanceMetric_ = DistanceMetric::AreaWeight;
	BoundaryFaceClassify::FaceGroups bcGroups_; //旧接口，后期将被bcGroups_byTag_替代
    BoundaryFaceClassify_byTag::FaceGroups bcGroups_byTag_;
    size_t boundaryCount_ = 0;

    /** * @brief 基岩自由度总数 (Number of Matrix Degrees of Freedom)
     * @note 同时也作为裂缝单元全局编号的起始偏移量 (Offset)
     */
    int n_matrix_dofs_ = 0;

    // =========================================================
    // 全局查找表
    // =========================================================
    /**
     * @brief 全局 SolverIndex 到 FractureElement 指针的映射表
     * @details
     * Index: SolverIndex (0 ~ TotalDOFs - 1)
     * Value: FractureElement* (若是基岩DOF，则为 nullptr)
     * @note 仅在 BuildGlobalSystemIndexing 后有效。若 frNet_ 发生重分配，此指针可能失效。
     */
    std::vector<FractureElement*> globalSolverIndexToElementPtr_;

};