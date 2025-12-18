#pragma once
#include "Mesh.h"
#include "FractureNetwork.h"
#include "BoundaryFaceClassify.h"

/**
 * @class MeshManager
 * @brief 统一管理基岩 & 裂缝网格及几何CI计算
 */

class MeshManager
{
public:
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

    // —— 网格构建 & 预处理 ——————————————————————————————————
    void BuildSolidMatrixGrid(NormalVectorCorrectionMethod corr = NormalVectorCorrectionMethod::OrthogonalCorrection);


    // —— 裂缝几何 & CI 计算 —————————————————————————————
    /// 添加裂缝（起点→终点）
    void addFracture(const Vector& start, const Vector& end); //手动插入裂缝用于调试

    //DFN裂缝模型
    void setDFNRandomSeed(unsigned seed);  // 设置随机种子
    void generateDFN(int N,
        const Vector& minPoint,
        const Vector& maxPoint,
        double Lmin,
        double Lmax,
        double alpha,
        double kappa,
        bool avoidOverlap); //生成随机裂缝网络


    /// 处理裂缝几何：交点→排序→subdivide
    void DetectAndSubdivideFractures(bool useAABBAccelerationctrl);
    void setDistanceMetric(DistanceMetric m) { distanceMetric_ = m; }
    DistanceMetric distanceMetric() const { return distanceMetric_; }

    /// 只算几何耦合（调用 computeGeometryCoupling）
    void ComputeFractureGeometryCouplingCoefficient();

    ///访问边界面
    const BoundaryFaceClassify::FaceGroups& boundaryFaces() const { return bcGroups_; }
    // 统计所有边界面总数（自检用）
    size_t countBoundaryFaces() const { return boundaryCount_; }

    // —— 输出 & 调试 ——————————————————————————————————
    void exportMesh(const std::string& prefix) const;
    void exportFractures(const std::string& prefix) const;
    void printFractureInfo() const;
    void printCISourceTerms();

    // —— 访问底层 Mesh & FractureNetwork ——————————————————
    Mesh& mesh() { return mesh_; }
    FractureNetwork& fracture_network() { return frNet_; }



private:
    Mesh mesh_; // 网格成员对象
    FractureNetwork frNet_; // 裂缝网络成员对象
    double lx_, ly_, lz_; // 区域长度成员对象
    int nx_, ny_, nz_; // 网格数量成员对象
    bool usePrism_, useQuadBase_; // 网格单元类型成员对象
    DistanceMetric distanceMetric_ = DistanceMetric::AreaWeight;
    BoundaryFaceClassify::FaceGroups bcGroups_;
    size_t boundaryCount_ = 0;

};