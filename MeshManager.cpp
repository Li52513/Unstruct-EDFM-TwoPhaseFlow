#include "MeshManager.h"
#include "FractureNetwork.h"

MeshManager::MeshManager(double lx, double ly, double lz,
    int nx, int ny, int nz,
    bool usePrism, bool useQuadBase)
    : lx_(lx), ly_(ly), lz_(lz),
    nx_(nx), ny_(ny), nz_(nz),
    usePrism_(usePrism), useQuadBase_(useQuadBase)
{
    cout<<"调用MeshManager构造函数并初始化"<< "\n";
}

void MeshManager::BuildSolidMatrixGrid(NormalVectorCorrectionMethod corr) 
{
    mesh_.BuildMesh(lx_, ly_, lz_, nx_, ny_, nz_, usePrism_, useQuadBase_);
    mesh_.ClassifySolidMatrixCells();
    BoundaryClassify::ClassifySolidMatrixMeshBoundaryFaces(mesh_, lx_, ly_);
    mesh_.CreateSolidMatrixGhostCells();
    mesh_.ComputeSolidMatrixMeshFaceGeometricInfor(corr);
}

void MeshManager::addFracture(const Vector& s, const Vector& e) 
{
    frNet_.addFracture(s, e);
}

void MeshManager::DetectAndSubdivideFractures() 
{
    // 1) 先找裂缝–裂缝的交点
    frNet_.DetectFracturetoFractureIntersections();
	// 2) 然后找裂缝–网格面的交点
    for (auto& F : frNet_.fractures)
        F.DetectFracturetoMeshFaceIntersections(mesh_.getFaces(), mesh_.getCells(), mesh_.getNodesMap());
   
    // 3) 给全局裂缝–裂缝交点去重并编号
    frNet_.DeduplicateAndRenumberFractureToFractureIntersections();
    
    // 4) 把全局裂缝–裂缝交点也插入到每条裂缝的 intersections 中
    frNet_.DistributeFracture_FractureIntersectionsToGlobalInersections();

    // 5) 重新排序 & 编号，然后划分
    for (auto& F : frNet_.fractures)
    {
        F.sortAndRenumberIntersections();
        F.subdivide(mesh_.getCells(), mesh_.getNodesMap());
    }
}

void MeshManager::ComputeFractureGeometryCouplingCoefficient()
{
    for (auto& F : frNet_.fractures)
        F.computeGeometryCouplingCoefficientgeomCIandgeomAlpha();
}

void MeshManager::setDFNRandomSeed(unsigned seed)
{
	frNet_.setRandomSeed(seed);
}

void MeshManager::generateDFN(int N,
    const Vector& minPoint,
    const Vector& maxPoint,
    double Lmin,
    double Lmax,
    double alpha,
    double kappa,
    bool avoidOverlap)
{
	frNet_.generateDFN(N, minPoint, maxPoint, Lmin, Lmax, alpha, kappa, avoidOverlap);
}

void MeshManager::exportMesh(const std::string& pref) const {
    mesh_.exportToTxt(pref);
}

void MeshManager::exportFractures(const std::string& pref) const {
    frNet_.exportToTxt(pref);
}

void MeshManager::printFractureInfo() const {
    frNet_.printFractureInfo();
}

void MeshManager::printCISourceTerms()  
{
    std::cout << "\n=== 检查基岩–裂缝 CI 源项 ===\n";
    for (auto& cell : mesh_.getCells())
    {
        if (cell.id < 0) continue;               // 跳过 Ghost
        if (cell.CI.empty())     continue;       // 没有 fracture 交互就跳过

        std::cout << "Cell " << cell.id << ":\n";
        for (size_t k = 0; k < cell.CI.size(); ++k) {
            std::cout
                << "  SegID=" << cell.CI_belongFraction[k]
                << "  CI=" << cell.CI[k]
                << "  p_fr(cell)=" << cell.fracturePressure[k]
                << "\n";
        }
    }
}
