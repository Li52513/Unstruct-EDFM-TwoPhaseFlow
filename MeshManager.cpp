#include "MeshManager.h"
#include "FractureNetwork.h"
#include <iostream>

using namespace std;

// ========================= CONSTRUCTION =========================

MeshManager::MeshManager(double lx, double ly, double lz,
    int nx, int ny, int nz,
    bool usePrism, bool useQuadBase)
    : lx_(lx), ly_(ly), lz_(lz),
    nx_(nx), ny_(ny), nz_(nz),
    usePrism_(usePrism), useQuadBase_(useQuadBase)
{
    cout << "Initializing MeshManager constructor" << "\n";
}

// ========================= MESH OPERATIONS =========================

void MeshManager::BuildSolidMatrixGrid(NormalVectorCorrectionMethod corr) 
{
    mesh_.BuildMesh(lx_, ly_, lz_, nx_, ny_, nz_, usePrism_, useQuadBase_);
    mesh_.ClassifySolidMatrixCells();
    BoundaryClassify::ClassifySolidMatrixMeshBoundaryFaces(mesh_, lx_, ly_);
    mesh_.CreateSolidMatrixGhostCells();
    mesh_.ComputeSolidMatrixMeshFaceGeometricInfor(corr);
}

// ========================= FRACTURE OPERATIONS =========================

void MeshManager::addFracture(const Vector& s, const Vector& e) 
{
    frNet_.addFracture(s, e);
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

void MeshManager::DetectAndSubdivideFractures() 
{
    // 1) Detect fracture-fracture intersections
    frNet_.DetectFracturetoFractureIntersections();
    
    // 2) Detect fracture-mesh intersections
    for (auto& F : frNet_.fractures)
        F.DetectFracturetoMeshFaceIntersections(mesh_.getFaces(), mesh_.getCells(), mesh_.getNodesMap());
   
    // 3) Deduplicate and renumber global fracture-fracture intersections
    frNet_.DeduplicateAndRenumberFractureToFractureIntersections();
    
    // 4) Distribute fracture-fracture intersections to individual fractures
    frNet_.DistributeFracture_FractureIntersectionsToGlobalInersections();

    // 5) Sort intersections and subdivide fractures
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

// ========================= EXPORT & DEBUG OPERATIONS =========================

void MeshManager::exportMesh(const std::string& pref) const 
{
    mesh_.exportToTxt(pref);
}

void MeshManager::exportFractures(const std::string& pref) const 
{
    frNet_.exportToTxt(pref);
}

void MeshManager::printFractureInfo() const 
{
    frNet_.printFractureInfo();
}

void MeshManager::printCISourceTerms()  
{
    std::cout << "\n=== Matrix-Fracture CI Source Terms ===\n";
    for (auto& cell : mesh_.getCells())
    {
        if (cell.id < 0) continue;               // Skip Ghost cells
        if (cell.CI.empty())     continue;       // Skip cells without fracture connections

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
