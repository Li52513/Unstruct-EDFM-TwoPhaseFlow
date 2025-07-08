#pragma once
#include "mesh.h"
#include "FractureNetwork.h"
#include "BoundaryClassify.h"

/**
 * @class MeshManager
 * @brief Unified mesh management for solid matrix and fracture network including CI calculation and property assignment
 */
class MeshManager 
{
public:
    // ========================= CONSTRUCTION =========================
    /**
     * @brief Constructor
     * @param lx Domain length in X direction
     * @param ly Domain length in Y direction  
     * @param lz Domain length in Z direction (0 for 2D cases)
     * @param nx Number of divisions in X direction
     * @param ny Number of divisions in Y direction
     * @param nz Number of divisions in Z direction
     * @param usePrism Whether to use prism elements (true) or tetrahedron elements (false)
     * @param useQuadBase Whether to use quadrilateral base (true) or triangular base (false)
     */
    MeshManager(double lx, double ly, double lz, int nx, int ny, int nz, bool usePrism, bool useQuadBase);

    // ========================= MESH OPERATIONS =========================
    /**
     * @brief Build solid matrix grid and preprocessing
     * @param corr Normal vector correction method
     */
    void BuildSolidMatrixGrid(NormalVectorCorrectionMethod corr = NormalVectorCorrectionMethod::OrthogonalCorrection);

    // ========================= FRACTURE OPERATIONS =========================
    /**
     * @brief Add single fracture manually
     * @param start Fracture start point
     * @param end Fracture end point
     */
    void addFracture(const Vector& start, const Vector& end);

    /**
     * @brief Set random seed for DFN fracture generation
     * @param seed Random seed value
     */
    void setDFNRandomSeed(unsigned seed);
    
    /**
     * @brief Generate fractures using DFN (Discrete Fracture Network) method
     * @param N Number of fractures to generate
     * @param minPoint Minimum bounds for fracture centers
     * @param maxPoint Maximum bounds for fracture centers  
     * @param Lmin Minimum fracture length
     * @param Lmax Maximum fracture length
     * @param alpha Length power law exponent
     * @param kappa von Mises concentration parameter
     * @param avoidOverlap Whether to avoid overlap with existing fractures
     */
    void generateDFN(int N,
        const Vector& minPoint,
        const Vector& maxPoint,
        double Lmin,
        double Lmax,
        double alpha,
        double kappa,
        bool avoidOverlap);

    /**
     * @brief Detect fracture intersections and subdivide fractures
     */
    void DetectAndSubdivideFractures();

    /**
     * @brief Compute geometric coupling coefficients for fractures
     */
    void ComputeFractureGeometryCouplingCoefficient();
   
    // ========================= EXPORT & DEBUG OPERATIONS =========================
    /**
     * @brief Export mesh data to text files
     * @param prefix File prefix for output files
     */
    void exportMesh(const std::string& prefix) const;
    
    /**
     * @brief Export fracture data to text files
     * @param prefix File prefix for output files
     */
    void exportFractures(const std::string& prefix) const;
    
    /**
     * @brief Print fracture information for debugging
     */
    void printFractureInfo() const;
    
    /**
     * @brief Print CI source terms for debugging
     */
    void printCISourceTerms();

    // ========================= ACCESSORS =========================
    /**
     * @brief Get reference to underlying mesh
     * @return Reference to Mesh object
     */
    Mesh& mesh() { return mesh_; }
    
    /**
     * @brief Get reference to fracture network
     * @return Reference to FractureNetwork object
     */
    FractureNetwork& fracture_network() { return frNet_; }

private:
    // ========================= MEMBER VARIABLES =========================
    Mesh mesh_;                    ///< Solid matrix mesh
    FractureNetwork frNet_;        ///< Fracture network
    double lx_, ly_, lz_;          ///< Domain dimensions
    int nx_, ny_, nz_;             ///< Grid divisions
    bool usePrism_, useQuadBase_;  ///< Element type options
};

