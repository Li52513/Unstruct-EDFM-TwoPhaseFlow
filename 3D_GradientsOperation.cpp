#include "3D_GradientsOperation.h"

// 标准库
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// =========================================================
// 引入核心类的完整定义
// =========================================================
#include "3D_MeshManager.h"      
#include "3D_FieldManager.h"     
#include "Mesh.h"                
#include "Cell.h"                
#include "Face.h"                
#include "VolField.h"            
#include "2D_FractureNetwork.h"  
#include "2D_Fracture.h"         
#include "2D_FractureEdge.h"     

namespace FVM
{
    // =========================================================
    // 内部私有结构体实现: 3x3 矩阵 (用于 LSQ)
    // =========================================================
    struct EDFM_GradientsOperation_3D::Mat3
    {
        double a[3][3];

        Mat3() {
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    a[i][j] = 0.0;
        }

        void addOuter(const Vector& r, double w) {
            const double x = r.m_x;
            const double y = r.m_y;
            const double z = r.m_z;
            a[0][0] += w * x * x; a[0][1] += w * x * y; a[0][2] += w * x * z;
            a[1][0] += w * y * x; a[1][1] += w * y * y; a[1][2] += w * y * z;
            a[2][0] += w * z * x; a[2][1] += w * z * y; a[2][2] += w * z * z;
        }

        bool solve(const Vector& b, Vector& g) const {
            double A00 = a[0][0], A01 = a[0][1], A02 = a[0][2];
            double A10 = a[1][0], A11 = a[1][1], A12 = a[1][2];
            double A20 = a[2][0], A21 = a[2][1], A22 = a[2][2];

            double det = A00 * (A11 * A22 - A12 * A21)
                - A01 * (A10 * A22 - A12 * A20)
                + A02 * (A10 * A21 - A11 * A20);

            if (std::abs(det) < 1e-20) return false;

            double invDet = 1.0 / det;

            g.m_x = (b.m_x * (A11 * A22 - A12 * A21) + b.m_y * (A02 * A21 - A01 * A22) + b.m_z * (A01 * A12 - A02 * A11)) * invDet;
            g.m_y = (b.m_x * (A12 * A20 - A10 * A22) + b.m_y * (A00 * A22 - A02 * A20) + b.m_z * (A02 * A10 - A00 * A12)) * invDet;
            g.m_z = (b.m_x * (A10 * A21 - A11 * A20) + b.m_y * (A01 * A20 - A00 * A21) + b.m_z * (A00 * A11 - A01 * A10)) * invDet;

            return true;
        }
    };

    // =========================================================
    // 1. 基岩梯度计算实现
    // =========================================================
    std::vector<Vector> EDFM_GradientsOperation_3D::ComputeMatrixGradients(
        const ::MeshManager_3D& meshMgr,
        const ::FieldManager_3D& fieldMgr,
        const std::string& scalarName,
        int smoothIters)
    {
        const ::Mesh& mesh = meshMgr.mesh();
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        auto phiPtr = fieldMgr.getMatrixScalar(scalarName);

        if (!phiPtr) {
            std::cerr << "[Warning] ComputeMatrixGradients: Field " << scalarName << " not found. Returning zero gradients.\n";
            return std::vector<Vector>(fieldMgr.numMatrixCells, Vector(0, 0, 0));
        }

        // [修正 C2440] 使用 const volScalarField&，不再强转 std::vector
        const volScalarField& phi = *phiPtr;

        std::vector<Vector> grad(fieldMgr.numMatrixCells, Vector(0, 0, 0));

        // --- LSQ + GG 循环 ---
        for (size_t i = 0; i < cells.size(); ++i)
        {
            const Cell& C = cells[i];
            double phiP = phi[i]; // volScalarField 必须支持 operator[]

            Mat3 G;
            Vector b(0, 0, 0);
            int neighborCount = 0;

            for (int globalFaceID : C.CellFaceIDs)
            {
                int fIdx = mesh.getFaceIndex(globalFaceID);
                if (fIdx < 0) continue;

                const Face& F = mesh.getFaces()[fIdx];
                int nID = (F.ownerCell == C.id) ? F.neighborCell : F.ownerCell;

                if (nID >= 0) {
                    auto it = id2idx.find(nID);
                    if (it == id2idx.end()) continue;

                    size_t nIdx = it->second;
                    Vector r = cells[nIdx].center - C.center;
                    double distSq = r.Mag() * r.Mag();

                    if (distSq > 1e-20) {
                        double w = 1.0 / distSq;
                        double dPhi = phi[nIdx] - phiP;
                        G.addOuter(r, w);
                        b = b + (r * (w * dPhi));
                        neighborCount++;
                    }
                }
            }

            bool solved = false;
            if (neighborCount >= 3) {
                G.a[0][0] += 1e-18; G.a[1][1] += 1e-18; G.a[2][2] += 1e-18;
                solved = G.solve(b, grad[i]);
            }

            if (!solved) {
                grad[i] = _computeGG_Matrix(mesh, phi, C.id);
            }
        }

        // --- 平滑 ---
        if (smoothIters > 0) {
            std::vector<Vector> tempGrad = grad;
            for (int iter = 0; iter < smoothIters; ++iter) {
                _performSmoothing(mesh, grad, tempGrad);
                grad = tempGrad;
            }
        }

        return grad;
    }

    // =========================================================
    // 2. 裂缝梯度计算实现
    // =========================================================
    std::vector<Vector> EDFM_GradientsOperation_3D::ComputeFractureGradients(
        const ::MeshManager_3D& meshMgr,
        const ::FieldManager_3D& fieldMgr,
        const std::string& scalarName)
    {
        const auto& frNet = meshMgr.fracture_network();
        const auto& edges = frNet.getAllFractureEdges();

        auto phiPtr = fieldMgr.getFractureScalar(scalarName);
        if (!phiPtr) {
            std::cerr << "[Warning] ComputeFractureGradients: Field " << scalarName << " not found.\n";
            return std::vector<Vector>(fieldMgr.numFracCells, Vector(0, 0, 0));
        }

        //
        const volScalarField& phi = *phiPtr;

        std::vector<Vector> grads(fieldMgr.numFracCells, Vector(0, 0, 0));

        // 获取尺寸的最安全方法是使用 fieldMgr，避免调用 phi.size()
        // 因为 VolField 可能是一个结构体，size 可能是成员变量而非函数
        int numPhi = static_cast<int>(fieldMgr.numFracCells);

        // --- Green-Gauss 累加 ---
        for (const auto& edge : edges)
        {
            int idO = edge.ownerCellID;
            int idN = edge.neighborCellID;

            // 移除 phi.size() 调用，使用安全的 numPhi
            if (idO < 0 || idO >= numPhi) continue;

            double valFace = phi[idO];
            if (idN >= 0 && idN < numPhi) {
                valFace = phi[idO] * (1.0 - edge.interpolationCoef) + phi[idN] * edge.interpolationCoef;
            }

            Vector Sf = edge.normal * edge.length;
            Vector flux = Sf * valFace;

            grads[idO] = grads[idO] + flux;

            if (idN >= 0 && idN < (int)grads.size()) {
                grads[idN] = grads[idN] - flux;
            }
        }

        // --- 除以体积 (Area) ---
        const auto& fractures = frNet.getFractures();
        for (const auto& frac : fractures) {
            for (const auto& cell : frac.fracCells) {
                int gid = cell.id;
                if (gid >= 0 && gid < (int)grads.size()) {
                    double area = cell.area;
                    if (area > 1e-16) {
                        grads[gid] = grads[gid] / area;
                    }
                    else {
                        grads[gid] = Vector(0, 0, 0);
                    }
                }
            }
        }

        return grads;
    }

    // =========================================================
    // 内部私有函数实现
    // =========================================================
    Vector EDFM_GradientsOperation_3D::_computeGG_Matrix(
        const ::Mesh& mesh,
        const volScalarField& phi, // [修正] 类型同步
        int cellID)
    {
        const auto& id2idx = mesh.getCellId2Index();
        size_t idxP = id2idx.at(cellID);
        const Cell& CP = mesh.getCells()[idxP];
        double phiP = phi[idxP];

        if (CP.volume <= 1e-14) return Vector(0, 0, 0);

        Vector fluxSum(0, 0, 0);

        for (int globalFaceID : CP.CellFaceIDs)
        {
            int fIdx = mesh.getFaceIndex(globalFaceID);
            if (fIdx < 0) continue;

            const Face& F = mesh.getFaces()[fIdx];
            double area = F.length; // 3D Area
            Vector Sf = F.normal * area;

            double phi_f = phiP;

            if (!F.isBoundary()) {
                int neighborID = (F.ownerCell == cellID) ? F.neighborCell : F.ownerCell;
                if (neighborID >= 0) {
                    auto itN = id2idx.find(neighborID);
                    if (itN != id2idx.end()) {
                        double phiN = phi[itN->second];
                        double w = F.f_linearInterpolationCoef;
                        phi_f = phiP * (1.0 - w) + phiN * w;
                    }
                }
            }

            if (F.ownerCell == cellID) {
                fluxSum = fluxSum + (Sf * phi_f);
            }
            else {
                fluxSum = fluxSum - (Sf * phi_f);
            }
        }

        return fluxSum * (1.0 / CP.volume);
    }

    void EDFM_GradientsOperation_3D::_performSmoothing(
        const ::Mesh& mesh,
        const std::vector<Vector>& inputGrad,
        std::vector<Vector>& outputGrad)
    {
        const auto& cells = mesh.getCells();
        const auto& id2idx = mesh.getCellId2Index();

        for (size_t i = 0; i < cells.size(); ++i) {
            const Cell& C = cells[i];
            Vector sumGrad = inputGrad[i];
            int count = 1;

            for (int globalFaceID : C.CellFaceIDs) {
                int fIdx = mesh.getFaceIndex(globalFaceID);
                if (fIdx < 0) continue;

                const Face& F = mesh.getFaces()[fIdx];
                int nID = (F.ownerCell == C.id) ? F.neighborCell : F.ownerCell;

                if (nID >= 0) {
                    auto it = id2idx.find(nID);
                    if (it != id2idx.end()) {
                        size_t nIdx = it->second;
                        sumGrad = sumGrad + inputGrad[nIdx];
                        count++;
                    }
                }
            }
            outputGrad[i] = sumGrad * (1.0 / static_cast<double>(count));
        }
    }
}