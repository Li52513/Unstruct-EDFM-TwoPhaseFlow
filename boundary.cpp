//#include "Bound.h"
//#include <stdexcept>
//#include <cmath>
//using namespace std;
//
//void applyBoundaryConditions
//(
//    Mesh& mesh,
//    const vector<BoundaryCondition>& conditions,
//    const Fluid& fluid,
//    const Matrix& matrix
//) 
//{
//    for (auto const& bc : conditions) {
//        int fid = bc.faceId;
//        if (fid < 1 || fid >(int)mesh.faces.size())
//            throw std::out_of_range("faceId out of range");
//
//        auto& face = mesh.faces[fid - 1];
//        int ownerCID = face.ownerCell;
//        if (ownerCID < 0) continue;
//
//        auto it = mesh.cellId2index.find(ownerCID);
//        if (it == mesh.cellId2index.end())
//            throw std::runtime_error("No cellId2index for cell");
//
//        auto& cell = mesh.cells[it->second];
//        double A_mag = face.length;
//
//        switch (bc.type) {
//        case BoundaryType::Dirichlet:
//            cell.pressure = bc.value;
//            break;
//
//        case BoundaryType::Neumann: {
//            double λ = matrix.matrix_Permeability / fluid.fluid_viscosity;
//            cell.sourceTerm += fluid.fluid_rho * λ * bc.value * A_mag;
//            break;
//        }
//
//        case BoundaryType::Robin: {
//            double a = bc.A, b = bc.B, c = bc.C;
//            double d = face.ownerToNeighbor.Mag();
//            Vector E = face.vectorE, T = face.vectorT;
//            double λ = matrix.matrix_Permeability / fluid.fluid_viscosity;
//            double BEd = b * E.Mag() / d;
//            double denom = a * A_mag + BEd;
//            double α = BEd / denom;
//            double gradT = cell.pressureGradient * T;
//            double β = (c * A_mag - b * gradT) / denom;
//            cell.faceDiscreCoef += fluid.fluid_rho * λ * (E.Mag() / d) * (1 - α);
//            cell.sourceTerm += fluid.fluid_rho * λ * (E.Mag() / d * β + gradT);
//            break;
//        }
//        }
//    }
//}

