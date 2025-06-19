//#pragma once
//#include <vector>
//#include "Mesh.h"
//#include "Fluid.h"   // <-- 新增
//#include "Matrix.h"  // <-- 新增
//using namespace std;
//
//enum class BoundaryType { Dirichlet, Neumann, Robin };
//
//struct BoundaryCondition 
//{
//    BoundaryType  type;
//    int           faceId;
//    double        value;
//    double        A, B, C;
//    BoundaryCondition
//    (
//        BoundaryType t,
//        int          f,
//        double       v,
//        double       a = 0.0,
//        double       b = 0.0,
//        double       c = 0.0
//    ) : type(t), faceId(f), value(v), A(a), B(b), C(c) 
//    {
//    }
//};
//
//
//void applyBoundaryConditions
//(
//    Mesh& mesh,
//    const vector<BoundaryCondition>& conditions,
//    const Fluid& fluid,
//    const Matrix& matrix
//);

