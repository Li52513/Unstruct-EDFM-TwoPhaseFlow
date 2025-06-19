//#include <iostream> 
//using namespace std;
//#include <string>
//#include <vector>
//#include <algorithm>
//#include"mesh.h"
//#include"fluid.h"
//#include"fraction.h"
//#include"Matrix.h"
//#include<time.h>
//#include <fstream>
//#include <iterator>
//#include"BoundaryValuesofMatrix.h"
//#include"matrixMesheGeneration.h"
//
//const double pi = 3.14159265358979323846;
//
//
//
//struct DiscreteCoef
//{
//	vector<pair<double, int>> coef;
//	double source;
//};
//
//
//int main()
//{
//	int begintime, endtime;
//	begintime = clock();
//	////////基岩网格信息//////////////
//	int sectionNumX = 11;
//	int sectionNumY = 11;
//	double lengthX = 3.0;
//	double lengthY = 3.0;
//	double theta = pi / 2.0;
//	
//
//	
//	/////////基岩网格离散////////////
//	MatrixMesh mesh = matrixMesheGeneration(lengthX, lengthY, theta, sectionNumX, sectionNumY);
//	////////基岩边界条件设置/////////
//	BoundaryValuesofMatrix LeftBoundaryValues{ 1.0, 0.0,2e5 };
//	BoundaryValuesofMatrix RightBoundaryValues{ 1.0, 0.0,2e5 };
//	BoundaryValuesofMatrix DownBoundaryValues{ 1.0, 0.0,2e5 };
//	BoundaryValuesofMatrix UpBoundaryValues{ 1.0, 0.0,1e5 };
//}