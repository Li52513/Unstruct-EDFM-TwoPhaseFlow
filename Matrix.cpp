#include "Matrix.h"

// 构造函数，初始化基岩物性参数
Matrix::Matrix(double porosity, double permeability, double beta)
    : matrix_Porosity(porosity), matrix_Permeability(permeability), matrix_Beta(beta) 
{
}

// 设置基岩物性参数
void Matrix::setProperties(double porosity, double permeability, double beta) 
{
    matrix_Porosity = porosity;
    matrix_Permeability = permeability;
    matrix_Beta = beta;
}




