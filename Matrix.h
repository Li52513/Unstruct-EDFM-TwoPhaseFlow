#pragma once
#include"PhysicalPropertiesManager.h"
class Matrix 
{
public:
    double matrix_Porosity;   // 基岩的孔隙度
    double matrix_Permeability; // 基岩的渗透率
    double matrix_Beta;        // 基岩的某个物理常数（比如与流动相关的常数）

    // 构造函数，初始化基岩物性参数
    Matrix(double porosity = 0.2, double permeability = 1e-12, double beta = 1.0);

    // 设置基岩物性参数
    void setProperties(double porosity, double permeability, double beta);
   

};
