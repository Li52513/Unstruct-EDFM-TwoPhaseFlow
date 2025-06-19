#pragma once
#include "Fracture.h"
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <unordered_map>

class FractureNetwork 
{
public:

    //储存裂缝与裂缝交点的信息
    struct GlobalFFPoint
    {
        int id;               // 全局编号
        Vector point;         // 交点坐标
        int fracA, fracB;     // 发生交点的两条裂缝索引
        double paramA, paramB;// 在线段 [0,1] 上的归一化位置
    };

	vector<Fracture> fractures;  // 裂缝集合
	vector<GlobalFFPoint> globalFFPts; // 裂缝-裂缝交点
   /* void computeFractureFractureTI(const Fluid& fluid);*/
    // 添加一条裂缝
    void addFracture(const Vector& start, const Vector& end);
    // 主处理函数：检测裂缝-裂缝交点 + 裂缝-边界交点 + 分段
    void processFractures(const vector<Face>& meshFaces, const vector<Cell>& meshCells, const map<int, Node>& meshNodes, const Fluid& fluid, const Matrix& matrix, Mesh& mesh);

    ///@brief 将 globalFFPts 中的每个裂缝–裂缝交点分发至对应两条裂缝的 intersections 中
    void DistributeFracture_FractureIntersectionsToGlobalInersections();

    // 输出所有裂缝和交点信息
    void printFractureInfo() const;

    // 为所有裂缝设置物性参数
    /*void setFractureProperties(int fractureID, double porosity, double permeability,
        double compressibility, double aperture)*/;
    
    void exportToTxt(const string& prefix) const;

    void DetectFracturetoFractureIntersections();
    void DeduplicateAndRenumberFractureToFractureIntersections();
    bool isClose(const Vector& a, const Vector& b, double tol = 1e-6) const;

private:
    int nextFracID_{ 0 };          ///< 递增的裂缝 ID
};