#pragma once
#include "Fracture.h"
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <unordered_map>
#include <cstdlib>   // rand, srand
#include <ctime>     // time
#include <cmath>     // cos, sin, M_PI


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
	

    void setRandomSeed(unsigned seed);  //unsigned 是什么数据类型？

    //@brief 生成随机 DFN 裂缝网络
    /**
     * @brief 基于 DFN 方法随机生成裂缝
     * @param N            要生成的裂缝数量
     * @param minPoint     裂缝中心坐标下限 (x,y,z)
     * @param maxPoint     裂缝中心坐标上限 (x,y,z)
     * @param Lmin         裂缝最小长度
     * @param Lmax         裂缝最大长度
     * @param alpha        长度幂律指数 (p(L) ∝ L^{-α})
     * @param kappa        von Mises 浓度 (κ≈0→均匀取向)
     * @param avoidOverlap 是否简单避让已有裂缝重叠
     */
    void generateDFN(int N, const Vector& minPoint, const Vector& maxPoint, double Lmin, double Lmax, double alpha, double kappa, bool avoidOverlap);
    
    //@brief 添加裂缝
    void addFracture(const Vector& start, const Vector& end);  
    
    //@brief 主处理函数：检测裂缝-裂缝交点 + 裂缝-边界交点 + 分段
    void processFractures(const vector<Face>& meshFaces, const vector<Cell>& meshCells, const unordered_map<int, Node>& meshNodes, const Fluid& fluid, const Matrix& matrix, Mesh& mesh);
    
    //@brief 将 globalFFPts 中的每个裂缝–裂缝交点分发至对应两条裂缝的 intersections 中
    void DistributeFracture_FractureIntersectionsToGlobalInersections();  
    
    //@brief 输出所有裂缝和交点信息
    void printFractureInfo() const; 
	
    //@brief 检测裂缝–裂缝交点
    void DetectFracturetoFractureIntersections(); 
	
    //@brief 去重并重新编号裂缝–裂缝交点
    void DeduplicateAndRenumberFractureToFractureIntersections(); 
    
    //@brief 导出裂缝网络信息到文本文件
    void exportToTxt(const string& prefix) const;  
	
    //@brief 判断两个点是否在给定的公差范围内相等
    bool isClose(const Vector& a, const Vector& b, double tol = 1e-6) const; 

private:
    int nextFracID_{ 0 };          ///< 递增的裂缝 ID
};