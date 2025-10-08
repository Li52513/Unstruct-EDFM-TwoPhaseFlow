#pragma once
#include <vector>
#include <map>
#include "Node.h"
#include "Face.h"
#include "Cell.h"
#include <fstream>
#include <unordered_map>
#include "VectorHash.h"

class Mesh 
{
public:
 
    //------------------构造函数----------------------//
    Mesh();
    // --------------------- 只读参数 ------------------- //
	const vector<Node>& getNodes() const { return nodes_; } // 提供只读访问
	const vector<Face>& getFaces() const { return faces_; } // 提供只读访问
	vector<Cell>& getCells() { return cells_; } // 提供访问  考虑到物性会更新，所以是可修改的
    const std::vector<Cell>& getCells() const { return cells_; }  // 新增
	const unordered_map<int, Node>& getNodesMap() const { return nodesMap_; } // 提供只读访问
	const unordered_map<int, int>& getCellId2Index() const { return cellId2index_; } // 提供只读访问

	//==============基岩网格空间划分函数声明=================//
	void buildFaceBins(); // 构建面空间划分
	
    vector<int>getCandidateFacesFromBins(const AABB& box) const;

    
	int getGridCount() const { return gridCount_; } // 获取网格总数

    vector<reference_wrapper<const Cell>> getInnerCells() const;
    vector<reference_wrapper<const Cell>> getBoundaryCells() const;
    
    //------------------其他接口------------------------------//
    double getCellArea(int cellID)  const;

	//-------------------成员函数--------------------------//
   
    void BuildMesh(double lengthX, double lengthY, double lengthZ, int nx, int ny, int nz,bool usePrism, bool useQuadBase);
    void ClassifySolidMatrixCells();
    void printMeshInfo();
    void exportToTxt(const string& prefix) const;     
    void ComputeSolidMatrixMeshFaceGeometricInfor(NormalVectorCorrectionMethod method);   
   
    // （可选）Dirichlet 边界压力  //待优化

    //-----------------------网格边界识别--------------------//
  

private:
    vector<Node> nodes_;
	vector<Face> faces_;
	vector<Cell> cells_;
    unordered_map<int, Node> nodesMap_;
    unordered_map<int, int> cellId2index_;

	int binCountX_=0; // X 方向网格划分数
	int binCountY_ = 0; // Y 方向网格划分数
    double binSizeX_ = 0;   // X 方向实际 bin 宽度
    double binSizeY_ = 0;   // Y 方向实际 bin 高度

	unordered_map<int, vector<int>> faceBins_; // 面空间划分的 bin
 
	int gridCount_ = 0; // 网格总数

};

