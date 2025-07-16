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
	vector<Cell>& getCells() { return cells_; } // 提供科协访问  考虑到物性会更新，所以是可修改的
	const unordered_map<int, Node>& getNodesMap() const { return nodesMap_; } // 提供只读访问
	const unordered_map<int, int>& getCellId2Index() const { return cellId2index_; } // 提供只读访问
	int getGridCount() const { return gridCount_; } // 获取网格总数

    vector<reference_wrapper<const Cell>> getInnerCells() const;
    vector<reference_wrapper<const Cell>> getBoundaryCells() const;
	int getGhostStartIndex() const { return ghostStartIndex_; } // 获取 Ghost Cell 起始索引   
    
    //------------------其他接口------------------------------//
    double getCellArea(int cellID)  const;

	//-------------------成员函数--------------------------//
   
    void BuildMesh(double lengthX, double lengthY, double lengthZ, int nx, int ny, int nz,bool usePrism, bool useQuadBase);
    void ClassifySolidMatrixCells();
    void printMeshInfo();
    void exportToTxt(const string& prefix) const;     
    void ComputeSolidMatrixMeshFaceGeometricInfor(NormalVectorCorrectionMethod method);   
   
    // （可选）Dirichlet 边界压力  //待优化
    void CreateSolidMatrixGhostCells();    // 新增：生成 Ghost Cell
    void exportGhostInfo(const std::string& prefix) const;

private:
    vector<Node> nodes_;
	vector<Face> faces_;
	vector<Cell> cells_;
    unordered_map<int, Node> nodesMap_;
    unordered_map<int, int> cellId2index_;
	int gridCount_ = 0; // 网格总数
    int ghostStartIndex_ = 0;

};

