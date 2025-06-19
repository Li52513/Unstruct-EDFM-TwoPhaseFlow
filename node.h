#pragma once
#include <vector>
#include <algorithm>
#include "UserDefineVarType.h"
using namespace std;

class Node 
{
public:
    int id;         // 节点编号
    Vector coord;   // 节点坐标

	Node(int id, const Vector& coord) : id(id), coord(coord) {}// 有参构造函数

    Node() : id(-1), coord(0.0, 0.0, 0.0) {}// 初始化列表

 
};

//class MatrixNode
//{
//public:
//	///the Geometry&Mesh information of the node///////
//	int n_ID;					// node ID
//	Vector n_position;			// node position
//	vector<int> n_inCellID;		// cell ID included in the node
//	vector<int> n_inFaceID;		// face ID included in the node
//
//	///the Physical information of the node///////
//	double node_pressureValue;	// pressure value at the node 
//
//public:
//	MatrixNode();
//};


