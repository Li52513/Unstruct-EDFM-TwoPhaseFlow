/*---------------------------------------------------------------------------*\
*                                 Node 数据类型
* 作用：    读取gmsh生成的网格节点，并存储节点编号和坐标信息。
* 
* 成员变量： 节点编号 （id）;  节点坐标（coord）
* 
* 构造
\*---------------------------------------------------------------------------*/
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



