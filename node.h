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
    Node(int id = -1, const Vector& coord = Vector(0.0, 0.0, 0.0)) : id(id), coord(coord) {}//构造函数及初始化列表 
};



