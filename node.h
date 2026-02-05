#pragma once
#include <vector>
#include <algorithm>
#include "UserDefineVarType.h"
using namespace std;

class Node 
{
public:
    // =========================================================
    // 1. 成员属性 (Member Properties)
    // =========================================================

    int id;         // 节点编号
    Vector coord;   // 节点坐标

    // =========================================================
    // 2. 成员行为 (Member Behaviors)
    // =========================================================

    ///有参构造函数
    Node(int id = -1, const Vector& coord = Vector(0.0, 0.0, 0.0));
};



