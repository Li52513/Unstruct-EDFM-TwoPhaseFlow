#pragma once
#include "UserDefineVarType.h"
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
using namespace std;
struct BaseField
{
	string name;
	size_t size; //存储该物理场变量的个数
	explicit BaseField (string n, size_t nCells)
		: name(move(n)), size(nCells) {}
	virtual ~BaseField() = default; //虚析构函数
}; 

template<typename T>
struct VolField : BaseField
{
	vector<T> data; //所储存的物理量
	VolField(const string& n, size_t nCells, const T& init = T{}) : BaseField(n, nCells), data(nCells, init) {}
	inline  T& operator[] (size_t i) { return data[i]; } //重载下标运算符
	inline const T& operator[] (size_t i) const { return data[i]; } //重载下标运算符

};
using volScalarField = VolField<double>;
using volVectorField = VolField<Vector>;
