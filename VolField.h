#pragma once
#include "UserDefineVarType.h"
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include "FractureNetwork.h"

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

// 辅助函数：遍历裂缝网络中的所有裂缝微元 适用于2D-EDFM
template<class Fn>
inline void forEachFracElem(const FractureNetwork& frNet,
    const FracElemIndex& idx,
    Fn&& fn)
{
    for (std::size_t f = 0; f < frNet.fractures.size(); ++f)
    {
        const auto& F = frNet.fractures[f];
        const std::size_t base = idx.offset[f];
        for (std::size_t e = 0; e < F.elements.size(); ++e)
        {
            const std::size_t g = base + e;
            fn(f, e, g, F.elements[e]);
        }
    }
}
