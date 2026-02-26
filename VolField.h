#pragma once
#include "UserDefineVarType.h"
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include "FractureNetwork.h"
#include "ADVar.hpp"

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

// =========================================================
// FIM 自动微分场类型别名 (AD Field Type Aliases)
// =========================================================
/**
 * @brief 存储自动微分变量的体心物理场
 * @tparam N 独立自变量的数量 (如单相流 N=2, 两相流 N=3)
 * @details 用于 FIM 框架下未知量及其衍生状态方程物性 (P, T, Sw, Rho, Mu 等) 的存储
 */
template<int N>
using volADField = VolField<ADVar<N>>;

/**
 * @brief 存储自动微分变量的面心/边心物理场
 * @tparam N 独立自变量的数量
 * @details 用于 FIM 框架下通量 (Flux)、迎风格式流度等面心物理量的存储
 */
template<int N>
using faceADField = VolField<ADVar<N>>;

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
