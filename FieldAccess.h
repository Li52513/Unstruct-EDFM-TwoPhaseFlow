#pragma once
#include <limits>
#include <string>
#include "FieldRegistry.h"
#include "VolField.h"   

namespace FieldAccess
{
	inline double getScalar(const FieldRegistry& R, const std::string& name, size_t i)
	{
		auto field = R.get<volScalarField>(name);
		if (!field) return std::numeric_limits<double>::quiet_NaN();
		const auto& v = field->data;  // 与 VolField.h 对齐（data）
		if (i >= v.size()) return std::numeric_limits<double>::quiet_NaN();
		return v[i];

	}


	//对计算时间步长为时间步长得到的物理场变量进行插值得到以输出时间步长为时间步长得到的物理场变量场：val = (1-α)*prev + α*curr
	inline double getScalarLerp(const FieldRegistry& Prev, const FieldRegistry& Curr, const std::string& name, size_t i, double alpha)
	{
		double a = getScalar(Prev, name, i);
		double b = getScalar(Curr, name, i);
		if (!(a == a) || !(b == b)) return std::numeric_limits<double>::quiet_NaN(); //检测是否为NaN
		return a + (b - a) * alpha;

	}

}