#pragma once
#ifndef _VariableTensor_
#define _VariableTensor_
#include "Variable3D.hpp"

template<class variableType>
class VariableTensor
{
public:
	variableType m_xx;
	variableType m_xy;
	variableType m_xz;
	variableType m_yx;
	variableType m_yy;
	variableType m_yz;
	variableType m_zx;
	variableType m_zy;
	variableType m_zz;

public:
	VariableTensor()
	{
		m_xx = 0.0;
		m_xy = 0.0;
		m_xz = 0.0;
		m_yx = 0.0;
		m_yy = 0.0;
		m_yz = 0.0;
		m_zx = 0.0;
		m_zy = 0.0;
		m_zz = 0.0;
	}
	VariableTensor(variableType xx, variableType xy, variableType xz,
		variableType yx, variableType yy, variableType yz,
		variableType zx, variableType zy, variableType zz)
	{
		m_xx = xx;
		m_xy = xy;
		m_xz = xz;
		m_yx = yx;
		m_yy = yy;
		m_yz = yz;
		m_zx = zx;
		m_zy = zy;
		m_zz = zz;
	}
	VariableTensor(variableType x)
	{
		m_xx = x;
		m_xy = x;
		m_xz = x;
		m_yx = x;
		m_yy = x;
		m_yz = x;
		m_zx = x;
		m_zy = x;
		m_zz = x;
	}

	friend VariableTensor<variableType> operator + (const VariableTensor<variableType>& var0, const VariableTensor<variableType>& var1)
	{
		VariableTensor<variableType> object;
		object.m_xx = var0.m_xx + var1.m_xx;
		object.m_xy = var0.m_xy + var1.m_xy;
		object.m_xz = var0.m_xz + var1.m_xz;
		object.m_yx = var0.m_yx + var1.m_yx;
		object.m_yy = var0.m_yy + var1.m_yy;
		object.m_yz = var0.m_yz + var1.m_yz;
		object.m_zx = var0.m_zx + var1.m_zx;
		object.m_zy = var0.m_zy + var1.m_zy;
		object.m_zz = var0.m_zz + var1.m_zz;
		return object;
	}
	//template<class variableType>
	friend VariableTensor<variableType> operator - (const VariableTensor<variableType>& var0, const VariableTensor<variableType>& var1)
	{
		VariableTensor<variableType> object;
		object.m_xx = var0.m_xx - var1.m_xx;
		object.m_xy = var0.m_xy - var1.m_xy;
		object.m_xz = var0.m_xz - var1.m_xz;
		object.m_yx = var0.m_yx - var1.m_yx;
		object.m_yy = var0.m_yy - var1.m_yy;
		object.m_yz = var0.m_yz - var1.m_yz;
		object.m_zx = var0.m_zx - var1.m_zx;
		object.m_zy = var0.m_zy - var1.m_zy;
		object.m_zz = var0.m_zz - var1.m_zz;
		return object;
	}
	//template<class variableType>
	friend VariableTensor<variableType> operator / (const VariableTensor<variableType>& var0, const double var1)
	{
		if (var1 == 0)
		{
			std::cout << "the denominator is zero, please check it!" << std::endl;
			system("pause");
			exit(0);
		}
		VariableTensor<variableType> object;
		object.m_xx = var0.m_xx / (var1 + 1.0e-30);
		object.m_xy = var0.m_xy / (var1 + 1.0e-30);
		object.m_xz = var0.m_xz / (var1 + 1.0e-30);
		object.m_yx = var0.m_yx / (var1 + 1.0e-30);
		object.m_yy = var0.m_yy / (var1 + 1.0e-30);
		object.m_yz = var0.m_yz / (var1 + 1.0e-30);
		object.m_zx = var0.m_zx / (var1 + 1.0e-30);
		object.m_zy = var0.m_zy / (var1 + 1.0e-30);
		object.m_zz = var0.m_zz / (var1 + 1.0e-30);
		return object;
	}
	friend VariableTensor<variableType> operator / (const double var1, const VariableTensor<variableType>& var0)
	{
		VariableTensor<variableType> object;
		object.m_xx = var1 / (var0.m_xx + 1.0e-30);
		object.m_xy = var1 / (var0.m_xy + 1.0e-30);
		object.m_xz = var1 / (var0.m_xz + 1.0e-30);
		object.m_yx = var1 / (var0.m_yx + 1.0e-30);
		object.m_yy = var1 / (var0.m_yy + 1.0e-30);
		object.m_yz = var1 / (var0.m_yz + 1.0e-30);
		object.m_zx = var1 / (var0.m_zx + 1.0e-30);
		object.m_zy = var1 / (var0.m_zy + 1.0e-30);
		object.m_zz = var1 / (var0.m_zz + 1.0e-30);
		return object;
	}

	//template<class variableType>
	friend Variable3D<variableType> operator * (const Variable3D<variableType>& var, const VariableTensor<variableType> tensor)
	{
		Variable3D<variableType> object;
		object.m_x = var.m_x * tensor.m_xx + var.m_y * tensor.m_yx + var.m_z * tensor.m_zx;
		object.m_y = var.m_x * tensor.m_xy + var.m_y * tensor.m_yy + var.m_z * tensor.m_zy;
		object.m_z = var.m_x * tensor.m_xz + var.m_y * tensor.m_yz + var.m_z * tensor.m_zz;
		return object;
	}
	//template<class variableType>
	friend Variable3D<variableType> operator * (const VariableTensor<variableType> tensor, const Variable3D<variableType>& var)
	{
		Variable3D<variableType> object;
		object.m_x = var.m_x * tensor.m_xx + var.m_y * tensor.m_xy + var.m_z * tensor.m_xz;
		object.m_y = var.m_x * tensor.m_yx + var.m_y * tensor.m_yy + var.m_z * tensor.m_yz;
		object.m_z = var.m_x * tensor.m_zx + var.m_y * tensor.m_zy + var.m_z * tensor.m_zz;
		return object;
	}

	//template<class variableType>
	friend VariableTensor<variableType> operator * (const VariableTensor<variableType>& var0, const VariableTensor<variableType>& var1)
	{
		VariableTensor<variableType> object;
		object.m_xx = var0.m_xx * var1.m_xx + var0.m_xy * var1.m_yx + var0.m_xz * var1.m_zx;
		object.m_xy = var0.m_xx * var1.m_xy + var0.m_xy * var1.m_yy + var0.m_xz * var1.m_zy;
		object.m_xz = var0.m_xx * var1.m_xz + var0.m_xy * var1.m_yz + var0.m_xz * var1.m_zz;
		object.m_yx = var0.m_yx * var1.m_xx + var0.m_yy * var1.m_yx + var0.m_yz * var1.m_zx;
		object.m_yy = var0.m_yx * var1.m_xy + var0.m_yy * var1.m_yy + var0.m_yz * var1.m_zy;
		object.m_yz = var0.m_yx * var1.m_xz + var0.m_yy * var1.m_yz + var0.m_yz * var1.m_zz;
		object.m_zx = var0.m_zx * var1.m_xx + var0.m_zy * var1.m_yx + var0.m_zz * var1.m_zx;
		object.m_zy = var0.m_zx * var1.m_xy + var0.m_zy * var1.m_yy + var0.m_zz * var1.m_zy;
		object.m_zz = var0.m_zx * var1.m_xz + var0.m_zy * var1.m_yz + var0.m_zz * var1.m_zz;
		return object;
	}

	friend VariableTensor<variableType> operator * (const double var0, const VariableTensor<variableType>& var1)
	{
		VariableTensor<variableType> object;
		object.m_xx = var1.m_xx * var0;
		object.m_xy = var1.m_xy * var0;
		object.m_xz = var1.m_xz * var0;
		object.m_yx = var1.m_yx * var0;
		object.m_yy = var1.m_yy * var0;
		object.m_yz = var1.m_yz * var0;
		object.m_zx = var1.m_zx * var0;
		object.m_zy = var1.m_zy * var0;
		object.m_zz = var1.m_zz * var0;
		return object;
	}

	friend VariableTensor<variableType> operator * (const VariableTensor<variableType>& var1, const double var0)
	{
		VariableTensor<variableType> object;
		object.m_xx = var1.m_xx * var0;
		object.m_xy = var1.m_xy * var0;
		object.m_xz = var1.m_xz * var0;
		object.m_yx = var1.m_yx * var0;
		object.m_yy = var1.m_yy * var0;
		object.m_yz = var1.m_yz * var0;
		object.m_zx = var1.m_zx * var0;
		object.m_zy = var1.m_zy * var0;
		object.m_zz = var1.m_zz * var0;
		return object;
	}
};

#endif