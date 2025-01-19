#pragma once
#ifndef _Variable3D_
#define _Variable3D_
#include <math.h>
#include "GeneralMethods.h"
#include <iostream>
#include <array>
template<class variableType>
class VariableTensor;

template<class variableType>
class Variable3D
{
public:
	variableType m_x;
	variableType m_y;
	variableType m_z;

public:
	Variable3D()
	{
		m_x = 0.0;
		m_y = 0.0;
		m_z = 0.0;
	}
	Variable3D(variableType x, variableType y, variableType z)
	{
		m_x = x;
		m_y = y;
		m_z = z;
	}
	Variable3D(variableType x)
	{
		m_x = x;
		m_y = x;
		m_z = x;
	}
	Variable3D operator -()
	{
		Variable3D<variableType> object;
		object.m_x = -m_x;
		object.m_y = -m_y;
		object.m_z = -m_z;
		return object;
	}

	double Mag() const
	{
		return sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
	}
	void Normalize()
	{
		double mag = Mag();
		if (GeneralMethods::isZero(mag))
		{
			m_x = m_x / mag;
			m_y = m_y / mag;
			m_z = m_z / mag;
		}
		else
		{
			std::cout << "The module of Vector is too small, please check it" << std::endl;
			exit(0);
		}
	}
	Variable3D getNormal() const
	{
		Variable3D<variableType> object;
		object.m_x = m_x;
		object.m_y = m_y;
		object.m_z = m_z;
		object.Normalize();
		return object;
	}
	Variable3D getAbs()
	{
		Variable3D<variableType> object;
		object.m_x = std::abs(m_x);
		object.m_y = std::abs(m_y);
		object.m_z = std::abs(m_z);
		return object;
	}
	Variable3D fabs_variable3D(Variable3D& var)
	{
		Variable3D<variableType> object;
		object.m_x = std::abs(var.m_x);
		object.m_y = std::abs(var.m_y);
		object.m_z = std::abs(var.m_z);
		return object;
	}

	friend Variable3D<variableType> operator + (const Variable3D<variableType>& var0, const Variable3D<variableType>& var1)
	{
		Variable3D<variableType> object;
		object.m_x = var0.m_x + var1.m_x;
		object.m_y = var0.m_y + var1.m_y;
		object.m_z = var0.m_z + var1.m_z;
		return object;
	}
	//template<class variableType>
	friend Variable3D<variableType> operator + (double var0, const Variable3D<variableType>& var1)
	{
		Variable3D<variableType> object;
		object.m_x = var0 + var1.m_x;
		object.m_y = var0 + var1.m_y;
		object.m_z = var0 + var1.m_z;
		return object;
	}
	//template<class variableType>
	friend Variable3D<variableType> operator + (const Variable3D<variableType>& var0, double var1)
	{
		Variable3D<variableType> object;
		object.m_x = var0.m_x + var1;
		object.m_y = var0.m_y + var1;
		object.m_z = var0.m_z + var1;
		return object;
	}
	//template<class variableType>
	friend Variable3D<variableType> operator - (const Variable3D<variableType>& var0, const Variable3D<variableType>& var1)
	{
		Variable3D<variableType> object;
		object.m_x = var0.m_x - var1.m_x;
		object.m_y = var0.m_y - var1.m_y;
		object.m_z = var0.m_z - var1.m_z;
		return object;
	}
	//template<class variableType>
	friend Variable3D<variableType> operator - (double var0, const Variable3D<variableType>& var1)
	{
		Variable3D<variableType> object;
		object.m_x = var0 - var1.m_x;
		object.m_y = var0 - var1.m_y;
		object.m_z = var0 - var1.m_z;
		return object;
	}
	//template<class variableType>
	friend Variable3D<variableType> operator - (const Variable3D<variableType>& var0, double var1)
	{
		Variable3D<variableType> object;
		object.m_x = var0.m_x - var1;
		object.m_y = var0.m_y - var1;
		object.m_z = var0.m_z - var1;
		return object;
	}
	//template<class variableType>
	friend Variable3D<variableType> operator * (const double coef, const Variable3D<variableType>& var)
	{
		Variable3D<variableType> object;
		object.m_x = coef * var.m_x;
		object.m_y = coef * var.m_y;
		object.m_z = coef * var.m_z;
		return object;
	}
	//template<class variableType>
	friend Variable3D<variableType> operator * (const Variable3D<variableType>& var, const double coef)
	{
		Variable3D<variableType> object;
		object.m_x = coef * var.m_x;
		object.m_y = coef * var.m_y;
		object.m_z = coef * var.m_z;
		return object;
	}

	//template<class variableType>
	friend double operator * (const Variable3D<variableType>& var0, const Variable3D<variableType>& var1)
	{
		return var0.m_x * var1.m_x + var0.m_y * var1.m_y + var0.m_z * var1.m_z;
	}
	
	friend Variable3D<variableType> operator & (const Variable3D<variableType>& var0, const Variable3D<variableType>& var1)
	{
		Variable3D<variableType> object;
		object.m_x = var0.m_y * var1.m_z - var1.m_y * var0.m_z;
		object.m_y = var0.m_z * var1.m_x - var0.m_x * var1.m_z;
		object.m_z = var0.m_x * var1.m_y - var0.m_y * var1.m_x;
		return object;
	}

	//template<class variableType>
	friend Variable3D<variableType> operator / (const Variable3D<variableType>& var, const double coef)
	{
		Variable3D<variableType> object;
		object.m_x = var.m_x / (coef + 1.0e-30);
		object.m_y = var.m_y / (coef + 1.0e-30);
		object.m_z = var.m_z / (coef + 1.0e-30);
		return object;
	}

	friend Variable3D<variableType> operator / (const Variable3D<variableType>& var1, const Variable3D<variableType>& var0)
	{
		Variable3D<variableType> object;
		if (var0.m_x != 0.0 && var0.m_y != 0.0 && var0.m_z != 0.0)
		{
			object.m_x = var1.m_x / var0.m_x;
			object.m_y = var1.m_y / var0.m_y;
			object.m_z = var1.m_z / var0.m_z;
		}
		return object;
	}

	//template<class variableType>
	friend Variable3D<variableType> operator / (const double coef, const Variable3D<variableType>& var)
	{
		Variable3D<variableType> object;
		object.m_x = coef / (var.m_x + 1.0e-30);
		object.m_y = coef / (var.m_y + 1.0e-30);
		object.m_z = coef / (var.m_z + 1.0e-30);
		return object;
	}

	//template<class variableType>
	friend bool operator == (const Variable3D<variableType>& var0, const Variable3D<variableType>& var1)
	{
		if (var0.m_x == var1.m_x && var0.m_y == var1.m_y && var0.m_z == var1.m_z)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	//template<class variableType>
	friend bool operator > (const Variable3D<variableType>& var0, const Variable3D<variableType>& var1)
	{
		if (var0.m_x > var1.m_x && var0.m_y > var1.m_y && var0.m_z > var1.m_z)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	//template<class variableType>
	friend bool operator < (const Variable3D<variableType>& var0, const Variable3D<variableType>& var1)
	{
		if (var0.m_x < var1.m_x && var0.m_y < var1.m_y && var0.m_z < var1.m_z)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}

	friend VariableTensor<variableType> operator && (Variable3D<variableType>& var0, Variable3D<variableType>& var1)
	{
		VariableTensor<variableType> object;
		object.m_xx = var0.m_x * var1.m_x;
		object.m_xy = var0.m_x * var1.m_y;
		object.m_xz = var0.m_x * var1.m_z;
		object.m_yx = var0.m_y * var1.m_x;
		object.m_yy = var0.m_y * var1.m_y;
		object.m_yz = var0.m_y * var1.m_z;
		object.m_zx = var0.m_z * var1.m_x;
		object.m_zy = var0.m_z * var1.m_y;
		object.m_zz = var0.m_z * var1.m_z;
		return object;
	}

};



#endif