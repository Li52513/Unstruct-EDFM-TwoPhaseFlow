#pragma once
#include <vector>
#include <algorithm>
#include "UserDefineVarType.h"

struct AABB
{
	Vector min;// вСоб╫г
	Vector max;//срио╫г
	AABB() = default;
	AABB(const Vector& p1, const Vector& p2)
	{
		min.m_x = std::min(p1.m_x, p2.m_x);
		min.m_y = std::min(p1.m_y, p2.m_y);
		min.m_z = std::min(p1.m_z, p2.m_z);
		max.m_x = std::max(p1.m_x, p2.m_x);
		max.m_y = std::max(p1.m_y, p2.m_y);
		max.m_z = std::max(p1.m_z, p2.m_z);
	}
	bool overlaps(const AABB& other) const
	{
		return (min.m_x <= other.max.m_x && max.m_x >= other.min.m_x &&
			min.m_y <= other.max.m_y && max.m_y >= other.min.m_y &&
			min.m_z <= other.max.m_z && max.m_z >= other.min.m_z);
	}
};