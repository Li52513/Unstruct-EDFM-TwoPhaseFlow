#pragma once
#include "FieldRegistry.h"
#include "Mesh.h" 

template<class FieldT>
struct CellFieldView
{
	FieldT& field; // 引用的物理场变量
	const Mesh& mesh; // 引用的网格对象
	
	//按内部下标访问
	inline auto& byindex(size_t idx)
	{
		return field[idx];
	}

	inline const auto& byindex(size_t idx) const
	{
		return field[idx];
	}


	//按网格编号访问
	inline auto& byId(int cellId)
	{
		return field[mesh.getCellId2Index.at(cellId)]
	}

	inline const auto& byId(int cellId) const
	{
		return field[mesh.getCellId2Index.at(cellId)]
	}

};