#pragma once
#include "FieldRegistry.h"


//实现 *_old 场的存在性检查和创建
//若不存在 old 场就按 curr 拷贝一份；存在则不动
// 用在每个时间步开始（或第一次 t=0）确保 *_old 场可用
inline void ensureOldLayer(
	FieldRegistry& reg,
	const std::string& curr_name,
	const std::string& old_name
)
{
	auto curr = reg.get<volScalarField>(curr_name);
	if (!curr) return; // 当前时层不存在
	auto old = reg.get<volScalarField>(old_name);

	if (!old) {
		auto created = reg.getOrCreate<volScalarField>(old_name, curr->data.size(), 0.0);
		*created = *curr; // 拷贝
	}
}