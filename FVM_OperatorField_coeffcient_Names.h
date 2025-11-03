#pragma once
#include <string>

// 你的命名规则：aPP_conv_<tag> / aPN_conv_<tag> / bP_conv_<tag>
struct OperatorFieldNames 
{
    std::string a_f_diff, s_f_diff;             //扩散项
	std::string aPP_conv, aPN_conv, bP_conv;    //对流项
	std::string a_time, b_time;                 //时间项

};