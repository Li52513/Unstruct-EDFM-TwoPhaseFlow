#pragma once
#ifndef GeneralMethods_HEADER
#define GeneralMethods_HEADER
#include <cstdlib>
#include <string>      // ? 引入 std::string
#include <iostream>    // ? 引入 std::cerr
class GeneralMethods
{
public:
	GeneralMethods();
	static bool isZero(double var);
	
	//0423改动
	static void warn(const std::string& msg) 
	{
		std::cerr << "**[EDFM-Warn]** " << msg << '\n';
	}

};

#endif

