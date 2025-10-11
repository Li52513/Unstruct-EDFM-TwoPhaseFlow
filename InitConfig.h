#pragma once
#include <string>
#include <vector>
#include <cmath>

struct Gravity 
{
	double gx = 0.0;
	double gy = 0.0;
	double gz = -9.81; // 默认重力加速度指向负Z方向，单位 m/s²
};

//////////VG模型参数与相对渗透率参数//////////
struct VGParams  //VG模型参数 Se=(1+（alpha*P）^n)^(-m) 表征毛细压力-饱和度关系
{
	double alpha = 1.0 / 5e4;  // [1/Pa] 例：5e4 Pa 尺度
	double n = 2.0;            // 无量纲
	double m() const { return 1.0 - 1.0 / n; }
	double Swr = 0.2;       // 残余水
	double Sgr = 0.0;       // 残余气（CO2）

};

struct RelPermParams
{
	double L = 0.5;  //无量纲Mualem 指数
};



///////////////基岩物性参数//////////
struct RockDefaults  //当前已经进行了分区处理
{
	double phi = 0.15;          // 孔隙度
	double k_iso = 1e-13;      //基岩绝对渗透率，单位 m²
	double rho_r = 2650.0;      //基岩密度，单位 kg/m³
	double cp_r = 1000.0;      //基岩比热容，单位 J/(kg·K)
	double lambda_r = 2.5;   //基岩导热系数，单位 W/(m·K)
};

struct InitFields
{
	//基础场：uniform+线性梯度（可都为0）
	double p0 = 6e6; //基岩内部初始水相压力，单位 Pa
	double dpdx = 1000.0; //基岩压力梯度，单位 Pa/m
	double dpdy = 0.0; //基岩压力梯度，单位 Pa/m
	double dpdz = 0.0; //基岩压力梯度，单位 Pa/m

	double T0 = 373.15; //基岩初始温度，单位 K
	double dTdx = 0.0; //基岩温度梯度，单位 K/m
	double dTdy = 0.0; //基岩温度梯度，单位 K/m
	double dTdz = 0.0; //基岩温度梯度，单位 K/m

	double sw0 = 0.9; //基岩初始水相饱和度 (最好等于1-Sgr)

};

/////////////////运行选项开关//////////////////
struct RuntypeSwitch
{
	bool enableWater = true;
	bool enableCO2 = true;
	bool enableEnergy = true;
	bool cacheFluidPropsAsFields = true; //是否将流体物性参数作为场变量缓存

};


