#pragma once
#include <string>
#include <vector>
#include <cmath>

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
	double p_w0 = 8e6; //基岩内部初始水相压力，单位 Pa
	double dp_wdx = 0.0; //基岩压力梯度，单位 Pa/m
	double dp_wdy = 0.0; //基岩压力梯度，单位 Pa/m
	double dp_wdz = 0.0; //基岩压力梯度，单位 Pa/m

	double T0 = 573.15; //基岩初始温度，单位 K
	double dTdx = 0.0; //基岩温度梯度，单位 K/m
	double dTdy = 0.0; //基岩温度梯度，单位 K/m
	double dTdz = 0.0; //基岩温度梯度，单位 K/m

	double s_w = 0.05; //基岩初始水相饱和度 (最好等于1-Sgr)

	double x0 = 0;
	double x_dx = 0;
	double x_dy = 0;
	double x_dz = 0;
};


