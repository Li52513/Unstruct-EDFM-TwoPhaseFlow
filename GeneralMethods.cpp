#include "GeneralMethods.h"
#include <stdlib.h>
#include <iostream> 

GeneralMethods::GeneralMethods()
{

}
bool GeneralMethods::isZero(double var)
{
	double minValue = 1.0e-20;
	if (abs(var) < minValue)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}