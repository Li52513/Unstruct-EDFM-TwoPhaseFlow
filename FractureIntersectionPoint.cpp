#include "FractureIntersectionPoint.h"

// 有参构造函数的实现
FractureIntersectionPoint::FractureIntersectionPoint(
	int _id,
	const Vector& _point,
	int _edgeID,
	double _param,
	bool _isFF,
	int _globalFFID,
	IntersectionOrigin _origin)
	: id(_id),
	point(_point),
	edgeID(_edgeID),
	param(_param),
	isFF(_isFF),
	globalFFID(_globalFFID),
	origin(_origin)
{

}