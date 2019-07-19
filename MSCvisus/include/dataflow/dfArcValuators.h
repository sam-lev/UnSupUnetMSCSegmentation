#ifndef DF_ARC_VALUATORS
#define DF_ARC_VALUATORS

#include "mscBasicMSC.h"

using namespace std;

template <typename dtype>
class dfArcValuator {
public:
	virtual float Value(arc<dtype>* a) { return 0; }
};

template <typename dtype>
class dfArcValuatorEndPoint {
public:
	enum ENDPOINTTYPE = { LOWER, UPPER };
	ENDPOINTTYPE mEPType;
	dfArcValuatorEndPoint() : mEPType(LOWER) {}
	void SetEPType(ENDPOINTTYPE type) { mEPType = type; }
	
	virtual float Value(arc<dtype>* a) { 
		if (mEPType == LOWER) {
			return (float) a->lower->value;
		} else {
			return (float) a->upper->value;
		}
	}
};




#endif
