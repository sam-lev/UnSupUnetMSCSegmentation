#ifndef DF_ARC_TESTS
#define DF_ARC_TESTS

#include <vector>
#include "mscBasicMSC.h"
#include "dfArcValuators.h"

using namespace std;

template <typename dtype>
class dfArcTest {
public:
	virtual bool Test(arc<dtype>* a) { return false; }
};

template <typename dtype>
class dfArcRangeTest {

public:


};




#endif
