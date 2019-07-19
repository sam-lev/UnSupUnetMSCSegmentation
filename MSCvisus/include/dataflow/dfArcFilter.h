#ifndef DF_ARC_FILTER
#define DF_ARC_FILTER

#include <vector>
#include dfProducerTypes

using namespace std;

template <typename dtype>
class dfArcTest {
public:
	virtual bool Test(arc<dtype>* a) { return false; }
};


template <typename dtype>
class dfArcFilter : public dfArcProducer {
public:
	dfArcTest* mTest;
	dfArcFilter() : mTest(NULL){}

	// set the test function
	void SetTest(dfArcTest* t) { mTest = t; }

	virtual void Compute() {
		this->clear();
		if (mTest == NULL) {
			printf("RUNTIME ERROR: dfArcFilter Test == NULL\n");
			return;
		}
		mValid = true;
		for (int i = 0; i < upstream.size(); i++) {
			if (! upstream[i]->Valid())
				upstream[i]->Compute();
			dfArcProducer* ap = (dfArcProducer*) upstream[i];
			
			for (dfArcProducer::iterator it = ap->begin(); it != ap->end(); it++) {
				arc<dtype>* a = *it;
				if (mTest(a)) this->push_back(a);
			}
		}
	}
};


#endif
