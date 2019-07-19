#ifndef DF_PRODUCER_TYPES
#define DF_PRODUCER_TYPES

#include <vector>

using namespace std;

template <typename dtype>
class dfNodeProducer : public dfDataFlowNode, public vector<node<dtype>*> {
};

template <typename dtype>
class dfArcProducer : public dfDataFlowNode, public vector<arc<dtype>*> {
};

class dfGenericCellProducer : public dfDataFlowNode, public vector<mscGenericCell*> {
};



#endif
