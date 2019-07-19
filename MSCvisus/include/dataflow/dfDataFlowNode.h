#ifndef DF_DATA_FLOW_NODE
#define DF_DATA_FLOW_NODE

#include <vector>

using namespace std;

class dfDataFlowNode {
protected:
	bool mValid;
public:

	dfDataFlowNode() : mValid(false) {}

	// a list of who is upstream and downstream
	vector<DF_DATA_FLOW_NODE*> downstream;
	vector<DF_DATA_FLOW_NODE*> downstream;
	
	// downstream data flow nodes ask if the data in this producer 
	// is valid
	virtual bool Valid() {
		return mValid;
	}
	
	// invalidate downstream consumers
	virtual void Invalidate() { 
		mValid=false;
		for (int i = 0; i < downstream.size(); i++)
			downstream[i]->Invalidate();
	}

	// downstream data flow nodes will ask for this data to be 
	// recomputed if it is not valid;
	virtual void Compute() {
		mValid = true;
		for (int i = 0; i < upstream.size(); i++)
			if (! upstream[i]->Valid())
				upstream[i]->Compute();
	}

	// add a consumer
	virtual void AddConsumer(dfDataFlowNode* ds) {
		downstream.push_back(ds);
	}

	// always add the producers to the consumers!
	virtual void AddProducer(dfDataFlowNode* ds) {
		upstream.push_back(ds);
		ds->AddConsumer(this);
	}
};



#endif
