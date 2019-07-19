#ifndef MSC_NEGATING_MESH_FUNCTION
#define MSC_NEGATING_MESH_FUNCTION

#include "mscBasicMeshFunction.h"

template<typename dtype>
class mscNegatingMeshFunction : public mscBasicMeshFunction<dtype> {
   mscBasicMeshFunction<dtype>* mf;
public:
	mscNegatingMeshFunction(mscBasicMeshFunction<dtype>* mf) {
		this->mf = mf;
	}
	virtual ~mscNegatingMeshFunction(){
			printf("delete: mscNegatingMeshFunction \n");
	}
	virtual void initialize() {};
	virtual dtype cell_value(CELL_INDEX_TYPE cellid) {

		return -1 * mf->cell_value(cellid);
	}

	virtual bool less_than(CELL_INDEX_TYPE a, CELL_INDEX_TYPE b) {
		dtype av = cell_value(a);
		dtype bv = cell_value(b);
		if (av < bv) return true;
		if (bv < av) return false;
		return a > b; // NOTE THAT WE FLIP INDEX COMPARISON TOO!!!!!!
	}

};

#endif
