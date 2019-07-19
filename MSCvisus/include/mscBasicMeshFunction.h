#ifndef MSC_BASIC_MESH_FUNCTION
#define MSC_BASIC_MESH_FUNCTION

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"

template<typename dtype>
class mscBasicMeshFunction {
   
public:
	virtual ~mscBasicMeshFunction(){
			printf("delete: mscBasicMeshFunction \n");
}
	virtual void initialize() = 0;
	virtual dtype cell_value(CELL_INDEX_TYPE cellid)  = 0;
	virtual bool less_than(CELL_INDEX_TYPE a, CELL_INDEX_TYPE b) {
		dtype av = cell_value(a);
		dtype bv = cell_value(b);
		if (av < bv) return true;
		if (bv < av) return false;
		return a < b;
	}
};

#endif
