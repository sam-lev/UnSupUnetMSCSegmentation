#ifndef MSC_BASIC_SUB_MESH_HANDLER
#define MSC_BASIC_SUB_MESH_HANDLER

#include "mscBasicMeshHandler.h"

class mscBasicSubMeshHandler : mscBasicMeshHandler {

public:

	virtual ~mscBasicSubMeshHandler(){
			printf("delete: mscBasicSubMeshHandler \n");
	}
	// General mesh information
	virtual DIM_TYPE share_count(CELL_INDEX_TYPE cellid) = 0;

	virtual CELL_INDEX_TYPE local_2_global(CELL_INDEX_TYPE cellid) = 0;

	virtual CELL_INDEX_TYPE global_2_local(CELL_INDEX_TYPE cellid) = 0;

};

#endif
