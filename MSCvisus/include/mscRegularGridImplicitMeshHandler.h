#ifndef MSC_REGULAR_GRID_IMPLICIT_MESH_HANDLER
#define MSC_REGULAR_GRID_IMPLICIT_MESH_HANDLER

#include "mscIndexTypes.h"
#include "mscBasicIterator.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicIterator.h"

#include <map>

using namespace std;


class mscRegularGridImplicitMeshHandler : public mscBasicMeshHandler {

public:

	mscRegularGridImplicitMeshHandler() {}


	 virtual  ~mscRegularGridImplicitMeshHandler() {
		  		printf("delete: mscRegularGrid2DImplicitMeshHandler \n");
	  }

	 virtual CELL_INDEX_TYPE offset_2_pair(unsigned char offset) {return 0;}

	 virtual unsigned char pair_2_offset(CELL_INDEX_TYPE diff) { return 0 ;}

	
	  // Queries regarding individual cells
	 virtual inline void cellid_2_coords(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE* coords) {}  
	  
};

#endif