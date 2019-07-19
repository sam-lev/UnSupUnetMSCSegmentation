#ifndef MSC_REGULAR_3D_GRID_IMPLICIT_SUB_MESH_HANDLER
#define MSC_REGULAR_3D_GRID_IMPLICIT_SUB_MESH_HANDLER

#include "mscRegularGrid3DImplicitMeshHandler.h"
#include "mscBasicSubMeshHandler.h"


using namespace std;


class mscRegularGrid3DImplicitSubMeshHandler : public mscRegularGrid3DImplicitMeshHandler, 
	public mscBasicSubMeshHandler {

protected:
	
	CELL_INDEX_TYPE gX, gY, gZ; // global sizes
	CELL_INDEX_TYPE sX, sY, sZ; // the starts of the current block

	mscRegularGrid3DImplicitMeshHandler* mGlobalMesh;



public:

	mscRegularGrid3DImplicitSubMeshHandler(CELL_INDEX_TYPE x, CELL_INDEX_TYPE y, CELL_INDEX_TYPE z,
		CELL_INDEX_TYPE gx, CELL_INDEX_TYPE gy, CELL_INDEX_TYPE gz,
		CELL_INDEX_TYPE sx, CELL_INDEX_TYPE sy, CELL_INDEX_TYPE sz) :
		gX(gx), gY(gy), gZ(gz), sX(sx), sY(sy), sZ(sz), mscRegularGrid3DImplicitMeshHandler(x,y,z) {

	  }


	 virtual  ~mscRegularGrid3DImplicitSubMeshHandler() {
		  		printf("delete: mscRegularGrid3DImplicitSubMeshHandler \n");
	  }

	 
	 virtual DIM_TYPE share_count(CELL_INDEX_TYPE cellid) {
		 return mscRegularGrid3DImplicitMeshHandler::boundary_value(cellid);
	 }

	 virtual CELL_INDEX_TYPE local_2_global(CELL_INDEX_TYPE cellid) {
		return cellid;
	 }
	 
	 virtual CELL_INDEX_TYPE global_2_local(CELL_INDEX_TYPE cellid) {
		return cellid;
	 }

	  friend class mscRegularGrid3DGradientField;
	  
};

#endif