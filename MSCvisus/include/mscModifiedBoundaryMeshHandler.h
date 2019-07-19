#ifndef MSC_MODIFIED_BOUNDARY_MESH_HANDLER
#define MSC_MODIFIED_BOUNDARY_MESH_HANDLER

#include "mscBasicMeshHandler.h"
#include "mscBasicGradientField.h"

class mscModifiedBoundaryMeshHandler : public mscBasicMeshHandler {
protected:
	mscBasicMeshHandler* mh;
	mscBasicGradientField* gf;
	DIM_TYPE md;

public:

	mscModifiedBoundaryMeshHandler(mscBasicMeshHandler* mh, mscBasicGradientField* gf) {
		this->mh = mh;
		md = mh->max_dim();
		this->gf = gf;
	}
	virtual ~mscModifiedBoundaryMeshHandler(){
			printf("delete: mscComplementMeshHandler \n");
	}
	
	
	// Queries regarding individual cells
	virtual BOUNDARY_TYPE boundary_value(CELL_INDEX_TYPE cellid) {
		BOUNDARY_TYPE bd;// = mh->boundary_value(cellid);
		//if (bd > 0) bd += (md + 1);

		// UNCOMMENT THE FOLLOWING PART TO ENABLE CONVERGENT STUFF
		bd = (/*md -*/ gf->get_dim_asc_man(cellid)) * mh->max_dim() + mh->boundary_value(cellid);
		return bd;
	}	
	
	
	
	// General mesh information
	virtual CELL_INDEX_TYPE num_cells() {
		return mh->num_cells();
	}
	virtual CELL_INDEX_TYPE num_cells(DIM_TYPE dim) {
		return mh->num_cells(dim);
	}

	virtual DIM_TYPE max_dim() {
		return md;
	}

	virtual iteratorOperator& all_cells_iterator(cellIterator& it) {
		return mh->all_cells_iterator(it);
	}

	virtual iteratorOperator& d_cells_iterator(DIM_TYPE dim, cellIterator& it) {
		return mh->d_cells_iterator(dim, it);
	}


	
	virtual DIM_TYPE dimension(const CELL_INDEX_TYPE& cellid) {
		return mh->dimension(cellid);
	}

	virtual iteratorOperator& facets(CELL_INDEX_TYPE cellid, cellIterator& it) {
		return mh->facets(cellid, it);
	}

	virtual iteratorOperator& cofacets(CELL_INDEX_TYPE cellid, cellIterator& it) {
		return mh->cofacets(cellid, it);
	}

	virtual iteratorOperator& neighbor_vertices(CELL_INDEX_TYPE cellid, cellIterator& it) {
		return mh->neighbor_vertices(cellid, it);
	}
	//virtual iteratorOperator& vertices(CELL_INDEX_TYPE cellid, cellIterator& it) {
	//	return mh->vertices(cellid, it);
	//}
	//virtual iteratorOperator& dcofaces(CELL_INDEX_TYPE cellid, cellIterator& it) {
	//	return mh->dcofaces(cellid, it);
	//}
	virtual void centroid(CELL_INDEX_TYPE cellid, float* coords) {
		mh->centroid(cellid, coords);
	}
};

#endif
