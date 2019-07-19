#ifndef MSC_COMPLEMENT_MESH_HANDLER
#define MSC_COMPLEMENT_MESH_HANDLER

#include "mscBasicMeshHandler.h"


class mscComplementMeshHandler : public mscBasicMeshHandler {
protected:
	mscBasicMeshHandler* mh;
	DIM_TYPE md;

public:

	mscComplementMeshHandler(mscBasicMeshHandler* mh) {
		this->mh = mh;
		md = mh->max_dim();
	}
	virtual ~mscComplementMeshHandler(){
			printf("delete: mscComplementMeshHandler \n");
	}
	// General mesh information
	virtual CELL_INDEX_TYPE num_cells() {
		return mh->num_cells();
	}
	virtual CELL_INDEX_TYPE num_cells(DIM_TYPE dim) {
		return mh->num_cells(md - dim);
	}

	virtual DIM_TYPE max_dim() {
		return md;
	}

	virtual iteratorOperator& all_cells_iterator(cellIterator& it) {
		return mh->all_cells_iterator(it);
	}

	virtual iteratorOperator& d_cells_iterator(DIM_TYPE dim, cellIterator& it) {
		return mh->d_cells_iterator(md - dim, it);
	}

	// Queries regarding individual cells
	virtual BOUNDARY_TYPE boundary_value(CELL_INDEX_TYPE cellid) {
		return mh->boundary_value(cellid);
	}
	
	virtual DIM_TYPE dimension(const CELL_INDEX_TYPE& cellid) {
		return md - mh->dimension(cellid);
	}

	virtual iteratorOperator& facets(CELL_INDEX_TYPE cellid, cellIterator& it) {
		return mh->cofacets(cellid, it);
	}

	virtual iteratorOperator& cofacets(CELL_INDEX_TYPE cellid, cellIterator& it) {
		return mh->facets(cellid, it);
	}
	virtual iteratorOperator& neighbor_vertices(CELL_INDEX_TYPE cellid, cellIterator& it) {
		return mh->neighbor_vertices(cellid, it);
	}
	//virtual iteratorOperator& vertices(CELL_INDEX_TYPE cellid, cellIterator& it) {
	//	return mh->dcofaces(cellid, it);
	//}
	//virtual iteratorOperator& dcofaces(CELL_INDEX_TYPE cellid, cellIterator& it) {
	//	return mh->vertices(cellid, it);
	//}
	virtual void centroid(CELL_INDEX_TYPE cellid, float* coords) {
		mh->centroid(cellid, coords);
	}

 
};

#endif
