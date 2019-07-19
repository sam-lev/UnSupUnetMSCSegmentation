#ifndef MSC_BASIC_MESH_HANDLER
#define MSC_BASIC_MESH_HANDLER

#include "mscIndexTypes.h"
#include "mscBasicIterator.h"

class mscBasicMeshHandler;

struct cellIterator {

	DIM_TYPE my_dim;
	mscBasicMeshHandler* my_mesh_handler;
	CELL_INDEX_TYPE my_base;
	CELL_INDEX_TYPE my_location;
	CELL_INDEX_TYPE my_neighbor;
	CELL_INDEX_TYPE my_total_size;
	CELL_INDEX_TYPE my_coords[3];

};

class iteratorOperator {

public:

	virtual void begin(cellIterator& it) {
		it.my_location = 0;
		it.my_total_size = 0;
	}

	virtual bool valid(cellIterator& it) {
		return it.my_location < it.my_total_size;
	}
	virtual void advance(cellIterator& it) {
		it.my_location++;
	}
	virtual void advance(const CELL_INDEX_TYPE num, cellIterator& it) {
		for (int i = 0; this->valid(it) && i < num; i++) this->advance(it);
	}
	virtual CELL_INDEX_TYPE value(cellIterator& it) {
		return it.my_location;
	}
	virtual CELL_INDEX_TYPE index(cellIterator& it) {
		return it.my_location;
	}
};

class mscBasicMeshHandler {

public:

	virtual ~mscBasicMeshHandler(){
			printf("delete: mscBasicMeshHandler \n");
}
	// General mesh information
	virtual CELL_INDEX_TYPE num_cells() = 0;
	virtual CELL_INDEX_TYPE num_cells(DIM_TYPE dim) = 0;
	virtual DIM_TYPE max_dim() = 0;
	virtual iteratorOperator& all_cells_iterator(cellIterator& it) = 0;
	virtual iteratorOperator& d_cells_iterator(DIM_TYPE dim, cellIterator& it) = 0;
	
	// Queries regarding individual cells
	virtual BOUNDARY_TYPE boundary_value(CELL_INDEX_TYPE cellid) = 0;
	virtual DIM_TYPE dimension(const CELL_INDEX_TYPE& cellid) = 0;
	virtual iteratorOperator& facets(CELL_INDEX_TYPE cellid, cellIterator& it) = 0;
	virtual iteratorOperator& cofacets(CELL_INDEX_TYPE cellid, cellIterator& it) = 0;
	virtual iteratorOperator& neighbor_vertices(CELL_INDEX_TYPE cellid, cellIterator& it) = 0;
	//virtual iteratorOperator& vertices(CELL_INDEX_TYPE cellid, cellIterator& it) = 0;
	//virtual iteratorOperator& dcofaces(CELL_INDEX_TYPE cellid, cellIterator& it) = 0;
	virtual void centroid(CELL_INDEX_TYPE cellid, float* coords) = 0; 




};

#endif
