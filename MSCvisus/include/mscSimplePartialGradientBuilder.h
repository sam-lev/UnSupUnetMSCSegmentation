#ifndef MSC_SIMPLE_PARTIAL_GRADIENT_BUILDER
#define MSC_SIMPLE_PARTIAL_GRADIENT_BUILDER

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#include "mscBasicGradientField.h"
#include "mscArrayFactory.h"

#include <vector>
#include <queue>

using namespace std;

template<typename dtype>
class mscSimplePartialGradientBuilder {
protected:
	mscBasicMeshFunction<dtype>* my_mesh_function;
    mscBasicMeshHandler* my_mesh_handler;
    mscBasicGradientField* my_grad_field;
	mscArrayFactory* my_array_factory;

	// set all cells to unassigned
	virtual void init_assigned() {
		cellIterator it;
		iteratorOperator& all_cells = my_mesh_handler->all_cells_iterator(it);
		for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
			CELL_INDEX_TYPE cellid = all_cells.value(it);
			my_grad_field->set_assigned(cellid, 0);
			my_grad_field->set_mark(cellid, 0);
		}
	}



	CELL_INDEX_TYPE doublecount;

	virtual void init_all() {
		init_assigned();
		//init_number_unpaired_facets();
	}


	void doZeroCellsSteepest() {
		minima.clear();
		cellIterator it;
		iteratorOperator& d_cells = my_mesh_handler->d_cells_iterator(0, it);
		for(d_cells.begin(it); d_cells.valid(it); d_cells.advance(it)) {
			CELL_INDEX_TYPE cellid = d_cells.value(it);
			BOUNDARY_TYPE boundary = my_mesh_handler->boundary_value(cellid);
			
			cellIterator it2;
			iteratorOperator& cofacets = my_mesh_handler->cofacets(cellid, it2);
			cofacets.begin(it2);

			CELL_INDEX_TYPE pairneg = cellid;
			CELL_INDEX_TYPE pairid = cellid;

			while(cofacets.valid(it2)) {
				CELL_INDEX_TYPE cofacetid = cofacets.value(it2);
				CELL_INDEX_TYPE neg = cellid + 2 * (cofacetid - cellid);
				if (boundary <= my_mesh_handler->boundary_value(neg)) {
					if (my_mesh_function->less_than(neg, pairneg)) {
						pairneg = neg;
						pairid = cofacetid;
					}
				}
				cofacets.advance(it2);
			}

			if (cellid == pairid) {
				//make critical
				my_grad_field->set_assigned(cellid, true);
				my_grad_field->set_critical(cellid, true);
				my_grad_field->set_dim_asc_man(cellid, 3);
				minima.push_back(cellid);
		
			} else  {
				//pair
			//printf("p(%d->%d)\n", tail, head);
				my_grad_field->set_assigned(cellid, true);
				my_grad_field->set_assigned(pairid, true);
				

				my_grad_field->set_critical(pairid, false);
				my_grad_field->set_critical(cellid, false);
				my_grad_field->set_dim_asc_man(pairid, 3);
				my_grad_field->set_dim_asc_man(cellid, 3);
				my_grad_field->set_pair(pairid, cellid);
				my_grad_field->set_pair(cellid, pairid);

			}
		}
		
	}



public:
	vector<CELL_INDEX_TYPE> minima;

	mscSimplePartialGradientBuilder(
		mscBasicMeshFunction<dtype>* mesh_function,
		mscBasicMeshHandler* mesh_handler,
		mscBasicGradientField* grad_field,
		mscArrayFactory* array_factory) : 
			my_mesh_function(mesh_function),
			my_mesh_handler(mesh_handler), 
			my_grad_field(grad_field),
			my_array_factory(array_factory) {

	}

	virtual void computeGradient() {
		init_all();
		doZeroCellsSteepest();
	}

};




#endif
