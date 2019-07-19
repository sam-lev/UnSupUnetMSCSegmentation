#ifndef MSC_SIMPLE_LINE_CONSTRAINED_GRADIENT_BUILDER
#define MSC_SIMPLE_LINE_CONSTRAINED_GRADIENT_BUILDER

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#include "mscBasicGradientField.h"
#include "mscArrayFactory.h"
#include "mscSimpleGradientBuilder.h"
#include "mscPreClassifier.h"

#include <vector>
#include <queue>

using namespace std;

template<typename dtype>
class mscSimpleLineConstrainedGradientBuilder {
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
			my_grad_field->set_dim_asc_man(cellid, 0);
			//INT_TYPE temp_num = 0;
			//cellIterator fit;
			//iteratorOperator& facets = my_mesh_handler->facets(cellid, fit);
			//for (facets.begin(fit); facets.valid(fit); facets.advance(fit)) {
			//	temp_num++;
			//}
			//my_grad_field->set_num_unpaired_facets(cellid, temp_num);
		}
	}

	
	virtual void init_dsc() {
		cellIterator it;

		iteratorOperator& zero_cells = my_mesh_handler->d_cells_iterator(0, it);
		for (zero_cells.begin(it); zero_cells.valid(it); zero_cells.advance(it)) {
			
			// the id of the 0-cell
			CELL_INDEX_TYPE cellid = zero_cells.value(it);
			int value =  (mClassifier->getId(cellid) > 0 ? 1 : 0);
			my_grad_field->set_dim_asc_man(cellid, mClassifier->getId(cellid));
		}

		iteratorOperator& one_cells = my_mesh_handler->d_cells_iterator(1, it);
		for (one_cells.begin(it); one_cells.valid(it); one_cells.advance(it)) {
			
			// the id of the 1-cell
			CELL_INDEX_TYPE cellid = one_cells.value(it);
			
			// the "vertices" of the 1-cell
			cellIterator fit;
			iteratorOperator& facets = my_mesh_handler->facets(cellid, fit);
			facets.begin(fit);	//set to first vertex
			CELL_INDEX_TYPE vertex1 = facets.value(fit);
			facets.advance(fit);  // next vertex
			CELL_INDEX_TYPE vertex2 = facets.value(fit);

			// look up value
			
			int equals = 0;
			if (mClassifier->getId(vertex1) != 0 && mClassifier->getId(vertex2) != 0) equals = 1;

			my_grad_field->set_dim_asc_man(cellid, equals);
		}
	}
	
	void set_cell_dsc_count(const CELL_INDEX_TYPE& cellid, const int& value) {
		int v1 = my_grad_field->get_dim_asc_man(cellid);
		if (value > 0) {
			v1 = v1 + 2; // equivalent to v1 OR value >> 1
		}
		my_grad_field->set_dim_asc_man(cellid, v1);
	}
	
	int get_cell_dsc_count(const CELL_INDEX_TYPE& cellid) {
		int v1 = my_grad_field->get_dim_asc_man(cellid);
		return v1 >> 1;
	}


	virtual void init_asc() {
		
		cellIterator it;
		iteratorOperator& n_cells = my_mesh_handler->d_cells_iterator(my_mesh_handler->max_dim()-1, it);
		for (n_cells.begin(it); n_cells.valid(it); n_cells.advance(it)) {
			CELL_INDEX_TYPE cellid = n_cells.value(it);		
			if (my_mesh_handler->boundary_value(cellid) > 0) {
				set_cell_dsc_count(cellid, 1);
				continue;
			}
			int value =  (mClassifier->getId(cellid) > 0 ? 1 : 0);
			my_grad_field->set_dim_asc_man(cellid, value);
		}

		iteratorOperator& n2_cells = my_mesh_handler->d_cells_iterator(my_mesh_handler->max_dim()-2, it);
		for (n2_cells.begin(it); n2_cells.valid(it); n2_cells.advance(it)) {
			CELL_INDEX_TYPE cellid = n2_cells.value(it);		

			bool all_b = true;
			cellIterator fit;
			iteratorOperator& cofacets = my_mesh_handler->cofacets(cellid, fit);
		
			for (cofacets.begin(fit); cofacets.valid(fit); cofacets.advance(fit)) {
				CELL_INDEX_TYPE cofacet = cofacets.value(fit);
				if (mClassifier->getId(cofacet) == 0) {
					all_b = false;
				}
			}

			if (all_b) {
				int v = my_grad_field->get_dim_asc_man(cellid);
				my_grad_field->set_dim_asc_man(cellid, v+1);
			}
		}

		// now the 1 cells are 1 if boundary, and 0 if not
		// so for each 2 and (3) cell only check if facets are 1 or 0 
		for (int i = my_mesh_handler->max_dim()-2; i >= 0; i--) {
			cellIterator it2;
			iteratorOperator& d_cells = my_mesh_handler->d_cells_iterator(i, it2);
			for (d_cells.begin(it2); d_cells.valid(it2); d_cells.advance(it2)) {
				CELL_INDEX_TYPE cellid = d_cells.value(it2);
				cellIterator fit;
				iteratorOperator& cofacets = my_mesh_handler->cofacets(cellid, fit);
				int value = 0;
				for (cofacets.begin(fit); cofacets.valid(fit); cofacets.advance(fit)) {
					CELL_INDEX_TYPE cofacet = cofacets.value(fit);
					if (get_cell_dsc_count(cofacet) != 0) {
						value = 1;
						break;
					}
				}
				set_cell_dsc_count(cellid, value);
			}
		}
	}

	// set asc/dsc boundary for all cells
	virtual void init_asc_dsc() {
		init_asc();
		// NOW ASCENDING INDICES ARE DONE!!!
		init_dsc();
	}

	virtual void init_all() {
		init_assigned();
		init_asc_dsc();
	}

	struct vertex_record {
		CELL_INDEX_TYPE cellid;
		dtype value;
	};

	struct vertex_comparator {
		bool operator()(const vertex_record& a, const vertex_record& b) {
			if (a.value < b.value) return false;
			if (a.value > b.value) return true;
			return a.cellid > b.cellid;
		}
	};

	priority_queue<vertex_record, vector<vertex_record>, vertex_comparator> sorted_vertices;

	void enqueue_vertices() {
		cellIterator it;
		iteratorOperator& verts = my_mesh_handler->d_cells_iterator(0, it);
		for (verts.begin(it); verts.valid(it); verts.advance(it)) {
			CELL_INDEX_TYPE cellid = verts.value(it);
			vertex_record r;
			r.cellid = cellid;
			r.value = my_mesh_function->cell_value(cellid);
			sorted_vertices.push(r);
		}

	}


	void process_lower_star(const CELL_INDEX_TYPE& vert) {



	}



	void assign_gradients() {

		while(! sorted_vertices.empty()) {
			CELL_INDEX_TYPE vert = sorted_vertices.top().cellid;
			sorted_vertices.pop();
			process_lower_star(vert);
		}

	}

	mscPreClassifier* mClassifier;
public:

	mscSimpleLineConstrainedGradientBuilder(
		mscPreClassifier* classifier,
		mscBasicMeshFunction<dtype>* mesh_function,
		mscBasicMeshHandler* mesh_handler,
		mscBasicGradientField* grad_field,
		mscArrayFactory* array_factory) :			
			my_mesh_function(mesh_function),
			my_mesh_handler(mesh_handler), 
			my_grad_field(grad_field),
			my_array_factory(array_factory),
			mClassifier(classifier) { 		
	}

	virtual void computeGradient() {
		init_all();
		//enqueue_vertices();
		//assign_gradients();	
	}

};




#endif
