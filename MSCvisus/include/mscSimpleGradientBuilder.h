#ifndef MSC_SIMPLE_GRADIENT_BUILDER
#define MSC_SIMPLE_GRADIENT_BUILDER

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
class mscSimpleGradientBuilder {
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
	// set number of unpaired facets for all cells
	virtual void init_number_unpaired_facets() {
		cellIterator it;
		iteratorOperator& all_cells = my_mesh_handler->all_cells_iterator(it);
		for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
			CELL_INDEX_TYPE cellid = all_cells.value(it);
			INT_TYPE temp_num = 0;
			cellIterator fit;
			iteratorOperator& facets = my_mesh_handler->facets(cellid, fit);
			for (facets.begin(fit); facets.valid(fit); facets.advance(fit)) {
				temp_num++;
			}
			my_grad_field->set_num_unpaired_facets(cellid, temp_num);
		}
	}

	CELL_INDEX_TYPE doublecount;

	virtual void init_all() {
		init_assigned();
		init_number_unpaired_facets();
	}

	// comparison object for queue
	struct comparison_element {
		CELL_INDEX_TYPE cellid;
		CELL_INDEX_TYPE insertion_time;
		dtype value;
		BOUNDARY_TYPE boundary;
	};

	// comparator of objects 
	struct element_comparator {
		// return true if b comes first, then a
		bool operator() (const comparison_element& a, const comparison_element& b) {

			// then values
			if (a.value < b.value) return false;
			if (a.value > b.value) return true;

			// first boundaries
			if (a.boundary < b.boundary) return true;
			if (a.boundary > b.boundary) return false;
			
			// then insertion time, later first
			if (a.insertion_time > b.insertion_time) return true;
			if (a.insertion_time < b.insertion_time) return false;

			return a.cellid > b.cellid;
		}
	};

	struct oneleft_element {
		CELL_INDEX_TYPE cellid;
		BOUNDARY_TYPE boundary;
		DIM_TYPE dim;
		dtype value;
		INT_TYPE weight;
	};

	struct oneleft_comparator {
		bool operator() (const oneleft_element& a, const oneleft_element& b) {
			// first boundaries
			if (a.boundary < b.boundary) return true;
			if (a.boundary > b.boundary) return false;
			
			// then dim
			if (a.dim < b.dim) return false;
			if (a.dim > b.dim) return true;
			
			// then values
			if (a.value < b.value) return false;
			if (a.value > b.value) return true;
			
			// then weight, lower first
			if (a.weight < b.weight) return false;
			if (a.weight > b.weight) return true;

			return a.cellid > b.cellid;			

		}
	};

	// simple simulation of simplicity comparator
	//virtual bool less_than(const CELL_INDEX_TYPE& a, const CELL_INDEX_TYPE& b) {
	//	if (my_mesh_function->cell_value(a) < my_mesh_function->cell_value(b)) return true;
	//	if (my_mesh_function->cell_value(a) > my_mesh_function->cell_value(b)) return false;
	//	printf("INDEX BREAKING TIES\n");
	//	return a < b;
	//}

	virtual bool f_greater_than(const CELL_INDEX_TYPE& a, const CELL_INDEX_TYPE& b) {
		return my_mesh_function->cell_value(a) > my_mesh_function->cell_value(b);
	}

	virtual bool less_than_all_unassigned_cofacets(const CELL_INDEX_TYPE& cellid) {
	
		BOUNDARY_TYPE boundary = my_mesh_handler->boundary_value(cellid);
		cellIterator it;
		iteratorOperator& cofacets = my_mesh_handler->cofacets(cellid, it);
		cofacets.begin(it);

		while(cofacets.valid(it)) {
			CELL_INDEX_TYPE cofacetid = cofacets.value(it);
			if (my_grad_field->get_assigned(cofacetid) == 0 && 
				boundary == my_mesh_handler->boundary_value(cofacetid)) {
				//if (f_greater_than(cofacetid, cellid)) {
				//	cofacets.advance(it);
				//	continue;
				//} else {

					cellIterator it2;
					iteratorOperator& facets = my_mesh_handler->facets( cofacetid, it2);
					facets.begin(it2);
					while(facets.valid(it2)) {
						CELL_INDEX_TYPE facetid = facets.value(it2);
						if (facetid != cellid &&
							my_grad_field->get_assigned(facetid) == 0 &&
							! my_mesh_function->less_than(cellid, facetid))
							return false;
						facets.advance(it2);
					}
				//}
			}
			cofacets.advance(it);
		}
		return true;
	}

	priority_queue<comparison_element, vector< comparison_element > , element_comparator> sorted_cell_queue;
	priority_queue<oneleft_element, vector< oneleft_element > , oneleft_comparator> oneleft_cell_queue;

	virtual void enqueue_oneleft_element(const CELL_INDEX_TYPE& cellid) {
		oneleft_element to_insert;
		// basic optimization: 1-cells will never have to be "zipped up"
		to_insert.dim = my_mesh_handler->dimension(cellid);	
		if (to_insert.dim == 1) return;

		// compute weight
		int weight = 0;
		int SANITY_COUNT = 0;
		cellIterator it;
		iteratorOperator& facets = my_mesh_handler->facets(cellid, it);
		for (facets.begin(it); facets.valid(it); facets.advance(it)) {
			CELL_INDEX_TYPE temp_id = facets.value(it);
			if (my_grad_field->get_assigned(temp_id)) {
				if (my_grad_field->get_critical(temp_id)) {
					weight += my_mesh_handler->dimension(temp_id);
				} else {
					weight += 
						my_mesh_handler->dimension(my_grad_field->get_pair(temp_id));
				}
			} else {
				SANITY_COUNT++;
			}
		}

		//if (SANITY_COUNT != 1) printf("INSANE IN THE MAINFRAME\n");

		to_insert.cellid = cellid;
		to_insert.boundary = my_mesh_handler->boundary_value(cellid);
		to_insert.dim = my_mesh_handler->dimension(cellid);
		to_insert.value = my_mesh_function->cell_value(cellid);
		to_insert.weight = weight;
		oneleft_cell_queue.push(to_insert);

	}

	CELL_INDEX_TYPE insert_time;
	virtual void enqueue_sorted_element(const CELL_INDEX_TYPE& cellid, CELL_INDEX_TYPE itime, unsigned char mark) {
		comparison_element to_insert;
		to_insert.cellid = cellid;
		to_insert.insertion_time = itime++;
		to_insert.value = my_mesh_function->cell_value(cellid);
		to_insert.boundary = my_mesh_handler->boundary_value(cellid);
		//printf("insert %d -> fval=%.2f\n", cellid, to_insert.value);
		my_grad_field->set_mark(cellid, mark);
		sorted_cell_queue.push(to_insert);
	}

	// seed my queue
	virtual void seedQueueWithDMinima(const DIM_TYPE& dim) {
		int counter = 0;
		cellIterator it;
		iteratorOperator& d_cells = my_mesh_handler->d_cells_iterator(dim, it);
		for(d_cells.begin(it); d_cells.valid(it); d_cells.advance(it)) {
			CELL_INDEX_TYPE cellid = d_cells.value(it);

			if (! my_grad_field->get_assigned(cellid) &&
				less_than_all_unassigned_cofacets(cellid)) {
				// potential critical point, so enqueue
					//printf("potential minimum %d\n", cellid);
				enqueue_sorted_element(cellid, doublecount, 0);
				counter++;
			}
		}
		printf("seeding queue with %d %d-cells\n", counter, dim);
	}

	// we don't need to add d+1 cells, only d+2 cells
	virtual void decrement_cofacets_num_unpaired_facets(const CELL_INDEX_TYPE& cellid, 
		const bool& add_to_oneleft) {
		cellIterator it;
		iteratorOperator& cofacets = my_mesh_handler->cofacets(cellid, it);
		for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
			CELL_INDEX_TYPE temp_cell = cofacets.value(it);
			if (my_grad_field->get_assigned(temp_cell) == false) {
				CELL_INDEX_TYPE num_unpaired = 
					my_grad_field->get_num_unpaired_facets(temp_cell);
				//if (num_unpaired == 0) continue;
				num_unpaired -= 1;
				//if (num_unpaired < 0 || num_unpaired > 7) printf("UNP=%lld\n", num_unpaired);
				my_grad_field->set_num_unpaired_facets(temp_cell, num_unpaired);
				if (num_unpaired == 1 && add_to_oneleft) {
					enqueue_oneleft_element(temp_cell);
				}
			}
		}
	}
	

	virtual void make_critical(const CELL_INDEX_TYPE& cellid, 
		const DIM_TYPE& dim) {
			//printf("making critical %d %d\n", cellid, dim);
		my_grad_field->set_assigned(cellid, true);
		my_grad_field->set_critical(cellid, true);
		my_grad_field->set_dim_asc_man(cellid, my_mesh_handler->max_dim() - dim);

		decrement_cofacets_num_unpaired_facets(cellid, false);
	}


	virtual void pair(const CELL_INDEX_TYPE& tail,
		const CELL_INDEX_TYPE& head, 
		const bool& add_tail_to_oneleft) {
			//printf("p(%d->%d)\n", tail, head);
		my_grad_field->set_assigned(tail, true);
		my_grad_field->set_assigned(head, true);

		DIM_TYPE mindim = my_mesh_handler->max_dim();
		cellIterator it;
		iteratorOperator& facets = my_mesh_handler->facets(head, it);
		for (facets.begin(it); facets.valid(it); facets.advance(it)) {
			CELL_INDEX_TYPE temp_id = facets.value(it);
			if (temp_id == tail) continue;
			DIM_TYPE otherdim = my_grad_field->get_dim_asc_man(temp_id);
			if (otherdim < mindim) mindim = otherdim;
		}

		my_grad_field->set_critical(head, false);
		my_grad_field->set_critical(tail, false);
		my_grad_field->set_dim_asc_man(head, mindim);
		my_grad_field->set_dim_asc_man(tail, mindim);
		my_grad_field->set_pair(head, tail);
		my_grad_field->set_pair(tail, head);

		decrement_cofacets_num_unpaired_facets(tail, add_tail_to_oneleft);
		decrement_cofacets_num_unpaired_facets(head, true);

	}
	////////////////////////
	/////// ADD NEIGHBORS - facets of co-facets that are not assigned and not marked
	////////////

	virtual void add_neighbors_to_sort(const CELL_INDEX_TYPE& cellid) {
		
		cellIterator it;
		iteratorOperator& cofacets = my_mesh_handler->cofacets(cellid, it);
		for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
			CELL_INDEX_TYPE temp_cell = cofacets.value(it);

			cellIterator it2;
			iteratorOperator& facets = my_mesh_handler->facets(temp_cell, it2);
			for (facets.begin(it2); facets.valid(it2); facets.advance(it2)) {
				CELL_INDEX_TYPE temp_neg = facets.value(it2);

				if (temp_neg != cellid &&
					! my_grad_field->get_assigned(temp_neg) &&
					! my_grad_field->get_mark(temp_neg)
					) {
						//printf("enqueuing %d\n", temp_neg);
					enqueue_sorted_element(temp_neg, ++insert_time, 1);
				}
			}
		}
	}



	// return index in candidates of pair
	virtual int pick_from_candidates(const CELL_INDEX_TYPE& cellid, 
		const vector<CELL_INDEX_TYPE>& candidates) {
		int minloc = 0; 
		dtype minval = lowest_facet_value(candidates[minloc]);
		//printf("%d=%.2f\n", candidates[0], minval);
		for (int i = 1; i < candidates.size(); i++) {
			
			dtype otherval = lowest_facet_value(candidates[i]);
			//printf("%d=%.2f\n", candidates[i], otherval);
			if (otherval < minval) {
				minval = otherval;
				minloc = i;
			}
		}
		return minloc;
	}

	virtual void pick_and_pair(const comparison_element& element,
		const DIM_TYPE& dim) {
		// assume it's unassigned
			//if (my_grad_field->get_assigned(element.cellid)) 
			//	printf("WHOA assigned already!!\n");

		vector<CELL_INDEX_TYPE> candidates;
		cellIterator it;
		iteratorOperator& cofacets = my_mesh_handler->cofacets(element.cellid, it);
		for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
			CELL_INDEX_TYPE temp_cell = cofacets.value(it);
			// it's not assigned and is "lower"
			if (my_grad_field->get_assigned(temp_cell) == false &&
				my_grad_field->get_num_unpaired_facets(temp_cell) == 1 &&
				element.value >= my_mesh_function->cell_value(temp_cell) &&
				element.boundary == my_mesh_handler->boundary_value(temp_cell)) {
					candidates.push_back(temp_cell);
			}
		}

		add_neighbors_to_sort(element.cellid);

		// if no candidates, we have a critical point
		if (candidates.size() == 0) {
			make_critical(element.cellid, dim);
			return;
		}

		// so we have candidates
		if (candidates.size() == 1) {
			pair(element.cellid, candidates[0], false);
			return;
		}

		// else we have a choice.
		// for now pick lowest value
		//printf("have coice:\n");
		int minloc = pick_from_candidates(element.cellid, candidates);
		pair(element.cellid, candidates[minloc], false);

	}

	virtual void assignDArrows(const DIM_TYPE& dim) {
		// initialize queue
		seedQueueWithDMinima(dim);

		while (! sorted_cell_queue.empty()) {
			comparison_element top_element = sorted_cell_queue.top();
			sorted_cell_queue.pop();

			//if (! my_grad_field->get_assigned(top_element.cellid) &&
			//	! my_grad_field->get_mark(top_element.cellid)) {
			//		printf("WHOA THERE NELLY\n");
			//}
			//intf("considering %d\n", top_element.cellid);
			if (! my_grad_field->get_assigned(top_element.cellid)) {
				//intf("pairing->%d\n", top_element.cellid);
				pick_and_pair(top_element, dim);
			} 
			//else {
			//	printf("whoa\n");
			//}
		}
	}

	virtual void find_only_and_pair(const oneleft_element& element) {
		
		int SANITY_COUNT = 0;
		cellIterator it;
		iteratorOperator& facets = my_mesh_handler->facets(element.cellid, it);
		for (facets.begin(it); facets.valid(it); facets.advance(it)) {
			CELL_INDEX_TYPE temp_id = facets.value(it);

			if(my_grad_field->get_assigned(temp_id)) continue;
			
			SANITY_COUNT++;
			// now we have the unassigned one
			// TEST FOR DECREASING
			if (my_mesh_function->cell_value(temp_id) == 
				my_mesh_function->cell_value(element.cellid) &&
				element.boundary == my_mesh_handler->boundary_value(temp_id)) 
				pair(temp_id, element.cellid, true);

		}
		//if (SANITY_COUNT != 1) printf("UUUUBER INSAAAANE! %d %d\n", SANITY_COUNT,
		//	my_mesh_handler->dimension(element.cellid));


	}

	virtual void zipUp() {

		while (! oneleft_cell_queue.empty()) {
			oneleft_element top_element = oneleft_cell_queue.top();
			oneleft_cell_queue.pop();

			if(my_grad_field->get_assigned(top_element.cellid)) continue;
			find_only_and_pair(top_element);
		}

	}

public:

	mscSimpleGradientBuilder(
		mscBasicMeshFunction<dtype>* mesh_function,
		mscBasicMeshHandler* mesh_handler,
		mscBasicGradientField* grad_field,
		mscArrayFactory* array_factory) : 
			my_mesh_function(mesh_function),
			my_mesh_handler(mesh_handler), 
			my_grad_field(grad_field),
			my_array_factory(array_factory) {
				insert_time = 0;
				doublecount = mesh_handler->num_cells() *2;

	}

	virtual void computeGradient() {
		init_all();
		for (DIM_TYPE i = 0; i <= my_mesh_handler->max_dim(); i++) {
			assignDArrows(i);
			zipUp();
		}
	}

		virtual dtype lowest_facet_value(const CELL_INDEX_TYPE& cellid) {
		
			if (my_mesh_handler->dimension(cellid) == 0) return my_mesh_function->cell_value(cellid);
			
			bool init = false;
		dtype value;
		cellIterator it2;
		iteratorOperator& facets = my_mesh_handler->facets(cellid, it2);
		for (facets.begin(it2); facets.valid(it2); facets.advance(it2)) {
			CELL_INDEX_TYPE temp_neg = facets.value(it2);
			//dtype nvalue = my_mesh_function->cell_value(temp_neg);
			dtype nvalue = lowest_facet_value(temp_neg);
			if (! init) {
				value = nvalue;	init = true;
			} else if (nvalue < value) {
				value = nvalue;
			}
		}		
		return value;
	}
};




#endif
