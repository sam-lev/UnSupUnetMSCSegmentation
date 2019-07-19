#ifndef MSC_ASSISTED_GRADIENT_BUILDER
#define MSC_ASSISTED_GRADIENT_BUILDER

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscRegularGridImplicitMeshHandler.h"
#include "mscRegularGrid3DMeshFunction.h"
#include "mscRegularGrid3DGradientField.h"
#include "mscArrayFactory.h"

using namespace std;

template<typename dtype>
class mscAssistedGradientBuilder {
protected:
	mscRegularGrid3DMeshFunction<dtype>* my_mesh_function;
	mscBasicMeshHandler* my_mesh_handler;
	mscRegularGrid3DGradientField* my_grad_field;
	mscArrayFactory* my_array_factory;
	const char* my_fname;

	virtual void init_all() {

		// load "assist" information from file
		char tfname[2048];
		sprintf(tfname, "%s.pregrad", my_fname);
		my_grad_field->load_from_file(tfname);
	}


	// throw in simulation of simplicity
	virtual bool f_greater_than(const CELL_INDEX_TYPE& a, const CELL_INDEX_TYPE& b) {
		dtype av = my_mesh_function->cell_value(a);
		dtype bv = my_mesh_function->cell_value(b);
		if (av == bv) return a > b;
		return  av > bv;
	}

	//bool lowest_facet(const CELL_INDEX_TYPE& cellid, 
	//	const CELL_INDEX_TYPE& cofacetid, CELL_INDEX_TYPE& lowest) {
	//	
	//		bool haslowest = false;

	//		cellIterator it2;
	//		iteratorOperator& facets = my_mesh_handler->facets( cofacetid, it2);

	//		for(facets.begin(it2);facets.valid(it2);facets.advance(it2)) {
	//			CELL_INDEX_TYPE facetid = facets.value(it2);
	//			if (facetid != cellid &&
	//				my_grad_field->get_num_unpaired_facets(cofacetid) ==
	//				my_grad_field->get_num_unpaired_facets(facetid) &&
	//				f_greater_than(cellid, facetid)) {
	//					if (haslowest) {
	//						if (f_greater_than(facetid, lowest) {
	//							lowest = facetid;
	//						}
	//					} else {
	//						lowest = facetid;
	//						haslowest = true;
	//					}
	//			}
	//		}
	//		return haslowest;
	//}

	CELL_INDEX_TYPE other_vertex(const CELL_INDEX_TYPE& cellid, 
		const CELL_INDEX_TYPE& edgeid){

			cellIterator it2;
			iteratorOperator& facets = my_mesh_handler->facets( edgeid, it2);

			for(facets.begin(it2);facets.valid(it2);facets.advance(it2)) {
				CELL_INDEX_TYPE facetid = facets.value(it2);
				if (facetid != cellid ) {
					return facetid;
				}
			}
			printf("ERROR this should never happen, no other vertex found\n");
			return 0;
	}

	virtual void do_restricted_vertex(const CELL_INDEX_TYPE& cellid) {

		BOUNDARY_TYPE boundary = my_mesh_handler->boundary_value(cellid);

		bool haslower = false;
		CELL_INDEX_TYPE lowestcofacet;
		CELL_INDEX_TYPE lowestface;
		std::vector<CELL_INDEX_TYPE> lowerstar;

		cellIterator it;
		iteratorOperator& cofacets = my_mesh_handler->cofacets(cellid, it);
		for(cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {

			CELL_INDEX_TYPE cofacetid = cofacets.value(it);

			if (my_grad_field->get_num_unpaired_facets(cofacetid) == 
				my_grad_field->get_num_unpaired_facets(cellid) && 
				boundary == my_mesh_handler->boundary_value(cofacetid)) {
					CELL_INDEX_TYPE otherv = other_vertex(cellid, cofacetid);
					// now test
					if (! f_greater_than(cellid, otherv)) continue;
					lowerstar.push_back(cofacetid);
					if (haslower) {
						if (f_greater_than(lowestface, otherv)) {
							lowestface = otherv;
							lowestcofacet = cofacetid;
						}
					} else {
						lowestface = otherv;
						lowestcofacet = cofacetid;
						haslower = true;
					}
			}
		}
		//

		if (lowerstar.size() == 0) {
			my_grad_field->set_critical(cellid, true);
			my_grad_field->set_assigned(cellid, true);
		} else if (lowerstar.size() == 1) {
			my_grad_field->set_pair(cellid, lowestcofacet);
			my_grad_field->set_assigned(cellid, true);
			my_grad_field->set_pair(lowestcofacet, cellid);
			my_grad_field->set_assigned(lowestcofacet, true);
			// for sanity
			if (lowestcofacet != lowerstar[0]) {
				printf("ERROR, lowestface != lowerstar[0]\n");
			}
		} else {
			my_grad_field->set_pair(cellid, lowestcofacet);
			my_grad_field->set_assigned(cellid, true);
			my_grad_field->set_pair(lowestcofacet, cellid);
			my_grad_field->set_assigned(lowestcofacet, true);
			for (int i = 0; i < lowerstar.size(); i++) {
				if (lowerstar[i] == lowestcofacet) continue;
				my_grad_field->set_critical(lowerstar[i], true);
				my_grad_field->set_assigned(lowerstar[i],true);
			}
		}
	}



	virtual void do_01s() {
		cellIterator it;
		iteratorOperator& verts = my_mesh_handler->d_cells_iterator(0, it);
		for (verts.begin(it); verts.valid(it); verts.advance(it)) {
			do_restricted_vertex(verts.value(it));
		}
	}
	CELL_INDEX_TYPE other_quad(const CELL_INDEX_TYPE& cellid, 
		const CELL_INDEX_TYPE& edgeid){

			cellIterator it2;
			iteratorOperator& cofacets = my_mesh_handler->cofacets( edgeid, it2);

			for(cofacets.begin(it2);cofacets.valid(it2);cofacets.advance(it2)) {
				CELL_INDEX_TYPE cofacetid = cofacets.value(it2);
				if (cofacetid != cellid ) {
					return cofacetid;
				}
			}
			printf("ERROR this should never happen, no other quad found\n");
			return 0;
	}

	virtual void do_restricted_quads(const CELL_INDEX_TYPE& cellid) {

		BOUNDARY_TYPE boundary = my_mesh_handler->boundary_value(cellid);

		bool hashigher = false;
		CELL_INDEX_TYPE highestfacet;
		CELL_INDEX_TYPE highestcofacet;
		std::vector<CELL_INDEX_TYPE> upperstar;

		cellIterator it;
		iteratorOperator& facets = my_mesh_handler->facets(cellid, it);
		for(facets.begin(it); facets.valid(it); facets.advance(it)) {

			CELL_INDEX_TYPE facetid = facets.value(it);

			if (my_grad_field->get_num_unpaired_facets(facetid) == 
				my_grad_field->get_num_unpaired_facets(cellid) && 
				boundary == my_mesh_handler->boundary_value(facetid)) {
					CELL_INDEX_TYPE otherq = other_quad(cellid, facetid);
					//CELL_INDEX_TYPE otherv = other_vertex(cellid, cofacetid);
					// now test
					if (! f_greater_than(otherq, cellid)) continue;
					upperstar.push_back(facetid);
					if (hashigher) {
						if (f_greater_than(otherq, highestfacet)) {
							highestfacet = facetid;
							highestcofacet = otherq;
						}
					} else {
						highestcofacet = otherq;
						highestfacet = facetid;
						hashigher = true;
					}
			}
		}
		//

		if (upperstar.size() == 0) {
			my_grad_field->set_critical(cellid, true);
			my_grad_field->set_assigned(cellid, true);
		} else if (upperstar.size() == 1) {
			my_grad_field->set_pair(cellid, highestfacet);
			my_grad_field->set_assigned(cellid, true);
			my_grad_field->set_pair(highestfacet, cellid);
			my_grad_field->set_assigned(highestfacet, true);
			// for sanity
			if (highestfacet != upperstar[0]) {
				printf("ERROR, highestcoface != lowerstar[0]\n");
			}
		} else {
			my_grad_field->set_pair(cellid, highestfacet);
			my_grad_field->set_assigned(cellid, true);
			my_grad_field->set_pair(highestfacet, cellid);
			my_grad_field->set_assigned(highestfacet, true);
			for (int i = 0; i < upperstar.size(); i++) {
				if (upperstar[i] == highestfacet) continue;
				my_grad_field->set_critical(upperstar[i], true);
				my_grad_field->set_assigned(upperstar[i],true);
			}
		}
	}



	virtual void do_21s() {
		cellIterator it;
		iteratorOperator& quads = my_mesh_handler->d_cells_iterator(2, it);
		for (quads.begin(it); quads.valid(it); quads.advance(it)) {
			do_restricted_quads(quads.value(it));
		}
	}





public:

	mscAssistedGradientBuilder(
		mscRegularGrid3DMeshFunction<dtype>* mesh_function,
		mscRegularGridImplicitMeshHandler* mesh_handler,
		mscRegularGrid3DGradientField* grad_field,
		const char* fname,
		mscArrayFactory* array_factory) : 
	my_mesh_function(mesh_function),
		my_mesh_handler(mesh_handler), 
		my_grad_field(grad_field),
		my_fname(fname),
		my_array_factory(array_factory) {
	}

	virtual void computeGradient() {
		init_all();
		do_01s();
		do_21s();
		cellIterator it;
		iteratorOperator& all = my_mesh_handler->all_cells_iterator(it);
		for (all.begin(it); all.valid(it); all.advance(it)) {
			CELL_INDEX_TYPE cid = all.value(it);
			if (! my_grad_field->get_assigned(cid)) {
				float coords[3];
				my_mesh_handler->centroid(cid, coords);
				printf("ERROR %d not assigned, %d, (%f, %f, %f)\n", cid,
					my_mesh_handler->dimension(cid), coords[0], coords[1], coords[2]);
			}
		}
	}

};




#endif
