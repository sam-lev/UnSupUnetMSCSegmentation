#ifndef MSC_SHORTCUTTING_CONVERGENT_GRADIENT_BUILDER
#define MSC_SHORTCUTTING_CONVERGENT_GRADIENT_BUILDER

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#include "mscBasicGradientField.h"
#include "mscArrayFactory.h"
#include "mscSimpleGradientBuilder.h"
#include "mscConvergentGradientBuilder.h"

#include <vector>
#include <queue>
#include <map>
#ifndef WIN32
#include <cmath>
#endif


using namespace std;
template<typename dtype>
class mscShortcuttingConvergentGradientBuilder : public mscConvergentGradientBuilder<dtype> {

protected:


	virtual bool less_than_all_unassigned_neighbors(const CELL_INDEX_TYPE& cellid) {

		BOUNDARY_TYPE boundary = this->my_mesh_handler->boundary_value(cellid);
		cellIterator it;
		iteratorOperator& neighbors = this->my_mesh_handler->neighbor_vertices(cellid, it);
		for (neighbors.begin(it);neighbors.valid(it); neighbors.advance(it)) {
			CELL_INDEX_TYPE nid = neighbors.value(it);
			if (this->my_grad_field->get_assigned(nid) == 0 && 
				boundary <= this->my_mesh_handler->boundary_value(nid) &&
				! this->my_mesh_function->less_than(cellid, nid)) {
					return false;
			}
		}
		return true;
	}


	// seed my queue
	virtual void seedQueueWithDMinima(const DIM_TYPE& dim) {
		int counter = 0;
		cellIterator it;
		iteratorOperator& d_cells = this->my_mesh_handler->d_cells_iterator(dim, it);
		// start specialized stuff for low dim stuff
		if (dim == 0) {

			for(d_cells.begin(it); d_cells.valid(it); d_cells.advance(it)) {
				CELL_INDEX_TYPE cellid = d_cells.value(it);

				if (! this->my_grad_field->get_assigned(cellid) &&
					this->less_than_all_unassigned_neighbors(cellid)) {
						// potential critical point, so enqueue
						//printf("potential minimum %d\n", cellid);
						this->enqueue_sorted_element(cellid, this->doublecount, 0);
						counter++;

				}
			}

		} else {

			for(d_cells.begin(it); d_cells.valid(it); d_cells.advance(it)) {
				CELL_INDEX_TYPE cellid = d_cells.value(it);

				if (! this->my_grad_field->get_assigned(cellid) &&
					this->less_than_all_unassigned_cofacets(cellid)) {
						// potential critical point, so enqueue
						//printf("potential minimum %d\n", cellid);
						this->enqueue_sorted_element(cellid, this->doublecount, 0);
						counter++;
				}
			}
		}
		printf("seeding queue with %d %d-cells\n", counter, dim);
	}

#ifdef PROBABILITY_CUTOFF_VALUE
	void hack_reduce(MemberDist& md) {

		vector<idfpair> res;
		if (md.pairs.size() < 2) return;
		for (int i = 0; i < md.pairs.size(); i++) {
			if (md.pairs[i].prob > PROBABILITY_CUTOFF_VALUE) res.push_back(md.pairs[i]);
		}
		md.pairs = res;

	}
#endif

	// return index in candidates of pair
	virtual int pick_from_candidatesV(const CELL_INDEX_TYPE& cellid, 
		const vector<CELL_INDEX_TYPE>& candidates,
		const vector<CELL_INDEX_TYPE>& candidatesV,
		const vector<MemberDist>& mdcand,
		int countDepMD) {

			// find weights for computing my probability
			dtype cell_value = this->my_mesh_function->cell_value(cellid);

			int result = 0; 
			vector<float> values;
			float temp_sum = 0;

			for (int i = 0; i < candidates.size(); i++) {
				float diff_val = (float) (cell_value - this->my_mesh_function->cell_value(candidatesV[i]));
				temp_sum += diff_val;
				values.push_back(diff_val);
			}
			if (temp_sum == 0.0f) {
				for (int i=0;i<values.size();i++) values[i] = 1.0f / (float)values.size();
			} else {
				for (int i=0;i<values.size();i++) values[i] = values[i] / temp_sum;
			}

			// compute my local probabilities
			MemberDist my_dist;

			for (int i = 0; i<mdcand.size(); i++) {
				this->mdCombine(values[i], mdcand[i], my_dist);
			}
			// reduce
#ifdef PROBABILITY_CUTOFF_VALUE
			this->hack_reduce(my_dist);
#endif
			// yay, have my distribution now!
			my_dist.count = countDepMD;

#ifdef DEBUG_SIZE_OF_MD
			this->mdcounts[my_dist.pairs.size()]++;
#endif
#ifdef FANCY_PROB_OUTPUT
			// THIS CAN BE REMOVED- ONLY USED FOR FANCY RENDERING
			idfpair p; p.id = cellid;
			p.prob = this->mdMax(my_dist).prob;
			maxvals.push_back(p);
#endif
#ifdef DEBUG_SIZE_OF_MD
			if (mdcounter % 200 == 0) {
				mdsizes.push_back(my_dists.size());
			}
#endif
			// compute my local weights

			int maxloc = 0; 
			vector<float> tress;
			float maxval = this->mdDot3(mdcand[0], my_dist);
			tress.push_back(maxval);
			//printf("%d=%.2f\n", candidates[0], minval);
			for (int i = 1; i < mdcand.size(); i++) {

				float otherval = this->mdDot3(mdcand[i], my_dist);
				//printf("%d=%.2f\n", candidates[i], otherval);
				tress.push_back(otherval);
				if (otherval > maxval) {
					maxval = otherval;
					maxloc = i;
				}
			}

			// check if lowest
			for(int i = 0; i < mdcand.size(); i++) {
				if (i == maxloc) continue;
				if (tress[maxloc] == tress[i]) {
					// pick lower value
					if (this->my_mesh_function->cell_value(candidatesV[i]) <= 
						this->my_mesh_function->cell_value(candidatesV[maxloc])){
							//printf("EHHEHEHE\n");
							maxloc = i;
					}
				}
			}


			// maxloc has my pair
			for (int i =0; i < my_dist.pairs.size(); i++) {
				my_dist.pairs[i].picked = false;
			}
			for (int i =0; i < my_dist.pairs.size(); i++) {
				for (int j = 0; j < mdcand[maxloc].pairs.size(); j++) {
					if (my_dist.pairs[i].id == mdcand[maxloc].pairs[j].id) {
						my_dist.pairs[i].picked |= mdcand[maxloc].pairs[j].picked;
					}
				}
			}

			if (my_dist.count > 0) this->my_dists[cellid] = my_dist;
			return maxloc;			


	}

	virtual bool mdDecrementVertex(const CELL_INDEX_TYPE& cellid, const DIM_TYPE& dim) {

		DIM_TYPE mydim = this->my_mesh_handler->max_dim() - dim;
		if (this->my_dists[cellid].count == 1 && this->my_erase) {
			//printf("erasing %d\n", candidates[i]);
			this->my_dists.erase(cellid);

		} else {
			this->my_dists[cellid].count--;
		}

		return true;
	}

	virtual bool mdCopyNeighborMDAndDec(const CELL_INDEX_TYPE& cellid,
		const DIM_TYPE& dim, MemberDist& result) {
			//if(my_dists.count(cellid) == 0) 
			//	printf("ERROR neighbor not in thingy!!\m");
			MemberDist& neb = this->my_dists[cellid];

			for (int i = 0; i < neb.pairs.size(); i++) {
				//idfpair p; p.id = neb.pairs[i].id;
				//p.picked = neb.pairs[i].picked;
				//p.prob = neb.pairs[i].prob;
				//result.pairs.push_back(p);
				result.pairs.push_back(neb.pairs[i]);
			}
			// now decrement 
			if (this->my_dists[cellid].count == 1 && this->my_erase) {
				//printf("erasing %d\n", candidates[i]);
				this->my_dists.erase(cellid);

			} else {
				this->my_dists[cellid].count--;
			}
			return true;
	}


	virtual void add_neighbors_to_sortV(const CELL_INDEX_TYPE& cellid) {

		cellIterator it;
		iteratorOperator& cofacets = this->my_mesh_handler->cofacets(cellid, it);
		for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
			CELL_INDEX_TYPE temp_cell = cofacets.value(it);

			cellIterator it2;
			iteratorOperator& facets = this->my_mesh_handler->facets(temp_cell, it2);
			for (facets.begin(it2); facets.valid(it2); facets.advance(it2)) {
				CELL_INDEX_TYPE temp_neg = facets.value(it2);

				if (temp_neg != cellid &&
					! this->my_grad_field->get_assigned(temp_neg) &&
					! this->my_grad_field->get_mark(temp_neg)
					) {
						//printf("enqueuing %d\n", temp_neg);
						this->enqueue_sorted_element(temp_neg, ++this->insert_time, 1);
				}
			}
		}
	}


	virtual void pairV(const CELL_INDEX_TYPE& tail,
		const CELL_INDEX_TYPE& head, 
		const bool& add_tail_to_oneleft) {
			//printf("p(%d->%d)\n", tail, head);
			this->my_grad_field->set_assigned(tail, true);
			this->my_grad_field->set_assigned(head, true);

			DIM_TYPE maxdim = this->my_mesh_handler->max_dim();

			this->my_grad_field->set_critical(head, false);
			this->my_grad_field->set_critical(tail, false);
			this->my_grad_field->set_dim_asc_man(head, maxdim);
			this->my_grad_field->set_dim_asc_man(tail, maxdim);
			this->my_grad_field->set_pair(head, tail);
			this->my_grad_field->set_pair(tail, head);

			//-- replace....
	//	decrement_cofacets_num_unpaired_facets(tail, false);
			this->decrement_cofacets_num_unpaired_facets(head, true);

	}

	virtual void make_critical_vertex(const CELL_INDEX_TYPE& cellid, 
		const DIM_TYPE& dim, int countDepMD) {
			//printf("making critical %d %d\n", cellid, dim);
			this->my_grad_field->set_assigned(cellid, true);
			this->my_grad_field->set_critical(cellid, true);
			this->my_grad_field->set_dim_asc_man(cellid, this->my_mesh_handler->max_dim() - dim);
			
			this->decrement_cofacets_num_unpaired_facets(cellid, false);

			idfpair p; p.id = cellid; p.prob = 1.0f; p.picked = true;
#ifdef FANCY_PROB_OUTPUT
			maxvals.push_back(p);
#endif
			//now add the distribution!!
			if (countDepMD > 0) {
				MemberDist md;
				md.count = countDepMD;
				md.mydest = cellid;
				md.pairs.push_back(p);
				this->my_dists[cellid] = md;
			}
	}

	virtual void pick_and_pair(const typename mscSimpleGradientBuilder<dtype>::comparison_element& element,
		const DIM_TYPE& dim) {

			if (dim == 0) {



				vector<CELL_INDEX_TYPE> candidates;
				vector<CELL_INDEX_TYPE> candidatesV;
				vector<MemberDist> mdcand;

#ifdef EXPECTED_NEIGHBORHOOD_SIZE
				candidates.reserve(EXPECTED_NEIGHBORHOOD_SIZE);
				candidatesV.reserve(EXPECTED_NEIGHBORHOOD_SIZE);
				mdcand.reserve(EXPECTED_NEIGHBORHOOD_SIZE);
#endif

				int countDepMD = 0;
				cellIterator it;
				iteratorOperator& neighbors = this->my_mesh_handler->neighbor_vertices(element.cellid, it);
				for (neighbors.begin(it);
					neighbors.valid(it); 
					neighbors.advance(it)) {
						CELL_INDEX_TYPE nid = neighbors.value(it);
						CELL_INDEX_TYPE edgeid = (element.cellid + nid) / 2; // edge id is always average
						// find candidates - it's not assigned and is "lower"
						if (this->my_grad_field->get_assigned(edgeid) == false &&
							this->my_grad_field->get_assigned(nid) == true ) {
								//if (element.value < this->my_mesh_function->cell_value(nid)) {
								//	printf("ERROR: %f < %f\n",
								//		element.value, 
								//		this->my_mesh_function->cell_value(nid));
								//}
								if (element.value >= this->my_mesh_function->cell_value(edgeid) &&
									element.boundary == this->my_mesh_handler->boundary_value(edgeid)) {
									// THEN THIS IS A CANDIDATE
									MemberDist res;
									if (this->mdCopyNeighborMDAndDec(nid, dim, res)) {
										mdcand.push_back(res);
										candidates.push_back(edgeid);
										candidatesV.push_back(nid);
									}
								} else {
									this->mdDecrementVertex(nid, dim);
								}
						}

						// enqueue neighbors
						if (! this->my_grad_field->get_assigned(nid) && ! this->my_grad_field->get_mark(nid)) {
							//printf("enqueuing %d\n", temp_neg);
							this->enqueue_sorted_element(nid, ++this->insert_time, 1);
						}

						// count number of times this membership distribution will be looked up
						if (! this->my_grad_field->get_assigned(nid)) countDepMD++;

				}

				//printf("%d cd %d\n", element.cellid, countDepMD);
				// if no candidates, we have a critical point
				if (candidates.size() == 0) {
					
					this->make_critical_vertex(element.cellid, dim, countDepMD);
					return;
				}

				// so we have candidates.
				// for now pick lowest value
				//printf("have coice:\n");
				int minloc = this->pick_from_candidatesV(element.cellid, candidates, candidatesV, mdcand, countDepMD);

				this->pairV(element.cellid, candidates[minloc], false);

			} else {
				this->mscConvergentGradientBuilder<dtype>::pick_and_pair(element, dim);
			}


	}

public:

	void set_eraser(bool val) { this->my_erase = val; }

	mscShortcuttingConvergentGradientBuilder(
		mscBasicMeshFunction<dtype>* mesh_function,
		mscBasicMeshHandler* mesh_handler,
		mscBasicGradientField* grad_field,
		mscArrayFactory* array_factory) : 
	mscConvergentGradientBuilder<dtype>(mesh_function, mesh_handler, grad_field, array_factory) {
		this->my_erase = true;

	}



#ifdef FANCY_PROB_OUTPUT
	vector<idfpair>& getmaxvals() {
		return maxvals;
	}
#endif





	virtual idfpair mdMax(MemberDist& a) {
		int size = a.pairs.size();
		if (size == 0) printf("ERRROROROROROR: mdMax has no elements\n");
		int tempi = 0;
		float sum = a.pairs[tempi].prob;
		for (int i = 1; i < size; i++) {
			sum += a.pairs[i].prob;
			if (a.pairs[i].prob > a.pairs[tempi].prob) tempi = i;
		}
		if (sum < .99f || sum > 1.001f) printf("WHOA Sum of probs = %f\n", sum);
		return a.pairs[tempi];
	}

	// res += scale*a
	virtual void mdCombine(float scale, const MemberDist& a, MemberDist& res) {
		for (int i = 0; i < a.pairs.size(); i++) {
			bool has = false;
			for (int j = 0; j < res.pairs.size(); j++) {
				if (a.pairs[i].id == res.pairs[j].id) {
					has = true;
					res.pairs[j].prob += a.pairs[i].prob * scale;
					res.pairs[j].picked = res.pairs[j].picked || a.pairs[i].picked;
					break;
				}
			}
			if (! has) {
				idfpair p;
				p.id = a.pairs[i].id;
				p.prob = scale * a.pairs[i].prob;
				p.picked = a.pairs[i].picked;
				res.pairs.push_back(p);
			}
		}
	}

};




#endif
