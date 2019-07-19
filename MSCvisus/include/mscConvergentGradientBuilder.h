#ifndef MSC_CONVERGENT_GRADIENT_BUILDER
#define MSC_CONVERGENT_GRADIENT_BUILDER

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#include "mscBasicGradientField.h"
#include "mscArrayFactory.h"
#include "mscSimpleGradientBuilder.h"

#include <vector>
#include <queue>
#include <map>
#ifndef WIN32
#include <cmath>
#endif

using namespace std;
struct idfpair {
	CELL_INDEX_TYPE id;
	float prob;
	bool picked;
};
	struct MemberDist {
		INT_TYPE count;
		CELL_INDEX_TYPE mydest;
		vector<idfpair> pairs;
	};
template<typename dtype>
class mscConvergentGradientBuilder : public mscSimpleGradientBuilder<dtype> {
public:
	//struct MemberDist;
protected:


	bool my_erase;

#ifdef FANCY_PROB_OUTPUT
	vector<idfpair> maxvals;
#endif

	virtual bool mdDecrement(const CELL_INDEX_TYPE& cellid, const DIM_TYPE& dim) {

		DIM_TYPE mydim = this->my_mesh_handler->max_dim() - dim;

		cellIterator it;
		iteratorOperator& facets = this->my_mesh_handler->facets(cellid, it);
		for (facets.begin(it); facets.valid(it); facets.advance(it)) {
			CELL_INDEX_TYPE temp_cell = facets.value(it);
			// it's not assigned and is "lower"
			if (this->my_grad_field->get_assigned(temp_cell) == true &&
				this->my_grad_field->get_dim_asc_man(temp_cell) == mydim) {
					// now decrement 
					if (my_dists[temp_cell].count == 1 && my_erase) {
						//printf("erasing %d\n", candidates[i]);
						my_dists.erase(temp_cell);

					} else {
						my_dists[temp_cell].count--;
					}
			}
		}

		return true;
	}

	virtual bool mdCombineAllFacetsAndDecrement(const CELL_INDEX_TYPE& cellid,
		const DIM_TYPE& dim, MemberDist& result) {

			vector<CELL_INDEX_TYPE> candidates;
			DIM_TYPE mydim = this->my_mesh_handler->max_dim() - dim;

			cellIterator it;
			iteratorOperator& facets = this->my_mesh_handler->facets(cellid, it);
			for (facets.begin(it); facets.valid(it); facets.advance(it)) {
				CELL_INDEX_TYPE temp_cell = facets.value(it);
				// it's not assigned and is "lower"
				if (this->my_grad_field->get_assigned(temp_cell) == true &&
					this->my_grad_field->get_dim_asc_man(temp_cell) == mydim) {
						candidates.push_back(temp_cell);
				}
			}

			if (candidates.size() == 0) {
				//printf("ERROR: mdCombineAllFacets candidates.size == 0, %d, %d\n", dim, mydim);
				return false;
			}

			////if (candidates.size() == 1) {
			////	if (my_dists.count(candidates[0]) == 0) {
			////		printf("ERROR: mdCombineAllFacets, my_dists does not have candidate\n");
			////	}
			////	return my_dists[candidates[0]];
			////}

			float oneoversize = 1.0f / (float) candidates.size();
			for (int i = 0; i < candidates.size(); i++) {
				if (my_dists.count(candidates[i]) == 0) {
					printf("ERROR: mdCombineAllFacets, my_dists does not have candidate %d\n", i);
				}
				mdCombine(oneoversize, my_dists[candidates[i]], result);

				// now decrement 
				if (my_dists[candidates[i]].count == 1 && my_erase) {
					//printf("erasing %d\n", candidates[i]);
					my_dists.erase(candidates[i]);

				} else {
					my_dists[candidates[i]].count--;
				}

			}
			return true;
	}




	map<CELL_INDEX_TYPE, MemberDist> my_dists;





	// set all cells to unassigned
	virtual void init_assigned() {
		cellIterator it;
		iteratorOperator& all_cells = this->my_mesh_handler->all_cells_iterator(it);
		for (all_cells.begin(it); all_cells.valid(it); all_cells.advance(it)) {
			CELL_INDEX_TYPE cellid = all_cells.value(it);
			this->my_grad_field->set_assigned(cellid, 0);
			this->my_grad_field->set_mark(cellid, 0);
		}
	}

	virtual void init_all() {
		init_assigned();
		this->init_number_unpaired_facets();
	}









	// seed my queue
	virtual void seedQueueWithDMinima(const DIM_TYPE& dim) {
		int counter = 0;
		cellIterator it;
		iteratorOperator& d_cells = this->my_mesh_handler->d_cells_iterator(dim, it);
		for(d_cells.begin(it); d_cells.valid(it); d_cells.advance(it)) {
			CELL_INDEX_TYPE cellid = d_cells.value(it);

			if (! this->my_grad_field->get_assigned(cellid) &&
				this->less_than_all_unassigned_cofacets(cellid)) {
					// potential critical point, so enqueue
					//if (dim ==0)printf("potential minimum %d\n", cellid);
					this->enqueue_sorted_element(cellid, this->doublecount, 0);
					counter++;
			}
		}
		printf("seeding queue with %d %d-cells\n", counter, dim);
	}

	////// we don't need to add d+1 cells, only d+2 cells
	////virtual void decrement_cofacets_num_unpaired_facets(const CELL_INDEX_TYPE& cellid, 
	////	const bool& add_to_oneleft, int& count) {
	////		count = 0;
	////	cellIterator it;
	////	iteratorOperator& cofacets = my_mesh_handler->cofacets(cellid, it);
	////	for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
	////		CELL_INDEX_TYPE temp_cell = cofacets.value(it);
	////		if (this->my_grad_field->get_assigned(temp_cell) == false) {
	////			count++;
	////			CELL_INDEX_TYPE num_unpaired = 
	////				this->my_grad_field->get_num_unpaired_facets(temp_cell);
	////			//if (num_unpaired == 0) continue;
	////			num_unpaired -= 1;
	////			//if (num_unpaired < 0 || num_unpaired > 7) printf("UNP=%lld\n", num_unpaired);
	////			this->my_grad_field->set_num_unpaired_facets(temp_cell, num_unpaired);
	////			if (num_unpaired == 1 && add_to_oneleft) {
	////				enqueue_oneleft_element(temp_cell);
	////			}
	////		}
	////	}
	////}

	// number of cofacets that will look up this distribution via a call to
	// mdCombineAllFacetsAndDecrement
	int count_unpaired_facets(const CELL_INDEX_TYPE& cellid) {
		int counter= 0;
		cellIterator it;
		iteratorOperator& facets = this->my_mesh_handler->facets(cellid, it);
		for (facets.begin(it); facets.valid(it); facets.advance(it)) {
			CELL_INDEX_TYPE temp_cell = facets.value(it);
			// it's not assigned and is "lower"
			if (this->my_grad_field->get_assigned(temp_cell) == false) {
					counter++;
			}
		}

		return counter;
	}

	virtual int countMDCofacets(const CELL_INDEX_TYPE& cellid) {
		int count = 0;
		cellIterator it;
		iteratorOperator& cofacets = this->my_mesh_handler->cofacets(cellid, it);
		for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
			CELL_INDEX_TYPE temp_cell = cofacets.value(it);
			//if (this->my_grad_field->get_assigned(temp_cell) == false &&
			//	count_unpaired_facets(temp_cell) > 1) {
			//		count++;
			//}
						if (this->my_grad_field->get_assigned(temp_cell) == false &&
				this->my_grad_field->get_num_unpaired_facets(temp_cell) > 1 /*&&
				this->my_mesh_function->cell_value(cellid) >= this->my_mesh_function->cell_value(temp_cell)*/) {
					count++;
			}
		}
		return count;
	}

	virtual void make_critical(const CELL_INDEX_TYPE& cellid, 
		const DIM_TYPE& dim) {
			//printf("making critical %d %d\n", cellid, dim);
			this->my_grad_field->set_assigned(cellid, true);
			this->my_grad_field->set_critical(cellid, true);
			this->my_grad_field->set_dim_asc_man(cellid, this->my_mesh_handler->max_dim() - dim);

			int count = countMDCofacets(cellid);

			this->decrement_cofacets_num_unpaired_facets(cellid, false);

			idfpair p; p.id = cellid; p.prob = 1.0f; p.picked = true;
#ifdef FANCY_PROB_OUTPUT
			maxvals.push_back(p);
#endif
			//now add the distribution!!
			if (count > 0) {
				MemberDist md;
				md.count = count;
				md.mydest = cellid;
				md.pairs.push_back(p);
				my_dists[cellid] = md;
			}
	}


	// return index in candidates of pair
	virtual int pick_from_candidates(const CELL_INDEX_TYPE& cellid, 
		const vector<CELL_INDEX_TYPE>& candidates,
		const vector<MemberDist>& mdcand) {

			// find weights for computing my probability
			dtype cell_value = this->my_mesh_function->cell_value(cellid);

			int result = 0; 
			vector<float> values;
			float temp_sum = 0;

			for (int i = 0; i < candidates.size(); i++) {
				float diff_val = (float) (cell_value - this->lowest_facet_value(candidates[i]));
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
			// yay, have my distribution now!
			my_dist.count = countMDCofacets(cellid);

#ifdef DEBUG_SIZE_OF_MD
			this->mdcounts[my_dist.pairs.size()]++;
#endif

#ifdef FANCY_PROB_OUTPUT
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
			vector<float> tress; // the similarity of candidate[i]
			float maxval = mdDot3(mdcand[0], my_dist);
			tress.push_back(maxval);
			//printf("%d=%.2f\n", candidates[0], minval);
			for (int i = 1; i < mdcand.size(); i++) {

				float otherval = mdDot3(mdcand[i], my_dist);
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
					if (this->lowest_facet_value(candidates[i]) <= this->lowest_facet_value(candidates[maxloc])){
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

	virtual void pick_and_pair(const typename mscSimpleGradientBuilder<dtype>::comparison_element& element,
		const DIM_TYPE& dim) {
			// assume it's unassigned
			//if (this->my_grad_field->get_assigned(element.cellid)) 
			//	printf("WHOA assigned already!!\n");

			vector<CELL_INDEX_TYPE> candidates;
			vector<MemberDist> mdcand;
			cellIterator it;
			iteratorOperator& cofacets = this->my_mesh_handler->cofacets(element.cellid, it);
			for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
				CELL_INDEX_TYPE temp_cell = cofacets.value(it);
				// it's not assigned and is "lower"
				if (this->my_grad_field->get_assigned(temp_cell) == false &&
					this->my_grad_field->get_num_unpaired_facets(temp_cell) == 1 ) {
						if (element.value >= this->my_mesh_function->cell_value(temp_cell) && 
							element.boundary == this->my_mesh_handler->boundary_value(temp_cell)) {
							MemberDist res;
							if (mdCombineAllFacetsAndDecrement(temp_cell, dim, res)) {
								mdcand.push_back(res);
								candidates.push_back(temp_cell);
							}
						} else {
							mdDecrement(temp_cell, dim);
						}
				}
			}

			this->add_neighbors_to_sort(element.cellid);

			// if no candidates, we have a critical point
			if (candidates.size() == 0) {
				make_critical(element.cellid, dim);
				return;
			}

			// so we have candidates.
			// for now pick lowest value
			//printf("have coice:\n");
			int minloc = pick_from_candidates(element.cellid, candidates, mdcand);

			this->pair(element.cellid, candidates[minloc], false);

	}

#ifdef DEBUG_SIZE_OF_MD
	vector<int> mdcounts;
	vector<int> mdsizes;
	int mdcounter;
#endif

public:

	void set_eraser(bool val) { my_erase = val; }

	mscConvergentGradientBuilder(
		mscBasicMeshFunction<dtype>* mesh_function,
		mscBasicMeshHandler* mesh_handler,
		mscBasicGradientField* grad_field,
		mscArrayFactory* array_factory) : 
	mscSimpleGradientBuilder<dtype>(mesh_function, mesh_handler, grad_field, array_factory) {
		my_erase = true;


#ifdef DEBUG_SIZE_OF_MD
				mdcounts.resize(1024, 0);
#endif
	}

	void check_ldir_vs_uafacetcount(int stage0, int stage, int stage2) {
		cellIterator it;
		iteratorOperator& all = this->my_mesh_handler->all_cells_iterator(it);
		for (all.begin(it); all.valid(it); all.advance(it)) {
			CELL_INDEX_TYPE cellid = all.value(it);
			if (! this->my_grad_field->get_assigned(cellid)) {
				// check values
				int actual = 0;
				cellIterator fit;
				iteratorOperator& facets = this->my_mesh_handler->facets(cellid, fit);
				for (facets.begin(fit); facets.valid(fit); facets.advance(fit)) {
					if (! this->my_grad_field->get_assigned(facets.value(fit))) actual++;
				}

				int recorded = this->my_grad_field->get_num_unpaired_facets(cellid);

				if (recorded != actual) 
					printf("WHOA-%d-%d-%d, %d's recorded = %d, actual %d num unpaired facets dim = %d\n", 
							stage0, stage, stage2, (int) cellid, recorded,actual, 
							this->my_mesh_handler->dimension(cellid));
			}
		}
	}



#ifdef DEBUG_SIZE_OF_MD
	virtual void computeGradient(int stage=0) {
		mdcounter = 0;
#else
	virtual void computeGradient() {
#endif
		init_all();
		for (DIM_TYPE i = 0; i <= this->my_mesh_handler->max_dim(); i++) {
			this->my_dists.clear();
			this->assignDArrows(i);
			printf("Have %d remaining dists\n", my_dists.size());

#ifdef DEBUG_SIZE_OF_MD
			//check_ldir_vs_uafacetcount(stage,i,0);


			char fname[1024];
			sprintf(fname, "mdcount_%d_%d.txt", stage, i);
			FILE* fout = fopen(fname, "w");
			for (int i = 0; i < mdcounts.size(); i++) {
				fprintf(fout, "%d %d\n", i, mdcounts[i]);
				mdcounts[i] = 0;
			}
			fclose(fout);
			sprintf(fname, "mdsizes_%d_%d.txt", stage, i);
			fout = fopen(fname, "w");
			for (int i = 0; i < mdsizes.size(); i++) {
				fprintf(fout, "%d %d\n", i, mdsizes[i]);
			}
			fclose(fout);
			mdsizes.clear();
			mdcounter = 0;
#endif
			this->my_dists.clear();
			this->zipUp();
#ifdef DEBUG_SIZE_OF_MD
			//check_ldir_vs_uafacetcount(stage, i, 1);
#endif
		}
	}

#ifdef FANCY_PROB_OUTPUT
	vector<idfpair>& getmaxvals() {
		return maxvals;
	}
#endif





	virtual float mdDot(const MemberDist& a, const MemberDist& b) {
		float result = 0.0f;
		for (int i = 0; i < a.pairs.size(); i++) {
			for (int j = 0; j < b.pairs.size(); j++) {
				if (a.pairs[i].id == b.pairs[j].id) {
					result += a.pairs[i].prob * b.pairs[j].prob;
					//break;
				}
			}
		}
		return result;
	}

	virtual float mdDot3(const MemberDist& a, const MemberDist& b) {
			float result = 0.0f;
		for (int i = 0; i < a.pairs.size(); i++) {
			for (int j = 0; j < b.pairs.size(); j++) {
				if (a.pairs[i].id == b.pairs[j].id) {
					result += a.pairs[i].prob * b.pairs[j].prob;
					//break;
				}
			}
		}
		return result;
		
		
		//float result = 0.0f;
		//for (int i = 0; i < a.pairs.size(); i++) {
		//	for (int j = 0; j < b.pairs.size(); j++) {
		//		if (a.pairs[i].id == b.pairs[j].id && a.pairs[i].picked) {
		//			result += /*a.pairs[i].prob **/ b.pairs[j].prob;
		//			//break;
		//		}
		//	}
		//}
		//return result;
	}

	virtual float mdDot2(const MemberDist& a, const MemberDist& b) {
		int tempi = 0;
		float mv = 0;
		//float sum = a.pairs[tempi].prob;
		for (int i = 0; i < a.pairs.size(); i++) {
			for (int j = 0; j < b.pairs.size(); j++) {
				if (a.pairs[i].id == b.pairs[j].id) {
					float v = a.pairs[i].prob * b.pairs[j].prob;
					if (v > mv) { tempi = i; mv = v; }
					//sum += a.pairs[i].prob;
					//if (a.pairs[i].prob > a.pairs[tempi].prob) tempi = i;
				}
			}
		}
		//float result = a.pairs[tempi].prob;
		return mv; //result;
	}
	virtual float mdDist(const MemberDist& a, const MemberDist& b) {
		int tempi = 0;
		float mv = 0;
		//float sum = a.pairs[tempi].prob;
		for (int i = 0; i < a.pairs.size(); i++) {
			for (int j = 0; j < b.pairs.size(); j++) {
				if (a.pairs[i].id == b.pairs[j].id) {
					mv += (a.pairs[i].prob - b.pairs[j].prob)*(a.pairs[i].prob - b.pairs[j].prob);
					//float v = a.pairs[i].prob * b.pairs[j].prob;
					//if (v > mv) { tempi = i; mv = v; }
					//sum += a.pairs[i].prob;
					//if (a.pairs[i].prob > a.pairs[tempi].prob) tempi = i;
				}
			}
		}
		//float result = a.pairs[tempi].prob;
		return sqrt(mv); //result;
	}
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
