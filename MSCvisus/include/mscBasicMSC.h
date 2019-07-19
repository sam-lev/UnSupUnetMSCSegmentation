#ifndef MSC_BASIC_MSC
#define MSC_BASIC_MSC


#include <stdio.h>
#include <vector>
#include <map>
#include <queue>
#include <set>

#include "mscBasicGradientField.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#define MAX_CELLID_VALUE 999999999

using namespace std;


template <typename dtype>
class node;
template <typename dtype>
class arc;

template<typename dtype>
class manifold {
public:
	CELL_INDEX_TYPE created;
	CELL_INDEX_TYPE destroyed;
	node<dtype>* basenode;

	manifold<dtype>* merge[2];

	manifold(node<dtype>* n) : created(0), destroyed(MAX_CELLID_VALUE), basenode(n) {
		merge[0] = NULL; merge[1] = NULL;
	}

	manifold(manifold<dtype>* base, manifold<dtype>* other, CELL_INDEX_TYPE created) :
		created(created), destroyed(MAX_CELLID_VALUE), basenode(base->basenode) {
			merge[0] = base; merge[1] = other;
	}
};

template <typename dtype>
class node {
public:
	arc<dtype>* firstarc;
	dtype value;
	DIM_TYPE index;
	CELL_INDEX_TYPE cellid;
	BOUNDARY_TYPE boundary;
	CELL_INDEX_TYPE destroyed;
	manifold<dtype>* dsc_man;
	manifold<dtype>* asc_man;

	node(CELL_INDEX_TYPE id, DIM_TYPE index, dtype value, BOUNDARY_TYPE boundary) :
		value(value), cellid(id), index(index), firstarc(NULL), boundary(boundary),
		destroyed(MAX_CELLID_VALUE) {
			asc_man = new manifold<dtype>(this);
			dsc_man = new manifold<dtype>(this);
	}
};



template <typename dtype>
class arc {
public:

	node<dtype>* lower;
	node<dtype>* upper;
	arc<dtype>* lower_next;
	arc<dtype>* upper_next;

	dtype persistence;

	CELL_INDEX_TYPE created;
	CELL_INDEX_TYPE destroyed;

	arc<dtype>* merge[3];
	vector<CELL_INDEX_TYPE> base_geometry;

	int cancellation_weight;

	arc(node<dtype>* upper, node<dtype>* lower, 
		vector<CELL_INDEX_TYPE> geometry) : lower(lower), upper(upper),
		lower_next(NULL), upper_next(NULL), created(0), destroyed(MAX_CELLID_VALUE) {
			base_geometry.insert(base_geometry.begin(), geometry.begin(), geometry.end());
			persistence = upper->value - lower->value;
			lower_next = lower->firstarc;
			lower->firstarc = this;
			upper_next = upper->firstarc;
			upper->firstarc = this;
			cancellation_weight = 0;
	}

	arc(arc<dtype>* l_upper, arc<dtype>* a, arc<dtype>* u_lower, CELL_INDEX_TYPE created) :
		created(created), destroyed(MAX_CELLID_VALUE) {
			//printf("creating %d-%d->%d-%d->%d-%d\n",
			//	l_upper->upper, l_upper->lower, a->lower, a->upper, 
			//	u_lower->upper, u_lower->lower);

			merge[0] = l_upper; merge[1] = a; merge[2] = u_lower;
			upper = l_upper->upper;
			lower = u_lower->lower;
			lower_next = lower->firstarc;
			upper_next = upper->firstarc;
			lower->firstarc = this;
			upper->firstarc = this;
			cancellation_weight = 0;

			persistence = upper->value - lower->value;
	}


};



template<typename dtype>
class BasicMSC {
protected:
	struct arcCompare {
		bool operator()(const arc<dtype>* left, const arc<dtype>* right) {
			if (left->persistence < right->persistence) return false;
			if (left->persistence > right->persistence ) return true;

			if (left->cancellation_weight < right->cancellation_weight) return false;
			if (left->cancellation_weight > right->cancellation_weight) return true;

			return left < right;
		}
	};
	priority_queue<arc<dtype>*, vector< arc<dtype>* >, arcCompare > edges_to_cancel;


	float m_temp_perc_pers;
	CELL_INDEX_TYPE m_sel_persistence;
	CELL_INDEX_TYPE m_num_destroyed;
	mscBasicGradientField* m_grad;
	mscBasicMeshHandler* m_mesh;
	mscBasicMeshFunction<dtype>* m_func;

	void rec_tdcr(const CELL_INDEX_TYPE& cellid, DIM_TYPE& temp_dim, vector<CELL_INDEX_TYPE>& res) {
		res.push_back(cellid);
		CELL_INDEX_TYPE current = cellid;
		cellIterator it;
		iteratorOperator& facets = m_mesh->facets(current, it);
		for (facets.begin(it); facets.valid(it); facets.advance(it)) {
			CELL_INDEX_TYPE temp_id = facets.value(it);
			if (m_grad->get_critical(temp_id) ) { 
				res.push_back(temp_id);
				//printf("adding arc: %llu\n", temp_id);
				add_arc(nodes[res[0]], nodes[temp_id], res);
				res.pop_back();
			} else /*if (m_grad->get_dim_asc_man(temp_id) == temp_dim)*/ {
				CELL_INDEX_TYPE pair = m_grad->get_pair(temp_id);
				if(pair != cellid && m_mesh->dimension(pair) == m_mesh->dimension(cellid)) {
					res.push_back(temp_id);
					//result.push_back(pair);
					rec_tdcr(pair, temp_dim, res);
					res.pop_back();
				}
			}		

		}
		res.pop_back();
	}

	void set_cancellation_weight(arc<dtype>* a) {
		int res1 = 0;
		arc<dtype>* t = a->lower->firstarc;
		while(t != NULL) {
			//if (isAlive(t) && a->lower == t->lower) res++;
			if (a->lower == t->lower) {
				if (isAlive(t)) res1++;
				t = t->lower_next;
			} else {
				t = t->upper_next;
			}
		} 
		int res2 = 0;
		t = a->upper->firstarc;
		while(t != NULL) {
			//if (isAlive(t) && a->lower == t->lower) res++;
			if (a->lower == t->lower) {
				t = t->lower_next;
			} else {
				if (isAlive(t)) res2++;
				t = t->upper_next;
			}
		} 

		a->cancellation_weight = res1 * res2;
	}


	void rec_man_trace_down(CELL_INDEX_TYPE& cellid, set<CELL_INDEX_TYPE>& res) {
		res.insert(cellid);
		CELL_INDEX_TYPE current = cellid;
		DIM_TYPE cdim = m_mesh->dimension(cellid);
		cellIterator it;
		iteratorOperator& facets = m_mesh->facets(current, it);
		for (facets.begin(it); facets.valid(it); facets.advance(it)) {
			CELL_INDEX_TYPE temp_id = facets.value(it);
			if (m_grad->get_critical(temp_id) ) continue; 

			CELL_INDEX_TYPE temp_pair = m_grad->get_pair(temp_id);

			if (temp_pair == cellid) continue;

			if (m_mesh->dimension(temp_pair) != cdim) continue;

			rec_man_trace_down(temp_pair, res);
		}
	}

	void rec_man_trace_down_restricted(CELL_INDEX_TYPE& cellid, set<CELL_INDEX_TYPE>& res, set <CELL_INDEX_TYPE>& ids) {
		res.insert(cellid);
		CELL_INDEX_TYPE current = cellid;
		DIM_TYPE cdim = m_mesh->dimension(cellid);
		cellIterator it;
		iteratorOperator& facets = m_mesh->facets(current, it);
		for (facets.begin(it); facets.valid(it); facets.advance(it)) {
			CELL_INDEX_TYPE temp_id = facets.value(it);
			if (m_grad->get_critical(temp_id) ) continue; 

			CELL_INDEX_TYPE temp_pair = m_grad->get_pair(temp_id);

			if (temp_pair == cellid) continue;

			if (m_mesh->dimension(temp_pair) != cdim) continue;

			if (ids.count(temp_pair) == 0) continue;

			rec_man_trace_down(temp_pair, res);
		}
	}

	// this function adds EVERY cell in the ascending manifold of cellid - including boundary and 
	// inbetween stuff. quite literally anything that can go boundary->vector repeated to the 
	// critical point.
	void rec_man_trace_up_fill(CELL_INDEX_TYPE& cellid, set<CELL_INDEX_TYPE>& res) {
		if (res.count(cellid) != 0) return;
		
		res.insert(cellid);
		CELL_INDEX_TYPE current = cellid;
		DIM_TYPE cdim = m_mesh->dimension(cellid);
		cellIterator it;
		iteratorOperator& cofacets = m_mesh->cofacets(current, it);
		for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
			CELL_INDEX_TYPE temp_id = cofacets.value(it);
			rec_man_trace_up_fill(temp_id, res);
		}
		if  (m_grad->get_critical(cellid)) return;
		CELL_INDEX_TYPE nid = m_grad->get_pair(cellid);
		if (m_mesh->dimension(nid) != cdim - 1) return;
		rec_man_trace_up_fill(nid, res);
	}

	void rec_res_man_trace_down_fill(CELL_INDEX_TYPE& cellid, set<CELL_INDEX_TYPE>& ids, set<CELL_INDEX_TYPE>& res) {
		if (res.count(cellid) != 0) return;
		if (ids.count(cellid) == 0) return;

		res.insert(cellid);
		CELL_INDEX_TYPE current = cellid;
		DIM_TYPE cdim = m_mesh->dimension(cellid);
		cellIterator it;
		iteratorOperator& facets = m_mesh->facets(current, it);
		for (facets.begin(it); facets.valid(it); facets.advance(it)) {
			CELL_INDEX_TYPE temp_id = facets.value(it);
			rec_res_man_trace_down_fill(temp_id, ids, res);
		}
		if  (m_grad->get_critical(cellid)) return;
		CELL_INDEX_TYPE nid = m_grad->get_pair(cellid);
		if (m_mesh->dimension(nid) != cdim + 1) return;
		rec_res_man_trace_down_fill(nid, ids, res);
	}


	void rec_man_trace_up(CELL_INDEX_TYPE& cellid, set<CELL_INDEX_TYPE>& res) {
		res.insert(cellid);
		CELL_INDEX_TYPE current = cellid;
		DIM_TYPE cdim = m_mesh->dimension(cellid);
		cellIterator it;
		iteratorOperator& cofacets = m_mesh->cofacets(current, it);
		for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
			CELL_INDEX_TYPE temp_id = cofacets.value(it);
			if (m_grad->get_critical(temp_id) ) continue; 

			CELL_INDEX_TYPE temp_pair = m_grad->get_pair(temp_id);

			if (temp_pair == cellid) continue;

			if (m_mesh->dimension(temp_pair) != cdim) continue;

			rec_man_trace_up(temp_pair, res);
		}
	}
	void rec_man_trace_up_restricted(CELL_INDEX_TYPE& cellid, set<CELL_INDEX_TYPE>& res, set <CELL_INDEX_TYPE>& ids) {
		res.insert(cellid);
		CELL_INDEX_TYPE current = cellid;
		DIM_TYPE cdim = m_mesh->dimension(cellid);
		cellIterator it;
		iteratorOperator& cofacets = m_mesh->cofacets(current, it);
		for (cofacets.begin(it); cofacets.valid(it); cofacets.advance(it)) {
			CELL_INDEX_TYPE temp_id = cofacets.value(it);
			if (m_grad->get_critical(temp_id) ) continue; 

			CELL_INDEX_TYPE temp_pair = m_grad->get_pair(temp_id);

			if (temp_pair == cellid) continue;

			if (m_mesh->dimension(temp_pair) != cdim) continue;

			if (ids.count(temp_pair) == 0) continue;

			rec_man_trace_up(temp_pair, res);
		}
	}

	void trace_down_cells_restricted(const CELL_INDEX_TYPE& cellid) {
		//printf("tracing down from %llu\n", cellid);
		vector<CELL_INDEX_TYPE> result;

		DIM_TYPE temp_dim = m_grad->get_dim_asc_man(cellid) + 1;
		result.clear();
		rec_tdcr(cellid, temp_dim, result);
	}



	int countMultiplicity(arc<dtype>* a) {
		int res = 0;
		arc<dtype>* t = a->lower->firstarc;
		while(t != NULL) {
			if (isAlive(t) && a->lower == t->lower &&
				a->upper == t->upper) res++;
			if (a->lower == t->lower) {
				t = t->lower_next;
			} else {
				t = t->upper_next;
			}
		} 
		return res;
	}

	bool validToCancel(arc<dtype>* a) {
		//printf("valid...\n");
		if (! isAlive(a)) { 
			//printf("done\n"); 
			return false; 
		}

		if (a->lower->boundary != a->upper->boundary) { 
			//printf("done\n"); 
			return false; 
		}
		if (countMultiplicity(a) != 1)  { 
			//printf("done\n"); 
			return false; 
		}
		//if (a->upper->value > 118 && a->lower->value <= 118) return false;

		int old_cancellation_weight = a->cancellation_weight;
		this->set_cancellation_weight(a);
		if (a->cancellation_weight > 10000)  {
			//printf("done\n"); 
			return false; }
		if (old_cancellation_weight < a->cancellation_weight) {
			edges_to_cancel.push(a);
			{ 
				//printf("done\n"); 
				return false; }
		}
		//printf("done1\n");
		return true;
	}


	//float persistenceFunction(arc<dtype>* a) {
	//	return a->upper->value - a->lower->value;

	//}

	arc<dtype>* nextArc(node<dtype>* n, arc<dtype>* a) {
		if (a->lower == n) return a->lower_next;
		return a->upper_next;
	}

	virtual arc<dtype>* NewArc( arc<dtype>* l_upper, arc<dtype>* a, arc<dtype>* u_lower, CELL_INDEX_TYPE created) {

		return new arc<dtype>(l_upper, a, u_lower, created);

	}

	void cancel(arc<dtype>* a) {

		//printf("cancelling %d-%d, %f\n", a->lower->index, a->upper->index, a->persistence);

		node<dtype>* lower = a->lower;
		node<dtype>* upper = a->upper;
		this->m_num_destroyed++;

		vector<arc<dtype>* > to_insert;

		// first create new ones
		arc<dtype>* lower_arc = lower->firstarc;
		while (lower_arc != NULL) {

			if (! isAlive(lower_arc) || lower_arc->lower != lower || lower_arc == a) {
				lower_arc = nextArc(lower, lower_arc);
				continue;
			}
			//printf("get here a\n");

			arc<dtype>* upper_arc = upper->firstarc;
			while (upper_arc != NULL ) {
				if (! isAlive(upper_arc)|| upper_arc->upper != upper || upper_arc == a) {
					upper_arc = nextArc(upper, upper_arc);
					continue;
				}
				////////////////
				// so both are alive, connect upper of lower to lower of upper
				//printf("get here\n");
				arc<dtype>* newarc = NewArc(lower_arc, a, upper_arc, m_num_destroyed);
				this->arcs.push_back(newarc);

				to_insert.push_back(newarc);
				// INSERT INTO ARC SORTER!

				//this->edges_to_cancel.push(newarc);
				////////////////
				upper_arc = nextArc(upper, upper_arc);
			}
			lower_arc = nextArc(lower, lower_arc);

		}
		// now destroy old ones
		lower_arc = lower->firstarc;
		while (lower_arc != NULL) {
			if (isAlive(lower_arc) && lower_arc != a) {
				if (lower_arc->lower == lower) {
					lower_arc->upper->dsc_man = 
						new manifold<dtype>(lower_arc->upper->dsc_man, upper->dsc_man, m_num_destroyed);
				}
				lower_arc->destroyed = this->m_num_destroyed;
			}
			lower_arc = nextArc(lower, lower_arc);
		}
		arc<dtype>* upper_arc = upper->firstarc;
		while (upper_arc != NULL) {
			if (isAlive(upper_arc) && upper_arc != a) {
				if (upper_arc->upper == upper) {
					upper_arc->lower->asc_man = 
						new manifold<dtype>(upper_arc->lower->asc_man, lower->asc_man, m_num_destroyed);
				}
				upper_arc->destroyed = this->m_num_destroyed;
			}
			upper_arc = nextArc(upper, upper_arc);
		}
		lower->destroyed = this->m_num_destroyed;
		upper->destroyed = this->m_num_destroyed;
		a->destroyed = this->m_num_destroyed;

		this->m_sel_persistence = this->m_num_destroyed;

		for (int i = 0; i < to_insert.size(); i++) {
			set_cancellation_weight(to_insert[i]);
			edges_to_cancel.push(to_insert[i]);
		}
	}

	vector<arc<dtype>*> cancel_history;

	// NOTE: this only works for 2d and 3d data sets.
	void recCollectMans(manifold<dtype>* m, set<CELL_INDEX_TYPE>& res) {
		if (m->basenode->index > this->m_mesh->max_dim()/2) {
			recCollectMans(m, res, 1);
		} else {
			recCollectMans(m, res, 0);
		}

	}
	void recCollectMans(manifold<dtype>* m, set<CELL_INDEX_TYPE>& res, set<CELL_INDEX_TYPE>& ids) {
		if (m->basenode->index > this->m_mesh->max_dim()/2) {
			recCollectMans(m, res, ids, 1);
		} else {
			recCollectMans(m, res, ids, 0);
		}

	}
	void recCounLeaftManifolds(manifold<dtype>* m, map< manifold<dtype>*, int >& counter) {
		if (m->merge[0] != NULL) {
			recCounLeaftManifolds(m->merge[0], counter);
			recCounLeaftManifolds(m->merge[1], counter);
		} else {
			if (counter.count(m) == 0) {
				counter[m] = 1;
			} else {
				counter[m]++;
			}
		}
	}

	void recCollectAscMansFill(manifold<dtype>* m, set<CELL_INDEX_TYPE>& res) {

		map<manifold<dtype>*,int> counter;
		recCounLeaftManifolds(m, counter);
		printf("found %d leaves\n", counter.size());
		for (typename map<manifold<dtype>*,int>::iterator it = counter.begin(); it != counter.end(); it++) {
			if ((*it).second % 2)
				this->rec_man_trace_up_fill((*it).first->basenode->cellid, res);				
		}
	}


	void recCollectDscMansFillRes(manifold<dtype>* m,set<CELL_INDEX_TYPE>& ids, set<CELL_INDEX_TYPE>& res) {

		map<manifold<dtype>*,int> counter;
		recCounLeaftManifolds(m, counter);
		for (typename map<manifold<dtype>*,int>::iterator it = counter.begin(); it != counter.end(); it++) {
			if ((*it).second % 2)
				this->rec_res_man_trace_down_fill((*it).first->basenode->cellid, ids, res);				
		}
	}


	void recCollectMans(manifold<dtype>* m, set<CELL_INDEX_TYPE>& res, int direction) {

		map<manifold<dtype>*,int> counter;
		recCounLeaftManifolds(m, counter);
		if (direction == 1) {
			for (typename map<manifold<dtype>*, int>::iterator it = counter.begin(); it != counter.end(); it++) {
				if ((*it).second % 2)
					this->rec_man_trace_down((*it).first->basenode->cellid, res);
			}
		} else {
			for (typename map<manifold<dtype>*,int>::iterator it = counter.begin(); it != counter.end(); it++) {
				if ((*it).second % 2)
					this->rec_man_trace_up((*it).first->basenode->cellid, res);
			}				
		}
	}	

	void recCollectMans(manifold<dtype>* m, set<CELL_INDEX_TYPE>& res, set<CELL_INDEX_TYPE>& ids, int direction) {

		//printf("new code gets called!\n");
		map< manifold<dtype>*, int> counter;
		recCounLeaftManifolds(m, counter);
		if (direction == 1) {
			for (typename map<manifold<dtype>*,int>::iterator it = counter.begin(); it != counter.end(); it++) {
				if ((*it).second % 2)
					this->rec_man_trace_down_restricted((*it).first->basenode->cellid, res, ids);
			}
		} else {
			for (typename map<manifold<dtype>*,int>::iterator it = counter.begin(); it != counter.end(); it++) {
				if ((*it).second % 2)
					this->rec_man_trace_up_restricted((*it).first->basenode->cellid, res, ids);
			}				
		}
	}

	//void recCollectMansRestricted(manifold<dtype>* m, set<CELL_INDEX_TYPE>& res, set<CELL_INDEX_TYPE>& ids) {
	//	if (m->basenode->index > 1) {
	//		recCollectMans(m, res, ids, 1);
	//	} else {
	//		recCollectMans(m, res, ids, 0);
	//	}

	//}

	//void recCollectMansRestricted(manifold<dtype>* m, set<CELL_INDEX_TYPE>& res, set<CELL_INDEX_TYPE>& ids, int direction) {

	//	if (m->merge[0] != NULL) {
	//		recCollectMansRestricted(m->merge[0], res, direction);
	//		recCollectMansRestricted(m->merge[1], res, direction);
	//	} else {
	//		if (direction == 1) {
	//			this->rec_man_trace_down_restricted(m->basenode->cellid, res);
	//		} else {
	//			this->rec_man_trace_down_restricted(m->basenode->cellid, res);
	//		}
	//	}
	//}
	void recCountGeometry(arc<dtype>* a, map<arc<dtype>*, int>& count) {

		// base geometry
		if (count.count(a) == 0) {
			count[a] = 1;
		} else {
			count[a]++;
		}
		if (a->base_geometry.size() != 0) {

			return;
		}
		recCountGeometry(a->merge[0],  count);
		recCountGeometry(a->merge[1],  count);
		recCountGeometry(a->merge[2],  count);
		//recFillGeometry(a, dir, res);
	}


	void recFillGeometry(arc<dtype>* a, bool dir, vector<CELL_INDEX_TYPE>& res) {

		if (a->base_geometry.size() != 0) {
			// base geometry
			if (dir) {

				int i = 0;
				while (res.size() > 0 &&
					i < a->base_geometry.size()-1 &&
					res[res.size()-1] == a->base_geometry[i]) {
						res.pop_back();
						i++;
				}
				if (i > 0) i--;
				while (i <= a->base_geometry.size() - 1) {
					res.push_back(a->base_geometry[i]);
					i++;
				}


				////res.insert(res.end(), a->base_geometry.begin(), a->base_geometry.end());
				//for (int i=0;i<a->base_geometry.size(); i++) {
				//	if (//i > 0 && i < a->base_geometry.size() - 1 && 
				//		res.size() > 0 && res[res.size()-1] == a->base_geometry[i]) {
				//			if (i < a->base_geometry.size()-1 && 
				//				res.size() > 1 && 
				//				res[res.size()-2] == a->base_geometry[i+1]) {
				//		res.pop_back();
				//			} else {
				//			}
				//	} else {
				//		res.push_back(a->base_geometry[i]);
				//	}
				//}
			} else {
				//res.insert(res.end(), a->base_geometry.rend(), a->base_geometry.rbegin());
				int i = a->base_geometry.size()-1;
				while (res.size() > 0 &&
					i >= 0 &&
					res[res.size()-1] == a->base_geometry[i]) {
						res.pop_back();
						i--;
				}
				if (i < a->base_geometry.size()-1) i++;
				while (i >= 0 ) {
					res.push_back(a->base_geometry[i]);
					i--;
				}

				//for (int i=a->base_geometry.size()-1;i>=0; i--) {
				//	if (//i > 0 && i < a->base_geometry.size() - 1 && 
				//		res.size() > 0 && res[res.size()-1] == a->base_geometry[i]) {
				//			if (i > 1 && res.size() > 1 && 
				//				res[res.size()-2] == a->base_geometry[i-1]) {
				//				res.pop_back();
				//			} else {
				//			}
				//	} else {
				//		res.push_back(a->base_geometry[i]);
				//	}
				//}
			}
			return;
		}
		if (dir) {
			recFillGeometry(a->merge[0], true, res);
			recFillGeometry(a->merge[1], false, res);
			recFillGeometry(a->merge[2], true, res);
			//recFillGeometry(a, dir, res);
		} else {
			recFillGeometry(a->merge[2], false, res);
			recFillGeometry(a->merge[1], true, res);
			recFillGeometry(a->merge[0], false, res);
			//recFillGeometry(a, dir, res);
		}

	}
	void recFillGeometry(arc<dtype>* a, bool dir, vector<CELL_INDEX_TYPE>& res,map<arc<dtype>*, int>& count) {

		if (count[a] % 2 == 0) return;

		if (a->base_geometry.size() != 0) {
			// base geometry
			if (dir) {

				int i = 0;
				while (res.size() > 0 &&
					i < a->base_geometry.size()-1 &&
					res[res.size()-1] == a->base_geometry[i]) {
						res.pop_back();
						i++;
				}
				if (i > 0) i--;
				while (i <= a->base_geometry.size() - 1) {
					res.push_back(a->base_geometry[i]);
					i++;
				}


				////res.insert(res.end(), a->base_geometry.begin(), a->base_geometry.end());
				//for (int i=0;i<a->base_geometry.size(); i++) {
				//	if (//i > 0 && i < a->base_geometry.size() - 1 && 
				//		res.size() > 0 && res[res.size()-1] == a->base_geometry[i]) {
				//			if (i < a->base_geometry.size()-1 && 
				//				res.size() > 1 && 
				//				res[res.size()-2] == a->base_geometry[i+1]) {
				//		res.pop_back();
				//			} else {
				//			}
				//	} else {
				//		res.push_back(a->base_geometry[i]);
				//	}
				//}
			} else {
				//res.insert(res.end(), a->base_geometry.rend(), a->base_geometry.rbegin());
				int i = a->base_geometry.size()-1;
				while (res.size() > 0 &&
					i >= 0 &&
					res[res.size()-1] == a->base_geometry[i]) {
						res.pop_back();
						i--;
				}
				if (i < a->base_geometry.size()-1) i++;
				while (i >= 0 ) {
					res.push_back(a->base_geometry[i]);
					i--;
				}

				//for (int i=a->base_geometry.size()-1;i>=0; i--) {
				//	if (//i > 0 && i < a->base_geometry.size() - 1 && 
				//		res.size() > 0 && res[res.size()-1] == a->base_geometry[i]) {
				//			if (i > 1 && res.size() > 1 && 
				//				res[res.size()-2] == a->base_geometry[i-1]) {
				//				res.pop_back();
				//			} else {
				//			}
				//	} else {
				//		res.push_back(a->base_geometry[i]);
				//	}
				//}
			}
			return;
		}
		if (dir) {
			recFillGeometry(a->merge[0], true, res,count);
			recFillGeometry(a->merge[1], false, res,count);
			recFillGeometry(a->merge[2], true, res,count);
			//recFillGeometry(a, dir, res);
		} else {
			recFillGeometry(a->merge[2], false, res,count);
			recFillGeometry(a->merge[1], true, res,count);
			recFillGeometry(a->merge[0], false, res,count);
			//recFillGeometry(a, dir, res);
		}

	}
public:
	virtual void setFirstNonzeroPersistence(float perc) {
		int pos = 1;
		while (this->cancel_history[pos]->persistence <= 0  && pos < this->cancel_history.size()) pos++;
		this->m_sel_persistence = pos-1;
		printf("pos is %d\n",pos-1);

	}
	void setPercentPersistence(float perc) {
		float maxval = this->cancel_history[this->cancel_history.size()-1]->persistence;
		float looky = maxval * 0.01 * perc;

		printf("maxval = %f, looky = %f, num = %d\n", maxval, looky, this->cancel_history.size());
		int pos = 0;
		while (pos < this->cancel_history.size() && this->cancel_history[pos]->persistence < looky) pos++;
		this->m_sel_persistence = pos;
		this->m_temp_perc_pers = perc;
		printf("pos is %d\n",pos);

	}
	void setAbsolutePersistence(float perc) {
		float maxval = this->cancel_history[this->cancel_history.size()-1]->persistence;
		float looky = perc;

		printf("maxval = %f, looky = %f, num = %d\n", maxval, looky, this->cancel_history.size());
		int pos = 0;
		while (pos < this->cancel_history.size() && this->cancel_history[pos]->persistence < looky) pos++;
		this->m_sel_persistence = pos;
		this->m_temp_perc_pers = perc;
		printf("pos is %d\n",pos);

	}
	float getPercentPersistence() {
		return this->m_temp_perc_pers;
	}
	void addPersistence(int val) {
		this->m_sel_persistence += val;
		if (m_sel_persistence < 0)m_sel_persistence=0;
		//printf("pos is %d\n",pos);

	}
	void fillGeometry(arc<dtype>* a, vector<CELL_INDEX_TYPE>& res) {

		map<arc<dtype>*, int> count;
		//printf("counting...\n");

		recCountGeometry(a, count);
		//printf("done\nfilling...\n");

		recFillGeometry(a, true, res, count);
		//printf("done\n");
	}

	manifold<dtype>* getActiveMan(manifold<dtype>* m) {
		while (m->created > this->m_sel_persistence) m = m->merge[0];
		return m;
	}

	void fillGeometry(node<dtype>* n, set<CELL_INDEX_TYPE>& res) {
		if (n->destroyed < this->m_sel_persistence) return;

		manifold<dtype>* man;
		if (n->index <= this->m_mesh->max_dim()/2) {
			man = n->asc_man;
		} else {
			man = n->dsc_man;
		}
		man = getActiveMan(man);

		recCollectMans(man, res);
	}

	void fillGeometryRestricted(node<dtype>* n, set<CELL_INDEX_TYPE>& res, set<CELL_INDEX_TYPE>& ids) {
		if (n->destroyed < this->m_sel_persistence) return;

		manifold<dtype>* man = n->man;
		while (man->created > this->m_sel_persistence) man = man->merge[0];

		recCollectMans(man, res, ids);
	}

	bool isAlive(arc<dtype>* a) {
		return a->created <= this->m_sel_persistence &&
			a->destroyed > this->m_sel_persistence;
	}
	bool isAlive(node<dtype>* n) {
		return n->destroyed > this->m_sel_persistence;
	}
	map<CELL_INDEX_TYPE, node<dtype>* > nodes;
	vector<arc<dtype>* > arcs;

	BasicMSC(mscBasicGradientField* grad, 
		mscBasicMeshHandler* mesh,
		mscBasicMeshFunction<dtype>* func) :
	m_grad(grad), m_mesh(mesh), m_func(func)
	{
		m_num_destroyed = 0;
		m_sel_persistence = 0;
	}

	void add_arc(node<dtype>* u, node<dtype>* l, 
		vector<CELL_INDEX_TYPE>& res) {
			arc<dtype>* a =
				new arc<dtype>(u, l, res);
			this->arcs.push_back(a);
	}

	//vector<node<dtype> > node_vec;

	void add_node(CELL_INDEX_TYPE id) {

		//node_vec.push_back(node<dtype>(id, m_mesh->dimension(id),
		//m_func->cell_value(id), m_mesh->boundary_value(id)));

		node<dtype>* n = new node<dtype>(id, m_mesh->dimension(id),
			m_func->cell_value(id), m_mesh->boundary_value(id));
		nodes[id] = n;

	}

	void ComputeFromGrad(bool restricted=false) {
		cellIterator t_it;
		// first add nodes
		iteratorOperator& t_cells = m_mesh->all_cells_iterator(t_it);
		for (t_cells.begin(t_it); t_cells.valid(t_it); t_cells.advance(t_it)) {
			CELL_INDEX_TYPE t_id = t_cells.value(t_it);
			if (m_grad->get_critical(t_id)) 
				add_node(t_id);
		}

		// then add arcs
		if (! restricted){
			for (typename map< CELL_INDEX_TYPE, node<dtype>* >::iterator nit = nodes.begin();
				nit != nodes.end(); nit++) {
					CELL_INDEX_TYPE t_id = (*nit).first;
					this->trace_down_cells_restricted(t_id);
			}
		} else {
			for (typename map< CELL_INDEX_TYPE, node<dtype>* >::iterator nit = nodes.begin();
				nit != nodes.end(); nit++) {
					CELL_INDEX_TYPE t_id = (*nit).first;
					if ((*nit).second->index == 1)
						this->trace_down_cells_restricted(t_id);

			}
		}
		ValidateComplex();
	}
	void WriteComplex(const char* filenamebase) {
		char outname[1024];

		/// nodes text file: each node on separate line
		// id dimension functionval boundaryflag X Y Z
		 
		printf("writing nodes\n");
		sprintf(outname, "%s.nodes.txt", filenamebase);
		FILE* fout = fopen(outname, "w");
		for (auto np : this->nodes) {

			node<dtype>* n = np.second;

			if (! isAlive(n)) continue;
			//printf("node %d\n", n->cellid);

			float coo[3];
			this->m_mesh->centroid(n->cellid, coo);

			fprintf(fout, "%d %d %f %d %.4f %.4f %.4f\n", n->cellid, n->index, n->value, n->boundary, coo[0] * 0.5, coo[1] * 0.5, coo[2] * 0.5);
		}
		fclose(fout);

		printf("writing arcs\n");
		/// arcs text file: each arc on separate line
		// kind lowernodeID uppernodeid lowerfunctionval upperfunctionval numpoints [x y z * numpoints]

		sprintf(outname, "%s.arcs.txt", filenamebase);
		fout = fopen(outname, "w");
		for (auto a : arcs) {
			//node<dtype>* b = this->nodes[i];
			if (!isAlive(a)) continue;

			vector<CELL_INDEX_TYPE> geom;
			this->fillGeometry(a, geom);

			fprintf(fout, "%d %d %d %f %f %d ", a->lower->index, a->lower->cellid, a->upper->cellid, a->lower->value, a->upper->value, geom.size());

			for (auto id : geom) {
				float coo[3];
				this->m_mesh->centroid(id, coo);
				fprintf(fout, "%.4f %.4f %.4f ", coo[0] * 0.5, coo[1] * 0.5, coo[2] * 0.5);
			}
			fprintf(fout, "\n");
	
		}
		fclose(fout);

		/// ascman text file: each minimum's ascending manifold on separate line
		// nodeid numpoints numpoints [x y z * numpoints]
		
		//sprintf(outname, "%s.ascman.txt", filenamebase);
		//fout = fopen(outname, "w");
		//for (auto a : arcs) {
		//	//node<dtype>* b = this->nodes[i];
		//	if (!isAlive(a)) continue;

		//	vector<CELL_INDEX_TYPE> geom;
		//	this->fillGeometry(a, geom);

		//	fprintf(fout, "%d %d %d %f %f %d ", a->lower->index, a->lower->cellid, a->upper->cellid, a->lower->value, a->upper->value, geom.size());

		//	for (auto id : geom) {
		//		float coo[3];
		//		this->m_mesh->centroid(id, coo);
		//		fprintf(fout, "%.4f %.4f %.4f ", coo[0] * 0.5, coo[1] * 0.5, coo[2] * 0.5);
		//	}
		//	fprintf(fout, "\n");

		//}
		//fclose(fout);


	}
	void WriteCancelHistory(const char* fnamebase) {
		int counts[4]; for (int i = 0; i < 4; i++) counts[i] = 0;
		for (auto n : this->nodes) {
			counts[n.second->index]++;
		}
		char fname[1024];
		sprintf(fname, "%s.pers", fnamebase);
		FILE* fout = fopen(fname, "w");
		for (int i = 0; i < this->cancel_history.size(); i++) {
			arc<dtype>* a = this->cancel_history[i];
			counts[a->lower->index]--;
			counts[a->upper->index]--;
			fprintf(fout, "%d %f %f %f %d %d %d %d\n", a->lower->index, a->lower->value,
				a->upper->value, a->persistence, counts[0], counts[1], counts[2], counts[3]);

		}
		fclose(fout);
	}
	void WriteCancelHistory() {
		WriteCancelHistory("PERSISTENCE.txt");
	}
	virtual void ComputeHeirarchy() {
		for (int i = 0; i < this->arcs.size(); i++) {
			this->set_cancellation_weight(this->arcs[i]);
			edges_to_cancel.push(this->arcs[i]);
		}

		while(! edges_to_cancel.empty()) {
			arc<dtype>* a = edges_to_cancel.top();
			edges_to_cancel.pop();
			//printf("isvalid?\n");
			if (validToCancel(a)) {
				//printf("%d -> %d, %f\n", a->lower->index, a->upper->index, a->persistence);
				cancel(a);
				this->cancel_history.push_back(a);
			}
		}
		//ValidateComplex();
	}
	void ComputeHeirarchy(float value) {
		for (int i = 0; i < this->arcs.size(); i++) {
			if ( this->arcs[i]->persistence <= value){
				this->set_cancellation_weight(this->arcs[i]);
				edges_to_cancel.push(this->arcs[i]);
			}
		}
		int counter  =0 ;
		while(! edges_to_cancel.empty()) {
			arc<dtype>* a = edges_to_cancel.top();
			if (a->persistence > value) return;
			edges_to_cancel.pop();
			//printf("isvalid?\n");
			if (validToCancel(a)) {
				//printf("%d-%d: %d -> %d, %f, %d, %d\n", counter++, nodes.size() - 2*counter, a->lower->index, a->upper->index, a->persistence, a->cancellation_weight, edges_to_cancel.size());
				cancel(a);
				this->cancel_history.push_back(a);
				//printf("done cancel\n");
			}
		}
		//ValidateComplex();
	}
	void ValidateComplex() {
		typename map<CELL_INDEX_TYPE, node<float>*>::iterator nit = nodes.begin();
		while (nit != nodes.end()) {
			node<float>* n = (*nit).second;
			nit++;	

			arc<dtype>* a = n->firstarc;
			while( a != NULL)  {
				if (a->lower != n && a->upper != n)
					printf("FUCKFUKCUF\n");

				a = this->nextArc(n, a);
			}
		}

		for (int i = 0; i < this->arcs.size(); i++) {
			arc<dtype>* a = this->arcs[i];
			if (a->lower->index != a->upper->index-1) {
				printf("ERROR: arc's lower id %d, upper %d\n", a->lower->index,
					a->upper->index);
			}
			arc<dtype>* la = a->lower->firstarc;
			bool has = false;
			while (la != NULL) {
				if (la == a) {
					has = true;
					break;
				}
				la = nextArc(a->lower, la);
			}
			if (! has) printf("ERROR lower node's  list does not have arc\n");
			la = a->upper->firstarc;
			has = false;
			while (la != NULL) {
				if (la == a) {
					has = true;
					break;
				}
				la = nextArc(a->upper, la);
			}
			if (! has) printf("ERROR upper node's list does not have arc\n");
		}



	}
};



template<typename dtype> 
class CrystalExtractingMSC : public BasicMSC<dtype> {


	void FillInCofacets(const set<CELL_INDEX_TYPE>& ids, set<CELL_INDEX_TYPE>& res) {
		for (set<CELL_INDEX_TYPE>::iterator it = ids.begin(); it != ids.end(); it++) {
			CELL_INDEX_TYPE cellid = *it;
			cellIterator cit;
			iteratorOperator& cofacets = this->m_mesh->cofacets(cellid, cit);
			for (cofacets.begin(cit); cofacets.valid(cit); cofacets.advance(cit)) {
				res.insert(cofacets.value(cit));
			}
		}
	}

	void FillInCofaces(set<CELL_INDEX_TYPE>& ids) {


	}

public:

	struct DCell {
		CELL_INDEX_TYPE upper;
		CELL_INDEX_TYPE lower;
		int D;
		set<DCell*> neighbors;	
	};

	typedef pair<CELL_INDEX_TYPE, CELL_INDEX_TYPE> CPair;
	map<CPair, DCell*> mCellRegistry;
	typedef pair<node<dtype>*, node<dtype>*>  CellPair;

	//set<CellPair> mCellRegistry;

	void Connect(DCell* A, DCell* B) {
		A->neighbors.insert(B);
		B->neighbors.insert(A);
	}

	void AddUpCrystalCells(node<dtype>* n, vector<node<dtype>*>& history) {

		history.push_back(n);

		int last = history.size()-1;
		for (int i = last; i >= 0; i--) {




			DCell* dc;
			//CellPair cell(history[i], n);
			CPair cell(history[i]->cellid, history[last]->cellid);
			if (mCellRegistry.count(cell) == 0) {
				dc = new DCell();
				//fill in dc
				dc->lower = cell.first;
				dc->upper = cell.second;
				dc->D = last - i;
				mCellRegistry[cell] = dc;
			} else {
				dc = mCellRegistry[cell];
			}

			if (i < last) {
				// connect to i + 1
				DCell* A = mCellRegistry[CPair(history[i+1]->cellid, history[last]->cellid)];
				DCell* B = mCellRegistry[CPair(history[i]->cellid, history[last-1]->cellid)];

				Connect(dc,A);
				Connect(dc,B);

			}




		}
		arc<dtype>* a = n->firstarc;
		while (a != NULL) {
			if (this->isAlive(a) && a->lower == n) {
				AddUpCrystalCells(a->upper, history);			
			}
			a = this->nextArc(n, a);
		}
		history.pop_back();
	}

	void AddFaces(set<CELL_INDEX_TYPE>& res) {
		set<CELL_INDEX_TYPE> source;

		source.insert(res.begin(), res.end());

		set<CELL_INDEX_TYPE> dest;

		int DIM = this->m_mesh->dimension(*(res.begin()));

		while(DIM > 0) {

			for (set<CELL_INDEX_TYPE>::iterator it = source.begin(); it != source.end(); it++) {
				cellIterator cit;
				iteratorOperator& facets = this->m_mesh->facets(*it, cit);
				for (facets.begin(cit); facets.valid(cit); facets.advance(cit)) {

					CELL_INDEX_TYPE fid = facets.value(cit);
					// check if all cofacets are also in shit
					bool allin = true;
					int count = 0;
					cellIterator cit2;
					iteratorOperator& cofacets = this->m_mesh->cofacets(fid, cit2);
					for (cofacets.begin(cit2); cofacets.valid(cit2); cofacets.advance(cit2)) {
						CELL_INDEX_TYPE fid2 = cofacets.value(cit2);
						if (res.count(fid2) == 0) {
							allin = false;
							break;
						}
						count++;
					}
					if (allin && count > 1) dest.insert(fid);
				}
			}

			res.insert(dest.begin(), dest.end());

			source.clear();
			source.insert(dest.begin(), dest.end());
			dest.clear();
			DIM--;
		}
	}
	void AddCoFaces(set<CELL_INDEX_TYPE>& res) {
		set<CELL_INDEX_TYPE> source;

		source.insert(res.begin(), res.end());

		set<CELL_INDEX_TYPE> dest;

		int DIM = this->m_mesh->dimension(*(res.begin()));

		while(DIM < this->m_mesh->max_dim()) {

			for (set<CELL_INDEX_TYPE>::iterator it = source.begin(); it != source.end(); it++) {
				cellIterator cit;
				iteratorOperator& cofacets = this->m_mesh->cofacets(*it, cit);
				for (cofacets.begin(cit); cofacets.valid(cit); cofacets.advance(cit)) {
					CELL_INDEX_TYPE fid = cofacets.value(cit);
					// check if all cofacets are also in shit
					bool allin = true;
					int count = 0;
					cellIterator cit2;
					iteratorOperator& facets = this->m_mesh->facets(fid, cit2);
					for (facets.begin(cit2); facets.valid(cit2); facets.advance(cit2)) {
						CELL_INDEX_TYPE fid2 = facets.value(cit2);
						if (res.count(fid2) == 0) {
							allin = false;
							break;
						}
						count++;
					}
					if (allin && count > 1) dest.insert(fid);
				}
			}

			res.insert(dest.begin(), dest.end());

			source.clear();
			source.insert(dest.begin(), dest.end());
			dest.clear();
			DIM++;

		}
	}

	void AddUpBoundaryNodes(node<dtype>* nl, set<node<dtype>*>& nodeplusboundary) {
		nodeplusboundary.insert(nl);

		for (int i = nl->index; i < this->m_mesh->max_dim(); i++) {
			set<node<dtype>*> temp;
			for (typename set<node<dtype>*>::iterator it = nodeplusboundary.begin();
				it != nodeplusboundary.end(); it++) {
					if ((*it)->index != i) continue;

					arc<dtype>* a = (*it)->firstarc;
					while (a != NULL) {
						if (isAlive(a) && a->lower == (*it)) {
							temp.insert(a->upper);
						}
						a = nextArc((*it), a);
					}
			}
			nodeplusboundary.insert(temp.begin(), temp.end());
		}

	}

	void AddDownBoundaryNodes(node<dtype>* nl, set<node<dtype>*>& nodeplusboundary) {
		nodeplusboundary.insert(nl);

		for (int i = nl->index; i > 0; i--) {
			set<node<dtype>*> temp;
			for (typename set<node<dtype>*>::iterator it = nodeplusboundary.begin();
				it != nodeplusboundary.end(); it++) {
					if ((*it)->index != i) continue;

					arc<dtype>* a = (*it)->firstarc;
					while (a != NULL) {
						if (isAlive(a) && a->upper == (*it)) {
							temp.insert(a->lower);
						}
						a = nextArc((*it), a);
					}
			}
			nodeplusboundary.insert(temp.begin(), temp.end());
		}

	}

	void CellGeometry(DCell* d, set<CELL_INDEX_TYPE>& result) {

		// first find the right merge node

		// gather the nodes in the boundary of the descending manifold

		printf("CellGeometry called: %d-%d\n", d->lower, d->upper);
		node<dtype>* nl = this->nodes[d->lower];
		set<node<dtype>*> nodeplusboundary;
		//AddUpBoundaryNodes(nl, nodeplusboundary);
		nodeplusboundary.insert(nl);

		printf("-- found %d nodes in up boundary nodes\n", nodeplusboundary.size());

		// now nodeplusboundary has the list of all nodes belonging to the boundary. 
		// so collect their descending manifolds.
		set<CELL_INDEX_TYPE> upids;
		for(typename set<node<dtype>*>::iterator nit = nodeplusboundary.begin();
			nit != nodeplusboundary.end(); nit++) {
				node<dtype>* n = *nit;
				manifold<dtype>* man = n->asc_man;
				while (man->created > this->m_sel_persistence) man = man->merge[0];
				// collect upward 
				this->recCollectAscMansFill(man, upids);
				printf("-- id size: %d\n", upids.size());
		}

		node<dtype>* nu = this->nodes[d->upper];
		set<node<dtype>*> nodeplusboundary2;
		nodeplusboundary2.insert(nu);
		//AddDownBoundaryNodes(nu, nodeplusboundary2);
		printf("-- found %d nodes in down boundary nodes\n", nodeplusboundary2.size());
		//set<CELL_INDEX_TYPE> bothids;
		for(typename set<node<dtype>*>::iterator nit = nodeplusboundary2.begin();
			nit != nodeplusboundary2.end(); nit++) {
				node<dtype>* n = *nit;
				manifold<dtype>* man = n->dsc_man;
				while (man->created > this->m_sel_persistence) man = man->merge[0];
				// collect downward 
				this->recCollectDscMansFillRes(man, upids, result);
				printf("-- res size: %d\n", result.size());
		}

		//// now intersect
		//set<CELL_INDEX_TYPE>::iterator it1 = upids.begin();
		//set<CELL_INDEX_TYPE>::iterator it2 = downids.begin();
		//// sets are sorted so mergesort
		//while (it1 != upids.end() && it2!= downids.end()) {
		//	if (*it1 < *it2) ++it1;
		//	else if (*it2 < *it1) ++it2;
		//	else {
		//		result.insert(result.end(), *it1);
		//		++it1; ++it2;
		//	}
		//}


	}


	void RegisterCells() {
		mCellRegistry.clear();

		vector<node<dtype>*> history;
		for (typename map<CELL_INDEX_TYPE, node<dtype>* >::iterator it = this->nodes.begin();
            it != this->nodes.end(); it++) {
			node<dtype>* n = (*it).second;
			if(n->index == 0 && this->isAlive(n)) {
				AddUpCrystalCells(n, history);
			}
		}

		printf("Cell Registry Size: %d\n", mCellRegistry.size());
	}


	CrystalExtractingMSC(mscBasicGradientField* grad, 
		mscBasicMeshHandler* mesh,
		mscBasicMeshFunction<dtype>* func) : BasicMSC<dtype>(grad, mesh, func)
	{
	}



};

template<typename dtype>
class RestrictedMSC : public BasicMSC<dtype> {
public:
	mscBasicGradientField* mBoundary;
	RestrictedMSC(mscBasicGradientField* grad, 
		mscBasicMeshHandler* mesh,
		mscBasicMeshFunction<dtype>* func, mscBasicGradientField* boundarymesh) : BasicMSC<dtype>(grad, mesh, func)
	{
		mBoundary = boundarymesh;
	}

	virtual void setFirstNonzeroPersistence(float perc) {
		int pos = 1;
		while (this->cancel_history[pos]->persistence < 0  && pos < this->cancel_history.size()) pos++;
		this->m_sel_persistence = pos-1;
		printf("pos is %d\n",pos-1);

	}


	virtual arc<dtype>* NewArc( arc<dtype>* l_upper, arc<dtype>* a, arc<dtype>* u_lower, CELL_INDEX_TYPE created) {
		arc<dtype>* na = new arc<dtype>(l_upper, a, u_lower, created);
		if (mBoundary->get_dim_asc_man(na->lower->cellid) == 
			mBoundary->get_dim_asc_man(na->upper->cellid)){
				a->persistence = -0.01;
		} else {
			//a->persistence += 0.1;
		}
		return na;		
	}

	virtual void ComputeHeirarchy() {
		for (int i = 0; i < this->arcs.size(); i++) {
			arc<dtype>* a = this->arcs[i];
			this->set_cancellation_weight(a);
			if (mBoundary->get_dim_asc_man(a->lower->cellid) == 
				mBoundary->get_dim_asc_man(a->upper->cellid)) {
					a->persistence = -0.01;
			} else {
				//a->persistence += 0.1;
			}

			this->edges_to_cancel.push(this->arcs[i]);
		}

		while(! this->edges_to_cancel.empty()) {
			arc<dtype>* a = this->edges_to_cancel.top();
			this->edges_to_cancel.pop();
			//printf("isvalid?\n");
			if (this->validToCancel(a)) {
				//printf("%d -> %d, %f\n", a->lower->index, a->upper->index, a->persistence);
				this->cancel(a);
				this->cancel_history.push_back(a);
			}
		}
		//ValidateComplex();
	}


};


template<typename dtype>
class TalassMSC : public BasicMSC<dtype> {
public:
	TalassMSC(mscBasicGradientField* grad,
		mscBasicMeshHandler* mesh,
		mscBasicMeshFunction<dtype>* func) : BasicMSC<dtype>(grad, mesh, func)
	{
	}


	virtual void ComputeHeirarchy() {
		for (int i = 0; i < this->arcs.size(); i++) {
			arc<dtype>* a = this->arcs[i];
			this->set_cancellation_weight(a);
			if (a->lower->index != 0) {
				continue;
			}
			else {
				//a->persistence += 0.1;
			}

			this->edges_to_cancel.push(this->arcs[i]);
		}

		while (!this->edges_to_cancel.empty()) {
			arc<dtype>* a = this->edges_to_cancel.top();
			this->edges_to_cancel.pop();
			//printf("isvalid?\n");
			if (validToCancel(a)) {
				//printf("%d -> %d, %f\n", a->lower->index, a->upper->index, a->persistence);
				cancel(a);
				this->cancel_history.push_back(a);
			}
		}
		//ValidateComplex();
	}

	//void CreateTALASSFormatFromMorseComplex() {
	//	vector<vector<CELL_INDEX_TYPE> > segments;
	//	vector<node<dtype>*> node_vec;
	//	for (node<dtype>* n : nodes) {
	//		if (n->index == 0 && isAlive(n)) {
	//			node_vec.push_back(n);
	//			set<CELL_INDEX_TYPE> geom;
	//			this->fillGeometry(n, geom);
	//			vector<CELL_INDEX_TYPE> geomv;
	//			for (auto id : geom) {
	//				if (m_mesh->dimension(id) == 0) {
	//					geomv.push_back(id); // want to convert
	//	}





	//}
};

#endif


