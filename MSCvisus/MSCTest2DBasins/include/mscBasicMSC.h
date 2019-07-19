#ifndef MSC_BASIC_MSC
#define MSC_BASIC_MSC


#include <stdio.h>
#include <vector>
#include <map>
#include <queue>

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
	manifold<dtype>* man;
	
	node(CELL_INDEX_TYPE id, DIM_TYPE index, dtype value, BOUNDARY_TYPE boundary) :
		value(value), cellid(id), index(index), firstarc(NULL), boundary(boundary),
		destroyed(MAX_CELLID_VALUE) {
			man = new manifold<dtype>(this);
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

	arc(node<dtype>* upper, node<dtype>* lower, 
		vector<CELL_INDEX_TYPE> geometry) : lower(lower), upper(upper),
	lower_next(NULL), upper_next(NULL), created(0), destroyed(MAX_CELLID_VALUE) {
		base_geometry.insert(base_geometry.begin(), geometry.begin(), geometry.end());
		persistence = upper->value - lower->value;
		lower_next = lower->firstarc;
		lower->firstarc = this;
		upper_next = upper->firstarc;
		upper->firstarc = this;
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

			persistence = upper->value - lower->value;
	}


};


template<typename dtype>
class BasicMSC {

	struct arcCompare {
		bool operator()(const arc<dtype>* left, const arc<dtype>* right) {
			if (left->persistence < right->persistence) return false;
			if (left->persistence > right->persistence ) return true;
			return left < right;
		}
	};
	priority_queue<arc<dtype>*, vector< arc<dtype>* >, arcCompare > edges_to_cancel;

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
		if (! isAlive(a)) return false;
		if (a->lower->boundary != a->upper->boundary) return false;
		if (countMultiplicity(a) != 1) return false;
		return true;
	}
		



	arc<dtype>* nextArc(node<dtype>* n, arc<dtype>* a) {
		if (a->lower == n) return a->lower_next;
		return a->upper_next;
	}

	void cancel(arc<dtype>* a) {

		//printf("cancelling %d-%d, %f\n", a->lower->index, a->upper->index, a->persistence);

		node<dtype>* lower = a->lower;
		node<dtype>* upper = a->upper;
		this->m_num_destroyed++;

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
				arc<dtype>* newarc = new arc<dtype>(lower_arc, a, upper_arc, m_num_destroyed);
				arcs.push_back(newarc);

		// INSERT INTO ARC SORTER!

				this->edges_to_cancel.push(newarc);
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
					lower_arc->upper->man = 
						new manifold<dtype>(lower_arc->upper->man, upper->man, m_num_destroyed);
				}
				lower_arc->destroyed = this->m_num_destroyed;
			}
			lower_arc = nextArc(lower, lower_arc);
		}
		arc<dtype>* upper_arc = upper->firstarc;
		while (upper_arc != NULL) {
			if (isAlive(upper_arc) && upper_arc != a) {
				if (upper_arc->upper == upper) {
					upper_arc->lower->man = 
						new manifold<dtype>(upper_arc->lower->man, lower->man, m_num_destroyed);
				}
				upper_arc->destroyed = this->m_num_destroyed;
			}
			upper_arc = nextArc(upper, upper_arc);
		}
		lower->destroyed = this->m_num_destroyed;
		upper->destroyed = this->m_num_destroyed;
		a->destroyed = this->m_num_destroyed;

		this->m_sel_persistence = this->m_num_destroyed;
	}

	vector<arc<dtype>*> cancel_history;


	void recCollectMans(manifold<dtype>* m, set<CELL_INDEX_TYPE>& res) {
	
		if (m->merge[0] != NULL) {
			recCollectMans(m->merge[0], res);
			recCollectMans(m->merge[1], res);
		} else {
			if (m->basenode->index == 2) {
				this->rec_man_trace_down(m->basenode->cellid, res);
			} else {
				this->rec_man_trace_up(m->basenode->cellid, res);
			}
		}
	}
	




	void recFillGeometry(arc<dtype>* a, bool dir, vector<CELL_INDEX_TYPE>& res) {

		if (a->base_geometry.size() != 0) {
			// base geometry
			if (dir) {
				for (int i=0;i<a->base_geometry.size(); i++) {
					if (res.size() > 0 && res[res.size()-1] == a->base_geometry[i]) {
						res.pop_back();
					} else {
						res.push_back(a->base_geometry[i]);
					}
				}
			} else {
				for (int i=a->base_geometry.size()-1;i>=0; i--) {
					if (res.size() > 0 && res[res.size()-1] == a->base_geometry[i]) {
						res.pop_back();
					} else {
						res.push_back(a->base_geometry[i]);
					}
				}
			}
			return;
		}
		if (dir) {
			recFillGeometry(a->merge[0], true, res);
			recFillGeometry(a->merge[1], false, res);
			recFillGeometry(a->merge[2], true, res);
		} else {
			recFillGeometry(a->merge[2], false, res);
			recFillGeometry(a->merge[1], true, res);
			recFillGeometry(a->merge[0], false, res);
		
		}

	}

public:
	void setFirstNonzeroPersistence(float perc) {
		int pos = 1;
		while (cancel_history[pos]->persistence <= 0  && pos < cancel_history.size()) pos++;
		this->m_sel_persistence = pos-1;
		printf("pos is %d\n",pos-1);

	}
	void setPercentPersistence(float perc) {
		float maxval = cancel_history[cancel_history.size()-1]->persistence;
		float looky = maxval * 0.01 * perc;

		printf("maxval = %f, looky = %f, num = %d\n", maxval, looky, cancel_history.size());
		int pos = 1;
		while (cancel_history[pos]->persistence < looky && pos < cancel_history.size()) pos++;
		this->m_sel_persistence = pos-1;
		printf("pos is %d\n",pos);

	}
	void addPersistence(int val) {
		this->m_sel_persistence += val;
		if (m_sel_persistence < 0)m_sel_persistence=0;
		//printf("pos is %d\n",pos);

	}
	void fillGeometry(arc<dtype>* a, vector<CELL_INDEX_TYPE>& res) {
		recFillGeometry(a, true, res);
	}

		void fillGeometry(node<dtype>* n, set<CELL_INDEX_TYPE>& res) {
		if (n->destroyed < this->m_sel_persistence) return;

		manifold<dtype>* man = n->man;
		while (man->created > this->m_sel_persistence) man = man->merge[0];

		recCollectMans(man, res);
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
			arcs.push_back(a);
	}
	void add_node(CELL_INDEX_TYPE id) {

		node<dtype>* n = new node<dtype>(id, m_mesh->dimension(id),
			m_func->cell_value(id), m_mesh->boundary_value(id));
		nodes[id] = n;

	}

	void ComputeFromGrad() {
		cellIterator t_it;
		// first add nodes
		iteratorOperator& t_cells = m_mesh->all_cells_iterator(t_it);
		for (t_cells.begin(t_it); t_cells.valid(t_it); t_cells.advance(t_it)) {
			CELL_INDEX_TYPE t_id = t_cells.value(t_it);
			if (m_grad->get_critical(t_id)) 
				add_node(t_id);
		}

		// then add arcs
		for (typename map< CELL_INDEX_TYPE, node<dtype>* >::iterator nit = nodes.begin();
			nit != nodes.end(); nit++) {
			CELL_INDEX_TYPE t_id = (*nit).first;

			this->trace_down_cells_restricted(t_id);

		}
		ValidateComplex();
	}
	
	void WriteCancelHistory() {
		FILE* fout = fopen("PERSISTENCE.txt", "w");
		for (int i = 0; i < cancel_history.size(); i++) {
			arc<dtype>* a = cancel_history[i];

			fprintf(fout, "%f %f %f\n", a->lower->value, 
				a->upper->value, a->persistence);

		}
		fclose(fout);
	}
	void ComputeHeirarchy() {
		for (int i = 0; i < arcs.size(); i++) {
			edges_to_cancel.push(arcs[i]);
		}

		while(! edges_to_cancel.empty()) {
			arc<dtype>* a = edges_to_cancel.top();
			edges_to_cancel.pop();
			//printf("isvalid?\n");
			if (validToCancel(a)) {
				//printf("%d -> %d, %f\n", a->lower->index, a->upper->index, a->persistence);
				cancel(a);
				cancel_history.push_back(a);
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

		for (int i = 0; i < arcs.size(); i++) {
			arc<dtype>* a = arcs[i];

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





#endif


