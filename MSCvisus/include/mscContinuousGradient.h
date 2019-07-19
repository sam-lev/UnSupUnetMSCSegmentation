#ifndef MSC_CONTINUOUS_GRADIENT
#define MSC_CONTINUOUS_GRADIENT

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscArrayFactory.h"
#include "mscConvergentGradientBuilder.h"
#include "mscRegularRawDataHandler.h"


#include <vector>
#include <math.h>

#include "ap.h"
#include "interpolation.h"

using namespace alglib_impl;
using namespace std;

// assume all stuff is doubles
struct dvect2d {
	double x;
	double y;
};

inline double dot2d(dvect2d& a, dvect2d& b) {
	return a.x * b.x + a.y * b.y;
};

inline double magnitude2d(dvect2d& a) {
	return sqrt(a.x * a.x + a.y * a.y);
};

inline void normalize2d(dvect2d& a) {
	double div = 1.0 / magnitude2d(a);
	a.x *= div;
	a.y *= div;
};

inline dvect2d minus2d(dvect2d& a, dvect2d& b) {
	dvect2d res;
	res.x = a.x - b.x;
	res.y = a.y - b.y;
	return res;
};

//class mscGeneralContinousGradient


class mscContinuousGradient2D {
protected:
	double* mValues;
	alglib::real_1d_array mArray;
	CELL_INDEX_TYPE mNumElements;
	alglib::real_1d_array mXs;
	alglib::real_1d_array mYs;
	alglib::ae_int_t mX;
	alglib::ae_int_t mY;

	alglib::spline2dinterpolant mSpline;
public:




	bool load_data(char* filename, CELL_INDEX_TYPE x, CELL_INDEX_TYPE y) {
		CELL_INDEX_TYPE num_elements = x*y;
		this->mNumElements = num_elements;
		FILE* fin = fopen(filename, "rb");
		if (fin == NULL) return false;

		double* tValues = new double[num_elements];

		for (CELL_INDEX_TYPE i = 0; i < num_elements; i++) {
			float ftmp;
			fread(&ftmp, sizeof(float), 1, fin);
			tValues[i] = (double) ftmp;
		}
		fclose(fin);

		this->mArray.setcontent(num_elements, tValues);

		delete[] tValues;

		tValues = new double[x];
		for (int i = 0; i < x; i++) tValues[i] = (double) i;
		this->mXs.setcontent(x, tValues);
		delete[] tValues;

		tValues = new double[y];
		for (int i = 0; i < y; i++) tValues[i] = (double) i;
		this->mYs.setcontent(y, tValues);
		delete[] tValues;

		this->mX = x;
		this->mY = y;

		alglib::spline2dbuildbicubicv(this->mXs, this->mX,
			this->mYs, this->mY, 
			this->mArray, 1, this->mSpline);
		printf("done!\n");

		double v, dx, dy, dxy;

		//double V2[2] = {.2, .4};
		//double V1[2];
		//for (int i = 0; i < 50; i++) {
		// alglib::spline2ddiff(this->mSpline, (double) i, (double) i, v, dx, dy, dxy); 
		//V1[0] = dx; V1[1] = dy;
		//printf("(%d, %d) = %.4f, (%.4f, %.4f), %.4f\n", i, i, v, dx, dy, alglib::vdotproduct(V1, V2, 2));
		//}

		return true;
	}

};



// for each cell of each dimension, list "downwards" directions
typedef  unsigned char dependency;
class mscVectorToScalarGridFunc2D : public mscBasicDataHandler<float> {

protected:

	class graphedgenode;

	class graphnode {
	public:
		CELL_INDEX_TYPE id;
		float value;
		graphedgenode* nodes;
		int edgecount;
		graphedgenode* depends;
		int flag;

		graphnode(CELL_INDEX_TYPE mid) : flag(0), id(mid), value(0.0f), edgecount(0),nodes(NULL), depends(NULL){}
		void add_node(graphnode* n, double magnitude) {
			edgecount++;
			graphedgenode* m = new graphedgenode( n, this, nodes, magnitude);
			nodes = m;
			graphedgenode* ss = new graphedgenode(this, n, n->depends, magnitude);
			n->depends = ss;
		}
	};

	class graphedgenode {
	public:
		graphnode* node;
		graphnode* prev;
		graphedgenode* next;
		double magnitude;
		bool visited;
		graphedgenode(graphnode* n, graphnode* p, graphedgenode* e, double mag) : 
		prev(p), node(n), next(e), magnitude(mag), visited(false) {}


	};

	dependency* cells;

	graphnode** nodes;
	vector<graphedgenode*> edges;


	int ddim(CELL_INDEX_TYPE cellid) {
		return (cellid % DX) +
			((cellid / DX) % DY);
	}

	CELL_INDEX_TYPE g2d(CELL_INDEX_TYPE cellid) {
		return ((cellid%DX) / 2) + (((cellid / DX) % DY) / 2) * XMAX;
	}
	queue<graphnode*> nnn;

public:

	vector<CELL_INDEX_TYPE> barfedges;
	vector<CELL_INDEX_TYPE> barfproc;
	inline virtual float value(CELL_INDEX_TYPE index) { return nodes[index]->value;}

	int XMAX, YMAX, DX, DY;
	mscBasicMeshHandler* mesh;
	dvect2d* vectors;
	void coordinates(CELL_INDEX_TYPE cellid, dvect2d& c){
		c.x = (cellid % DX);
		c.y = ((cellid / DX) % DY);
	}

	virtual dvect2d getVector(CELL_INDEX_TYPE cid) {
		return vectors[g2d(cid)];
	}


	void remove(graphedgenode* e) {

		printf("called\n");
		graphnode* nodeE = e->node;
		graphnode* nodeS = e->prev;
		printf("e=%d, s=%d, ee=%d\n", nodeE->id, nodeS->id, nodeS->nodes);

		graphedgenode* ee = nodeS->nodes;
		if (ee == e) { 
			nodeS->nodes = nodeS->nodes->next; 
		} else { 
			graphedgenode* pp = ee;
			while (ee != e) { pp =ee; ee = ee->next; }
			pp->next = ee->next;
		}
		nodeS->edgecount--;
		if (nodeS->edgecount == 0) nnn.push(nodeS);
		printf("2\n");

		//graphedgenode* ee2 = nodeE->depends;
		//if (ee2->node == nodeS) {
		//	nodeE->depends = nodeE->depends->next;
		//} else {
		//	graphedgenode* pp = ee2;
		//	while (ee2->node != nodeS) { pp == ee2; ee2 = ee2->next; }
		//	pp->next = ee2->next;
		//}
	}

	mscVectorToScalarGridFunc2D(double* vin, int x, int y, mscBasicMeshHandler* mh) {
		XMAX = x;
		YMAX = y;
		DX = x*2-1;
		DY = y*2-1;
		mesh = mh;
		nodes = new graphnode*[x*y];
		vectors = new dvect2d[x*y];
		for (int i = 0; i < x*y; i++) {
			vectors[i].x = vin[i*2];
			vectors[i].y = vin[i*2+1];
		}
	}

	vector<CELL_INDEX_TYPE> barfpoop;

	bool destroy_loop(graphnode* n, vector<graphedgenode*>& stack) {

		n->flag = 2;

		graphedgenode* ee = n->nodes;
		while (ee != NULL) {
			if (ee->node->flag == 0 /*&& ! ee->visited*/) {
				stack.push_back(ee);
				//printf("adding %d\n", ee->node->id);
				if (! destroy_loop(ee->node, stack)) {
					stack.pop_back();
					n->flag = 0;
					return false;
				}
				//ee->visited = true;
				stack.pop_back();
			} else if (ee->node->flag == 2) {
				printf("detected cycle...\n");
				double minmag = ee->magnitude;
				graphedgenode* mm = ee;
				int clok = stack.size()-1;
				while (clok >=0 && stack[clok]->node != ee->node) {
					printf("-%d\n", mm->node->id);
					if (stack[clok]->magnitude < minmag) {
						minmag = stack[clok]->magnitude;
						printf("----\n");
						//mm2 = stack[clok-1];
						mm = stack[clok];
					}
					clok--;
				}
				printf("nnnnn %d %d\n", mm->node->id, mm->prev->id);
				barfpoop.push_back(mm->prev->id);
				barfpoop.push_back(mm->node->id);
				// so mm is culprit!, remove it
				remove(mm);
				n->flag = 0;
				return false;
			}
			ee = ee->next;
		}
		n->flag = 3;
		return true;
	}

	void region_grow_assign_vals() {


		float counter = 1.0f;

		// initialize queue
		for (int i =0; i < XMAX*YMAX; i++) {
			if (nodes[i]->edgecount == 0) {
				nnn.push(nodes[i]);
			}
		}
		// traverse queue
		while (! nnn.empty()) {
			graphnode* nn = nnn.front(); nnn.pop();
			barfproc.push_back(nn->id);
			nn->flag = 1;
			nn->value = counter;
			counter += 1.0f;
			graphedgenode* ee = nn->depends;
			while (ee != NULL) {
				int count = ee->node->edgecount;
				ee->node->edgecount--;
				if (count == 1) {
					nnn.push(ee->node);
				}
				ee = ee->next;
			}

		}
		printf (" count = %f\n", counter); 
		// now fix any missing loops
		for (int i =0; i < XMAX*YMAX; i++) {
			// not processed. 
			if (nodes[i]->flag == 0) {
				vector<graphedgenode*> stack;
				printf("doing xxxx%d\n",i);
				while (! destroy_loop(nodes[i], stack)) {}
			}
		}
		while (! nnn.empty()) {
			graphnode* nn = nnn.front(); nnn.pop();
			barfproc.push_back(nn->id);
			nn->flag = 1;
			nn->value = counter;
			counter += 1.0f;
			graphedgenode* ee = nn->depends;
			while (ee != NULL) {
				int count = ee->node->edgecount;
				ee->node->edgecount--;
				if (count == 1) {
					nnn.push(ee->node);
				}
				ee = ee->next;
			}

		}
		printf (" count = %f\n", counter); 
	}

	void init() {
		printf("HERE IN INIT\n");
		// create all the nodes
		for (int i =0; i < XMAX*YMAX; i++) {
			nodes[i] = new graphnode(i);
		}
		barfedges.clear();
		// iterate over each edge in dataset and assign to one endpoint
		cellIterator it;
		iteratorOperator& dcells = mesh->d_cells_iterator(1, it); 
		for (dcells.begin(it); dcells.valid(it); dcells.advance(it)) {
			CELL_INDEX_TYPE cid = dcells.value(it);	


			struct t_pair { double mag; CELL_INDEX_TYPE id; };
			t_pair tbarf[2];
			int counter= 0;

			dvect2d edgeC;
			coordinates(cid, edgeC);
			//printf("doing %d\n", cid);

			cellIterator fit;
			iteratorOperator& facets = mesh->facets(cid, fit);
			for (facets.begin(fit); facets.valid(fit); facets.advance(fit)) {
				CELL_INDEX_TYPE vid = facets.value(fit);

				//printf("a\n");
				// now cid is cell id of edge, and vid is cell id of vertex
				dvect2d vertC;
				//printf("a\n");
				coordinates(vid, vertC);
				//printf("a\n");

				// reverse direction to get negative gradient
				dvect2d diffC = minus2d(vertC, edgeC);
				//printf("a\n");

				CELL_INDEX_TYPE vdid = g2d(vid);
				dvect2d gradV = vectors[vdid];
				//printf("a\n");

				double magnitude = dot2d(gradV, diffC);

				tbarf[counter].mag = magnitude;
				tbarf[counter].id = vdid;
				//printf("e\n");
				counter++;
			}

			// now tbarf has magnitude of either side
			// add the appropriate edges to dependency graph
			if (tbarf[0].mag > tbarf[1].mag && tbarf[0].mag > 0) {
				nodes[tbarf[0].id]->add_node(nodes[tbarf[1].id], tbarf[0].mag);

				barfedges.push_back(tbarf[0].id);
				barfedges.push_back(tbarf[1].id);

			} else if (tbarf[1].mag > tbarf[0].mag && tbarf[1].mag > 0){
				nodes[tbarf[1].id]->add_node(nodes[tbarf[0].id], tbarf[1].mag);

				barfedges.push_back(tbarf[1].id);
				barfedges.push_back(tbarf[0].id);


			} 

		}

		// NOW IN THEORY WE HAVE WHOLE DEPENDENCY GRAPH
		printf("OUT OF INIT\n");
		region_grow_assign_vals();
	}
};











class mscConvergentVFBuilder : public mscConvergentGradientBuilder<float> {
protected:
	// return index in candidates of pair
	virtual int pick_from_candidates(const CELL_INDEX_TYPE& cellid, 
		const vector<CELL_INDEX_TYPE>& candidates,
		const vector<MemberDist>& mdcand) {

			// find weights for computing my probability
			float cell_value = this->my_mesh_function->cell_value(cellid);
			//-->
			dvect2d cell_grad = this->mVF->getVector(cellid);

			int result = 0; 
			vector<float> values;
			float temp_sum = 0;

			dvect2d oV;
			this->mVF->coordinates(cellid, oV);

			for (int i = 0; i < candidates.size(); i++) {
				//float diff_val = (float) (cell_value - this->lowest_facet_value(candidates[i]));
				dvect2d dV;
				this->mVF->coordinates(candidates[i], dV);

				dvect2d ddV = minus2d(oV, dV);

				float diff_val = (float) ( dot2d(cell_grad, ddV));
				if (diff_val < 0.0f) diff_val = 0.0f;
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

#ifdef FANCY_PROB_OUTPUT
			idfpair p; p.id = cellid;
			p.prob = this->mdMax(my_dist).prob;
			maxvals.push_back(p);
#endif
			// compute my local weights

			int maxloc = 0; 
			vector<float> tress;
			float maxval = mdDot3(mdcand[0], my_dist);
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
			//for(int i = 0; i < mdcand.size(); i++) {
			//	if (i == maxloc) continue;
			//	if (tress[maxloc] == tress[i]) {
			//		// pick lower value
			//		if (this->lowest_facet_value(candidates[i]) <= this->lowest_facet_value(candidates[maxloc])){
			//			//printf("EHHEHEHE\n");
			//			maxloc = i;
			//		}
			//	}
			//}


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

			this->my_dists[cellid] = my_dist;
			return maxloc;			


	}

	virtual void pick_and_pair(comparison_element& element,
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
					this->my_grad_field->get_num_unpaired_facets(temp_cell) == 1 &&
					element.value >= this->my_mesh_function->cell_value(temp_cell) &&
					element.boundary == this->my_mesh_handler->boundary_value(temp_cell)) {
						MemberDist res;
						if (mdCombineAllFacetsAndDecrement(temp_cell, dim, res)) {
							mdcand.push_back(res);
							candidates.push_back(temp_cell);
						}
				}
			}

			add_neighbors_to_sort(element.cellid);

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

	mscVectorToScalarGridFunc2D* mVF;
public:
	mscConvergentVFBuilder (mscVectorToScalarGridFunc2D* vects,
		mscBasicMeshFunction<float>* mesh_function,
		mscBasicMeshHandler* mesh_handler,
		mscBasicGradientField* grad_field,
		mscArrayFactory* array_factory) : 
	mscConvergentGradientBuilder<float>(mesh_function, mesh_handler, grad_field, array_factory) {
		my_erase = true;
		mVF = vects;
	}
};


//class mscContGradMeshFunction : public mscBasicMeshFunction<float> {
//protected:
//	mscContinuousGradient2D* mGrad;
//	
//
//
//
//
//
//
//
//public:
//
//
//	mscContGradMeshFunction(mscContinuousGradient2D* grad) {
//		this->mGrad = grad;
//	}
//
//	virtual ~mscContGradMeshFunction(){
//			printf("delete: mscContGradMeshFunction \n");
//	}
//	virtual void initialize() = 0;
//	virtual float cell_value(CELL_INDEX_TYPE cellid)  = 0;
//	virtual bool less_than(CELL_INDEX_TYPE a, CELL_INDEX_TYPE b) {
//		float av = cell_value(a);
//		float bv = cell_value(b);
//		if (av < bv) return true;
//		if (bv < av) return false;
//		return a < b;
//	}
//
//};





class mscVectorNegatorGridFunc2d: public mscVectorToScalarGridFunc2D
{

public:
	virtual dvect2d getVector(CELL_INDEX_TYPE cid) {
		dvect2d tmp = this->vectors[g2d(cid)];
		tmp.x *= -1.0;
		tmp.y *= -1.0;
		return tmp;
	}

		mscVectorNegatorGridFunc2d(double* vin, int x, int y, mscBasicMeshHandler* mh) :
			mscVectorToScalarGridFunc2D(vin, x, y, mh) {
			}


};




//
//template<typename dtype>
//class mscRegularRawDataHandler : public mscBasicDataHandler<dtype> {
//protected:
//	mscBasicArray<dtype>* values;
//
//	public: 
//
//	  mscRegularRawDataHandler() {}
//   virtual ~mscRegularRawDataHandler() {
//		  		printf("delete: mscRegularRawDataHandler \n");
//
//		 delete values;
//	  }
//
//	  bool load_data(char* filename, CELL_INDEX_TYPE num_elements, mscArrayFactory* array_factory) {
//	
//		  FILE* fin = fopen(filename, "rb");
//		  if (fin == NULL) return false;
//
//		  values = array_factory->create<dtype>(num_elements);
//		  mscBasicArray<dtype>& values_r = *(values);
//		  
//		  for (CELL_INDEX_TYPE i = 0; i < num_elements; i++) {
//			  fread(&(values_r[i]), sizeof(dtype), 1, fin);
//		  }
//		  fclose(fin);
//		  return true;
//	  }
//
//	  void logify() {
//		  mscBasicArray<dtype>& values_r = *(values);
//		  CELL_INDEX_TYPE num_elements = values->size();
//		  for (CELL_INDEX_TYPE i = 0; i < num_elements; i++) {
//			  values_r[i] = log(values_r[i]);
//		  }
//	  }
//
//	  void hack_cut() {
//		  mscBasicArray<dtype>& values_r = *(values);
//		  CELL_INDEX_TYPE num_elements = values->size();
//		  for (CELL_INDEX_TYPE i = 0; i < num_elements; i++) {
//			  if (values_r[i] < .2) values_r[i] = .2;
//		  }
//	  }
//
//	  void negate() {
//		  mscBasicArray<dtype>& values_r = *(values);
//		  CELL_INDEX_TYPE num_elements = values->size();
//		  for (CELL_INDEX_TYPE i = 0; i < num_elements; i++) {
//			  values_r[i] = -1.0* (values_r[i]);
//		  }
//	  }
//      inline dtype value(CELL_INDEX_TYPE index) {
//		  mscBasicArray<dtype>& values_r = *values;
//		  return values_r[index];
//	  }
//
//
//
//	  void dump_vals(char* filename, int X, int Y, int Z, vector<idfpair>& v) {
//
//		  int dX = 2*X-1;
//		  int dY = 2*Y-1;
//		  int dZ = 2*Z-1;
//
//		  float* result = new float[X*Y*Z];
//
//		  for (int i = 0; i < X*Y*Z; i++) result[i] = 1.0f;
//
//		  for (int i =0; i < v.size(); i++) {
//
//			  CELL_INDEX_TYPE id = v[i].id;
//			  CELL_INDEX_TYPE x = id % dX;
//			  CELL_INDEX_TYPE y = (id / dX) % dY;
//			  CELL_INDEX_TYPE z = id / (dX*dY);
//
//			  //if (x%2 + y%2 + z%2 == 0) {
//				  result[(x/2)+(y/2)*X+(z/2)*Y*X] = min(result[(x/2)+(y/2)*X+(z/2)*Y*X], v[i].prob);
//			  //}
//		  }
//
//		  char newname[1024];
//		  sprintf(newname, "%s.prob", filename);
//		  FILE* fout = fopen(newname, "wb");
//		  fwrite(result, sizeof(float), X*Y*Z, fout);
//		fclose(fout);
//	  }
//
//};

#endif