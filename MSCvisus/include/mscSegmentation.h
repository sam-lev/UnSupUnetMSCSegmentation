#ifndef MSC_SEGMENTATION
#define MSC_SEGMENTATION

#include "mscIndexTypes.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMSC.h"
#include <set>
#include <map>

#ifndef WIN32
#include <cstdlib>
#endif

using namespace std;

class mscSegmentation {
protected:
public:
	mscSegmentation() {}

	virtual void Set(CELL_INDEX_TYPE at, CELL_INDEX_TYPE id) {};
	virtual CELL_INDEX_TYPE SegID(CELL_INDEX_TYPE cellid) {
		return 0;
	}


};


class mscDumbStoringSegmentation : public mscSegmentation {
protected:
	CELL_INDEX_TYPE* mSeg;
	CELL_INDEX_TYPE mNumIds;

public:

	mscDumbStoringSegmentation(mscBasicMeshHandler* mesh) {
		mSeg = new CELL_INDEX_TYPE[mesh->num_cells()];
		mNumIds = mesh->num_cells();
	}

	void Set(CELL_INDEX_TYPE at, CELL_INDEX_TYPE id) {
		mSeg[at] = id;
	}
	virtual CELL_INDEX_TYPE SegID(CELL_INDEX_TYPE cellid) {
		return mSeg[cellid];
	}
};

template<typename dtype>
class mscComputeBasinSegmentation {
protected:
	mscDumbStoringSegmentation* mSeg;
	BasicMSC<dtype>* mMSC;
	mscBasicMeshHandler* mMesh;

	struct SetLess {
		bool operator() (const set<CELL_INDEX_TYPE>& a, const set<CELL_INDEX_TYPE>& b) const {
			//printf("a.size=%d, b.size=%d\n", a.size(), b.size());
			if (a.size() < b.size()) return true;
			if (a.size() > b.size()) return false;

			set<CELL_INDEX_TYPE>::iterator ait = a.begin();
			set<CELL_INDEX_TYPE>::iterator bit = b.begin();

			while( ait != a.end()) {
				//printf("a=%d, b=%d\n", *ait, *bit);
				if (*ait != *bit && *ait < *bit) return true;
				if (*ait != *bit && *ait > *bit) return false;
				ait++; bit++;
			}
			return false;
		}
	};

	map< set<CELL_INDEX_TYPE>, CELL_INDEX_TYPE, SetLess> mIDMAP;

public:

	mscComputeBasinSegmentation(BasicMSC<dtype>* msc, mscBasicMeshHandler* mesh) :
		mMSC(msc), mMesh(mesh) {
			mSeg = new mscDumbStoringSegmentation(mesh);
			for (CELL_INDEX_TYPE i =0; i < mesh->num_cells(); i++) {
				mSeg->Set(i, -1);		
			}
	}

	mscSegmentation* GetSeg() {
		return this->mSeg;
	}


	void Compute() {
		// zero everything out
		//printf("gothere2\n");
		mIDMAP.clear();
		for (CELL_INDEX_TYPE i =0; i < mMesh->num_cells(); i++) {
				mSeg->Set(i, -1);
		}
		// first set the vertices
		for (typename map<CELL_INDEX_TYPE, node<dtype>* >::iterator it = mMSC->nodes.begin(); it != mMSC->nodes.end(); it++) {
			node<dtype>* n = (*it).second;
			if (mMSC->isAlive(n) && n->index == 0) {
				set<CELL_INDEX_TYPE> res;
				mMSC->fillGeometry(n, res);
				for (typename set<CELL_INDEX_TYPE>::iterator sit = res.begin(); sit != res.end(); sit++) {
					if (mMesh->dimension(*sit) == 0) {
						mSeg->Set(*sit, n->cellid);
					}
				}
			}
		}
		//printf("gothere3\n");

		// now higher dim stuff
		for (int d = 1; d <= mMesh->max_dim(); d++) {
			//printf("dim = %d\n", d);
			cellIterator it;
			iteratorOperator& dcells = mMesh->d_cells_iterator(d, it);
			for (dcells.begin(it); dcells.valid(it); dcells.advance(it)) {
				CELL_INDEX_TYPE id = dcells.value(it);

				set<CELL_INDEX_TYPE> ids;
				CELL_INDEX_TYPE ids2;
				cellIterator fit;
				iteratorOperator& facets = mMesh->facets(id, fit);
				for (facets.begin(fit); facets.valid(fit); facets.advance(fit)) {
					CELL_INDEX_TYPE f = facets.value(fit);
					CELL_INDEX_TYPE target = mSeg->SegID(f);
					DIM_TYPE tdim = mMesh->dimension(target);
					if (tdim == d - 1) {
						ids.insert(target);
					} else {
						ids2 = target;
					}
				}
				//printf("%d -> %d\n", id, ids.size());
				if (ids.size() == 0) {
					mSeg->Set(id, ids2); 
				} else if (ids.size() == 1) {
					mSeg->Set(id, *(ids.begin()));
				} else {
					//printf("a\n");
					if (mIDMAP.count(ids) == 0) {
						//printf("b\n");
						mIDMAP[ids] = id;
						//printf("c\n");
					}
					mSeg->Set(id, mIDMAP[ids]);
				}
			}
		}
	}
};

#endif
