#ifndef MSC_GENERIC_COMPLEX
#define MSC_GENERIC_COMPLEX

#include "mscIndexTypes.h"
#include <set>
#include <vector>
#include <map>

#ifndef WIN32
#include <cstdlib>
#endif

using namespace std;

class mscGenericCell {
public:

	mscGenericCell(CELL_INDEX_TYPE n, DIM_TYPE d, set<CELL_INDEX_TYPE>& i) {
		name = n; dim = d; ids.insert(i.begin(), i.end());
	}
	mscGenericCell(CELL_INDEX_TYPE n, DIM_TYPE d) {
		name = n; dim = d; 
	}
	~mscGenericCell() {
		neighbors.clear();
		ids.clear();
	}

	CELL_INDEX_TYPE name;
	DIM_TYPE dim;
	set<mscGenericCell*> neighbors;
	set<CELL_INDEX_TYPE> ids;	
};

class mscGenericComplex {
protected:

	map<CELL_INDEX_TYPE, mscGenericCell*> mCells[4]; // cells of each dimension


public:

	mscGenericComplex() {}
	
	map<CELL_INDEX_TYPE, mscGenericCell*>& GetCells(DIM_TYPE dim) {
		return mCells[dim];
	}

	~mscGenericComplex() {
		Clear();
	}

	void Clear() {
		for (int i = 0; i < 4; i++) {
			map<CELL_INDEX_TYPE, mscGenericCell*>::iterator it;
			for (it = mCells[i].begin(); it != mCells[i].end(); it++) 
				delete (*it).second;
			mCells[i].clear();
		}
	}

	bool InsertCell(CELL_INDEX_TYPE name, DIM_TYPE dim, set<CELL_INDEX_TYPE> ids) {
		if (mCells[dim].count(name) != 0) return false;
		mscGenericCell* c = new mscGenericCell(name, dim, ids);
		mCells[dim][name] = c;
		return true;
	}
	bool InsertCell(CELL_INDEX_TYPE name, DIM_TYPE dim) {
		if (mCells[dim].count(name) != 0) return false;
		mscGenericCell* c = new mscGenericCell(name, dim);
		mCells[dim][name] = c;
		return true;
	}
	bool InsertCell(mscGenericCell* c) {
		if (mCells[c->dim].count(c->name) != 0) return false;
		mCells[c->dim][c->name] = c;
		return true;
	}

	mscGenericCell* GetCell(CELL_INDEX_TYPE name, DIM_TYPE dim) {
		return mCells[dim][name];
	}
	mscGenericCell* GetCell(CELL_INDEX_TYPE name) {
		for (int i =0; i < 4; i ++) {
			if (mCells[i].count(name) != 0) return mCells[i][name];
		}
		return NULL;
	}
	
	void ConnectCells(mscGenericCell* a, mscGenericCell* b) {
		a->neighbors.insert(b);
		b->neighbors.insert(a);
	}

};


#endif
