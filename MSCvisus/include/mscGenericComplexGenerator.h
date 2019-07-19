#ifndef MSC_GENERIC_COMPLEX_GENERATOR
#define MSC_GENERIC_COMPLEX_GENERATOR

#include "mscGenericComplex.h"
#include "mscSegmentation.h"
#include "mscBasicMeshHandler.h"
#include <map>

#ifndef WIN32
#include <cstdlib>
#endif

using namespace std;

class mscGenericComplexGenerator {
protected:
	mscBasicMeshHandler* mMesh;
	mscSegmentation* mSeg;
	mscGenericComplex* mCom;

	map<CELL_INDEX_TYPE, mscGenericCell*> mCellMap;

public:

	mscGenericComplexGenerator(mscBasicMeshHandler* m, mscSegmentation* s) :
		mMesh(m), mSeg(s) {
		mCom = new mscGenericComplex();

	}

	mscGenericComplex* Complex() {
		return mCom;
	}

	void Generate() {
		mCom->Clear();
		
		// now go through and populate map
		cellIterator it;
		iteratorOperator& allcells = mMesh->all_cells_iterator(it);
		for (allcells.begin(it); allcells.valid(it); allcells.advance(it)){
			CELL_INDEX_TYPE id = allcells.value(it);
			CELL_INDEX_TYPE target = mSeg->SegID(id);
			DIM_TYPE d = mMesh->dimension(target);

			if (mCom->GetCells(d).count(target) == 0) {
				mCom->InsertCell(target, d);
			}
			// add the current cell to the big cell
			mCom->GetCell(target, d)->ids.insert(id);
		}
		

		for (int i =1; i <4; i++) {
			map<CELL_INDEX_TYPE, mscGenericCell*>& cells = mCom->GetCells(i);

			for (map<CELL_INDEX_TYPE, mscGenericCell*>::iterator it = cells.begin();
				it != cells.end(); it++) {
					mscGenericCell* m = (*it).second;

					cellIterator cit;
					iteratorOperator& facets = mMesh->facets(m->name, cit);
					for (facets.begin(cit); facets.valid(cit); facets.advance(cit)) {
						CELL_INDEX_TYPE target = mSeg->SegID(facets.value(cit));
						//if (mMesh->dimension(target) != m->dim-1) printf("WHREOSDLFKJSDL:FKJSD:LFKJS\n");
						mscGenericCell* other = mCom->GetCell(target);
						mCom->ConnectCells(m, other);
					}

			}


		}

		printf("done generating new complex\n");

	}
};


#endif
