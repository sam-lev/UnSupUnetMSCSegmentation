#ifndef MSC_PRECLASSIFIER
#define MSC_PRECLASSIFIER

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#include "mscBasicGradientField.h"
#include "mscArrayFactory.h"
#include "mscSimpleGradientBuilder.h"

#include <vector>
#include <queue>

using namespace std;


class mscPreClassifier {
protected:
	int* mIDVol;
	mscBasicMeshHandler* mMesh;
public:

	mscPreClassifier(int* ascID, int* dscID, mscBasicMeshHandler* mesh) {
		
		mIDVol = new int[mesh->num_cells()];
	   cellIterator it;
	   iteratorOperator& it_all = mesh->all_cells_iterator(it);
	   it_all.begin(it);
	   while (it_all.valid(it)) {
		   mIDVol[it_all.value(it)] = -1;
		   it_all.advance(it);
	   }
	   iteratorOperator& it_ver = mesh->d_cells_iterator(0,it);
	   it_ver.begin(it);
	   CELL_INDEX_TYPE dataid = 0;
	   while (it_ver.valid(it)) {
		   mIDVol[it_ver.value(it)] = ascID[dataid];
		   dataid++;
		   it_ver.advance(it);
	   }
	   iteratorOperator& it_vox = mesh->d_cells_iterator(mesh->max_dim(),it);
	   it_vox.begin(it);
	   dataid = 0;
	   while (it_vox.valid(it)) {
		   mIDVol[it_vox.value(it)] = dscID[dataid];
		   dataid++;
		   it_vox.advance(it);
	   }
	   
	}

	
	int getId(const CELL_INDEX_TYPE& cellid) {
		return mIDVol[cellid];
	}

};




#endif
