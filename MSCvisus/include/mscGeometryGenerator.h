#ifndef MSC_GEOMETRY_GENERATOR
#define MSC_GEOMETRY_GENERATOR

#include "mscIndexTypes.h"
#include "mscGenericComplex.h"
#include "mscBasicMeshHandler.h"
#include "mscRegularGrid3DImplicitMeshHandler.h"
#include "mscRegularGrid2DImplicitMeshHandler.h"


#include <set>
#include <vector>
#include <map>

#ifndef WIN32
#include <cstdlib>
#endif

using namespace std;

const int onelink[26][4] =  {
/*0*/	{ 1, 0, 0, 0 },
/*1*/	{ -1, 0, 0, 0 },
/*2*/	{ 0, 1, 0, 0 },
/*3*/	{ 1, 1, 0, 0 },
/*4*/	{ -1, 1, 0, 0 },
/*5*/	{ 0, -1, 0, 0 },
/*6*/	{ 1, -1, 0, 0 },
/*7*/	{ -1, -1, 0, 0 },
/*8*/	{ 0, 0, 1, 0 },
/*9*/	{ 1, 0, 1, 0 },
/*0*/	{ -1, 0, 1, 0 },
/*1*/	{ 0, 1, 1, 0 },
/*2*/	{ 1, 1, 1, 0 },
/*3*/	{ -1, 1, 1, 0 },
/*4*/	{ 0, -1, 1, 0 },
/*5*/	{ 1, -1, 1, 0 },
/*6*/	{ -1, -1, 1, 0 },
/*7*/	{ 0, 0, -1, 0 },
/*8*/	{ 1, 0, -1, 0 },
/*9*/	{ -1, 0, -1, 0 },
/*0*/	{ 0, 1, -1, 0 },
/*1*/	{ 1, 1, -1, 0 },
/*2*/	{ -1, 1, -1, 0 },
/*3*/	{ 0, -1, -1, 0 },
/*4*/	{ 1, -1, -1, 0 },
/*5*/	{ -1, -1, -1, 0 }
};

class mscRegular3dGridSurfExtractor {
protected:
	mscRegularGrid3DImplicitMeshHandler* mMesh;
	
	
	inline void setOLneighbor(CELL_INDEX_TYPE* cell, CELL_INDEX_TYPE* neighb, int idx) {
		  
		  for (int i = 0; i < 3; i++) {
			  int v = ((int)cell[i]) + onelink[idx][i];
			  if (v < 0 || v >= mMesh->extent(i)) {
				neighb[i] = v - onelink[idx][i];
			  } else {
				neighb[i] = v;
			  }
		  }
	  }

	 int _av(CELL_INDEX_TYPE* coords2, int ox,mscGeomQuadSurface& s) {
		CELL_INDEX_TYPE negs[3];
		float fcoords[3];
		setOLneighbor(coords2, negs, ox); /// set neighbor coords in neg
		CELL_INDEX_TYPE id = mMesh->coords_2_cellid(negs);  /// get the index of neg
		for(int i = 0; i < 3; i++) {
			fcoords[i] = (float) negs[i];
		}
		return s.add_vertex(fcoords, id);
	  }

  bool edgeToQuad(CELL_INDEX_TYPE p, mscGeomQuadSurface& s) {
		CELL_INDEX_TYPE c2[3];
		mMesh->cellid_2_coords(p, c2);
		
		if (c2[0] % 2 == 1) {
			s.add_quad(_av(c2, 23,s),_av(c2, 20,s),_av(c2, 11,s),_av(c2, 14,s));
		} else if (c2[1] % 2 == 1) {
			s.add_quad(_av(c2, 19,s),_av(c2, 18,s),_av(c2, 9,s),_av(c2, 10,s));
		} else if (c2[2] % 2 == 1) {
			s.add_quad(_av(c2, 7,s),_av(c2, 6,s),_av(c2, 3,s),_av(c2, 4,s));
		} 
		return true;
	  }
  	 
  bool quadToQuad(CELL_INDEX_TYPE p, mscGeomQuadSurface& s) {
		CELL_INDEX_TYPE c2[3];
		mMesh->cellid_2_coords(p, c2);
		if (c2[0] % 2 == 0) {
			s.add_quad(_av(c2, 23,s),_av(c2, 20,s),_av(c2, 11,s),_av(c2, 14,s));
		} else if (c2[1] % 2 == 0) {
			s.add_quad(_av(c2, 19,s),_av(c2, 18,s),_av(c2, 9,s),_av(c2, 10,s));
		} else if (c2[2] % 2 == 0) {
			s.add_quad(_av(c2, 7,s),_av(c2, 6,s),_av(c2, 3,s),_av(c2, 4,s));
		} 
		return true;
	  }

public:
	mscRegular3dGridSurfExtractor (mscRegularGrid3DImplicitMeshHandler* mesh) :
		mMesh(mesh) {

	}
	// same as createsurface, but restrict do d-dimensional cells
	void CreateSurface(set<CELL_INDEX_TYPE>& ids, mscGeomQuadSurface& res, DIM_TYPE i) {

		if (ids.size() == 0) return;
		for (set<CELL_INDEX_TYPE>::iterator it = ids.begin(); it != ids.end(); it++) {
				CELL_INDEX_TYPE id = *it;
				
				DIM_TYPE d = mMesh->dimension(id);
				if (d == i && d == 1) {
					edgeToQuad(id, res);
				} else if (d == i && d == 2) {
					quadToQuad(id, res);
				}

		}


	}
	void CreateSurface(set<CELL_INDEX_TYPE>& ids, mscGeomQuadSurface& res) {
		// assume they are all the same dimension
		if (ids.size() == 0) return;
		for (set<CELL_INDEX_TYPE>::iterator it = ids.begin(); it != ids.end(); it++) {
				CELL_INDEX_TYPE id = *it;
				
				DIM_TYPE d = mMesh->dimension(id);
				if (d == 1) {
					edgeToQuad(id, res);
				} else if (d == 2) {
					quadToQuad(id, res);
				}

		}	
		
		//DIM_TYPE element_dimension = mMesh->dimension(*ids.begin());

		//if (element_dimension == 1) {
		//	for (set<CELL_INDEX_TYPE>::iterator it = ids.begin(); it != ids.end(); it++) {
		//		CELL_INDEX_TYPE id = *it;
		//		edgeToQuad(id, res);
		//	}
		//
		//} else if (element_dimension == 2) {
		//	for (set<CELL_INDEX_TYPE>::iterator it = ids.begin(); it != ids.end(); it++) {
		//		CELL_INDEX_TYPE id = *it;
		//		quadToQuad(id, res);
		//	}
		//}

	}




};



#endif
