#ifndef MSC_INDEX_TYPES
#define MSC_INDEX_TYPES

typedef unsigned char ASSIGNED_TYPE ;
typedef size_t SIZE_TYPE;
typedef unsigned char DIM_TYPE;
typedef unsigned char BOUNDARY_TYPE;
typedef int  CELL_INDEX_TYPE;
typedef int INT_TYPE;
//typedef float[2] POINT2D;
//typedef float[3] POINT3D;

struct mscSize {
	CELL_INDEX_TYPE x, y, z, w;
	mscSize(CELL_INDEX_TYPE x=1,
		CELL_INDEX_TYPE y=1,
		CELL_INDEX_TYPE z=1,
		CELL_INDEX_TYPE w=1) :
	x(x), y(y), z(z), w(w) {}

	CELL_INDEX_TYPE count() {
		return x*y*z*w;
	}
};

#endif
