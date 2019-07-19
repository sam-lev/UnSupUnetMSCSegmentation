#ifndef MSC_INDEXED_GRAPH
#define MSC_INDEXED_GRAPH

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#include "mscBasicGradientField.h"
#include "mscArrayFactory.h"

#include <vector>
#include <queue>
#include <map>
#include <set>




using namespace std;


template <class Element>
class mscIndexedGraph {


public:

map<CELL_INDEX_TYPE, Pair<INT_TYPE, Element> > vertices;
 










};








#endif