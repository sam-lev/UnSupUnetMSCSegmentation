#ifndef MSC_GRAPH_INTERFACE
#define MSC_GRAPH_INTERFACE


//DEFINE A VALUATOR
typedef int IDTYPE;





class mscGraph {

public:
	mscGraph() {};

	//node iterator
	virtual IDTYPE FirstNode() {return 0;}
	virtual IDTYPE NextNode(IDTYPE id) {return 0;}
	virtual bool ValidNode(IDTYPE id) {return false;}

	//arc iterator
	virtual IDTYPE FirstArc() {return 0;}
	virtual IDTYPE NextArc(IDTYPE id) {return 0;}
	virtual bool ValidArc(IDTYPE id) {return false;}

	//node connection arc iterator
	virtual IDTYPE FirstNodeArc(IDTYPE nodeid) {return 0;}
	virtual IDTYPE NextNodeArc(IDTYPE nodeid, IDTYPE arcid) {return 0;}
	virtual bool ValidNodeArc(IDTYPE arcid) {return false;}
	
	//nodes connected to an arc
	virtual IDTYPE ArcNode(int position, IDTYPE arcid) {return 0;}

};


#endif
