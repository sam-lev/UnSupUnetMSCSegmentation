#ifndef MSC_BASIC_GRADIENT_FIELD
#define MSC_BASIC_GRADIENT_FIELD

#include "mscIndexTypes.h"
#include "mscBasicMeshHandler.h"

class mscBasicGradientField {
    
   public:

	   virtual ~mscBasicGradientField() {		
		   printf("delete: mscBasicGradientField \n");
	}
	   // constructor takes a pointer to an mscBasicMeshFunction

	   //mscBasicGradientField(mscBasicMeshFunction<dtype>* basicMeshFunction) {};

	   // get set state in gradient.
	   // note the assumptions on when each call is valid!

	   // ALWAYS VALID 
	   virtual ASSIGNED_TYPE get_assigned(CELL_INDEX_TYPE cellid) = 0;
	   virtual void set_assigned(CELL_INDEX_TYPE cellid, ASSIGNED_TYPE value) = 0;

	   // VALID ONLY AFTER FIRST ASSIGNMENT
	   virtual CELL_INDEX_TYPE get_pair(CELL_INDEX_TYPE cellid) = 0;
	   virtual void set_pair(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE value)  = 0;

	   virtual bool get_critical(CELL_INDEX_TYPE cellid)  = 0;
	   virtual void set_critical(CELL_INDEX_TYPE cellid, bool value) = 0;

	   virtual DIM_TYPE get_dim_asc_man(CELL_INDEX_TYPE cellid) = 0;
	   virtual void set_dim_asc_man(CELL_INDEX_TYPE cellid, DIM_TYPE value) = 0;

	   // VALID ONLY BEFORE FIRST ASSIGNEMT
	   virtual CELL_INDEX_TYPE get_num_unpaired_facets(CELL_INDEX_TYPE cellid) = 0;
	   virtual void set_num_unpaired_facets(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE value) = 0;

	   virtual DIM_TYPE get_mark(CELL_INDEX_TYPE cellid) = 0;
	   virtual void set_mark(CELL_INDEX_TYPE cellid, DIM_TYPE value) = 0;

	   
      
};

#endif
