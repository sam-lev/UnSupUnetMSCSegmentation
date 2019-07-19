#ifndef MSC_COMBINATION_GRADIENT_FIELD
#define MSC_COMBINATION_GRADIENT_FIELD

#include "mscBasicGradientField.h"

class mscCombinationGradientField : public mscBasicGradientField {
    
	mscBasicGradientField* f[2];
	int my_assigned, my_pair, my_critical, my_dim_asc_man, my_mark, my_num_unpaired_facets;

   public:

	   mscCombinationGradientField(mscBasicGradientField* f1,
		   mscBasicGradientField* f2,
		   int assigned, int pair,
		   int critical, int dim_asc_man,
		   int num_unpaired_facets,
		   int mark) {
			   my_assigned = assigned;
			   my_pair = pair;
			   my_critical = critical;
			   my_dim_asc_man = dim_asc_man;
			   my_mark = mark;
			   my_num_unpaired_facets = num_unpaired_facets;

	   }

	   virtual ~mscCombinationGradientField() {		
		   printf("delete: mscComboGradientField \n");
		}
	   // constructor takes a pointer to an mscBasicMeshFunction

	   //mscBasicGradientField(mscBasicMeshFunction<dtype>* basicMeshFunction) {};

	   // get set state in gradient.
	   // note the assumptions on when each call is valid!

	   // ALWAYS VALID 
	   virtual ASSIGNED_TYPE get_assigned(CELL_INDEX_TYPE cellid){
		   return f[my_assigned]->get_assigned(cellid);
	   }

	   virtual void set_assigned(CELL_INDEX_TYPE cellid, ASSIGNED_TYPE value) {
		   f[my_assigned]->set_assigned(cellid, value);
	   }

	   // VALID ONLY AFTER FIRST ASSIGNMENT
	   virtual CELL_INDEX_TYPE get_pair(CELL_INDEX_TYPE cellid) {
		   return f[my_pair]->get_pair(cellid);
	   }

	   virtual void set_pair(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE value) {
		   f[my_pair]->set_assigned(cellid, value);
	   }

	   virtual bool get_critical(CELL_INDEX_TYPE cellid) {
		   return f[my_critical]->get_critical(cellid);
	   }
	   
	   virtual void set_critical(CELL_INDEX_TYPE cellid, bool value) {
		   f[my_critical]->set_critical(cellid, value);
	   }

	   virtual DIM_TYPE get_dim_asc_man(CELL_INDEX_TYPE cellid) {
		   return f[my_dim_asc_man]->get_dim_asc_man(cellid);
	   }

	   virtual void set_dim_asc_man(CELL_INDEX_TYPE cellid, DIM_TYPE value) {
		   f[my_dim_asc_man]->set_assigned(cellid, value);
	   }


	   // VALID ONLY BEFORE FIRST ASSIGNEMT
	   virtual CELL_INDEX_TYPE get_num_unpaired_facets(CELL_INDEX_TYPE cellid) {
		   return f[my_num_unpaired_facets]->get_num_unpaired_facets(cellid);
	   }
	   
	   virtual void set_num_unpaired_facets(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE value){
		   f[my_num_unpaired_facets]->set_num_unpaired_facets(cellid, value);
	   }

	   virtual DIM_TYPE get_mark(CELL_INDEX_TYPE cellid) {
		   return f[my_mark]->get_mark(cellid);
	   }

	   virtual void set_mark(CELL_INDEX_TYPE cellid, DIM_TYPE value) {
		   f[my_mark]->set_mark(cellid, value);   
	   }

	   
      
};

#endif
