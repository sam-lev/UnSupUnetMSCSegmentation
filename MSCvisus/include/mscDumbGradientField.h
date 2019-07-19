#ifndef MSC_DUMB_GRADIENT_FIELD
#define MSC_DUMB_GRADIENT_FIELD

#include "mscIndexTypes.h"
#include "mscBasicGradientField.h"

// THIS IS A DUMB CLASS FOR STORING THE GRADIENT, i.e. it is agnostic
// of any values. Uses 1+8+1+1+4 = approx 15 bytes per cell

class mscDumbGradientField : public mscBasicGradientField {

protected:
	// store everything in arrays in least efficient way possible
	mscBasicArray<ASSIGNED_TYPE>* dumb_assigned;
	mscBasicArray<CELL_INDEX_TYPE>* dumb_pair;
	mscBasicArray<bool>* dumb_critical;
	mscBasicArray<DIM_TYPE>* dumb_dim_asc_man;

	CELL_INDEX_TYPE dumb_num_cells;

   public:

	   mscDumbGradientField(mscBasicMeshHandler* basic_mesh_handler, 
		   mscArrayFactory* array_factory=NULL) {
			   if (basic_mesh_handler == NULL) {
				   printf("Error: mscDumbGradientField constructor, basicMeshFunction == NULL\n");
				   return;
			   }
			if (array_factory == NULL) array_factory = new mscArrayFactory(REGULAR_ARRAY);

			dumb_num_cells = basic_mesh_handler->num_cells();
			
			dumb_assigned = array_factory->create<ASSIGNED_TYPE>(dumb_num_cells);
			dumb_pair = array_factory->create<CELL_INDEX_TYPE>(dumb_num_cells);
			dumb_critical = array_factory->create<bool>(dumb_num_cells);
			dumb_dim_asc_man = array_factory->create<DIM_TYPE>(dumb_num_cells);

	   };
	   virtual ~mscDumbGradientField() {
		   		printf("delete: mscDumbGradientField \n");

		   delete dumb_assigned;
		   delete dumb_pair;
		   delete dumb_critical;
		   delete dumb_dim_asc_man;
	   }

	   // get set state in gradient.
	   // note the assumptions on when each call is valid!

	   // ALWAYS VALID 
	   virtual ASSIGNED_TYPE get_assigned(CELL_INDEX_TYPE cellid) {
		   return (*dumb_assigned)[cellid];
	   }

	   virtual void set_assigned(CELL_INDEX_TYPE cellid, ASSIGNED_TYPE value) {
		   (*dumb_assigned)[cellid] = value;
	   }

	   // VALID ONLY AFTER FIRST ASSIGNMENT
	   virtual CELL_INDEX_TYPE get_pair(CELL_INDEX_TYPE cellid){
		   return (*dumb_pair)[cellid];
	   }
	   virtual void set_pair(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE value) {
		   (*dumb_pair)[cellid] = value;
	   }

	   virtual bool get_critical(CELL_INDEX_TYPE cellid)  {
		   return (*dumb_critical)[cellid];
	   }
	   virtual void set_critical(CELL_INDEX_TYPE cellid, bool value) {
		   (*dumb_critical)[cellid] = value;
	   }

	   virtual DIM_TYPE get_dim_asc_man(CELL_INDEX_TYPE cellid){
		   return (*dumb_dim_asc_man)[cellid];
	   }
	   virtual void set_dim_asc_man(CELL_INDEX_TYPE cellid, DIM_TYPE value){
		   (*dumb_dim_asc_man)[cellid] = value;
	   }

	   // VALID ONLY BEFORE FIRST ASSIGNEMT
	   virtual CELL_INDEX_TYPE get_num_unpaired_facets(CELL_INDEX_TYPE cellid) {
		   return (*dumb_pair)[cellid];
	   }
	   virtual void set_num_unpaired_facets(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE value) {
		   (*dumb_pair)[cellid] = value;
	   }

	   virtual DIM_TYPE get_mark(CELL_INDEX_TYPE cellid){
		   return (*dumb_dim_asc_man)[cellid];
	   }
	   virtual void set_mark(CELL_INDEX_TYPE cellid, DIM_TYPE value){
		   (*dumb_dim_asc_man)[cellid] = value;
	   }
      
};

#endif
