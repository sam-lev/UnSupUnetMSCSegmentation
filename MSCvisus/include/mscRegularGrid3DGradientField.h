#ifndef MSC_REGULAR_GRID_3D_GRADIENT_FIELD
#define MSC_REGULAR_GRID_3D_GRADIENT_FIELD

#include <stdio.h>
#include "mscArrayFactory.h"
#include "mscRegularGrid3DImplicitMeshHandler.h"
#include "mscBasicGradientField.h"
//#include "mscRegularGrid3DDecomposition.h"


class mscRegularGrid3DGradientField : public mscBasicGradientField {
protected:
	struct bitfield {
		unsigned char assigned : 1;
		unsigned char flag : 1;
		//unsigned char critical : 1;
		//unsigned char insorter : 1;
		//unsigned char dimA : 3;
		unsigned char pair : 3;
		unsigned char ldir : 3;
		//unsigned char empty_flag : 1;
	}; 
   
	mscBasicArray<bitfield>* cells;
	mscRegularGridImplicitMeshHandler* my_mesh_handler;
	

   public:

	   void resetMeshHandler(mscRegularGridImplicitMeshHandler* basic_mesh_handler) {
		   my_mesh_handler = basic_mesh_handler;
	   }
	   mscRegularGrid3DGradientField(mscRegularGridImplicitMeshHandler* basic_mesh_handler, 
		   mscArrayFactory* array_factory=NULL) {
			   if (basic_mesh_handler == NULL) {
				   printf("Error: mscRegularGrid3dGradientField constructor, basicMeshFunction == NULL\n");
				   return;
			   }
			   my_mesh_handler = basic_mesh_handler;
			   if (array_factory == NULL) array_factory = new mscArrayFactory(REGULAR_ARRAY);

			   cells = array_factory->create<bitfield>(basic_mesh_handler->num_cells());
	   }

	 virtual   ~mscRegularGrid3DGradientField() {
		   printf("delete: mscRegularGrid3DGradientField \n");

		   delete cells;
	   }
	   // get set state in gradient.
	   // note the assumptions on when each call is valid!

	   // ALWAYS VALID 
		   virtual ASSIGNED_TYPE get_assigned(CELL_INDEX_TYPE cellid) {
		   return (*cells)[cellid].assigned;
	   }

	   virtual void set_assigned(CELL_INDEX_TYPE cellid, ASSIGNED_TYPE value) {
		   (*cells)[cellid].assigned = value;
	   }

	   // VALID ONLY AFTER FIRST ASSIGNMENT
	   virtual CELL_INDEX_TYPE get_pair(CELL_INDEX_TYPE cellid){
		   return cellid + my_mesh_handler->offset_2_pair((*cells)[cellid].pair);
	   }

	   virtual void set_pair(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE value) {
		   (*cells)[cellid].pair = my_mesh_handler->pair_2_offset(value - cellid);
	   }

	   virtual bool get_critical(CELL_INDEX_TYPE cellid)  {
		   return (*cells)[cellid].pair == 7;
	   }
	   virtual void set_critical(CELL_INDEX_TYPE cellid, bool value) {
		   if (value) (*cells)[cellid].pair = 7;
	   }

	   virtual DIM_TYPE get_dim_asc_man(CELL_INDEX_TYPE cellid){
		   return (*cells)[cellid].ldir;
	   }
	   virtual void set_dim_asc_man(CELL_INDEX_TYPE cellid, DIM_TYPE value){
		   (*cells)[cellid].ldir = value;
	   }

	   // VALID ONLY BEFORE FIRST ASSIGNEMT
	   virtual CELL_INDEX_TYPE get_num_unpaired_facets(CELL_INDEX_TYPE cellid) {
		   return (*cells)[cellid].ldir;
	   }
	   virtual void set_num_unpaired_facets(CELL_INDEX_TYPE cellid, CELL_INDEX_TYPE value) {
		   (*cells)[cellid].ldir = value;
	   }

	   virtual DIM_TYPE get_mark(CELL_INDEX_TYPE cellid){
		   return (*cells)[cellid].pair;
	   }
	   virtual void set_mark(CELL_INDEX_TYPE cellid, DIM_TYPE value){
		   (*cells)[cellid].pair = value;
	   }
	   
  

	 void output_to_renderer(const char* filenamebase) {
		
		char grad_name[1024];
		sprintf(grad_name, "%s.grad", filenamebase);

		FILE* fgrad = fopen(grad_name, "wb");

		for (int i = 0; i < my_mesh_handler->num_cells(); i++) {

			(*cells)[i].flag = (bool) my_mesh_handler->boundary_value(i);
			//if (i % 1000 == 0)
			//	printf("cell[%d]:(%d, %d, %d, %d)\n", i,
			//	(*cells)[i].pair,
			//	(*cells)[i].ldir,
			//	(*cells)[i].flag,
			//	(*cells)[i].assigned);

		
			fwrite(&(*cells)[i], sizeof(bitfield), 1, fgrad);		
		}
		fclose(fgrad);
	}

	 void output_to_file(const char* filenamebase) {
		

		FILE* fgrad = fopen(filenamebase, "wb");

		for (int i = 0; i < my_mesh_handler->num_cells(); i++) {

			(*cells)[i].flag = (bool) my_mesh_handler->boundary_value(i);

		
			fwrite(&(*cells)[i], sizeof(bitfield), 1, fgrad);		
		}
		fclose(fgrad);
	}


	 bool load_from_file(const char* filename) {
	
		
		mscBasicArray<bitfield>& cells_r = *(cells);

		FILE* fdat = fopen(filename, "rb");
		if (fdat == NULL) {
			return false;

		}

		for (int i = 0; i < my_mesh_handler->num_cells(); i++) {
			
			bitfield val;
			fread(&val, sizeof(bitfield), 1, fdat);
			cells_r[i]=val;

		}

		fclose(fdat);
		return true;
	}
};

#endif
