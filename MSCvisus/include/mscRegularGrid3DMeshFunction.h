#ifndef MSC_REGULAR_GRID_3D_MESH_FUNCTION
#define MSC_REGULAR_GRID_3D_MESH_FUNCTION

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#include "mscArrayFactory.h"
#include "mscDumbStoringMeshFunction.h"
#include "mscRegularGrid3DImplicitMeshHandler.h"

template<typename dtype>
class mscRegularGrid3DMeshFunction : public mscDumbStoringMeshFunction<dtype> {

public:

	mscRegularGrid3DMeshFunction (mscBasicDataHandler<dtype>* data_handler,
		mscBasicMeshHandler* mesh_handler,
		mscArrayFactory* array_factory) : 
		mscDumbStoringMeshFunction<dtype>(data_handler, mesh_handler, array_factory) {

	}
		virtual ~mscRegularGrid3DMeshFunction() {}

	void output_to_renderer(const char* filenamebase) {
		
		char dat_name[1024];
		sprintf(dat_name, "%s.dat", filenamebase);

		FILE* fdat = fopen(dat_name, "wb");

		for (int i = 0; i < this->my_mesh_handler->num_cells(); i++) {
			
			dtype val = this->cell_value(i);
			fwrite(&val, sizeof(dtype), 1, fdat);
		}

		fclose(fdat);
	}

	void load_from_file(const char* filename) {
	
		this->my_values = this->my_array_factory->create(this->my_mesh_handler->num_cells());
		mscBasicArray<dtype>& values_r = *(this->my_values);

		FILE* fdat = fopen(filename, "rb");

		for (int i = 0; i < this->my_mesh_handler->num_cells(); i++) {
			
			dtype val;
			fread(&val, sizeof(dtype), 1, fdat);
			values_r[i]=val;

		}

		fclose(fdat);
	}

};

#endif
