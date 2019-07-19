#ifndef MSC_DUMB_STORING_MIN_FUNCTION
#define MSC_DUMB_STORING_MIN_FUNCTION

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicMeshFunction.h"
#include "mscArrayFactory.h"

template<typename dtype>
class mscDumbStoringMinFunction : public mscBasicMeshFunction<dtype> {
protected:
   mscBasicDataHandler<dtype>* my_data_handler;
   mscBasicMeshHandler* my_mesh_handler;
   mscArrayFactory* my_array_factory;

   mscBasicArray<dtype>* my_values;

   // assume all facet values are correct!!!
   dtype min_of_facets(CELL_INDEX_TYPE cellid) {
	   
	   dtype result;
	   bool has_result = false;

	   cellIterator it;
	   iteratorOperator& ito = my_mesh_handler->facets(cellid, it);
	   ito.begin(it);
	   
	   while (ito.valid(it)) {
		   if (! has_result) {
			   result = cell_value(ito.value(it));
			   has_result = true;
		   } else {
			   dtype temp_other = cell_value(ito.value(it)); 
			   result = (result < temp_other? result : temp_other);
		   }
		   ito.advance(it);
	   }
	   return result;
   }
      // assume all facet values are correct!!!
   dtype ave_of_facets(CELL_INDEX_TYPE cellid) {
	   
	   dtype result = 0;
	   dtype count = 0;
	  
	   cellIterator it;
	   iteratorOperator& ito = my_mesh_handler->facets(cellid, it);
	   ito.begin(it);
	   
	   while (ito.valid(it)) {
		   result += cell_value(ito.value(it)); 
		   count += 1;
	   	   ito.advance(it);
	   }
	   return result / count;
   }
   // first do 0's. Assume that the i'th 0-cell in the mesh is the i'th vertex in data
   void set_vertex_values() {
	   mscBasicArray<dtype>& values_r = *(my_values);

	   cellIterator it;
	   iteratorOperator& ito = my_mesh_handler->d_cells_iterator(0, it);
	   ito.begin(it);
	   CELL_INDEX_TYPE temp_grad_2_data = 0;
	   
	   while (ito.valid(it)) {
		   values_r[ito.value(it)] = my_data_handler->value(temp_grad_2_data);
		   temp_grad_2_data++;
		   ito.advance(it);
	   }  
   }

   void set_cell_values(DIM_TYPE dim) {
	   mscBasicArray<dtype>& values_r = *(my_values);

	   cellIterator it;
	   iteratorOperator& ito = my_mesh_handler->d_cells_iterator(dim, it);
	   ito.begin(it);
	   
	   while (ito.valid(it)) {
		   CELL_INDEX_TYPE temp_id = ito.value(it);
		   values_r[temp_id] = min_of_facets(temp_id);
		   ito.advance(it);
	   }  
   }

public:
	mscDumbStoringMinFunction(mscBasicDataHandler<dtype>* data_handler,
		mscBasicMeshHandler* mesh_handler,
		mscArrayFactory* array_factory) : my_data_handler(data_handler),
		my_mesh_handler(mesh_handler), my_array_factory(array_factory) {
			my_values = NULL;

	}
	virtual ~mscDumbStoringMinFunction() {
				printf("delete: mscDumbStoringMinFunction \n");

		delete my_values;
	}
	virtual void initialize() {
		 my_values = my_array_factory->create<dtype>(my_mesh_handler->num_cells());

		 // set the vertex values to data values
		 set_vertex_values();

		 // now set values for higher dim cells as max of facets
		 for (DIM_TYPE i = 1; i <= my_mesh_handler->max_dim(); i++) {
			 set_cell_values(i);
		 }
	}

	virtual dtype cell_value(CELL_INDEX_TYPE cellid) {
	   mscBasicArray<dtype>& values_r = *(my_values);
		return values_r[cellid];
	}



};

#endif