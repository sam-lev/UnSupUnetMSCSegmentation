#ifndef MSC_REGULAR_RAW_DATA_HANDLER
#define MSC_REGULAR_RAW_DATA_HANDLER

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"
#include "mscArrayFactory.h"
#include "mscConvergentGradientBuilder.h"


#include <vector>
#include <math.h>

using namespace std;

template<typename dtype>
class mscRegularRawDataHandler : public mscBasicDataHandler<dtype> {
protected:
	mscBasicArray<dtype>* values;

	public: 

	  mscRegularRawDataHandler() {}
   virtual ~mscRegularRawDataHandler() {
		  		printf("delete: mscRegularRawDataHandler \n");

		 delete values;
	  }

	  bool load_data(char* filename, CELL_INDEX_TYPE num_elements, mscArrayFactory* array_factory) {
	
		  FILE* fin = fopen(filename, "rb");
		  if (fin == NULL) return false;

		  values = array_factory->create<dtype>(num_elements);
		  mscBasicArray<dtype>& values_r = *(values);
		  
		  for (CELL_INDEX_TYPE i = 0; i < num_elements; i++) {
			  fread(&(values_r[i]), sizeof(dtype), 1, fin);
		  }
		  fclose(fin);




		  return true;
	  }

	  void logify() {
		  mscBasicArray<dtype>& values_r = *(values);
		  CELL_INDEX_TYPE num_elements = values->size();
		  for (CELL_INDEX_TYPE i = 0; i < num_elements; i++) {
			  values_r[i] = log(values_r[i]);
		  }
	  }

	  void hack_cut() {
		  mscBasicArray<dtype>& values_r = *(values);
		  CELL_INDEX_TYPE num_elements = values->size();
		  for (CELL_INDEX_TYPE i = 0; i < num_elements; i++) {
			  if (values_r[i] < .2) values_r[i] = .2;
		  }
	  }

	  void negate() {
		  mscBasicArray<dtype>& values_r = *(values);
		  CELL_INDEX_TYPE num_elements = values->size();
		  for (CELL_INDEX_TYPE i = 0; i < num_elements; i++) {
			  values_r[i] = -1.0* (values_r[i]);
		  }
	  }
      inline dtype value(CELL_INDEX_TYPE index) {
		  mscBasicArray<dtype>& values_r = *values;
		  return values_r[index];
	  }

	  //void testfart() {
		 // for (int i =0; i < 5; i++) 
			//  printf("%d = %.4f\n", i, value(i + i*100));
	  //}

	  void dump_vals(char* filename, int X, int Y, int Z, vector<idfpair>& v) {

		  int dX = 2*X-1;
		  int dY = 2*Y-1;
		  int dZ = 2*Z-1;

		  float* result = new float[X*Y*Z];

		  for (int i = 0; i < X*Y*Z; i++) result[i] = 1.0f;

		  for (int i =0; i < v.size(); i++) {

			  CELL_INDEX_TYPE id = v[i].id;
			  CELL_INDEX_TYPE x = id % dX;
			  CELL_INDEX_TYPE y = (id / dX) % dY;
			  CELL_INDEX_TYPE z = id / (dX*dY);

			  //if (x%2 + y%2 + z%2 == 0) {
				  result[(x/2)+(y/2)*X+(z/2)*Y*X] = min(result[(x/2)+(y/2)*X+(z/2)*Y*X], v[i].prob);
			  //}
		  }

		  char newname[1024];
		  sprintf(newname, "%s.prob", filename);
		  FILE* fout = fopen(newname, "wb");
		  fwrite(result, sizeof(float), X*Y*Z, fout);
		fclose(fout);
	  }

};

#endif