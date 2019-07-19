#ifndef MSC_SIMPLE_POINTER_DATA_HANDLER
#define MSC_SIMPLE_POINTER_DATA_HANDLER

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"

#include <math.h>

template<typename dtype>
class mscSimplePointerDataHandler : public mscBasicDataHandler<dtype> {
protected:
	dtype* mValues;

	public: 

	  mscSimplePointerDataHandler() {}
	  virtual ~mscSimplePointerDataHandler() {
		  		printf("delete: mscSimplePointerDataHandler \n");

		delete[] mValues;
	  }

	  bool set_data(dtype* values) {
		  mValues = values;
		  return true;
	  }

      inline dtype value(CELL_INDEX_TYPE index) {
		  return mValues[index];
	  }



};

#endif