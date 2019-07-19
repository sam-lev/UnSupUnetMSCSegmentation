#ifndef MSC_BASIC_DATA_HANDLER
#define MSC_BASIC_DATA_HANDLER

#include "mscIndexTypes.h"

template<typename dtype>
class mscBasicDataHandler {
   public: 

	   virtual ~mscBasicDataHandler() {
		printf("delete: mscBasicDataHandler \n");
	   }
      inline virtual dtype value(CELL_INDEX_TYPE index) = 0;

};

#endif
