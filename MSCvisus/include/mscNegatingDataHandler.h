#ifndef MSC_NEGATING_DATA_HANDLER
#define MSC_NEGATING_DATA_HANDLER

#include "mscIndexTypes.h"
#include "mscBasicDataHandler.h"

template<typename dtype>
class mscNegatingDataHandler : public mscBasicDataHandler<dtype> {
protected:
	mscBasicDataHandler<dtype>* my_data_handler;
public: 
	   mscNegatingDataHandler(mscBasicDataHandler<dtype>* data_handler) :
		 my_data_handler(data_handler) {
		 }
	   virtual ~mscNegatingDataHandler() {
		printf("delete: mscNegatingDataHandler \n");
	   }
      inline virtual dtype value(CELL_INDEX_TYPE index) {
		  return -1 * my_data_handler->value(index);
	  }

};

#endif
