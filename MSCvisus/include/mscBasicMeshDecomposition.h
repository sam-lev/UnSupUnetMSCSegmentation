#ifndef MSC_BASIC_MESH_DECOMPOSITION
#define MSC_BASIC_MESH_DECOMPOSITION

#include "mscIndexTypes.h"
#include "mscBasicMeshFunction.h"
#include "mscBasicDataHandler.h"
#include "mscBasicMeshHandler.h"
#include "mscBasicGradientField.h"


template<typename dtype>
class mscBasicMeshDecomposition {

   public:

	   virtual ~mscBasicMeshDecomposition(){
	   		printf("delete: mscBasicDataHandler \n");
}
	   virtual void decompose() = 0;
	   virtual int numBlocks() = 0;
	   
	   virtual void loadBlock(INT_TYPE block_id) = 0;
	   virtual void unloadBlock(INT_TYPE block_id) = 0;

	   virtual mscBasicMeshFunction<dtype>* getMeshFunction(INT_TYPE block_id) = 0;
	   virtual mscBasicMeshHandler* getMeshHandler(INT_TYPE block_id) = 0;
	   virtual mscBasicGradientField* getGradientField(INT_TYPE block_id) = 0;

	   virtual void writeToDisk(char* filename) = 0;

};


#endif
