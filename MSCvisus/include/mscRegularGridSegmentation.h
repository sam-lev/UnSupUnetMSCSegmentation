#ifndef MSC_SIMPLE_RANDOM_GRADIENT_BUILDER
#define MSC_SIMPLE_RANDOM_GRADIENT_BUILDER

#include "mscSimpleGradientBuilder.h"

#ifndef WIN32
#include <cstdlib>
#endif

using namespace std;

template<typename dtype>
class mscSimpleRandomGradientBuilder : public mscSimpleGradientBuilder<dtype> {
protected:

	// return index in candidates of pair
	virtual int pick_from_candidates(const CELL_INDEX_TYPE& cellid, 
		const vector<CELL_INDEX_TYPE>& candidates) {
		
		dtype cell_value = this->my_mesh_function->cell_value(cellid);

		int result = 0; 
		vector<dtype> values;
		dtype temp_sum = 0;

		//printf("%d=%.2f\n", candidates[0], minval);
		for (int i = 0; i < candidates.size(); i++) {
			dtype diff_val = cell_value - this->lowest_facet_value(candidates[i]);
			temp_sum += diff_val;
			values.push_back(temp_sum);
		}
		float random = ((rand() % 1000) * 0.001f) * temp_sum;
		//printf("tempsum:%f, random=%f\n",temp_sum, random);
		if (temp_sum == 0.0f) {
			return rand() % candidates.size();
		}

		for (int i = 0; i < candidates.size(); i++) {
			if (values[i] > random) return i;
		}
		return candidates.size() - 1;
	}


public:

	void set_random_seed(int number) {
		srand(number);
	}

	mscSimpleRandomGradientBuilder(
		mscBasicMeshFunction<dtype>* mesh_function,
		mscBasicMeshHandler* mesh_handler,
		mscBasicGradientField* grad_field,
		mscArrayFactory* array_factory) : 
			mscSimpleGradientBuilder<dtype>(mesh_function, mesh_handler, grad_field, array_factory) {
				
	}

};


#endif
