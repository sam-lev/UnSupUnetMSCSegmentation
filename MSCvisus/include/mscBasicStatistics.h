#ifndef MSC_BASIC_STATISTICS
#define MSC_BASIC_STATISTICS

#include <vector>
#include <numeric>

#undef min
#undef max

namespace mscStatistics {

template <typename dtype> 
dtype mean(vector<dtype>& v) {
	dtype sum = std::accumulate(std::begin(v), std::end(v), 0.0);
	double m =  sum / v.size();
	return m;
};

template <typename dtype> 
dtype variance(vector<dtype>& v, dtype m) {
	dtype accum = 0.0;
	std::for_each (std::begin(v), std::end(v), [&](const dtype d) {
	  accum += (d - m) * (d - m);
	});

	return (accum / (v.size()-1));
};

template <typename dtype>
dtype min(vector<dtype>& v) {
	dtype mv = v[0];
	for (int i = 1; i < v.size(); i++) {
		if (v[i] < mv) mv = v[i];
	}
	return mv;
};

template <typename dtype>
dtype max(vector<dtype>& v) {
	dtype mv = v[0];
	for (int i = 1; i < v.size(); i++) {
		if (v[i] > mv) mv = v[i];
	}
	return mv;
};

};

#endif