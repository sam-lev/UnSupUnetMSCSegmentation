#ifndef MSC_ARRAY_FACTORY
#define MSC_ARRAY_FACTORY

#include "mscBasicArray.h"
#include "mscRegularArray.h"

enum MSC_ARRAY_KIND {
	REGULAR_ARRAY
};


class mscArrayFactory {
private:
	MSC_ARRAY_KIND my_array_kind;
public:

	mscArrayFactory(MSC_ARRAY_KIND kind=REGULAR_ARRAY) : my_array_kind(kind) {}

	void set_kind(MSC_ARRAY_KIND kind) {
		my_array_kind = kind;
	}

	template<class T>
	mscBasicArray<T>* create(size_t size) {
		switch(my_array_kind) {
		case REGULAR_ARRAY:
			return new mscRegularArray<T>(size);
			break;
		default:
			return new mscRegularArray<T>(size);
		}
	}

};


#endif