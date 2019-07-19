#ifndef MSC_BASIC_ITERATOR
#define MSC_BASIC_ITERATOR


template <class T>
class mscBasicIterator {
private:
	T* val_p;
public:

	virtual bool valid() { return false; }
	virtual void operator++() {}
	virtual void operator+=(const SIZE_TYPE num) {
		SIZE_TYPE temp = num;
		while (temp > 0) this->operator++();
	}
	virtual T& operator*() { 
		return *val_p; }

};

class mscSimpleCellCountIterator : public mscBasicIterator<CELL_INDEX_TYPE> {

protected:

	CELL_INDEX_TYPE my_location;
	CELL_INDEX_TYPE my_total_size;

public:

	mscSimpleCellCountIterator(CELL_INDEX_TYPE total_size) 
		: my_total_size(total_size), my_location(0) {}

	bool valid() { 
		return my_location < my_total_size; 
	}

	virtual void operator++() {
		my_location++;
	}
	virtual void operator+=(const CELL_INDEX_TYPE num) {
		my_location+=num;
	}
	virtual CELL_INDEX_TYPE& operator*() {
		return my_location;
	}

};


#endif
