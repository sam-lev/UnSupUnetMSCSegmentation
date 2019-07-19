#ifndef MSC_VALUATOR_INTERFACE
#define MSC_VALUATOR_INTERFACE



template<class RETURN_TYPE>
class mscValuatorR {

public:
	virtual RETURN_TYPE eval();
};

template<class RETURN_TYPE, class INPUT1_TYPE>
class mscValuatorRI : public mscValuatorR<RETURN_TYPE> 
{
protected:
	mscValuatorR* input1
public:
	virtual RETURN_TYPE eval(INPUT1_TYPE in1);
};



#endif
