#ifndef MSC_BASIC_ARRAY
#define MSC_BASIC_ARRAY

template<class T>
class mscBasicArray {

   public:
      
	   virtual ~mscBasicArray() {
	   		printf("delete: mscBasicArray \n");

	   };
      virtual size_t size() = 0;
      inline virtual T &operator[](size_t index) = 0;
      inline virtual const T &operator[](size_t index) const  = 0;

};


#endif
