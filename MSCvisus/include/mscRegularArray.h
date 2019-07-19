#ifndef MSC_REGULAR_ARRAY
#define MSC_REGULAR_ARRAY

#include "mscBasicArray.h"

template<class T>
class mscRegularArray : public mscBasicArray<T> {
   
   protected:
      T* my_array;
      size_t my_size;
   public:

      mscRegularArray(size_t size) : my_size(size) {
         my_array = new T[size];
      }

      virtual ~mscRegularArray() {
		  		printf("delete: mscRegularArray \n");

         delete [] my_array;
      }

      size_t size() { return my_size; }

      inline T &operator[](size_t index) { return my_array[index]; }
      
      inline const T &operator[](size_t index) const { return my_array[index]; }

};

#endif
