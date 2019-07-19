#ifndef MSC_TOPO_ARRAY
#define MSC_TOPO_ARRAY
#include <vector>


template<class T>
class mscTopoArray {
protected:
	vector<T> elements;

   public:
      
	   virtual ~mscTopoArray() {
	   		printf("delete: mscTopoArray \n");

	   };
	   virtual size_t size() { return elements.size(); }
	   inline virtual T operator[](size_t index) { return elements[index]; }
	   inline virtual bool set(T element, size_t index) {
		   if (index > elements.size()) return false;
			// USE LOCK HERE
		   elements[index] = element;
		   // RELEASE LOCK
		   return true;
	   }
	   inline virtual size_t append(T element) {
		   // vector does bad things since it reallocs, and not thread safe
		   elements.push_back(element);
		   return elements.size();
	   }
};


#endif