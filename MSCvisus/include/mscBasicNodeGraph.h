#ifndef MSC_BASIC_NODE_GRAPH
#define MSC_BASIC_NODE_GRAPH

template<class T>
class mscBasicNodeGraph {

   public:
      
	   virtual ~mscBasicArray() {
	   		printf("delete: mscBasicArray \n");

	   };
      virtual size_t size() = 0;
      inline virtual T &operator[](size_t index) = 0;
      inline virtual const T &operator[](size_t index) const  = 0;

};


#endif
