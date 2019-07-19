#ifndef FILTER_EXPRESSION
#define FILTER_EXPRESSION


enum ExpressionType {
	INT8, INT16, INT32, INT64,
	UINT8, UINT16, UINT32, UINT64,
	FLOAT, DOUBLE,
	BOOL,
	GRAPH, SET, MESH
};



template<class T>
class Expression {
protected:
	
	vector< Expression<T>* > operands;

   public:

	void addOperand(Expression<T>* e) { operands.push_back(e); }
    virtual T evaluate() ;

};






/////////////////// ARITHMETIC EXPRESSIONS

template<typename T>
class Sum : public Expression<T> {
   public:
	   T evaluate() {
		  T result = 0;
		  for (int i = 0; i < this->operands.size(); i++) {
			  T = T + this->operands[i]->evaluate();
		  }
		  return T;
	  }
};

template<typename T>
class Product : public Expression<T>{

   public:
		  
	   T evaluate() {
		  T result = 1;
		  for (int i = 0; i < this->operands.size(); i++) {
			  T = T * this->operands[i]->evaluate();
		  }
		  return T;
	  }
};

template<typename T>
class Divide : public Expression<T>{

   public:
	   
	   T evaluate() {
		  return this->operands[0]->evaluate() / this->operands[1]->evaluate();
	  }
};


template<typename T>
class Negate : public Expression<T>{

   public:
	   
	   T evaluate() {
		  return -this->operands[0]->evaluate();
	  }
};


/////////////////// BOOLEAN EXPRESSIONS


class Not : public Expression<bool>{

   public:
	   
	   bool evaluate() {
		  return -this->operands[0]->evaluate();
	  }
};

class Or : public Expression<bool>{

   public:
	   
	   bool evaluate() {
		   bool result = false;
		   for (int i = 0; i < this->operands.size(); i++) {
			   result = result || this->operands[0].evaluate();
		   }
		  return result;
	  }
};


class And : public Expression<bool>{

   public:
	   
	   bool evaluate() {
		   bool result = false;
		   for (int i = 0; i < this->operands.size(); i++) {
			   result = result &&| this->operands[0].evaluate();
		   }
		  return result;
	  }
};

















#endif
