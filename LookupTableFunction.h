#ifndef FunctionLookupTable_h
#define FunctionLookupTable_h

#include <vector>

template <typename FuncType>
class FunctionLookupTable  {
  public:
  typedef typename FuncType::value_type T;

  FunctionLookupTable(
    const T cubeSize,
    const size_t lutSize,
    const FuncType *func
  ) :
  cubeSize(cubeSize),
  lut(lutSize),
  lutResolution(
    ( (T) cubeSize) / 
    ( (T) (2 * (lutSize - 1) ) ) ),
  invLUTResolution(((T) 1.0) / lutResolution) {
    populateLUT(func); 
  }
  
  FunctionLookupTable(
    const T cubeSize,
    const size_t lutSize
  ) :
  cubeSize(cubeSize),
  lut(lutSize),
  lutResolution(
    ( (T) cubeSize) / 
    ( (T) (2 * (lutSize - 1) ) ) ),
  invLUTResolution(((T) 1.0) / lutResolution) {
    FuncType func(cubeSize);
    populateLUT(&func); 
  }


  public:
  T operator() (T z, T y, T x) const {
    return (*this)( std::sqrt(z*z + y*y + x*x) );   
  }

  T operator() (T radius) const {
    size_t lookupIndex = radius * invLUTResolution;

    if(lookupIndex >= lut.size()) {
      return 0; 
    }
    else {
      return lut[lookupIndex];
    }
  }

  protected:
  void populateLUT(const FuncType *func) {
    for(size_t i = 0; i < lut.size(); i++) {
      lut[i] = (*func)(i * lutResolution);
    }
  }

  protected:
  const T cubeSize;
  const T lutResolution;
  const T invLUTResolution;
  std::vector<T> lut;
};

#endif
