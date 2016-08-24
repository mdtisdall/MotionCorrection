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
  func(func),
  lut(lutSize),
  lutResolution(
    ( (T) cubeSize) / 
    ( (T) (2 * (lutSize - 1) ) ) ),
  invLUTResolution(((T) 1.0) / lutResolution) {
    populateLUT(); 
  }


  public:
  T operator() (T radius) const {
    return lut[(size_t) (radius * invLUTResolution)];
  }

  protected:
  void populateLUT() {
    for(size_t i = 0; i < lut.size(); i++) {
      lut[i] = (*func)(i * lutResolution);
    }
  }

  protected:
  const T cubeSize;
  const T lutResolution;
  const T invLUTResolution;
  const FuncType *func;
  std::vector<T> lut;
};

#endif
