#ifndef DerivWeightFunction_h
#define DerivWeightFunction_h

#include <Eigen/Dense>

#include <cmath>

template <typename T>
class DerivWeightFunction {
  public:
    typedef T value_type;

    DerivWeightFunction(const T cubeSize) :
      invRadius(((T) 2.0) / ((T) cubeSize))
      {}

    T operator() (T z, T y, T x) const {
      return maskValue( std::sqrt(z*z + y*y + x*x) );   
    }
    
    T operator() (T radius) const {
        T r = radius * invRadius;

        if(r < (T) 0.75) {
          return ((T) 1.0); 
        }
        else if(r > (T) 1.0) {
          return ((T) 0.0); 
        }
        
        return wcos(r * (T) 2.0 - (T) 1.5);
    }

  protected:
    T wcos(T t) const {
      if(t < -((T) 0.5) || t > ((T) 0.5)) {
        return 0;
      }
      else {
        return - sin(t * (T) M_PI) * (T) M_PI; 
      }
    }

  protected:
    T invRadius;    
};

#endif
