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
      return (*this)( std::sqrt(z*z + y*y + x*x) );   
    }
    
    T operator() (T radius) const {
        T r = radius * invRadius;

        if(r < (T) 0.75) {
          return ((T) 0.0); 
        }
        else if(r > (T) 1.0) {
          return ((T) 0.0); 
        }
        
        return DerivHannWindow(r * (T) 2.0 - (T) 1.5);
    }

  protected:
    T DerivHannWindow(T t) const {
      if(t < -((T) 0.5) || t > ((T) 0.5)) {
        return 0;
      }
      else {
        const T Pi2 = 2.0 * M_PI;
        return - M_PI * sin(t * Pi2); 
      }
    }

  protected:
    T invRadius;    
};

#endif
