#ifndef WeightFunction_h
#define WeightFunction_h

#include <Eigen/Dense>

#include <cmath>

template <typename T>
class WeightFunction {
  public:
    typedef T value_type;

    WeightFunction(const T cubeSize) :
      invRadius(((T) 2.0) / ((T) cubeSize))
      {}

    T operator() (T z, T y, T x) const {
      return (*this)( std::sqrt(z*z + y*y + x*x) );   
    }
    
    T operator() (T radius) const {
        T r = radius * invRadius;

        if(r < 0.75) {
          return ((T) 1.0); 
        }
        else if(r > 1) {
          return ((T) 0.0); 
        }
        
        return HannWindow(r * (T) 2.0 - (T) 1.5);
    }

  protected:
    T HannWindow(T t) const {
      if(t < -((T) 0.5) || t > ((T) 0.5)) {
        return 0;
      }
      else {
        const T cos2PiT = cos(t * (T) (2.0 * M_PI));
        return 0.5 + 0.5 * cos2PiT; 
      }
    }

  protected:
    T invRadius;    
};

#endif
