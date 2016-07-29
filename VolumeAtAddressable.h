#ifndef VolumeAtAddressable_h
#define VolumeAtAddressable_h

#include <cmath>
#include <algorithm>

#ifdef LINUX
#include <cstddef>
#endif

template < typename AtAddressableT >
class VolumeAtAddressable : public AtAddressableT {
  public:
    typedef typename AtAddressableT::value_type value_type;

    VolumeAtAddressable(
      const size_t cubeSize
    ) :
    AtAddressableT(cubeSize * cubeSize * cubeSize), 
    cubeSize(cubeSize),
    cubeSizeInt(cubeSize),
    cubeSizeFloat(cubeSize),
    cubeSizeDouble(cubeSize),
    totalPoints(cubeSize * cubeSize * cubeSize),
    maxCubeIndex(cubeSize - 1),
    buffer(&(asAtAddressable().at(0)))
    {}
    
    VolumeAtAddressable(
      const size_t cubeSize,
      const size_t atAddressableSize
    ) :
    AtAddressableT(atAddressableSize),
    cubeSize(cubeSize),
    cubeSizeInt(cubeSize),
    cubeSizeFloat(cubeSize),
    cubeSizeDouble(cubeSize),
    totalPoints(cubeSize * cubeSize * cubeSize),
    maxCubeIndex(cubeSize - 1),
    buffer(&(asAtAddressable().at(0)))
    {}
    
    VolumeAtAddressable(
      const size_t cubeSize,
      const AtAddressableT atAddressable
    ) :
    AtAddressableT(atAddressable),
    cubeSize(cubeSize),
    cubeSizeInt(cubeSize),
    cubeSizeFloat(cubeSize),
    cubeSizeDouble(cubeSize),
    totalPoints(cubeSize * cubeSize * cubeSize),
    maxCubeIndex(cubeSize - 1),
    buffer(&(asAtAddressable().at(0)))
    {}

    const value_type& at(const size_t z, const size_t y, const size_t x) const {
        return asAtAddressable().at((z * cubeSize + y) * cubeSize + x);
    }
    
    value_type& at(const size_t z, const size_t y, const size_t x) {
        return asAtAddressable().at((z * cubeSize + y) * cubeSize + x);
    } 
    
    const value_type& at(const size_t offset) const {
        return asAtAddressable().at(offset);
    }
    
    value_type& at(const size_t offset) {
        return asAtAddressable().at(offset);
    } 
   
    value_type max() {
      value_type max_ret = asAtAddressable().at(0);

      for(size_t offset = 1; offset < totalPoints; offset++) {
        max_ret = std::max(max_ret, asAtAddressable().at(offset));
      }

      return max_ret;
    }

    int clipIndex(const int index) const {
        if(index > maxCubeIndex) {
            return maxCubeIndex;
        }
        
        if (index < 0) {
            return 0;
        }
        
        return index;
    }
   
    int wrapIndex(int index) const {
        while(index < 0) {
          index += cubeSizeInt; 
        }

        while(index >= cubeSizeInt) {
          index -= cubeSizeInt;  
        }

        return index;

        //return (index + cubeSizeInt) % cubeSizeInt;
    }
    
    size_t wrapIndex(size_t index) const {
        // index can't be negative

        while(index >= cubeSize) {
          index -= cubeSize;  
        }

        return index;

        //return (index + cubeSize) % cubeSize;
    }
    
    float wrapIndex(float index) const {
        while(index < 0) {
          index += cubeSizeFloat; 
        }

        while(index >= cubeSizeFloat) {
          index -= cubeSizeFloat;  
        }

        return index;

        //return std::fmod(index + cubeSizeFloat, cubeSizeFloat);
    }
    
    double wrapIndex(double index) const {
        while(index < 0) {
          index += cubeSizeDouble; 
        }

        while(index >= cubeSizeDouble) {
          index -= cubeSizeDouble;  
        }

        return index;

        //return std::fmod(index + cubeSizeDouble, cubeSizeDouble);
    }
    

  protected:
    AtAddressableT& asAtAddressable() {  
      return static_cast<AtAddressableT&>(*this);
    }
    
    const AtAddressableT& asAtAddressable() const {  
      return static_cast<const AtAddressableT&>(*this);
    }

  public:
    const size_t cubeSize; 
    const int cubeSizeInt; 
    const float cubeSizeFloat; 
    const double cubeSizeDouble; 
    const size_t totalPoints;
    const size_t maxCubeIndex; 
    value_type* buffer;
};
    

#endif
