#ifndef CircularMaskOp_h
#define CircularMaskOp_h

#include "SymmetricHalfVolumeAtAddressable.h"

#include <Eigen/Dense>

#include <cmath>

template <typename T>
class CircularMaskOpFunction {
  public:

    CircularMaskOpFunction(T cubeSize) :
      invRadius(((T) 2.0) / ((T) cubeSize)),
      eightOverThree(((T) 8.0) / ((T) 3.0))
      {}

    T maskValue(T z, T y, T x) {
      return maskValue( std::sqrt(z*z + y*y + x*x) );   
    }
    
    T maskValue(T radius) {
        T r = radius * invRadius;

        if(r < 0.75) {
          return ((T) 1.0); 
        }
        else if(r > 1) {
          return ((T) 0.0); 
        }
        
        return wcos(r * eightOverThree - ((T) 2.0));
    }

  protected:
    T wcos(T t) {
      if(t < -((T) 0.5) || t > ((T) 0.5)) {
        return 0;
      }
      else {
        return cos(((T) M_PI) * t); 
      }
    }

  protected:
    T invRadius;    
    const T eightOverThree;    
};

template < typename DataVolumeT, typename MaskVolumeT >
class CircularMaskOp_Base {
  protected:
    typedef typename DataVolumeT::value_type DataT;
    typedef typename MaskVolumeT::value_type MaskDataT;
    typedef Eigen::Map< Eigen::Matrix<MaskDataT, Eigen::Dynamic, 1> >
      MaskVolumeMapT;
    typedef Eigen::Map< Eigen::Matrix<DataT, Eigen::Dynamic, 1> >
      DataVolumeMapT;

  public:
    CircularMaskOp_Base(const size_t cubeSize) :
      cubeSize(cubeSize),
      volMask(cubeSize),
      volMaskMap(volMask.buffer, volMask.totalPoints)
      {}

    void applyMask(DataVolumeT *vol) const {
      DataVolumeMapT volMap(vol->buffer, vol->totalPoints);

      volMap = volMap.array() * volMaskMap.array();
    }
    
    void applyMask(Eigen::Matrix<DataT, Eigen::Dynamic, 1> *vol) const {
      vol->array() = vol->array() * volMaskMap.array();
    }
    
    void applyMask(const DataVolumeT *inVol, DataVolumeT *outVol) const {
      DataVolumeMapT inVolMap(inVol->buffer, inVol->totalPoints);
      DataVolumeMapT outVolMap(outVol->buffer, outVol->totalPoints);

      outVolMap.array() = inVolMap.array() * volMaskMap.array();
    }

  protected:
    size_t cubeSize;
    MaskVolumeT volMask;
    MaskVolumeMapT volMaskMap;
};

template < typename AtAddressableT, typename T >
class CircularMaskOp {
};


template < typename VolAtAddressableT, typename MaskAtAddressableT>
class CircularMaskOp<
    SymmetricHalfVolumeAtAddressable<VolAtAddressableT>,
    SymmetricHalfVolumeAtAddressable<MaskAtAddressableT>
  > :
  public CircularMaskOp_Base<
    SymmetricHalfVolumeAtAddressable<VolAtAddressableT>,
    SymmetricHalfVolumeAtAddressable<MaskAtAddressableT>
    >
    {
  public:
    typedef typename VolAtAddressableT::value_type VolDataT;
    typedef typename MaskAtAddressableT::value_type MaskDataT;
    typedef SymmetricHalfVolumeAtAddressable<VolAtAddressableT>
      SymmetricHalfVolT;

    typedef CircularMaskOp_Base<
        SymmetricHalfVolumeAtAddressable<VolAtAddressableT>,
        SymmetricHalfVolumeAtAddressable<MaskAtAddressableT>
      > Parent;

    CircularMaskOp(
      const size_t cubeSize,
      const MaskDataT maskVal = ((MaskDataT) 1.0)) :
      Parent(cubeSize),
      maskFunction(cubeSize)
    {
      const size_t lastDim = SymmetricHalfVolT::getLastFourierDimension(cubeSize);
      //const int twoLastDim = 2 * lastDim; 
      size_t offset = 0;
      for(int z = 0; z < cubeSize; z++) {
        int zIndex = z;
        if(zIndex >= lastDim - 1) {
          zIndex -= cubeSize; 
        }

        for(int y = 0; y < cubeSize; y++) {
          int yIndex = y;
          if(yIndex >= lastDim - 1) {
            yIndex -= cubeSize; 
          }
          
          for(int x = 0; x < lastDim; x++, offset++) {
            int xIndex = x;
            if(xIndex >= lastDim - 1) {
              xIndex -= cubeSize; 
            }

            this->volMask.at(offset) = maskVal *
              maskFunction.maskValue(zIndex, yIndex, xIndex);
          }
        }
      }
    }

  protected:
    CircularMaskOpFunction<MaskDataT> maskFunction;
};

template < typename VolAtAddressableT, typename MaskAtAddressableT>
class CircularMaskOp<
    VolumeAtAddressable<VolAtAddressableT>,
    VolumeAtAddressable<MaskAtAddressableT>
  > :
  public CircularMaskOp_Base<
    VolumeAtAddressable<VolAtAddressableT>,
    VolumeAtAddressable<MaskAtAddressableT>
    >
    {
  public:
    typedef typename VolAtAddressableT::value_type VolDataT;
    typedef typename MaskAtAddressableT::value_type MaskDataT;
    typedef VolumeAtAddressable<VolAtAddressableT>
      VolumeT;

    typedef CircularMaskOp_Base<
        VolumeAtAddressable<VolAtAddressableT>,
        VolumeAtAddressable<MaskAtAddressableT>
      > Parent;

    CircularMaskOp(
      const size_t cubeSize,
      const MaskDataT maskVal = ((MaskDataT) 1.0)) :
      Parent(cubeSize),
      maskFunction(cubeSize)
    {
      size_t offset = 0;

      MaskDataT startIndex = - ((MaskDataT) (cubeSize - 1)) / ((MaskDataT) 2.0);
      MaskDataT endIndex = startIndex + cubeSize; 

      for(MaskDataT z = startIndex; z < endIndex; z += ((MaskDataT) 1.0)) {

        for(MaskDataT y = startIndex; y < endIndex; y += ((MaskDataT) 1.0)) {
          
          for(MaskDataT x = startIndex; x < endIndex; x += ((MaskDataT) 1.0), offset++) {
            this->volMask.at(offset) = maskVal * 
              maskFunction.maskValue(z, y, x);
          }
        }
      }
    }

  protected:
    CircularMaskOpFunction<MaskDataT> maskFunction;
};

#endif
