#ifndef UpsampledTrilinearInterpolator_h
#define UpsampledTrilinearInterpolator_h

#include "TrilinearInterpolator.h"

#include <Eigen/Dense>
#include <vector>
#include <cmath>

template <
  typename UpsamplingInterpolatorT,
  int upsamplingFactor,
  typename VolumeT,
  typename CoordT
  >
class UpsampledTrilinearInterpolator :
  public TrilinearInterpolator<VolumeT, CoordT> {

  public:
    typedef typename VolumeT::value_type T;
    typedef TrilinearInterpolator<VolumeT, CoordT> ParentT; 

  public:
    UpsampledTrilinearInterpolator(
      const VolumeT *volume,
      const UpsamplingInterpolatorT *upsamplingInterpolator) :
        ParentT(
          createUpsampledVolume(volume, upsamplingInterpolator)),
        upsamplingFactorCoordT(upsamplingFactor) {
        upsampledVolume = this->volume; 
      }
    
    T interp(
      const CoordT z,
      const CoordT y,
      const CoordT x) const {
      return ParentT::interp(
        z * upsamplingFactorCoordT,
        y * upsamplingFactorCoordT,
        x * upsamplingFactorCoordT);
    }

    ~UpsampledTrilinearInterpolator() {
      delete this->volume; 
    }


  protected:
    const VolumeT *upsampledVolume;
    CoordT upsamplingFactorCoordT;

    static VolumeT* createUpsampledVolume(
      const VolumeT *volume,
      const UpsamplingInterpolatorT *upsamplingInterpolator) {
      
      size_t upsampledCubeSize = volume->cubeSize * upsamplingFactor;
      VolumeT *upsampledVolume = new VolumeT(upsampledCubeSize);

//      CoordT cubeCenter = ((CoordT) volume->cubeSize)/(CoordT)2.0 - (CoordT)0.5;

      CoordT invUpsamplingFactorCoordT = ((CoordT) 1.0) / ((CoordT) upsamplingFactor);

      size_t offset = 0;

      for(size_t z = 0; z < upsampledCubeSize; z++) {
//        CoordT zCoord = ((CoordT) z) * invUpsamplingFactorCoordT - cubeCenter; 
        CoordT zCoord = ((CoordT) z) * invUpsamplingFactorCoordT; 
      
        for(size_t y = 0; y < upsampledCubeSize; y++) {
//          CoordT yCoord = ((CoordT) y) * invUpsamplingFactorCoordT- cubeCenter; 
          CoordT yCoord = ((CoordT) y) * invUpsamplingFactorCoordT; 
          
          for(size_t x = 0; x < upsampledCubeSize; x++, offset++) {
//            CoordT xCoord = ((CoordT) x) * invUpsamplingFactorCoordT- cubeCenter;
            CoordT xCoord = ((CoordT) x) * invUpsamplingFactorCoordT;

            upsampledVolume->at(offset) =
              upsamplingInterpolator->interp(zCoord, yCoord, xCoord);
          }
        }
      }

      return upsampledVolume;
    }
};

#endif
