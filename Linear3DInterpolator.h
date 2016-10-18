#ifndef Linear3DInterpolator_h
#define Linear3DInterpolator_h

#include "Interpolator3D.h"

#include <Eigen/Dense>
#include <cmath>

#include <vector>

template <typename VolumeT, typename CoordT>
class Linear3DInterpolator : public Interpolator3D<VolumeT, CoordT> {
  public:
    typedef typename VolumeT::value_type T;
    typedef Eigen::Matrix< T, 8, 8 >  Matrix_8_8_T;
    typedef Eigen::Matrix< T, 8, Eigen::Dynamic >  Matrix_8_X_T;
    
    T interp(
      const CoordT z,
      const CoordT y,
      const CoordT x) const {

      CoordT xInt, yInt, zInt;
      CoordT xFrac, yFrac, zFrac;

      xFrac = std::modf(x, &xInt);
      yFrac = std::modf(y, &yInt);
      zFrac = std::modf(z, &zInt);
   
      // ideally we don't want to allocate these on every call to interp
      // but for right now this ensures the method is thread-safe
      CoordT target_YArr[8];
      const Eigen::Map< Eigen::Matrix<T, 8, 1> > target_YVec(target_YArr);

      fill_target_Y(target_YArr, zFrac, yFrac, xFrac);
      return target_YVec.dot(
        coefficients.col(
          ( zInt * cubeSizeCoordT + yInt) *
            cubeSizeCoordT + xInt
          )
//        coefficients.col(
//          (((size_t) zInt) * cubeSize + ((size_t) yInt)) *
//            cubeSize + ((size_t) xInt)
//          )
        );
    }

  protected:
    Linear3DInterpolator(const VolumeT *volume) :
      Interpolator3D<VolumeT, CoordT>(volume),
      coefficientsInner(volume->totalPoints * 8),
      coefficients(&(coefficientsInner[0]), 8, volume->totalPoints),
      cubeSize(volume->cubeSize),
      cubeSizeCoordT(volume->cubeSize) {}

    void fill_target_Y(T* target_YArr, const T z, const T y, const T x) const {
        target_YArr[0] = (T) 1.0;
        target_YArr[1] = x; //x
        target_YArr[2] = y; //y
        target_YArr[3] = z; //z
        target_YArr[4] = target_YArr[1] * y; // x*y
        target_YArr[5] = target_YArr[2] * z; // y*z
        target_YArr[6] = target_YArr[3] * x; // x*z
        target_YArr[7] = target_YArr[4] * z; // xyz

    }


  protected:
    std::vector<T> coefficientsInner;
    Eigen::Map< Matrix_8_X_T > coefficients;
    const size_t cubeSize;
    const CoordT cubeSizeCoordT;

};
#endif
