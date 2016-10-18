#ifndef Cubic3DInterpolator_h
#define Cubic3DInterpolator_h

#include "Interpolator3D.h"

#include <Eigen/Dense>
#include <cmath>

#include <vector>

template <typename VolumeT, typename CoordT>
class Cubic3DInterpolator : public Interpolator3D<VolumeT, CoordT> {
  public:
    typedef typename VolumeT::value_type T;
    typedef Eigen::Matrix< T, 64, 64 >  Matrix_64_64_T;
    typedef Eigen::Matrix< T, 64, Eigen::Dynamic >  Matrix_64_X_T;
    
    T interp(
      const CoordT z,
      const CoordT y,
      const CoordT x) const {

      CoordT xInt, yInt, zInt;
      CoordT xFrac, yFrac, zFrac;

      xFrac = std::modf(x, &xInt);
      yFrac = std::modf(y, &yInt);
      zFrac = std::modf(z, &zInt);

      const size_t xIntSizeT = xInt;
      const size_t yIntSizeT = yInt;
      const size_t zIntSizeT = zInt;

      // ideally we don't want to allocate these on every call to interp
      // but for right now this ensures the method is thread-safe
      CoordT target_YArr[64];
      const Eigen::Map< Eigen::Matrix<T, 64, 1> > target_YVec(target_YArr);

      fill_target_Y(target_YArr, zFrac, yFrac, xFrac);
      return target_YVec.dot(
        coefficients.col(
          (zIntSizeT * cubeSize + yIntSizeT) * cubeSize + xIntSizeT
          )
        );
    }

  protected:
    Cubic3DInterpolator(const VolumeT *volume) :
      Interpolator3D<VolumeT, CoordT>(volume),
      coefficientsInner(volume->totalPoints * 64),
      coefficients(&(coefficientsInner[0]), 64, volume->totalPoints),
      cubeSize(volume->cubeSize) {}

    void fill_target_Y(T* target_YArr, const T z, const T y, const T x) const {
        T xSq = x*x;
        T xCubed = xSq * x;

        target_YArr[0] = (T) 1.0;
        target_YArr[1] = x;
        target_YArr[2] = xSq;
        target_YArr[3] = xCubed; 
        target_YArr[4] = target_YArr[0] * y; // y
        target_YArr[5] = target_YArr[1] * y; // x*y
        target_YArr[6] = target_YArr[2] * y; // x*x*y
        target_YArr[7] = target_YArr[3] * y; // x*x*x*y
        target_YArr[8] = target_YArr[4] * y; // y*y
        target_YArr[9] = target_YArr[5] * y; // x*y*y
        target_YArr[10] = target_YArr[6] * y; // x*x*y*y
        target_YArr[11] = target_YArr[7] * y; // x*x*x*y*y
        target_YArr[12] = target_YArr[8] * y; // y*y*y
        target_YArr[13] = target_YArr[9] * y; // x*y*y*y
        target_YArr[14] = target_YArr[10] * y; // x*x*y*y*y
        target_YArr[15] = target_YArr[11] * y; // x*x*x*y*y*y
        
        target_YArr[16] = target_YArr[0] * z; // z
        target_YArr[17] = target_YArr[1] * z; // x*z
        target_YArr[18] = target_YArr[2] * z; // x*x*z
        target_YArr[19] = target_YArr[3] * z; // x*x*x*z
        target_YArr[20] = target_YArr[4] * z; // y*z
        target_YArr[21] = target_YArr[5] * z; // x*y*z
        target_YArr[22] = target_YArr[6] * z; // x*x*y*z
        target_YArr[23] = target_YArr[7] * z; // x*x*x*y*z
        target_YArr[24] = target_YArr[8] * z; // y*y*z
        target_YArr[25] = target_YArr[9] * z; // x*y*y*z
        target_YArr[26] = target_YArr[10] * z; // x*x*y*y*z
        target_YArr[27] = target_YArr[11] * z; // x*x*x*y*y*z
        target_YArr[28] = target_YArr[12] * z; // y*y*y*z
        target_YArr[29] = target_YArr[13] * z; // x*y*y*y*z
        target_YArr[30] = target_YArr[14] * z; // x*x*y*y*y*z
        target_YArr[31] = target_YArr[15] * z; // x*x*x*y*y*y*z

        target_YArr[32] = target_YArr[16] * z; // z*z
        target_YArr[33] = target_YArr[17] * z; // x*z*z
        target_YArr[34] = target_YArr[18] * z; // x*x*z*z
        target_YArr[35] = target_YArr[19] * z; // x*x*x*z*z
        target_YArr[36] = target_YArr[20] * z; // y*z*z
        target_YArr[37] = target_YArr[21] * z; // x*y*z*z
        target_YArr[38] = target_YArr[22] * z; // x*x*y*z*z
        target_YArr[39] = target_YArr[23] * z; // x*x*x*y*z*z
        target_YArr[40] = target_YArr[24] * z; // y*y*z*z
        target_YArr[41] = target_YArr[25] * z; // x*y*y*z*z
        target_YArr[42] = target_YArr[26] * z; // x*x*y*y*z*z
        target_YArr[43] = target_YArr[27] * z; // x*x*x*y*y*z*z
        target_YArr[44] = target_YArr[28] * z; // y*y*y*z*z
        target_YArr[45] = target_YArr[29] * z; // x*y*y*y*z*z
        target_YArr[46] = target_YArr[30] * z; // x*x*y*y*y*z*z
        target_YArr[47] = target_YArr[31] * z; // x*x*x*y*y*y*z*z

        target_YArr[48] = target_YArr[32] * z; // z*z*z;
        target_YArr[49] = target_YArr[33] * z; // x*z*z*z;
        target_YArr[50] = target_YArr[34] * z; // x*x*z*z*z;
        target_YArr[51] = target_YArr[35] * z; // x*x*x*z*z*z;
        target_YArr[52] = target_YArr[36] * z; // y*z*z*z;
        target_YArr[53] = target_YArr[37] * z; // x*y*z*z*z;
        target_YArr[54] = target_YArr[38] * z; // x*x*y*z*z*z;
        target_YArr[55] = target_YArr[39] * z; // x*x*x*y*z*z*z;
        target_YArr[56] = target_YArr[40] * z; // y*y*z*z*z;
        target_YArr[57] = target_YArr[41] * z; // x*y*y*z*z*z;
        target_YArr[58] = target_YArr[42] * z; // x*x*y*y*z*z*z;
        target_YArr[59] = target_YArr[43] * z; // x*x*x*y*y*z*z*z;
        target_YArr[60] = target_YArr[44] * z; // y*y*y*z*z*z;
        target_YArr[61] = target_YArr[45] * z; // x*y*y*y*z*z*z;
        target_YArr[62] = target_YArr[46] * z; // x*x*y*y*y*z*z*z;
        target_YArr[63] = target_YArr[47] * z; // x*x*x*y*y*y*z*z*z;
    }

  protected:
    std::vector<T> coefficientsInner;
    Eigen::Map< Matrix_64_X_T > coefficients;
    const size_t cubeSize;

};
#endif
