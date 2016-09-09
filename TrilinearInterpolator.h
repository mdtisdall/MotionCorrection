#ifndef TrilinearInterpolator_h
#define TrilinearInterpolator_h

#include "Linear3DInterpolator.h"

#include <Eigen/Dense>
#include <vector>
#include <cmath>

template <
    typename VolumeT,
    typename CoordT
  >
class TrilinearInterpolator :
  public Linear3DInterpolator<VolumeT, CoordT> {

  public:
    typedef typename VolumeT::value_type T;
    typedef Eigen::Matrix< T, 8, 8 >  Matrix_8_8_T;
    typedef Eigen::Matrix< T, 8, Eigen::Dynamic >  Matrix_8_X_T;

    Matrix_8_8_T generate_X_inv(){
        Matrix_8_8_T X_inv(8,8);
        X_inv << 1,0,0,0,0,0,0,0,
            		-1,1,0,0,0,0,0,0,
            		-1,0,1,0,0,0,0,0,
            		-1,0,0,1,0,0,0,0,
                 1,-1,-1,0,1,0,0,0,
                 1,0,-1,-1,0,1,0,0,
                 1,-1,0,-1,0,0,1,0,
                -1,1,1,1,-1,-1,-1,1;
        return X_inv;
    }

  protected:
    static void constructCoeffVectorSubpart(
      T *temp,
      const VolumeT *volume,
      const int z,
      const int y,
      const int x) {
        size_t xPlus1 = volume->wrapIndex(x + 1);
        size_t yPlus1 = volume->wrapIndex(y + 1);
        size_t zPlus1 = volume->wrapIndex(z + 1);

        temp[0]= volume->at(z, y, x);
        temp[1]= volume->at(z, y, xPlus1);
        temp[2]= volume->at(z, yPlus1, x);
        temp[3]= volume->at(zPlus1, y, x);
        temp[4]= volume->at(z, yPlus1, xPlus1);
        temp[5]= volume->at(zPlus1, yPlus1, x);
        temp[6]= volume->at(zPlus1, y, xPlus1);
        temp[7]= volume->at(zPlus1, yPlus1, xPlus1); 
    }

    void computeCoefficients(
      const VolumeT *volume) {
          T temp[8];

          Eigen::Map< Eigen::Matrix<T, 8, 1> > tempVector(temp);

          const size_t totalPoints =
            volume->cubeSize * volume->cubeSize * volume->cubeSize;

          Eigen::Matrix<T, 8, Eigen::Dynamic> tempMat(8, totalPoints);

          size_t tempMatrixOffset = 0;

          for(int z = 0; z < volume->cubeSize; ++z){
            for(int y = 0; y < volume->cubeSize; ++y){
              for(int x = 0; x < volume->cubeSize; ++x, tempMatrixOffset++){
                //values of f(x,y,z) at each corner.
                constructCoeffVectorSubpart(temp, volume, z, y, x);

                // store the new coefficients in the temp matrix
                tempMat.col(tempMatrixOffset) = tempVector;
              }
            }
          }

          // having made the argument vectors for all the points,
          // now apply the X_inv matrix to all of them to get the
          // coefficient vectors used for interpolation
          this->coefficients.noalias() = X_inv * tempMat;
    }
 

  public:
    TrilinearInterpolator(
      const VolumeT *volume) :
        Linear3DInterpolator<VolumeT, CoordT>(volume),
        X_inv(generate_X_inv())
        {
          computeCoefficients(volume);
        }


  protected:
    const Matrix_8_8_T X_inv;
};

#endif
