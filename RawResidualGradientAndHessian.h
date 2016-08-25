#ifndef RawResidualGradientAndHessian_h
#define RawResidualGradientAndHessian_h

#include "ResidualGradientAndHessianBase.h"

template < typename _VolumeT, typename CoordT >
class RawResidualGradientAndHessian :
  public ResidualGradientAndHessianBase <
    typename _VolumeT::value_type,
    CoordT > {
public:
  typedef _VolumeT VolumeT;
  typedef typename VolumeT::value_type T;
  typedef ResidualGradientAndHessianBase <
    typename _VolumeT::value_type,
    CoordT > ParentT;
  typedef typename ParentT::PointListT PointListT;
  typedef Eigen::Matrix< T, 6, Eigen::Dynamic > ResidualGradientT;
  typedef Eigen::Matrix< T, 6, 6 > ResidualHessianT;
  typedef Eigen::LDLT< ResidualHessianT, Eigen::Upper > ResidualHessianLDLT;
 
  RawResidualGradientAndHessian(const size_t cubeSize) :
    ParentT(cubeSize)
    {}

  void generateResidualGradientAndApproxHessian(
    const PointListT *pointList,
    const VolumeT *refdz,
    const VolumeT *refdy,
    const VolumeT *refdx,
    ResidualGradientT *residualGradient,
    ResidualHessianT *approxResidualHessian,
    ResidualHessianLDLT *residualHessianLDL,
    double *elapsedTime = NULL 
    ) const {
 
    struct timeval timeBefore, timeAfter;

    if(NULL != elapsedTime) {
      gettimeofday(&timeBefore, NULL);
    }
  
    typedef Eigen::Map< Eigen::Matrix<T, 1, Eigen::Dynamic > > RefDMatT;

    RefDMatT refDzMat(refdz->buffer, refdz->totalPoints);      
    RefDMatT refDyMat(refdy->buffer, refdy->totalPoints);      
    RefDMatT refDxMat(refdx->buffer, refdx->totalPoints);      

    // the first three colums of the M matrix are just negative copies
    // of the spatial gradients
    residualGradient->row(0).noalias() = -refDzMat;
    residualGradient->row(1).noalias() = -refDyMat;
    residualGradient->row(2).noalias() = -refDxMat;

    // the last three colums of the M matrix are element-wise products 
    // of the point-lists and the spatial gradients  
    residualGradient->row(3).array() =
      pointList->row(2).array() * refDyMat.array();
    
    residualGradient->row(3).array() -=
      pointList->row(1).array() * refDxMat.array();
   
    residualGradient->row(4).array() =
      pointList->row(0).array() * refDxMat.array();
    
    residualGradient->row(4).array() -=
      pointList->row(2).array() * refDzMat.array();
    
    residualGradient->row(5).array() =
      pointList->row(1).array() * refDzMat.array();
   
    residualGradient->row(5).array() -=
      pointList->row(0).array() * refDyMat.array();

    //std::cout << "residualGradient->col(0): " <<
    //  residualGradient->col(0).transpose() << std::endl;

    // now we can compute the Hessian
    approxResidualHessian->setZero(6, 6);
   
    approxResidualHessian->template selfadjointView<Eigen::Upper>().rankUpdate(
      *(residualGradient));

    residualHessianLDL->compute(*approxResidualHessian);
    
    if(NULL != elapsedTime) { 
      gettimeofday(&timeAfter, NULL);
  
      *elapsedTime =
        ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
        ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
    }

//    std::cout << "approxResidualHessian: " << *approxResidualHessian << std::endl;
  }
};

#endif

