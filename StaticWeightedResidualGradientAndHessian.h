#ifndef StaticWeightedResidualGradientAndHessian_h
#define StaticWeightedResidualGradientAndHessian_h

#include "ResidualGradientAndHessianBase.h"

template <
  typename _VolumeT,
  typename CoordT,
  typename _WeightFuncT>
class StaticWeightedResidualGradientAndHessian :
  public ResidualGradientAndHessianBase <
    typename _VolumeT::value_type,
    CoordT > {
public:
  typedef _VolumeT VolumeT;
  typedef _WeightFuncT WeightFuncT;
  typedef typename VolumeT::value_type T;
  typedef ResidualGradientAndHessianBase <
    typename _VolumeT::value_type,
    CoordT > ParentT;
  typedef typename ParentT::PointListT PointListT;
  typedef typename ParentT::PointT PointT;
  typedef typename ParentT::ParamT ParamT; 
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> ResidualT;
  typedef Eigen::Matrix< T, 6, Eigen::Dynamic > ResidualGradientT;
  typedef Eigen::Matrix< T, 6, 6 > ResidualHessianT;
  typedef Eigen::LDLT< ResidualHessianT, Eigen::Upper > ResidualHessianLDLT;
 
  StaticWeightedResidualGradientAndHessian(
    const size_t cubeSize,
    const WeightFuncT *weightFunc) :
    ParentT(cubeSize),
    weightFunc(weightFunc),
    weightPoints(cubeSize * cubeSize * cubeSize, 1),
    unweightedResidualGradient(6, cubeSize * cubeSize * cubeSize)
    {}


  
  
  void initializeResidualGradientAndApproxHessian(
    const PointListT *pointList,
    const ParamT *curParam,
    const VolumeT *refdz,
    const VolumeT *refdy,
    const VolumeT *refdx,
    ResidualGradientT *residualGradient,
    ResidualHessianT *approxResidualHessian,
    ResidualHessianLDLT *residualHessianLDL,
    double *elapsedTime = NULL 
    ) {
 
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
    unweightedResidualGradient.row(0).noalias() = -refDzMat;
    unweightedResidualGradient.row(1).noalias() = -refDyMat;
    unweightedResidualGradient.row(2).noalias() = -refDxMat;

    // the last three colums of the M matrix are element-wise products 
    // of the point-lists and the spatial gradients  
    unweightedResidualGradient.row(3).array() =
      pointList->row(2).array() * refDyMat.array();
    
    unweightedResidualGradient.row(3).array() -=
      pointList->row(1).array() * refDxMat.array();
   
    unweightedResidualGradient.row(4).array() =
      pointList->row(0).array() * refDxMat.array();
    
    unweightedResidualGradient.row(4).array() -=
      pointList->row(2).array() * refDzMat.array();
    
    unweightedResidualGradient.row(5).array() =
      pointList->row(1).array() * refDzMat.array();
   
    unweightedResidualGradient.row(5).array() -=
      pointList->row(0).array() * refDyMat.array();

    //std::cout << "residualGradient->col(0): " <<
    //  residualGradient->col(0).transpose() << std::endl;
   
    updateResidualGradientAndApproxHessian(
      pointList, curParam, refdz, refdy, refdx,
      residualGradient, approxResidualHessian, residualHessianLDL);

    if(NULL != elapsedTime) { 
      gettimeofday(&timeAfter, NULL);
  
      *elapsedTime =
        ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
        ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
    }

//    std::cout << "approxResidualHessian: " << *approxResidualHessian << std::endl;
  }
  
  void updateResidualGradientAndApproxHessian(
    const PointListT *pointList,
    const ParamT *curParam,
    const VolumeT *refdz,
    const VolumeT *refdy,
    const VolumeT *refdx,
    ResidualGradientT *residualGradient,
    ResidualHessianT *approxResidualHessian,
    ResidualHessianLDLT *residualHessianLDL,
    double *elapsedTime = NULL 
    ) {
    
    const size_t pointListLength = refdz->totalPoints;

    for(size_t offset = 0; offset < pointListLength; offset++) {
      PointT curPoint =
        pointList->col(offset) + this->cubeCenterPoint;
      
      PointT wrappedPoint;

      wrappedPoint(0) = refdz->wrapIndex(curPoint(0));
      wrappedPoint(1) = refdz->wrapIndex(curPoint(1));
      wrappedPoint(2) = refdz->wrapIndex(curPoint(2));
      
      PointT weightPoint = wrappedPoint - this->cubeCenterPoint;
      
      weightPoints(offset, 0) =
        (*weightFunc)(
          weightPoint(0),
          weightPoint(1),
          weightPoint(2)
        ); 
    }
    
    // weight the residual gradients based on where our sample
    // points are now located
    for(unsigned int i = 0; i < 6; i++) {
      residualGradient->row(i).array() = 
        weightPoints.transpose().array() * unweightedResidualGradient.row(i).array();
    }

    // now we can compute the Hessian
    approxResidualHessian->setZero(6, 6);
   
    approxResidualHessian->template selfadjointView<Eigen::Upper>().rankUpdate(
      *(residualGradient));

    residualHessianLDL->compute(*approxResidualHessian);
  }
  
protected:  
  const WeightFuncT *weightFunc;
  ResidualT weightPoints;
  ResidualGradientT unweightedResidualGradient;
};

#endif

