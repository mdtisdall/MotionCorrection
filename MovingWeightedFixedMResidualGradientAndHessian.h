#ifndef MovingWeightedFixedMResidualGradientAndHessian_h
#define MovingWeightedFixedMResidualGradientAndHessian_h

#include "MovingWeightedResidualGradientAndHessian.h"

template <
  typename _VolumeT,
  typename CoordT,
  typename _WeightFuncT,
  typename _WeightGradientFuncT,
  typename _ResidualOpT>
class MovingWeightedFixedMResidualGradientAndHessian :
  public MovingWeightedResidualGradientAndHessian <
    _VolumeT,
    CoordT,
    _WeightFuncT,
    _WeightGradientFuncT,
    _ResidualOpT> {
public:
  typedef _VolumeT VolumeT;
  typedef _WeightFuncT WeightFuncT;
  typedef _WeightGradientFuncT WeightGradientFuncT;
  typedef _ResidualOpT ResidualOpT;
  typedef typename VolumeT::value_type T;
  typedef MovingWeightedResidualGradientAndHessian <
    _VolumeT,
    CoordT,
    _WeightFuncT,
    _WeightGradientFuncT,
    _ResidualOpT> ParentT;
  typedef typename ParentT::PointListT PointListT;
  typedef typename ParentT::PointT PointT;
  typedef typename ParentT::ParamT ParamT; 
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> ResidualT;
  typedef Eigen::Matrix< T, 6, Eigen::Dynamic > ResidualGradientT;
  typedef Eigen::Matrix< T, 6, 6 > ResidualHessianT;
  typedef Eigen::LDLT< ResidualHessianT, Eigen::Upper > ResidualHessianLDLT;
 
  MovingWeightedFixedMResidualGradientAndHessian(
    const size_t cubeSize,
    const WeightFuncT *weightFunc,
    const WeightGradientFuncT *weightGradientFunc,
    ResidualOpT *residualOp,
    PointListT *initialPointList) :
    ParentT(cubeSize),
    weightFunc(weightFunc),
    weightGradientFunc(weightGradientFunc),
    residualOp(residualOp),
    initialPointList(initialPointList),
    weightFuncPoints(cubeSize * cubeSize * cubeSize, 1),
    unweightedResidualGradient(6, cubeSize * cubeSize * cubeSize),
    weightGradientPointList(3, cubeSize * cubeSize * cubeSize)
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

    this->refdz = refdz; 
    this->refdy = refdy; 
    this->refdx = refdx; 

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

    PointT transPart = curParam->block(0,0,3,1);

    for(size_t offset = 0; offset < pointListLength; offset++) {
      PointT curPoint =
        initialPointList->col(offset) - transPart;
     
      PointT weightPoint;

      weightPoint(0) = refdz->wrapCoord(curPoint(0));
      weightPoint(1) = refdz->wrapCoord(curPoint(1));
      weightPoint(2) = refdz->wrapCoord(curPoint(2));
       
      weightGradientPointList.col(offset) = weightPoint * ( 
        (*weightGradientFunc)(
          weightPoint(0),
          weightPoint(1),
          weightPoint(2)) / weightPoint.norm());

//      if(16832 == offset) {
//        std::cout << "curPoint(16832): " << curPoint << std::endl;
//        std::cout << "weightPoint(16832): " << weightPoint << std::endl;
//        std::cout << "weightPoint(16832).norm(): " << weightPoint.norm() << std::endl; 
//        std::cout << "weightGradientFunc(16832): " << (*weightGradientFunc)(
//          weightPoint(0),
//          weightPoint(1),
//          weightPoint(2)) << std::endl; 
//      }

      weightFuncPoints(offset, 0) =
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
        weightFuncPoints.transpose().array() * unweightedResidualGradient.row(i).array();
    }

//    std::cout << "residualGradient.col(16832) just first term: " <<
//      residualGradient->col(16832).transpose() << std::endl;
//
//    std::cout << "residualOp->getUnweightedResidual()->row(16832): " << 
//      residualOp->getUnweightedResidual()->row(16832) << std::endl;
//    
//    std::cout << "weightGradientPointList.col(16832): " << 
//      weightGradientPointList.col(16832).transpose() << std::endl;

    residualGradient->row(0).array() -= 
      weightGradientPointList.row(0).array() *
        residualOp->getUnweightedResidual()->transpose().array();
    
    residualGradient->row(1).array() -= 
      weightGradientPointList.row(1).array() *
        residualOp->getUnweightedResidual()->transpose().array();
    
    residualGradient->row(2).array() -= 
      weightGradientPointList.row(2).array() *
        residualOp->getUnweightedResidual()->transpose().array();
    
//    std::cout << "residualGradient.col(16832): " <<
//      residualGradient->col(16832).transpose() << std::endl;
//
//        {
//          std::string filePath("Moving_Weighted_Gauss_Newton_Fixed_M_New_Grad_tests/residualGradient");
//          int outputFile = open(filePath.c_str(), O_WRONLY | O_CREAT, 0600);
//    
//          if(-1 == outputFile) {
//            std::cerr << "Could not open file: " << filePath << std::endl;
//            exit(1);
//          }
//
//          int bytesWritten = 
//              ::write(outputFile, residualGradient->data(),
//                sizeof(T) * 6 * this->cubeSize * this->cubeSize * this->cubeSize);
//    
//          close(outputFile);
//        }
//
//    exit(0);

    // now we can compute the Hessian
    approxResidualHessian->setZero(6, 6);
   
    approxResidualHessian->template selfadjointView<Eigen::Upper>().rankUpdate(
      *(residualGradient));

    residualHessianLDL->compute(*approxResidualHessian);
  }

  void setResidualOp(ResidualOpT *newResidualOp) {
    residualOp = newResidualOp; 
  }
  
  void setInitialPointList(PointListT *newInitialPointList) {
    initialPointList = newInitialPointList; 
  }
  
  
protected:  
  const GradientT *refdz, *refdy, *refdx; 
  const WeightFuncT *weightFunc;
  const WeightGradientFuncT *weightGradientFunc;
  ResidualOpT *residualOp;
  PointListT *initialPointList;
  ResidualT weightFuncPoints;
  ResidualGradientT unweightedResidualGradient;
  PointListT weightGradientPointList;
};

#endif

