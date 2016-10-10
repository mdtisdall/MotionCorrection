#ifndef MovingWeightedResidualGradientAndHessian_h
#define MovingWeightedResidualGradientAndHessian_h

#include "ResidualGradientAndHessianBase.h"

template <
  typename _VolumeT,
  typename CoordT,
  typename _WeightFuncT,
  typename _WeightGradientFuncT,
  typename _ResidualOpT>
class MovingWeightedResidualGradientAndHessian :
  public ResidualGradientAndHessianBase <
    typename _VolumeT::value_type,
    CoordT > {
public:
  typedef _VolumeT VolumeT;
  typedef _WeightFuncT WeightFuncT;
  typedef _WeightGradientFuncT WeightGradientFuncT;
  typedef _ResidualOpT ResidualOpT;
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
 
  MovingWeightedResidualGradientAndHessian(
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
    weightGradientPointList(3, cubeSize * cubeSize * cubeSize),
    jointGradDz(cubeSize),
    jointGradDy(cubeSize),
    jointGradDx(cubeSize)
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
      PointT weightPoint;

      weightPoint(0) = refdz->wrapCoord((*pointList)(0,offset));
      weightPoint(1) = refdz->wrapCoord((*pointList)(1,offset));
      weightPoint(2) = refdz->wrapCoord((*pointList)(2,offset));
       
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
    
    typedef Eigen::Map< Eigen::Matrix<T, 1, Eigen::Dynamic > > RefDMatT;
    
    RefDMatT jointGradDzMat(jointGradDz.buffer, jointGradDz.totalPoints);      
    RefDMatT jointGradDyMat(jointGradDy.buffer, jointGradDy.totalPoints);      
    RefDMatT jointGradDxMat(jointGradDx.buffer, jointGradDx.totalPoints);      
    
    RefDMatT refdzMat(refdz->buffer, refdz->totalPoints);      
    RefDMatT refdyMat(refdy->buffer, refdy->totalPoints);      
    RefDMatT refdxMat(refdx->buffer, refdx->totalPoints);      

    jointGradDzMat.array() = weightFuncPoints.transpose().array() * refdzMat.array();
    jointGradDyMat.array() = weightFuncPoints.transpose().array() * refdyMat.array();
    jointGradDxMat.array() = weightFuncPoints.transpose().array() * refdxMat.array();

    jointGradDzMat.array() += residualOp->getUnweightedResidual()->transpose().array() * weightGradientPointList.row(0).array();
    jointGradDyMat.array() += residualOp->getUnweightedResidual()->transpose().array() * weightGradientPointList.row(1).array();
    jointGradDxMat.array() += residualOp->getUnweightedResidual()->transpose().array() * weightGradientPointList.row(2).array();

    // the first three colums of the M matrix are just negative copies
    // of the spatial gradients
    residualGradient->row(0).noalias() = -jointGradDzMat;
    residualGradient->row(1).noalias() = -jointGradDyMat;
    residualGradient->row(2).noalias() = -jointGradDxMat;

    // the last three colums of the M matrix are element-wise products 
    // of the point-lists and the spatial gradients  
    residualGradient->row(3).array() =
      pointList->row(2).array() * jointGradDyMat.array()
      - pointList->row(1).array() * jointGradDxMat.array();
   
    residualGradient->row(4).array() =
      pointList->row(0).array() * jointGradDxMat.array()
      - pointList->row(2).array() * jointGradDzMat.array();
    
    residualGradient->row(5).array() =
      pointList->row(1).array() * jointGradDzMat.array()
      - pointList->row(0).array() * jointGradDyMat.array();

    //std::cout << "residualGradient->col(0): " <<
    //  residualGradient->col(0).transpose() << std::endl;
    

//    std::cout << "residualGradient.col(16832) just first term: " <<
//      residualGradient->col(16832).transpose() << std::endl;
//
//    std::cout << "residualOp->getUnweightedResidual()->row(16832): " << 
//      residualOp->getUnweightedResidual()->row(16832) << std::endl;
//    
//    std::cout << "weightGradientPointList.col(16832): " << 
//      weightGradientPointList.col(16832).transpose() << std::endl;

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
  const WeightFuncT *weightFunc;
  const WeightGradientFuncT *weightGradientFunc;
  ResidualOpT *residualOp;
  PointListT *initialPointList;
  ResidualT weightFuncPoints;
  ResidualGradientT unweightedResidualGradient;
  PointListT weightGradientPointList;
  VolumeT jointGradDz, jointGradDy, jointGradDx;
};

#endif

