#ifndef StaticWeightedResidualOp_h
#define StaticWeightedResidualOp_h

#include "ResidualOpBase.h"

#include "WeightFunction.h"

template <
  typename _InterpolatorT
  >
class StaticWeightedResidualOp :
  public ResidualOpBase<
    typename _InterpolatorT::VolumeT::value_type,
    typename _InterpolatorT::CoordT> {

public:
  typedef _InterpolatorT InterpolatorT;
  typedef typename InterpolatorT::VolumeT VolumeT;
  typedef typename InterpolatorT::CoordT CoordT;
  typedef typename VolumeT::value_type T;

protected:
  typedef ResidualOpBase<
    typename _InterpolatorT::VolumeT::value_type,
    typename _InterpolatorT::CoordT> ParentT;

public:
  typedef typename ParentT::PointT PointT;
  typedef typename ParentT::ResidualT ResidualT;
  typedef typename ParentT::ParamT ParamT;
  typedef typename ParentT::NewVolVecT NewVolVecT;
  typedef typename ParentT::PointListT PointListT;

public:
  StaticWeightedResidualOp(
    const size_t cubeSize,
    const InterpolatorT *interpRef) :
    ParentT(cubeSize),
    interpRef(interpRef),
    weightFunc(cubeSize),
    interpPoints(cubeSize * cubeSize * cubeSize, 1),
    weightPoints(cubeSize * cubeSize * cubeSize, 1),
    unweightedResidual(cubeSize * cubeSize * cubeSize, 1)
    {}

  void operator() (
    const VolumeT *newVol,
    const NewVolVecT *newVolVec,
    const PointListT *points,
    const ParamT *param,
    ResidualT *residual) {

    const size_t pointListLength = newVol->totalPoints;

    for(size_t offset = 0; offset < pointListLength; offset++) {
      PointT curPoint =
        points->col(offset) ;

      PointT wrappedPoint = curPoint + this->cubeCenterPoint;

      for(size_t i = 0; i < 3; i++) {
        wrappedPoint(i) = newVol->wrapIndex(wrappedPoint(i));
      }

      interpPoints(offset, 0) =
        interpRef->interp(
          wrappedPoint(0),
          wrappedPoint(1),
          wrappedPoint(2)
        ); 
     
      PointT weightPoint;

      for(size_t i = 0; i < 3; i++) {
        weightPoint(i) = newVol->wrapIndex(curPoint(i) + this->cubeSizeDiv2)
          - this->cubeSizeDiv2;
      }

      weightPoints(offset, 0) =
        weightFunc(
          weightPoint(0),
          weightPoint(1),
          weightPoint(2)
        ); 
      
//      if(12305 == offset) {
//        std::cout << "points[12305]: " << points->col(offset).transpose() << std::endl;
//        std::cout << "wrappedPoint[12305]: " << wrappedPoint.transpose() << std::endl;
//        std::cout << "weightPoint[12305]: " << weightPoint.transpose() << std::endl;
//      }

    }
 
//    std::cout << "weightPoints[12305]: " << weightPoints(12305) << std::endl;
//    std::cout << "interpPoints[12305]: " << interpPoints(12305) << std::endl;
//    std::cout << "newVolVec[12305]: " << (*newVolVec)(12305) << std::endl;

    unweightedResidual.array() = interpPoints.array() - (*newVolVec).array();

    residual->array() = weightPoints.array() * unweightedResidual.array();
  }

  const ResidualT* getUnweightedResidual() {
    return &unweightedResidual; 
  }

protected:
  const InterpolatorT *interpRef;
  const WeightFunction<T> weightFunc;
  ResidualT interpPoints;
  ResidualT weightPoints;
  ResidualT unweightedResidual;

};

#endif

