#ifndef RawResidualOp_h
#define RawResidualOp_h

#include "ResidualOpBase.h"

template <
  typename _InterpolatorT
  >
class RawResidualOp :
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
  RawResidualOp(
    const size_t cubeSize,
    const InterpolatorT *interpRef) :
    ParentT(cubeSize),
    interpRef(interpRef),
    interpPoints(cubeSize * cubeSize * cubeSize, 1)
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
        points->col(offset) + this->cubeCenterPoint;

      interpPoints(offset, 0) =
        interpRef->interp(
          newVol->wrapIndex(curPoint(0)),
          newVol->wrapIndex(curPoint(1)),
          newVol->wrapIndex(curPoint(2))
        ); 
    }
 
//    std::cout << "interpPoints[0]: " << interpPoints(0) << std::endl;
//    std::cout << "newVolVec[0]: " << (*newVolVec)(0) << std::endl;

    residual->noalias() = interpPoints - (*newVolVec);
  }

protected:
  const InterpolatorT *interpRef;
  ResidualT interpPoints;

};

#endif

