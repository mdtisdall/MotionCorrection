#ifndef ResidualGradientAndHessianBase_h
#define ResidualGradientAndHessianBase_h


template < typename T, typename CoordT >
class ResidualGradientAndHessianBase {
public:
  typedef Eigen::Matrix< T, 3, 1 > PointT;
  typedef Eigen::Matrix< T, 3, Eigen::Dynamic > PointListT;

protected:

  ResidualGradientAndHessianBase(const size_t cubeSize) :
  cubeSize(cubeSize),
  cubeCenter(cubeCenterFromCubeSize(cubeSize)),
  cubeCenterPoint(cubeCenter, cubeCenter, cubeCenter) {}
    
  static CoordT cubeCenterFromCubeSize(const size_t cubeSize) {
      return ((CoordT) cubeSize)/(CoordT)2.0 - (CoordT)0.5;
  }

  const size_t cubeSize;
  const CoordT cubeCenter;
  const PointT cubeCenterPoint;
};

#endif
