#ifndef ResidualOpBase_h
#define ResidualOpBase_h

template < typename T, typename CoordT >
class ResidualOpBase {
public:
  typedef Eigen::Matrix< T, 3, 1 > PointT;
  typedef Eigen::Matrix<T, Eigen::Dynamic, 1> ResidualT;
  typedef typename Eigen::Matrix<T, 6, 1> ParamT;
  typedef Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, 1 > > NewVolVecT;
  typedef Eigen::Matrix< T, 3, Eigen::Dynamic > PointListT;

protected:
  ResidualOpBase(const size_t cubeSize) :
    cubeSize(cubeSize),
    cubeSizeDiv2( ((CoordT) cubeSize) / (CoordT) 2.0),
    cubeCenter(cubeCenterFromCubeSize(cubeSize)),
    cubeCenterPoint(cubeCenter, cubeCenter, cubeCenter) {}

  static CoordT cubeCenterFromCubeSize(const size_t cubeSize) {
      return ((CoordT) cubeSize)/(CoordT)2.0 - (CoordT)0.5;
  }

  const size_t cubeSize;
  const CoordT cubeSizeDiv2;
  const CoordT cubeCenter;
  const PointT cubeCenterPoint;
};

#endif

