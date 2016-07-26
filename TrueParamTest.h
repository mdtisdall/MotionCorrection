#ifndef TrueParamTest_h
#define TrueParamTest_h


template< typename T >
class TrueParamTest {
  public:
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    TrueParamTest() {}

    bool operator()(const ParamT *paramUpdate) const {
      return true;
    }
};

#endif
