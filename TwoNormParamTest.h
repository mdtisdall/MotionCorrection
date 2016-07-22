#ifndef TwoNormParamTest_h
#define TwoNormParamTest_h


template< typename T >
class TwoNormParamTest {
  public:
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    TwoNormParamTest(const T paramUpdate2NormLimit) :
      paramUpdate2NormLimit(paramUpdate2NormLimit) {}

    bool operator()(ParamT *paramUpdate) const {
      return (paramUpdate->norm() < paramUpdate2NormLimit);
    }

  protected:
    T paramUpdate2NormLimit;
};

#endif
