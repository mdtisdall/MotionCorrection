#ifndef SumParamAccumulator_h
#define SumParamAccumulator_h


template< typename T >
class SumParamAccumulator {
  public:
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    static ParamT accumulate(const ParamT *curParam, const ParamT *deltaParam) {
      return (*curParam + *deltaParam); 
    }
};

#endif

