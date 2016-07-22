#ifndef MMParamTest
#define MMParamTest


template< typename T >
class MMParamTest {
  public:
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    MMParamTest(
      const T paramUpdateMMLimit,
      const T paramUpdateTransScaleMM,
      const T paramUpdateRotScaleMM
      ) :
      paramUpdateMMLimit(paramUpdateMMLimit),
      paramUpdateTransScaleMM(paramUpdateTransScaleMM),
      paramUpdateRotScaleMM(paramUpdateRotScaleMM) {}

    bool operator()(ParamT *paramUpdate) const {
      T paramUpdateMMScore = negParamUpdate.head(3).norm() * 
        paramUpdateTransScaleMM;

      // paramUpdateMMScore += sqrt(2.0 * (1.0 -
      //   cos(negParamUpdate.tail(3).norm()))) * 
      //   paramUpdateRotScaleMM;

      // While the above is correct, both the cos and sqrt are
      // expensive. We note that 
      // sqrt(2.0 * (1.0 - cos(x))) ~= x
      // which allows us to approximate the above via

      paramUpdateMMScore += negParamUpdate.tail(3).norm() * 
        paramUpdateRotScaleMM;

      return (paramUpdateMMScore < paramUpdateMMLimit);
    }

  protected:
    T paramUpdateMMLimit;
    T paramUpdateTransScaleMM;
    T paramUpdateRotScaleMM;
};

#endif
