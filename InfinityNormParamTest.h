#ifndef InfinityNormParamTest
#define InfinityNormParamTest


template< typename T >
class InfinityNormParamTest {
  public:
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    InfinityNormParamTest(const T paramUpdateInfinityNormLimit) :
      paramUpdateInfinityNormLimit(paramUpdateInfinityNormLimit) {}

    bool operator()(ParamT *paramUpdate) const {
      for (int i = 0; i < 6; i++) {
        if(std::abs(negParamUpdate(i)) >  paramUpdateInfinityNormLimit) {
          return false; 
        }
      }
      
      return true; 
    }

  protected:
    T paramUpdateInfinityNormLimit;
};

#endif
