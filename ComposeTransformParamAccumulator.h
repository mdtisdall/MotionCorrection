#ifndef ComposeTransformParamAccumulator_h
#define ComposeTransformParamAccumulator_h


template< typename T >
class ComposeTransformParamAccumulator {
  public:
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    static ParamT accumulate(const ParamT *curParam, const ParamT *deltaParam) {
      typedef Eigen::AngleAxis<T> RotationT;

      ParamT retParam;

      const Eigen::Matrix<T, 3, 1> curRotVec = curParam->tail(3);

      const T curAngle = curRotVec.norm();

      Eigen::Matrix<T, 3, 1> curRotAxis;
      if(0 == curAngle) {
        curRotAxis << 1, 0, 0;
      }
      else {
        curRotAxis = curRotVec.normalized();
      }

      RotationT curRotation(curAngle, curRotAxis);
      
      const Eigen::Matrix<T, 3, 1> deltaRotVec = deltaParam->tail(3);

      const T deltaAngle = deltaRotVec.norm();

      Eigen::Matrix<T, 3, 1> deltaRotAxis;
      if(0 == deltaAngle) {
        deltaRotAxis << 1, 0, 0;
      }
      else {
        deltaRotAxis = deltaRotVec.normalized();
      }
      
      RotationT deltaRotation(deltaAngle, deltaRotAxis);
      
      // To update points, we want to apply cur, then delta. However,
      // we update points with the inverse of the transform we have
      // parameterized. Thus, updating points with a single new transform is
      // equivalent to applying:
      // new^-1 = delta^-1 . cur^-1
      // and so
      // new = cur . delta

      // combine the translations
      retParam.head(3) = curParam->head(3) + (curRotation * deltaParam->head(3));
      RotationT retRotation(curRotation * deltaRotation);
      retParam.tail(3) = retRotation.axis() * retRotation.angle();

      return retParam;
    }
};

#endif

