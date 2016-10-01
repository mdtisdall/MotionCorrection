#ifndef Moving_Weighted_Gauss_Newton_Moving_M_New_Grad_h
#define Moving_Weighted_Gauss_Newton_Moving_M_New_Grad_h

#include "Moving_Weighted_Gauss_Newton_New_Grad.h"

#include "MovingWeightedMovingMResidualGradientAndHessian.h"


template <
  typename _InterpolatorT,
  typename _ParamAccumulatorT,
  typename _WeightFuncT,
  typename _WeightGradientFuncT,
  typename _ConvergenceTestT = void
  >
class Moving_Weighted_Gauss_Newton_Moving_M_New_Grad : 
  public Moving_Weighted_Gauss_Newton_New_Grad <
    MovingWeightedMovingMResidualGradientAndHessian<
      typename _InterpolatorT::VolumeT,
      typename _InterpolatorT::CoordT,
      _WeightFuncT,
      _WeightGradientFuncT,
      StaticWeightedResidualOp<_InterpolatorT>
      >,
    _InterpolatorT, 
    _ParamAccumulatorT,
    _WeightFuncT,
    _WeightGradientFuncT,
    _ConvergenceTestT >{
  public:
    typedef Moving_Weighted_Gauss_Newton_New_Grad <
      MovingWeightedMovingMResidualGradientAndHessian<
        typename _InterpolatorT::VolumeT,
        typename _InterpolatorT::CoordT,
        _WeightFuncT,
        _WeightGradientFuncT,
        StaticWeightedResidualOp<_InterpolatorT>
        >,
      _InterpolatorT,
      _ParamAccumulatorT,
      _WeightFuncT,
      _WeightGradientFuncT,
      _ConvergenceTestT > Parent;
    typedef typename Parent::ResidualOpT ResidualOpT;
    typedef typename Parent::ResidualGradientAndHessianT
      ResidualGradientAndHessianT;
    typedef _WeightFuncT WeightFuncT;
    typedef _WeightGradientFuncT WeightGradientFuncT;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::ConvergenceTestT ConvergenceTestT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;

    Moving_Weighted_Gauss_Newton_Moving_M_New_Grad(
      const InterpolatorT *interpRef,
      const size_t cubeSize,
      WeightFuncT *weightFunc,
      WeightGradientFuncT *weightGradientFunc
      ) :
      Parent(interpRef, cubeSize, weightFunc, weightGradientFunc){}

  protected:
    typedef typename Parent::NewVolVecT NewVolVecT;
    typedef typename Parent::PointListT PointListT;

};


#endif
