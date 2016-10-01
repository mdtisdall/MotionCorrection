#ifndef Moving_Weighted_Gauss_Newton_Fixed_M_Ref_Grad_h
#define Moving_Weighted_Gauss_Newton_Fixed_M_Ref_Grad_h

#include "Moving_Weighted_Gauss_Newton_Ref_Grad.h"
#include "MovingWeightedFixedMResidualGradientAndHessian.h"

#include <fcntl.h>
#include <unistd.h>
#ifdef LINUX
#include <stdlib.h>
#endif 

#include <iostream>

#include <cfloat>

template <
  typename _InterpolatorT,
  typename _ParamAccumulatorT,
  typename _WeightFuncT,
  typename _WeightGradientFuncT,
  typename _ConvergenceTestT = void
  >
class Moving_Weighted_Gauss_Newton_Fixed_M_Ref_Grad : 
  public Moving_Weighted_Gauss_Newton_Ref_Grad <
    MovingWeightedFixedMResidualGradientAndHessian<
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
    typedef Moving_Weighted_Gauss_Newton_Ref_Grad <
      MovingWeightedFixedMResidualGradientAndHessian<
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

    Moving_Weighted_Gauss_Newton_Fixed_M_Ref_Grad(
      const InterpolatorT *interpRef,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      WeightFuncT *weightFunc,
      WeightGradientFuncT *weightGradientFunc,
      double *gradientAndHessianComputeTime = NULL
      ) :
      Parent(
        interpRef, refdz, refdy, refdx,
        weightFunc, weightGradientFunc,
        gradientAndHessianComputeTime) {}

  protected:
    typedef typename Parent::NewVolVecT NewVolVecT;
    typedef typename Parent::PointListT PointListT;
    
};


#endif
