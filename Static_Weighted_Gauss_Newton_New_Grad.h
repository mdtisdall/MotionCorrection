#ifndef Static_Weighted_Gauss_Newton_New_Grad_h
#define Static_Weighted_Gauss_Newton_New_Grad_h

#include "Gauss_Newton.h"

#include "StaticWeightedResidualOp.h"
#include "StaticWeightedResidualGradientAndHessian.h"

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
  typename _ConvergenceTestT = void
  >
class Static_Weighted_Gauss_Newton_New_Grad : 
  Gauss_Newton <
    StaticWeightedResidualOp<_InterpolatorT>,
    StaticWeightedResidualGradientAndHessian<
      typename _InterpolatorT::VolumeT,
      typename _InterpolatorT::CoordT,
      _WeightFuncT
      >,
    _InterpolatorT, 
    _ParamAccumulatorT,
    _ConvergenceTestT >{
  public:
    typedef Gauss_Newton <
      StaticWeightedResidualOp<_InterpolatorT>,
      StaticWeightedResidualGradientAndHessian<
        typename _InterpolatorT::VolumeT,
        typename _InterpolatorT::CoordT,
        _WeightFuncT
        >,
      _InterpolatorT,
      _ParamAccumulatorT,
      _ConvergenceTestT > Parent;
    typedef typename Parent::ResidualOpT ResidualOpT;
    typedef typename Parent::ResidualGradientAndHessianT
      ResidualGradientAndHessianT;
    typedef _WeightFuncT WeightFuncT;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::ConvergenceTestT ConvergenceTestT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;

    Static_Weighted_Gauss_Newton_New_Grad(
      const InterpolatorT *interpRef,
      const size_t cubeSize,
      WeightFuncT *weightFunc
      ) :
      Parent(
        interpRef,
        new ResidualOpT(cubeSize, interpRef),
        new ResidualGradientAndHessianT(cubeSize, weightFunc),
        cubeSize) {}

  protected:
    typedef typename Parent::NewVolVecT NewVolVecT;
    typedef typename Parent::PointListT PointListT;

  public:
    
    void minimize(
      const VolumeT *newVolume,
      const VolumeT *newdz,
      const VolumeT *newdy,
      const VolumeT *newdx,
      const ParamT *initialParam,
      ParamT *finalParam,
      const size_t maxSteps = 20,
      const T stepSizeScale = 0.25,
      const T stepSizeLimit = 0,
      ConvergenceTestT *convergenceTest = NULL, 
      size_t *elapsedSteps = NULL, 
      double *elapsedTime = NULL,
      double *gradientAndHessianComputeTime = NULL
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }
     
      this->residualGradientAndHessian->
        initializeResidualGradientAndApproxHessian(
          &(this->pointList),
          newdz, newdy, newdx,
          &(this->residualGradient), &(this->approxResidualHessian),
          &(this->residualHessianLDL),
          gradientAndHessianComputeTime);

      Parent::minimize(newVolume, newdz, newdy, newdx,
        initialParam, finalParam,
        maxSteps, stepSizeScale, stepSizeLimit,
        convergenceTest, 
        elapsedSteps);

      if(NULL != elapsedTime) { 
        gettimeofday(&timeAfter, NULL);
  
        *elapsedTime =
          ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
          ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
      }
    }

};


#endif
