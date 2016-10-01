#ifndef Moving_Weighted_Gauss_Newton_Ref_Grad_h
#define Moving_Weighted_Gauss_Newton_Ref_Grad_h

#include "Gauss_Newton.h"

#include "StaticWeightedResidualOp.h"

#include <fcntl.h>
#include <unistd.h>
#ifdef LINUX
#include <stdlib.h>
#endif 

#include <iostream>

#include <cfloat>

template <
  typename _ResidualGradientAndHessianT,
  typename _InterpolatorT,
  typename _ParamAccumulatorT,
  typename _WeightFuncT,
  typename _WeightGradientFuncT,
  typename _ConvergenceTestT = void
  >
class Moving_Weighted_Gauss_Newton_Ref_Grad : 
  Gauss_Newton <
    StaticWeightedResidualOp<_InterpolatorT>,
    _ResidualGradientAndHessianT,
    _InterpolatorT, 
    _ParamAccumulatorT,
    _ConvergenceTestT >{
  public:
    typedef Gauss_Newton <
      StaticWeightedResidualOp<_InterpolatorT>,
      _ResidualGradientAndHessianT,
      _InterpolatorT,
      _ParamAccumulatorT,
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

    Moving_Weighted_Gauss_Newton_Ref_Grad(
      const InterpolatorT *interpRef,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      WeightFuncT *weightFunc,
      WeightGradientFuncT *weightGradientFunc,
      double *gradientAndHessianComputeTime = NULL
      ) :
      Parent(
        interpRef,
        new ResidualOpT(refdz->cubeSize, interpRef),
        new ResidualGradientAndHessianT(
          refdz->cubeSize, weightFunc, weightGradientFunc, NULL, NULL),
        refdz->cubeSize),
      refdz(refdz),
      refdy(refdy),
      refdx(refdx){
        this->residualGradientAndHessian->setResidualOp(this->residualOp);
        this->residualGradientAndHessian->setInitialPointList(
          &(this->pointList)); 
      }

  protected:
    typedef typename Parent::NewVolVecT NewVolVecT;
    typedef typename Parent::PointListT PointListT;
    
    const VolumeT *refdz;
    const VolumeT *refdy;
    const VolumeT *refdx;

  public:
    
    void minimize(
      const VolumeT *newVolume,
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

      this->computeResidual(newVolume, initialParam);

      this->residualGradientAndHessian->
        initializeResidualGradientAndApproxHessian(
          &(this->pointList),
          initialParam,
          refdz, refdy, refdx, 
          &(this->residualGradient), &(this->approxResidualHessian),
          &(this->residualHessianLDL),
          gradientAndHessianComputeTime);
      
      Parent::minimize(newVolume, this->refdz, this->refdy, this->refdx,
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
