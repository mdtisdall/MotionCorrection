#ifndef Weighted_Gauss_Newton_Ref_Grad_h
#define Weighted_Gauss_Newton_Ref_Grad_h

#include "Gauss_Newton_Base.h"

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
  typename _CircularMaskOpT,
  typename _ConvergenceTestT = void,
  typename _GradientUpdateTestT = void
  >
class Weighted_Gauss_Newton_Ref_Grad :
  Gauss_Newton_Base <
    _InterpolatorT, 
    _ParamAccumulatorT,
    _ConvergenceTestT, 
    _GradientUpdateTestT >{
  public:
    typedef Gauss_Newton_Base <
      _InterpolatorT,
      _ParamAccumulatorT,
      _ConvergenceTestT,
      _GradientUpdateTestT > Parent;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::ConvergenceTestT ConvergenceTestT;
    typedef typename Parent::GradientUpdateTestT GradientUpdateTestT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;
    typedef _CircularMaskOpT CircularMaskOpT;

    Weighted_Gauss_Newton_Ref_Grad(
      const InterpolatorT *interpRef,
      const CircularMaskOpT *circularMaskOp,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      double *gradientAndHessianComputeTime
      ) :
      Parent(interpRef, refdz->cubeSize),
      circularMaskOp(circularMaskOp),
      weightedRefdz(computeWeightedDerivative(refdz, circularMaskOp)),
      weightedRefdy(computeWeightedDerivative(refdy, circularMaskOp)),
      weightedRefdx(computeWeightedDerivative(refdx, circularMaskOp)) {
     
      this->generateResidualGradientAndApproxHessian(
        &(this->residualGradient), &(this->approxResidualHessian),
        &(this->residualHessianLDL),
        &(this->pointList),
        &weightedRefdz, &weightedRefdy, &weightedRefdx,
        this->cubeSize, this->cubeCenter,
        gradientAndHessianComputeTime);
     
//      std::cout << "approxResidualHessian:" << std::endl <<
//        approxResidualHessian << std::endl;
    }

  protected:
    typedef typename Parent::NewVolVecT NewVolVecT;
    typedef typename Parent::PointListT PointListT;

    static VolumeT computeWeightedDerivative(
      const VolumeT *refD,
      const CircularMaskOpT *circularMaskOp) {
      VolumeT weightedRefD(refD->cubeSize);
      circularMaskOp->applyMask(refD, &weightedRefD);
      return weightedRefD;
    }

    virtual void computeResidual(
      const VolumeT *newVol,
      const NewVolVecT *newVolVec,
      const PointListT *initialPoints,
      const ParamT *param) {
        Parent::computeResidual(newVol, newVolVec, initialPoints, param);
        
        circularMaskOp->applyMask(&this->residual);
    }

  public:
    
    void minimize(
      const VolumeT *newVolume,
      const ParamT *initialParam,
      ParamT *finalParam,
      const size_t maxSteps = 20,
      const T stepSizeScale = 0.25,
      const T stepSizeLimit = 0,
      const ConvergenceTestT *convergenceTest = NULL, 
      const GradientUpdateTestT *gradientUpdateTest = NULL, 
      size_t *elapsedSteps = NULL, 
      double *elapsedTime = NULL 
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }

      Parent::minimize(newVolume, &weightedRefdz, &weightedRefdy, &weightedRefdx,
        initialParam, finalParam,
        maxSteps, stepSizeScale, stepSizeLimit,
        convergenceTest, gradientUpdateTest, 
        elapsedSteps);

      if(NULL != elapsedTime) { 
        gettimeofday(&timeAfter, NULL);
  
        *elapsedTime =
          ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
          ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
      }
    }
    
  protected: 
    const CircularMaskOpT *circularMaskOp;
    const VolumeT weightedRefdz;
    const VolumeT weightedRefdy;
    const VolumeT weightedRefdx;

};


#endif
