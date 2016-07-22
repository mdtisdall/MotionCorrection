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
  typename _ConvergenceTestT,
  typename _CircularMaskOpT
  >
class Weighted_Gauss_Newton_Ref_Grad :
  Gauss_Newton_Base<_InterpolatorT, _ConvergenceTestT> {
  public:
    typedef Gauss_Newton_Base<_InterpolatorT, _ConvergenceTestT> Parent;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::ConvergenceTestT ConvergenceTestT;
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
      circularMaskOp(circularMaskOp) {
     
      VolumeT weightedRefdz(this->cubeSize);
      VolumeT weightedRefdy(this->cubeSize);
      VolumeT weightedRefdx(this->cubeSize);

      circularMaskOp->applyMask(refdz, &weightedRefdz); 
      circularMaskOp->applyMask(refdy, &weightedRefdy); 
      circularMaskOp->applyMask(refdx, &weightedRefdx);

      this->generateResidualGradientAndApproxHessian(
        &(this->residualGradient), &(this->approxResidualHessian),
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
      ConvergenceTestT *convergenceTest = NULL, 
      size_t *elapsedSteps = NULL, 
      double *elapsedTime = NULL 
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }

      Parent::minimize(newVolume, initialParam, finalParam,
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
    
  protected: 
    const CircularMaskOpT *circularMaskOp;

};


#endif
