#ifndef Weighted_Gauss_Newton_New_Grad_h
#define Weighted_Gauss_Newton_New_Grad_h

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
  typename _CircularMaskOpT
  >
class Weighted_Gauss_Newton_New_Grad : Gauss_Newton_Base<_InterpolatorT>{
  public:
    typedef Gauss_Newton_Base<_InterpolatorT> Parent;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;
    typedef _CircularMaskOpT CircularMaskOpT;

    Weighted_Gauss_Newton_New_Grad(
      const InterpolatorT *interpRef,
      const size_t cubeSize, 
      const CircularMaskOpT *circularMaskOp
      ) :
      Parent(interpRef, cubeSize),
      circularMaskOp(circularMaskOp) {
    }

  protected:
    typedef typename Parent::NewVolVecT NewVolVecT;
    
    virtual void computeResidual(
      const VolumeT *newVol,
      const NewVolVecT *newVolVec,
      const ParamT *param) {
        Parent::computeResidual(newVol, newVolVec, param);

        circularMaskOp->applyMask(&this->residual);
    }

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
      const T paramUpdate2NormLimit = 0,
      const T paramUpdateInfinityNormLimit = 0,
      const T paramUpdateMMLimit = 0,
      const T paramUpdateTransScaleMM = 0,
      const T paramUpdateRotScaleMM = 0,
      size_t *elapsedSteps = NULL, 
      double *elapsedTime = NULL,
      double *gradientAndHessianComputeTime = NULL
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }
     
      VolumeT weightedNewdz(this->cubeSize);
      VolumeT weightedNewdy(this->cubeSize);
      VolumeT weightedNewdx(this->cubeSize);

      circularMaskOp->applyMask(newdz, &weightedNewdz); 
      circularMaskOp->applyMask(newdy, &weightedNewdy); 
      circularMaskOp->applyMask(newdx, &weightedNewdx);

      this->generateResidualGradientAndApproxHessian(
        &(this->residualGradient), &(this->approxResidualHessian),
        &(this->pointList),
        &weightedNewdz, &weightedNewdy, &weightedNewdx,
        this->cubeSize, this->cubeCenter,
        gradientAndHessianComputeTime);

      this->residualHessianLDL.compute(this->approxResidualHessian);

      Parent::minimize(newVolume, initialParam, finalParam,
        maxSteps, stepSizeScale, stepSizeLimit,
        paramUpdate2NormLimit, paramUpdateInfinityNormLimit,
        paramUpdateMMLimit, paramUpdateTransScaleMM, paramUpdateRotScaleMM,
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
