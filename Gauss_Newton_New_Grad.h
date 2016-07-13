#ifndef Gauss_Newton_New_Grad_h
#define Gauss_Newton_New_Grad_h

#include "Gauss_Newton_Base.h"

template <
  typename _InterpolatorT,
  typename _ConvergenceTestT 
  >
class Gauss_Newton_New_Grad : 
  Gauss_Newton_Base<_InterpolatorT, _ConvergenceTestT>{
  public:
    typedef Gauss_Newton_Base<_InterpolatorT, _ConvergenceTestT> Parent;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::ConvergenceTestT ConvergenceTestT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;

    Gauss_Newton_New_Grad(
      const InterpolatorT *interpRef, 
      const size_t cubeSize
      ) :
      Parent(interpRef, cubeSize) {}
    
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
      const ConvergenceTestT *convergenceTest = NULL, 
      size_t *elapsedSteps = NULL, 
      double *elapsedTime = NULL,
      double *gradientAndHessianComputeTime = NULL
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }
      
      this->generateResidualGradientAndApproxHessian(
        &(this->residualGradient), &(this->approxResidualHessian),
        &this->pointList,  
        newdz, newdy, newdx, this->cubeSize, this->cubeCenter,
        gradientAndHessianComputeTime);
     
      this->residualHessianLDL.compute(this->approxResidualHessian);

      Parent::minimize(newVolume, initialParam, finalParam,
        maxSteps, stepSizeScale, stepSizeLimit,
        convergenceTest, 
        elapsedSteps, NULL);

      if(NULL != elapsedTime) { 
        gettimeofday(&timeAfter, NULL);
  
        *elapsedTime =
          ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
          ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
      }
    }
    
  protected: 
    typedef typename Parent::NewVolVecT NewVolVecT;

};


#endif
