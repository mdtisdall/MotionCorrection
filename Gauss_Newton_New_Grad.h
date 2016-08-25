#ifndef Gauss_Newton_New_Grad_h
#define Gauss_Newton_New_Grad_h

#include "Gauss_Newton.h"

#include "RawResidualOp.h"
#include "RawResidualGradientAndHessian.h"

template <
  typename _InterpolatorT,
  typename _ParamAccumulatorT,
  typename _ConvergenceTestT = void,
  typename _GradientUpdateTestT = void
  >
class Gauss_Newton_New_Grad : 
  Gauss_Newton <
    RawResidualOp<_InterpolatorT>,
    RawResidualGradientAndHessian<
      typename _InterpolatorT::VolumeT,
      typename _InterpolatorT::CoordT
      >,
    _InterpolatorT,
    _ParamAccumulatorT,
    _ConvergenceTestT,
    _GradientUpdateTestT>{
  public:
    typedef Gauss_Newton <
      RawResidualOp<_InterpolatorT>,
      RawResidualGradientAndHessian<
        typename _InterpolatorT::VolumeT,
        typename _InterpolatorT::CoordT
        >,
      _InterpolatorT,
      _ParamAccumulatorT,
      _ConvergenceTestT,
      _GradientUpdateTestT > Parent;
    typedef typename Parent::ResidualOpT ResidualOpT;
    typedef typename Parent::ResidualGradientAndHessianT
      ResidualGradientAndHessianT;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::ConvergenceTestT ConvergenceTestT;
    typedef typename Parent::GradientUpdateTestT GradientUpdateTestT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;

    Gauss_Newton_New_Grad(
      const InterpolatorT *interpRef, 
      const size_t cubeSize
      ) :
      Parent(
        interpRef,
        new ResidualOpT(cubeSize, interpRef),
        new ResidualGradientAndHessianT(cubeSize),
        cubeSize) {}
    
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
      const GradientUpdateTestT *gradientUpdateTest = NULL, 
      size_t *elapsedSteps = NULL, 
      double *elapsedTime = NULL,
      double *gradientAndHessianComputeTime = NULL
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }
      
      this->residualGradientAndHessian->
        generateResidualGradientAndApproxHessian(
          &this->pointList, 
          newdz, newdy, newdx,
          &(this->residualGradient), &(this->approxResidualHessian),
          &(this->residualHessianLDL),  
          gradientAndHessianComputeTime);
     
      Parent::minimize(newVolume, newdz, newdy, newdx,
        initialParam, finalParam,
        maxSteps, stepSizeScale, stepSizeLimit,
        convergenceTest, gradientUpdateTest,
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
