#ifndef Gauss_Newton_Ref_Grad_New_h
#define Gauss_Newton_Ref_Grad_New_h

#include "Gauss_Newton.h"

#include "RawResidualOp.h"
#include "RawResidualGradientAndHessian.h"

template <
  typename _InterpolatorT,
  typename _ParamAccumulatorT,
  typename _ConvergenceTestT = void
  >
class Gauss_Newton_Ref_Grad :
  Gauss_Newton <
    RawResidualOp<_InterpolatorT>,
    RawResidualGradientAndHessian<
      typename _InterpolatorT::VolumeT,
      typename _InterpolatorT::CoordT
      >,
    _InterpolatorT,
    _ParamAccumulatorT,
    _ConvergenceTestT >{
  public:
    typedef Gauss_Newton <
      RawResidualOp<_InterpolatorT>,
      RawResidualGradientAndHessian<
        typename _InterpolatorT::VolumeT,
        typename _InterpolatorT::CoordT
        >,
      _InterpolatorT,
      _ParamAccumulatorT,
      _ConvergenceTestT > Parent;
    typedef typename Parent::ResidualOpT ResidualOpT;
    typedef typename Parent::ResidualGradientAndHessianT
      ResidualGradientAndHessianT;
    typedef typename Parent::InterpolatorT InterpolatorT;
    typedef typename Parent::ConvergenceTestT ConvergenceTestT;
    typedef typename Parent::VolumeT VolumeT;
    typedef typename Parent::CoordT CoordT;
    typedef typename Parent::T T;
    typedef typename Parent::ParamT ParamT;

    Gauss_Newton_Ref_Grad(
      const InterpolatorT *interpRef, 
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      double *gradientAndHessianComputeTime
      ) :
      Parent(
        interpRef,
        new ResidualOpT(refdz->cubeSize, interpRef),
        new ResidualGradientAndHessianT(refdz->cubeSize),
        refdz->cubeSize),
      refdz(refdz),
      refdy(refdy),
      refdx(refdx) {

      this->residualGradientAndHessian->
        initializeResidualGradientAndApproxHessian(
          &(this->pointList),
          refdz, refdy, refdx,
          &(this->residualGradient), &(this->approxResidualHessian),
          &(this->residualHessianLDL),
          gradientAndHessianComputeTime);
     
      //std::cout << "approxResidualHessian:" << std::endl <<
      //  this->approxResidualHessian << std::endl;
    }
    
    ~Gauss_Newton_Ref_Grad() {
      if(NULL != this->residualOp) {
        delete this->residualOp;
        this->residualOp = NULL;
      }
    } 
    
    void minimize(
      const VolumeT *newVolume,
      const ParamT *initialParam,
      ParamT *finalParam,
      const size_t maxSteps = 20,
      const T stepSizeScale = 0.25,
      const T stepSizeLimit = 0,
      const ConvergenceTestT *convergenceTest = NULL, 
      size_t *elapsedSteps = NULL,
      double *elapsedTime = NULL 
      ) {
      std::cout << "In Gauss_Newton_Ref_Grad_New::minimize()" << std::endl;

      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }

      Parent::minimize(newVolume, refdz, refdy, refdx,
        initialParam, finalParam,
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
    const VolumeT *refdz;
    const VolumeT *refdy;
    const VolumeT *refdx;

};


#endif
