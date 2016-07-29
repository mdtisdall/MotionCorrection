#ifndef Algorithms_h
#define Algorithms_h

#include "VolumeAtAddressable.h"

#include "MMParamTest.h"
#include "TrueParamTest.h"

#include "SumParamAccumulator.h"
#include "ComposeTransformParamAccumulator.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTOp.h"
#include "FFTWBuffer.h"

#include "CircularMaskOp.h"

template<typename InterpolatorT, typename dataT>
class Algorithm {
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef std::complex<dataT> complexT;
  typedef FFTOp<dataT> DataFFTOpT;
  typedef DataFFTOpT::spatialVolumeT DataVolumeT; 
  typedef DataFFTOpT::fourierVolumeT ComplexVolumeT; 
  typedef CircularMaskOp< DataVolumeT, DataVolumeT> DataCircularMaskOpT;
  typedef CircularMaskOp< ComplexVolumeT, 
    SymmetricHalfVolumeAtAddressable< FFTWBuffer<dataT> > >
    ComplexDataCircularMaskOpT;
  typedef SumParamAccumulator<dataT> SumParamAccumulatorT;
  typedef ComposeTransformParamAccumulator<dataT> ComposeParamAccumulatorT;
  typedef TrueParamTest<dataT> GradientUpdateTestT;
  typedef MMParamTest<dataT> ConvergenceTestT;
  typedef SumParamAccumulatorT::ParamT ParamT;  

  Algorithm(
    InterpolatorT *interpolator

  class Algorithm1 : public Weighted_Gauss_Newton_Ref_Grad<
    InterpolatorT,
    SumParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT,
    GradientUpdateTestT > {
    
    protected:
    
    typedef Weighted_Gauss_Newton_Ref_Grad<
      InterpolatorT,
      SumParamAccumulatorT,
      DataCircularMaskOpT,
      ConvergenceTestT,
      GradientUpdateTestT > Parent;
    
    public:
    
    Algorithm1 (
      const InterpolatorT *interpRef,
      const DataCircularMaskOpT *circularMaskOp,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      double *gradientAndHessianComputeTime
    ) : Parent(interpRef, circularMaskOp,
      refdz, refdy, refdx,
      gradientAndHessianComputeTime) {}

    void register(
      VolumeT *newVol,
      ParamT *p,
      double *elapsedTime,
      size_t *elapsedSteps) {
      
    }
};
  
  template<typename InterpolatorT>
  class Algorithm2 : public Weighted_Gauss_Newton_Ref_Grad<
    InterpolatorT,
    SumParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT > {
    
    protected:
    
    typedef Weighted_Gauss_Newton_Ref_Grad<
      InterpolatorT,
      SumParamAccumulatorT,
      DataCircularMaskOpT,
      ConvergenceTestT > Parent;
    
    public:
    
    Algorithm2 (
      const InterpolatorT *interpRef,
      const DataCircularMaskOpT *circularMaskOp,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      double *gradientAndHessianComputeTime
    ) : Parent(interpRef, circularMaskOp,
      refdz, refdy, refdx,
      gradientAndHessianComputeTime) {}
  };

  template<typename InterpolatorT>
  class Algorithm3 : public Weighted_Gauss_Newton_New_Grad<
    InterpolatorT,
    SumParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT,
    GradientUpdateTestT > {
    
    protected:
    
    typedef Weighted_Gauss_Newton_New_Grad<
    InterpolatorT,
    SumParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT,
    GradientUpdateTestT > Parent;
    
    public:
    
    Algorithm3 (
      const InterpolatorT *interpRef,
      const size_t cubeSize,
      const DataCircularMaskOpT *circularMaskOp
    ) : Parent(interpRef, cubeSize, circularMaskOp) {}
  };

  template<typename InterpolatorT>
  class Algorithm4 : public Weighted_Gauss_Newton_New_Grad<
    InterpolatorT,
    SumParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT > {
    
    protected:
    
    typedef Weighted_Gauss_Newton_New_Grad<
    InterpolatorT,
    SumParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT > Parent;
    
    public:
    
    Algorithm4 (
      const InterpolatorT *interpRef,
      const size_t cubeSize,
      const DataCircularMaskOpT *circularMaskOp
    ) : Parent(interpRef, cubeSize, circularMaskOp) {}
  };

  template<typename InterpolatorT>
  class Algorithm5 : public Weighted_Gauss_Newton_Ref_Grad<
    InterpolatorT,
    ComposeParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT,
    GradientUpdateTestT > {
    
    protected:
    
    typedef Weighted_Gauss_Newton_Ref_Grad<
    InterpolatorT,
    ComposeParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT,
    GradientUpdateTestT > Parent;
    
    public:
    
    Algorithm5 (
      const InterpolatorT *interpRef,
      const DataCircularMaskOpT *circularMaskOp,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      double *gradientAndHessianComputeTime
    ) : Parent(interpRef, circularMaskOp,
      refdz, refdy, refdx,
      gradientAndHessianComputeTime) {}
  };

  template<typename InterpolatorT>
  class Algorithm6 : public Weighted_Gauss_Newton_Ref_Grad<
    InterpolatorT,
    ComposeParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT > {
    
    protected:
    
    typedef Weighted_Gauss_Newton_Ref_Grad<
    InterpolatorT,
    ComposeParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT > Parent;
    
    public:
    
    Algorithm6 (
      const InterpolatorT *interpRef,
      const DataCircularMaskOpT *circularMaskOp,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      double *gradientAndHessianComputeTime
    ) : Parent(interpRef, circularMaskOp,
      refdz, refdy, refdx,
      gradientAndHessianComputeTime) {}
  };

  template<typename InterpolatorT>
  class Algorithm7 : public Weighted_Gauss_Newton_New_Grad<
    InterpolatorT,
    ComposeParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT,
    GradientUpdateTestT > {
    
    protected:
    
    typedef Weighted_Gauss_Newton_New_Grad<
    InterpolatorT,
    ComposeParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT,
    GradientUpdateTestT > Parent;
    
    public:
    
    Algorithm7 (
      const InterpolatorT *interpRef,
      const size_t cubeSize,
      const DataCircularMaskOpT *circularMaskOp
    ) : Parent(interpRef, cubeSize, circularMaskOp) {}
  };

  template<typename InterpolatorT>
  class Algorithm8 : public Weighted_Gauss_Newton_New_Grad<
    InterpolatorT,
    ComposeParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT > {
    
    protected:
    
    typedef Weighted_Gauss_Newton_New_Grad<
    InterpolatorT,
    ComposeParamAccumulatorT,
    DataCircularMaskOpT,
    ConvergenceTestT  > Parent;
    
    public:
    
    Algorithm8 (
      const InterpolatorT *interpRef,
      const size_t cubeSize,
      const DataCircularMaskOpT *circularMaskOp
    ) : Parent(interpRef, cubeSize, circularMaskOp) {}
  };
}

#endif
