#ifndef AlgorithmBase_h
#define AlgorithmBase_h

#include "FFTWBuffer.h"
#include "FFTOp.h"

#include "TrilinearInterpolator.h"

#include "TricubicInterpolator.h"
#include "CentralDifferenceDifferentiator.h"

#include "CubicBSplineInterpolator.h"

#include "CircularMaskOp.h"

#include "WeightFunction.h"
#include "DerivWeightFunction.h"

#include <sys/time.h>

template < typename AlgorithmT >
class AlgorithmInterpolatorFactory {
};

template <
  typename InterpolatorT,
  typename ConvergenceTestT >
class AlgorithmBase {
  public: 
  typedef typename InterpolatorT::VolumeT DataVolumeT;
  typedef DataVolumeT VolumeT;
  typedef typename DataVolumeT::value_type dataT;
  typedef typename Eigen::Matrix<dataT, 6, 1> ParamT;
  typedef WeightFunction<dataT> WeightFunctionT;
  typedef DerivWeightFunction<dataT> WeightGradientFunctionT;
  typedef FFTOp<dataT> DataFFTOpT;
  typedef typename DataFFTOpT::fourierVolumeT ComplexVolumeT; 
  typedef CircularMaskOp< ComplexVolumeT, 
    SymmetricHalfVolumeAtAddressable< FFTWBuffer<dataT> >,
    WeightFunctionT >
    ComplexDataCircularMaskOpT;

  
  AlgorithmBase(
    DataVolumeT *refVol,
    ConvergenceTestT *convergenceTest) :
    maxSteps(50),
    stepSizeScale(0.25),
    stepSizeLimit(1e-5),
    cubeSize(refVol->cubeSize),
    cubeVectorLength(refVol->totalPoints),
    interpolator(NULL),
    weightFunction(cubeSize),
    weightGradientFunction(cubeSize),
    fourierMaskedRefVol(cubeSize),
    fourierMaskedNewVol(cubeSize),
    fourierMaskOp(cubeSize,
      &weightFunction,
      1.0/(refVol->max() * ((dataT) cubeVectorLength))),
    fourierData(cubeSize),
    fftOp(cubeSize),
    convergenceTest(convergenceTest) {

    fourierMaskVolume(refVol, &fourierMaskedRefVol); 

    interpolator =
      AlgorithmInterpolatorFactory<InterpolatorT>::newInterpolator(
        &fourierMaskedRefVol); 
  }

  ~AlgorithmBase() {
    if(NULL != interpolator) {
      delete interpolator; 
      interpolator = NULL;
    }
  }

  void registerNewVolume(
    DataVolumeT *newVol,
    ParamT *p,
    double *elapsedTime,
    size_t *elapsedSteps) {
    struct timeval timeBefore, timeAfter;

    if(NULL != elapsedTime) {
      gettimeofday(&timeBefore, NULL);
    }
   
    this->fourierMaskVolume(newVol, &fourierMaskedNewVol);

    registerNewVolumeInner(&fourierMaskedNewVol, p, elapsedSteps);

    if(NULL != elapsedTime) { 
      gettimeofday(&timeAfter, NULL);
  
      *elapsedTime =
        ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
        ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
    } 
  }

  protected:

  virtual void registerNewVolumeInner( 
    DataVolumeT *newVol,
    ParamT *p, 
    size_t *elapsedSteps) = 0;

  void fourierMaskVolume(DataVolumeT *inVol, DataVolumeT *fourierMaskedVol) { 
    fftOp.forward(inVol, &fourierData); 
    
    fourierMaskOp.applyMask(&fourierData);  
        
    fftOp.backward(&fourierData, fourierMaskedVol);   
  }

  protected:

  const size_t maxSteps;
  const dataT stepSizeScale;
  const dataT stepSizeLimit;
  const size_t cubeSize;
  const size_t cubeVectorLength; 
  InterpolatorT *interpolator;
  WeightFunctionT weightFunction;
  WeightGradientFunctionT weightGradientFunction;
  DataVolumeT fourierMaskedRefVol;
  DataVolumeT fourierMaskedNewVol;
  ComplexDataCircularMaskOpT fourierMaskOp;
  ComplexVolumeT fourierData;
  DataFFTOpT fftOp;
  ConvergenceTestT *convergenceTest;
};


template <typename DataVolumeT>
class AlgorithmInterpolatorFactory<
    TrilinearInterpolator<
      DataVolumeT, typename DataVolumeT::value_type
      >
  > {
  protected: 
  typedef TrilinearInterpolator<
    DataVolumeT,
    typename DataVolumeT::value_type
    > InterpolatorT;

  public:

  static InterpolatorT* newInterpolator(
    DataVolumeT *fourierMaskedRefVol) {
    return new InterpolatorT(fourierMaskedRefVol);
  }
};

template <typename DataVolumeT>
class AlgorithmInterpolatorFactory<
    TricubicInterpolator<
      DataVolumeT, typename DataVolumeT::value_type
      >
  > {
  protected: 
  typedef TricubicInterpolator<
    DataVolumeT,
    typename DataVolumeT::value_type
    > InterpolatorT;

  public:

  static InterpolatorT* newInterpolator(
    DataVolumeT *fourierMaskedRefVol) {
   
    const size_t cubeSize = fourierMaskedRefVol->cubeSize;

    CentralDifferencesDifferentiator<DataVolumeT>
      volDiffer(fourierMaskedRefVol);
    DataVolumeT dx(cubeSize);
    volDiffer.xDerivative(&dx);

    DataVolumeT dy(cubeSize);
    volDiffer.yDerivative(&dy);

    DataVolumeT dz(cubeSize);
    volDiffer.zDerivative(&dz);

    CentralDifferencesDifferentiator<DataVolumeT> dxDiffer(&dx);
   
    DataVolumeT dxy(cubeSize);
    dxDiffer.yDerivative(&dxy);
    
    DataVolumeT dxz(cubeSize);
    dxDiffer.zDerivative(&dxz);

    CentralDifferencesDifferentiator<DataVolumeT> dyDiffer(&dy);

    DataVolumeT dyz(cubeSize);
    dyDiffer.zDerivative(&dyz);

    CentralDifferencesDifferentiator<DataVolumeT> dxyDiffer(&dxy);
    
    DataVolumeT dxyz(cubeSize);
    dxyDiffer.zDerivative(&dxyz);
    return new InterpolatorT(
      fourierMaskedRefVol, &dx, &dy, &dz, &dxy, &dxz, &dyz, &dxyz);
  }
};

template <typename DataVolumeT>
class AlgorithmInterpolatorFactory<
    CubicBSplineInterpolator<
      DataVolumeT, typename DataVolumeT::value_type
      >
  > {
  protected: 
  typedef CubicBSplineInterpolator<
    DataVolumeT,
    typename DataVolumeT::value_type
    > InterpolatorT;

  public:

  static InterpolatorT* newInterpolator(
    DataVolumeT *fourierMaskedRefVol) {
    return new InterpolatorT(fourierMaskedRefVol);
  }
};

#endif
