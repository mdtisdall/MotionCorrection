#ifndef Algorithm_h
#define Algorithm_h

#include "FFTWBuffer.h"
#include "FFTOp.h"

#include "TricubicInterpolator.h"
#include "CubicBSplineInterpolator.h"

#include <sys/time.h>



template < typename AlgorithmT >
class AlgorithmInterpolatorFactory {
};


template <
  typename dataT,
  typename VolumeT,
  typename InterpolatorT,
  typename MinimizerFactoryT >
class Algorithm {
  public: 

  typedef Algorithm<
    dataT, VolumeT, InterpolatorT, MinimizerFactoryT>
    AlgorithmT;

  typedef std::complex<dataT> complexT;
  typedef FFTOp<dataT> DataFFTOpT;
  typedef typename DataFFTOpT::spatialVolumeT DataVolumeT; 
  typedef typename DataFFTOpT::fourierVolumeT ComplexVolumeT; 
  typedef CircularMaskOp< DataVolumeT, DataVolumeT> DataCircularMaskOpT;
  typedef CircularMaskOp< ComplexVolumeT, 
    SymmetricHalfVolumeAtAddressable< FFTWBuffer<dataT> > >
    ComplexDataCircularMaskOpT;
  typedef typename MinimizerFactoryT::MinimizerT MinimizerT;
  typedef typename MinimizerT::ParamT ParamT;

  
  Algorithm(VolumeT *refVol) :
    maxSteps(50),
    stepSizeScale(0.25),
    stepSizeLimit(1e-5),
    cubeSize(refVol->cubeSize),
    cubeVectorLength(refVol->totalPoints),
    interpolator(NULL),
    fourierMaskedRefVol(cubeSize),
    fourierMaskedNewVol(cubeSize),
    fourierMaskOp(cubeSize,
      1.0/(refVol->max() * ((dataT) cubeVectorLength))),
    fourierData(cubeSize),
    fftOp(cubeSize),
    imageMaskOp(cubeSize) {

    fourierMaskVolume(refVol, &fourierMaskedRefVol); 

    interpolator =
      AlgorithmInterpolatorFactory<AlgorithmT>::newInterpolator(
        this, &fourierMaskedRefVol); 
    
    minimizerFactory = new MinimizerFactoryT(refVol);
    minimizer = minimizerFactory->newMinimizer(
      interpolator, &imageMaskOp);
  }

  ~Algorithm() {
    if(NULL != interpolator) {
      delete interpolator; 
      interpolator = NULL;
    }
    
    if(NULL != minimizer) {
      delete minimizer;
      minimizer = NULL;
    }
    
    if(NULL != minimizerFactory) {
      delete minimizerFactory;
      minimizerFactory = NULL;
    }
  }

  void registerNewVolume(
    VolumeT *newVol,
    ParamT *p,
    double *elapsedTime,
    size_t *elapsedSteps) {
    struct timeval timeBefore, timeAfter;

    if(NULL != elapsedTime) {
      gettimeofday(&timeBefore, NULL);
    }
      
    registerNewVolumeInner(newVol, p, elapsedSteps);

    if(NULL != elapsedTime) { 
      gettimeofday(&timeAfter, NULL);
  
      *elapsedTime =
        ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
        ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
    } 
  }

  protected:

  virtual void registerNewVolumeInner( 
    VolumeT *newVol,
    ParamT *p, 
    size_t *elapsedSteps) = 0;

  void fourierMaskVolume(VolumeT *inVol, VolumeT *fourierMaskedVol) { 
    fftOp.forward(inVol, &fourierData); 
    
    fourierMaskOp.applyMask(&fourierData);  
        
    fftOp.backward(&fourierData, fourierMaskedVol);   
  }

  public:

  const size_t maxSteps;
  const dataT stepSizeScale;
  const dataT stepSizeLimit;
  const size_t cubeSize;
  const size_t cubeVectorLength; 
  InterpolatorT *interpolator;
  VolumeT fourierMaskedRefVol;
  VolumeT fourierMaskedNewVol;
  ComplexDataCircularMaskOpT fourierMaskOp;
  ComplexVolumeT fourierData;
  DataFFTOpT fftOp;
  DataCircularMaskOpT imageMaskOp;
  MinimizerFactoryT *minimizerFactory;
  MinimizerT *minimizer;
};



template <
  typename dataT,
  typename VolumeT,
  typename MinimizerFactoryT>
class AlgorithmInterpolatorFactory<
  Algorithm < 
    dataT,
    VolumeT,
    CubicBSplineInterpolator<VolumeT, dataT>,
    MinimizerFactoryT >
  > {
  protected: 
  typedef CubicBSplineInterpolator<VolumeT, dataT>
    InterpolatorT;
  typedef Algorithm < 
    dataT,
    VolumeT,
    InterpolatorT,
    MinimizerFactoryT > AlgorithmT;

  public:

  static InterpolatorT* newInterpolator(
    AlgorithmT *algorithm,
    VolumeT *fourierMaskedRefVol) {
    return new InterpolatorT(&(algorithm->fourierMaskedRefVol));
  }
};

#endif

