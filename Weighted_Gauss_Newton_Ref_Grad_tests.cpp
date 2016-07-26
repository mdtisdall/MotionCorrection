#include "catch.hpp"

#include "Weighted_Gauss_Newton_Ref_Grad.h"

#include "CubicBSplineInterpolator.h"
#include "VolumeAtAddressable.h"

#include "TwoNormParamTest.h"
#include "TrueParamTest.h"

#include "SumParamAccumulator.h"
#include "ComposeTransformParamAccumulator.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTOp.h"
#include "FFTWBuffer.h"

#include "CircularMaskOp.h"

#include "BinaryFile.h"

TEST_CASE("a weighted Gauss-Newton minimizer using reference-image gradients can be instantiated") {
  typedef float dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT; 

  const size_t cubeSize = 32;
  const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;
 
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > MaskVolumeT; 
  typedef CircularMaskOp< VolumeT, MaskVolumeT> CircularMaskOpT;
  typedef SumParamAccumulator<dataT> ParamAccumulatorT;
  typedef TwoNormParamTest<dataT> ConvergenceTestT;
  typedef Weighted_Gauss_Newton_Ref_Grad<
    InterpolatorT,
    ParamAccumulatorT,
    CircularMaskOpT,
    ConvergenceTestT> MinimizerT; 
  typedef MinimizerT::ParamT ParamT;
 
  VolumeT volume(cubeSize);
  
  for(size_t i = 0; i < cubeVectorLength; i++) {
      volume.at(i) = ((dataT) i) / ((dataT) cubeVectorLength); 
  }

  InterpolatorT interpolator(&volume);
  
  double gradientAndHessianComputeTime;

  CentralDifferencesDifferentiator<VolumeT> volDiffer(&volume);
  VolumeT dx(cubeSize, cubeVectorLength);
  volDiffer.xDerivative(&dx);

  VolumeT dy(cubeSize, cubeVectorLength);
  volDiffer.yDerivative(&dy);

  VolumeT dz(cubeSize, cubeVectorLength);
  volDiffer.zDerivative(&dz);

  CircularMaskOpT imageMaskOp(cubeSize);
  
  MinimizerT minimizer(&interpolator, &imageMaskOp,
    &dz, &dy, &dx,
    &gradientAndHessianComputeTime);
    
  WARN("elapsed time computing gradient and Hessian: "
    << gradientAndHessianComputeTime << " ms");

  SECTION("and registering an image with itself produces 0 transformation") {
    ParamT initialParam;
    initialParam << 0, 0, 0, 0, 0, 0;

    ParamT finalParam;

    size_t maxSteps = 20;
    const dataT stepSizeScale = 0.25;
    const dataT stepSizeLimit = 1e-5;

    const dataT paramUpdate2NormLimit = 1e-6;

    ConvergenceTestT convergenceTest(paramUpdate2NormLimit);

    double elapsedTime;
    size_t elapsedSteps;

    minimizer.minimize(&volume, &initialParam, &finalParam,
      maxSteps, stepSizeScale, stepSizeLimit,
      &convergenceTest,  NULL,
      &elapsedSteps, &elapsedTime);

    WARN("elapsed time: " << elapsedTime << " ms");
    WARN("elapsed steps: " << elapsedSteps);
    WARN("finalParam: " << finalParam.transpose());

    for(int i = 0; i < 6; i++) {
      REQUIRE(0 == Approx(finalParam(i)));
    }
  } 
}


TEST_CASE("a weighted Gauss-Newton minimizer using reference-image gradients can be instantiated from image data") {
  typedef double dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef CubicBSplineInterpolator<VolumeT, dataT> InterpolatorT; 

  const size_t cubeSize = 32;
  const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

  typedef std::complex<dataT> complexT;
  typedef FFTOp<dataT> DataFFTOpT;
  typedef DataFFTOpT::spatialVolumeT DataVolumeT; 
  typedef DataFFTOpT::fourierVolumeT ComplexVolumeT; 
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > ImageMaskVolumeT; 
  typedef CircularMaskOp< DataVolumeT, ImageMaskVolumeT> DataCircularMaskOpT;
  typedef SymmetricHalfVolumeAtAddressable< FFTWBuffer<dataT> >
    FourierMaskVolumeT; 
  typedef CircularMaskOp< ComplexVolumeT, FourierMaskVolumeT>
    ComplexDataCircularMaskOpT;
  
  const dataT maskScale =
    1.0/((dataT) cubeVectorLength);
  
  ComplexDataCircularMaskOpT fourierMaskOp(cubeSize, maskScale);
  DataVolumeT maskedRefVolume(cubeSize);
  DataVolumeT maskedNewVolume(cubeSize);

  {
    VolumeT refVolume(cubeSize);
    VolumeT newVolume(cubeSize);
  
    ComplexVolumeT fourierData(cubeSize);
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&refVolume,
                "Weighted_Gauss_Newton_Ref_Grad_tests/refVolInput.dat"));
    
    DataFFTOpT fftOp(cubeSize);
  
    fftOp.forward(&refVolume, &fourierData); 
    
    fourierMaskOp.applyMask(&fourierData);  
        
    fftOp.backward(&fourierData, &maskedRefVolume); 
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&newVolume,
                "Weighted_Gauss_Newton_Ref_Grad_tests/newVolInput.dat"));
    
    fftOp.forward(&newVolume, &fourierData);
    
    fourierMaskOp.applyMask(&fourierData); 
        
    fftOp.backward(&fourierData, &maskedNewVolume);
  } 

  InterpolatorT interpolator(&maskedRefVolume);

  CentralDifferencesDifferentiator<VolumeT> volDiffer(&maskedRefVolume);
  VolumeT dx(cubeSize, cubeVectorLength);
  volDiffer.xDerivative(&dx);

  VolumeT dy(cubeSize, cubeVectorLength);
  volDiffer.yDerivative(&dy);

  VolumeT dz(cubeSize, cubeVectorLength);
  volDiffer.zDerivative(&dz); 

  double gradientAndHessianComputeTime;
  
  DataCircularMaskOpT imageMaskOp(cubeSize);


  SECTION(std::string("and registering two images ") +
    std::string("without updating gradients ") + 
    std::string("returns identical result to Mathematica")) {
  
    typedef SumParamAccumulator<dataT> ParamAccumulatorT;
    
    typedef Weighted_Gauss_Newton_Ref_Grad<
      InterpolatorT, ParamAccumulatorT, DataCircularMaskOpT > MinimizerT; 
    typedef MinimizerT::ParamT ParamT;

    MinimizerT minimizer(&interpolator, &imageMaskOp,
      &dz, &dy, &dx,
      &gradientAndHessianComputeTime);
          
    WARN("elapsed time computing gradient and Hessian: "
      << gradientAndHessianComputeTime << " ms");
    
    ParamT initialParam;
    initialParam << 0, 0, 0, 0, 0, 0;
  
    ParamT finalParam;

    // We force this to stop as soon as we take a bad step, so that we get to
    // the same point as the Mathematica code
    size_t maxSteps = 20;
    const dataT stepSizeScale = 0.25;
    const dataT stepSizeLimit = 1.0;
 
    double elapsedTime;
    size_t elapsedSteps;
  
    minimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
      maxSteps, stepSizeScale, stepSizeLimit,
      NULL,  NULL,
      &elapsedSteps, &elapsedTime);

    std::vector<dataT> paramSolution(6);

    WARN("elapsed time: " << elapsedTime << " ms");
    WARN("elapsed steps: " << elapsedSteps);
    WARN("finalParam: " << finalParam.transpose());

    REQUIRE(6 * sizeof(dataT)
          == BinaryFile< std::vector<dataT> >::read(&paramSolution,
              "Weighted_Gauss_Newton_Ref_Grad_tests/parameterOutput.dat"));

    for(int i = 0; i < 6; i++) {
      REQUIRE(paramSolution[i] == Approx(finalParam(i)));
    }
  }


  SECTION(std::string("and registering two images ") +
    std::string("updating gradients every step ") + 
    std::string("returns identical result to Mathematica")) {
  
    typedef void ConvergenceTestT;
    typedef TrueParamTest<dataT> GradientUpdateTestT;
    typedef SumParamAccumulator<dataT> ParamAccumulatorT;
    
    typedef Weighted_Gauss_Newton_Ref_Grad<
      InterpolatorT,
      ParamAccumulatorT,
      DataCircularMaskOpT,
      ConvergenceTestT,
      GradientUpdateTestT> MinimizerT; 
    typedef MinimizerT::ParamT ParamT;

    MinimizerT minimizer(&interpolator, &imageMaskOp,
      &dz, &dy, &dx,
      &gradientAndHessianComputeTime);
          
    WARN("elapsed time computing gradient and Hessian: "
      << gradientAndHessianComputeTime << " ms");
    
    ParamT initialParam;
    initialParam << 0, 0, 0, 0, 0, 0;
  
    ParamT finalParam;

    // We force this to stop as soon as we take a bad step, so that we get to
    // the same point as the Mathematica code
    size_t maxSteps = 20;
    const dataT stepSizeScale = 0.25;
    const dataT stepSizeLimit = 1.0;
 
    double elapsedTime;
    size_t elapsedSteps;
  
    GradientUpdateTestT gradientUpdateTest;

    minimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
      maxSteps, stepSizeScale, stepSizeLimit,
      NULL,  &gradientUpdateTest,
      &elapsedSteps, &elapsedTime);

    std::vector<dataT> paramSolution(6);

    WARN("elapsed time: " << elapsedTime << " ms");
    WARN("elapsed steps: " << elapsedSteps);
    WARN("finalParam: " << finalParam.transpose());

    REQUIRE(6 * sizeof(dataT)
          == BinaryFile< std::vector<dataT> >::read(&paramSolution,
              "Weighted_Gauss_Newton_Ref_Grad_tests/gradientUpdateParameterOutput.dat"));

    for(int i = 0; i < 6; i++) {
      REQUIRE(paramSolution[i] == Approx(finalParam(i)));
    }
  }

  SECTION(std::string("and registering two images ") +
    std::string("updating gradients every step ") + 
    std::string("and using compose accumulator ") + 
    std::string("returns identical result to Mathematica")) {
  
    typedef void ConvergenceTestT;
    typedef TrueParamTest<dataT> GradientUpdateTestT;
    typedef ComposeTransformParamAccumulator<dataT> ParamAccumulatorT;
    
    typedef Weighted_Gauss_Newton_Ref_Grad<
      InterpolatorT,
      ParamAccumulatorT,
      DataCircularMaskOpT,
      ConvergenceTestT,
      GradientUpdateTestT> MinimizerT; 
    typedef MinimizerT::ParamT ParamT;

    MinimizerT minimizer(&interpolator, &imageMaskOp,
      &dz, &dy, &dx,
      &gradientAndHessianComputeTime);
          
    WARN("elapsed time computing gradient and Hessian: "
      << gradientAndHessianComputeTime << " ms");
    
    ParamT initialParam;
    initialParam << 0, 0, 0, 0, 0, 0;
  
    ParamT finalParam;

    // We force this to stop as soon as we take a bad step, so that we get to
    // the same point as the Mathematica code
    size_t maxSteps = 20;
    const dataT stepSizeScale = 0.25;
    const dataT stepSizeLimit = 1.0;
 
    double elapsedTime;
    size_t elapsedSteps;
  
    GradientUpdateTestT gradientUpdateTest;

    minimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
      maxSteps, stepSizeScale, stepSizeLimit,
      NULL,  &gradientUpdateTest,
      &elapsedSteps, &elapsedTime);

    std::vector<dataT> paramSolution(6);

    WARN("elapsed time: " << elapsedTime << " ms");
    WARN("elapsed steps: " << elapsedSteps);
    WARN("finalParam: " << finalParam.transpose());

    REQUIRE(6 * sizeof(dataT)
          == BinaryFile< std::vector<dataT> >::read(&paramSolution,
              "Weighted_Gauss_Newton_Ref_Grad_tests/gradientUpdateComposeAccumulateParameterOutput.dat"));

    for(int i = 0; i < 6; i++) {
      REQUIRE(paramSolution[i] == Approx(finalParam(i)));
    }
  }
}


