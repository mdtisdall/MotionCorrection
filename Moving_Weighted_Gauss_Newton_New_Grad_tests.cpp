#include "catch.hpp"

#include "Moving_Weighted_Gauss_Newton_New_Grad.h"

#include "CubicBSplineInterpolator.h"
#include "VolumeAtAddressable.h"

#include "TwoNormParamTest.h"
#include "TrueParamTest.h"

#include "SumParamAccumulator.h"
#include "ComposeTransformParamAccumulator.h"

#include "WeightFunction.h"
#include "DerivWeightFunction.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTOp.h"
#include "FFTWBuffer.h"

#include "CircularMaskOp.h"

#include "BinaryFile.h"


TEST_CASE("a Gauss-Newton minimizer using moving weights and new-image gradients can be instantiated") {
  typedef float dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef CubicBSplineInterpolator<VolumeT, float> InterpolatorT; 

  const size_t cubeSize = 32;
  const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

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
  
  typedef TwoNormParamTest<dataT> ConvergenceTestT;
  typedef SumParamAccumulator<dataT> ParamAccumulatorT;
  typedef WeightFunction<dataT> WeightFuncT;
  typedef DerivWeightFunction<dataT> WeightGradientFuncT;
  typedef Moving_Weighted_Gauss_Newton_New_Grad <
    InterpolatorT,
    ParamAccumulatorT,
    WeightFuncT,
    WeightGradientFuncT,
    ConvergenceTestT > MinimizerT; 
  typedef MinimizerT::ParamT ParamT;
  
  WeightFuncT weightFunc(cubeSize);
  WeightGradientFuncT weightGradientFunc(cubeSize);

  MinimizerT minimizer(
    &interpolator, cubeSize, &weightFunc, &weightGradientFunc);

  SECTION("and registering an image with itself produces 0 transformation") {
    ParamT initialParam;
    initialParam << 0, 0, 0, 0, 0, 0;
  
    ParamT finalParam;
 
    size_t maxSteps = 20;
    float stepSizeScale = 0.25;
    float stepSizeLimit = 1.0;


    const dataT paramUpdate2NormLimit = 1e-6;
    
    ConvergenceTestT convergenceTest(paramUpdate2NormLimit);
  
    double elapsedTime;
    size_t elapsedSteps;
  
    minimizer.minimize(&volume, &dz, &dy, &dx,
      &initialParam, &finalParam,
      maxSteps, stepSizeScale, stepSizeLimit,
      &convergenceTest,
      &elapsedSteps, &elapsedTime, &gradientAndHessianComputeTime);
    
    WARN("elapsed time computing gradient and Hessian: "
      << gradientAndHessianComputeTime << " ms");

 
      for(int i = 0; i < 6; i++) {
        REQUIRE(0 == Approx(finalParam(i)));
      }

    WARN("elapsed time: " << elapsedTime << " ms");
    WARN("elapsed steps: " << elapsedSteps);
    WARN("finalParam: " << finalParam.transpose());
  }


}

TEST_CASE("a Gauss-Newton minimizer using moving weights, moving M, and new-image gradients can be instantiated from image data") {
    typedef double dataT;
    typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
    typedef CubicBSplineInterpolator<VolumeT, dataT> InterpolatorT; 

  const size_t cubeSize = 32;
  const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

  typedef std::complex<dataT> complexT;
  typedef FFTOp<dataT> DataFFTOpT;
  typedef DataFFTOpT::spatialVolumeT DataVolumeT; 
  typedef DataFFTOpT::fourierVolumeT ComplexVolumeT; 
  typedef SymmetricHalfVolumeAtAddressable< FFTWBuffer<dataT> >
    FourierMaskVolumeT;
  typedef WeightFunction<dataT> WeightFunctionT;
  typedef CircularMaskOp< ComplexVolumeT, FourierMaskVolumeT, WeightFunctionT>
    ComplexDataCircularMaskOpT;
  
  const dataT maskScale =
    1.0/((dataT) cubeVectorLength); 
 
  WeightFunctionT weightFunction(cubeSize);

  ComplexDataCircularMaskOpT
    fourierMaskOp(cubeSize, &weightFunction, maskScale);
  DataVolumeT maskedRefVolume(cubeSize);
  DataVolumeT maskedNewVolume(cubeSize);

  {
    VolumeT refVolume(cubeSize);
    VolumeT newVolume(cubeSize);
  
    ComplexVolumeT fourierData(cubeSize);
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&refVolume,
                "Moving_Weighted_Gauss_Newton_New_Grad_tests/refVolInput.dat"));
    
    DataFFTOpT fftOp(cubeSize);
  
    fftOp.forward(&refVolume, &fourierData);
    
    fourierMaskOp.applyMask(&fourierData); 
        
    fftOp.backward(&fourierData, &maskedRefVolume);
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&newVolume,
                "Moving_Weighted_Gauss_Newton_New_Grad_tests/newVolInput.dat"));
    
    fftOp.forward(&newVolume, &fourierData);
    
    fourierMaskOp.applyMask(&fourierData); 
        
    fftOp.backward(&fourierData, &maskedNewVolume);
  }

//  std::cout << "maskedNewVolume.at(0): " << maskedNewVolume.at(0) << std::endl;
//  std::cout << "maskedRefVolume.at(0): " << maskedRefVolume.at(0) << std::endl;

  InterpolatorT interpolator(&maskedRefVolume);

  CentralDifferencesDifferentiator<VolumeT> volDiffer(&maskedNewVolume);
  VolumeT dx(cubeSize, cubeVectorLength);
  volDiffer.xDerivative(&dx);

  VolumeT dy(cubeSize, cubeVectorLength);
  volDiffer.yDerivative(&dy);

  VolumeT dz(cubeSize, cubeVectorLength);
  volDiffer.zDerivative(&dz);
  
  double gradientAndHessianComputeTime;
  
  WARN("elapsed time computing gradient and Hessian: "
    << gradientAndHessianComputeTime << " ms");

        
      SECTION(std::string("and registering two images ") +
        std::string("returns identical result to Mathematica")) {
    
        typedef SumParamAccumulator<dataT> ParamAccumulatorT;
        typedef WeightFunction<dataT> WeightFuncT;
        typedef DerivWeightFunction<dataT> WeightGradientFuncT;
        typedef Moving_Weighted_Gauss_Newton_New_Grad<
          InterpolatorT, ParamAccumulatorT, WeightFuncT, WeightGradientFuncT>
          MinimizerT;
        typedef MinimizerT::ParamT ParamT;
      
        WeightFuncT weightFunc(cubeSize);
        WeightGradientFuncT weightGradientFunc(cubeSize);

        MinimizerT minimizer(
          &interpolator, cubeSize, &weightFunc, &weightGradientFunc);
        
        ParamT initialParam;
        initialParam << 0, 0, 0, 0, 0, 0;
  
        ParamT finalParam;

        size_t maxSteps = 20;
        dataT stepSizeScale = 0.25;
        dataT stepSizeLimit = 1.0;
 
        double elapsedTime;
        size_t elapsedSteps;
  
        minimizer.minimize(&maskedNewVolume, &dz, &dy, &dx,
          &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          NULL,
          &elapsedSteps, &elapsedTime, &gradientAndHessianComputeTime);
      
        WARN("elapsed time computing gradient and Hessian: "
          << gradientAndHessianComputeTime << " ms");


        std::vector<dataT> paramSolution(6);
        
        INFO("finalParam: " << finalParam.transpose());


        REQUIRE(6 * sizeof(dataT)
              == BinaryFile< std::vector<dataT> >::read(&paramSolution,
                  "Moving_Weighted_Gauss_Newton_New_Grad_tests/parameterOutput.dat"));

        for(int i = 0; i < 6; i++) {
          REQUIRE(paramSolution[i] == Approx(finalParam(i)));
        }

        WARN("elapsed time: " << elapsedTime << " ms");
        WARN("elapsed steps: " << elapsedSteps);
      }

      SECTION(std::string("and registering two images ") +
        std::string("using transform composition accumulator ") +
        std::string("returns identical result to Mathematica")) {
    
        typedef ComposeTransformParamAccumulator<dataT> ParamAccumulatorT;
        typedef WeightFunction<dataT> WeightFuncT;
        typedef DerivWeightFunction<dataT> WeightGradientFuncT;
        typedef Moving_Weighted_Gauss_Newton_New_Grad<
          InterpolatorT, ParamAccumulatorT, WeightFuncT, WeightGradientFuncT>
          MinimizerT; 
        typedef MinimizerT::ParamT ParamT;
      
        WeightFuncT weightFunc(cubeSize);
        WeightGradientFuncT weightGradientFunc(cubeSize);

        MinimizerT minimizer(
          &interpolator, cubeSize, &weightFunc, &weightGradientFunc);
        
        ParamT initialParam;
        initialParam << 0, 0, 0, 0, 0, 0;
  
        ParamT finalParam;

        // We force this to go exactly 20 steps, so that we get to the same
        // point as the Mathematica code
        size_t maxSteps = 20;
        dataT stepSizeScale = 0.25;
        dataT stepSizeLimit = 1.0;
 
        double elapsedTime;
        size_t elapsedSteps;
  
        minimizer.minimize(&maskedNewVolume, &dz, &dy, &dx,
          &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          NULL,
          &elapsedSteps, &elapsedTime, &gradientAndHessianComputeTime);
      
        WARN("elapsed time computing gradient and Hessian: "
          << gradientAndHessianComputeTime << " ms");


        std::vector<dataT> paramSolution(6);

        INFO("finalParam: " << finalParam.transpose());

        REQUIRE(6 * sizeof(dataT)
              == BinaryFile< std::vector<dataT> >::read(&paramSolution,
                  "Moving_Weighted_Gauss_Newton_New_Grad_tests/accumulateParameterOutput.dat"));

        for(int i = 0; i < 6; i++) {
          REQUIRE(paramSolution[i] == Approx(finalParam(i)));
        }

        WARN("elapsed time: " << elapsedTime << " ms");
        WARN("elapsed steps: " << elapsedSteps);
      }

}

