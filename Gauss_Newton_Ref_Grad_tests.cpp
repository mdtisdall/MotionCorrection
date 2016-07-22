#include "catch.hpp"

#include "Gauss_Newton_Ref_Grad.h"
//#include "Gauss_Newton_New_Grad.h"

#include "CubicBSplineInterpolator.h"
#include "VolumeAtAddressable.h"

#include "TwoNormParamTest.h"
#include "TrueParamTest.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTWBuffer.h"

#include "BinaryFile.h"

TEST_CASE("a Gauss-Newton minimizer using reference-image gradients can be instantiated") {
  typedef float dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef CubicBSplineInterpolator<VolumeT, dataT> InterpolatorT; 

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
    typedef TwoNormParamTest<dataT> GradientUpdateTestT;
    typedef Gauss_Newton_Ref_Grad<
      InterpolatorT, ConvergenceTestT, GradientUpdateTestT > MinimizerT; 
    typedef MinimizerT::ParamT ParamT;
    
    MinimizerT minimizer(&interpolator, &dz, &dy, &dx,
      &gradientAndHessianComputeTime);
      
    WARN("elapsed time computing gradient and Hessian: "
      << gradientAndHessianComputeTime << " ms");

    SECTION("and registering an image with itself produces 0 transformation") {
      ParamT initialParam;
      initialParam << 0, 0, 0, 0, 0, 0;
  
      ParamT finalParam;
 
      size_t maxSteps = 20;
      dataT stepSizeScale = 0.25;
      dataT stepSizeLimit = 1.0;

      ConvergenceTestT convergenceTest(1e-6);

      double elapsedTime;
      size_t elapsedSteps;
  
      minimizer.minimize(&volume, &initialParam, &finalParam,
        maxSteps, stepSizeScale, stepSizeLimit,
        &convergenceTest, NULL,
        &elapsedSteps, &elapsedTime);
 
        for(int i = 0; i < 6; i++) {
          REQUIRE(0 == Approx(finalParam(i)));
        }

      WARN("elapsed time: " << elapsedTime << " ms");
      WARN("elapsed steps: " << elapsedSteps);
      WARN("finalParam: " << finalParam.transpose());
    } 
}


TEST_CASE("a Gauss-Newton minimizer using reference-image gradients can be instantiated from image data") {
    typedef double dataT;
    typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
    typedef CubicBSplineInterpolator<VolumeT, dataT> InterpolatorT; 
    typedef TwoNormParamTest<dataT> ConvergenceTestT;
    typedef TrueParamTest<dataT> GradientUpdateTestT;
    typedef Gauss_Newton_Ref_Grad<
      InterpolatorT, ConvergenceTestT, GradientUpdateTestT > MinimizerT; 
    typedef MinimizerT::ParamT ParamT;

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    VolumeT refVolume(cubeSize);
    VolumeT newVolume(cubeSize);

    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&refVolume,
                "Gauss_Newton_Ref_Grad_tests/refVolInput.dat"));
    
    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<VolumeT>::read(&newVolume,
                "Gauss_Newton_Ref_Grad_tests/newVolInput.dat"));

    InterpolatorT interpolator(&refVolume);


      CentralDifferencesDifferentiator<VolumeT> volDiffer(&refVolume);
      VolumeT dx(cubeSize, cubeVectorLength);
      volDiffer.xDerivative(&dx);

      VolumeT dy(cubeSize, cubeVectorLength);
      volDiffer.yDerivative(&dy);

      VolumeT dz(cubeSize, cubeVectorLength);
      volDiffer.zDerivative(&dz);
     
      double gradientAndHessianComputeTime;

      MinimizerT minimizer(&interpolator, &dz, &dy, &dx,
        &gradientAndHessianComputeTime);
        
      WARN("elapsed time computing gradient and Hessian: "
        << gradientAndHessianComputeTime << " ms");

      SECTION(std::string("and registering two images ") +
        std::string("without updating gradients ") + 
        std::string("returns identical result to Mathematica")) {
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
  
        minimizer.minimize(&newVolume, &initialParam, &finalParam,
        maxSteps, stepSizeScale, stepSizeLimit,
        NULL,  NULL,
        &elapsedSteps, &elapsedTime);

        std::vector<dataT> paramSolution(6);

        REQUIRE(6 * sizeof(dataT)
              == BinaryFile< std::vector<dataT> >::read(&paramSolution,
                  "Gauss_Newton_Ref_Grad_tests/parameterOutput.dat"));

        for(int i = 0; i < 6; i++) {
          REQUIRE(paramSolution[i] == Approx(finalParam(i)));
        }

        WARN("elapsed time: " << elapsedTime << " ms");
        WARN("elapsed steps: " << elapsedSteps);
      }

      SECTION(std::string("and registering two images ") +
        std::string("updating gradients every step ") + 
        std::string("returns identical result to Mathematica")) {
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
 
        GradientUpdateTestT gradientUpdateTest;

        minimizer.minimize(&newVolume, &initialParam, &finalParam,
        maxSteps, stepSizeScale, stepSizeLimit,
        NULL,  &gradientUpdateTest,
        &elapsedSteps, &elapsedTime);

        std::vector<dataT> paramSolution(6);

        REQUIRE(6 * sizeof(dataT)
              == BinaryFile< std::vector<dataT> >::read(&paramSolution,
                  "Gauss_Newton_Ref_Grad_tests/gradientUpdateParameterOutput.dat"));

        for(int i = 0; i < 6; i++) {
          REQUIRE(paramSolution[i] == Approx(finalParam(i)));
        }

        WARN("elapsed time: " << elapsedTime << " ms");
        WARN("elapsed steps: " << elapsedSteps);
      }
}

