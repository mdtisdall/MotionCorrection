#include "Weighted_Gauss_Newton_Ref_Grad.h"
#include "Weighted_Gauss_Newton_New_Grad.h"

#include "Algorithms.h"

#include "TricubicInterpolator.h"
#include "CubicBSplineInterpolator.h"

#include "BinaryFile.h"

#include <iostream>
#include <fstream>

#include <string>
#include <sstream>

#include <tclap/CmdLine.h>

using namespace Algorithms;

int main(int argc, char* argv[]) {
  typedef TricubicInterpolator<VolumeT, dataT> TricubicInterpolatorT; 
  typedef CubicBSplineInterpolator<VolumeT, dataT> CubicBSplineInterpolatorT; 
  
  size_t cubeSize;
  std::string basePath;
  std::string outputPath;
  dataT translationScaleMM;
  dataT rotationScaleMM;

  try {
    TCLAP::CmdLine cmd("Registering vNav volumes", ' ', "dev");

    TCLAP::ValueArg<std::string> inputFileArg("i", "input", "Path to a volume data, including file name, with the final \"xx.dat\" removed.", true, "", "path", cmd);
    TCLAP::ValueArg<std::string> outputFileArg("o", "output", "Path to output file that will be written.", true, "", "path", cmd);
    TCLAP::ValueArg<unsigned int> widthArg("", "width", "Number of voxels along the side of the vNav", true, 32, "integer", cmd);
    TCLAP::ValueArg<dataT> resolutionArg("", "res", "Size of a voxel's side in mm", true, 8, "mm", cmd);

    cmd.parse(argc, argv);

    basePath = inputFileArg.getValue();
    outputPath = outputFileArg.getValue();

    cubeSize = widthArg.getValue();

    translationScaleMM = resolutionArg.getValue();
    rotationScaleMM = 100.0;

  }   catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    exit(1);
  } 

  const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

  
  ComplexDataCircularMaskOpT fourierMaskOp(cubeSize);
  DataVolumeT maskedRefVolume(cubeSize);
  ComplexVolumeT fourierData(cubeSize);

  
  DataFFTOpT fftOp(cubeSize);
  std::ofstream outputFile(outputPath);

  size_t step = 0;

  for(int baseIndex = 0; baseIndex <= 36; baseIndex += 36) {

    {
      VolumeT refVolume(cubeSize);
  
      std::stringstream refPath;
      refPath << basePath << baseIndex << ".dat";
      
      if(cubeVectorLength * sizeof(dataT)
              != BinaryFile<VolumeT>::read(&refVolume, refPath.str())) {
        std::cerr << "read incorrect length from file: " << refPath.str()
          << std::endl;
        exit(1);
      }
       
      fftOp.forward(&refVolume, &fourierData);

      fourierMaskOp.applyMask(&fourierData); 
          
      fftOp.backward(&fourierData, &maskedRefVolume);
    }

    BinaryFile<VolumeT>::write(&maskedRefVolume,
      "debug_filtered_ref_volume.dat");

    CubicBSplineInterpolatorT cubicBSplineInterpolator(&maskedRefVolume);
  
    CentralDifferencesDifferentiator<VolumeT> volDiffer(&maskedRefVolume);
    VolumeT dx(cubeSize, cubeVectorLength);
    volDiffer.xDerivative(&dx);
  
    VolumeT dy(cubeSize, cubeVectorLength);
    volDiffer.yDerivative(&dy);
  
    VolumeT dz(cubeSize, cubeVectorLength);
    volDiffer.zDerivative(&dz);
    
    CentralDifferencesDifferentiator<VolumeT> dxDiffer(&dx);
   
    VolumeT dxy(cubeSize, cubeVectorLength);
    dxDiffer.yDerivative(&dxy);
    
    VolumeT dxz(cubeSize, cubeVectorLength);
    dxDiffer.zDerivative(&dxz);

    CentralDifferencesDifferentiator<VolumeT> dyDiffer(&dy);

    VolumeT dyz(cubeSize, cubeVectorLength);
    dyDiffer.zDerivative(&dyz);

    CentralDifferencesDifferentiator<VolumeT> dxyDiffer(&dxy);
    
    VolumeT dxyz(cubeSize, cubeVectorLength);
    dxyDiffer.zDerivative(&dxyz);
     
    TricubicInterpolatorT tricubicInterpolator(&maskedRefVolume,
      &dx, &dy, &dz, &dxy, &dxz, &dyz, &dxyz);

    DataCircularMaskOpT imageMaskOp(cubeSize);
  
    Algorithm1<CubicBSplineInterpolatorT> algo1CubicBSplineMinimizer(
      &cubicBSplineInterpolator, &imageMaskOp,
      &dz, &dy, &dx, NULL);
    
    Algorithm2<CubicBSplineInterpolatorT> algo2CubicBSplineMinimizer(
      &cubicBSplineInterpolator, &imageMaskOp,
      &dz, &dy, &dx, NULL);
    
    Algorithm3<CubicBSplineInterpolatorT> algo3CubicBSplineMinimizer(
      &cubicBSplineInterpolator, cubeSize, &imageMaskOp);
    
    Algorithm4<CubicBSplineInterpolatorT> algo4CubicBSplineMinimizer(
      &cubicBSplineInterpolator, cubeSize, &imageMaskOp);
    
    Algorithm5<CubicBSplineInterpolatorT> algo5CubicBSplineMinimizer(
      &cubicBSplineInterpolator, &imageMaskOp,
      &dz, &dy, &dx, NULL);
    
    Algorithm6<CubicBSplineInterpolatorT> algo6CubicBSplineMinimizer(
      &cubicBSplineInterpolator, &imageMaskOp,
      &dz, &dy, &dx, NULL);
    
    Algorithm7<CubicBSplineInterpolatorT> algo7CubicBSplineMinimizer(
      &cubicBSplineInterpolator, cubeSize, &imageMaskOp);
    
    Algorithm8<CubicBSplineInterpolatorT> algo8CubicBSplineMinimizer(
      &cubicBSplineInterpolator, cubeSize, &imageMaskOp);
    
          
    DataVolumeT maskedNewVolume(cubeSize);
    VolumeT newVolume(cubeSize);
    ParamT initialParam;
    ParamT finalParam;
 
    for(unsigned int i = baseIndex + 1; i < baseIndex + 36; i++, step++) { 
      std::stringstream newPath;
      newPath << basePath << i << ".dat";
      
      if(cubeVectorLength * sizeof(dataT)
              != BinaryFile<VolumeT>::read(&newVolume, newPath.str())) {
        std::cerr << "read incorrect length from file: " << newPath.str()
          << std::endl;
        exit(1);
      }
      
      fftOp.forward(&newVolume, &fourierData);
      
      fourierMaskOp.applyMask(&fourierData); 
          
      fftOp.backward(&fourierData, &maskedNewVolume);
      
      initialParam << 0, 0, 0, 0, 0, 0;
      
      const size_t maxSteps = 50;
      const dataT stepSizeScale = 0.25;
      const dataT stepSizeLimit = 1e-5;
  
      const dataT paramUpdateMMLimit = 0.01;
      
      ConvergenceTestT convergenceTest(
        paramUpdateMMLimit,
        translationScaleMM,
        rotationScaleMM);
     
      GradientUpdateTestT gradientUpdateTest;

      std::cout << "----------" << std::endl;
      std::cout << "step " << step << std::endl;
    
      { 
        double elapsedTime;
        size_t elapsedSteps;
       
        algo1CubicBSplineMinimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          &convergenceTest, NULL,
          &elapsedSteps, &elapsedTime);
        
        outputFile << elapsedTime << " " << elapsedSteps
          << " " << finalParam.transpose() << std::endl;
      
        std::cout << "algo 1 cubic b-spline elapsed time: " << elapsedTime << " ms" << std::endl;
        std::cout << "algo 1 cubic b-spline elapsed steps: " << elapsedSteps << std::endl;
      }
      
      { 
        double elapsedTime;
        size_t elapsedSteps;
       

        algo2CubicBSplineMinimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          &convergenceTest, &gradientUpdateTest,
          &elapsedSteps, &elapsedTime);
        
        outputFile << elapsedTime << " " << elapsedSteps
          << " " << finalParam.transpose() << std::endl;
      
        std::cout << "algo 2 cubic b-spline elapsed time: " << elapsedTime << " ms" << std::endl;
        std::cout << "algo 2 cubic b-spline elapsed steps: " << elapsedSteps << std::endl;
      }
       
      CentralDifferencesDifferentiator<VolumeT> newVolDiffer(&maskedNewVolume);
      VolumeT newVolDx(cubeSize, cubeVectorLength);
      newVolDiffer.xDerivative(&newVolDx);
  
      VolumeT newVolDy(cubeSize, cubeVectorLength);
      newVolDiffer.yDerivative(&newVolDy);
  
      VolumeT newVolDz(cubeSize, cubeVectorLength);
      newVolDiffer.zDerivative(&newVolDz);
      
      { 
        double elapsedTime;
        size_t elapsedSteps;
        
        algo3CubicBSplineMinimizer.minimize(&maskedNewVolume,
          &newVolDz, &newVolDy, &newVolDx, 
          &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          &convergenceTest, &gradientUpdateTest,
          &elapsedSteps, &elapsedTime);
        
        outputFile << elapsedTime << " " << elapsedSteps
          << " " << finalParam.transpose() << std::endl;
      
        std::cout << "algo 3 cubic b-spline elapsed time: " << elapsedTime << " ms" << std::endl;
        std::cout << "algo 3 cubic b-spline elapsed steps: " << elapsedSteps << std::endl;
      }

      { 
        double elapsedTime;
        size_t elapsedSteps;
        
        algo4CubicBSplineMinimizer.minimize(&maskedNewVolume,
          &newVolDz, &newVolDy, &newVolDx, 
          &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          &convergenceTest, &gradientUpdateTest,
          &elapsedSteps, &elapsedTime);
        
        outputFile << elapsedTime << " " << elapsedSteps
          << " " << finalParam.transpose() << std::endl;
      
        std::cout << "algo 4 cubic b-spline elapsed time: " << elapsedTime << " ms" << std::endl;
        std::cout << "algo 4 cubic b-spline elapsed steps: " << elapsedSteps << std::endl;
      }

      { 
        double elapsedTime;
        size_t elapsedSteps;
       
        algo5CubicBSplineMinimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          &convergenceTest, NULL,
          &elapsedSteps, &elapsedTime);
        
        outputFile << elapsedTime << " " << elapsedSteps
          << " " << finalParam.transpose() << std::endl;
      
        std::cout << "algo 5 cubic b-spline elapsed time: " << elapsedTime << " ms" << std::endl;
        std::cout << "algo 5 cubic b-spline elapsed steps: " << elapsedSteps << std::endl;
      }
      
      { 
        double elapsedTime;
        size_t elapsedSteps;
       

        algo6CubicBSplineMinimizer.minimize(&maskedNewVolume, &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          &convergenceTest, &gradientUpdateTest,
          &elapsedSteps, &elapsedTime);
        
        outputFile << elapsedTime << " " << elapsedSteps
          << " " << finalParam.transpose() << std::endl;
      
        std::cout << "algo 6 cubic b-spline elapsed time: " << elapsedTime << " ms" << std::endl;
        std::cout << "algo 6 cubic b-spline elapsed steps: " << elapsedSteps << std::endl;
      }

      { 
        double elapsedTime;
        size_t elapsedSteps;
        
        algo7CubicBSplineMinimizer.minimize(&maskedNewVolume,
          &newVolDz, &newVolDy, &newVolDx, 
          &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          &convergenceTest, &gradientUpdateTest,
          &elapsedSteps, &elapsedTime);
        
        outputFile << elapsedTime << " " << elapsedSteps
          << " " << finalParam.transpose() << std::endl;
      
        std::cout << "algo 7 cubic b-spline elapsed time: " << elapsedTime << " ms" << std::endl;
        std::cout << "algo 7 cubic b-spline elapsed steps: " << elapsedSteps << std::endl;
      }

      { 
        double elapsedTime;
        size_t elapsedSteps;
        
        algo8CubicBSplineMinimizer.minimize(&maskedNewVolume,
          &newVolDz, &newVolDy, &newVolDx, 
          &initialParam, &finalParam,
          maxSteps, stepSizeScale, stepSizeLimit,
          &convergenceTest, &gradientUpdateTest,
          &elapsedSteps, &elapsedTime);
        
        outputFile << elapsedTime << " " << elapsedSteps
          << " " << finalParam.transpose() << std::endl;
      
        std::cout << "algo 8 cubic b-spline elapsed time: " << elapsedTime << " ms" << std::endl;
        std::cout << "algo 8 cubic b-spline elapsed steps: " << elapsedSteps << std::endl;
      }
    }
  }

  return 0;
}
