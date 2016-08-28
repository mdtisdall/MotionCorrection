#include "Weighted_Gauss_Newton_Ref_Grad.h"
#include "Weighted_Gauss_Newton_New_Grad.h"

#include "Algorithms.h"

#include "BinaryFile.h"
#include "ReadVolume.h"

#include <iostream>
#include <fstream>

#include <string>
#include <sstream>

#include <tclap/CmdLine.h>

using namespace Algorithms;

template<typename Algo>
void runAlgo(
  Algo *algo,
  VolumeT *newVolume,
  std::ofstream *outputFile,
  std::string algoName) 
{
  ParamT finalParam;
  double elapsedTime;
  size_t elapsedSteps;

  // make a local copy since some of the operations we'll perform may be
  // destructive and we want to keep the original buffer clean.
  VolumeT localNewVolume(newVolume->cubeSize);

  memcpy(localNewVolume.buffer, newVolume->buffer,
    newVolume->totalPoints * sizeof(dataT));

  algo->registerNewVolume(&localNewVolume, &finalParam,
    &elapsedTime, &elapsedSteps);
  
  (*outputFile) << elapsedTime << " " << elapsedSteps
    << " " << finalParam.transpose() << std::endl;

  std::cout << algoName << " cubic b-spline elapsed time: "
    << elapsedTime << " ms" << std::endl;
  std::cout << algoName << " cubic b-spline elapsed steps: "
    << elapsedSteps << std::endl;
}

int main(int argc, char* argv[]) {
  typedef TricubicInterpolator<VolumeT, dataT> TricubicInterpolatorT; 
  typedef CubicBSplineInterpolator<VolumeT, dataT> CubicBSplineInterpolatorT; 
  typedef CentralDifferencesDifferentiator<VolumeT> CentralDiffDifferentiatorT;
  
  size_t cubeSize;
  std::string basePath;
  std::string outputPath;
  dataT translationScaleMM;
  dataT rotationScaleMM;

  try {
    TCLAP::CmdLine cmd("Registering vNav volumes", ' ', "dev");

    TCLAP::ValueArg<std::string> inputFileArg("i", "input", "Path to the first slice of a volume data, including file name, with the final \"_rep_x_slice_x.dat\" removed.", true, "", "path", cmd);
    TCLAP::ValueArg<std::string> outputFileArg("o", "output", "Path to output file that will be written. Output file name should follow \"8mm_bspline_x_rot_3_0_to_5_0_deg_z_trans\"", true, "", "path", cmd);
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
 
  VolumeT refVolume(cubeSize);

  std::ofstream outputFile(outputPath);

  size_t step = 1;

  for(int baseIndex = 0; baseIndex < 36; baseIndex += 36) {

    {
      std::stringstream refPath;
      refPath << basePath << "_rep_" << baseIndex;
      
      if(cubeVectorLength * sizeof(dataT)
              != ReadVolume<VolumeT>::read_volume(&refVolume, refPath.str(), cubeSize)) {
        std::cerr << "read incorrect length from file: " << refPath.str()
          << std::endl;
        exit(1);
      }
    }

    const dataT paramUpdateMMLimit = 0.01;

    MMParamTestT convergenceTest(
      paramUpdateMMLimit,
      translationScaleMM,
      rotationScaleMM);
   
    Algorithm1<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT 
      > algo1CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm2<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT 
      > algo2CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm3<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT 
      > algo3CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm4<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT 
      > algo4CubicBSplineMinimizer( &refVolume, &convergenceTest);

    Algorithm5<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT 
      > algo5CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm6<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT 
      > algo6CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm7<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT 
      > algo7CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm8<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT 
      > algo8CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    VolumeT newVolume(cubeSize);
    ParamT finalParam;

    for(unsigned int i = baseIndex + 1; i < baseIndex + 36; i++, step++) { 
      std::stringstream newPath;
      newPath << basePath << "_rep_" << i;
      
      if(cubeVectorLength * sizeof(dataT)
              != ReadVolume<VolumeT>::read_volume(&newVolume, newPath.str(), cubeSize)) {
        std::cerr << "read incorrect length from file: " << newPath.str()
          << std::endl;
        exit(1);
      }
      std::cout << "----------" << std::endl;
      std::cout << "step " << step << std::endl;
      
      runAlgo(&algo1CubicBSplineMinimizer, &newVolume, &outputFile, "algo1");
      runAlgo(&algo2CubicBSplineMinimizer, &newVolume, &outputFile, "algo2");
      runAlgo(&algo3CubicBSplineMinimizer, &newVolume, &outputFile, "algo3");
      runAlgo(&algo4CubicBSplineMinimizer, &newVolume, &outputFile, "algo4");
      runAlgo(&algo5CubicBSplineMinimizer, &newVolume, &outputFile, "algo5");
      runAlgo(&algo6CubicBSplineMinimizer, &newVolume, &outputFile, "algo6");
      runAlgo(&algo7CubicBSplineMinimizer, &newVolume, &outputFile, "algo7");
      runAlgo(&algo8CubicBSplineMinimizer, &newVolume, &outputFile, "algo8");
    }
  }

  return 0;
}
