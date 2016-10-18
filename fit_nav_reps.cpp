//#include "Weighted_Gauss_Newton_Ref_Grad.h"
//#include "Weighted_Gauss_Newton_New_Grad.h"

#include "Algorithms.h"

#include "CentralDifferenceDifferentiator.h"
#include "MMParamTest.h"

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
  typename Algo::VolumeT *newVolume,
  std::ofstream *outputFile,
  std::string algoName,
  std::string interpName) 
{
  typedef typename Algo::ParamT ParamT;
  typedef typename Algo::VolumeT VolumeT;
  typedef typename VolumeT::value_type dataT;

  ParamT finalParam;
  double elapsedTime;
  size_t elapsedSteps;
  size_t elapsedSearchSteps;

  // make a local copy since some of the operations we'll perform may be
  // destructive and we want to keep the original buffer clean.
  VolumeT localNewVolume(newVolume->cubeSize);

  memcpy(localNewVolume.buffer, newVolume->buffer,
    newVolume->totalPoints * sizeof(dataT));

  algo->registerNewVolume(&localNewVolume, &finalParam,
    &elapsedTime, &elapsedSteps, &elapsedSearchSteps);
  
  (*outputFile) << elapsedTime << " " << elapsedSteps
    << " " << finalParam.transpose() << std::endl;

  std::cout << algoName << " " << interpName << " elapsed time: "
    << elapsedTime << " ms" << std::endl;
  std::cout << algoName << " " << interpName << " elapsed steps: "
    << elapsedSteps << std::endl;
  std::cout << algoName << " " << interpName << " elapsed search steps: "
    << elapsedSearchSteps << std::endl;
}

int main(int argc, char* argv[]) {
  typedef float dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef std::complex<dataT> complexT;
   
  typedef MMParamTest<dataT> MMParamTestT;
  
  typedef TrilinearInterpolator<VolumeT, dataT> TrilinearInterpolatorT; 
  typedef TricubicInterpolator<VolumeT, dataT> TricubicInterpolatorT; 
  typedef CubicBSplineInterpolator<VolumeT, dataT> CubicBSplineInterpolatorT; 
  typedef UpsampledTrilinearInterpolator<
    CubicBSplineInterpolatorT, 8, VolumeT, dataT>
    CubicBSplineUpsampledTrilinearInterpolatorT; 
  
  size_t cubeSize;
  std::string basePath;
  std::string outputPath;
  dataT translationScaleMM;
  dataT rotationScaleMM;

  bool isTrilinear = false;
  bool isTricubic = false;
  bool isCubicBSpline = false;
  bool isUpsampledTrilinear = false;

  try {
    TCLAP::CmdLine cmd("Registering vNav volumes", ' ', "dev");

    TCLAP::ValueArg<std::string> inputFileArg("i", "input", "Path to the first slice of a volume data, including file name, with the final \"_rep_x_slice_x.dat\" removed.", true, "", "path", cmd);
    TCLAP::ValueArg<std::string> outputFileArg("o", "output", "Path to output file that will be written. Output file name should follow \"8mm_bspline_x_rot_3_0_to_5_0_deg_z_trans\"", true, "", "path", cmd);
    TCLAP::ValueArg<unsigned int> widthArg("", "width", "Number of voxels along the side of the vNav", true, 32, "integer", cmd);
    TCLAP::ValueArg<dataT> resolutionArg("", "res", "Size of a voxel's side in mm", true, 8, "mm", cmd);
    TCLAP::SwitchArg linearInterpArg("", "linear", "Use trilinear interpoloation", false);
    TCLAP::SwitchArg cubicInterpArg("", "cubic", "Use tricubic interpoloation", false);
    TCLAP::SwitchArg cubicBSplineInterpArg("", "cubicBSpline", "Use cubic B-spline interpoloation", false);
    TCLAP::SwitchArg upsampledLinearInterpArg("", "upsampled-linear", "Use cubic B-spline upsampled data with trilinear interpoloation", false);

    std::vector<TCLAP::Arg*>  xorlist;
    xorlist.push_back(&linearInterpArg);
    xorlist.push_back(&cubicInterpArg);
    xorlist.push_back(&cubicBSplineInterpArg);
    xorlist.push_back(&upsampledLinearInterpArg);

    cmd.xorAdd( xorlist );

    cmd.parse(argc, argv);

    basePath = inputFileArg.getValue();
    outputPath = outputFileArg.getValue();

    cubeSize = widthArg.getValue();

    translationScaleMM = resolutionArg.getValue();
    rotationScaleMM = 100.0;

    isTrilinear = linearInterpArg.getValue();
    isTricubic = cubicInterpArg.getValue();
    isCubicBSpline = cubicBSplineInterpArg.getValue();
    isUpsampledTrilinear = upsampledLinearInterpArg.getValue();

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
    
    typedef CentralDifferencesDifferentiator<VolumeT>
      CentralDiffDifferentiatorT;
    
    Algorithm1<
      TrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo1TrilinearMinimizer( &refVolume, &convergenceTest);
     
    Algorithm1<
      TricubicInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo1TricubicMinimizer( &refVolume, &convergenceTest);
    
    Algorithm1<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo1CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm1<
      CubicBSplineUpsampledTrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo1UpsampledTrilinearMinimizer( &refVolume, &convergenceTest);
     
    Algorithm2<
      TrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo2TrilinearMinimizer( &refVolume, &convergenceTest);
     
    Algorithm2<
      TricubicInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo2TricubicMinimizer( &refVolume, &convergenceTest);
    
    Algorithm2<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo2CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm2<
      CubicBSplineUpsampledTrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo2UpsampledTrilinearMinimizer( &refVolume, &convergenceTest);
     
    Algorithm3<
      TrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo3TrilinearMinimizer( &refVolume, &convergenceTest);
     
    Algorithm3<
      TricubicInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo3TricubicMinimizer( &refVolume, &convergenceTest);
    
    Algorithm3<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo3CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm3<
      CubicBSplineUpsampledTrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo3UpsampledTrilinearMinimizer( &refVolume, &convergenceTest);
    
    Algorithm4<
      TrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo4TrilinearMinimizer( &refVolume, &convergenceTest);
     
    Algorithm4<
      TricubicInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo4TricubicMinimizer( &refVolume, &convergenceTest);
    
    Algorithm4<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo4CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm4<
      CubicBSplineUpsampledTrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo4UpsampledTrilinearMinimizer( &refVolume, &convergenceTest);
     
    Algorithm5<
      TrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo5TrilinearMinimizer( &refVolume, &convergenceTest);
     
    Algorithm5<
      TricubicInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo5TricubicMinimizer( &refVolume, &convergenceTest);
    
    Algorithm5<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo5CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm5<
      CubicBSplineUpsampledTrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo5UpsampledTrilinearMinimizer( &refVolume, &convergenceTest);
    
    Algorithm6<
      TrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo6TrilinearMinimizer( &refVolume, &convergenceTest);
    
    Algorithm6<
      TricubicInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo6TricubicMinimizer( &refVolume, &convergenceTest);
    
    Algorithm6<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo6CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
    Algorithm6<
      CubicBSplineUpsampledTrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo6UpsampledTrilinearMinimizer( &refVolume, &convergenceTest);
    
    VolumeT newVolume(cubeSize);

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

      if(isTrilinear) {
        runAlgo(&algo1TrilinearMinimizer, &newVolume, &outputFile, "algo1", "trilinear");
        runAlgo(&algo2TrilinearMinimizer, &newVolume, &outputFile, "algo2", "trilinear");
        runAlgo(&algo3TrilinearMinimizer, &newVolume, &outputFile, "algo3", "trilinear");
        runAlgo(&algo4TrilinearMinimizer, &newVolume, &outputFile, "algo4", "trilinear");
        runAlgo(&algo5TrilinearMinimizer, &newVolume, &outputFile, "algo5", "trilinear");
        runAlgo(&algo6TrilinearMinimizer, &newVolume, &outputFile, "algo6", "trilinear");
      }
      if(isTricubic) {
        runAlgo(&algo1TricubicMinimizer, &newVolume, &outputFile, "algo1", "tricubic");
        runAlgo(&algo2TricubicMinimizer, &newVolume, &outputFile, "algo2", "tricubic");
        runAlgo(&algo3TricubicMinimizer, &newVolume, &outputFile, "algo3", "tricubic");
        runAlgo(&algo4TricubicMinimizer, &newVolume, &outputFile, "algo4", "tricubic");
        runAlgo(&algo5TricubicMinimizer, &newVolume, &outputFile, "algo5", "tricubic");
        runAlgo(&algo6TricubicMinimizer, &newVolume, &outputFile, "algo6", "tricubic");
      }
      if(isCubicBSpline) {
        runAlgo(&algo1CubicBSplineMinimizer, &newVolume, &outputFile, "algo1", "cubic b-spline");
        runAlgo(&algo2CubicBSplineMinimizer, &newVolume, &outputFile, "algo2", "cubic b-spline");
        runAlgo(&algo3CubicBSplineMinimizer, &newVolume, &outputFile, "algo3", "cubic b-spline");
        runAlgo(&algo4CubicBSplineMinimizer, &newVolume, &outputFile, "algo4", "cubic b-spline");
        runAlgo(&algo5CubicBSplineMinimizer, &newVolume, &outputFile, "algo5", "cubic b-spline");
        runAlgo(&algo6CubicBSplineMinimizer, &newVolume, &outputFile, "algo6", "cubic b-spline");
      }
      if(isUpsampledTrilinear) {
        runAlgo(&algo1UpsampledTrilinearMinimizer, &newVolume, &outputFile, "algo1", "upsampled trilinear");
        runAlgo(&algo2UpsampledTrilinearMinimizer, &newVolume, &outputFile, "algo2", "upsampled trilinear");
        runAlgo(&algo3UpsampledTrilinearMinimizer, &newVolume, &outputFile, "algo3", "upsampled trilinear");
        runAlgo(&algo4UpsampledTrilinearMinimizer, &newVolume, &outputFile, "algo4", "upsampled trilinear");
        runAlgo(&algo5UpsampledTrilinearMinimizer, &newVolume, &outputFile, "algo5", "upsampled trilinear");
        runAlgo(&algo6UpsampledTrilinearMinimizer, &newVolume, &outputFile, "algo6", "upsampled trilinear");
      }
    }
  }

  return 0;
}
