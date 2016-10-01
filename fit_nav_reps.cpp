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
  std::string algoName) 
{
  typedef typename Algo::ParamT ParamT;
  typedef typename Algo::VolumeT VolumeT;
  typedef typename VolumeT::value_type dataT;

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
  typedef float dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef std::complex<dataT> complexT;
   
  typedef MMParamTest<dataT> MMParamTestT;
  
  typedef TrilinearInterpolator<VolumeT, dataT> TrilinearInterpolatorT; 
  typedef TricubicInterpolator<VolumeT, dataT> TricubicInterpolatorT; 
  typedef CubicBSplineInterpolator<VolumeT, dataT> CubicBSplineInterpolatorT; 
  
  size_t cubeSize;
  std::string basePath;
  std::string outputPath;
  dataT translationScaleMM;
  dataT rotationScaleMM;

  bool isTrilinear = false;
  bool isTricubic = false;
  bool isCubicBSpline = false;

  try {
    TCLAP::CmdLine cmd("Registering vNav volumes", ' ', "dev");

    TCLAP::ValueArg<std::string> inputFileArg("i", "input", "Path to the first slice of a volume data, including file name, with the final \"_rep_x_slice_x.dat\" removed.", true, "", "path", cmd);
    TCLAP::ValueArg<std::string> outputFileArg("o", "output", "Path to output file that will be written. Output file name should follow \"8mm_bspline_x_rot_3_0_to_5_0_deg_z_trans\"", true, "", "path", cmd);
    TCLAP::ValueArg<unsigned int> widthArg("", "width", "Number of voxels along the side of the vNav", true, 32, "integer", cmd);
    TCLAP::ValueArg<dataT> resolutionArg("", "res", "Size of a voxel's side in mm", true, 8, "mm", cmd);
    TCLAP::SwitchArg linearInterpArg("", "linear", "Use trilinear interpoloation", false);
    TCLAP::SwitchArg cubicInterpArg("", "cubic", "Use tricubic interpoloation", false);
    TCLAP::SwitchArg cubicBSplineInterpArg("", "cubicBSpline", "Use cubic B-spline interpoloation", false);

    std::vector<TCLAP::Arg*>  xorlist;
    xorlist.push_back(&linearInterpArg);
    xorlist.push_back(&cubicInterpArg);
    xorlist.push_back(&cubicBSplineInterpArg);

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
/*   
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
*/    
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
    
    Algorithm10<
      TrilinearInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo10TrilinearMinimizer( &refVolume, &convergenceTest);
    
    Algorithm10<
      TricubicInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo10TricubicMinimizer( &refVolume, &convergenceTest);
    
    Algorithm10<
      CubicBSplineInterpolatorT,
      CentralDiffDifferentiatorT,
      MMParamTestT
      > algo10CubicBSplineMinimizer( &refVolume, &convergenceTest);
    
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
        runAlgo(&algo1TrilinearMinimizer, &newVolume, &outputFile, "algo1");
        runAlgo(&algo2TrilinearMinimizer, &newVolume, &outputFile, "algo2");
        runAlgo(&algo3TrilinearMinimizer, &newVolume, &outputFile, "algo3");
        runAlgo(&algo4TrilinearMinimizer, &newVolume, &outputFile, "algo4");
        runAlgo(&algo5TrilinearMinimizer, &newVolume, &outputFile, "algo5");
        runAlgo(&algo10TrilinearMinimizer, &newVolume, &outputFile, "algo10");
      }
      if(isTricubic) {
        runAlgo(&algo1TricubicMinimizer, &newVolume, &outputFile, "algo1");
        runAlgo(&algo2TricubicMinimizer, &newVolume, &outputFile, "algo2");
        runAlgo(&algo3TricubicMinimizer, &newVolume, &outputFile, "algo3");
        runAlgo(&algo4TricubicMinimizer, &newVolume, &outputFile, "algo4");
        runAlgo(&algo5TricubicMinimizer, &newVolume, &outputFile, "algo5");
        runAlgo(&algo10TricubicMinimizer, &newVolume, &outputFile, "algo10");
      }
      if(isCubicBSpline) {
        runAlgo(&algo1CubicBSplineMinimizer, &newVolume, &outputFile, "algo1");
        runAlgo(&algo2CubicBSplineMinimizer, &newVolume, &outputFile, "algo2");
        runAlgo(&algo3CubicBSplineMinimizer, &newVolume, &outputFile, "algo3");
        runAlgo(&algo4CubicBSplineMinimizer, &newVolume, &outputFile, "algo4");
        runAlgo(&algo5CubicBSplineMinimizer, &newVolume, &outputFile, "algo5");
        runAlgo(&algo10CubicBSplineMinimizer, &newVolume, &outputFile, "algo10");
      }
    }
  }

  return 0;
}
