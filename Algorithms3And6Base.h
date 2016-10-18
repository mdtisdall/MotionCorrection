#ifndef Algorithms3And6Base_h
#define Algorithms3And6Base_h

namespace Algorithms {

template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT,
  typename ParamAccumulatorT
  >
class Algorithm3and6Base : public AlgorithmBase <
  InterpolatorT,
  ConvergenceTestT
  > {
public:
  typedef AlgorithmBase < InterpolatorT,
    ConvergenceTestT
    > ParentT;
  typedef typename ParentT::VolumeT VolumeT;
  typedef typename ParentT::dataT dataT;
  typedef typename ParentT::ParamT ParamT;
  typedef typename ParentT::WeightFunctionT WeightFunctionT;
  
  typedef Static_Weighted_Gauss_Newton_New_Grad <
    InterpolatorT,
    ParamAccumulatorT,
    WeightFunctionT,
    ConvergenceTestT > MinimizerT; 

  Algorithm3and6Base(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest),
    minimizer(this->interpolator, this->cubeSize, &(this->weightFunction)) {
    
  }
protected:
  virtual void registerNewVolumeInner( 
    VolumeT *newVol,
    ParamT *p, 
    size_t *elapsedSteps,
    size_t *elapsedSearchSteps) {
    
    ParamT initialParam; 
    initialParam << 0, 0, 0, 0, 0, 0;
    
    DifferentiatorT newVolDiffer(newVol);
    VolumeT newVolDx(this->cubeSize);
    newVolDiffer.xDerivative(&newVolDx);
  
    VolumeT newVolDy(this->cubeSize);
    newVolDiffer.yDerivative(&newVolDy);
  
    VolumeT newVolDz(this->cubeSize);
    newVolDiffer.zDerivative(&newVolDz);
    
    minimizer.minimize(
      newVol,
      &newVolDz, &newVolDy, &newVolDx, 
      &initialParam, p,
      this->maxSteps, this->stepSizeScale, this->stepSizeLimit,
      this->convergenceTest,
      elapsedSteps, elapsedSearchSteps, NULL);
  }

  MinimizerT minimizer;

};

} // Algorithms namespace

#endif

