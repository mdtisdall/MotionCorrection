#ifndef Algorithms1And4Base_h
#define Algorithms1And4Base_h

namespace Algorithms {

template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT,
  typename ParamAccumulatorT
  >
class Algorithm1and4Base : public AlgorithmBase <
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
  typedef typename ParentT::WeightGradientFunctionT WeightGradientFunctionT;
  
  typedef Moving_Weighted_Gauss_Newton_Ref_Grad <
    InterpolatorT,
    ParamAccumulatorT,
    WeightFunctionT,
    WeightGradientFunctionT,
    ConvergenceTestT > MinimizerT; 

  Algorithm1and4Base(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest),
    refVolDx(refVol->cubeSize),
    refVolDy(refVol->cubeSize),
    refVolDz(refVol->cubeSize),
    minimizer(
      createMinimizer(refVol,
      &refVolDz,
      &refVolDy,
      &refVolDx,
      this->interpolator,
      &(this->weightFunction),
      &(this->weightGradientFunction))) {}

protected:
  virtual void registerNewVolumeInner( 
    VolumeT *newVol,
    ParamT *p, 
    size_t *elapsedSteps,
    size_t *elapsedSearchSteps) {
    
    ParamT initialParam; 
    initialParam << 0, 0, 0, 0, 0, 0;
     
    minimizer.minimize(
      newVol,
      &initialParam, p,
      this->maxSteps, this->stepSizeScale,
      this->convergenceTest,
      elapsedSteps, elapsedSearchSteps, NULL);
  }

  static MinimizerT createMinimizer(
    VolumeT *refVol,
    VolumeT *refVolDz,
    VolumeT *refVolDy,
    VolumeT *refVolDx,
    InterpolatorT *interpolator,
    WeightFunctionT *weightFunction,
    WeightGradientFunctionT *weightGradientFunction) {
    DifferentiatorT refVolDiffer(refVol);
    refVolDiffer.xDerivative(refVolDx);
    refVolDiffer.yDerivative(refVolDy);
    refVolDiffer.zDerivative(refVolDz);
    
    return MinimizerT(
      interpolator,
      refVolDz, refVolDy, refVolDx,
      weightFunction, weightGradientFunction);
  }

  VolumeT refVolDx;
  VolumeT refVolDy;
  VolumeT refVolDz;
  MinimizerT minimizer;

};

} // Algorithms namespace

#endif

