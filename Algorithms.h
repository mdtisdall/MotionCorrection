#ifndef Algorithms_h
#define Algorithms_h

#include "AlgorithmBase.h"

#include "Static_Weighted_Gauss_Newton_New_Grad.h"
#include "Moving_Weighted_Gauss_Newton_Ref_Grad.h"
#include "Moving_Weighted_Gauss_Newton_New_Grad.h"

#include "SumParamAccumulator.h"
#include "ComposeTransformParamAccumulator.h"

#include "Algorithms3And6Base.h"

namespace Algorithms {


template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm3 : public Algorithm3and6Base <
  InterpolatorT,
  DifferentiatorT,
  ConvergenceTestT,
  SumParamAccumulator< typename InterpolatorT::VolumeT::value_type > > {

public:
  typedef Algorithm3and6Base <
    InterpolatorT,
    DifferentiatorT,
    ConvergenceTestT,
    SumParamAccumulator< typename InterpolatorT::VolumeT::value_type >
    > ParentT;
  typedef typename ParentT::VolumeT VolumeT;
  typedef typename ParentT::dataT dataT;
  typedef typename ParentT::ParamT ParamT;
  typedef typename ParentT::WeightFunctionT WeightFunctionT;

  Algorithm3(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest) {}
};


template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm6 : public Algorithm3and6Base <
  InterpolatorT,
  DifferentiatorT,
  ConvergenceTestT,
  ComposeTransformParamAccumulator<
    typename InterpolatorT::VolumeT::value_type > > {

public:
  typedef Algorithm3and6Base <
    InterpolatorT,
    DifferentiatorT,
    ConvergenceTestT,
    ComposeTransformParamAccumulator<
      typename InterpolatorT::VolumeT::value_type >
    > ParentT;
  typedef typename ParentT::VolumeT VolumeT;
  typedef typename ParentT::dataT dataT;
  typedef typename ParentT::ParamT ParamT;
  typedef typename ParentT::WeightFunctionT WeightFunctionT;
  
  Algorithm6(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest) {}
};


template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm1 : public AlgorithmBase <
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
  typedef SumParamAccumulator<dataT> ParamAccumulatorT; 

  typedef Moving_Weighted_Gauss_Newton_Ref_Grad <
    InterpolatorT,
    ParamAccumulatorT,
    WeightFunctionT,
    WeightGradientFunctionT,
    ConvergenceTestT > MinimizerT; 

  Algorithm1(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest),
    minimizer(
      createMinimizer(refVol,
      this->interpolator,
      &(this->weightFunction),
      &(this->weightGradientFunction))) {}

protected:
  virtual void registerNewVolumeInner( 
    VolumeT *newVol,
    ParamT *p, 
    size_t *elapsedSteps) {
    
    ParamT initialParam; 
    initialParam << 0, 0, 0, 0, 0, 0;
     
    minimizer.minimize(
      newVol,
      &initialParam, p,
      this->maxSteps, this->stepSizeScale, this->stepSizeLimit,
      this->convergenceTest,
      elapsedSteps, NULL);
  }

  static MinimizerT createMinimizer(
    VolumeT *refVol,
    InterpolatorT *interpolator,
    WeightFunctionT *weightFunction,
    WeightGradientFunctionT *weightGradientFunction) {
    DifferentiatorT refVolDiffer(refVol);
    VolumeT refVolDx(refVol->cubeSize);
    refVolDiffer.xDerivative(&refVolDx);
  
    VolumeT refVolDy(refVol->cubeSize);
    refVolDiffer.yDerivative(&refVolDy);
  
    VolumeT refVolDz(refVol->cubeSize);
    refVolDiffer.zDerivative(&refVolDz);
    
    return MinimizerT(
      interpolator,
      &refVolDz, &refVolDy, &refVolDx,
      weightFunction, weightGradientFunction);
  }

  MinimizerT minimizer;

};




template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm2 : public AlgorithmBase <
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
  typedef SumParamAccumulator<dataT> ParamAccumulatorT; 

  typedef Moving_Weighted_Gauss_Newton_New_Grad <
    InterpolatorT,
    ParamAccumulatorT,
    WeightFunctionT,
    WeightGradientFunctionT,
    ConvergenceTestT > MinimizerT; 

  Algorithm2(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest),
    minimizer( 
      this->interpolator,
      this->cubeSize,
      &(this->weightFunction),
      &(this->weightGradientFunction)) {}

protected:
  virtual void registerNewVolumeInner( 
    VolumeT *newVol,
    ParamT *p, 
    size_t *elapsedSteps) {
    
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
      elapsedSteps, NULL);
  }

  MinimizerT minimizer;

};


template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm4 : public AlgorithmBase <
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
  typedef ComposeTransformParamAccumulator<
    typename InterpolatorT::VolumeT::value_type > ParamAccumulatorT; 

  typedef Moving_Weighted_Gauss_Newton_Ref_Grad <
    InterpolatorT,
    ParamAccumulatorT,
    WeightFunctionT,
    WeightGradientFunctionT,
    ConvergenceTestT > MinimizerT; 

  Algorithm4(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest),
    minimizer(
      createMinimizer(refVol,
      this->interpolator,
      &(this->weightFunction),
      &(this->weightGradientFunction))) {}

protected:
  virtual void registerNewVolumeInner( 
    VolumeT *newVol,
    ParamT *p, 
    size_t *elapsedSteps) {
    
    ParamT initialParam; 
    initialParam << 0, 0, 0, 0, 0, 0;
     
    minimizer.minimize(
      newVol,
      &initialParam, p,
      this->maxSteps, this->stepSizeScale, this->stepSizeLimit,
      this->convergenceTest,
      elapsedSteps, NULL);
  }

  static MinimizerT createMinimizer(
    VolumeT *refVol,
    InterpolatorT *interpolator,
    WeightFunctionT *weightFunction,
    WeightGradientFunctionT *weightGradientFunction) {
    DifferentiatorT refVolDiffer(refVol);
    VolumeT refVolDx(refVol->cubeSize);
    refVolDiffer.xDerivative(&refVolDx);
  
    VolumeT refVolDy(refVol->cubeSize);
    refVolDiffer.yDerivative(&refVolDy);
  
    VolumeT refVolDz(refVol->cubeSize);
    refVolDiffer.zDerivative(&refVolDz);
    
    return MinimizerT(
      interpolator,
      &refVolDz, &refVolDy, &refVolDx,
      weightFunction, weightGradientFunction);
  }

  MinimizerT minimizer;

};




template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm5 : public AlgorithmBase <
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
  typedef ComposeTransformParamAccumulator<
    typename InterpolatorT::VolumeT::value_type > ParamAccumulatorT; 

  typedef Moving_Weighted_Gauss_Newton_New_Grad <
    InterpolatorT,
    ParamAccumulatorT,
    WeightFunctionT,
    WeightGradientFunctionT,
    ConvergenceTestT > MinimizerT; 

  Algorithm5(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest),
    minimizer( 
      this->interpolator,
      this->cubeSize,
      &(this->weightFunction),
      &(this->weightGradientFunction)) {}

protected:
  virtual void registerNewVolumeInner( 
    VolumeT *newVol,
    ParamT *p, 
    size_t *elapsedSteps) {
    
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
      elapsedSteps, NULL);
  }

  MinimizerT minimizer;

};



} // end of Algorithms namespace
#endif
