#ifndef Algorithms_h
#define Algorithms_h

#include "AlgorithmBase.h"

#include "Static_Weighted_Gauss_Newton_New_Grad.h"
#include "Moving_Weighted_Gauss_Newton_Ref_Grad.h"
#include "Moving_Weighted_Gauss_Newton_New_Grad.h"

#include "SumParamAccumulator.h"
#include "ComposeTransformParamAccumulator.h"

#include "Algorithms1And4Base.h"
#include "Algorithms2And5Base.h"
#include "Algorithms3And6Base.h"

namespace Algorithms {


template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm1 : public Algorithm1and4Base <
  InterpolatorT,
  DifferentiatorT,
  ConvergenceTestT,
  SumParamAccumulator< typename InterpolatorT::VolumeT::value_type >
   > {

public:
  typedef Algorithm1and4Base <
    InterpolatorT,
    DifferentiatorT,
    ConvergenceTestT,
    SumParamAccumulator< typename InterpolatorT::VolumeT::value_type >
    > ParentT;
  typedef typename ParentT::VolumeT VolumeT;
  typedef typename ParentT::dataT dataT;
  typedef typename ParentT::ParamT ParamT;
  typedef typename ParentT::WeightFunctionT WeightFunctionT;
  
  Algorithm1(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest) {}
};

template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm4 : public Algorithm1and4Base <
  InterpolatorT,
  DifferentiatorT,
  ConvergenceTestT,
  ComposeTransformParamAccumulator<
    typename InterpolatorT::VolumeT::value_type > > {

public:
  typedef Algorithm1and4Base <
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
  
  Algorithm4(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest) {}
};


template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm2 : public Algorithm2and5Base <
  InterpolatorT,
  DifferentiatorT,
  ConvergenceTestT,
  SumParamAccumulator< typename InterpolatorT::VolumeT::value_type >
   > {

public:
  typedef Algorithm2and5Base <
    InterpolatorT,
    DifferentiatorT,
    ConvergenceTestT,
    SumParamAccumulator< typename InterpolatorT::VolumeT::value_type >
    > ParentT;
  typedef typename ParentT::VolumeT VolumeT;
  typedef typename ParentT::dataT dataT;
  typedef typename ParentT::ParamT ParamT;
  typedef typename ParentT::WeightFunctionT WeightFunctionT;
  
  Algorithm2(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest) {}
};

template <
  typename InterpolatorT,
  typename DifferentiatorT,
  typename ConvergenceTestT
  >
class Algorithm5 : public Algorithm2and5Base <
  InterpolatorT,
  DifferentiatorT,
  ConvergenceTestT,
  ComposeTransformParamAccumulator<
    typename InterpolatorT::VolumeT::value_type > > {

public:
  typedef Algorithm2and5Base <
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
  
  Algorithm5(VolumeT *refVol, ConvergenceTestT *convergenceTest) :
    ParentT(refVol, convergenceTest) {}
};



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



} // end of Algorithms namespace
#endif
