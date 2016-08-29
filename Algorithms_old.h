#ifndef Algorithms_h
#define Algorithms_h

#include "VolumeAtAddressable.h"

#include "CircularMaskOp.h"

#include "SumParamAccumulator.h"
#include "ComposeTransformParamAccumulator.h"

#include "TrueParamTest.h"
#include "MMParamTest.h"

#include "CentralDifferenceDifferentiator.h"

#include "FFTOp.h"

#include "Algorithm.h"

#include <Eigen/Dense>

namespace Algorithms {
  typedef float dataT;
  typedef VolumeAtAddressable< FFTWBuffer<dataT> > VolumeT; 
  typedef std::complex<dataT> complexT;
  
  typedef SumParamAccumulator<dataT> SumParamAccumulatorT;
  typedef ComposeTransformParamAccumulator<dataT>
    ComposeTransformParamAccumulatorT;
  
  typedef TrueParamTest<dataT> TrueParamTestT;
  typedef MMParamTest<dataT> MMParamTestT;
  
  typedef CircularMaskOp< VolumeT, VolumeT> DataCircularMaskOpT;
  
  typedef Eigen::Matrix<dataT, 6, 1> ParamT;

  template < typename DifferentiatorT >
  class AlgorithmDifferentiatorFactory {
  };

  template < typename VolumeT >
  class AlgorithmDifferentiatorFactory<
    CentralDifferencesDifferentiator<VolumeT> > {
    protected: 
    typedef CentralDifferencesDifferentiator<VolumeT> DifferentiatorT;
  
    public:
  
    static DifferentiatorT* newVolDiffer(
      VolumeT *maskedRefVol) {
      return new DifferentiatorT(maskedRefVol);
    }
  };
  
  template <
    typename InterpolatorT,
    typename DifferentiatorT,
    typename ParamAccumulatorT,
    typename ConvergenceTestT,
    typename GradientUpdateTestT = void>
  class WGNNGMinimizerFactory {
    public:
    typedef Weighted_Gauss_Newton_New_Grad<
      InterpolatorT,
      ParamAccumulatorT,
      DataCircularMaskOpT,
      ConvergenceTestT,
      GradientUpdateTestT > MinimizerT;
   
    typedef Algorithm <
      dataT,
      VolumeT,
      InterpolatorT,
      WGNNGMinimizerFactory<
        InterpolatorT,
        DifferentiatorT,
        ParamAccumulatorT,
        ConvergenceTestT,
        GradientUpdateTestT >
      > AlgorithmT;
   
    WGNNGMinimizerFactory(VolumeT* fourierMaskedRefVol) :
      cubeSize(fourierMaskedRefVol->cubeSize) {}
   
    MinimizerT* newMinimizer(
      InterpolatorT *interpolator,
      typename AlgorithmT::DataCircularMaskOpT *imageMaskOp) {
      return new MinimizerT(
        interpolator,
        cubeSize,
        imageMaskOp);
    }

    protected:
    size_t cubeSize;
  };


  template<
    typename InterpolatorT,
    typename DifferentiatorT,
    typename ParamAccumulatorT,
    typename ConvergenceTestT,
    typename GradientUpdateTestT = void>
  class AlgorithmWGNNG : public Algorithm <
    dataT,
    VolumeT,
    InterpolatorT,
    WGNNGMinimizerFactory<
      InterpolatorT,
      DifferentiatorT,
      ParamAccumulatorT,
      ConvergenceTestT,
      GradientUpdateTestT >
    > {

    protected:
    typedef WGNNGMinimizerFactory<
      InterpolatorT,
      DifferentiatorT,
      ParamAccumulatorT,
      ConvergenceTestT,
      GradientUpdateTestT > MinimizerFactoryT;

    typedef Algorithm <
      dataT,
      VolumeT,
      InterpolatorT,
      MinimizerFactoryT > Parent;

    virtual void registerNewVolumeInner( 
      VolumeT *newVol,
      ParamT *p, 
      size_t *elapsedSteps) {
      
      this->fourierMaskVolume(newVol, &(this->fourierMaskedNewVol));

      ParamT initialParam; 
      initialParam << 0, 0, 0, 0, 0, 0;
      
      DifferentiatorT newVolDiffer(&(this->fourierMaskedNewVol));
      VolumeT newVolDx(this->cubeSize);
      newVolDiffer.xDerivative(&newVolDx);
  
      VolumeT newVolDy(this->cubeSize);
      newVolDiffer.yDerivative(&newVolDy);
  
      VolumeT newVolDz(this->cubeSize);
      newVolDiffer.zDerivative(&newVolDz);
      
      this->minimizer->minimize(
        &(this->fourierMaskedNewVol),
        &newVolDz, &newVolDy, &newVolDx, 
        &initialParam, p,
        this->maxSteps, this->stepSizeScale, this->stepSizeLimit,
        convergenceTest, gradientUpdateTest,
        elapsedSteps, NULL);
    
    }

    public:

    AlgorithmWGNNG(
      VolumeT *refVol,
      ConvergenceTestT *convergenceTest,
      GradientUpdateTestT *gradientUpdateTest) :
      Parent(refVol),
      convergenceTest(convergenceTest),
      gradientUpdateTest(gradientUpdateTest){}

    protected:
      ConvergenceTestT *convergenceTest;
      GradientUpdateTestT *gradientUpdateTest;

  };
  
  template <
    typename InterpolatorT,
    typename DifferentiatorT,
    typename ParamAccumulatorT,
    typename ConvergenceTestT,
    typename GradientUpdateTestT = void>
  class WGNRGMinimizerFactory {
    public:
    typedef Weighted_Gauss_Newton_Ref_Grad<
      InterpolatorT,
      ParamAccumulatorT,
      DataCircularMaskOpT,
      ConvergenceTestT,
      GradientUpdateTestT > MinimizerT;
   
    typedef Algorithm <
      dataT,
      VolumeT,
      InterpolatorT,
      WGNRGMinimizerFactory<
        InterpolatorT,
        DifferentiatorT,
        ParamAccumulatorT,
        ConvergenceTestT,
        GradientUpdateTestT >
      > AlgorithmT;

    WGNRGMinimizerFactory(VolumeT *fourierMaskedRefVol) :
      cubeSize(fourierMaskedRefVol->cubeSize),
      fourierMaskedRefdz(cubeSize),
      fourierMaskedRefdy(cubeSize),
      fourierMaskedRefdx(cubeSize)
      { 
        DifferentiatorT *volDiffer =
          AlgorithmDifferentiatorFactory<DifferentiatorT>::newVolDiffer(
            fourierMaskedRefVol); 

        volDiffer->xDerivative(&fourierMaskedRefdx);
        volDiffer->yDerivative(&fourierMaskedRefdy);
        volDiffer->zDerivative(&fourierMaskedRefdz);

        delete volDiffer;
      }

    MinimizerT* newMinimizer(
      InterpolatorT *interpolator,
      typename AlgorithmT::DataCircularMaskOpT *imageMaskOp) {
      return new MinimizerT(
        interpolator,
        imageMaskOp,
        &fourierMaskedRefdz,
        &fourierMaskedRefdy,
        &fourierMaskedRefdx,
        NULL);
    }

    protected:
    size_t cubeSize;
    VolumeT fourierMaskedRefdz;
    VolumeT fourierMaskedRefdy;
    VolumeT fourierMaskedRefdx;
  };


  template<
    typename InterpolatorT,
    typename DifferentiatorT,
    typename ParamAccumulatorT,
    typename ConvergenceTestT,
    typename GradientUpdateTestT = void>
  class AlgorithmWGNRG : public Algorithm <
    dataT,
    VolumeT,
    InterpolatorT,
    WGNRGMinimizerFactory<
      InterpolatorT,
      DifferentiatorT,
      ParamAccumulatorT,
      ConvergenceTestT,
      GradientUpdateTestT >
    > {

    protected:
    typedef WGNRGMinimizerFactory<
      InterpolatorT,
      DifferentiatorT,
      ParamAccumulatorT,
      ConvergenceTestT,
      GradientUpdateTestT > MinimizerFactoryT;

    typedef Algorithm <
      dataT,
      VolumeT,
      InterpolatorT,
      MinimizerFactoryT > Parent;

    virtual void registerNewVolumeInner( 
      VolumeT *newVol,
      ParamT *p, 
      size_t *elapsedSteps) {
      
      this->fourierMaskVolume(newVol, &(this->fourierMaskedNewVol));
     
      ParamT initialParam; 
      initialParam << 0, 0, 0, 0, 0, 0;
      
      this->minimizer->minimize(
        &(this->fourierMaskedNewVol), &initialParam, p,
        this->maxSteps, this->stepSizeScale, this->stepSizeLimit,
        convergenceTest, gradientUpdateTest,
        elapsedSteps, NULL);
    
    }

    public:

    AlgorithmWGNRG(
      VolumeT *refVol,
      ConvergenceTestT *convergenceTest,
      GradientUpdateTestT *gradientUpdateTest) :
      Parent(refVol),
      convergenceTest(convergenceTest),
      gradientUpdateTest(gradientUpdateTest) {
    }

    protected:
      ConvergenceTestT *convergenceTest;
      GradientUpdateTestT *gradientUpdateTest;
  };

  template<typename InterpolatorT, typename DifferentiatorT>
  class Algorithm1 : public AlgorithmWGNRG <
    InterpolatorT,
    DifferentiatorT,
    SumParamAccumulatorT,
    MMParamTestT,
    TrueParamTestT
    > {

    protected:
    typedef AlgorithmWGNRG <
      InterpolatorT,
      DifferentiatorT,
      SumParamAccumulatorT,
      MMParamTestT,
      TrueParamTestT
    > Parent; 

    public:
    Algorithm1(VolumeT *refVol, MMParamTestT *convergenceTest) :
      Parent(refVol, convergenceTest, &gradientUpdateTest) {}

    protected:
    TrueParamTestT gradientUpdateTest;

  };
  
  template<typename InterpolatorT, typename DifferentiatorT>
  class Algorithm2 : public AlgorithmWGNRG <
    InterpolatorT,
    DifferentiatorT,
    SumParamAccumulatorT,
    MMParamTestT
    > {

    protected:
    typedef AlgorithmWGNRG <
      InterpolatorT,
      DifferentiatorT,
      SumParamAccumulatorT,
      MMParamTestT
    > Parent; 

    public:
    Algorithm2(VolumeT *refVol, MMParamTestT *convergenceTest) :
      Parent(refVol, convergenceTest, NULL) {}

  };


  template<typename InterpolatorT, typename DifferentiatorT>
  class Algorithm3 : public AlgorithmWGNNG <
    InterpolatorT,
    DifferentiatorT,
    SumParamAccumulatorT,
    MMParamTestT,
    TrueParamTestT
    > {

    protected:
    typedef AlgorithmWGNNG <
      InterpolatorT,
      DifferentiatorT,
      SumParamAccumulatorT,
      MMParamTestT,
      TrueParamTestT
    > Parent; 

    public:
    Algorithm3(VolumeT *refVol, MMParamTestT *convergenceTest) :
      Parent(refVol, convergenceTest, &gradientUpdateTest) {}

    protected:
    TrueParamTestT gradientUpdateTest;
  };
  
  template<typename InterpolatorT, typename DifferentiatorT>
  class Algorithm4 : public AlgorithmWGNNG <
    InterpolatorT,
    DifferentiatorT,
    SumParamAccumulatorT,
    MMParamTestT
    > {

    protected:
    typedef AlgorithmWGNNG <
      InterpolatorT,
      DifferentiatorT,
      SumParamAccumulatorT,
      MMParamTestT
    > Parent; 

    public:
    Algorithm4(VolumeT *refVol, MMParamTestT *convergenceTest) :
      Parent(refVol, convergenceTest, NULL) {}

  };

  template<typename InterpolatorT, typename DifferentiatorT>
  class Algorithm5 : public AlgorithmWGNRG <
    InterpolatorT,
    DifferentiatorT,
    ComposeTransformParamAccumulatorT,
    MMParamTestT,
    TrueParamTestT
    > {

    protected:
    typedef AlgorithmWGNRG <
      InterpolatorT,
      DifferentiatorT,
      ComposeTransformParamAccumulatorT,
      MMParamTestT,
      TrueParamTestT
    > Parent; 

    public:
    Algorithm5(VolumeT *refVol, MMParamTestT *convergenceTest) :
      Parent(refVol, convergenceTest, &gradientUpdateTest) {}

    protected:
    TrueParamTestT gradientUpdateTest;

  };

  template<typename InterpolatorT, typename DifferentiatorT>
  class Algorithm6 : public AlgorithmWGNRG <
    InterpolatorT,
    DifferentiatorT,
    ComposeTransformParamAccumulatorT,
    MMParamTestT
    > {

    protected:
    typedef AlgorithmWGNRG <
      InterpolatorT,
      DifferentiatorT,
      ComposeTransformParamAccumulatorT,
      MMParamTestT
    > Parent; 

    public:
    Algorithm6(VolumeT *refVol, MMParamTestT *convergenceTest) :
      Parent(refVol, convergenceTest, NULL) {}

  };


  template<typename InterpolatorT, typename DifferentiatorT>
  class Algorithm7 : public AlgorithmWGNNG <
    InterpolatorT,
    DifferentiatorT,
    ComposeTransformParamAccumulatorT,
    MMParamTestT,
    TrueParamTestT
    > {

    protected:
    typedef AlgorithmWGNNG <
      InterpolatorT,
      DifferentiatorT,
      ComposeTransformParamAccumulatorT,
      MMParamTestT,
      TrueParamTestT
    > Parent; 

    public:
    Algorithm7(VolumeT *refVol, MMParamTestT *convergenceTest) :
      Parent(refVol, convergenceTest, &gradientUpdateTest) {}

    protected:
    TrueParamTestT gradientUpdateTest;
  };
  
  template<typename InterpolatorT, typename DifferentiatorT>
  class Algorithm8 : public AlgorithmWGNNG <
    InterpolatorT,
    DifferentiatorT,
    ComposeTransformParamAccumulatorT,
    MMParamTestT
    > {

    protected:
    typedef AlgorithmWGNNG <
      InterpolatorT,
      DifferentiatorT,
      ComposeTransformParamAccumulatorT,
      MMParamTestT
    > Parent; 

    public:
    Algorithm8(VolumeT *refVol, MMParamTestT *convergenceTest) :
      Parent(refVol, convergenceTest, NULL) {}

  };
}

#endif
