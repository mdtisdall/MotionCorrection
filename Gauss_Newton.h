#ifndef Gauss_Newton_h
#define Gauss_Newton_h

#include <iostream>

#include <sys/time.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>


#include <cfloat>

#include <string>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>

#include "TestHelper.h"

template <
  typename _ResidualOpT,
  typename _ResidualGradientAndHessianT,
  typename _InterpolatorT,
  typename _ParamAccumulatorT,
  typename _ConvergenceTestT = void
  >
class Gauss_Newton {
  public:
    typedef _ResidualOpT ResidualOpT;
    typedef _ResidualGradientAndHessianT ResidualGradientAndHessianT;
    typedef _InterpolatorT InterpolatorT;
    typedef _ConvergenceTestT ConvergenceTestT;
    typedef _ParamAccumulatorT ParamAccumulatorT;
    typedef typename InterpolatorT::VolumeT VolumeT;
    typedef typename InterpolatorT::CoordT CoordT;
    typedef typename VolumeT::value_type T;
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    Gauss_Newton(
      const InterpolatorT *interpRef,
      ResidualOpT *residualOp,
      ResidualGradientAndHessianT *residualGradientAndHessian,
      const size_t cubeSize 
      ) :
      interpRef(interpRef),
      residualOp(residualOp),
      residualGradientAndHessian(residualGradientAndHessian),
      cubeSize(cubeSize),
      cubeCenter(cubeCenterFromCubeSize(cubeSize)),
      residualGradient(6, cubeSize * cubeSize * cubeSize),
      pointList(3, cubeSize * cubeSize * cubeSize),
      transformedPointList(3, cubeSize * cubeSize * cubeSize),
      residual(cubeSize * cubeSize * cubeSize, 1) {
      
      generatePointList(&pointList, cubeSize, cubeCenter);
    }
    

  protected:
    typedef typename ResidualOpT::ResidualT ResidualT;
    typedef typename ResidualOpT::NewVolVecT NewVolVecT;
    typedef typename ResidualOpT::PointListT PointListT;
    typedef typename ResidualGradientAndHessianT::ResidualGradientT ResidualGradientT;
    typedef typename ResidualGradientAndHessianT::ResidualHessianT ResidualHessianT;
    typedef typename ResidualGradientAndHessianT::ResidualHessianLDLT ResidualHessianLDLT;
    typedef Eigen::Matrix< T, 3, 1 > PointT;
    
    static CoordT cubeCenterFromCubeSize(const size_t cubeSize) {
        return ((CoordT) cubeSize)/(CoordT)2.0 - (CoordT)0.5;
    }
    


    static void generatePointList(
      PointListT *pointList, const size_t cubeSize, const CoordT cubeCenter) {

      size_t offset = 0;

      for(size_t z = 0; z < cubeSize; z++) {
        CoordT zCoord = ((CoordT) z) - cubeCenter; 
        
        for(size_t y = 0; y < cubeSize; y++) {
          CoordT yCoord = ((CoordT) y) - cubeCenter; 
          
          for(size_t x = 0; x < cubeSize; x++, offset++) {
            CoordT xCoord = ((CoordT) x) - cubeCenter;

            pointList->col(offset) = PointT(zCoord, yCoord, xCoord);
          }
        }
      }
    }


    static void transformPointListWithParam(
      const ParamT *param,
      const PointListT *pointList,
      PointListT *transformedPointList) {

      typedef Eigen::AngleAxis<T> RotationT;
      typedef Eigen::Translation<T, 3> TranslationT;

      const Eigen::Matrix<T, 3, 1> rotVec = param->tail(3);

      const T angle = rotVec.norm();

      Eigen::Matrix<T, 3, 1> rotAxis;
      if(0 == angle) {
        rotAxis << 1, 0, 0;
      }
      else {
        rotAxis = rotVec.normalized();
      }
      
      transformedPointList->noalias() =
        ( RotationT(-angle, rotAxis) * TranslationT(-param->head(3)) ) *
        (*pointList);
    }
   
    void computeResidual(
      const VolumeT *newVolume,
      const ParamT *curParam) {

      transformPointListWithParam(
        curParam, &pointList, &transformedPointList);
      
      NewVolVecT newVolVec(
        newVolume->buffer,
        this->cubeSize * this->cubeSize * this->cubeSize, 1);
      
      (*residualOp)(newVolume, &newVolVec, &transformedPointList,
        curParam, &residual);
    }

    void minimize(
      const VolumeT *newVolume,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      const ParamT *initialParam,
      ParamT *finalParam,
      const size_t maxSteps = 20,
      const T stepSizeScale = 0.25,
      const T stepSizeLimit = 10e-6,
      const ConvergenceTestT *convergenceTest = NULL, 
      size_t *elapsedSteps = NULL,
      size_t *elapsedSearchSteps = NULL,
      double *elapsedTime = NULL 
      ) {
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }

      //
      // establish point coordinates based on initialParam guess
      //
      

      ParamT curParam = *initialParam;
      ParamT prevParam = curParam;


      ParamT reducedResidual;
      size_t step = 0;
      size_t searchStep = 0;
    
      computeResidual(newVolume, &curParam);

      T prevResidualNorm = this->residual.norm();
        
//      std::cout << "prevResidualNorm " << prevResidualNorm << std::endl; 

      
      for(; step < maxSteps; step++) {
//        std::cout << "-------" << std::endl; 
//        std::cout << "step " << step << std::endl; 
//        std::cout << "curParam: " << std::endl <<
//          curParam.transpose() << std::endl; 
//        std::cout << "residualNorm: " << prevResidualNorm << std::endl; 
//        std::cout << "residual(12305): " << residual(16911) << std::endl;


        //
        // compute the direction of the next step
        //
        
//        std::cout << "residualGradient(12305): " << std::endl <<
//          residualGradient.col(12305).transpose() << std::endl;

//        if(step == 1) 
//        {
//          std::string filePath("Moving_Weighted_Gauss_Newton_Fixed_M_New_Grad_tests/residualGradient");
//          int outputFile = open(filePath.c_str(), O_WRONLY | O_CREAT, 0600);
//    
//          if(-1 == outputFile) {
//            std::cerr << "Could not open file: " << filePath << std::endl;
//            exit(1);
//          }
//
//          int bytesWritten = 
//              ::write(outputFile, residualGradient.data(),
//                sizeof(T) * 6 * this->cubeSize * this->cubeSize * this->cubeSize);
//    
//          close(outputFile);
//          
//          exit(0);
//        }


//        {
//          std::string filePath("Static_Weighted_Gauss_Newton_New_Grad_tests/residual_" + std::to_string(step));
//          int outputFile = open(filePath.c_str(), O_WRONLY | O_CREAT, 0600);
//    
//          if(-1 == outputFile) {
//            std::cerr << "Could not open file: " << filePath << std::endl;
//            exit(1);
//          }
//
//          int bytesWritten = 
//              ::write(outputFile, residual.data(),
//                sizeof(T) * this->cubeSize * this->cubeSize * this->cubeSize);
//    
//          close(outputFile);
//        }

        reducedResidual.noalias() = this->residualGradient * this->residual;
    
//        std::cout << "approxResidualHessian: " << std::endl <<
//          this->approxResidualHessian << std::endl;
//        
//        std::cout << "reducedResidual: " << std::endl <<
//          reducedResidual.transpose() << std::endl;

        // This equation solves the parameter update, but with signs negated
        ParamT negParamUpdateDirection;
        negParamUpdateDirection.noalias() =
          this->residualHessianLDL.solve(reducedResidual);
       
//        std::cout << "negParamUpdateDirection: " << std::endl <<
//          negParamUpdateDirection.transpose() << std::endl;

        //
        // perform backtracking line search in the direction of the next step 
        //
     
        T stepSize = 1.0;

        bool improved = false;

        ParamT paramUpdate;

        while(!improved && stepSize >= stepSizeLimit) {
          searchStep++;

          // negate this to get the correct sign for the update
          paramUpdate = -negParamUpdateDirection * stepSize;
          
          ParamT newParam =
            ParamAccumulatorT::accumulate(&curParam, &paramUpdate);

          computeResidual(newVolume, &newParam);

          T newResidualNorm = this->residual.norm();

          // If the residual has become smaller since the last step, take step 
          if(newResidualNorm <= prevResidualNorm) {
            prevResidualNorm = newResidualNorm;
            curParam = newParam;
            improved = true;
          }
          // otherwise, make the step smaller
          else {
            stepSize *= stepSizeScale; 
          }
        }

        // if we can't make an improvement by stepping in this direction
        // then we should take no step and stop the minimization
        if(!improved) {
          break;
        }

        //
        // now check to see if the steps have converged
        //
        
        if(
          TestHelper<ConvergenceTestT>::eval(convergenceTest, &paramUpdate)
          ) {
            step++; 
            break; 
        }

        residualGradientAndHessian->updateResidualGradientAndApproxHessian(
          &transformedPointList,
          &curParam,
          refdz, refdy, refdx,
          &residualGradient, &approxResidualHessian, 
          &residualHessianLDL);
      }

      if(NULL != elapsedSteps) {
        *elapsedSteps = step; 
      }
      
      if(NULL != elapsedSearchSteps) {
        *elapsedSearchSteps = searchStep; 
      }

      *finalParam = curParam;

    }

  protected:
    const InterpolatorT *interpRef;
    ResidualOpT *residualOp;
    ResidualGradientAndHessianT *residualGradientAndHessian;
    const size_t cubeSize;
    const CoordT cubeCenter;
    ResidualGradientT residualGradient;
    ResidualHessianT approxResidualHessian;
    ResidualHessianLDLT residualHessianLDL;
    PointListT pointList;
    PointListT transformedPointList;
    ResidualT residual;
};


#endif
