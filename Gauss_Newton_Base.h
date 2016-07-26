#ifndef Gauss_Newton_Base_h
#define Gauss_Newton_Base_h

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
  typename _InterpolatorT,
  typename _ParamAccumulatorT,
  typename _ConvergenceTestT = void,
  typename _GradientUpdateTestT = void
  >
class Gauss_Newton_Base{
  public:
    typedef _InterpolatorT InterpolatorT;
    typedef _ConvergenceTestT ConvergenceTestT;
    typedef _GradientUpdateTestT GradientUpdateTestT;
    typedef _ParamAccumulatorT ParamAccumulatorT;
    typedef typename InterpolatorT::VolumeT VolumeT;
    typedef typename InterpolatorT::CoordT CoordT;
    typedef typename VolumeT::value_type T;
    typedef typename Eigen::Matrix<T, 6, 1> ParamT;

    Gauss_Newton_Base(
      const InterpolatorT *interpRef, 
      const size_t cubeSize 
      ) :
      interpRef(interpRef),
      cubeSize(cubeSize),
      cubeCenter(cubeCenterFromCubeSize(cubeSize)),
      residualGradient(6, cubeSize * cubeSize * cubeSize),
      pointList(3, cubeSize * cubeSize * cubeSize),
      transformedPointList(3, cubeSize * cubeSize * cubeSize),
      interpPoints(cubeSize * cubeSize * cubeSize, 1),
      residual(cubeSize * cubeSize * cubeSize, 1) {
      
      generatePointList(&pointList, cubeSize, cubeCenter);
    }
    

  protected:
    typedef Eigen::Matrix< T, 6, Eigen::Dynamic > ResidualGradientT;
    typedef Eigen::Matrix< T, 6, 6 > ResidualHessianT;
    typedef Eigen::LDLT< ResidualHessianT, Eigen::Upper > ResidualHessianLDLT;
    typedef Eigen::Matrix< T, 3, Eigen::Dynamic > PointListT;
    typedef Eigen::Matrix< T, 3, 1 > PointT;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> ResidualT;
    typedef Eigen::Map< Eigen::Matrix<T, Eigen::Dynamic, 1 > > NewVolVecT;
    
    static CoordT cubeCenterFromCubeSize(const size_t cubeSize) {
        return ((CoordT) cubeSize)/(CoordT)2.0 - (CoordT)0.5;
    }
    
    static void generateResidualGradientAndApproxHessian(
      ResidualGradientT *residualGradient,
      ResidualHessianT *approxResidualHessian,
      ResidualHessianLDLT *residualHessianLDL,
      const PointListT *pointList,
      const VolumeT *refdz,
      const VolumeT *refdy,
      const VolumeT *refdx,
      const size_t cubeSize,
      const CoordT cubeCenter,
      double *elapsedTime = NULL 
      ) {
 
      struct timeval timeBefore, timeAfter;

      if(NULL != elapsedTime) {
        gettimeofday(&timeBefore, NULL);
      }
    
      typedef Eigen::Map< Eigen::Matrix<T, 1, Eigen::Dynamic > > RefDMatT;

      RefDMatT refDzMat(refdz->buffer, refdz->totalPoints);      
      RefDMatT refDyMat(refdy->buffer, refdy->totalPoints);      
      RefDMatT refDxMat(refdx->buffer, refdx->totalPoints);      

      // the first three colums of the M matrix are just negative copies
      // of the spatial gradients
      residualGradient->row(0).noalias() = -refDzMat;
      residualGradient->row(1).noalias() = -refDyMat;
      residualGradient->row(2).noalias() = -refDxMat;

      // the last three colums of the M matrix are element-wise products 
      // of the point-lists and the spatial gradients  
      residualGradient->row(3).array() =
        pointList->row(2).array() * refDyMat.array();
      
      residualGradient->row(3).array() -=
        pointList->row(1).array() * refDxMat.array();
     
      residualGradient->row(4).array() =
        pointList->row(0).array() * refDxMat.array();
      
      residualGradient->row(4).array() -=
        pointList->row(2).array() * refDzMat.array();
      
      residualGradient->row(5).array() =
        pointList->row(1).array() * refDzMat.array();
     
      residualGradient->row(5).array() -=
        pointList->row(0).array() * refDyMat.array();

      //std::cout << "residualGradient->col(0): " <<
      //  residualGradient->col(0).transpose() << std::endl;

      // now we can compute the Hessian
      approxResidualHessian->setZero(6, 6);
     
      approxResidualHessian->template selfadjointView<Eigen::Upper>().rankUpdate(
        *(residualGradient));

      residualHessianLDL->compute(*approxResidualHessian);
      
      if(NULL != elapsedTime) { 
        gettimeofday(&timeAfter, NULL);
  
        *elapsedTime =
          ((double) (timeAfter.tv_sec - timeBefore.tv_sec)) * 1000.0 +
          ((double) (timeAfter.tv_usec - timeBefore.tv_usec)) * 0.001;
      }

//      std::cout << "approxResidualHessian: " << *approxResidualHessian << std::endl;
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
    
    virtual void computeResidual(
      const VolumeT *newVol,
      const NewVolVecT *newVolVec,
      const ParamT *param) {

      computeResidual(newVol, newVolVec, &pointList, param);
    }

    virtual void computeResidual(
      const VolumeT *newVol,
      const NewVolVecT *newVolVec,
      const PointListT *initialPoints,
      const ParamT *param) {

      transformPointListWithParam(
        param, initialPoints, &transformedPointList);

      const size_t pointListLength = newVol->totalPoints;

      PointT cubeCenterPoint(cubeCenter, cubeCenter, cubeCenter);

      for(size_t offset = 0; offset < pointListLength; offset++) {
        PointT curPoint =
          transformedPointList.col(offset) + cubeCenterPoint;

        interpPoints(offset, 0) =
          interpRef->interp(
            newVol->wrapIndex(curPoint(0)),
            newVol->wrapIndex(curPoint(1)),
            newVol->wrapIndex(curPoint(2))
          ); 
      }
 
//      std::cout << "interpPoints[0]: " << interpPoints(0) << std::endl;
//      std::cout << "newVolVec[0]: " << (*newVolVec)(0) << std::endl;

      residual.noalias() = interpPoints - (*newVolVec);
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
      const GradientUpdateTestT *gradientUpdateTest = NULL, 
      size_t *elapsedSteps = NULL,
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

      NewVolVecT newVolVec(
        newVolume->buffer,
        this->cubeSize * this->cubeSize * this->cubeSize, 1);

      ParamT reducedResidual;
      size_t step = 0;
      
      this->computeResidual(newVolume, &newVolVec,
        &pointList, &curParam);
      
      T prevResidualNorm = this->residual.norm();
        
//      std::cout << "prevResidualNorm " << prevResidualNorm << std::endl; 

      
      for(; step < maxSteps; step++) {
        std::cout << "-------" << std::endl; 
        std::cout << "step " << step << std::endl; 
        std::cout << "curParam: " << std::endl <<
          curParam.transpose() << std::endl; 
        std::cout << "residualNorm: " << prevResidualNorm << std::endl; 

        //
        // compute the direction of the next step
        //
        
        reducedResidual.noalias() = this->residualGradient * this->residual;
    
        //std::cout << "reducedResidual: " << std::endl <<
        //  reducedResidual << std::endl;

        // This equation solves the parameter update, but with signs negated
        ParamT negParamUpdateDirection;
        negParamUpdateDirection.noalias() =
          this->residualHessianLDL.solve(reducedResidual);
       
        //std::cout << "negParamUpdateDirection: " << std::endl <<
        //  negParamUpdateDirection << std::endl;

        //
        // perform backtracking line search in the direction of the next step 
        //
     
        T stepSize = 1.0;

        bool improved = false;

        ParamT paramUpdate;

        while(!improved && stepSize >= stepSizeLimit) {
          // negate this to get the correct sign for the update
          paramUpdate = -negParamUpdateDirection * stepSize;
          
          ParamT newParam =
            ParamAccumulatorT::accumulate(&curParam, &paramUpdate);
      
          this->computeResidual(newVolume, &newVolVec,
            &pointList, &newParam);

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

        //
        // now check to see if we need to update our gradients
        //

        if(
          TestHelper<GradientUpdateTestT>::eval(
            gradientUpdateTest, &paramUpdate) ) {
          
          PointListT accumulatedPoints;

          transformPointListWithParam(&curParam,
            &pointList, &accumulatedPoints);
    
          generateResidualGradientAndApproxHessian(
            &residualGradient, &approxResidualHessian, 
            &residualHessianLDL, &accumulatedPoints,
            refdz, refdy, refdx, cubeSize, cubeCenter);

          //std::cout << "approxResidualHessian:" << std::endl <<
          //  approxResidualHessian << std::endl;
        }
      }

      if(NULL != elapsedSteps) {
        *elapsedSteps = step; 
      }

      *finalParam = curParam;

    }

  protected:
    const InterpolatorT *interpRef;
    const size_t cubeSize;
    const CoordT cubeCenter;
    ResidualGradientT residualGradient;
    ResidualHessianT approxResidualHessian;
    ResidualHessianLDLT residualHessianLDL;
    PointListT pointList;
    PointListT transformedPointList;
    ResidualT interpPoints;
    ResidualT residual;
};


#endif
