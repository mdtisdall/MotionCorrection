#include "catch.hpp"

#include <complex>

#include "FFTOp.h"
#include "FFTWBuffer.h"
#include "VolumeAtAddressable.h"
#include "SymmetricHalfVolumeAtAddressable.h"

#include "CircularMaskOp.h"

#include "WeightFunction.h"

#include "BinaryFile.h"

TEST_CASE("A circular mask can be instantiated") {
    typedef double dataT;
    typedef std::complex<dataT> complexT;
    typedef FFTOp<dataT> DataFFTOpT;
    typedef DataFFTOpT::spatialVolumeT DataVolumeT; 
    typedef DataFFTOpT::fourierVolumeT ComplexVolumeT;
    typedef WeightFunction<dataT> WeightFunctionT;
    typedef CircularMaskOp<
      DataVolumeT,
      DataVolumeT,
      WeightFunctionT >
      DataCircularMaskOpT;
    typedef CircularMaskOp< ComplexVolumeT,
      SymmetricHalfVolumeAtAddressable< FFTWBuffer<dataT> >,
      WeightFunctionT >
      ComplexDataCircularMaskOpT;

    const size_t cubeSize = 32;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    DataVolumeT initialData(cubeSize);

    REQUIRE(cubeVectorLength * sizeof(dataT)
            == BinaryFile<DataVolumeT>::read(&initialData,
                "CircularMaskOp_tests/volInput.dat"));


    SECTION("and masking in the Fourier domain gives the precomputed answers") {
      WeightFunctionT weightFunc(cubeSize);
      ComplexDataCircularMaskOpT fourierMaskOp(cubeSize, &weightFunc);
      ComplexVolumeT fourierData(cubeSize);
      ComplexVolumeT maskedFourierData(cubeSize);
      DataFFTOpT fftOp(cubeSize);

      fftOp.forward(&initialData, &fourierData);
      
      fourierMaskOp.applyMask(&fourierData, &maskedFourierData); 
      
      SECTION("in the Fourier domain") { 
        size_t pointsLength = fourierData.totalPoints;
        std::vector<complexT> solutionData(pointsLength); 
        
        REQUIRE(pointsLength * sizeof(complexT)
              == BinaryFile< std::vector<complexT> >::read(&solutionData,
                  "CircularMaskOp_tests/CircularMaskOpFourierDomainOutput.dat"));

        for(size_t i = 0; i < pointsLength; i++) {
          INFO("i: " << i);
          INFO("fourierData.at(i): " << fourierData.at(i));
          INFO("maskedFourierData.at(i): " << maskedFourierData.at(i));
          INFO("solutionData[i]: " << solutionData[i]);
          REQUIRE(maskedFourierData.at(i).real() == Approx(solutionData[i].real()));
          REQUIRE(maskedFourierData.at(i).imag() == Approx(solutionData[i].imag()));
        }
      }
       
    

      SECTION("in the image domain") { 
        DataVolumeT resultData(cubeSize);
     
        fftOp.backward(&maskedFourierData, &resultData);

        size_t pointsLength = resultData.totalPoints;
        std::vector<dataT> solutionData(pointsLength); 
        
        REQUIRE(pointsLength * sizeof(dataT)
              == BinaryFile< std::vector<dataT> >::read(&solutionData,
                  "CircularMaskOp_tests/CircularMaskOpFourierOutput.dat"));

        for(size_t i = 0; i < pointsLength; i++) {
          INFO("i: " << i);
          REQUIRE(resultData.at(i) == Approx(solutionData[i]).epsilon(0.0005));
        }
      }
    }
    
    SECTION("and masking in the image domain gives the precomputed answers") {
      WeightFunctionT weightFunc(cubeSize);
      DataCircularMaskOpT imageMaskOp(cubeSize, &weightFunc);
      
      imageMaskOp.applyMask(&initialData); 
    
      size_t pointsLength = initialData.totalPoints;
      std::vector<dataT> solutionData(pointsLength);
    
      REQUIRE(pointsLength * sizeof(dataT)
            == BinaryFile< std::vector<dataT> >::read(&solutionData,
                "CircularMaskOp_tests/CircularMaskOpImageOutput.dat"));

      for(size_t i = 0; i < pointsLength; i++) {
        INFO("i: " << i);
        REQUIRE(initialData.at(i) == Approx(solutionData[i]));
      }
    }
}
