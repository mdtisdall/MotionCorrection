#include "catch.hpp"

#include "TrilinearInterpolator.h"

#include "VolumeAtAddressable.h"

#include "Interpolator3D_tests.h"

#include <vector>
#include <complex>

TEST_CASE("a trilinear interpolator can be created from a volume") {
    typedef std::complex<float> dataT;
    typedef VolumeAtAddressable< std::vector<dataT> > VolumeT; 
    typedef TrilinearInterpolator<VolumeT, float> InterpolatorT;

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = i; 
    }

    VolumeT volume(cubeSize, initialData);

    InterpolatorT interpolator(&volume);

    InterpolatorTests<InterpolatorT>::identity_tests(&interpolator, &volume); 

    SECTION("and all the points are averaged when interpolated in x") { 
        for(size_t z = 0; z < cubeSize; z++) {
            for(size_t y = 0; y < cubeSize; y++) {
                for(size_t x = 0; x < cubeSize; x++) {
                    dataT avg = (
                      volume.at(z, y, x) +
                      volume.at(z, y, volume.wrapIndex(x + 1))
                    ) * (float) 0.5;
                    REQUIRE(interpolator.interp(z, y, x + 0.5) == avg);
                }
            }
        }
    }
    
    SECTION("and all the points are averaged when interpolated in y") { 
        for(size_t z = 0; z < cubeSize; z++) {
            for(size_t y = 0; y < cubeSize; y++) {
                for(size_t x = 0; x < cubeSize; x++) {
                    dataT avg = (
                      volume.at(z, y, x) +
                      volume.at(z, volume.wrapIndex(y + 1), x)
                    ) * (float) 0.5;
                    REQUIRE(interpolator.interp(z, y + 0.5, x) == avg);
                }
            }
        }
    }
    
    SECTION("and all the points are averaged when interpolated in z") { 
        for(size_t z = 0; z < cubeSize; z++) {
            for(size_t y = 0; y < cubeSize; y++) {
                for(size_t x = 0; x < cubeSize; x++) {
                    dataT avg = (
                      volume.at(z, y, x) +
                      volume.at(volume.wrapIndex(z + 1), y, x)
                    ) * (float) 0.5;
                    REQUIRE(interpolator.interp(z + 0.5, y, x) == avg);
                }
            }
        }
    }
}

TEST_CASE(
  "a trilinear interpolator can be created from a constant-valued volume") {
    typedef std::complex<float> dataT;
    typedef VolumeAtAddressable< std::vector<dataT> > VolumeT; 
    typedef TrilinearInterpolator<VolumeT, float> InterpolatorT;

    const size_t cubeSize = 10;
    const size_t cubeVectorLength = cubeSize * cubeSize * cubeSize;

    std::vector<dataT> initialData(cubeVectorLength);

    const dataT constValue = 1.0;

    for(size_t i = 0; i < cubeVectorLength; i++) {
        initialData[i] = constValue; 
    }

    VolumeT volume(cubeSize, initialData);

    InterpolatorT interpolator(&volume);

    InterpolatorTests<InterpolatorT>::identity_tests(&interpolator, &volume); 
    InterpolatorTests<InterpolatorT>::constant_tests(
      &interpolator, &volume, constValue); 
}
