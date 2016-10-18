#include "catch.hpp"

#include "BinaryFile.h"
#include "Interpolation_3D.h"

#include <vector>
#include <complex>

TEST_CASE("a vector of complex floats can be written") {
    typedef std::complex<float> dataT;
    typedef std::vector< dataT > vectorT;

    vectorT initialData(100);
    
    for(size_t i = 0; i < 100; i++) {
        initialData[i] = dataT(i, -i); 
    }

    REQUIRE(100 * sizeof(dataT)
        == BinaryFile<vectorT>::write(&initialData, "temp.dat"));

    SECTION("and read back into a new vector") {
        vectorT readData(100);

        REQUIRE(100 * sizeof(dataT)
            == BinaryFile<vectorT>::read(&readData, "temp.dat"));

        SECTION("and the data is identical") {
            for(size_t i = 0; i < 100; i++) {
                REQUIRE(initialData[i] == readData[i]); 
            }
        }
    }
}
