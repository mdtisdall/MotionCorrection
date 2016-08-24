#include "catch.hpp"

#include "DerivWeightFunction.h"
#include "LookUpTableFunction.h"

TEST_CASE("A derivative weight function lookup table can be instantiated") {
  typedef float dataT;
  typedef DerivWeightFunction<dataT> DerivWeightFunctionT;
  typedef FunctionLookupTable<DerivWeightFunctionT> LookupT;

  const size_t cubeSize = 32;
  const dataT radius = ((dataT) cubeSize) / (dataT) 2.0;


  SECTION(" and using the lookup table gives errors approximately the same as direct evaluation") {
    DerivWeightFunctionT func(cubeSize);
    LookupT lut(cubeSize, 10e6, &func);
  
    for(dataT i = 0; i <= radius; i += 0.001) {
      INFO("i: " << i);
      INFO("func(i): " << func(i));
      INFO("lut(i): " << lut(i));
      REQUIRE(func(i) == Approx(lut(i)).epsilon(0.00001));
    }
  }
}
