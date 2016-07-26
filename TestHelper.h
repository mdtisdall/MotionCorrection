#ifndef TestHelper_h
#define TestHelper_h

template <typename ParamTestT>
class TestHelper {
  public:
    static bool eval(
      const ParamTestT *test, const typename ParamTestT::ParamT *param) {
      return (NULL != test && (*test)(param));
    }
};

template <>
class TestHelper<void> {
  public:
    static bool eval(const void *test, const void *param) {
      return false;
    }
};


#endif

