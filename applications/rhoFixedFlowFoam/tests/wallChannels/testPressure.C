#include "../../phaseChangeKinetics.H"
#include <cfenv>
#include <iostream>

int main() {
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  for (double product : {0., std::numeric_limits<double>::denorm_min(), 1e-319, 1e-100, 1., 101325.}) {
    const long double logBar = (LESTO::logBarPressure(product)-std::log(1e30L))/3;
    const double actual = LESTO::pressureFromLogBar(logBar);
    // Independent cube-root reference evaluated in extended precision.
    const long double expected = 1e5L*std::cbrt(static_cast<long double>(product)/1e35L);
    if (product == 0 ? actual != 0 : std::fabs(actual/expected-1)>1e-14L) return 1;
  }
  if (LESTO::pressureFromLogBar(-1e4L)!=0) return 1;
  try {LESTO::pressureFromLogBar(1e4L);return 1;} catch(const std::overflow_error&) {}
  try {LESTO::logBarPressure(-1);return 1;} catch(const std::domain_error&) {}
  std::cout << "Reactive wall pressure: zero/subnormal/ordinary inputs, independent cube roots <=1e-14, overflow refused; FPE traps enabled\n";
  return 0;
}
