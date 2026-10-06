#ifndef LESTO_ChapmanEnskogDiffusivity_H
#define LESTO_ChapmanEnskogDiffusivity_H
#include <cmath>

namespace LESTO {
// Binary dilute-gas diffusion [m2/s]; T [K], p [Pa], masses [kg/mol],
// Lennard-Jones sigma [angstrom], epsilon/k [K]. Lorentz-Berthelot mixing.
// Neufeld, Janzen & Aziz, J. Chem. Phys. 57, 1100 (1972),
// doi:10.1063/1.1678363. See also doi:10.5194/acp-20-3669-2020, Eq. 4-6.
inline double chapmanEnskogDiffusivity (
  double T, double p, double M, double carrierM,
  double sigma, double carrierSigma, double epsilon, double carrierEpsilon
) {
  const double s = (sigma + carrierSigma)/2;
  const double t = T/std::sqrt(epsilon*carrierEpsilon);
  const double omega = 1.06036/std::pow(t, 0.15610)
    + 0.19300*std::exp(-0.47635*t) + 1.03587*std::exp(-1.52996*t)
    + 1.76474*std::exp(-3.89411*t);
  return 1.858e-7*std::pow(T, 1.5)
    *std::sqrt(1/(1000*M) + 1/(1000*carrierM))
    /((p/101325)*s*s*omega);
}
}
#endif
