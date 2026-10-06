#include "../../ChapmanEnskogDiffusivity.H"
#include <cassert>
#include <iostream>
int main() {
  using LESTO::chapmanEnskogDiffusivity;
  const double d = chapmanEnskogDiffusivity(300, 101325, .0280134, .0319988,
    3.798, 3.467, 71.4, 106.7);
  // N2-O2 at 300 K, 1 atm: independent kinetic-theory SI expression
  // 3/(16 p sigma^2 Omega) sqrt(2(kT)^3/(pi mu)); rounded CE prefactor.
  const double k = 1.380649e-23, N = 6.02214076e23;
  const double mu = (.0280134*.0319988/(.0280134+.0319988))/N;
  const double t = 300/std::sqrt(71.4*106.7);
  const double omega = 1.06036/std::pow(t,.15610)+.193*std::exp(-.47635*t)
    +1.03587*std::exp(-1.52996*t)+1.76474*std::exp(-3.89411*t);
  const double sigma = (3.798+3.467)/2*1e-10;
  const double reference = 3/(16*101325*sigma*sigma*omega)
    *std::sqrt(2*std::pow(k*300,3)/(std::acos(-1.)*mu));
  assert(std::abs(d/reference-1) < .001);
  assert(d > 2e-5 && d < 2.2e-5);
  assert(chapmanEnskogDiffusivity(300, 202650, .0280134, .0319988,
    3.798, 3.467, 71.4, 106.7) == d/2);
  assert(chapmanEnskogDiffusivity(300, 101325, .0319988, .0280134,
    3.467, 3.798, 106.7, 71.4) == d);
  std::cout << "Chapman-Enskog SI reference, inverse pressure and symmetry: " << d << '\n';
}
