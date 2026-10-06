// SPDX-License-Identifier: GPL-3.0-or-later
#include "../../gasSpeciation.H"
#include <cfenv>
#include <iostream>

int main() {
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  const std::vector<LESTO::gasSpecies> formulas = {
    {1,0,0},{1,0,1},{1,0,2},{0,1,0},{0,2,0},
    {0,1,1},{0,1,3},{0,0,1},{0,0,2}
  };
  LESTO::gasSpeciation kernel(formulas);
  const double k[] = {0,30,60,0,10,35,70,0,20};
  // At these dilute bulks, molecular concentrations are below double's
  // representable range. Independent monatomic balances are the reference.
  const double bulks[] = {1e-280,1e-300,3.0847282814072739e-310,1e-315,
                         1e-320,std::numeric_limits<double>::denorm_min()};
  for (double bulk : bulks) for (double scale : {1.,20.,32.94,1000.}) {
    std::array<double,3> b = {bulk,bulk,bulk}, pi = {NAN,NAN,NAN};
    double c[9];
    if (kernel.solve(k,scale,b,pi,c)<0) return 1;
    for (int i=0;i<9;++i) {
      const bool mono = i==0 || i==3 || i==7;
      if (!std::isfinite(c[i]) || c[i] != (mono ? bulk : 0.)) {
        std::cerr << "Tiny bulk mismatch: " << bulk << ", scale " << scale
                  << ", species " << i << ", result " << c[i] << '\n';
        return 1;
      }
    }
  }
  std::cout << "24 tiny-bulk equilibria: independent monatomic balances exact, FPE traps enabled\n";
}
