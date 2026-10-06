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
  // Synthetic formation constants, independent of external thermodynamic
  // records. The first cold solve formerly exhausted 200 iterations:
  // subtraction in the iodine derivative could make it negative.
  const double k[] = {0,161.13759369955756,10.003964640827878,0,
                     35.850403165015415,100.42985376482571,
                     323.55159582233881,0,80.965703031464329};
  int count = 0;
  for (double scale : {30.,20.,1000.})
    for (double pb : {2e-159,2e-15})
      for (double ratio : {1.,2.,3.})
        for (double offset : {0.,-1e-14,1e-14}) {
          const std::array<double,3> b = {pb,5e-7,ratio*5e-7*(1+offset)};
          std::array<double,3> pi = {NAN,NAN,NAN};
          double c[9], again[9];
          if (kernel.solve(k,scale,b,pi,c)<0) {
            std::cerr << "Strong-binding cold solve did not converge\n";
            return 1;
          }
          std::array<long double,3> totals = {0,0,0};
          for (int j=0;j<9;++j) {
            if (!std::isfinite(c[j]) || c[j]<0) return 1;
            totals[0] += formulas[j].nPb*static_cast<long double>(c[j]);
            totals[1] += formulas[j].nBi*static_cast<long double>(c[j]);
            totals[2] += formulas[j].nI*static_cast<long double>(c[j]);
          }
          for (int e=0;e<3;++e) if (std::fabs(totals[e]/b[e]-1)>1e-14L) {
            std::cerr << "Strong-binding element balance failed: " << e << '\n';
            return 1;
          }
          if (kernel.solve(k,scale,b,pi,again)<0) return 1;
          for (int j=0;j<9;++j) if (std::fabs(again[j]-c[j])>1e-15*b[1]) {
            std::cerr << "Strong-binding warm/cold equilibrium differs\n";
            return 1;
          }
          ++count;
        }
  std::cout << count << " strong-binding equilibria: convergence, all-element "
               "balance and warm/cold agreement, FPE traps enabled\n";
}
