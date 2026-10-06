/* M10c optional bridge owner. Bridge calls live in gemsEquilibrium.C so the
 * existing build toggle and header checksum recompile every bridge call. */
#ifndef gemsGasSpeciation_H
#define gemsGasSpeciation_H
#include "fvMesh.H"
#include "dictionary.H"
#include "gasSpeciation.H"
#include <memory>
namespace LESTO {
class gemsGasSpeciation {
  struct implementation;
  std::unique_ptr<implementation> impl_;
public:
  gemsGasSpeciation(const Foam::fvMesh&, const Foam::dictionary&,
                   const Foam::wordList&, const std::vector<gasSpecies>&,
                   const Foam::scalarList&);
  ~gemsGasSpeciation();
  bool local() const;
  Foam::string signature() const;
  bool constants(double T, double P, std::vector<double>& lnK) const;
  // A successful guarded IPM call supplies the kernel's initial potentials.
  // The kernel always enforces the physical (unfloored) element amounts.
  bool warm(Foam::label cell, double T, double P,
            const std::array<double,3>& b, double trace,
            std::array<double,3>& pi, int& iterations, bool& guarded);
};
}
#endif
