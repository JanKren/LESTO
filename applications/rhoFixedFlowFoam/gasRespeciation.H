/* M10b: cell-local gas equilibrium and persistent species/element accounts.
 * Default engine none creates no objects or output. Formation pressures use
 * the 1 bar standard state; the frozen inert carrier remains unchanged. */
#ifndef gasRespeciationAdapter_H
#define gasRespeciationAdapter_H
#include "gasSpeciation.H"
#include "gemsGasSpeciation.H"
#include "interfaceExchange.H"
#include "binaryIOList.H"
#include "compensatedSum.H"

namespace LESTO {
class gasRespeciation {
  enum entry { INITIAL, RELEASED, TRANSPORT, EXCHANGED, DEFECT, RESTART,
               REACTION, COPRODUCT, HELD, REACTION_VARIATION, NENTRY };
  const Foam::fvMesh& mesh_;
  const Foam::fluidThermo& thermo_;
  const Foam::volScalarField& rho_;
  Foam::PtrList<Foam::volScalarField>& fields_;
  interfaceExchange& exchange_;
  bool active_ = false;
  Foam::label interval_ = 1;
  Foam::scalar minT_ = 350;
  Foam::labelList ids_, index_;
  Foam::scalarList masses_;
  std::vector<gasSpecies> formulas_;
  std::vector<std::vector<double>> constants_;
  Foam::autoPtr<gasSpeciation> kernel_;
  Foam::autoPtr<gemsGasSpeciation> gems_;
  Foam::autoPtr<binaryGlobalIOList<Foam::scalar>> values_, elements_;
  Foam::autoPtr<binaryGlobalIOList<Foam::scalar>> condensates_;
  Foam::scalarList flux_, defect_;
  Foam::PtrList<Foam::OFstream> files_, elementFiles_;
  Foam::label offset(Foam::label i, entry e) const { return i*NENTRY+e; }
  Foam::scalar value(Foam::label i, entry e) const;
  void add(Foam::label i, entry e, Foam::scalar x);
  void set(Foam::label i, entry e, Foam::scalar x);
  Foam::scalar inventory(Foam::label i) const;
public:
  gasRespeciation(const Foam::fvMesh&, const Foam::fluidThermo&,
                 const Foam::volScalarField&, Foam::PtrList<Foam::volScalarField>&,
                 const Foam::wordList&, const Foam::List<Foam::word>&,
                 const Foam::List<Foam::scalar>&, interfaceExchange&);
  bool active() const { return active_; }
  bool gasOnly(Foam::label si) const { return active_ && index_[si] >= 0
                                      && exchange_.pairOfGas(si) < 0; }
  void recordSolve(Foam::label si, const Foam::fvScalarMatrix&);
  void checkTransport(Foam::label si) const;
  void apply();
  void book();
  void writeElements();
  void addGasProfiles(const Foam::word&, Foam::scalarField&, Foam::scalar&) const;
};
}
#endif
