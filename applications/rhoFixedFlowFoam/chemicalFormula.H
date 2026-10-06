#ifndef LESTO_chemicalFormula_H
#define LESTO_chemicalFormula_H

#include "dictionary.H"
#include <map>
#include <cmath>

namespace LESTO {
using chemicalFormula = std::map<Foam::word, Foam::scalar>;

// Atomic masses [kg/mol]. Keep the convention explicit: formula validation
// must not silently accept a misspelled element or a different mass convention.
inline chemicalFormula readFormula(const Foam::dictionary& species, double M) {
  chemicalFormula result;
  if (!species.found("formula")) return result;
  const Foam::dictionary& d = species.subDict("formula");
  const chemicalFormula atomicMass {
    {"Pb", 0.2072}, {"Bi", 0.2089804}, {"I", 0.12690447},
    {"He", 0.004002602}, {"Ar", 0.039948}, {"N", 0.014007},
    {"O", 0.015999}, {"H", 0.001008}, {"K", 0.0390983}
  };
  double calculated = 0;
  for (const Foam::entry& e : d) {
    const auto mass = atomicMass.find(e.keyword());
    const double nu = d.get<Foam::scalar>(e.keyword());
    if (mass == atomicMass.end() || !std::isfinite(nu) || nu <= 0
      || nu != std::floor(nu)) {
      FatalIOErrorInFunction(d) << "formula requires known elements and "
        << "positive integer atom counts; got " << e.keyword() << " " << nu
        << Foam::exit(Foam::FatalIOError);
    }
    result.emplace(e.keyword(), nu);
    calculated += mass->second*nu;
  }
  if (result.empty() || !std::isfinite(M) || M <= 0
    || std::abs(M/calculated - 1) > 1e-6) {
    FatalIOErrorInFunction(species) << "formula molarMass is " << calculated
      << " kg/mol; specified molarMass " << M
      << " must agree within relative 1e-6." << Foam::exit(Foam::FatalIOError);
  }
  return result;
}
}
#endif
