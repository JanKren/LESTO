/*------------------------------------------------------------------------------

phaseChangeLedger.C

PURPOSE

Implementation of the cumulative phase-change mass ledger; see
phaseChangeLedger.H for the entries, the closure and the file formats.

------------------------------------------------------------------------------*/

#include "phaseChangeLedger.H"
#include "compensatedSum.H"
#include "emptyFvPatch.H"
#include "Time.H"
#include "OSspecific.H"

#include <cstdlib>
#include <sstream>

using namespace Foam;


Foam::word LESTO::phaseChangeLedger::exactDecimal(const scalar x) {
  std::string text;
  for (int digits = 15; digits <= 17; ++digits) {
    std::ostringstream os;
    os.precision(digits);
    os << x;
    text = os.str();
    if (std::strtod(text.c_str(), nullptr) == x) {
      break;
    }
  }
  return word(text, false);
}


LESTO::phaseChangeLedger::phaseChangeLedger (
  const fvMesh&   mesh,
  const wordList& pairNames,
  const wordList& condensateNames
)
:
  mesh_(mesh),
  pairNames_(pairNames),
  patchIds_(),
  nEntries_(0),
  values_ (
    IOobject("phaseChangeBalance",
             mesh.time().timeName(),
             "uniform",
             mesh.time(),
             IOobject::READ_IF_PRESENT,
             IOobject::AUTO_WRITE)
  ),
  layout_ (
    IOobject("phaseChangeLayout",
             mesh.time().timeName(),
             "uniform",
             mesh.time(),
             IOobject::READ_IF_PRESENT,
             IOobject::AUTO_WRITE)
  ),
  restarted_(false),
  files_()
{
  /*--------------------------------------------------------------------------
  Transport entries: the non-coupled, non-empty patches among the
  non-processor patches, which are the same on every rank (B1).
  --------------------------------------------------------------------------*/
  const fvBoundaryMesh& boundary = mesh.boundary();
  DynamicList<label> ids;
  for (label patchi = 0; patchi < boundary.size(); ++patchi) {
    if (patchi >= mesh.boundaryMesh().nNonProcessor()) {
      break;
    }
    if (!boundary[patchi].coupled() && !isA<emptyFvPatch>(boundary[patchi])) {
      ids.append(patchi);
    }
  }
  patchIds_.transfer(ids);

  /*--------------------------------------------------------------------------
  The entries, followed by their compensations.
  --------------------------------------------------------------------------*/
  nEntries_ =
    nGlobalEntries + (nBaseEntries + patchIds_.size())*pairNames_.size();
  const label size = 2*nEntries_;

  restarted_ = returnReduceOr(!values_.empty());

  /*--------------------------------------------------------------------------
  Layout: the names behind the positions of the entries.  A restart must
  have been written for the same pairs, condensates and transport patches,
  in the same order; the ledger would otherwise book its entries to the
  wrong pair or patch.
  --------------------------------------------------------------------------*/
  wordList patchNames(patchIds_.size());
  forAll(patchIds_, k) {
    patchNames[k] = boundary[patchIds_[k]].name();
  }

  if (restarted_) {

    const fileName layoutPath(layout_.objectRelPath());

    if (!returnReduceAnd(layout_.found("pairs"))) {
      FatalErrorInFunction
        << values_.objectRelPath() << " exists, but not " << layoutPath
        << ", which names its pairs and transport patches." << nl
        << "Copy the complete uniform/ directory of the run that wrote the "
        << "ledger." << exit(FatalError);
    }

    const wordList oldPairs(layout_.get<wordList>("pairs"));
    const wordList oldCondensates(layout_.get<wordList>("condensates"));
    const wordList oldPatches(layout_.get<wordList>("transportPatches"));

    if (
         oldPairs != pairNames_
      || oldCondensates != condensateNames
      || oldPatches != patchNames
    ) {
      FatalErrorInFunction
        << "The ledger " << values_.objectRelPath() << " was written for"
        << nl << "    pairs " << oldPairs << ", condensates "
        << oldCondensates << ", transport patches " << oldPatches << nl
        << "but this run has" << nl
        << "    pairs " << pairNames_ << ", condensates " << condensateNames
        << ", transport patches " << patchNames << " (" << layoutPath << ")."
        << nl << "Its entries are addressed by position; restart with the "
        << "pairs of constant/thermochemistryProperties and the patches in "
        << "the order of that run." << exit(FatalError);
    }
  }

  layout_.set("pairs", pairNames_);
  layout_.set("condensates", condensateNames);
  layout_.set("transportPatches", patchNames);

  if (!restarted_) {
    values_.resize(size, 0.0);
  } else if (returnReduceOr(values_.size() != size)) {
    FatalErrorInFunction
      << values_.objectPath() << " holds " << values_.size()
      << " values, but " << pairNames_.size() << " pair(s) with "
      << patchIds_.size() << " non-coupled patch(es) need " << size
      << " = 2 x (" << label(nGlobalEntries) << " + ("
      << label(nBaseEntries) << " + " << patchIds_.size() << ") per pair)"
      << ", the entries and their compensations.  Were the pairs or the"
      << " patches changed since the ledger was written, or was it written"
      << " by another version of the solver?" << exit(FatalError);
  }
}


void LESTO::phaseChangeLedger::addTo(const label i, const scalar x) {
  compensatedAdd(values_[i], values_[nEntries_ + i], x);
}


Foam::scalar LESTO::phaseChangeLedger::totalTransport (
  const label pairi
) const {
  compensatedSum sum;
  forAll(patchIds_, k) {
    const label i = index(pairi, nBaseEntries + k);
    sum.add(values_[i]);
    sum.add(values_[nEntries_ + i]);
  }
  return sum.value();
}


Foam::scalar LESTO::phaseChangeLedger::supplied(const label pairi) const {

  /* sign*(value + compensation) of every term, added compensated */
  compensatedSum sum;
  auto term = [&](const label i, const scalar sign) {
    sum.add(sign*values_[i]);
    sum.add(sign*values_[nEntries_ + i]);
  };

  term(index(pairi, INITIAL), 1);
  term(index(pairi, RELEASED), 1);
  forAll(patchIds_, k) {
    term(index(pairi, nBaseEntries + k), -1);
  }
  term(index(pairi, REMOVED), -1);
  term(index(pairi, CLAMPED), -1);
  term(index(pairi, SOLVER_DEFECT), -1);
  term(index(pairi, RESTART), 1);
  if (!reaction_.empty()) sum.add(reaction_[pairi]);

  return sum.value();
}


Foam::scalar LESTO::phaseChangeLedger::inflow(const label pairi) const {
  scalar sum = 0;
  forAll(patchIds_, k) {
    sum += max(-transport(pairi, k), scalar(0));
  }
  return sum;
}


Foam::scalar LESTO::phaseChangeLedger::reference(const label pairi) const {
  const phaseChangeLedger& l = *this;
  return max (
    max(mag(l(pairi, INITIAL) + l(pairi, RELEASED)), mag(l(pairi, HELD))),
    inflow(pairi)
  );
}


void LESTO::phaseChangeLedger::openFiles() {

  if (!Pstream::master()) {
    return;
  }

  const Time& runTime = mesh_.time();
  const fileName dir (
    runTime.globalPath()/"postProcessing"/"phaseChangeBalance"
   /runTime.timeName()
  );
  mkDir(dir);

  files_.resize(pairNames_.size());

  forAll(pairNames_, pairi) {
    files_.set(pairi, new OFstream(dir/(pairNames_[pairi] + ".dat")));
    OFstream& os = files_[pairi];
    os.precision(17);

    os << "# phase-change mass balance of " << pairNames_[pairi]
       << " [kg]; closure = gas + wall + sample - (initial + released"
       << " - transport - removed - clamped - solverDefect + restart)" << nl
       << "# the time column is the exact time at the end of the step [s],"
       << " as many digits as needed to read it back exactly"
       << nl
       << "# time gas wall sample initial released";
    forAll(patchIds_, k) {
      os << " transport_" << mesh_.boundary()[patchIds_[k]].name();
    }
    os << " removed clamped solverDefect restart projected"
       << " evaporatedAboveWarning closure" << endl;
  }
}


void LESTO::phaseChangeLedger::writeLine (
  const label  pairi,
  const scalar gas,
  const scalar wall,
  const scalar sample,
  const scalar closure
) {

  if (!Pstream::master() || files_.empty()) {
    return;
  }

  const phaseChangeLedger& l = *this;
  OFstream& os = files_[pairi];

  os << exactDecimal(time()) << ' ' << gas << ' ' << wall << ' '
     << sample << ' ' << l(pairi, INITIAL) << ' ' << l(pairi, RELEASED);
  forAll(patchIds_, k) {
    os << ' ' << transport(pairi, k);
  }
  os << ' ' << l(pairi, REMOVED) << ' ' << l(pairi, CLAMPED) << ' '
     << l(pairi, SOLVER_DEFECT) << ' ' << l(pairi, RESTART) << ' '
     << l(pairi, PROJECTED) << ' ' << l(pairi, EVAPORATED_ABOVE_WARNING)
     << ' ' << closure << endl;
}
