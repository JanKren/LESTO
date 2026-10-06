/*------------------------------------------------------------------------------

axialProfiles.C

PURPOSE

Implementation of the axial profiles of rhoFixedFlowFoam; see axialProfiles.H
for the definitions, the dictionary and the file format, and profileBinning.H
for the overlap weights and the onset search.

Parallel hygiene: every reduction is called on all ranks (the constructor
and write() are collective), on lists sized by the dictionary only; files
are written by the master only.

------------------------------------------------------------------------------*/

#include "axialProfiles.H"
#include "profileBinning.H"
#include "surfaceFields.H"
#include "compensatedSum.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "StringStream.H"
#include "Time.H"

#include <algorithm>

using namespace Foam;


/*------------------------------------------------------------------------------
Local helpers
------------------------------------------------------------------------------*/

namespace {

  /* the number x as text, as the log prints it */
  string text(const scalar x) {
    OStringStream os;
    os << x;
    return os.str();
  }

  /* the file's status word of an onset search */
  const char* statusWord(const int status) {
    return status == LESTO::ONSET_FOUND ? "found"
      : status == LESTO::ONSET_NO_PEAK ? "noPeak"
      : status == LESTO::ONSET_INSIGNIFICANT ? "insignificant"
      : "notBracketed";
  }
}


/*------------------------------------------------------------------------------
Exact overlap weights of the intervals [lo_i, hi_i] with the bins
------------------------------------------------------------------------------*/

void LESTO::axialProfiles::binIntervals (
  const scalarField& edges,
  const scalarField& lo,
  const scalarField& hi,
  labelList&         start,
  labelList&         bin,
  scalarField&       weight,
  scalarField&       outside
) {

  const int nBins = int(edges.size() - 1);
  DynamicList<label>  bins(lo.size());
  DynamicList<scalar> weights(lo.size());

  start.resize(lo.size() + 1);
  outside.resize(lo.size());
  start[0] = 0;

  forAll(lo, i) {
    outside[i] = intervalOverlap (
      lo[i], hi[i], edges.cdata(), nBins,
      [&](const int b, const double w) {
        bins.append(label(b));
        weights.append(w);
      }
    );
    start[i + 1] = bins.size();
  }

  bin.transfer(bins);
  weight.transfer(weights);
}


/*------------------------------------------------------------------------------
Construction
------------------------------------------------------------------------------*/

LESTO::axialProfiles::axialProfiles (
  const fvMesh&         mesh,
  const dictionary&     dict,
  const labelList&      elementPatch,
  const labelList&      elementFace,
  const labelList&      elementCell,
  const labelList&      wallPatch,
  const labelList&      wallFace,
  const volScalarField& T,
  const boolList&       sampleCells,
  const scalar          significance,
  const scalar          tStart,
  const scalar          deltaT,
  const scalarList*     written,
  const bool            restarted
)
:
  mesh_(mesh),
  axis_(dict.getOrDefault<vector>("axis", vector(1, 0, 0))),
  origin_(dict.getOrDefault<point>("origin", Zero)),
  sets_(),
  sampleCell_(sampleCells),
  times_(),
  done_(),
  written_(),
  fraction_(0.01),
  significance_(significance),
  fromInlet_(true),
  fromPeak_(true),
  haveSearchFrom_(false),
  searchFrom_(0),
  wallPlusGas_(true)
{
  const scalar axisLength = mag(axis_);
  if (!(axisLength > 0)) {
    FatalIOErrorInFunction(dict) << "profiles: the axis " << axis_
      << " has no length." << exit(FatalIOError);
  }
  axis_ /= axisLength;

  /*--------------------------------------------------------------------------
  Axial extent [lo, hi] of every face from its points, of every cell from
  its faces, of every interface element from its wall face (layer faceCells)
  or from its cell (layer distance, elementPatch -1: the deposit of a layer
  cell is binned like the gas of that cell), and of every wall face of the
  frozen wall temperature and area per bin.
  --------------------------------------------------------------------------*/
  const pointField& points = mesh_.points();
  const faceList& faces = mesh_.faces();

  scalarField faceLo(faces.size()), faceHi(faces.size());
  forAll(faces, facei) {
    scalar lo = GREAT, hi = -GREAT;
    for (const label pointi : faces[facei]) {
      const scalar x = (points[pointi] - origin_) & axis_;
      lo = min(lo, x);
      hi = max(hi, x);
    }
    faceLo[facei] = lo;
    faceHi[facei] = hi;
  }

  scalarField cellLo(mesh_.nCells(), GREAT), cellHi(mesh_.nCells(), -GREAT);
  const labelUList& owner = mesh_.faceOwner();
  const labelUList& neighbour = mesh_.faceNeighbour();
  forAll(owner, facei) {
    cellLo[owner[facei]] = min(cellLo[owner[facei]], faceLo[facei]);
    cellHi[owner[facei]] = max(cellHi[owner[facei]], faceHi[facei]);
  }
  forAll(neighbour, facei) {
    cellLo[neighbour[facei]] = min(cellLo[neighbour[facei]], faceLo[facei]);
    cellHi[neighbour[facei]] = max(cellHi[neighbour[facei]], faceHi[facei]);
  }

  const scalar meshLo = gMin(cellLo);
  const scalar meshHi = gMax(cellHi);

  const label nElements = elementPatch.size();
  scalarField elementLo(nElements), elementHi(nElements);

  forAll(elementPatch, e) {
    const label patchi = elementPatch[e];
    if (patchi < 0) {
      elementLo[e] = cellLo[elementCell[e]];
      elementHi[e] = cellHi[elementCell[e]];
    } else {
      const label facei = mesh_.boundaryMesh()[patchi].start() + elementFace[e];
      elementLo[e] = faceLo[facei];
      elementHi[e] = faceHi[facei];
    }
  }

  const label nWallFaces = wallPatch.size();
  scalarField wallLo(nWallFaces), wallHi(nWallFaces);
  scalarField wallFaceArea(nWallFaces), wallFaceT(nWallFaces);
  const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();

  forAll(wallPatch, i) {
    const label patchi = wallPatch[i];
    const label facei = mesh_.boundaryMesh()[patchi].start() + wallFace[i];
    wallLo[i] = faceLo[facei];
    wallHi[i] = faceHi[facei];
    wallFaceArea[i] = magSfb[patchi][wallFace[i]];
    wallFaceT[i] = T.boundaryField()[patchi][wallFace[i]];
  }

  /*--------------------------------------------------------------------------
  Onset and observable
  --------------------------------------------------------------------------*/
  const dictionary onset(dict.subOrEmptyDict("onset"));
  fraction_ = onset.getOrDefault<scalar>("fraction", 0.01);
  if (!(fraction_ > 0 && fraction_ <= 1)) {
    FatalIOErrorInFunction(dict) << "profiles: onset fraction " << fraction_
      << " must lie in (0, 1]." << exit(FatalIOError);
  }

  wordList modes({"fromInlet", "fromPeak"});
  if (onset.found("search")) {
    ITstream& is = onset.lookup("search");
    if (is.size() == 1 && is.front().isWord()) {
      modes = wordList(1, is.front().wordToken());
    } else {
      modes = onset.get<wordList>("search");
    }
  }
  fromInlet_ = modes.found("fromInlet");
  fromPeak_ = modes.found("fromPeak");
  for (const word& mode : modes) {
    if (mode != "fromInlet" && mode != "fromPeak") {
      FatalIOErrorInFunction(dict) << "profiles: unknown onset search "
        << mode << ".  Choose fromInlet, fromPeak or both, e.g. "
        << "search (fromInlet fromPeak);" << exit(FatalIOError);
    }
  }
  if (!fromInlet_ && !fromPeak_) {
    FatalIOErrorInFunction(dict) << "profiles: no onset search given."
      << exit(FatalIOError);
  }
  haveSearchFrom_ = onset.readIfPresent("searchFrom", searchFrom_);
  significance_ = onset.getOrDefault<scalar>("significance", significance);
  if (!(significance_ >= 0)) {
    FatalIOErrorInFunction(dict) << "profiles: onset significance "
      << significance_ << " must be >= 0." << exit(FatalIOError);
  }

  const word observable(dict.getOrDefault<word>("observable", "wallPlusGas"));
  if (observable != "wall" && observable != "wallPlusGas") {
    FatalIOErrorInFunction(dict) << "profiles: unknown observable "
      << observable << ".  Choose wall or wallPlusGas." << exit(FatalIOError);
  }
  wallPlusGas_ = observable == "wallPlusGas";

  /*--------------------------------------------------------------------------
  Bin sets: edges, overlap weights, frozen wall area and temperature
  --------------------------------------------------------------------------*/
  const dictionary& binsDict = dict.subDict("bins");
  if (binsDict.empty()) {
    FatalIOErrorInFunction(binsDict) << "profiles: no bin set in bins."
      << exit(FatalIOError);
  }

  const scalarField& V = mesh_.V();
  const scalar meshVolume = gSum(V);
  const scalar interfaceArea = gSum(wallFaceArea);

  sets_.resize(binsDict.size());
  label seti = 0;
  OStringStream setsText;

  for (const entry& e : binsDict) {

    if (!e.isDict()) {
      FatalIOErrorInFunction(binsDict) << "profiles: bin set " << e.keyword()
        << " must be a dictionary." << exit(FatalIOError);
    }
    const dictionary& bd = e.dict();
    const scalar binsMin = bd.getOrDefault<scalar>("min", meshLo);
    const scalar binsMax = bd.getOrDefault<scalar>("max", meshHi);
    const bool haveNumber = bd.found("nBins");
    const bool haveWidth = bd.found("width");

    if (haveNumber == haveWidth || !(binsMax > binsMin)) {
      FatalIOErrorInFunction(bd) << "profiles: bin set " << e.keyword()
        << " needs max > min (defaults: the extent of the mesh along the "
        << "axis) and exactly one of nBins or width." << exit(FatalIOError);
    }

    /* at most maxBins bins (a width in the wrong unit would otherwise
       overflow the bin count); a width at or above the range gives the
       single bin [min, max] (profileBinning.H) */
    const scalar maxBins = 1e6;
    std::vector<double> edges;
    if (haveNumber) {
      const label nBins = bd.get<label>("nBins");
      if (nBins < 1 || nBins > maxBins) {
        FatalIOErrorInFunction(bd) << "profiles: bin set " << e.keyword()
          << " needs 1 <= nBins <= " << maxBins << "." << exit(FatalIOError);
      }
      edges = uniformBinEdges(binsMin, binsMax, int(nBins));
    } else {
      const scalar width = bd.get<scalar>("width");
      if (!(width > 0) || !((binsMax - binsMin)/width <= maxBins)) {
        FatalIOErrorInFunction(bd) << "profiles: bin set " << e.keyword()
          << " needs width > 0 and at most " << maxBins << " bins, i.e. "
          << "width >= " << (binsMax - binsMin)/maxBins << " m."
          << exit(FatalIOError);
      }
      edges = widthBinEdges(binsMin, binsMax, width);
    }

    sets_.set(seti, new binSet());
    binSet& set = sets_[seti];
    set.name = e.keyword();
    set.edges.resize(label(edges.size()));
    forAll(set.edges, k) {
      set.edges[k] = edges[k];
    }
    const label nBins = set.nBins();

    binIntervals (
      set.edges, cellLo, cellHi,
      set.cellStart, set.cellBin, set.cellWeight, set.cellOutside
    );
    binIntervals (
      set.edges, elementLo, elementHi,
      set.elementStart, set.elementBin, set.elementWeight, set.elementOutside
    );

    /* the overlap weights of the wall faces (for layer faceCells the same
       faces, in the same order, as the elements) */
    labelList wallStart, wallBin;
    scalarField wallWeight, wallOutside;
    binIntervals(set.edges, wallLo, wallHi, wallStart, wallBin, wallWeight,
                 wallOutside);

    /* frozen wall area and area-weighted T per bin, and the coverage */
    scalarField frozen(2*nBins + 2, 0.0);
    forAll(wallPatch, i) {
      for (label k = wallStart[i]; k < wallStart[i + 1]; ++k) {
        const scalar area = wallWeight[k]*wallFaceArea[i];
        frozen[wallBin[k]] += area;
        frozen[nBins + wallBin[k]] += area*wallFaceT[i];
      }
      frozen[2*nBins + 1] += wallOutside[i]*wallFaceArea[i];
    }
    forAll(V, celli) {
      frozen[2*nBins] += set.cellOutside[celli]*V[celli];
    }
    Pstream::listCombineReduce(frozen, plusEqOp<scalar>());

    set.wallArea.resize(nBins);
    set.wallT.resize(nBins);
    forAll(set.wallT, b) {
      set.wallArea[b] = frozen[b];
      set.wallT[b] = set.wallArea[b] > 0
        ? frozen[nBins + b]/set.wallArea[b]
        : 0;
    }

    /* first bin at or after searchFrom (fromInlet) */
    set.firstSearchBin = 0;
    if (haveSearchFrom_) {
      label& b = set.firstSearchBin;
      while (b < nBins && 0.5*(set.edges[b] + set.edges[b + 1]) < searchFrom_) {
        ++b;
      }
    }

    const scalar lastWidth = set.edges[nBins] - set.edges[nBins - 1];
    setsText << (seti ? "; " : "") << set.name << " " << nBins
      << " bins over [" << binsMin << ", " << binsMax << "] m, width "
      << set.edges[1] - set.edges[0] << " m";
    if (nBins > 1 && mag(lastWidth - (set.edges[1] - set.edges[0]))
                   > 1e-9*lastWidth) {
      setsText << " (last bin " << lastWidth << " m)";
    }
    const scalar volumeOutside = frozen[2*nBins]/meshVolume;
    const scalar areaOutside =
      interfaceArea > 0 ? frozen[2*nBins + 1]/interfaceArea : 0;
    if (volumeOutside > 0 || areaOutside > 0) {
      setsText << " (outside the bins: " << volumeOutside*100
        << " % of the mesh volume, " << areaOutside*100
        << " % of the interface area; their amounts are reported as "
        << "'outside')";
    }

    ++seti;
  }

  /*--------------------------------------------------------------------------
  Profile times: ascending and distinct.  On a restart with the record of
  the entries written by the runs before (uniform/phaseChangeLayout), those
  entries are skipped exactly, and so are the entries at or before the
  start that they did not write (e.g. added to 'times' since); every other
  entry is written by this run, whatever the deltaT of either run.  On a
  fresh start (or a restart from a layout of an earlier build, without the
  record) an entry whose nearest step end is the start state or earlier
  (t_i - deltaT/2 <= tStart) is skipped: the initial state has profiles
  only if Time writes it (so an entry up to deltaT/2 after a fresh start is
  skipped as well).
  --------------------------------------------------------------------------*/
  scalarList times(dict.getOrDefault<scalarList>("times", scalarList()));
  std::sort(times.begin(), times.end());
  DynamicList<scalar> distinct(times.size());
  for (const scalar t : times) {
    if (distinct.empty() || t != distinct.back()) {
      distinct.append(t);
    }
  }
  times_.transfer(distinct);
  done_.resize(times_.size(), false);

  label nWrittenBefore = 0, nBeforeStart = 0;
  if (written) {
    written_ = *written;
    forAll(times_, i) {
      if (written->found(times_[i])) {
        done_[i] = true;
        ++nWrittenBefore;
      } else if (times_[i] <= tStart) {
        done_[i] = true;
        ++nBeforeStart;
      }
    }
  } else {
    forAll(times_, i) {
      if (times_[i] - 0.5*deltaT <= tStart) {
        done_[i] = true;
        ++nBeforeStart;
      }
    }
  }
  const label nPending = times_.size() - nWrittenBefore - nBeforeStart;

  /*--------------------------------------------------------------------------
  Start-up line
  --------------------------------------------------------------------------*/
  Info<< "Profiles: axis " << axis_ << ", origin " << origin_
    << ", the mesh spans [" << meshLo << ", " << meshHi << "] m along it; "
    << setsText.str().c_str() << "; written at the write times and at "
    << nPending << " profile time(s)";
  if (written && nWrittenBefore + nBeforeStart > 0) {
    Info<< " (";
    if (nWrittenBefore > 0) {
      Info<< nWrittenBefore << " skipped: written by the run that wrote the "
        << "restart state" << (nBeforeStart > 0 ? "; " : "");
    }
    if (nBeforeStart > 0) {
      Info<< nBeforeStart << " skipped: at or before the start time, not "
        << "written by that run";
    }
    Info<< ")";
  } else if (nBeforeStart > 0) {
    Info<< " (" << nBeforeStart << " skipped: their nearest step end is not "
      << "after the start time"
      << (restarted ? "; the restart state has no record of the profile "
          "times written (an earlier build)" : "") << ")";
  }
  Info<< "; T_dep at " << fraction_ << " of the peak of " << observable
    << (wallPlusGas_ ? " (wall + gas outside the sample cells)" : "")
    << " (significant when the peak bin holds more than " << significance_
    << " x the pair inventory), search";
  if (fromInlet_) {
    Info<< " fromInlet";
    if (haveSearchFrom_) {
      Info<< " (from x = " << searchFrom_ << " m)";
    }
  }
  if (fromPeak_) {
    Info<< " fromPeak";
  }
  Info<< nl;
}


/*------------------------------------------------------------------------------
Profile times reached by the state at the end of a step
------------------------------------------------------------------------------*/

Foam::string LESTO::axialProfiles::reached (
  const scalar stepEnd,
  const scalar deltaT
) {
  string entries;

  forAll(times_, i) {
    if (!done_[i] && stepEnd >= times_[i] - 0.5*deltaT) {
      entries += string(entries.empty() ? "" : "; ") + "profile time "
        + text(times_[i]);
      done_[i] = true;
      written_.append(times_[i]);
    }
  }

  return entries;
}


Foam::scalarList LESTO::axialProfiles::written() const {
  scalarList out(written_);
  std::sort(out.begin(), out.end());
  return out;
}


/*------------------------------------------------------------------------------
Onset of one search mode, as text for the log or the file
------------------------------------------------------------------------------*/

Foam::string LESTO::axialProfiles::onsetText (
  const binSet&      set,
  const scalarField& observable,
  const bool         fromPeak,
  const bool         forFile,
  const scalar       minimumAmount
) const {

  const label nBins = set.nBins();
  const label first = fromPeak ? 0 : set.firstSearchBin;
  scalarField widths(nBins);
  forAll(widths, b) {
    widths[b] = set.edges[b + 1] - set.edges[b];
  }
  const onsetCrossing c = onsetSearch (
    observable.cdata(), int(nBins), int(first), fromPeak, fraction_,
    widths.cdata(), minimumAmount
  );

  auto centre = [&](const label b) {
    return 0.5*(set.edges[b] + set.edges[b + 1]);
  };

  scalar x = 0, T = 0;
  bool haveT = false;
  if (c.status == ONSET_FOUND) {
    const label lo = c.lower;
    x = centre(lo) + c.weight*(centre(lo + 1) - centre(lo));
    haveT = set.wallArea[lo] > 0 && set.wallArea[lo + 1] > 0;
    if (haveT) {
      T = set.wallT[lo] + c.weight*(set.wallT[lo + 1] - set.wallT[lo]);
    }
  }

  OStringStream os;
  const word mode(fromPeak ? "fromPeak" : "fromInlet");

  if (forFile) {
    os.precision(17);
    os << mode << " status " << statusWord(c.status);
    if (c.status == ONSET_FOUND) {
      os << " T ";
      if (haveT) {
        os << T;
      } else {
        os << "none";
      }
      os << " x " << x;
    }
    if (c.peakBin >= 0) {
      os << " peak " << c.peak << " xPeak " << centre(c.peakBin)
        << " threshold " << c.threshold;
    }
    if (!fromPeak && haveSearchFrom_) {
      os << " searchFrom " << searchFrom_;
    }
    return os.str();
  }

  os << mode;
  if (!fromPeak && haveSearchFrom_) {
    os << " (from x = " << searchFrom_ << " m)";
  }
  if (c.status == ONSET_NO_PEAK) {
    os << " none (no " << (wallPlusGas_ ? "deposit or gas" : "deposit") << ")";
  } else if (c.status == ONSET_INSIGNIFICANT) {
    os << " none (insignificant: the peak bin holds "
       << c.peak*widths[c.peakBin] << " mol, at most " << significance_
       << " x the pair inventory = " << minimumAmount << " mol; peak "
       << c.peak << " mol/m at x = " << centre(c.peakBin) << " m)";
  } else if (c.status == ONSET_NOT_BRACKETED) {
    os << " none (at or above the threshold from the first bin of the "
       << "search; peak " << c.peak << " mol/m at x = " << centre(c.peakBin)
       << " m)";
  } else {
    if (haveT) {
      os << " " << T << " K";
    } else {
      os << " (no wall face in the bins)";
    }
    os << " at x = " << x << " m (peak " << c.peak << " mol/m at x = "
       << centre(c.peakBin) << " m)";
  }
  return os.str();
}


/*------------------------------------------------------------------------------
Profiles of one pair
------------------------------------------------------------------------------*/

void LESTO::axialProfiles::write (
  const word&        timeName,
  const scalar       time,
  const string&      reason,
  const word&        gasName,
  const word&        condensateName,
  const scalar       molarMass,
  const scalarField& gasMass,
  const scalarField& wallMass,
  const scalarField& sampleMass,
  const FixedList<scalar, 3>& inventory,
  const bool writePeak,
  const boolList* sampleCells
) const {

  const scalar M = molarMass;
  const char* names[3] = {"gas", "wall", "sample"};

  for (const binSet& set : sets_) {

    const label nBins = set.nBins();

    /*------------------------------------------------------------------------
    Binned amounts [kg] of the gas, the wall, the sample and the gas of the
    sample cells per bin, followed by the four amounts outside the bins,
    summed over the ranks in one reduction of a list sized by the bin set
    only.
    ------------------------------------------------------------------------*/
    scalarField values(4*nBins + 4, 0.0);

    forAll(gasMass, celli) {
      const bool inSample = sampleCells ? (*sampleCells)[celli] : sampleCell_[celli];
      for (label k = set.cellStart[celli]; k < set.cellStart[celli + 1]; ++k) {
        values[set.cellBin[k]] += set.cellWeight[k]*gasMass[celli];
        if (inSample) {
          values[3*nBins + set.cellBin[k]] += set.cellWeight[k]*gasMass[celli];
        }
      }
      values[4*nBins] += set.cellOutside[celli]*gasMass[celli];
      if (inSample) {
        values[4*nBins + 3] += set.cellOutside[celli]*gasMass[celli];
      }
    }
    forAll(wallMass, e) {
      for (label k = set.elementStart[e]; k < set.elementStart[e + 1]; ++k) {
        values[nBins + set.elementBin[k]] += set.elementWeight[k]*wallMass[e];
      }
      values[4*nBins + 1] += set.elementOutside[e]*wallMass[e];
    }
    forAll(sampleMass, celli) {
      for (label k = set.cellStart[celli]; k < set.cellStart[celli + 1]; ++k) {
        values[2*nBins + set.cellBin[k]] += set.cellWeight[k]*sampleMass[celli];
      }
      values[4*nBins + 2] += set.cellOutside[celli]*sampleMass[celli];
    }

    Pstream::listCombineReduce(values, plusEqOp<scalar>());

    /*------------------------------------------------------------------------
    Line densities [mol/m], the observable, the totals [mol] in the bins and
    outside, and their relative difference from the inventory.
    ------------------------------------------------------------------------*/
    scalarField lineDensity(4*nBins), observable(nBins);
    FixedList<scalar, 3> binned, outside, relative;

    for (label q = 0; q < 4; ++q) {
      compensatedSum sum;
      for (label b = 0; b < nBins; ++b) {
        const scalar dx = set.edges[b + 1] - set.edges[b];
        lineDensity[q*nBins + b] = values[q*nBins + b]/(M*dx);
        sum.add(values[q*nBins + b]);
      }
      if (q == 3) {
        break;
      }
      binned[q] = sum.value()/M;
      outside[q] = values[4*nBins + q]/M;
      const scalar difference =
        sum.value() + values[4*nBins + q] - inventory[q];
      relative[q] = difference == 0
        ? 0
        : difference/max(mag(inventory[q]), VSMALL);
    }

    /* wallPlusGas: the deposit and the gas outside the sample cells (in a
       bin without sample cells the gas itself, exactly) */
    for (label b = 0; b < nBins; ++b) {
      const scalar sampleGas = lineDensity[3*nBins + b];
      observable[b] = lineDensity[nBins + b]
        + (
            wallPlusGas_
          ? (sampleGas != 0 ? lineDensity[b] - sampleGas : lineDensity[b])
          : 0
          );
    }

    /* the significance of an onset [mol]: of the pair inventory */
    const scalar minimumAmount =
      significance_*(inventory[0] + inventory[1] + inventory[2])/M;

    /*------------------------------------------------------------------------
    Log
    ------------------------------------------------------------------------*/
    Info<< "Profiles " << gasName << ", set " << set.name << " ("
      << reason.c_str() << "): in the bins";
    for (label q = 0; q < 3; ++q) {
      Info<< ' ' << names[q] << ' ' << binned[q];
    }
    Info<< " mol, outside " << outside[0] << ' ' << outside[1] << ' '
      << outside[2] << " mol; relative difference of bins + outside from "
      << "the inventory " << relative[0] << ' ' << relative[1] << ' '
      << relative[2] << nl
      << "  T_dep of " << (wallPlusGas_ ? "wallPlusGas" : "wall") << " at "
      << fraction_ << " of the peak:";
    if (fromInlet_) {
      Info<< ' '
        << onsetText(set, observable, false, false, minimumAmount).c_str()
        << (fromPeak_ ? ";" : "");
    }
    if (fromPeak_) {
      Info<< ' '
        << onsetText(set, observable, true, false, minimumAmount).c_str();
    }
    Info<< nl;

    /*------------------------------------------------------------------------
    File (master only)
    ------------------------------------------------------------------------*/
    if (!Pstream::master()) {
      continue;
    }

    const Time& runTime = mesh_.time();
    const fileName dir
    (
      runTime.globalPath()/"postProcessing"/"axialProfiles"/timeName
    );
    mkDir(dir);

    OFstream os(dir/(gasName + "_" + set.name + ".dat"));
    os.precision(17);

    os << "# axial profiles of the pair " << gasName << " -> "
       << condensateName << ", bin set " << set.name
       << " (rhoFixedFlowFoam, axialProfiles.H)" << nl
       << "# time " << time << " directory " << timeName << " reason "
       << reason.c_str() << nl
       << "# molarMass " << M << nl
       << "# axis " << axis_ << " origin " << origin_ << nl
       << "# nBins " << nBins << nl;

    const char* keys[3] = {"inventory", "binned", "outside"};
    for (label row = 0; row < 3; ++row) {
      os << "# " << keys[row];
      for (label q = 0; q < 3; ++q) {
        const scalar value = row == 0 ? inventory[q]/M
          : row == 1 ? binned[q] : outside[q];
        os << ' ' << names[q] << ' ' << value;
      }
      os << nl;
    }

    os << "# observable " << (wallPlusGas_ ? "wallPlusGas" : "wall") << nl
       << "# fraction " << fraction_ << nl
       << "# significance " << significance_ << " minimumAmount "
       << minimumAmount << nl;
    if (fromInlet_) {
      os << "# Tdep "
         << onsetText(set, observable, false, true, minimumAmount).c_str()
         << nl;
    }
    if (fromPeak_) {
      os << "# Tdep "
         << onsetText(set, observable, true, true, minimumAmount).c_str()
         << nl;
    }
    if (writePeak) {
      label peak = 0;
      for (label b = 1; b < nBins; ++b) {
        if (observable[b] > observable[peak]) peak = b;
      }
      if (observable[peak]*(set.edges[peak+1] - set.edges[peak]) > minimumAmount
        && set.wallArea[peak] > 0) {
        os << "# Tpeak " << set.wallT[peak] << " K at x "
          << 0.5*(set.edges[peak] + set.edges[peak+1])
          << " m; bin maximum of observable" << nl;
      } else os << "# Tpeak none (no significant peak with wall temperature)" << nl;
    }
    os << "# columns x_lo x_hi x gas wall sample observable Twall wallArea"
       << " sampleGas" << nl
       << "# units m m m mol/m mol/m mol/m mol/m K m2/m mol/m" << nl;

    for (label b = 0; b < nBins; ++b) {
      const scalar dx = set.edges[b + 1] - set.edges[b];
      os << set.edges[b] << ' ' << set.edges[b + 1] << ' '
         << 0.5*(set.edges[b] + set.edges[b + 1]) << ' '
         << lineDensity[b] << ' ' << lineDensity[nBins + b] << ' '
         << lineDensity[2*nBins + b] << ' ' << observable[b] << ' '
         << set.wallT[b] << ' ' << set.wallArea[b]/dx << ' '
         << lineDensity[3*nBins + b] << nl;
    }
  }
}
