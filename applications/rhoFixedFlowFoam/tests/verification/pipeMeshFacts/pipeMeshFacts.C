/*------------------------------------------------------------------------------

pipeMeshFacts.C

PURPOSE

Test helper of tests/verification (not part of the solver).  Prints the
facts of a pipe mesh that the M7 refinement studies of rhoFixedFlowFoam
compare (doc/phase-change-plan.md, milestone M7: the h vs h/2 first-layer
mesh of acceptance 2 and the axial x2 mesh; tests/verification/pipeMeshes
and meshChecks.py).  Serial; nothing is written.  One 'key value ...' line
per fact:

  - the cells, and the wall faces of the interface patches;
  - the interaction layer 'layer faceCells' of the solver: elements, cells,
    A_I = sum |S_f|, V_I = the volume of the face cells and A_I/V_I (plan
    E2), and the per-element |S_f|/V_c (areaPerVolume cell);
  - the first cell layer: the half-cell distance d_e = 1/deltaCoeffs of
    the wall face (what wallResistance halfCell uses), the wall-normal
    thickness of the cell from its points, V_c/|S_f|, and the centre
    distance R - r_c;
  - the axial cells: their axial extent from the points and the number of
    distinct axial positions of the cell centres;
  - the cell layers counted from the wall by face neighbours (layer 1 = the
    face cells, layer k+1 = the face neighbours of layer k that are in no
    earlier layer; on an O-grid these are the rings of cells), with their
    centre distances R - r_c and patchWave distances;
  - the cells within 'distance' of the wall by the centre distance R - r_c
    <= distance (T-Flows: 'centre distance <= 0.1 mm') and by the
    patchWave distance <= distance (corrected near the wall), with their
    volume and A_I/V_layer.  'layer distance' of milestone M4 uses the
    exact distance from the nearest wall face (plan section 30, item 1),
    which selects the same 30528 cells on the base mesh (criterion M4.2).

Usage

  pipeMeshFacts [-case dir] [-patches '(WALL)'] [-distance 1e-4]
                [-radius 2.4e-3] [-axis '(1 0 0)']

------------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "patchWave.H"
#include "IOmanip.H"

namespace {

  /* the minimum, mean and maximum of a sequence of values */
  struct stats {

    scalar min = GREAT;
    scalar max = -GREAT;
    scalar sum = 0;
    label n = 0;

    void add(const scalar v) {
      min = Foam::min(min, v);
      max = Foam::max(max, v);
      sum += v;
      ++n;
    }

    scalar mean() const {
      return n ? sum/n : 0;
    }
  };

  void print(const word& key, const stats& s) {
    Info<< key << " min " << s.min << " mean " << s.mean()
      << " max " << s.max << " n " << s.n << nl;
  }
}


int main(int argc, char *argv[]) {

  argList::addNote (
    "Facts of a pipe mesh for the M7 refinement studies"
    " (interaction layers, first-layer thickness, axial cells)"
  );
  argList::noParallel();
  argList::addOption("patches", "wordRes", "interface patches (default (WALL))");
  argList::addOption("distance", "scalar", "layer distance [m] (default 1e-4)");
  argList::addOption("radius", "scalar", "pipe radius [m] (default 2.4e-3)");
  argList::addOption (
    "axis", "vector", "pipe axis through the origin (default (1 0 0))"
  );

  #include "setRootCase.H"
  #include "createTime.H"
  #include "createMesh.H"

  wordRes patchNames(1, wordRe("WALL"));
  args.readListIfPresent<wordRe>("patches", patchNames);
  const scalar delta = args.getOrDefault<scalar>("distance", 1e-4);
  const scalar R = args.getOrDefault<scalar>("radius", 2.4e-3);
  const vector axis (
    normalised(args.getOrDefault<vector>("axis", vector(1, 0, 0)))
  );

  const labelHashSet patchIDs(mesh.boundaryMesh().patchSet(patchNames));
  if (patchIDs.empty()) {
    FatalErrorInFunction << "No patch matches " << patchNames
      << exit(FatalError);
  }

  const volVectorField& C = mesh.C();
  const scalarField& V = mesh.V();
  const pointField& points = mesh.points();
  const labelListList& cellPoints = mesh.cellPoints();
  const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();

  auto radial = [&axis](const vector& p) {
    return mag(p - (p & axis)*axis);
  };

  Info<< setprecision(10);
  Info<< nl << "cells " << mesh.nCells() << " points " << mesh.nPoints()
    << " faces " << mesh.nFaces() << nl;
  Info<< "axis " << axis << " radius " << R << " distance " << delta << nl;

  /*--------------------------------------------------------------------------
  The interaction layer 'faceCells' and the first cell layer
  --------------------------------------------------------------------------*/
  label nElements = 0;
  scalar A_I = 0;
  boolList faceCell(mesh.nCells(), false);
  stats areaPerVolume, halfCell, thickness, volumePerArea, centreDistance;

  for (const label patchi : patchIDs.sortedToc()) {

    const fvPatch& p = mesh.boundary()[patchi];
    const labelUList& fc = p.faceCells();
    const vectorField nf(p.nf());
    const vectorField& Cf = p.Cf();
    const scalarField& magSf = p.magSf();
    const scalarField& dc = deltaCoeffs.boundaryField()[patchi];

    Info<< "patch " << p.name() << " faces " << p.size()
      << " area " << sum(magSf) << nl;

    forAll(fc, facei) {
      const label celli = fc[facei];
      ++nElements;
      A_I += magSf[facei];
      faceCell[celli] = true;

      areaPerVolume.add(magSf[facei]/V[celli]);
      halfCell.add(1.0/dc[facei]);
      volumePerArea.add(V[celli]/magSf[facei]);
      centreDistance.add(R - radial(C[celli]));

      scalar extent = 0;
      for (const label pointi : cellPoints[celli]) {
        extent = max(extent, (Cf[facei] - points[pointi]) & nf[facei]);
      }
      thickness.add(extent);
    }
  }

  label nFaceCells = 0;
  scalar V_I = 0;
  forAll(faceCell, celli) {
    if (faceCell[celli]) {
      ++nFaceCells;
      V_I += V[celli];
    }
  }

  Info<< "faceCells elements " << nElements << " cells " << nFaceCells
    << " A_I " << A_I << " V_I " << V_I << " A_I/V_I " << A_I/V_I
    << " V_I/A_I " << V_I/A_I << nl;
  print("faceCells |S_f|/V_c", areaPerVolume);
  print("firstLayer halfCellDistance(1/deltaCoeffs)", halfCell);
  print("firstLayer wallNormalThickness", thickness);
  print("firstLayer V_c/|S_f|", volumePerArea);
  print("firstLayer centreDistance(R-r_c)", centreDistance);

  /*--------------------------------------------------------------------------
  The axial cells: the axial extent of every cell from its points, and the
  distinct axial positions of the cell centres
  --------------------------------------------------------------------------*/
  stats axialExtent;
  DynamicList<scalar> xc(mesh.nCells());
  forAll(C, celli) {
    scalar lo = GREAT, hi = -GREAT;
    for (const label pointi : cellPoints[celli]) {
      const scalar x = points[pointi] & axis;
      lo = min(lo, x);
      hi = max(hi, x);
    }
    axialExtent.add(hi - lo);
    xc.append(C[celli] & axis);
  }
  Foam::sort(xc);
  label nAxial = xc.size() ? 1 : 0;
  for (label i = 1; i < xc.size(); ++i) {
    if (xc[i] - xc[i-1] > 1e-3*axialExtent.min) {
      ++nAxial;
    }
  }
  print("axialExtent", axialExtent);
  Info<< "axialPositions " << nAxial << " cellsPerAxialPosition "
    << scalar(mesh.nCells())/max(nAxial, 1) << nl;

  /*--------------------------------------------------------------------------
  The patchWave distance (corrected near the wall) and the cell layers
  counted from the wall by face neighbours
  --------------------------------------------------------------------------*/
  const patchWave wave(mesh, patchIDs, true);
  const scalarField& dist = wave.distance();

  labelList layer(mesh.nCells(), -1);
  DynamicList<label> front(nFaceCells);
  forAll(faceCell, celli) {
    if (faceCell[celli]) {
      layer[celli] = 1;
      front.append(celli);
    }
  }
  const labelListList& cellCells = mesh.cellCells();
  const label nLayersShown = 6;
  for (label k = 1; k <= nLayersShown && front.size(); ++k) {

    stats cd, wd;
    scalar Vk = 0;
    for (const label celli : front) {
      cd.add(R - radial(C[celli]));
      wd.add(dist[celli]);
      Vk += V[celli];
    }
    Info<< "layer " << k << " cells " << front.size() << " volume " << Vk
      << " centreDistance min " << cd.min << " max " << cd.max
      << " patchWaveDistance min " << wd.min << " max " << wd.max
      << " withinDistance(centre) "
      << (cd.max <= delta ? "all" : cd.min <= delta ? "some" : "none") << nl;

    DynamicList<label> next;
    for (const label celli : front) {
      for (const label nbr : cellCells[celli]) {
        if (layer[nbr] < 0) {
          layer[nbr] = k + 1;
          next.append(nbr);
        }
      }
    }
    front.transfer(next);
  }

  /*--------------------------------------------------------------------------
  The cells within 'distance' of the wall, by both distances
  --------------------------------------------------------------------------*/
  label nCentre = 0, nWave = 0;
  scalar Vcentre = 0, Vwave = 0;
  labelList centreLayers(nLayersShown + 2, Zero);
  forAll(C, celli) {
    if (R - radial(C[celli]) <= delta) {
      ++nCentre;
      Vcentre += V[celli];
      const label k = layer[celli] < 0 ? nLayersShown + 1 : layer[celli];
      ++centreLayers[min(k, nLayersShown + 1)];
    }
    if (dist[celli] <= delta) {
      ++nWave;
      Vwave += V[celli];
    }
  }
  Info<< "distanceLayer centre(R-r_c<=" << delta << ") cells " << nCentre
    << " layers " << scalar(nCentre)/max(nElements, 1)
    << " volume " << Vcentre << " A_I/V_layer "
    << (Vcentre > 0 ? A_I/Vcentre : 0) << nl;
  Info<< "distanceLayer patchWave(<=" << delta << ") cells " << nWave
    << " layers " << scalar(nWave)/max(nElements, 1)
    << " volume " << Vwave << " A_I/V_layer "
    << (Vwave > 0 ? A_I/Vwave : 0) << nl;
  Info<< "distanceLayer centre cells per layer index (1.." << nLayersShown
    << ", beyond)";
  for (label k = 1; k < centreLayers.size(); ++k) {
    Info<< ' ' << centreLayers[k];
  }
  Info<< nl << nl << "End" << nl;

  return 0;
}
