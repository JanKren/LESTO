/*------------------------------------------------------------------------------

interfaceExchange.H

PURPOSE

Gas/condensate exchange of rhoFixedFlowFoam (doc/phase-change-plan.md,
sections 2 and 3).  One concrete class, LESTO::interfaceExchange, reads
constant/thermochemistryProperties and owns everything the exchange needs:

  - the gas/condensate PAIRS, e.g. PbI2_g -> PbI2_s.  Only a paired gas is
    coupled to a condensate; every other species keeps the unchanged solver
    path (solveGasSpecies.H);
  - the INTERFACE ELEMENTS of the interaction layer.  layer faceCells: one
    element per face of the listed wall patches, with its cell, its area
    A_e (times areaMultiplier; areaPerVolume layer rescales it to the layer
    mean A_I/V_I), its temperature T_e (cell or wall face) and its deposit
    m''_e [kg/m2].  layer distance (milestone M4): one element per cell
    whose centre lies within 'distance' of the wall patches, of the area
    A_e = V_c (A_I/V_I) areaMultiplier, with its nearest wall face for the
    output (see the fifth comment block);
  - the SAMPLES, keyed to a pair; their cells are selected in memory
    (cellBitSet) or taken from an existing cellZone.  A release sample
    injects a given amount over a finite time window (mode release); an
    inventory sample (mode inventory, model HKS) holds a condensate
    reservoir Cs [kg/m3] in its cells that evaporates by the HKS law and is
    removed at removeTime (see the third comment block);
  - per pair: the linearised gas source S = Su - Sp*Y (volScalarField::
    Internal, Su >= 0 and Sp >= 0), the regime of every element, the bounds
    of the step, and the cumulative mass ledger (phaseChangeLedger).

Models

  model mock;          no exchange: the solver runs exactly as the pristine
                       build of tag 0.1, the coupled mock of the main branch
                       (also without the file).  Nothing is read, created
                       or printed here.
  model temperature;   Lobresco et al. (2026) Eq. 5 in flux form,
                       mdot''_e = -rho_c v_d f(T_e) Y_c, irreversible.
  model HKS;           Hertz-Knudsen-Schrage (plan E2), SI:
                       mdot''_e = G_e (p_eq(T_e) - beta_c Y_c), with
                       G_e = kineticScale Ce 2 sigma/(2 - sigma)
                             sqrt(M/(2 pi R T_e))  [s/m],
                       beta_c = rho_c R T_c/M (p_i = beta_c Y_c from the
                       carrier's actual rho) and p_eq from a table
                       (vapourPressureTable.H) or from GEMS3K
                       (equilibrium GEMS, gemsEquilibrium.H; the sixth
                       comment block); both directions where there is a
                       condensate.

Every element carries the affine law q_e = a_e - g_e*Y [kg/s] and the bounds
lower_e <= q_e <= upper_e [kg/s] of the step (plan E3):

  wall element, temperature   a_e = 0, g_e = A_e rho_c v_d f(T_e);
                              upper = 0 (irreversible)
  wall element, HKS           a_e = A_e G_e p_eq(T_e), g_e = A_e G_e beta_c;
                              upper = m''_e A_e/deltaT (a bare wall only
                              deposits: upper = 0)
  sample element, HKS         one per cell of an inventory sample, of the
                              area A_e = areaPerVolume V_c and the sample's
                              own kineticScale and Ce; lower = 0 (no
                              deposition onto the sample), upper =
                              Cs_c V_c/deltaT
  all                         lower = -infinity on the wall

The wall elements come first (sorted by cell), then the sample elements,
grouped by exchange cell. Different pairs may share sample cells; samples
of the same pair must be disjoint. Each sample retains its own element map.

Step algorithm (plan section 3)

  beginStep()      the exact step [t, t + deltaT] on the clock of the ledger,
                   the removal of inventory samples (t >= removeTime -
                   deltaT/2), release rates (exact overlap of the step with
                   the release window), element bounds, per-step counters,
                   and a check
                   of the final solver dictionary, of the Euler ddt scheme
                   and of the conservative (Gauss) convection scheme
                   (fvSolution and fvSchemes may have been modified during
                   the run)
  linearise()      at every outer corrector: the local implicit predictor
                   (phaseChangeKinetics.H, cellPredictor) from the transport-
                   only matrix, starting from the regimes of the previous
                   corrector (a sample element never jumps between CAPPED
                   and NONE; it passes IMPLICIT), Su and Sp of every
                   exchange cell, and the regime changes (number, mass,
                   significant ones; reduced over the ranks)
  outerConverged() after an outer corrector: initial residual below the
                   tolerance and the regimes settled (see below); keeps the
                   regime changes of the deciding corrector for the log
  recordSolve()    after the final solve: transport through every
                   non-coupled patch (fvMatrix::flux()) and the per-cell
                   solver defect of the transport-only matrix in flux form
                   plus the realised element flows and the release
                   (matrixDefect.H)
  guard()          the bounded law at the solved Y: the realised flow of
                   every element against the admissible flow min(max(a - g
                   Y, lower), upper); if the violations exceed
                   guardTolerance*inventory in total, violating IMPLICIT
                   elements are switched to CAPPED or NONE and violating
                   CAPPED or NONE elements back to IMPLICIT, the largest
                   misses first until the misses left are within the
                   tolerance, and the equation is re-solved (at most
                   nGuard times; see below)
  realise()        conservative same-cell projection of what is left
                   (round-off), condensate update, CAPPED residue check,
                   complementarity check, monitors, log lines
  balance()        ledger and closure, balance file, derived Y_<condensate>,
                   the deposit fields mw_<condensate> and mDep_<condensate>,
                   the optional exchange field exch_<gas>, and the element
                   state of uniform/ (below); all before runTime.write()
  writeProfiles()  after runTime.write(): the axial profiles and T_dep
                   (axialProfiles.H) at the profile times reached by the
                   step, into postProcessing/axialProfiles only (M3); at a
                   write of Time runTime.write() writes them (below)

Settled regimes (a deviation from plan E4, recorded in its section 14).
The outer loop is converged when the initial residual is below the
tolerance and the regime changes of the latest predictor are settled: the
flow they change, times deltaT, is at most guardTolerance*inventory, the
threshold at which the guard accepts a bound violation.  Sign changes of Y
at round-off level (e.g. 1e-40 ahead of a gas front) flip elements between
IMPLICIT and NONE without carrying mass; counting them would never let the
outer loop converge.  After convergence the Final solve (Y_<gas>Final) runs
its own predictor; the log prints the changes of the deciding corrector
and those of this final predictor separately.

Significant counts.  Regime changes, projections and complementarity
violations are counted in total and as 'significant': an element whose
mass (changed flow, projected flow or missed flow, times deltaT) exceeds
guardTolerance*inventory on its own.  Noise-level elements (|Y| ~ 1e-40)
thus cannot hide a real event, and 'significant 0' is a meaningful check.
A complementarity violation is an element that is not IMPLICIT and whose
realised flow differs from the admissible flow of a - g Y at the final Y,
e.g. a CAPPED sample cell in supersaturated gas.  The masses of a step are
also compared in total: a projected or missed mass above guardTolerance*
inventory is named in the exchange line ('... kg > tolerance ... kg') and
the steps are counted in the end-of-run summary (review of M5, round 3;
several sub-threshold misses had added up to 1.43 times the tolerance
with '0 significant').

Settling a block of sample cells (review of M5, round 2).  The predictor
lags the neighbour cells by one corrector (H).  In the cells of an
inventory sample with a small reservoir, every cell switched at once: all
CAPPED (the neighbours are empty), then all NONE (the neighbours are
supersaturated), and so on at every corrector; the last corrector's parity
decided the step, and neither state obeyed the law.  Two rules remove it.
The predictor never switches a sample element directly between its two
bounds: CAPPED <-> NONE passes IMPLICIT, whose coupled solve resolves the
block (acrossBounds).  The guard enforces the law in both directions at the
solved Y, where only the adjacent switch back to IMPLICIT is safe (with
lambda deltaT of 1e3 to 1e4, the gas of a cell CAPPED and of the same cell
NONE differ by far more than the width of its bounds).  It is triggered by
the total violation and acts on the total (review of M5, round 3): the
CAPPED or NONE elements go back to IMPLICIT with the largest misses first,
until the misses left add up to at most guardTolerance*inventory
(switchBackThreshold, phaseChangeKinetics.H); round 2 switched back only an
element whose own miss was significant.  The outer loop then settles
within a few correctors whatever nCorr is, and a step is accepted with a
total violation above the tolerance only when the guard passes are
exhausted (reported in the exchange line and the end-of-run summary).

(continued in the next comment block)

------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

interfaceExchange.H, continued.  The description is split into several
comment blocks on purpose: wmkdepend of OpenFOAM v2412 reads a comment as
one token into a 16 KB buffer, and at a longer comment it stops scanning the
file ('wmkdepend: buffer full'), so the #include lines below would be
missing from the dependencies.  Keep every comment block below 16 KB.

State, restart and output (the state lives in plain lists; what is written
is refreshed in every step by balance(), before runTime.write())

  uniform/phaseChangeBalance   global: the ledger (phaseChangeLedger.H),
  uniform/phaseChangeLayout    global: the names of its pairs and patches,
                   and the restart record of the interface and the samples
                   (below)
  uniform/phaseChange/phaseChangeDeposits  per rank: the deposit m''_e
                   [kg/m2] of every element of every pair, the RESTART
                   STATE of the condensate
  uniform/phaseChange/phaseChangeRegimes  per rank: the regime of every
                   element of every pair (review item B5), so that a binary
                   restart repeats the continuous run bit for bit
  uniform/phaseChange/phaseChangeReservoirs  per rank: the reservoir Cs of
                   every sample element of every pair (inventory samples
                   only; third comment block)
  uniform/phaseChange/phaseChangeEquilibrium  per rank, equilibrium GEMS
                   in mode local only: the lagged state of the last update
                   (p_eq, chi and the source of every element of every
                   pair, and the counts of that update; the sixth comment
                   block), so that a restart with updateInterval > 1
                   repeats the continuous run up to its next update

  The ledger (binaryGlobalIOList) and the per-rank lists
  (binaryIOList; binaryIOList.H) are registered AUTO_WRITE, so that every
  write of Time stores them with the fields (write times, and
  Time::writeNow()/writeAndEnd() of function objects between two steps),
  and always written in BINARY, as raw blocks that read back bit for bit
  whatever writeFormat and writePrecision are; the layout is a global
  dictionary registered AUTO_WRITE as well.  The global objects are the same
  on every rank: with the collated file handler the master writes them as
  plain files, which reconstructPar and decomposePar carry to any number of
  ranks.  Each per-rank list starts with a fingerprint of its elements on
  the rank (their number and a 62-bit hash; elementsFingerprint()): the
  deposits list that of the wall elements (cells, patch faces and face
  centres, as exact scalars), the regimes list that of all elements, wall
  and sample (equal to the wall fingerprint without inventory samples; the
  equilibrium list too, as exact scalars), the reservoirs list that of the
  sample elements (cells, sample names and cell centres).  They are
  restored only onto the same elements, not after a re-decomposition or
  renumberMesh, and read only when the file handler can read them for this
  run (perRankFile(): not a collated container of another number of ranks,
  as reconstructPar or decomposePar leave it, and only by the ranks of the
  container's IO group); otherwise mw_ and Cs_ are re-read (below).  The
  per-rank lists lie in the subdirectory uniform/phaseChange/, never in
  uniform/ itself: the master-reading file handlers (masterUncollated,
  collated) take an object whose local directory is exactly "uniform" for
  global data and resolve it to the master's file on every rank, so every
  restart on the same decomposition lost the exact state (review of M5,
  round 4).  decomposePar and reconstructPar copy (or link) the
  subdirectory with uniform/.

  mw_<condensate>  [kg/m2] volScalarField, output: its values on the
                   interface patches are the element deposits m''_e, per
                   EFFECTIVE element area A_e = |S_f|*areaMultiplier (times
                   the scale of areaPerVolume layer), not per wall area;
                   the internal field is the area mean of the elements of a
                   cell.  surfaceFieldValue areaIntegrate of mw_ on the
                   wall is therefore the wall inventory only for
                   areaMultiplier 1 and areaPerVolume cell (in general it is
                   sum_e m''_e |S_f|, the wall inventory being sum_e m''_e
                   A_e); the deposit per wall area, sum_e m''_e A_e/|S_f|,
                   is mDep_ (below).  mw_ is refreshed every step
                   (function objects may sample it), its processor patches
                   are evaluated, and it is written with the fields.  It is
                   read back only when uniform/phaseChange/
                   phaseChangeDeposits does not fit the elements (another
                   decomposition, a renumbered mesh), with a warning.  The
                   deposits are not kept in mw_ itself because OpenFOAM
                   writes a patch (or internal field) whose values are all
                   equal as 'uniform <value>' in text at writePrecision,
                   also in a binary file: a rank with a single WALL face
                   lost the deposit's digits.  layer distance (M4, review
                   item B3) keeps its elements in the same per-rank lists;
                   its mw_ is a cell field (the deposit m''_e of every
                   layer cell; 0 on the patches), re-read from the cells.
  mDep_<condensate> [kg/m2] volScalarField, output: the deposit per WALL
                   area.  On the interface patches, face f holds
                   sum_{e -> f} m''_e A_e/|S_f| (each element of layer
                   faceCells maps to its own face, so mDep_ = m''_e
                   A_e/|S_f|, computed as m''_e (A_e/|S_f|): exactly
                   2 m''_e for areaMultiplier 2; an element of layer
                   distance maps to its nearest wall face, on any rank);
                   areaIntegrate of mDep_ on the interface patches is the
                   wall inventory for any layer, areaMultiplier and
                   areaPerVolume.  The internal field is the deposit per
                   wall area of a cell, (sum_e m''_e A_e)/(sum_e |S_f|)
                   (layer distance: per the wall area V_c A_I/V_I that
                   the cell represents).  Refreshed every step, processor
                   patches evaluated, written with the fields, never read.
  exch_<gas>       [kg/(m3 s)] volScalarField, output (balance {
                   writeExchangeField yes; }, default no): the realised
                   exchange of the last step per cell volume, sum_{e in c}
                   q_e/V_c, positive = evaporation into the gas (the release
                   of the samples is not included).  Refreshed every step,
                   written with the fields, never read.
  Y_<condensate>   derived output (sum_e m''_e A_e)/(rho V), recomputed every
                   step before runTime.write(); never solved.

  Axial profiles at the writes of Time.  profileWriter (object
  phaseChangeProfiles, created with a profiles dictionary) is registered
  AUTO_WRITE like the lists above, but its writeObject() writes the axial
  profiles of the state just written to postProcessing/axialProfiles
  instead of a file of the time directory.  Every write of Time thus has
  its profiles, also Time::writeNow()/writeAndEnd() of a function object
  between two steps (e.g. runTimeControl), which the solver does not see
  (review of M3, round 1).  All ranks write, so the reduction of the
  profiles is collective.  The profiles of one state are written once;
  a write of Time after the per-step profiles of a 'times' entry of the
  same state writes them again with the reason 'write time; profile time
  ...', as a write time of controlDict would.

  Fresh start = no ledger and no condensate state (uniform/phaseChange/
  phaseChangeDeposits or mw_, uniform/phaseChange/phaseChangeReservoirs or
  Cs_) at the start time.  The state without the
  ledger stops the run (review, non-blocking 2), and so does the ledger
  with neither deposits of these elements nor all mw_ fields.  A fresh start
  after the release window of a sample has begun stops the run too: that
  part of the amount would never be released.

  Restart booking.  The inventory re-read (gas from rho*Y_<gas>, wall from
  the deposits) minus the inventory held by the ledger is booked as
  'restart'.  When the inventory re-read is stored exactly it must be
  <= 1e-14 relative (FatalError otherwise: the files do not belong to the
  same run and time); otherwise a warning is printed above
  10^(1 - writePrecision).  The inventory is exact when the deposits come
  from uniform/phaseChange/phaseChangeDeposits, or from mw_ written in
  binary with no rank whose interface patch holds one repeated non-zero
  value; Y_<gas> is written in binary with no rank whose internal field is
  a rounded 'uniform' entry; and, when the carrier is read from the restart
  directory (writeFrozenFields yes, where this solver wrote it), the same
  holds for rho.  OpenFOAM writes a field or patch whose values are all
  equal as 'uniform <value>' in text at writePrecision, also in a binary
  file; a uniform carrier computed with more digits is then 'rounded'
  (review round 4).

  The verdict also covers what the next steps re-read beyond the inventory
  (review, final M2 pass): Y_<gas> on the non-coupled patches (a fixedValue
  or inletOutlet value) and, with writeFrozenFields yes, the patches of rho
  and T, p and phi (format, internal field, patches; the processor patches
  of phi, which cannot be evaluated).  A patch whose value is recomputed
  when it is read is not concerned: zeroGradient, fixedGradient and mixed
  (from the cells and the patch parameters), symmetry, symmetryPlane and
  wedge (from the cells), empty, and coupled patches such as cyclic.  If
  one of them is a rounded 'uniform' entry, the inventory bound still
  applies, but the restart does not repeat the continuous run bit for bit:
  the verdict is 'inventory exact, continuation rounded' and the values are
  listed ('rounded for the next steps: ...').  The verdict covers the
  inventory, the internal fields and the patch VALUES only: other boundary
  parameters that OpenFOAM writes as 'uniform' text as well (inletValue of
  inletOutlet, refValue/refGradient/valueFraction of mixed, gradient of
  fixedGradient) are not checked, and the restart line says so
  ('checked: the inventory, the internal fields and the patch values
  re-read, not other boundary parameters').
  Whether a 'uniform' entry lost digits is decided by the run that writes
  it (the text of the value at writePrecision reads back to it or not;
  a value re-read always reproduces its own text): recordUniformValues()
  stores the fields checked and their rounded entries, merged over the
  ranks, in uniform/phaseChangeLayout ('uniformValues').  A restart uses
  that record when it belongs to the files (the element deposits of
  uniform/ fit the elements, i.e. not after decomposePar, reconstructPar or
  renumberMesh); otherwise every repeated non-zero value counts as rounded
  ('possibly rounded for the next steps').

  The restart line names every source and the verdict ('gas: ...;
  carrier: ...; deposits: ...[; rounded for the next steps: ...]; exact |
  inventory exact, continuation rounded | rounded'), followed by what the
  verdict covers.  The format is read
  from the field headers ('data.format' of a collated decomposedBlockData
  container, whose own 'format' is always binary; writtenBinary()).  The
  inventories are compensated sums (compensatedSum.H), so the bound holds
  for any summation order, e.g. after decomposePar or renumberMesh.  A
  bit-identical continuation needs the carrier of the continuous run as
  well: writeFrozenFields no (the carrier is then read from its original
  directory in both runs), or a carrier without rounded 'uniform' entries.

(continued in the next comment block: restart record, processor patches,
clock, solver dictionaries, wall resistance, schemes)

------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

interfaceExchange.H, continued.

  Restart record (uniform/phaseChangeLayout).  The deposits of a restart
  are per effective element area A_e on the faces of the interface
  patches, and the ledger has booked the release up to the restart under
  the recorded schedule.  The layout therefore records the interface
  (layer, patches, areaPerVolume, areaMultiplier; for layer distance also
  the distance, which is not recorded for faceCells) and every release sample
  (pair, amount, startTime, duration; scalars exact, as strings of
  phaseChangeLedger::exactDecimal()).  A restart stops with a FatalIOError
  that names the change when the interface differs (the wall inventory
  would change by the ratio of the areas, silently booked as rounding in
  ascii), when a sample whose window has begun is changed, when a sample
  is removed before its window has ended, and when a new sample's window
  has begun (as on a fresh start).  Samples whose windows lie after the
  restart may be added, changed or removed.  A layout of an earlier build
  without these records gives a warning.

  Processor patches of Y_<gas>.  They are evaluated from the neighbour
  cells at start-up: the file holds a single-face (or all-equal) processor
  patch as 'uniform <value>' at writePrecision, and these values enter the
  first matrix of a restart (the SuperBee limiter, the explicit gradients),
  so a binary restart would not repeat the continuous run.  At the end of
  every step of the continuous run they are the neighbour cell values.

  Clock.  The release schedule uses the exact step times kept in the
  ledger (phaseChangeLedger.H), never the rounded time-directory names; a
  restart continues the ledger's clock and warns if the start-time name
  differs from it by more than 1e-6*deltaT.

Solver dictionaries.  A separate Final solve exists when Y_<gas>Final
selects its own entry of the fvSolution solvers (solverEntry(), as
solution::solverDict(): the literal keyword, else the last matching
pattern).  A single pattern "Y_.*" serves both names, so there is no
separate Final solve and "Y_.*" needs relTol 0; "Y_.*Final" must follow
"Y_.*" to be used.  The start-up line names a pattern entry.

wallResistance halfCell.  g_e = A_e rho v' with v' = 1/(1/v + R),
v = v_d f(T_e) and R = rho d_e/(rhoD)_e the resistance of the half cell,
d_e = 1/deltaCoeffs of the wall face; rhoD = 0 on the face (no diffusion)
gives v' = 0.  The delta coefficients are built collectively and are
therefore obtained on every rank, before the element loop.  The series law
is a Robin condition per physical wall area |S_f|, so halfCell needs
areaMultiplier 1 and areaPerVolume cell (FatalIOError otherwise).  With
A_e = m |S_f| the element law is |S_f| rho/(1/(m v) + R/m): the scale
divides the half-cell resistance by m as well, whereas the Robin condition
with the kinetic area m |S_f| keeps the half cell of the mesh face,
|S_f| rho/(1/(m v) + R); the two agree only in the kinetic limit
R << 1/(m v).  The model line then gives the kinetic rate; the rate applied
with the half cell, sum_e A_e v'_e/V_I, is printed per pair.

Convection scheme.  The ledger books the transport as face fluxes, so a
paired species needs 'Gauss <interpolation>' for div(phi,Y_<gas>)
(checked at start-up and every step); 'bounded Gauss' subtracts
fvm::Sp(div(phi),Y), a source the ledger would not book.

Keys of the other sample mode are refused instead of ignored: a release
sample stops at removeTime, areaPerVolume, kineticScale or Ce, an inventory
sample at startTime or duration.

WALL mask.  The mask of the cells next to the patch WALL (and its line
'Cells adjacent to WALL') belongs to model mock only (rhoFixedFlowFoam.C);
the interface patches of the other models may have any name.

(continued in the next comment block: HKS and inventory samples)

------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

interfaceExchange.H, continued: model HKS and inventory samples (milestone
M5; plan E2, E3, E6, E7, E9 and section 20).

HKS coefficients (HKSCoeffs: accommodation sigma, Ce, kineticScale; T is
frozen, so they are evaluated once).  A wall element e in cell c has

  G_e    = kineticScale Ce 2 sigma/(2 - sigma) sqrt(M/(2 pi R T_e))  [s/m]
  beta_c = rho_c R T_c/M  [Pa],   Y_eq,e = p_eq(T_e)/beta_c,
  a_e    = A_e G_e p_eq(T_e),   g_e = A_e G_e beta_c,

and with wallResistance halfCell (faceCells, areaMultiplier 1, areaPerVolume
cell; Robin-equivalent) the conductance in series with the half cell,
G'_e = G_e/(1 + G_e beta_c d_e/(rhoD)_e) (0 for rhoD = 0 on the face).
kineticScale 1 is SI; 1e-5/sqrt(1000) = 3.16228e-7 reproduces the T-Flows
expression (bar, g/mol).  p_eq comes from the table of the pair
(vapourPressure { file; phases; units; interpolation; outOfRange; },
vapourPressureTable.H): the minimum over the listed phases.  The start-up
prints the table (nodes, range, the temperatures where the stable phase
changes), and per pair the range of T_e, p_eq and Y_eq, the layer mean of
lambda = sum_e g_e/(rho_c V_c) [1/s], the elements outside the table and
those with p_eq > p (the condensate would boil; reported, not capped).

Inventory samples (mode inventory; model HKS only).  Each cell of the
sample is a sample element of its pair (above) with the reservoir Cs_c
[kg/m3]:

  - fresh start: Cs = n0 M/V_Z in every cell of the sample (plan E7), part
    of the initial inventory of the ledger;
  - evaporation only: q_e = a_e - g_e Y within [0, Cs V/deltaT]; CAPPED sets
    Cs := 0 exactly (the round-off residue is booked as clamped, like a
    wall element), NONE when the gas is above Y_eq (no deposition);
  - Cs^{n+1} = Cs^n - deltaT q_e/V_c (plan E6);
  - removal at the start of the first step whose start t^n on the ledger's
    exact clock satisfies t^n >= removeTime - deltaT/2: the reservoir
    sum Cs V is booked as 'removed', Cs := 0, once; the ledger counts the
    removals of a pair (REMOVAL_DONE) and the layout records the sample as
    removed (with the clock of the removal), so a restart never removes it
    again.

State of the reservoirs: uniform/phaseChange/phaseChangeReservoirs, a
per-rank binary list (binaryIOList.H) like the deposits: a fingerprint of
the sample elements of the rank (their number and a hash of their cells,
the names of their samples and the cell centres), then Cs of every sample
element of every pair.  Cs_<condensate> [kg/m3] is an output volScalarField (the
reservoirs in the sample cells of the pair, 0 elsewhere), refreshed every
step, and read back only when the list does not fit the elements (another
decomposition, a renumbered mesh, inventory samples reordered or dropped),
with a warning, as mw_ for the deposits; the restart verdict names the
source.  When the list fits, a Cs_ that differs from it by more than its
rounding is reported (a warning, as for mw_) and rewritten from the list.
When every inventory sample has been removed, neither is needed
(the reservoirs are 0).  The list, the fields and the restart record exist
only when some inventory sample exists.  Y_<condensate> (derived) includes
Cs/rho in the sample cells, the balance line and file hold 'sample' = sum
Cs V, and the profiles bin it as the sample line density.

Restart record of the samples (uniform/phaseChangeLayout): every sample has
its mode; an inventory sample its pair, amount, number of cells, volume and
volume-weighted centre (compensated sums over the ranks, as exact
decimals) with the error bounds that the rounding of this run's mesh
points implies (volumeTolerance, centreTolerance, meshPoints;
sampleGeometryTolerance()), removeTime and whether (and at which clock) it
was removed.  A restart stops with a FatalIOError when an inventory sample
is new (it would never be filled: the fill belongs to a fresh start),
changed (pair, amount and number of cells exactly; volume and centre to
1e-12, since their last bits depend on the decomposition, plus the bounds
of the points of both runs: e.g. after decomposePar with writeFormat
ascii, or reconstructPar of such a run onto the binary mesh; its reservoir
belongs to those cells), changes its mode, or is dropped before its
removal; removeTime and the kinetics may change.  A sample dropped after
its removal keeps its record for the rest of the run (the ledger's count
of removals includes it).  A fresh start at or after removeTime - deltaT/2
of an inventory sample stops (it would be filled and removed in the first
step).

------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

interfaceExchange.H, continued: layer distance (milestone M4; plan sections
1, 2 'layer distance', 12 B3 and items 9-10, 16 item 1, and 25).

  interface
  {
      patches   (WALL);
      layer     distance;
      distance  1e-4;              // [m]
      centreDistance               // optional: a printed comparison only
      {
          axis (1 0 0);  origin (0 0 0);  radius 2.4e-3;
      }
  }

Elements.  One element per cell c whose centre lies within 'distance' of the
interface patches, d_c <= distance, with d_c the exact distance of the cell
centre from the nearest face of the interface patches, over the faces of
all patches and all ranks (face::nearestPoint, exact for planar faces),
and f(e) that nearest face, identified by its global index (globalIndex
over the interface faces of all ranks); it may lie on another rank.  Every
rank receives the interface faces that can lie within 2 x distance of one
of its cell centres (bounding boxes, PstreamBuffers) and searches them with
an octree; ties (a centre equidistant from several faces) go to the face
with the lexicographically smallest centre.  The set, the nearest faces
and A_I/V_I therefore depend neither on the order of the interface patches
nor on the decomposition (but for the last bits of the cell centres next
to processor faces).  The first M4 round used patchDataWave<wallPointData<
label>>, whose near-wall correction gives the wrong faces with the combined
wall patch, lets the last of several patches overwrite the distance of a
corner cell, and sees only the wall faces of the cell's own rank (review of
M4, round 1; plan section 26).  The elements are sorted by cell, one per
exchange cell.

  A_e   = V_c (A_I/V_I) areaMultiplier, A_I = sum |S_f| of the interface
          faces and V_I = sum V_c of the layer cells (all ranks): every
          cell of the layer has the mean area per volume of the layer (the
          paper's Eq. 7: A_I/V_I applies uniformly to the interaction
          region).  areaPerVolume is 'layer' by definition; 'cell' is
          refused (a cell of the second layer has no wall face of its own);
  T_e   = the cell T (interfaceTemperature cell), or the T of the nearest
          wall face f(e) (interfaceTemperature wall, fetched by the map
          below);
  v_d   = A V_I/A_I of this layer when A is given (temperature model): the
          rate per volume of every layer cell is A f(T_e) times
          areaMultiplier, the volumetric form of T-Flows.
  halfCell is refused: no half cell of a wall face belongs to a cell of the
  second layer.

The start-up line gives the number of layer cells, A_I, V_I, A_I/V_I, the
distances of the selected cells and the nearest distance outside the layer
within 2 x distance (the margin of the rule), the number of elements whose
nearest face lies on another rank, and, with centreDistance, the number of
cells of the centre-distance rule R - r <= distance about the given axis
(the T-Flows test with the analytic radius; plan section 12, item 10) and
the number of cells on which the two rules disagree.  This is a comparison
only: the layer is always the set of the exact distance.

Nearest-face map (plan section 12, B3).  A mapDistributeBase built from the
global face index of every element (compact numbering: the local interface
faces first, then the remote faces that the elements of this rank need)
moves values between the faces and the elements:

  pull  the values of the nearest faces to the elements (interfaceTemperature
        wall);
  push  the masses m''_e A_e of the elements to their faces, summed:
        mDep_ on face f = sum_{e -> f} m''_e A_e/|S_f|, so that sum_f mDep_
        |S_f| is the wall inventory to round-off in serial and in parallel,
        also for an element whose nearest face lies on another rank.  The
        sums of the remote slots go in PstreamBuffers to the rank of the
        face, which adds them to its own: its own elements first, then the
        other ranks in ascending order, every contribution once
        (interfaceExchange.C, pushToFaces; review of M4, round 3: the
        scheduled reverse distribute of round 2 added every contribution
        twice between two ranks that map elements to each other's faces,
        and a nonBlocking one adds them in the order the receives
        complete).

The exchange of the faces, the global numbering, the map and every
distribute are collective and called on every rank, also on a rank without
interface faces or layer cells.

State.  The deposit m''_e of every element is kept in uniform/phaseChange/
phaseChangeDeposits like the deposits of layer faceCells, behind a
fingerprint that hashes the rule and the distance first, then the cell and
the exact cell centre of every element (the faceCells fingerprint, and so
its restart format, is unchanged).  The layout records 'layer distance' and
the distance (exact decimal), so a restart with another rule or distance
stops with a FatalIOError that names the change; a restart after
decomposePar, reconstructPar or renumberMesh (the fingerprint does not fit)
re-reads the deposits from the cells of mw_, with a warning.  The cell
value of mw_ is m''_e itself (not m''_e A_e/A_e, which differs from it by
an ulp for many values), so mw_ equals the deposits of uniform/ bit for bit
and the fallback re-reads them exactly (review of M4, round 1).

Profiles (axialProfiles.H).  The deposit of a layer cell is binned by the
axial overlap of its own cell, like the gas of that cell, not by the extent
of its nearest wall face (plan section 25): bins + outside = inventory
exactly, local to every rank, and the line density is that of T-Flows,
which sums the cells of an axial slice.  On a mesh extruded along the axis
(the pipe) a layer cell and its nearest wall face have the same axial
extent, and the two choices agree.  The wall temperature and the wall area
per bin come from all faces of the interface patches, as for layer
faceCells, so that T_dep is the wall temperature at the onset for both
rules.

------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

interfaceExchange.H, continued: the GEMS3K backend of p_eq (milestone M9;
plan E9 and section 32; the backend itself is gemsEquilibrium.H).

  HKSCoeffs { equilibrium GEMS; speciation none; }   // or lagged
  GEMSCoeffs { system ...; mode frozen; ... }         // or local
  pairs { <gas> { ...; vapourPressure { ... } gems { gas; condensates;
                                                       elements; } } }

The HKS law of an element keeps its form; only p_eq (and, with speciation
lagged, the partial pressure) changes:

  a_e = A_e G'_e p_eq,e,   g_e = A_e G'_e chi_e beta_c,
  G'_e = G_e, or with wallResistance halfCell G_e/(1 + G_e chi_e beta_c
         d_e/(rhoD)_e) (hksCoefficients()),

with G_e and the half cell stored per element at start-up (T is frozen).
equilibrium table: p_eq,e = the table at T_e and chi_e = 1, with the
arithmetic of M5 (the results of every table case are unchanged to the
last bit).  The table of every pair stays required: it is the fallback and
the validation of GEMS3K.

  frozen   setCoefficients() takes p_eq,e from the GEMS table of the pair
           (gemsEquilibrium::frozenPEq(), built at start-up), or from the
           pair's table below minTemperature or outside the GEMS table; a
           start-up line per pair counts the elements of each source and
           gives the largest |dlog10 p_eq,e/p_eq,table| of the GEMS
           elements.  Nothing changes per step.
  local    setCoefficients() starts from the table; beginStep() updates
           p_eq,e and chi_e of every element from GEMS3K at the start of
           the steps of the schedule of the absolute time index n
           (gemsEquilibrium::updateDue(): (n - 1) mod updateInterval = 0,
           and a first step without a lagged state), before the correctors
           (updateEquilibrium(): the cell's composition from Y at the start
           of the step), and a_e, g_e with them.  A step thus uses one law
           throughout its correctors, the guard and realise(); the regimes,
           the bounds and the ledger are unaffected (the ledger books the
           realised flows of whatever law).  The step report has a GEMS
           line per pair after the exchange line (counts of calls and of
           the sources, iterations, deviation; gemsEquilibrium::
           stepReport()), end() the summary.  The lagged state of the last
           update is per-rank restart state (uniform/phaseChange/
           phaseChangeEquilibrium, first comment block): a restart onto the
           same elements restores it (restoreEquilibrium()) and repeats
           the continuous run bit for bit up to its next update; without
           it the first step updates.

The start-up refuses mode local with a carrier that is not a gas at its p
and T (max |rho - p W/(R T)|/rho = max |T - p W/(R rho)|/T > 0.01, the
carrier line; paper mode), unless GEMSCoeffs allowNonPhysicalCarrier yes,
and speciation lagged without equilibrium GEMS in mode local.

Output (GEMSCoeffs writeEquilibriumField yes): pEq_<gas> [Pa], the p_eq of
the law, gemsSource_<gas>, the source of every element (0 GEMS3K, 1 below
minTemperature, 2 below minMoleFraction, 3 outside the GEMS grid or
table, 4 rejected, 5 failed; -1 the table before the first update), and
with speciation lagged chi_<gas>: on the interface patches (layer
faceCells) the value of the element of each face; in the cells of the
elements the mean (the largest source) of their elements; elsewhere 0 (-1
for the source).  Refreshed in balance(), written with the fields, never
read.

------------------------------------------------------------------------------*/

#ifndef interfaceExchange_H
#define interfaceExchange_H

#include "chemicalFormula.H"

#include "fvMesh.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fluidThermo.H"
#include "IOdictionary.H"
#include "PtrList.H"
#include "autoPtr.H"
#include "binaryIOList.H"
#include "mapDistributeBase.H"
#include "phaseChangeLedger.H"
#include "axialProfiles.H"
#include "phaseChangeKinetics.H"
#include "vapourPressureTable.H"
#include "gemsEquilibrium.H"

namespace LESTO {

class gasRespeciation;

class interfaceExchange {

public:

  enum class modelType { mock, temperature, HKS };

private:

  /*--------------------------------------------------------------------------
  One gas/condensate pair: its species, the per-element coefficients and
  state, the source fields of its gas equation and the counters of the
  current step and of the run.
  --------------------------------------------------------------------------*/
  class pairData {

  public:

    Foam::label  gas;                  /* species index of the gas */
    Foam::label  condensed;            /* species index of the condensate */
    Foam::scalar molarMass;            /* [kg/mol] */
    chemicalFormula formula;
    Foam::scalar elementGas = 0, elementWall = 0, elementSample = 0;
    Foam::scalar accommodation, Ce, kineticScale;
    bool pairHKS = false;
    Foam::word   solverName;           /* Y_<gas> */
    bool         finalSolver;          /* Y_<gas>Final selects its own
                                          solvers entry */

    /* model HKS: the p_eq table of the pair (vapourPressure) */
    Foam::autoPtr<vapourPressureTable> table;

    /* model HKS, per element (the law is a_e = A_e G'_e p_eq,e and g_e =
       A_e G'_e chi_e beta_c; hksCoefficients()): the kinetic conductance
       G_e [s/m], the half cell d_e/(rhoD)_e [m2 s/kg] of wallResistance
       halfCell (-1: no diffusion to the wall), p_eq of the table at T_e
       and the p_eq of the law [Pa] (the table, or GEMS3K), the speciation
       factor chi_e (1 but with speciation lagged) and the source of p_eq
       (gemsEquilibrium::sourceType; -1 for the table of equilibrium
       table, or before the first update of mode local) */
    Foam::scalarField conductance;
    Foam::scalarField resistance;
    Foam::scalarField pEqTable;
    Foam::scalarField pEq;
    Foam::scalarField chi;
    Foam::labelList   source;

    /* equilibrium GEMS with writeEquilibriumField: pEq_<gas> [Pa],
       gemsSource_<gas> and, with speciation lagged, chi_<gas> (output) */
    Foam::autoPtr<Foam::volScalarField> equilibriumField;
    Foam::autoPtr<Foam::volScalarField> sourceField;
    Foam::autoPtr<Foam::volScalarField> chiField;

    /* per element (wall elements, then sample elements) */
    Foam::scalarField a;               /* [kg/s] */
    Foam::scalarField g;               /* [kg/s] per unit Y */
    Foam::scalarField lower;           /* [kg/s] */
    Foam::scalarField upper;           /* [kg/s] */
    Foam::scalarField flow;            /* realised flow of the step [kg/s] */
    Foam::labelList   regime;          /* exchangeRegime */

    /* per wall element: the deposit m''_e [kg/m2]; per sample element: the
       reservoir Cs [kg/m3] (0 for the samples of other pairs) */
    Foam::scalarField deposit;
    Foam::scalarField reservoir;
    bool              inventory;       /* the pair has inventory samples
                                          (all ranks) */

    /* per cell [kg/(m3 s)] */
    Foam::volScalarField::Internal Su;     /* exchange + release, >= 0 */
    Foam::volScalarField::Internal Sp;     /* >= 0 */
    Foam::volScalarField::Internal Srel;   /* release only */

    /* output of the deposits: mw_<condensate> (per effective element
       area), mDep_<condensate> (per wall area) and, optionally, the
       realised exchange exch_<gas> */
    Foam::autoPtr<Foam::volScalarField> depositField;
    Foam::autoPtr<Foam::volScalarField> wallDepositField;
    Foam::autoPtr<Foam::volScalarField> exchangeField;

    /* output of the reservoirs of inventory samples: Cs_<condensate>
       [kg/m3] (only when the pair has inventory samples) */
    Foam::autoPtr<Foam::volScalarField> reservoirField;

    /* this step */
    Foam::scalar      releasedStep;    /* [kg] */
    Foam::scalar      inflowStep;      /* convective inflow [kg], estimate */
    Foam::label       nChanged;        /* regime changes, latest predictor */
    Foam::label       nChangedSignificant;
    Foam::scalarField transport;       /* [kg/s] per ledger patch */
    Foam::scalar      defect;          /* [kg/s] */
    Foam::scalar      residualOF;      /* gSum(fvMatrix::residual()) [kg/s] */
    Foam::scalar      explicitRemainder;  /* sum_c E_c of the transport-only
                                             matrix [kg/s] (matrixDefect.H;
                                             reportMatrixResidual only) */
    Foam::label       nExplicitRemainder; /* cells with E_c != 0 */
    Foam::scalar      explicitSource;     /* sum_c |b_c - s_c| [kg/s]: the
                                             explicit source of the matrix
                                             (e.g. its face-flux
                                             correction, in source form) */
    Foam::scalar      transportIdentity;  /* D = sum_c r_c + dM/dt + sum_b
                                             (F_b + C_b) - sum_c E_c over
                                             all ranks, summed exactly
                                             [kg/s]: 0 in the flux form
                                             (matrixDefect.H) */
    Foam::scalar      identityTerms;      /* the magnitude of its terms
                                             [kg/s] */
    Foam::scalar      identityRelative;   /* |D|/identityTerms */
    Foam::scalar      clampedStep;     /* [kg] */
    Foam::scalar      projectedStep;   /* [kg] */
    Foam::scalar      evaporatedAboveWarningStep;  /* [kg] */
    Foam::scalar      changedMass;     /* |flow change| of the regime changes
                                          of the latest predictor x dt [kg] */
    Foam::label       nChangedOuter;   /* the same, of the corrector that */
    Foam::label       nChangedSignificantOuter;  /* decided the outer loop */
    Foam::scalar      changedMassOuter;
    Foam::label       nGuardPasses;
    Foam::label       nGuardSwitched;

    /* whole run */
    Foam::label       nUnconverged;
    Foam::label       nGuardPassesRun;
    Foam::label       nProjectedRun;
    Foam::label       nProjectedSignificantRun;
    Foam::label       nComplementarityRun;
    Foam::label       nComplementaritySignificantRun;
    Foam::label       nFinalChangesSignificantRun;
    Foam::label       nProjectedAboveRun;         /* steps whose total is */
    Foam::label       nComplementarityAboveRun;   /* above massTolerance() */

    pairData (
      const Foam::fvMesh& mesh,
      const Foam::word&   gasFieldName,
      const Foam::label   nElements
    );

    /* sizes the element lists: nElements in total, of which nWall wall
       elements (deposits) and nSample sample elements (reservoirs) */
    void resizeElements (
      const Foam::label nElements,
      const Foam::label nWall,
      const Foam::label nSample
    );
  };

  /*--------------------------------------------------------------------------
  One sample: mode release, or mode inventory (model HKS).
  --------------------------------------------------------------------------*/
  class sampleData {

  public:

    Foam::word       name;
    Foam::label      pair;             /* index in pairs_ */
    bool             inventory;        /* mode inventory (else release) */
    Foam::labelList  cells;            /* local cells */
    Foam::label      nCells;           /* all ranks */
    Foam::scalar     volume;           /* V_Z, all ranks [m3] */

    /* the restart record: V_Z and the volume-weighted centre [m] as
       compensated sums over all ranks (to round-off independent of the
       decomposition), and the bounds [m3], [m] of their error that the
       rounding of the mesh points of this run implies (0 when they are
       binary; sampleGeometryTolerance()), with the format of the points;
       not recorded: the largest differences [m3], [m] that a restart may
       accept for the rounding of both runs, below what a selection
       changed by one cell can cause */
    Foam::scalar     recordVolume;
    Foam::vector     centre;
    Foam::scalar     volumeTolerance;
    Foam::scalar     centreTolerance;
    Foam::string     meshPoints;
    Foam::scalar     volumeLimit;
    Foam::scalar     centreLimit;
    Foam::scalar     amount;           /* n0 [mol] */

    /* release */
    Foam::scalar     startTime;        /* [s] */
    Foam::scalar     duration;         /* t_ev [s] */
    Foam::scalar     length;           /* the window as represented,
                                          releaseWindowLength() [s] */

    /* inventory: elements in cell order; shared cells group elements
       from different pairs, so one sample need not be contiguous */
    Foam::label      firstElement;
    Foam::labelList  elements;         /* inventory elements, in cell order */
    Foam::scalar     areaPerVolume;    /* a_s [1/m] */
    Foam::scalar     kineticScale;
    Foam::scalar     Ce;
    bool             removable;        /* removeTime given */
    Foam::scalar     removeTime;       /* [s] */
    bool             removed;          /* the removal has happened */
    Foam::scalar     removalClock;     /* t^n of the removal [s] */
  };

  /* solver objects */
  const Foam::fvMesh&                        mesh_;
  const Foam::fluidThermo&                   thermo_;
  const Foam::volScalarField&                rho_;
  const Foam::surfaceScalarField&            phi_;
  Foam::PtrList<Foam::volScalarField>&       species_;
  const Foam::wordList&                      speciesNames_;
  const Foam::PtrList<Foam::volScalarField>& rhoD_;
  const Foam::label                          nCorr_;
  const Foam::scalar                         tolerance_;

  /* constant/thermochemistryProperties */
  Foam::IOdictionary  dict_;
  modelType           model_;

  /* pairs */
  Foam::labelList               pairOfSpecies_;     /* -1: not a paired gas */
  Foam::PtrList<pairData>       pairs_;
  gasRespeciation* respeciation_ = nullptr;
  Foam::wordList elementNames_;
  Foam::PtrList<Foam::OFstream> elementFiles_;
  void recordFormulas();
  void writeElementBalance();

  /* model HKS (HKSCoeffs) */
  Foam::scalar      accommodation_;     /* sigma */
  Foam::scalar      Ce_;
  Foam::scalar      kineticScale_;

  /* model HKS, equilibrium GEMS: the GEMS3K backend (GEMSCoeffs; sixth
     comment block) */
  Foam::autoPtr<gemsEquilibrium> gems_;

  /* interface layer */
  Foam::labelList   interfacePatches_;
  Foam::word        layer_;             /* faceCells | distance */
  Foam::scalar      layerDistance_;     /* [m], layer distance */
  Foam::word        areaPerVolume_;
  Foam::scalar      areaMultiplier_;
  Foam::word        interfaceTemperature_;
  Foam::word        wallResistance_;
  Foam::scalar      layerArea_;         /* A_I, all ranks, physical [m2] */
  Foam::scalar      layerVolume_;       /* V_I, all ranks [m3] */
  Foam::label       nLayerCells_;       /* all ranks */
  Foam::label       nElementsTotal_;    /* all ranks */

  /* elements: the nWall_ wall elements, sorted by cell, then the sample
     elements, one per cell of an inventory sample; exchange cell k owns
     the elements cellStart_[k] ... cellStart_[k+1]-1 (the first
     nWallCells_ exchange cells are the cells of the wall elements) */
  Foam::labelList   elementCell_;
  Foam::labelList   elementPatch_;      /* wall elements of layer faceCells
                                           (-1 for layer distance) */
  Foam::labelList   elementFace_;       /* wall elements: face index within
                                           its patch (-1 for distance) */
  Foam::scalarField elementArea_;       /* A_e [m2]; sample: a_s V_c */
  Foam::scalarField elementWallArea_;   /* wall elements: |S_f| [m2] (layer
                                           distance: V_c A_I/V_I) */
  Foam::scalarField elementT_;          /* T_e [K]; sample: the cell T */
  Foam::labelList   elementSample_;     /* sample elements: the sample */
  Foam::labelList   exchangeCells_;
  Foam::labelList   cellStart_;
  Foam::label       nWall_;             /* local wall elements */
  Foam::label       nWallCells_;        /* local exchange cells of the wall */
  Foam::label       nSampleElementsTotal_;  /* all ranks */

  /* layer distance: the interface faces of this rank (patch, face index
     within the patch, in the order of interfacePatches_; their local
     index is their slot in the face map), the slot of the nearest wall
     face of every wall element in the compact numbering of the map, and
     the map (interfaceExchange.H, fifth comment block) */
  Foam::labelList   interfaceFacePatch_;
  Foam::labelList   interfaceFaceFace_;
  Foam::labelList   elementFaceSlot_;
  Foam::autoPtr<Foam::mapDistributeBase> faceMap_;

  /* temperature model */
  Foam::scalar      depositionVelocity_;   /* v_d [m/s] */
  Foam::scalar      k_;                    /* [1/K] */
  Foam::scalar      Tdep_;                 /* [K] */

  /* samples */
  Foam::List<sampleData> samples_;
  Foam::boolList         sampleCell_;   /* cells of every sample (release
                                           and inventory), local */

  /* balance controls */
  Foam::label       nGuard_;
  Foam::scalar      guardTolerance_;
  Foam::scalar      moleFractionWarning_;
  bool              writeFile_;
  bool              reportMatrixResidual_;
  bool              writeExchangeField_;

  /*--------------------------------------------------------------------------
  Hook of the axial profiles into every write of Time (object
  phaseChangeProfiles, registered AUTO_WRITE): writeObject() writes the
  profiles of the state just written (writeProfilesOfState(true)) and no
  file of the time directory.
  --------------------------------------------------------------------------*/
  class profileWriter : public Foam::regIOobject {

    interfaceExchange& exchange_;

  public:

    profileWriter(const Foam::IOobject& io, interfaceExchange& exchange)
    :
      Foam::regIOobject(io),
      exchange_(exchange)
    {}

    /* the other (deprecated) overload stays visible */
    using Foam::regIOobject::writeObject;

    virtual bool writeData(Foam::Ostream&) const override {
      return true;
    }

    /* called by Time (objectRegistry::writeObject) for every write */
    virtual bool writeObject (
      Foam::IOstreamOption streamOpt,
      const bool           writeOnProc
    ) const override;
  };

  /* axial profiles (profiles sub-dictionary; absent: none) */
  Foam::autoPtr<axialProfiles> profiles_;
  Foam::autoPtr<profileWriter> profileWriter_;

  /* the state of the last profiles written: its clock, whether Time wrote
     it (reason 'write time'), and the profile times it covered; the
     profile times reached by the last step, not yet written */
  bool              profilesWritten_;
  Foam::scalar      profilesClock_;
  bool              profilesAtWrite_;
  Foam::string      profilesTimes_;
  Foam::string      profilesPending_;

  /* carrier molar concentration p/(R T) [mol/m3] per cell, of the frozen
     p and T, for the mole-fraction monitor */
  Foam::scalarField carrierConcentration_;

  /* the exact step [stepStart_, stepEnd_] on the clock of the ledger, and
     its deltaT */
  Foam::scalar      stepStart_;
  Foam::scalar      stepEnd_;
  Foam::scalar      stepDeltaT_;

  /* ledger and the per-rank element state of uniform/: the fingerprints
     of the wall elements (deposits), of all elements (regimes, and the
     lagged state of equilibrium GEMS, mode local) and of the sample
     elements (reservoirs) */
  Foam::autoPtr<phaseChangeLedger>          ledger_;
  Foam::labelList                           fingerprint_;
  Foam::labelList                           elementsFingerprint_;
  Foam::labelList                           samplesFingerprint_;
  Foam::autoPtr<binaryIOList<Foam::label>>  regimesIO_;
  Foam::autoPtr<binaryIOList<Foam::scalar>> depositsIO_;
  Foam::autoPtr<binaryIOList<Foam::scalar>> reservoirsIO_;
  Foam::autoPtr<binaryIOList<Foam::scalar>> equilibriumIO_;


  /* construction steps */
  void readPairs (
    const Foam::List<Foam::word>&   state,
    const Foam::List<Foam::scalar>& molarMass
  );
  void readInterface();

  /*--------------------------------------------------------------------------
  layer distance: the elements (the cells whose centre lies within
  layerDistance_ of the interface patches: the exact distance from the
  nearest interface face of any patch and any rank, searched with an
  octree over the faces within reach of the rank), their nearest wall
  faces and the face map; interfaceFaceArea returns |S_f| of the local
  interface faces.  Returns the start-up lines of the layer: the distances
  of the selected cells and the nearest outside, and, when d holds
  centreDistance, the comparison with the centre-distance rule.
  Collective.
  --------------------------------------------------------------------------*/
  Foam::string readDistanceLayer (
    const Foam::dictionary& d,
    Foam::scalarField&      interfaceFaceArea
  );

  /* layer distance: the values of the local interface faces (in the order
     of interfaceFacePatch_) at the nearest face of every wall element, on
     any rank.  Collective. */
  Foam::tmp<Foam::scalarField> nearestFaceValues (
    const Foam::scalarField& faceValues
  ) const;

  /* layer distance: the sum over the wall elements whose nearest face is a
     local interface face of their values, per local interface face (the
     reverse of nearestFaceValues(): its own elements first, then those of
     the other ranks in ascending order, each once).  Collective. */
  Foam::tmp<Foam::scalarField> nearestFaceSums (
    const Foam::scalarField& elementValues
  ) const;

  void readTemperatureModel();
  void readHKSModel();
  void readVapourPressure(pairData& pair, const Foam::dictionary& pd);
  void readSamples();
  void checkPairedSpecies(pairData& pair);
  void setCoefficients(pairData& pair) const;
  void readState();

  /* model HKS: a_e and g_e of element e from pair.pEq[e], pair.chi[e],
     its conductance and its half cell (halfCell: wallResistance halfCell;
     the same arithmetic as before M9, so equilibrium table is unchanged
     to the last bit); returns the conductance G'_e with the half cell */
  Foam::scalar hksCoefficients (
    pairData&         pair,
    const Foam::label e,
    const bool        halfCell
  ) const;

  /* equilibrium GEMS, mode local: p_eq (and chi) of every element at the
     composition of the start of the step, and a_e, g_e (beginStep()) */
  void updateEquilibrium();

  /* pEq_<gas>, gemsSource_<gas> and chi_<gas> from the elements of a pair
     (writeEquilibriumField) */
  void updateEquilibriumFields(pairData& pair) const;

  /*--------------------------------------------------------------------------
  Restart record in uniform/phaseChangeLayout (see "Restart record" above):
  recordInterface() stops a restart whose interface (patches,
  areaPerVolume, areaMultiplier) differs from the recorded one;
  recordSamples() checks the release windows on a fresh start, and on a
  restart the schedule of the release samples, the inventory samples
  (new, changed, dropped before their removal, a changed mode) and the
  removal count of the ledger, and restores the removal flags.  Both then
  record the current definitions; both run before the restart booking.
  --------------------------------------------------------------------------*/
  void recordInterface();
  void recordSamples();

  /*--------------------------------------------------------------------------
  The error bounds of the recorded V_Z and centre of an inventory sample
  that the rounding of the mesh points implies (review of M5, round 3),
  and the limits of what a restart accepts (round 4).  Binary points are
  exact: 0.  Points read from ascii carry p significant digits: at least
  the most any coordinate of the sample's cells needs (significantDigits()),
  and at least min(6, writePrecision), the fewest digits OpenFOAM's
  utilities write points with (6 by default; blockMesh, gmshToFoam and
  decomposePar at least 10).  The digits a value needs are no bound on
  their own: a round coordinate ('0.011') needs 2, and a mesh of round
  coordinates gave bounds larger than the domain, which accepted a sample
  moved by 20 mm (round 4).  Every point is then off by at most eta =
  sqrt(3)*0.5*10^(e + 1 - p), e the decimal exponent of the largest
  coordinate.  To first order V_Z changes by at most eta times the surface
  of the sample, and the centre by at most eta times that surface times
  the largest distance of a point from the centre, over V_Z; the surface is
  bounded by the sum of the face areas of the sample cells, and the bounds
  are doubled (the face centres of the decomposition, second order).  A
  restart compares with 1e-12 (the summation order) plus the bounds of
  both runs, but never with more than the limits: a quarter of h_min
  V_min/V_Z for the centre (h_min the thinnest sample cell, V/its largest
  face; a selection with one cell replaced by another of the same volume
  moves the centre by about a cell distance times V_min/V_Z), and half of
  V_min for the volume.  So the restarts of the ascii-decomposed meshes
  of the tests (10 and 12 digits) pass in either direction, while a
  selection moved, or changed by a cell of the same volume, is refused
  whatever digits the points have; a restart whose points are too coarse
  to tell their rounding from such a change is refused too, naming the
  limit.  Collective.
  --------------------------------------------------------------------------*/
  void sampleGeometryTolerance(sampleData& s, const bool exactPoints) const;

  /*--------------------------------------------------------------------------
  The final solver dictionary of a pair: Y_<gas>Final if it selects its
  own solvers entry (solverEntry), and its relTol must be 0 (FatalIOError
  otherwise).  Called at start-up and at every step, because fvSolution
  may be modified during the run (runTimeModifiable).
  --------------------------------------------------------------------------*/
  void checkFinalSolver(pairData& pair, const bool startup);

  /* the solvers entry of fvSolution that OpenFOAM selects for a field name
     (literal keyword, else the last matching pattern); nullptr if none */
  const Foam::entry* solverEntry(const Foam::word& fieldName) const;

  /*--------------------------------------------------------------------------
  The ddt scheme of a paired species must be Euler (FatalError otherwise).
  Called at start-up and at every step, because fvSchemes may be modified
  during the run (runTimeModifiable).
  --------------------------------------------------------------------------*/
  void checkDdtScheme(const pairData& pair, const bool startup) const;

  /*--------------------------------------------------------------------------
  The convection scheme div(phi,Y_<gas>) of a paired species must be
  conservative: 'Gauss <interpolation>', a sum of face fluxes (FatalIOError
  otherwise, e.g. for 'bounded Gauss ...', which subtracts
  fvm::Sp(fvc::div(phi),Y), a source the ledger does not book).  Called at
  start-up and at every step, like checkDdtScheme().
  --------------------------------------------------------------------------*/
  void checkConvectionScheme(const pairData& pair, const bool startup) const;

  /* Su, Sp of exchange cell k from the regimes of its elements */
  void cellCoefficients(pairData& pair, const Foam::label k) const;

  /* inventories [kg], compensated sums (compensatedSum.H) */
  Foam::scalar gasInventory(const pairData& pair) const;
  Foam::scalar wallInventory(const pairData& pair) const;
  Foam::scalar sampleInventory(const pairData& pair) const;

  /* the removal of the inventory samples due at the start of the step
     [t0, t0 + dt] (t0 >= removeTime - dt/2), booked in the ledger */
  void removeSamples(const Foam::scalar t0, const Foam::scalar dt);

  /* the reservoir mass [kg] of sample element e (index among all
     elements) */
  Foam::scalar reservoirMass(const pairData& pair, const Foam::label e) const;

  /*--------------------------------------------------------------------------
  Mass [kg] below which a bound violation, a regime change, a projection
  or a complementarity violation counts as round-off: guardTolerance times
  the inventory of the pair in this step (held at the end of the previous
  step, plus the release and the convective inflow through the patches of
  this step, so that a pair fed only through an inlet has a scale from its
  first step on).
  --------------------------------------------------------------------------*/
  Foam::scalar massTolerance(const Foam::label pairi) const;

  /*--------------------------------------------------------------------------
  Fingerprints of the elements of this rank: which = 0 the wall elements
  (the header of uniform/phaseChange/phaseChangeDeposits), 1 all elements
  (of uniform/phaseChange/phaseChangeRegimes; equal to the wall one without
  inventory samples, so the files of a run without them are unchanged), 2
  the sample elements (of uniform/phaseChange/phaseChangeReservoirs).
  --------------------------------------------------------------------------*/
  Foam::labelList elementsFingerprint(const int which) const;

  /* mw_<condensate> and mDep_<condensate> from the element deposits of a
     pair, and Cs_<condensate> from its reservoirs */
  void updateDepositField(pairData& pair) const;
  void updateReservoirField(pairData& pair) const;

  /* uniform/phaseChange/phaseChangeRegimes, phaseChangeDeposits,
     phaseChangeReservoirs and phaseChangeEquilibrium from the element
     regimes, deposits, reservoirs and lagged GEMS state, and
     recordUniformValues() (every step, before runTime.write()) */
  void storeElementState();

  /* the fields a restart re-reads whose 'uniform' entries will be written
     rounded (in uniform/phaseChangeLayout, "uniformValues", merged over
     the ranks; used by readState() of a restart from these files) */
  void recordUniformValues();

  /*--------------------------------------------------------------------------
  The axial profiles of the state of the ledger's clock (the end of the
  last step, or the start), with the profile times it reaches.  timeWrite:
  Time writes this state (profileWriter), reason 'write time', written
  unless done for this state already; otherwise (writeProfiles()) only
  when profile times are reached and this state has no profiles yet.
  Collective: every rank calls it.
  --------------------------------------------------------------------------*/
  void writeProfilesOfState(const bool timeWrite);

  /* the profile times written so far, in uniform/phaseChangeLayout
     (profiles { written (...); }), for a restart */
  void recordProfileTimes();

  /* true when the file of io (at its instance; a field of type typeName)
     was written in binary */
  static bool writtenBinary (
    const Foam::IOobject& io,
    const Foam::word&     typeName = Foam::volScalarField::typeName
  );

  /*--------------------------------------------------------------------------
  Whether this run can read a per-rank list of uniform/phaseChange/ (io,
  of the list type typeName) at the start time: ABSENT (no rank has it),
  READABLE (a plain list, or a collated container written for the ranks
  that read it: processors<N>/ with one block per rank of this run, or
  processors<N>_<lo>-<hi>/ of several IO ranks with one block per rank of
  the group lo..hi, which holds this rank; the fingerprint then decides
  whether it belongs to these elements), or FOREIGN (a collated container
  of another number of ranks, e.g. copied by reconstructPar into a serial
  case or by decomposePar onto another number of ranks, or a file that
  some ranks lack), which the file handler cannot read for this run.  why
  describes a FOREIGN file.  Collective.
  --------------------------------------------------------------------------*/
  enum class stateFile { ABSENT, READABLE, FOREIGN };
  static stateFile perRankFile (
    const Foam::IOobject& io,
    const Foam::word&     typeName,
    Foam::string&         why
  );

  /*--------------------------------------------------------------------------
  equilibrium GEMS, mode local, on a restart: the lagged state of the last
  update from uniform/phaseChange/phaseChangeEquilibrium (file and why as
  perRankFile() found it at the start time startName), restored when it
  belongs to these elements, with the law of every element recomputed from
  it; otherwise a warning when the first step would not update by the
  schedule (it then updates; gemsEquilibrium.H, mode local).  Collective.
  --------------------------------------------------------------------------*/
  void restoreEquilibrium (
    const stateFile         file,
    const Foam::string&     why,
    const Foam::word&       startName
  );

public:

  interfaceExchange (
    const Foam::fvMesh&                        mesh,
    const Foam::fluidThermo&                   thermo,
    const Foam::volScalarField&                rho,
    const Foam::surfaceScalarField&            phi,
    Foam::PtrList<Foam::volScalarField>&       species,
    const Foam::wordList&                      speciesNames,
    const Foam::List<Foam::word>&              state,
    const Foam::List<Foam::scalar>&            molarMass,
    const Foam::PtrList<Foam::volScalarField>& rhoD,
    const Foam::label                          nCorr,
    const Foam::scalar                         tolerance
  );

  interfaceExchange(const interfaceExchange&) = delete;
  void operator=(const interfaceExchange&) = delete;

  /* true for model mock (or no constant/thermochemistryProperties) */
  bool mock() const { return model_ == modelType::mock; }
  void attachRespeciation(gasRespeciation& r) { respeciation_ = &r; }
  void includeElement(const Foam::word& e) {
    if (elementNames_.find(e)<0) elementNames_.append(e);
  }
  phaseChangeLedger& ledger() { return ledger_(); }
  Foam::scalar wallMass(Foam::label pi) const { return wallInventory(pairs_[pi]); }
  Foam::scalar sampleMass(Foam::label pi) const { return sampleInventory(pairs_[pi]); }
  Foam::scalar cachedWallMass(Foam::label pi) const { return pairs_[pi].elementWall; }
  Foam::scalar cachedSampleMass(Foam::label pi) const { return pairs_[pi].elementSample; }
  Foam::scalar stepReleased(Foam::label pi) const { return pairs_[pi].releasedStep; }
  Foam::scalar stepDefect(Foam::label pi) const { return pairs_[pi].defect; }
  const Foam::scalarField& stepTransport(Foam::label pi) const { return pairs_[pi].transport; }
  Foam::scalar stepExchange(Foam::label pi) const;

  /* true for layer distance (the elements are cells, not wall faces) */
  bool distanceLayer() const { return layer_ == "distance"; }

  /* pair index of species speciesi, or -1 if it is not a paired gas */
  Foam::label pairOfGas(const Foam::label speciesi) const {
    return pairOfSpecies_[speciesi];
  }

  const Foam::volScalarField::Internal& Su(const Foam::label pairi) const {
    return pairs_[pairi].Su;
  }

  const Foam::volScalarField::Internal& Sp(const Foam::label pairi) const {
    return pairs_[pairi].Sp;
  }

  /* true when <field>Final selects its own solvers entry of fvSolution */
  bool finalSolver(const Foam::label pairi) const {
    return pairs_[pairi].finalSolver;
  }

  /* solver dictionary of the next solve */
  const Foam::dictionary& solverDict (
    const Foam::label pairi,
    const bool        final
  ) const;

  /* start of a physical time step */
  void beginStep();

  /*--------------------------------------------------------------------------
  Local implicit predictor at an outer corrector.  transportEqn is the gas
  equation WITHOUT exchange (ddt + div - laplacian).  Sets Su and Sp and
  counts the elements whose regime changed since the previous corrector
  (at the first corrector of a step: since the previous step), summed over
  all ranks, with the mass they carry and the significant ones.
  --------------------------------------------------------------------------*/
  void linearise (
    const Foam::label           pairi,
    const Foam::fvScalarMatrix& transportEqn
  );

  /*--------------------------------------------------------------------------
  Decision of the outer loop after a corrector: true when the initial
  residual is below the tolerance and the regime changes of the latest
  predictor are settled, i.e. the flow they change, times dt, is at most
  massTolerance() (the threshold at which the guard accepts a bound
  violation; see 'Settled regimes' above).  The regime changes of this
  deciding corrector are kept for the log.
  --------------------------------------------------------------------------*/
  bool outerConverged (
    const Foam::label  pairi,
    const Foam::scalar initialResidual
  );

  /*--------------------------------------------------------------------------
  After the final solve (and every guard re-solve): the transport per
  non-coupled patch of the solved matrix eqn, and the solver defect R' =
  sum_c [defect of the transport-only matrix transportEqn in flux form +
  the realised element flows + Srel V]_c at the solved Y (matrixDefect.H),
  which uses the flows q_e that update the condensate; the flux form
  telescopes the face fluxes and measures the storage as gasInventory()
  does, so the double-precision imbalances of the operator are booked in
  solverDefect, not in the closure (plan section 27, item 1)
  --------------------------------------------------------------------------*/
  void recordSolve (
    const Foam::label           pairi,
    const Foam::fvScalarMatrix& eqn,
    const Foam::fvScalarMatrix& transportEqn,
    const Foam::volScalarField& Y
  );

  /* true when elements were switched and the equation must be re-solved */
  bool guard (
    const Foam::label           pairi,
    const Foam::volScalarField& Y
  );

  /* realised exchange, projection, condensate update, monitors, log */
  void realise (
    const Foam::label     pairi,
    Foam::volScalarField& Y,
    const Foam::label     nOuter,
    const bool            converged,
    const Foam::scalar    initialResidual
  );

  /* ledger, closure, derived Y_<condensate>, mw_<condensate>,
     mDep_<condensate>, exch_<gas> and the element state of uniform/ of
     every pair (written by runTime.write()) */
  void balance();

  /*--------------------------------------------------------------------------
  Axial profiles and T_dep of every pair (axialProfiles.H) after
  runTime.write(), into postProcessing/axialProfiles only: at the profile
  times reached by this step, unless Time has written this state with its
  profiles already.  Nothing without a profiles dictionary.  Collective:
  every rank calls it at every step.
  --------------------------------------------------------------------------*/
  void writeProfiles();

  /* end-of-run summary */
  void end() const;
};

}

#endif
