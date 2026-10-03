/*------------------------------------------------------------------------------

solvePairedGasSpecies.H

PURPOSE

Solves one physical time step of a gaseous species that exchanges mass with
a condensate (a paired species of constant/thermochemistryProperties).  It
is included from the species loop instead of solveGasSpecies.H, which stays
unchanged for every other species, so their logs and results are identical
to a build without phase change.

For the paired species i the equation is

  ddt(rho,Y_i) + div(phi,Y_i) - laplacian(rhoD_i,Y_i) = Su - Sp*Y_i,

with the linearised exchange-plus-release source Su >= 0, Sp >= 0 of
phaseChange (interfaceExchange.H).  Compared with solveGasSpecies.H:

  - at every outer corrector the local implicit predictor chooses the
    regime of every interface element from the transport-only matrix
    (phaseChange.linearise()) BEFORE Su and Sp are added, starting from
    the regimes of the previous corrector (no element jumps between
    CAPPED and NONE; phaseChangeKinetics.H, acrossBounds);
  - the matrix is NEVER relaxed (a relaxation factor in fvSolution is
    ignored with a start-up warning): the realised exchange must be that of
    the solved equation;
  - the outer loop is converged only when the initial residual is below
    tolerance AND the regimes are settled (phaseChange.outerConverged()):
    the elements that changed their regime at this corrector carry at most
    guardTolerance*inventory (the threshold of the guard; changes of sign
    of Y at round-off level carry no mass and would otherwise never let the
    loop converge; plan section 14).  After convergence, or after nCorr+1
    correctors, one more solve with the <field>Final solver dictionary
    follows if it exists (tight, relTol 0), with its own predictor;
    otherwise the last solve was the final one (and <field> must have
    relTol 0, checked at start-up and at every step);
  - after the final solve, the transport through every non-coupled patch
    and the solver defect are recorded, and the guard re-solves (with the
    final dictionary, without a new predictor) as long as the realised
    exchange violates the bounded law by more than guardTolerance*
    inventory in total (an IMPLICIT element outside its bounds, or a
    CAPPED or NONE element whose admissible flow at the solved Y differs
    from its flow; both are switched, the latter back to IMPLICIT with the
    largest misses first, until those left are within the tolerance), at
    most nGuard times.  The
    defect is that of the transport-only matrix (a copy kept before Su and
    Sp are added) plus the realised element flows and the release
    (matrixDefect.H: the rounding of the assembled stiff diagonal would
    otherwise enter the closure);
  - a guard re-solve solves the SAME transport matrix as the solve it
    corrects (the transport-only copy of the recorded solve), with only Su
    and Sp changed by the guard's switches.  It never re-assembles ddt +
    div - laplacian: the SuperBee limiter (and the explicit non-orthogonal
    correction) would be re-evaluated at the solved Y, so a switch decided
    at the solution of one matrix would be applied to another.  Near a
    bound a single element then cycled between two regimes until nGuard
    was used up (a bare HKS wall element at Y_eq on a moving gas front:
    IMPLICIT, NONE, IMPLICIT, NONE), and the step was accepted with a
    bounded-law violation far above the tolerance (review of M5, round 4).
    For a fixed linear matrix, switching a violating element to its bound
    leaves the implicit flow on the same side of it, so the switches of
    the guard are monotone.  Transport convergence comes from nCorr and
    the outer tolerance, never from the guard;
  - phaseChange.realise() then books the realised exchange, updates the
    condensate and prints the exchange, transport and monitor lines.  The
    outer-convergence status is part of that line; paired species print no
    per-step warning.

------------------------------------------------------------------------------*/

volScalarField& Y = species[speciesi];
const volScalarField& rhoDi = rhoD[speciesi];
const word& fieldName = Y.name();
const word divScheme("div(phi," + fieldName + ")");
const label pairi = phaseChange.pairOfGas(speciesi);

scalar initialResidual = GREAT;
label nOuter = 0;           /* outer correctors solved */
bool converged = false;
bool finalSolve = false;    /* the next solve is the final (or a guard) solve */
bool guardResolve = false;  /* the next solve re-solves after a guard switch */

/* the transport-only matrix of the recorded solve, which a guard re-solve
   solves again with the switched Su and Sp */
autoPtr<fvScalarMatrix> recordedTransport;

while (true) {

  tmp<fvScalarMatrix> tYEqn;
  if (guardResolve) {
    tYEqn = tmp<fvScalarMatrix>::New(recordedTransport());
  } else {
    tYEqn =
      fvm::ddt(rho, Y)
        + fvm::div(phi, Y, divScheme)
        - fvm::laplacian(rhoDi, Y);
  }
  fvScalarMatrix& YEqn = tYEqn.ref();

  /*---------------------------------------------------------------------------
  Regimes from the transport-only matrix, except in a guard re-solve, which
  keeps the coefficients of the guard (and the matrix of the recorded
  solve, above).
  ---------------------------------------------------------------------------*/
  if (!guardResolve) {
    phaseChange.linearise(pairi, YEqn);
  }

  /*---------------------------------------------------------------------------
  The transport-only matrix of a solve whose defect may be recorded: the
  final solve and the guard re-solves, and without a <field>Final
  dictionary every corrector (any of them may turn out to be the last).
  ---------------------------------------------------------------------------*/
  autoPtr<fvScalarMatrix> transportEqn;
  if (finalSolve || !phaseChange.finalSolver(pairi)) {
    transportEqn.reset(new fvScalarMatrix(YEqn));
  }

  YEqn += fvm::Sp(phaseChange.Sp(pairi), Y);
  YEqn -= phaseChange.Su(pairi);

  const scalar residual =
    YEqn.solve(phaseChange.solverDict(pairi, finalSolve)).initialResidual();

  if (!finalSolve) {

    ++nOuter;
    initialResidual = residual;
    converged = phaseChange.outerConverged(pairi, initialResidual);

    if (!converged && nOuter <= nCorr) {
      continue;
    }

    finalSolve = true;
    if (phaseChange.finalSolver(pairi)) {
      continue;
    }
  }

  phaseChange.recordSolve(pairi, YEqn, transportEqn(), Y);
  recordedTransport = std::move(transportEqn);

  guardResolve = phaseChange.guard(pairi, Y);
  if (!guardResolve) {
    break;
  }
}

phaseChange.realise(pairi, Y, nOuter, converged, initialResidual);

Info<< fieldName << " min/max = " << gMin(Y.primitiveField())
  << " " << gMax(Y.primitiveField()) << nl;
