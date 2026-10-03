/*------------------------------------------------------------------------------

matrixDefect.H

PURPOSE

Defect [kg/s] of a solved scalar transport matrix per cell (cellDefect), in
the FLUX FORM of the conservation law it discretises, with the storage term
measured as the phase-change ledger measures the gas inventory:

  r_c = E_c - (m_c - m0_c)/deltaT - sum_{f of c} (F_f(psi) + C_f),

  E_c = b_c - fl(s_c + c_c),

with b the source of the matrix, s the source of its ddt term (the old-time
storage rho V psi^n/deltaT as the Euler scheme computes it), C the
face-flux correction of the matrix (fvMatrix::faceFluxCorrectionPtr(): the
explicit part of its face fluxes, e.g. the non-orthogonal correction of a
corrected laplacian of a field that is fluxRequired, as interfaceExchange
sets for every paired species), c = -V div(C) the contribution of C to the
source as fvMatrix computes it (gaussLaplacianScheme: source -= V
fvc::div(C)), fl(s + c) the double sum by which the matrix was assembled,
m0 and m the old and the new mass of the cell as the inventories compute
them (rho psi V, interfaceExchange::gasInventory()), and F_f the flux out
of c through face f as the off-diagonal coefficients give it.  Every face
flux F_f + C_f enters the two cells of an internal face once, with
opposite signs, so the sum over the cells telescopes to the boundary
transport exactly, whatever the rounding of the coefficients and of the
correction, and the storage term telescopes over the steps with the
inventory itself.  Besides round-off it is the only term by which the
discrete species inventory of a solved step differs from "storage =
sources - boundary transport":

  (M^{n+1} - M^n)/deltaT = sum_c E_c + sum_c V_c (Su - Sp*Y)_c - sum_p F_p - R,

where F_p is the transport out of the non-coupled boundary patch p (the
values of fvMatrix::flux(), correction included).  The phase-change ledger
books R*deltaT as "solverDefect", so its closure is at round-off for ANY
linear-solver tolerance and for any run length, and it changes by exactly
deltaT sum_c E_c per step, while |solverDefect| itself measures the
physical error of an incompletely converged solve plus the double-precision
imbalances of the discrete operator (below).

The explicit remainder E.  When the matrix is ddt + div - laplacian with an
interpolation scheme that has no explicit correction (upwind, the limited
schemes such as SuperBee: the solver's cases), its source is exactly fl(s
+ c), assembled in that order (v2412 and v2606: fvMatrix::operator-=,
gaussLaplacianScheme), and recomputing c with the same operations on the
same correction gives the same doubles, so E_c = 0 in every cell and the
closure does not change at all.  The source form b - s of round 2 left the
correction in r_c as V_c (sum_f C_f/V_c) plus the rounding of the assembly
s + c, whose sum over the cells telescopes only to round-off: on the
non-orthogonal pipe of run/014 (up to about 30 degrees) the closure
changed by about 1e-25 kg per step (review of M4, round 3; plan section
29, item 2).  Now that residue, the difference between the source form of
the correction and its face values, is part of R, i.e. of solverDefect.
An interpolation scheme with an explicit correction (linearUpwind, LUST,
...) adds it to the source only (gaussConvectionScheme::fvmDiv: fvm +=
fvc::surfaceIntegrate(phi*correction)), without face values in the matrix:
it stays in E, in source form, and the rounding of its cell sums reaches
the closure (a few eps of the correction per step).  E is printed with
balance { reportMatrixResidual yes; } (interfaceExchange).

The faces:
  - internal face f (owner l, neighbour u): F_f + C_f with F_f = upper_f
    psi_u - lower_f psi_l (lduMatrix::faceH), subtracted from l and added
    to u;
  - non-coupled boundary face: ic psi_c - bc + C (internalCoeffs,
    boundaryCoeffs), the value of fvMatrix::flux(), so the defect and the
    per-patch transport of the ledger use the same face values (C is 0
    there for the non-orthogonal correction: OpenFOAM's correction vectors
    vanish on non-coupled patches);
  - processor and cyclic face: each side's own flux ic psi_c - bc psi_nbr
    + C is exchanged (syncTools::swapBoundaryFaceList), and both sides use
    the shared value F = (F_own - F_other)/2, which is exactly
    antisymmetric.  The coefficients and the correction of the two sides
    of a processor face are computed separately (the interpolation weights
    w and w' of the two sides do not add up to 1 exactly, nor do the
    limiters or the interpolated gradients of the two sides agree to the
    last bit), so their own fluxes are not exactly opposite.  Coupled
    patches that are not face-matched (cyclicAMI and the like, not used by
    this solver) keep their own flux.

Every cell is accumulated in long double (64-bit mantissa on x86-64), so
that its rounding, which repeats identically in every step of a steady
state, stays far below the closure.

Why this form (review of M4, round 2; plan section 27, item 1).  The defect
was the row form b - A psi of the matrix, the residual of the LINEAR
SYSTEM, which cannot see two double-precision imbalances of the discrete
operator itself:
  - the column sums.  For a conservative operator the column sums of A
    (besides the storage) vanish in exact arithmetic, but the assembled
    diagonal of ddt + div - laplacian is the rounded negSumDiag of the
    convection and laplacian coefficients plus the rounded sums of the
    three matrices, so delta_c = diag_c - D_c + sum_{f: l(f)=c} lower_f +
    sum_{f: u(f)=c} upper_f (D the storage diagonal) is a few ulp of the
    diagonal per cell, and the solved system creates sum_c delta_c psi_c
    per unit time.  It drifted the closure of the T5_60 runs of M4 (dt 10
    ms: D dt/h^2 of the wall cells, h = 57 um, is about 300) by -1.2e-16
    of the inventory per step, to -7.6e-13 at 60 s;
  - the storage products.  The Euler scheme computes the new-time storage
    as ((rDeltaT rho) V) psi and the old-time one as ((rDeltaT rho) psi^n)
    V, and the ledger the inventory as (rho psi) V: three roundings of the
    same mass, the first fixed per cell, coherent where rho and V are
    uniform.  At a bitwise steady state (the gas of the M8-scale channel
    of M2.2h during its release) they repeat in every step, and the closure
    drifted by about 1e-17 of the gas per step.
The flux form with the inventory's masses books both in solverDefect, as
the stiff exchange's imbalance is (plan section 21, item 1), and leaves the
closure its random walk of round-off.

Which matrix.  The ledger does not take R from the ASSEMBLED matrix of the
paired gas (transport plus fvm::Sp plus Su).  In an exchange cell with
lambda deltaT of 1e3 to 6e3 (HKS) the diagonal and the source are that many
times the storage term, and their rounding when assembled is about eps
lambda deltaT of the gas of the cell, while the condensate changes by
deltaT q_e with q_e = a_e - g_e Y, a different rounding of the same stiff
quantity.  Near equilibrium both repeat every step, so the closure drifted
linearly with the number of steps (the closed box at 700 K: +6.7e-14 of the
inventory per step, 2.7e-12 after 50 steps; review of M5).  interfaceExchange
::recordSolve() therefore books

  R' = sum_c [ cellDefect(transport-only matrix, Y)_c
               + sum_{e in c} q_e(Y) + Srel_c V_c ],

the transport-only defect plus the realised flows q_e of exactly the values
that update the condensate, plus the release.  In exact arithmetic R' = R,
the defect of the assembled matrix; in double precision it closes the
ledger at the round-off of the inventory sums, and a real imbalance (of the
stiff exchange: a fixed point where the linear solver no longer iterates
but fl(a - g Y) is not 0; of the transport operator: above) is booked in
solverDefect, where it belongs.  (The rounding of the stored condensate is
booked as clamped, interfaceExchange::realise().)

fvMatrix::residual() is NOT used.  In v2412 it lets lduMatrix::residual()
add boundaryCoeffs*patchNeighbourField of the coupled patches through the
interface update and then adds the same term again (addBoundarySource with
couples = true), so its global sum is wrong in parallel (plan section 7).

Collective: cellDefect exchanges the fluxes of the coupled faces (and
evaluates the divergence of the correction, whose processor patches
communicate), so it must be called on every rank.

------------------------------------------------------------------------------*/

#ifndef matrixDefect_H
#define matrixDefect_H

#include "fvMatrices.H"
#include "fvcDiv.H"
#include "syncTools.H"
#include "processorFvPatch.H"
#include "cyclicFvPatch.H"

namespace LESTO {

  /* the defect r_c of every local cell [kg/s], in long double, and the
     explicit remainder E_c [kg/s] (header): eqn the transport-only
     matrix, ddtSource the source of its ddt term, oldMass and mass the
     inventory masses of the cells before and after the step [kg].
     Collective. */
  inline void cellDefect (
    const Foam::fvScalarMatrix&  eqn,
    const Foam::scalarField&     ddtSource,
    const Foam::scalarField&     oldMass,
    const Foam::scalarField&     mass,
    const Foam::scalar           deltaT,
    const Foam::volScalarField&  psi,
    Foam::List<long double>&     r,
    Foam::scalarField&           remainder
  ) {

    using namespace Foam;
    typedef long double ld;

    const fvMesh& mesh = psi.mesh();
    const scalarField& psiI = psi.primitiveField();
    const scalarField& b = eqn.source();
    const ld dt = static_cast<ld>(deltaT);

    /*------------------------------------------------------------------------
    The face-flux correction C and its source contribution c = -V div(C),
    with the operations of gaussLaplacianScheme (V times fvc::div), so that
    c holds the doubles the matrix was assembled with.  The matrix offers
    no const access to C; it is only read here.
    ------------------------------------------------------------------------*/
    const surfaceScalarField* correction = nullptr;
    scalarField correctionSource(psiI.size(), Zero);
    if (eqn.hasFaceFluxCorrection()) {
      correction =
        const_cast<fvScalarMatrix&>(eqn).faceFluxCorrectionPtr().get();
      const tmp<volScalarField> tDiv(fvc::div(*correction));
      const scalarField& divC = tDiv().primitiveField();
      const scalarField& V = mesh.V();
      forAll(correctionSource, celli) {
        correctionSource[celli] = -(V[celli]*divC[celli]);
      }
    }

    r.resize_nocopy(psiI.size());
    remainder.resize(psiI.size());
    forAll(r, celli) {
      const scalar assembled = ddtSource[celli] + correctionSource[celli];
      const ld E = static_cast<ld>(b[celli]) - static_cast<ld>(assembled);
      remainder[celli] = static_cast<scalar>(E);
      r[celli] = E - (static_cast<ld>(mass[celli])
                      - static_cast<ld>(oldMass[celli]))/dt;
    }

    /* internal faces: the flux out of the owner, once per face */
    const lduAddressing& addr = eqn.lduAddr();
    const labelUList& l = addr.lowerAddr();
    const labelUList& u = addr.upperAddr();

    if (eqn.hasLower() || eqn.hasUpper()) {
      const scalarField& lower = eqn.lower();
      const scalarField& upper = eqn.upper();
      forAll(l, facei) {
        const ld F = static_cast<ld> (
          upper[facei]*psiI[u[facei]] - lower[facei]*psiI[l[facei]]
        );
        r[l[facei]] -= F;
        r[u[facei]] += F;
      }
    }
    if (correction) {
      const scalarField& Ci = correction->primitiveField();
      forAll(l, facei) {
        const ld C = static_cast<ld>(Ci[facei]);
        r[l[facei]] -= C;
        r[u[facei]] += C;
      }
    }

    /* boundary faces: non-coupled ones directly (the value of flux()),
       coupled ones collect the flux of their own side (indexed by boundary
       face) */
    const label nInternal = mesh.nInternalFaces();
    scalarList ownFlux(mesh.nBoundaryFaces(), Zero);

    forAll(psi.boundaryField(), patchi) {

      const fvPatchScalarField& pf = psi.boundaryField()[patchi];
      const labelUList& faceCells = pf.patch().faceCells();
      const scalarField& ic = eqn.internalCoeffs()[patchi];
      const scalarField& bc = eqn.boundaryCoeffs()[patchi];
      const scalarField* Cp =
        correction ? &correction->boundaryField()[patchi] : nullptr;

      if (pf.coupled()) {
        const scalarField pnf(pf.patchNeighbourField());
        const label start = pf.patch().start() - nInternal;
        forAll(faceCells, facei) {
          ownFlux[start + facei] = ic[facei]*psiI[faceCells[facei]]
                                 - bc[facei]*pnf[facei];
          if (Cp) {
            ownFlux[start + facei] += (*Cp)[facei];
          }
        }
      } else {
        forAll(faceCells, facei) {
          r[faceCells[facei]] -= static_cast<ld> (
            ic[facei]*psiI[faceCells[facei]] - bc[facei]
          );
          if (Cp) {
            r[faceCells[facei]] -= static_cast<ld>((*Cp)[facei]);
          }
        }
      }
    }

    /* the flux of the other side of every processor and cyclic face */
    scalarList otherFlux(ownFlux);
    syncTools::swapBoundaryFaceList(mesh, otherFlux);

    forAll(psi.boundaryField(), patchi) {

      const fvPatchScalarField& pf = psi.boundaryField()[patchi];
      if (!pf.coupled()) {
        continue;
      }
      const fvPatch& patch = pf.patch();
      const bool matched =
        isA<processorFvPatch>(patch) || isA<cyclicFvPatch>(patch);
      const labelUList& faceCells = patch.faceCells();
      const label start = patch.start() - nInternal;

      forAll(faceCells, facei) {
        const scalar own = ownFlux[start + facei];
        r[faceCells[facei]] -= static_cast<ld> (
          matched ? 0.5*(own - otherFlux[start + facei]) : own
        );
      }
    }
  }
}

#endif
