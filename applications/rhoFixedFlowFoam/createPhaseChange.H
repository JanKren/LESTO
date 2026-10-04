/*------------------------------------------------------------------------------

createPhaseChange.H

PURPOSE

Creates the gas/condensate exchange object phaseChange once, before the
physical-time loop (LESTO::interfaceExchange, interfaceExchange.H).

Without constant/thermochemistryProperties, or with "model mock;", nothing
else is read, created or printed: the mock thermochemistry path of the
solver runs exactly as in the pristine build of tag 0.1 (the coupled mock
of the main branch): the logs of 008-011 and the fields of 010 (write
interval 3) are reproduced bit for bit.

Otherwise the constructor reads the dictionary, builds the interface
elements of the interaction layer and the samples (with one sample element
per cell of an inventory sample), reads the p_eq tables of model HKS
(vapourPressureTable.H) and, with equilibrium GEMS, creates the GEMS3K
engine of every rank and checks the system (mode frozen: the GEMS table of
every pair; gemsEquilibrium.H), runs the start-up checks of every paired gas
species (Euler ddt, a conservative Gauss convection scheme, zeroGradient on
the interface patches, relTol 0 of the final solve, relaxation ignored, an
outer tolerance below the solver tolerance, molar masses) and prints the
layer report (elements, A_I, V_I, A_I/V_I, v_d or the HKS ranges of T_e,
p_eq, Y_eq and lambda, sample cells and volume), the carrier consistency
and the fresh-start or restart inventory of every pair.  A fresh start
fills the reservoirs of the inventory samples; a fresh start after the
release window of a sample has begun stops the run; a restart continues
the clock of the ledger and the reservoirs.

------------------------------------------------------------------------------*/

LESTO::interfaceExchange phaseChange (
  mesh,
  thermo,
  rho,
  phi,
  species,
  speciesNames,
  state,
  molarMass,
  rhoD,
  nCorr,
  tolerance
);
