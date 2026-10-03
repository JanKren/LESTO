/*---------------------------------------------------------------------------
  gemsbridge.h  --  plain-C ABI to the GEMS3K Gibbs-energy-minimisation kernel

  Purpose: let a code built with an old compiler (OpenFOAM v2412 + gcc 7.5,
  C++17) call GEMS3K 4.6.x (needs C++20 -> conda-forge gcc) without sharing a
  single C++ type across the boundary.  The shared library libgemsbridge.so
  contains GEMS3K, libstdc++ and libgcc statically and hidden; only the
  gemsb_* symbols below are exported.

  Units: T [K], P [Pa], amounts [mol], masses [kg].  Indices are 0-based and
  follow the GEMS3K DATABR ("xDB") ordering of the loaded system.

  Threading: an engine is NOT re-entrant.  Use one engine per MPI rank (or per
  thread).  Engines never communicate; there is no MPI inside the bridge.
---------------------------------------------------------------------------*/
#ifndef GEMSBRIDGE_H
#define GEMSBRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gemsb_engine gemsb_engine;         /* opaque handle */

/* return codes of gemsb_equilibrate() */
enum {
  GEMSB_OK            = 0,  /* converged (OK_GEM_AIA / OK_GEM_SIA)          */
  GEMSB_OK_RETRIED    = 1,  /* warm start failed, cold (AIA) retry is OK    */
  GEMSB_BAD_QUALITY   = 2,  /* BAD_GEM_*: result returned but not trusted   */
  GEMSB_ERR_NOCONV    = 3,  /* ERR_GEM_*: no result                          */
  GEMSB_ERR_FATAL     = 4,  /* T_ERROR_GEM: engine must be recreated         */
  GEMSB_ERR_INPUT     = 5   /* T/P outside lookup grid, bad b, NULL args    */
};

/* initial-approximation modes */
enum {
  GEMSB_COLD      = 0,  /* AIA: simplex/LPP initial approximation            */
  GEMSB_WARM_LAST = 1,  /* SIA from the engine's previous call (any cell)    */
  GEMSB_WARM_CELL = 2   /* SIA from caller-stored per-cell state             */
};

/* ---- global settings (call before the first create) -------------------- */
/* Directory for ipmlog.txt ("" = cwd).  Use one directory per MPI rank.     */
void gemsb_set_log_directory(const char* dir);
/* spdlog level for all GEMS3K loggers: 0 trace .. 3 warn, 4 error, 6 off.
   GEMS3K logs IPM warnings (with matrix dumps) to stdout by default; that
   costs time and floods OpenFOAM logs.  Use 4 in production.              */
void gemsb_set_log_level(int level);

/* ---- lifecycle ---------------------------------------------------------- */
/* dat_lst: GEM-Selektor export "<name>-dat.lst" (key-value or JSON files)  */
gemsb_engine* gemsb_create_from_lst(const char* dat_lst, char* errmsg, int errlen);
/* the same three documents as strings (xGEMS initializeFromJsonStrings)     */
gemsb_engine* gemsb_create_from_strings(const char* dch, const char* ipm,
                                        const char* dbr, char* errmsg, int errlen);
void gemsb_destroy(gemsb_engine* e);

/* ---- system description -------------------------------------------------- */
int         gemsb_num_elements(const gemsb_engine* e);
int         gemsb_num_species (const gemsb_engine* e);
int         gemsb_num_phases  (const gemsb_engine* e);
const char* gemsb_element_name(const gemsb_engine* e, int i);
const char* gemsb_species_name(const gemsb_engine* e, int j);
const char* gemsb_phase_name  (const gemsb_engine* e, int k);
int         gemsb_element_index(const gemsb_engine* e, const char* name); /* -1 */
int         gemsb_species_index(const gemsb_engine* e, const char* name); /* -1 */
int         gemsb_phase_index  (const gemsb_engine* e, const char* name); /* -1 */
int         gemsb_species_phase(const gemsb_engine* e, int j);
int         gemsb_species_is_gas(const gemsb_engine* e, int j);
double      gemsb_species_molar_mass(const gemsb_engine* e, int j); /* kg/mol */
double      gemsb_stoich(const gemsb_engine* e, int j, int i);  /* a_ji     */
void        gemsb_TP_range(const gemsb_engine* e, double* Tmin, double* Tmax,
                           double* Pmin, double* Pmax);
/* number of doubles the caller must store per cell for GEMSB_WARM_CELL      */
int         gemsb_state_size(const gemsb_engine* e);

/* Metastability bounds (GEMS3K dll/dul) on species j, persistent across calls.
   Units: mol per mol of total elements (the internal normalised basis).
   lower = upper = 0 suppresses a condensed phase: the call then returns the
   homogeneous gas speciation, and gemsb_log10_activities() gives the
   log10 saturation ratio of every suppressed pure condensate, so that
   p_sat,i = p_i / 10^log10(Omega) (verified in satind_test.cpp).          */
void gemsb_set_species_bounds(gemsb_engine* e, int j, double lower, double upper);

/* ---- the call ------------------------------------------------------------ */
/* b_mol[num_elements]: element amounts of the cell (any scale; internally
   normalised to 1 mol total, results are scaled back).
   state[state_size]: per-cell warm-start state, in/out, may be NULL.
   *state_valid (in/out): 0 on first use; set to 1 by a successful call.    */
int gemsb_equilibrate(gemsb_engine* e, double T_K, double P_Pa,
                      const double* b_mol, int mode,
                      double* state, int* state_valid);

/* ---- results of the last call that returned GEMSB_OK/OK_RETRIED/BAD ----- */
void   gemsb_species_amounts(const gemsb_engine* e, double* n_mol);   /* [nDC] */
void   gemsb_gas_partial_pressures(const gemsb_engine* e, double* p_Pa); /* [nDC], 0 if not gas */
void   gemsb_log10_activities(const gemsb_engine* e, double* log10a); /* [nDC]
          dual-thermodynamic: gas -> log10(f/1 bar); pure condensed phase ->
          log10 saturation ratio (0 when stable, <0 when undersaturated)  */
void   gemsb_phase_amounts(const gemsb_engine* e, double* n_mol);     /* [nPH] */
void   gemsb_phase_log10_saturation(const gemsb_engine* e, double* lg); /* [nPH] */
int    gemsb_last_iterations(const gemsb_engine* e);
double gemsb_last_seconds(const gemsb_engine* e);
const char* gemsb_last_error(const gemsb_engine* e);

#ifdef __cplusplus
}
#endif
#endif /* GEMSBRIDGE_H */
