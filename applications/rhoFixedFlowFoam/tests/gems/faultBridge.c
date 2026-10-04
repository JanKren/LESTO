/*------------------------------------------------------------------------------

tests/gems/faultBridge.c

PURPOSE

A fault-injecting GEMS3K bridge for the criteria M9.5, M9.8 and M9.13 of
tests/gems/Allrun (plan section 12, item 7: after a failure the engine is
re-created, the species bounds are re-applied and the warm states
invalidated; plan section 35, items A, E, F, V and Y).  It exports every
function of thermochemistry/gemsbridge/gemsbridge.h and forwards each call
to the real libgemsbridge.so, which it opens with dlopen; it

  - returns GEMSB_ERR_FATAL instead of calling the real gemsb_equilibrate()
    at the calls numbered in LESTO_TEST_FATAL_CALLS (e.g. "300 900",
    counted from 1 over the process; at most 64), as GEMS3K does for a
    T_ERROR_GEM;
  - at the calls numbered in LESTO_TEST_ZERO_PGAS_CALLS lets the real call
    run (it converges and balances) and then reports a partial pressure 0
    for every species (gemsb_gas_partial_pressures()), and at those of
    LESTO_TEST_NAN_OMEGA_CALLS a log10 activity NaN for every condensed
    species (gemsb_log10_activities()): results the solver must not use;
  - reports LESTO_TEST_ITERATIONS (e.g. 1000000000) as the IPM iterations
    of every call (gemsb_last_iterations()), so that the totals of a short
    run exceed a 32-bit integer;
  - multiplies every partial pressure it reports by LESTO_TEST_PGAS_SCALE
    (e.g. 1.0000000000000004, 2 ulp), a stand-in for a bridge built with
    another compiler, whose results differ from these in their last bits
    (M9.8);
  - logs to LESTO_TEST_BRIDGE_LOG every engine created and destroyed, every
    bound set (species, lower, upper), and every gemsb_equilibrate() with
    its number, the per-element state array it got (an address, the same
    for the same element), whether that state was valid on entry and the
    status of the real call, and the faults injected.  An engine is named
    by its creation number (1, 2, ...), not by its address: a re-created
    engine may get the address of the one destroyed before it.

The solver is not rebuilt for it: its run path is a RUNPATH, so a run with
LD_LIBRARY_PATH pointing at the directory of this library loads it instead
of the real one, which it opens by its absolute path from
LESTO_TEST_REAL_BRIDGE.  Serial runs only (one log per process).

Built by tests/gems/Allrun with the C compiler of OpenFOAM: -shared -fPIC
-ldl, next to a copy of the real gemsbridge.h.

------------------------------------------------------------------------------*/

#include "gemsbridge.h"

#include <dlfcn.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void* real = NULL;
static FILE* logFile = NULL;
static long  nEquilibrate = 0;
static long  fatalAt[64];
static int   nFatal = 0;
static long  zeroPgasAt[64];
static int   nZeroPgas = 0;
static long  nanOmegaAt[64];
static int   nNanOmega = 0;
static int   iterations = -1;    /* LESTO_TEST_ITERATIONS, -1: the real */
static double pgasScale = 1.0;   /* LESTO_TEST_PGAS_SCALE */
static int   zeroPgasNow = 0;    /* the last call is one of them */
static int   nanOmegaNow = 0;
static const gemsb_engine* engines[64];   /* by creation number - 1 */
static int   nEngines = 0;

/* the creation number of an engine (the last one created at its address;
   0 for one not created through this bridge) */
static int engineNumber(const gemsb_engine* e) {
  int i;
  for (i = nEngines - 1; i >= 0; --i) {
    if (engines[i] == e) {
      return i + 1;
    }
  }
  return 0;
}

/* the call numbers of the environment variable name into at (at most 64) */
static int callNumbers(const char* name, long* at) {
  const char* calls = getenv(name);
  int n = 0;
  while (calls && *calls && n < 64) {
    char* end = NULL;
    const long k = strtol(calls, &end, 10);
    if (end == calls) {
      break;
    }
    at[n++] = k;
    calls = end;
  }
  return n;
}

static int listed(const long* at, const int n, const long call) {
  int i;
  for (i = 0; i < n; ++i) {
    if (at[i] == call) {
      return 1;
    }
  }
  return 0;
}

static void openReal(void) {
  const char* path;
  if (real) {
    return;
  }
  path = getenv("LESTO_TEST_REAL_BRIDGE");
  real = path ? dlopen(path, RTLD_NOW | RTLD_LOCAL) : NULL;
  if (!real) {
    fprintf(stderr, "faultBridge: cannot open LESTO_TEST_REAL_BRIDGE=%s: %s\n",
            path ? path : "(unset)", dlerror());
    abort();
  }
  path = getenv("LESTO_TEST_BRIDGE_LOG");
  if (path) {
    logFile = fopen(path, "a");
  }
  nFatal = callNumbers("LESTO_TEST_FATAL_CALLS", fatalAt);
  nZeroPgas = callNumbers("LESTO_TEST_ZERO_PGAS_CALLS", zeroPgasAt);
  nNanOmega = callNumbers("LESTO_TEST_NAN_OMEGA_CALLS", nanOmegaAt);
  path = getenv("LESTO_TEST_ITERATIONS");
  if (path && *path) {
    iterations = atoi(path);
  }
  path = getenv("LESTO_TEST_PGAS_SCALE");
  if (path && *path) {
    pgasScale = strtod(path, NULL);
  }
}

static void* symbol(const char* name) {
  void* f;
  openReal();
  f = dlsym(real, name);
  if (!f) {
    fprintf(stderr, "faultBridge: no %s in the real bridge\n", name);
    abort();
  }
  return f;
}

static void note(const char* text) {
  if (logFile) {
    fprintf(logFile, "%s\n", text);
    fflush(logFile);
  }
}

/* the real function NAME of type TYPE, looked up once */
#define REAL(TYPE, NAME) \
  static TYPE NAME##_real = NULL; \
  if (!NAME##_real) { *(void**)(&NAME##_real) = symbol(#NAME); }

typedef void (*setDir_t)(const char*);
void gemsb_set_log_directory(const char* dir) {
  REAL(setDir_t, gemsb_set_log_directory)
  gemsb_set_log_directory_real(dir);
}

typedef void (*setLevel_t)(int);
void gemsb_set_log_level(int level) {
  REAL(setLevel_t, gemsb_set_log_level)
  gemsb_set_log_level_real(level);
}

typedef gemsb_engine* (*createLst_t)(const char*, char*, int);
gemsb_engine* gemsb_create_from_lst(const char* dat_lst, char* errmsg,
                                    int errlen) {
  gemsb_engine* e;
  char text[64];
  REAL(createLst_t, gemsb_create_from_lst)
  e = gemsb_create_from_lst_real(dat_lst, errmsg, errlen);
  if (e && nEngines < 64) {
    engines[nEngines++] = e;
  }
  snprintf(text, sizeof(text), "create %d", engineNumber(e));
  note(text);
  return e;
}

typedef gemsb_engine* (*createStrings_t)(const char*, const char*,
                                         const char*, char*, int);
gemsb_engine* gemsb_create_from_strings(const char* dch, const char* ipm,
                                        const char* dbr, char* errmsg,
                                        int errlen) {
  REAL(createStrings_t, gemsb_create_from_strings)
  note("create from strings");
  return gemsb_create_from_strings_real(dch, ipm, dbr, errmsg, errlen);
}

typedef void (*destroy_t)(gemsb_engine*);
void gemsb_destroy(gemsb_engine* e) {
  char text[64];
  REAL(destroy_t, gemsb_destroy)
  snprintf(text, sizeof(text), "destroy %d", engineNumber(e));
  note(text);
  gemsb_destroy_real(e);
}

typedef int (*count_t)(const gemsb_engine*);
int gemsb_num_elements(const gemsb_engine* e) {
  REAL(count_t, gemsb_num_elements)
  return gemsb_num_elements_real(e);
}
int gemsb_num_species(const gemsb_engine* e) {
  REAL(count_t, gemsb_num_species)
  return gemsb_num_species_real(e);
}
int gemsb_num_phases(const gemsb_engine* e) {
  REAL(count_t, gemsb_num_phases)
  return gemsb_num_phases_real(e);
}
int gemsb_state_size(const gemsb_engine* e) {
  REAL(count_t, gemsb_state_size)
  return gemsb_state_size_real(e);
}
int gemsb_last_iterations(const gemsb_engine* e) {
  REAL(count_t, gemsb_last_iterations)
  return iterations >= 0 ? iterations : gemsb_last_iterations_real(e);
}

typedef const char* (*name_t)(const gemsb_engine*, int);
const char* gemsb_element_name(const gemsb_engine* e, int i) {
  REAL(name_t, gemsb_element_name)
  return gemsb_element_name_real(e, i);
}
const char* gemsb_species_name(const gemsb_engine* e, int j) {
  REAL(name_t, gemsb_species_name)
  return gemsb_species_name_real(e, j);
}
const char* gemsb_phase_name(const gemsb_engine* e, int k) {
  REAL(name_t, gemsb_phase_name)
  return gemsb_phase_name_real(e, k);
}

typedef int (*index_t)(const gemsb_engine*, const char*);
int gemsb_element_index(const gemsb_engine* e, const char* name) {
  REAL(index_t, gemsb_element_index)
  return gemsb_element_index_real(e, name);
}
int gemsb_species_index(const gemsb_engine* e, const char* name) {
  REAL(index_t, gemsb_species_index)
  return gemsb_species_index_real(e, name);
}
int gemsb_phase_index(const gemsb_engine* e, const char* name) {
  REAL(index_t, gemsb_phase_index)
  return gemsb_phase_index_real(e, name);
}

typedef int (*speciesInt_t)(const gemsb_engine*, int);
int gemsb_species_phase(const gemsb_engine* e, int j) {
  REAL(speciesInt_t, gemsb_species_phase)
  return gemsb_species_phase_real(e, j);
}
int gemsb_species_is_gas(const gemsb_engine* e, int j) {
  REAL(speciesInt_t, gemsb_species_is_gas)
  return gemsb_species_is_gas_real(e, j);
}

typedef double (*speciesDouble_t)(const gemsb_engine*, int);
double gemsb_species_molar_mass(const gemsb_engine* e, int j) {
  REAL(speciesDouble_t, gemsb_species_molar_mass)
  return gemsb_species_molar_mass_real(e, j);
}

typedef double (*stoich_t)(const gemsb_engine*, int, int);
double gemsb_stoich(const gemsb_engine* e, int j, int i) {
  REAL(stoich_t, gemsb_stoich)
  return gemsb_stoich_real(e, j, i);
}

typedef void (*range_t)(const gemsb_engine*, double*, double*, double*,
                        double*);
void gemsb_TP_range(const gemsb_engine* e, double* Tmin, double* Tmax,
                    double* Pmin, double* Pmax) {
  REAL(range_t, gemsb_TP_range)
  gemsb_TP_range_real(e, Tmin, Tmax, Pmin, Pmax);
}

typedef void (*bounds_t)(gemsb_engine*, int, double, double);
void gemsb_set_species_bounds(gemsb_engine* e, int j, double lower,
                              double upper) {
  char text[96];
  REAL(bounds_t, gemsb_set_species_bounds)
  snprintf(text, sizeof(text), "bounds %d %d %g %g", engineNumber(e), j,
           lower, upper);
  note(text);
  gemsb_set_species_bounds_real(e, j, lower, upper);
}

typedef int (*equilibrate_t)(gemsb_engine*, double, double, const double*,
                             int, double*, int*);
int gemsb_equilibrate(gemsb_engine* e, double T_K, double P_Pa,
                      const double* b_mol, int mode, double* state,
                      int* state_valid) {
  char text[160];
  int rc;
  REAL(equilibrate_t, gemsb_equilibrate)
  ++nEquilibrate;
  snprintf(text, sizeof(text), "equilibrate %ld engine %d state %p valid %d",
           nEquilibrate, engineNumber(e), (void*)state,
           state_valid ? *state_valid : -1);
  note(text);
  zeroPgasNow = listed(zeroPgasAt, nZeroPgas, nEquilibrate);
  nanOmegaNow = listed(nanOmegaAt, nNanOmega, nEquilibrate);
  if (listed(fatalAt, nFatal, nEquilibrate)) {
    note("fatal injected");
    return GEMSB_ERR_FATAL;
  }
  rc = gemsb_equilibrate_real(e, T_K, P_Pa, b_mol, mode, state, state_valid);
  snprintf(text, sizeof(text), "status %d", rc);
  note(text);
  if (zeroPgasNow) {
    note("zero p_g injected");
  }
  if (nanOmegaNow) {
    note("nan omega injected");
  }
  return rc;
}

typedef void (*results_t)(const gemsb_engine*, double*);
void gemsb_species_amounts(const gemsb_engine* e, double* n_mol) {
  REAL(results_t, gemsb_species_amounts)
  gemsb_species_amounts_real(e, n_mol);
}
void gemsb_gas_partial_pressures(const gemsb_engine* e, double* p_Pa) {
  int j;
  REAL(results_t, gemsb_gas_partial_pressures)
  gemsb_gas_partial_pressures_real(e, p_Pa);
  for (j = 0; j < gemsb_num_species(e); ++j) {
    p_Pa[j] = zeroPgasNow ? 0.0 : pgasScale*p_Pa[j];
  }
}
void gemsb_log10_activities(const gemsb_engine* e, double* log10a) {
  int j;
  REAL(results_t, gemsb_log10_activities)
  gemsb_log10_activities_real(e, log10a);
  for (j = 0; nanOmegaNow && j < gemsb_num_species(e); ++j) {
    if (!gemsb_species_is_gas(e, j)) {
      log10a[j] = NAN;
    }
  }
}
void gemsb_phase_amounts(const gemsb_engine* e, double* n_mol) {
  REAL(results_t, gemsb_phase_amounts)
  gemsb_phase_amounts_real(e, n_mol);
}
void gemsb_phase_log10_saturation(const gemsb_engine* e, double* lg) {
  REAL(results_t, gemsb_phase_log10_saturation)
  gemsb_phase_log10_saturation_real(e, lg);
}

typedef double (*seconds_t)(const gemsb_engine*);
double gemsb_last_seconds(const gemsb_engine* e) {
  REAL(seconds_t, gemsb_last_seconds)
  return gemsb_last_seconds_real(e);
}

typedef const char* (*error_t)(const gemsb_engine*);
const char* gemsb_last_error(const gemsb_engine* e) {
  REAL(error_t, gemsb_last_error)
  return gemsb_last_error_real(e);
}
