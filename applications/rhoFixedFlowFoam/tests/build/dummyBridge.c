/*------------------------------------------------------------------------------

tests/build/dummyBridge.c

PURPOSE

A dummy GEMS3K bridge for the build tests of tests/build/Allrun (M0.3b, c,
f, h, i): every function of thermochemistry/gemsbridge/gemsbridge.h,
without GEMS3K.  Since milestone M9 the solver's gemsEquilibrium.C calls
the bridge's whole API, so a bridge build compiles and links against the
real header and needs every symbol; the plumbing (the stamp, wclean, the
generated gemsConfig.H, -lgemsbridge, the run path, ldd) is tested just as
well as with the real library, which a fresh clone does not hold (the .so
is git-ignored).  Nothing runs it: an engine cannot be created
(gemsb_create_from_lst returns NULL with a message), so a run that selects
equilibrium GEMS against this library stops at start-up.

Built by tests/build/Allrun with the C compiler of OpenFOAM into a
directory of the work area, next to a copy of the real gemsbridge.h.

------------------------------------------------------------------------------*/

#include "gemsbridge.h"

#include <stddef.h>
#include <string.h>

int gemsb_dummy_version(void) { return 0; }

static void message(char* errmsg, int errlen) {
  const char* text = "dummy GEMS3K bridge of tests/build: no GEMS3K inside";
  if (errmsg && errlen > 0) {
    strncpy(errmsg, text, (size_t)errlen - 1);
    errmsg[errlen - 1] = 0;
  }
}

void gemsb_set_log_directory(const char* dir) { (void)dir; }
void gemsb_set_log_level(int level) { (void)level; }

gemsb_engine* gemsb_create_from_lst(const char* dat_lst, char* errmsg,
                                    int errlen) {
  (void)dat_lst;
  message(errmsg, errlen);
  return NULL;
}

gemsb_engine* gemsb_create_from_strings(const char* dch, const char* ipm,
                                        const char* dbr, char* errmsg,
                                        int errlen) {
  (void)dch; (void)ipm; (void)dbr;
  message(errmsg, errlen);
  return NULL;
}

void gemsb_destroy(gemsb_engine* e) { (void)e; }

int gemsb_num_elements(const gemsb_engine* e) { (void)e; return 0; }
int gemsb_num_species(const gemsb_engine* e) { (void)e; return 0; }
int gemsb_num_phases(const gemsb_engine* e) { (void)e; return 0; }

const char* gemsb_element_name(const gemsb_engine* e, int i) {
  (void)e; (void)i; return "";
}
const char* gemsb_species_name(const gemsb_engine* e, int j) {
  (void)e; (void)j; return "";
}
const char* gemsb_phase_name(const gemsb_engine* e, int k) {
  (void)e; (void)k; return "";
}

int gemsb_element_index(const gemsb_engine* e, const char* name) {
  (void)e; (void)name; return -1;
}
int gemsb_species_index(const gemsb_engine* e, const char* name) {
  (void)e; (void)name; return -1;
}
int gemsb_phase_index(const gemsb_engine* e, const char* name) {
  (void)e; (void)name; return -1;
}

int gemsb_species_phase(const gemsb_engine* e, int j) {
  (void)e; (void)j; return -1;
}
int gemsb_species_is_gas(const gemsb_engine* e, int j) {
  (void)e; (void)j; return 0;
}
double gemsb_species_molar_mass(const gemsb_engine* e, int j) {
  (void)e; (void)j; return 0;
}
double gemsb_stoich(const gemsb_engine* e, int j, int i) {
  (void)e; (void)j; (void)i; return 0;
}

void gemsb_TP_range(const gemsb_engine* e, double* Tmin, double* Tmax,
                    double* Pmin, double* Pmax) {
  (void)e;
  if (Tmin) *Tmin = 0;
  if (Tmax) *Tmax = 0;
  if (Pmin) *Pmin = 0;
  if (Pmax) *Pmax = 0;
}

int gemsb_state_size(const gemsb_engine* e) { (void)e; return 0; }

void gemsb_set_species_bounds(gemsb_engine* e, int j, double lower,
                              double upper) {
  (void)e; (void)j; (void)lower; (void)upper;
}

int gemsb_equilibrate(gemsb_engine* e, double T_K, double P_Pa,
                      const double* b_mol, int mode, double* state,
                      int* state_valid) {
  (void)e; (void)T_K; (void)P_Pa; (void)b_mol; (void)mode; (void)state;
  (void)state_valid;
  return GEMSB_ERR_FATAL;
}

void gemsb_species_amounts(const gemsb_engine* e, double* n_mol) {
  (void)e; (void)n_mol;
}
void gemsb_gas_partial_pressures(const gemsb_engine* e, double* p_Pa) {
  (void)e; (void)p_Pa;
}
void gemsb_log10_activities(const gemsb_engine* e, double* log10a) {
  (void)e; (void)log10a;
}
void gemsb_phase_amounts(const gemsb_engine* e, double* n_mol) {
  (void)e; (void)n_mol;
}
void gemsb_phase_log10_saturation(const gemsb_engine* e, double* lg) {
  (void)e; (void)lg;
}

int gemsb_last_iterations(const gemsb_engine* e) { (void)e; return 0; }
double gemsb_last_seconds(const gemsb_engine* e) { (void)e; return 0; }
const char* gemsb_last_error(const gemsb_engine* e) {
  (void)e; return "dummy GEMS3K bridge";
}
