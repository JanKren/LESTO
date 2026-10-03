// gemsbridge.cpp -- implementation of the plain-C GEMS3K bridge (see gemsbridge.h)
// Compiled with the conda-forge g++ (C++20).  No C++ type or exception may
// cross the extern "C" boundary: every entry point catches everything.
#include "gemsbridge.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

#include "node.h"
#include "jsonconfig.h"

#define GEMSB_API extern "C" __attribute__((visibility("default")))

struct gemsb_engine {
  std::unique_ptr<TNode> node;
  long nIC = 0, nDC = 0, nPH = 0, nPS = 0;
  long gasPhase = -1;                 // DBR index of the gas phase
  std::vector<double> dul, dll, aPH;  // metastability limits, surface areas
  std::vector<double> bn;             // normalised element amounts
  std::vector<char> isGas;
  // cached, scaled-back results of the last successful call
  std::vector<double> n, pgas, lga, nph, lgom;
  double scale = 1.0;
  int iters = 0;
  double secs = 0.0;
  std::string err;
};

static void copy_err(const std::string& s, char* buf, int len) {
  if (buf && len > 0) { std::strncpy(buf, s.c_str(), len - 1); buf[len - 1] = 0; }
}

static void setup(gemsb_engine* e) {
  DATACH* ch = e->node->pCSD();
  DATABR* br = e->node->pCNode();
  e->nIC = ch->nICb; e->nDC = ch->nDCb; e->nPH = ch->nPHb; e->nPS = ch->nPSb;
  // G0(T) is re-interpolated only when |T - T_prev_call| >= Ttol (node_copy.cpp:547),
  // and T_prev is updated on every call, so a slow T drift between neighbouring
  // cells never triggers a reload.  Clamp Ttol (interpolation mode only).
  if (ch->mLook == 0 && ch->Ttol > 1e-3) ch->Ttol = 1e-3;
  e->dul.assign(br->dul, br->dul + e->nDC);
  e->dll.assign(br->dll, br->dll + e->nDC);
  e->aPH.assign(e->nPH, 0.0);
  if (ch->nAalp > 0) e->aPH.assign(br->aPH, br->aPH + e->nPH);
  e->bn.assign(e->nIC, 0.0);
  e->isGas.assign(e->nDC, 0);
  for (long k = 0; k < e->nPH; ++k)
    if (ch->ccPH[ch->xph[k]] == 'g') e->gasPhase = k;
  for (long j = 0; j < e->nDC; ++j)
    e->isGas[j] = (e->node->DCtoPh_DBR(j) == e->gasPhase) ? 1 : 0;
  e->n.assign(e->nDC, 0.0); e->pgas.assign(e->nDC, 0.0); e->lga.assign(e->nDC, 0.0);
  e->nph.assign(e->nPH, 0.0); e->lgom.assign(e->nPH, 0.0);
}

// ---- per-cell warm-start state: exactly the DATABR fields that
//      TNode::unpackDataBr(uPrimalSol=true) reads (node.cpp:399-507) -------------
static long state_size(const gemsb_engine* e) {
  return 2 * e->nDC + 2 * e->nPH + 2 * e->nIC + 3 * e->nPS + e->nPS * e->nIC + 3;
}
template <class F> static void state_walk(gemsb_engine* e, DATABR* b, F&& f) {
  f(b->xDC, e->nDC); f(b->gam, e->nDC); f(b->xPH, e->nPH); f(b->omPH, e->nPH);
  f(b->uIC, e->nIC); f(b->rMB, e->nIC); f(b->vPS, e->nPS); f(b->mPS, e->nPS);
  f(b->xPA, e->nPS); f(b->bPS, e->nPS * e->nIC);
  f(&b->Ms, 1); f(&b->Vs, 1); f(&b->IC, 1);
}
static void state_load(gemsb_engine* e, const double* s) {
  long o = 0;
  state_walk(e, e->node->pCNode(), [&](double* p, long n) { std::copy(s + o, s + o + n, p); o += n; });
}
static void state_save(gemsb_engine* e, double* s) {
  long o = 0;
  state_walk(e, e->node->pCNode(), [&](double* p, long n) { std::copy(p, p + n, s + o); o += n; });
}

GEMSB_API void gemsb_set_log_directory(const char* dir) {
  try {
    std::string d = dir ? dir : "";
    if (!d.empty() && d.back() != '/') d += '/';
    GemsSettings::data_logger_directory = d;
  } catch (...) {}
}

GEMSB_API void gemsb_set_log_level(int level) {
  try { gemsSettings().gems3k_update_log_level((size_t)level); } catch (...) {}
}

static gemsb_engine* finish_create(std::unique_ptr<gemsb_engine> e, long rc, char* errmsg, int errlen) {
  if (rc != 0) { copy_err("GEM_init failed, code " + std::to_string(rc), errmsg, errlen); return nullptr; }
  setup(e.get());
  return e.release();
}

GEMSB_API gemsb_engine* gemsb_create_from_lst(const char* lst, char* errmsg, int errlen) {
  try {
    auto e = std::make_unique<gemsb_engine>();
    e->node = std::make_unique<TNode>();
    long rc = e->node->GEM_init(lst);
    return finish_create(std::move(e), rc, errmsg, errlen);
  } catch (std::exception& x) { copy_err(x.what(), errmsg, errlen); }
  catch (...) { copy_err("unknown exception in gemsb_create_from_lst", errmsg, errlen); }
  return nullptr;
}

GEMSB_API gemsb_engine* gemsb_create_from_strings(const char* dch, const char* ipm, const char* dbr,
                                                  char* errmsg, int errlen) {
  try {
    auto e = std::make_unique<gemsb_engine>();
    e->node = std::make_unique<TNode>();
    long rc = e->node->GEM_init(std::string(dch), std::string(ipm), std::string(dbr));
    return finish_create(std::move(e), rc, errmsg, errlen);
  } catch (std::exception& x) { copy_err(x.what(), errmsg, errlen); }
  catch (...) { copy_err("unknown exception in gemsb_create_from_strings", errmsg, errlen); }
  return nullptr;
}

GEMSB_API void gemsb_destroy(gemsb_engine* e) { try { delete e; } catch (...) {} }

GEMSB_API int gemsb_num_elements(const gemsb_engine* e) { return (int)e->nIC; }
GEMSB_API int gemsb_num_species(const gemsb_engine* e)  { return (int)e->nDC; }
GEMSB_API int gemsb_num_phases(const gemsb_engine* e)   { return (int)e->nPH; }
GEMSB_API int gemsb_state_size(const gemsb_engine* e)   { return (int)state_size(e); }
GEMSB_API const char* gemsb_element_name(const gemsb_engine* e, int i) {
  return e->node->xCH_to_IC_name(e->node->IC_xDB_to_xCH(i)).c_str(); }
GEMSB_API const char* gemsb_species_name(const gemsb_engine* e, int j) {
  return e->node->xCH_to_DC_name(e->node->DC_xDB_to_xCH(j)).c_str(); }
GEMSB_API const char* gemsb_phase_name(const gemsb_engine* e, int k) {
  return e->node->xCH_to_PH_name(e->node->Ph_xDB_to_xCH(k)).c_str(); }
GEMSB_API int gemsb_element_index(const gemsb_engine* e, const char* s) {
  try { return (int)e->node->IC_name_to_xDB(s); } catch (...) { return -1; } }
GEMSB_API int gemsb_species_index(const gemsb_engine* e, const char* s) {
  try { return (int)e->node->DC_name_to_xDB(s); } catch (...) { return -1; } }
GEMSB_API int gemsb_phase_index(const gemsb_engine* e, const char* s) {
  try { return (int)e->node->Ph_name_to_xDB(s); } catch (...) { return -1; } }
GEMSB_API int gemsb_species_phase(const gemsb_engine* e, int j) { return (int)e->node->DCtoPh_DBR(j); }
GEMSB_API int gemsb_species_is_gas(const gemsb_engine* e, int j) { return e->isGas[j]; }
GEMSB_API double gemsb_species_molar_mass(const gemsb_engine* e, int j) { return e->node->DCmm(j); }
GEMSB_API double gemsb_stoich(const gemsb_engine* e, int j, int i) { return e->node->DCaJI(j, i); }
GEMSB_API void gemsb_TP_range(const gemsb_engine* e, double* T0, double* T1, double* P0, double* P1) {
  const DATACH* c = e->node->pCSD();
  if (T0) *T0 = c->TKval[0]; if (T1) *T1 = c->TKval[c->nTp - 1];
  if (P0) *P0 = c->Pval[0];  if (P1) *P1 = c->Pval[c->nPp - 1];
}

GEMSB_API void gemsb_set_species_bounds(gemsb_engine* e, int j, double lo, double up) {
  if (e && j >= 0 && j < e->nDC && lo <= up) { e->dll[j] = lo; e->dul[j] = up; }
}

static int map_status(long st) {
  switch (st) {
    case OK_GEM_AIA: case OK_GEM_SIA: return GEMSB_OK;
    case BAD_GEM_AIA: case BAD_GEM_SIA: return GEMSB_BAD_QUALITY;
    case ERR_GEM_AIA: case ERR_GEM_SIA: return GEMSB_ERR_NOCONV;
    default: return GEMSB_ERR_FATAL;
  }
}

static long run_once(gemsb_engine* e, double T, double P, int mode, const double* st) {
  TNode* nd = e->node.get();
  long want = (mode == GEMSB_COLD) ? NEED_GEM_AIA : NEED_GEM_SIA;
  if (mode == GEMSB_WARM_CELL && st) state_load(e, st);
  // (8c) minimal input: T, P, bIC, dul, dll   (node.h:355-364)
  nd->GEM_from_MT(0, want, T, P, e->bn.data(), e->dul.data(), e->dll.data());
  // uPrimalSol=true -> take speciation from DATABR (per-cell state);
  // false -> keep MULTI from the previous call (node.h:454-465)
  return nd->GEM_run(mode == GEMSB_WARM_CELL && st);
}

GEMSB_API int gemsb_equilibrate(gemsb_engine* e, double T, double P, const double* b, int mode,
                                double* state, int* state_valid) {
  if (!e || !b) return GEMSB_ERR_INPUT;
  try {
    TNode* nd = e->node.get();
    DATACH* ch = nd->pCSD();
    e->err.clear();
    if (!nd->check_TP(T, P)) { e->err = "T/P outside DCH lookup grid"; return GEMSB_ERR_INPUT; }
    // normalise the bulk composition to 1 mol of elements: GEMS3K cut-offs
    // (DS, DcMin, PhMin, DB) are absolute amounts in mol (ipm_main.cpp:1896-1922)
    double sum = 0.0;
    for (long i = 0; i < e->nIC; ++i)
      if (ch->ccIC[ch->xic[i]] != 'z') sum += std::max(0.0, b[i]);
    if (!(sum > 0.0) || !std::isfinite(sum)) { e->err = "non-positive bulk composition"; return GEMSB_ERR_INPUT; }
    e->scale = sum;
    const double floor_rel = 1e-15;   // > pa_DB (1e-17 mol); node.cpp unpackDataBr check
    for (long i = 0; i < e->nIC; ++i)
      e->bn[i] = (ch->ccIC[ch->xic[i]] == 'z') ? 0.0 : std::max(b[i] / sum, floor_rel);

    bool have_state = (mode == GEMSB_WARM_CELL && state && state_valid && *state_valid);
    int m = (mode == GEMSB_WARM_CELL && !have_state) ? GEMSB_COLD : mode;
    auto t0 = std::chrono::steady_clock::now();
    long st = run_once(e, T, P, m, have_state ? state : nullptr);
    int rc = map_status(st);
    if (rc == GEMSB_ERR_FATAL) {
      e->err = nd->ipmLogError();
      return rc;
    }
    if (rc != GEMSB_OK && m != GEMSB_COLD) {           // fall back to a cold start
      st = run_once(e, T, P, GEMSB_COLD, nullptr);
      rc = map_status(st);
      if (rc == GEMSB_OK) rc = GEMSB_OK_RETRIED;
    }
    auto t1 = std::chrono::steady_clock::now();
    e->secs = std::chrono::duration<double>(t1 - t0).count();
    e->iters = (int)nd->pCNode()->IterDone;
    if (rc == GEMSB_ERR_NOCONV || rc == GEMSB_ERR_FATAL) { e->err = nd->ipmLogError(); return rc; }

    // cache results, scaled back to the caller's amounts
    DATABR* br = nd->pCNode();
    double xgas = (e->gasPhase >= 0) ? br->xPH[e->gasPhase] : 0.0;
    for (long j = 0; j < e->nDC; ++j) {
      e->n[j] = br->xDC[j] * e->scale;
      e->pgas[j] = (e->isGas[j] && xgas > 0.0) ? br->xDC[j] / xgas * P : 0.0;
      e->lga[j] = nd->Get_aDC(j, false);
    }
    for (long k = 0; k < e->nPH; ++k) { e->nph[k] = br->xPH[k] * e->scale; e->lgom[k] = nd->Ph_SatInd(k); }
    if (state && state_valid && rc != GEMSB_BAD_QUALITY) { state_save(e, state); *state_valid = 1; }
    return rc;
  } catch (std::exception& x) { e->err = x.what(); return GEMSB_ERR_FATAL; }
  catch (...) { e->err = "unknown exception"; return GEMSB_ERR_FATAL; }
}

GEMSB_API void gemsb_species_amounts(const gemsb_engine* e, double* v) { std::copy(e->n.begin(), e->n.end(), v); }
GEMSB_API void gemsb_gas_partial_pressures(const gemsb_engine* e, double* v) { std::copy(e->pgas.begin(), e->pgas.end(), v); }
GEMSB_API void gemsb_log10_activities(const gemsb_engine* e, double* v) { std::copy(e->lga.begin(), e->lga.end(), v); }
GEMSB_API void gemsb_phase_amounts(const gemsb_engine* e, double* v) { std::copy(e->nph.begin(), e->nph.end(), v); }
GEMSB_API void gemsb_phase_log10_saturation(const gemsb_engine* e, double* v) { std::copy(e->lgom.begin(), e->lgom.end(), v); }
GEMSB_API int gemsb_last_iterations(const gemsb_engine* e) { return e->iters; }
GEMSB_API double gemsb_last_seconds(const gemsb_engine* e) { return e->secs; }
GEMSB_API const char* gemsb_last_error(const gemsb_engine* e) { return e->err.c_str(); }
