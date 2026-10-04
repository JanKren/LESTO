/*------------------------------------------------------------------------------

vapourPressureTable.H

PURPOSE

Equilibrium vapour pressure p_eq(T) [Pa] of the gas of a gas/condensate pair
from a tabulated csv file (plan section 2, E9), for the Hertz-Knudsen-Schrage
law of rhoFixedFlowFoam (interfaceExchange.H).  Standalone C++ in the style
of phaseChangeKinetics.H: no OpenFOAM types, errors are thrown as
std::runtime_error (interfaceExchange turns them into a FatalIOError), and
tests/testVapourPressureTable.C checks it with plain g++.

File format

  - lines starting with '#' are comments;
  - the '#' line right before the first data row names the columns,
    separated by commas (e.g. "# T_K,PbI2(cr),PbI2(l)");
  - every data row holds one value per column, separated by commas;
  - column 0 is the temperature T [K], strictly increasing;
  - the other columns are the vapour pressures of condensed PHASES of the
    same substance, in the unit of the dictionary (below).

  thermochemistry/systems/make_pv_table.py writes data/pv_PbI2_phases.csv
  (Gurvich via NASA-9, log10 Pa, one column per phase every 5 K plus the
  melting point) and data/pv_PbI2_TFlows.csv (the T-Flows table, bar).

Dictionary (constant/thermochemistryProperties, pairs { <gas> {
vapourPressure { ... } } })

  file           "<constant>/pv_PbI2_phases.csv";
  phases         ("PbI2(cr)" "PbI2(l)");  // columns; p_eq = their minimum
  units          log10Pa;       // log10Pa | Pa | bar
  interpolation  logInverseT;   // logInverseT (default) | linearT (T-Flows)
  outOfRange     extrapolate;   // extrapolate (default) | hold (T-Flows) |
                                // fatal

Phases.  Each listed column is the vapour pressure over one condensed phase
over the whole range (a phase is metastable where its curve lies above
another one).  The stable condensate has the lowest vapour pressure, so

  p_eq(T) = min over the listed phases of p_phase(T),

and the kink at a phase transition (683 K for PbI2, crystal -> liquid) is
exact whatever the node spacing: each phase is interpolated on its own
smooth curve.  crossings() returns the temperatures where the stable phase
changes (the melting point of the data).

Interpolation of one phase between the nodes T_i < T < T_i+1

  logInverseT  log10 p linear in 1/T (Clausius-Clapeyron: exact for a
               constant enthalpy of sublimation or vaporisation);
               the plan's default;
  linearT      p linear in T, as the T-Flows reference implementation
               (code-to-code comparison only).

At a node T = T_i the table returns the node's value p_i exactly, in either
mode and unit (the value converted to Pa once when the file is read).

Outside [T_min, T_max]

  extrapolate  the end interval's law continued (logInverseT: Clausius-
               Clapeyron through the two end nodes; linearT: the straight
               line, not below 0);
  hold         the end value (T-Flows);
  fatal        an error naming the temperature and the range.

Units of the values: log10Pa (log10 of p in Pa), Pa, bar (1e5 Pa).
logInverseT needs p > 0 at every node; linearT accepts p >= 0.

------------------------------------------------------------------------------*/

#ifndef vapourPressureTable_H
#define vapourPressureTable_H

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace LESTO {

class vapourPressureTable {

public:

  enum interpolationType { LOG_INVERSE_T, LINEAR_T };
  enum unitType { LOG10_PA, PA, BAR };
  enum rangeType { EXTRAPOLATE, HOLD, FATAL };

private:

  std::string                      file_;
  std::vector<std::string>         phases_;
  std::vector<double>              T_;       /* nodes [K], increasing */
  std::vector<std::vector<double>> p_;       /* per phase, per node [Pa] */
  std::vector<std::vector<double>> log10p_;  /* log10(p/Pa); empty for
                                                linearT */
  interpolationType                interpolation_;
  rangeType                        outOfRange_;

  static std::string trim(const std::string& s) {
    const std::size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) {
      return "";
    }
    const std::size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
  }

  static std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) {
      out.push_back(trim(item));
    }
    if (!line.empty() && line.back() == ',') {
      out.push_back("");
    }
    return out;
  }

  /* the value of a data field; the whole field must be a number */
  static double number(const std::string& text, const std::string& where) {
    const char* begin = text.c_str();
    char* end = nullptr;
    errno = 0;
    const double value = std::strtod(begin, &end);
    if (text.empty() || end != begin + text.size() || errno == ERANGE
      || !std::isfinite(value)) {
      throw std::runtime_error("vapourPressureTable: '" + text
        + "' is not a finite number (" + where + ")");
    }
    return value;
  }

  static std::string text(const double x) {
    std::ostringstream os;
    os.precision(10);
    os << x;
    return os.str();
  }

  /* the interval [T_i, T_i+1] used for T (the end intervals outside) */
  std::size_t interval(const double T) const {
    const std::size_t n = T_.size();
    if (!(T > T_[0])) {
      return 0;
    }
    if (T >= T_[n - 1]) {
      return n - 2;
    }
    std::size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      const std::size_t mid = (lo + hi)/2;
      if (T_[mid] > T) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    return lo;
  }

public:

  /*--------------------------------------------------------------------------
  From node arrays: T [K] strictly increasing, values[phase][node] in the
  given unit.  Throws on inconsistent input.
  --------------------------------------------------------------------------*/
  vapourPressureTable (
    const std::vector<double>&              T,
    const std::vector<std::vector<double>>& values,
    const std::vector<std::string>&         phases,
    const unitType                          units,
    const interpolationType                 interpolation,
    const rangeType                         outOfRange,
    const std::string&                      file = "(arrays)"
  )
  :
    file_(file),
    phases_(phases),
    T_(T),
    p_(),
    log10p_(),
    interpolation_(interpolation),
    outOfRange_(outOfRange)
  {
    if (T_.size() < 2) {
      throw std::runtime_error("vapourPressureTable: " + file_
        + " needs at least two temperatures");
    }
    if (phases_.empty() || values.size() != phases_.size()) {
      throw std::runtime_error("vapourPressureTable: " + file_
        + ": one column of values per listed phase is needed");
    }
    for (std::size_t i = 0; i < T_.size(); ++i) {
      if (!(T_[i] > 0) || (i > 0 && !(T_[i] > T_[i - 1]))) {
        throw std::runtime_error("vapourPressureTable: " + file_
          + ": the temperatures must be positive and strictly increasing"
          + " (row " + std::to_string(i + 1) + ", T = " + text(T_[i])
          + " K)");
      }
    }

    for (std::size_t k = 0; k < phases_.size(); ++k) {
      if (values[k].size() != T_.size()) {
        throw std::runtime_error("vapourPressureTable: " + file_
          + ": phase " + phases_[k] + " has " + std::to_string(
          values[k].size()) + " values for " + std::to_string(T_.size())
          + " temperatures");
      }
      std::vector<double> p(T_.size()), l(T_.size());
      for (std::size_t i = 0; i < T_.size(); ++i) {
        const double v = values[k][i];
        p[i] = units == LOG10_PA ? std::pow(10.0, v)
             : units == BAR ? v*1e5
             : v;
        l[i] = units == LOG10_PA ? v : std::log10(p[i]);
        const bool valid = interpolation_ == LOG_INVERSE_T
          ? p[i] > 0 && std::isfinite(l[i])
          : p[i] >= 0;
        if (!valid || !std::isfinite(p[i])) {
          throw std::runtime_error("vapourPressureTable: " + file_
            + ": phase " + phases_[k] + " at T = " + text(T_[i])
            + " K: the vapour pressure " + text(p[i]) + " Pa must be "
            + (interpolation_ == LOG_INVERSE_T
               ? "positive (interpolation logInverseT)" : "finite, >= 0"));
        }
      }
      p_.push_back(p);
      log10p_.push_back(l);
    }
  }

  /*--------------------------------------------------------------------------
  Reads the file (see PURPOSE) and selects the columns of the phases.
  --------------------------------------------------------------------------*/
  static vapourPressureTable read (
    const std::string&              file,
    const std::vector<std::string>& phases,
    const unitType                  units,
    const interpolationType         interpolation,
    const rangeType                 outOfRange
  ) {
    std::ifstream is(file);
    if (!is) {
      throw std::runtime_error("vapourPressureTable: cannot open " + file);
    }

    std::string line, header;
    std::vector<std::string> names;
    std::vector<std::size_t> columns;
    std::vector<double> T;
    std::vector<std::vector<double>> values(phases.size());
    long row = 0;

    while (std::getline(is, line)) {
      ++row;
      const std::string content = trim(line);
      if (content.empty()) {
        continue;
      }
      if (content[0] == '#') {
        if (names.empty()) {
          header = content.substr(1);
        }
        continue;
      }

      const std::string where = file + ", line " + std::to_string(row);

      /* the column names: the '#' line right before the first data row */
      if (names.empty()) {
        names = split(header);
        if (names.size() < 2) {
          throw std::runtime_error("vapourPressureTable: " + file
            + ": no column names (a '#' line with comma-separated names"
            + " right before the first data row)");
        }
        for (const std::string& phase : phases) {
          std::size_t found = 0;
          for (std::size_t c = 1; c < names.size(); ++c) {
            if (names[c] == phase) {
              found = c;
            }
          }
          if (found == 0) {
            std::string list;
            for (std::size_t c = 1; c < names.size(); ++c) {
              list += (c > 1 ? ", " : "") + names[c];
            }
            throw std::runtime_error("vapourPressureTable: " + file
              + ": no column '" + phase + "' (columns: " + list + ")");
          }
          columns.push_back(found);
        }
      }

      const std::vector<std::string> fields = split(content);
      if (fields.size() != names.size()) {
        throw std::runtime_error("vapourPressureTable: " + where + ": "
          + std::to_string(fields.size()) + " values for "
          + std::to_string(names.size()) + " columns");
      }
      T.push_back(number(fields[0], where));
      for (std::size_t k = 0; k < phases.size(); ++k) {
        values[k].push_back(number(fields[columns[k]], where));
      }
    }

    if (names.empty()) {
      throw std::runtime_error("vapourPressureTable: " + file
        + " holds no data row");
    }

    return vapourPressureTable (
      T, values, phases, units, interpolation, outOfRange, file
    );
  }

  /* dictionary words of the options; they throw on an unknown word */
  static unitType units(const std::string& word) {
    if (word == "log10Pa") return LOG10_PA;
    if (word == "Pa") return PA;
    if (word == "bar") return BAR;
    throw std::runtime_error("vapourPressure: unknown units " + word
      + ".  Choose log10Pa, Pa or bar.");
  }

  static interpolationType interpolation(const std::string& word) {
    if (word == "logInverseT") return LOG_INVERSE_T;
    if (word == "linearT") return LINEAR_T;
    throw std::runtime_error("vapourPressure: unknown interpolation " + word
      + ".  Choose logInverseT or linearT.");
  }

  static rangeType outOfRange(const std::string& word) {
    if (word == "extrapolate") return EXTRAPOLATE;
    if (word == "hold") return HOLD;
    if (word == "fatal") return FATAL;
    throw std::runtime_error("vapourPressure: unknown outOfRange " + word
      + ".  Choose extrapolate, hold or fatal.");
  }

  const std::string& file() const { return file_; }
  const std::vector<std::string>& phases() const { return phases_; }
  std::size_t nNodes() const { return T_.size(); }

  /* the nodes [K] and the interpolation: the GEMS table of mode frozen is
     built at the nodes of the pair's table and interpolated alike
     (gemsEquilibrium.H) */
  const std::vector<double>& nodes() const { return T_; }
  interpolationType interpolationMode() const { return interpolation_; }
  double Tmin() const { return T_.front(); }
  double Tmax() const { return T_.back(); }
  bool inRange(const double T) const { return T >= Tmin() && T <= Tmax(); }

  /*--------------------------------------------------------------------------
  Vapour pressure [Pa] over phase k at T [K].
  --------------------------------------------------------------------------*/
  double phase(const std::size_t k, const double T) const {

    if (!(T > 0) || !std::isfinite(T)) {
      throw std::runtime_error("vapourPressureTable: " + file_
        + ": temperature " + text(T) + " K");
    }

    double t = T;
    if (!inRange(T)) {
      if (outOfRange_ == FATAL) {
        throw std::runtime_error("vapourPressureTable: " + file_ + ": T = "
          + text(T) + " K lies outside the table [" + text(Tmin()) + ", "
          + text(Tmax()) + "] K (outOfRange fatal)");
      }
      if (outOfRange_ == HOLD) {
        t = T < Tmin() ? Tmin() : Tmax();
      }
    }

    const std::size_t i = interval(t);
    const std::vector<double>& p = p_[k];

    /* the nodes themselves: exactly the value of the file */
    if (t == T_[i]) {
      return p[i];
    }
    if (t == T_[i + 1]) {
      return p[i + 1];
    }

    if (interpolation_ == LOG_INVERSE_T) {
      const std::vector<double>& l = log10p_[k];
      const double w = (1.0/t - 1.0/T_[i])/(1.0/T_[i + 1] - 1.0/T_[i]);
      return std::pow(10.0, (1.0 - w)*l[i] + w*l[i + 1]);
    }

    const double w = (t - T_[i])/(T_[i + 1] - T_[i]);
    const double value = (1.0 - w)*p[i] + w*p[i + 1];
    return value > 0 ? value : 0.0;
  }

  /* index of the stable phase at T: the lowest vapour pressure (the first
     of equal ones) */
  std::size_t stablePhase(const double T) const {
    std::size_t best = 0;
    double pBest = phase(0, T);
    for (std::size_t k = 1; k < phases_.size(); ++k) {
      const double pk = phase(k, T);
      if (pk < pBest) {
        best = k;
        pBest = pk;
      }
    }
    return best;
  }

  /* p_eq [Pa]: the minimum over the phases */
  double operator()(const double T) const {
    return phase(stablePhase(T), T);
  }

  /*--------------------------------------------------------------------------
  The temperatures [K] inside [T_min, T_max] where the stable phase changes
  (e.g. the melting point), with the phases before and after: bisection to
  1e-9 K on the interpolated curves of the two phases in every interval
  whose ends have different stable phases.
  --------------------------------------------------------------------------*/
  struct crossing {
    double      T;
    std::size_t below;
    std::size_t above;
  };

  std::vector<crossing> crossings() const {
    std::vector<crossing> out;
    for (std::size_t i = 0; i + 1 < T_.size(); ++i) {
      const std::size_t a = stablePhase(T_[i]);
      const std::size_t b = stablePhase(T_[i + 1]);
      if (a == b) {
        continue;
      }
      double lo = T_[i], hi = T_[i + 1];
      for (int iter = 0; iter < 200 && hi - lo > 1e-9; ++iter) {
        const double mid = 0.5*(lo + hi);
        if (phase(a, mid) <= phase(b, mid)) {
          lo = mid;
        } else {
          hi = mid;
        }
      }
      out.push_back(crossing{0.5*(lo + hi), a, b});
    }
    return out;
  }
};

}

#endif
