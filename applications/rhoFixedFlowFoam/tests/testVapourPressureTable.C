/*------------------------------------------------------------------------------

tests/testVapourPressureTable.C

PURPOSE

Standalone test of vapourPressureTable.H (no OpenFOAM libraries needed),
criterion M5.4 of tests/Alltest M5 (plan section 6, M5 acceptance 4):

  - synthetic tables: every node is returned exactly in every unit and
    interpolation; logInverseT (log10 p linear in 1/T) and linearT (p
    linear in T) between the nodes; extrapolate (Clausius-Clapeyron through
    the end nodes, the straight line not below 0), hold and fatal outside;
    p_eq = the minimum over the phases with its crossing; the errors of
    malformed files and unknown options;
  - the phase-resolved Gurvich table data/pv_PbI2_phases.csv (argument 1):
    every node exact, i.e. p_eq(T_i) = min over the phases of the value of
    the file; the melting point from the crossing of the phases 683.0 +-
    0.1 K; the table (5 K nodes plus T_m, logInverseT) against NASA-9
    (argument 2: T and log10 p_eq every 0.25 K, make_pv_table.py
    --reference) below 2e-3 decades, not checked when argument 2 is '-'
    (a machine without thermo.inp: criterion M5.4a);
  - the T-Flows table data/pv_PbI2_TFlows.csv (argument 3, bar, linearT,
    hold): its 95 nodes 283.15 + 10 k K reproduced to 1e-12, the end values
    held outside, and the mean of two neighbouring nodes at the midpoint.

Build and run from the solver directory:

  g++ -std=c++17 -Wall -Wextra tests/testVapourPressureTable.C \
    -o <dir>/testVapourPressureTable
  cd <dir> && ./testVapourPressureTable pv_PbI2_phases.csv reference.csv \
    pv_PbI2_TFlows.csv

(reference.csv '-': without the comparison with NASA-9).

It writes its small synthetic csv files into the current directory.
The last line of the output is "<n> checks, <m> failed."; the exit status is
0 when nothing failed.

------------------------------------------------------------------------------*/

#include "../vapourPressureTable.H"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

  int nFailed = 0;
  int nChecks = 0;

  void check(const bool ok, const std::string& what) {
    ++nChecks;
    if (!ok) {
      ++nFailed;
      std::printf("FAILED: %s\n", what.c_str());
    }
  }

  bool close(const double a, const double b, const double tol) {
    return std::abs(a - b) <= tol*std::abs(b);
  }

  /* true when the call throws std::runtime_error with the given text */
  template<class F>
  bool throwsWith(F f, const std::string& text) {
    try {
      f();
    } catch (const std::runtime_error& e) {
      return std::string(e.what()).find(text) != std::string::npos;
    }
    return false;
  }

  std::string writeFile(const std::string& name, const std::string& text) {
    const std::string path = "testVapourPressureTable_" + name;
    std::ofstream(path) << text;
    return path;
  }

  /* the numeric rows of a csv file of this format (comments skipped) */
  std::vector<std::vector<double>> rows(const std::string& path) {
    std::vector<std::vector<double>> out;
    std::ifstream is(path);
    std::string line;
    while (std::getline(is, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }
      std::vector<double> row;
      std::stringstream ss(line);
      std::string item;
      while (std::getline(ss, item, ',')) {
        row.push_back(std::strtod(item.c_str(), nullptr));
      }
      out.push_back(row);
    }
    return out;
  }
}

int main(int argc, char* argv[]) {

  using LESTO::vapourPressureTable;
  typedef vapourPressureTable table;

  /*--------------------------------------------------------------------------
  Synthetic one-phase table: log10 p = 10 - 3000/T exactly (Clausius-
  Clapeyron), nodes 500, 600, 800 K.
  --------------------------------------------------------------------------*/
  const std::vector<double> T3 = {500, 600, 800};
  auto cc = [](const double T) { return 10.0 - 3000.0/T; };
  std::vector<double> l3, p3, b3;
  for (const double T : T3) {
    l3.push_back(cc(T));
    p3.push_back(std::pow(10.0, cc(T)));
    b3.push_back(std::pow(10.0, cc(T))/1e5);
  }

  const table logT(T3, {l3}, {"s"}, table::LOG10_PA, table::LOG_INVERSE_T,
                   table::EXTRAPOLATE);
  const table pascal(T3, {p3}, {"s"}, table::PA, table::LOG_INVERSE_T,
                     table::EXTRAPOLATE);
  const table bar(T3, {b3}, {"s"}, table::BAR, table::LINEAR_T,
                  table::HOLD);
  const table fatal(T3, {l3}, {"s"}, table::LOG10_PA, table::LINEAR_T,
                    table::FATAL);

  /* nodes: exactly the value of the file converted once */
  for (std::size_t i = 0; i < T3.size(); ++i) {
    check(logT(T3[i]) == std::pow(10.0, l3[i]), "node exact, log10Pa");
    check(pascal(T3[i]) == p3[i], "node exact, Pa with logInverseT");
    check(bar(T3[i]) == b3[i]*1e5, "node exact, bar with linearT");
    check(fatal(T3[i]) == std::pow(10.0, l3[i]), "node exact, linearT");
  }

  /* logInverseT is exact for Clausius-Clapeyron, inside and extrapolated */
  for (const double T : {510.0, 550.0, 599.0, 700.0, 790.0, 400.0, 1000.0}) {
    check(close(logT(T), std::pow(10.0, cc(T)), 1e-13),
          "logInverseT = Clausius-Clapeyron at " + std::to_string(T) + " K");
    check(close(pascal(T), std::pow(10.0, cc(T)), 1e-13),
          "logInverseT from Pa at " + std::to_string(T) + " K");
  }

  /* linearT: p linear in T; hold outside */
  check(close(bar(550), 0.5*(p3[0] + p3[1]), 1e-14), "linearT midpoint");
  check(close(bar(700), 0.5*(p3[1] + p3[2]), 1e-14), "linearT midpoint 2");
  check(bar(300) == p3[0] && bar(2000) == p3[2], "hold outside");
  check(throwsWith([&]() { fatal(499.9); }, "outside the table"),
        "fatal below the table");
  check(throwsWith([&]() { fatal(800.1); }, "outOfRange fatal"),
        "fatal above the table");

  /* linearT extrapolated: the straight line, not below 0 */
  const table linearOut(T3, {p3}, {"s"}, table::PA, table::LINEAR_T,
                        table::EXTRAPOLATE);
  const double slope = (p3[1] - p3[0])/100.0;
  check(close(linearOut(490), p3[0] - 10*slope, 1e-12),
        "linearT extrapolated below");
  check(linearOut(300) == 0, "linearT extrapolated, not below 0");

  /*--------------------------------------------------------------------------
  Two phases crossing at 650 K: log10 p_a = 12 - 4000/T, log10 p_b =
  12 - 4000/650 + 2000/650 - 2000/T (equal at 650 K; a below, b above).
  --------------------------------------------------------------------------*/
  {
    const std::vector<double> T = {500, 550, 600, 650, 700, 750};
    std::vector<double> a, b;
    for (const double t : T) {
      a.push_back(12 - 4000/t);
      b.push_back(12 - 4000.0/650 + 2000.0/650 - 2000/t);
    }
    const table two(T, {a, b}, {"a", "b"}, table::LOG10_PA,
                    table::LOG_INVERSE_T, table::EXTRAPOLATE);
    check(two.stablePhase(600) == 0 && two.stablePhase(700) == 1,
          "stable phase below and above the crossing");
    check(close(two(620), std::pow(10.0, 12 - 4000.0/620), 1e-13)
       && close(two(720), std::pow(10.0, 12 - 4000.0/650 + 2000.0/650
                                    - 2000.0/720), 1e-13),
          "p_eq = the minimum over the phases");
    const auto x = two.crossings();
    check(x.size() == 1 && std::abs(x[0].T - 650) < 1e-6
       && x[0].below == 0 && x[0].above == 1, "crossing at 650 K");
  }

  /*--------------------------------------------------------------------------
  Files and options
  --------------------------------------------------------------------------*/
  {
    const std::string good = writeFile("good.csv",
      "# a comment, with commas\n"
      "# T_K, s1 ,s2\n"
      "500,1.0,2.0\n"
      "600,2.0,1.5\r\n"
      "\n"
      "# a comment between the rows\n"
      "700,3.0,2.5\n");
    const table t = table::read(good, {"s2", "s1"}, table::LOG10_PA,
                                table::LOG_INVERSE_T, table::EXTRAPOLATE);
    check(t.nNodes() == 3 && t.Tmin() == 500 && t.Tmax() == 700,
          "read: nodes, CR LF and comments");
    check(t(600) == std::pow(10.0, 1.5) && t(500) == 10.0,
          "read: the phases selected by name, minimum");

    check(throwsWith([&]() { table::read(good, {"s3"}, table::LOG10_PA,
      table::LOG_INVERSE_T, table::EXTRAPOLATE); }, "no column 's3'"),
          "unknown column");
    check(throwsWith([&]() { table::read("/nonexistent.csv", {"s"},
      table::LOG10_PA, table::LOG_INVERSE_T, table::EXTRAPOLATE); },
      "cannot open"), "missing file");

    const std::string decreasing = writeFile("decreasing.csv",
      "# T_K,s\n500,1\n500,2\n");
    check(throwsWith([&]() { table::read(decreasing, {"s"},
      table::LOG10_PA, table::LOG_INVERSE_T, table::EXTRAPOLATE); },
      "strictly increasing"), "temperatures not increasing");

    const std::string short1 = writeFile("short.csv",
      "# T_K,s,u\n500,1,2\n600,1\n");
    check(throwsWith([&]() { table::read(short1, {"s"}, table::LOG10_PA,
      table::LOG_INVERSE_T, table::EXTRAPOLATE); }, "2 values for 3"),
          "row with a missing value");

    const std::string text = writeFile("text.csv",
      "# T_K,s\n500,1\n600,abc\n");
    check(throwsWith([&]() { table::read(text, {"s"}, table::LOG10_PA,
      table::LOG_INVERSE_T, table::EXTRAPOLATE); }, "not a finite number"),
          "value that is not a number");

    const std::string zero = writeFile("zero.csv",
      "# T_K,s\n500,0\n600,1\n");
    check(throwsWith([&]() { table::read(zero, {"s"}, table::PA,
      table::LOG_INVERSE_T, table::EXTRAPOLATE); }, "must be positive"),
          "p = 0 refused by logInverseT");
    check(table::read(zero, {"s"}, table::PA, table::LINEAR_T,
                      table::EXTRAPOLATE)(550) == 0.5,
          "p = 0 accepted by linearT");

    check(throwsWith([]() { table::units("mbar"); }, "unknown units"),
          "unknown units");
    check(throwsWith([]() { table::interpolation("spline"); },
                     "unknown interpolation"), "unknown interpolation");
    check(throwsWith([]() { table::outOfRange("clip"); },
                     "unknown outOfRange"), "unknown outOfRange");
    check(table::units("bar") == table::BAR
       && table::interpolation("linearT") == table::LINEAR_T
       && table::outOfRange("hold") == table::HOLD, "option words");
  }

  /*--------------------------------------------------------------------------
  The generated tables
  --------------------------------------------------------------------------*/
  if (argc == 4) {

    /* phase-resolved Gurvich table */
    const table gurvich = table::read(argv[1], {"PbI2(cr)", "PbI2(l)"},
      table::LOG10_PA, table::LOG_INVERSE_T, table::EXTRAPOLATE);
    const auto nodes = rows(argv[1]);
    bool exact = nodes.size() == gurvich.nNodes();
    for (const auto& r : nodes) {
      const double want = std::min(std::pow(10.0, r[1]), std::pow(10.0, r[2]));
      exact = exact && gurvich(r[0]) == want;
    }
    check(exact, "Gurvich table: every node exact");

    const auto x = gurvich.crossings();
    check(x.size() == 1 && std::abs(x[0].T - 683.0) <= 0.1
       && gurvich.phases()[x[0].below] == "PbI2(cr)"
       && gurvich.phases()[x[0].above] == "PbI2(l)",
          "Gurvich table: melting point 683.0 +- 0.1 K");

    const bool haveReference = std::string(argv[2]) != "-";
    double worst = 0, where = 0;
    std::size_t nReference = 0;
    if (haveReference) {
      for (const auto& r : rows(argv[2])) {
        const double err = std::abs(std::log10(gurvich(r[0])) - r[1]);
        if (err > worst) {
          worst = err;
          where = r[0];
        }
        ++nReference;
      }
      check(nReference > 3000 && worst < 2e-3,
            "Gurvich table vs NASA-9 below 2e-3 decades");
    }

    /* T-Flows table */
    const table tflows = table::read(argv[3], {"GEMS_VAPOR_PRESSURE"},
      table::BAR, table::LINEAR_T, table::HOLD);
    const auto tf = rows(argv[3]);
    bool reproduced = tf.size() == 95;
    double worstNode = 0;
    for (std::size_t k = 0; k < tf.size(); ++k) {
      const double T = 283.15 + 10*double(k);
      reproduced = reproduced && std::abs(tf[k][0] - T) < 1e-9;
      const double rel = std::abs(tflows(tf[k][0])/(tf[k][1]*1e5) - 1);
      worstNode = std::max(worstNode, rel);
    }
    check(reproduced && worstNode <= 1e-12,
          "T-Flows table: 95 nodes 283.15 + 10 k K reproduced to 1e-12");
    check(tflows(250) == tf.front()[1]*1e5
       && tflows(1300) == tf.back()[1]*1e5, "T-Flows table: ends held");
    check(close(tflows(703.15 + 5), 0.5e5*(tf[42][1] + tf[43][1]), 1e-14),
          "T-Flows table: p linear in T between nodes");

    std::printf("Gurvich table %s: %zu nodes %g-%g K, melting point "
                "(phase crossing) %.6f K; ", argv[1], gurvich.nNodes(),
                gurvich.Tmin(), gurvich.Tmax(), x.empty() ? 0.0 : x[0].T);
    if (haveReference) {
      std::printf("vs NASA-9 at %zu temperatures: max |dlog10 p| %.3e "
                  "decades at %.2f K\n", nReference, worst, where);
    } else {
      std::printf("vs NASA-9: not checked (no reference)\n");
    }
    std::printf("T-Flows table %s: %zu nodes, max relative difference at "
                "the nodes %.1e\n", argv[3], tf.size(), worstNode);
  } else {
    check(false, "usage: testVapourPressureTable <phases.csv> "
          "<reference.csv | -> <TFlows.csv>");
  }

  std::printf("%d checks, %d failed.\n", nChecks, nFailed);
  return nFailed == 0 ? 0 : 1;
}
