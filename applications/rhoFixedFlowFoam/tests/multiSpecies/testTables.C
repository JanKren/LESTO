#include "../../vapourPressureTable.H"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
int main(int argc, char** argv) {
  assert(argc == 2);
  const std::string gases[] = {"BiI3", "Bi", "Bi2", "Pb"};
  const double expected[] = {680.85, 544.51, 544.51, 600.65};
  for (int g=0; g<4; ++g) {
    const std::string phase = g == 2 ? "Bi" : gases[g];
    const std::string file = std::string(argv[1])+"/pv_"+gases[g]+"_phases.csv";
    using V = LESTO::vapourPressureTable;
    const V table=V::read(file,{phase+"(cr)",phase+"(l)"},
      V::LOG10_PA,V::LOG_INVERSE_T,V::FATAL);
    std::ifstream input(file);
    std::string line;
    std::size_t rows=0;
    while (std::getline(input,line)) {
      if (line.empty() || line[0]=='#') continue;
      std::replace(line.begin(),line.end(),',',' ');
      std::istringstream row(line);
      double T,a,b; row>>T>>a>>b; assert(row);
      assert(table.nodes()[rows++]==T);
      assert(table.phase(0,T)==std::pow(10.,a));
      assert(table.phase(1,T)==std::pow(10.,b));
      assert(table(T)==std::min(std::pow(10.,a),std::pow(10.,b)));
    }
    assert(rows==table.nodes().size());
    double lo=500,hi=750;
    for (int n=0;n<70;++n) {
      const double T=(lo+hi)/2;
      if (table.phase(0,T)<table.phase(1,T)) lo=T; else hi=T;
    }
    assert(std::abs((lo+hi)/2-expected[g])<.1);
    std::cout<<gases[g]<<": all "<<rows<<" table nodes exact, crossing "<<(lo+hi)/2<<" K\n";
  }
}
