#include "../../gasSpeciation.H"
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <algorithm>
using namespace LESTO;
int main(int argc,char** argv) {
  if(argc<2 || argc>3)return 2;
  std::ifstream input(argv[1]); std::string line;
  std::getline(input,line);
  const bool external=line.find("Pb2(g)")!=std::string::npos;
  const std::vector<gasSpecies> formulas=external?std::vector<gasSpecies>{
    {1,0,0},{2,0,0},{1,0,1},{1,0,2},{1,0,3},{1,0,4},
    {0,0,1},{0,0,2},{0,1,0},{0,2,0},{0,1,1},{0,1,3}}:
    std::vector<gasSpecies>{{1,0,0},{1,0,1},{1,0,2},{0,1,0},
    {0,2,0},{0,1,1},{0,1,3},{0,0,1},{0,0,2}};
  gasSpeciation kernel(formulas);const int ns=kernel.size();
  std::vector<std::vector<double>> rows;
  while(std::getline(input,line)) {
    if(line.empty() || line[0]=='#')continue;
    std::replace(line.begin(),line.end(),',',' ');std::istringstream stream(line);
    std::vector<double> r;double x;while(stream>>x)r.push_back(x);
    if(r.size()!=size_t(5+2*ns))return 2;
    rows.push_back(r);
  }
  if(rows.empty())return 2;
  double worst=0,scaledDifference=0,balance=0,idempotence=0;int failures=0; size_t rowIndex=0,worstRow=0;int worstSpecies=0;
  std::vector<double> k(ns),c(ns),second(ns);
  auto run=[&](const auto& r,bool check) {
    std::array<double,3>b={r[2],r[3],r[4]},pi={NAN,NAN,NAN};
    for(int j=0;j<ns;++j)k[j]=r[5+2*j];
    if(kernel.solve(k.data(),r[1],b,pi,c.data())<0)++failures;
    if(!check)return;
    const double maximum=std::max({b[0],b[1],b[2]});
    std::array<long double,3> totals={0,0,0};
    for(int j=0;j<ns;++j) {
      if(!std::isfinite(c[j]) || c[j]<0)++failures;
      if(r[6+2*j]>1e-10*maximum && std::fabs(c[j]/r[6+2*j]-1)>worst) {
        worst=std::fabs(c[j]/r[6+2*j]-1);worstRow=rowIndex;worstSpecies=j;
      }
      if(maximum>0)scaledDifference=std::max(scaledDifference,std::fabs(c[j]-r[6+2*j])/maximum);
      const auto& f=formulas[j];totals[0]+=f.nPb*static_cast<long double>(c[j]);
      totals[1]+=f.nBi*static_cast<long double>(c[j]);totals[2]+=f.nI*static_cast<long double>(c[j]);
    }
    for(int e=0;e<3;++e) {
      if(b[e]>0)balance=std::max(balance,double(std::fabs(totals[e]/b[e]-1)));
      else if(totals[e]!=0)++failures;
    }
    std::array<double,3> newB={double(totals[0]),double(totals[1]),double(totals[2])};
    pi={NAN,NAN,NAN};kernel.solve(k.data(),r[1],newB,pi,second.data());
    for(int j=0;j<ns;++j)if(maximum>0)idempotence=std::max(idempotence,std::fabs(second[j]-c[j])/maximum);
  };
  for(const auto&r:rows){run(r,true);++rowIndex;}
  const auto start=std::chrono::steady_clock::now();
  const int repetitions=rows.size()>10000 ? 1 : 20;
  for(int rep=0;rep<repetitions;++rep)for(const auto&r:rows)run(r,false);
  const double cost=1e6*std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/(repetitions*rows.size());
  std::cout<<rows.size()<<" points; species error "<<worst<<"; largest-element scaled difference "<<scaledDifference<<"; balance "<<balance
           <<"; idempotence (largest-element scale) "<<idempotence<<"; failures "<<failures
           <<"; cold us/call "<<cost<<'\n';
  std::cout<<"Worst row "<<worstRow<<" species index "<<worstSpecies<<"\n";
  // Accuracy gates are deterministic. Timing is reported separately because
  // the machine may have other running solver jobs.
  return failures || (argc==3 ? scaledDifference>1e-12 : worst>1e-9) || balance>1e-14 || idempotence>1e-15;
}
