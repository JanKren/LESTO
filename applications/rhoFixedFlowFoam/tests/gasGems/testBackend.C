#include "../../gasSpeciation.H"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
using namespace LESTO;
int main(int argc,char** argv) {
  if(argc!=2)return 2;
  const std::vector<gasSpecies> f={{1,0,0},{1,0,1},{1,0,2},{0,1,0},{0,2,0},{0,1,1},{0,1,3},{0,0,1},{0,0,2}};
  gasSpeciation kernel(f);std::ifstream file(argv[1]);std::string line;double worst=0,balance=0;int count=0,failed=0;int worstRow=0,worstMode=0,worstSpecies=0;
  while(std::getline(file,line)){
    if(line.empty() || line[0]=='#')continue;
    std::replace(line.begin(),line.end(),',',' ');std::istringstream input(line);
    std::vector<double> r;double x;while(input>>x)r.push_back(x);
    if(r.size()!=26)return 2;
    std::array<double,3>b={r[2],r[3],r[4]},pi={r[23],r[24],r[25]};
    std::vector<double> k(9),c(9);for(int j=0;j<9;++j)k[j]=r[5+2*j];
    const double largest=std::max({b[0],b[1],b[2]});
    for(int mode=0;mode<2;++mode){
      if(mode==0)pi={NAN,NAN,NAN};else pi={r[23],r[24],r[25]};
      if(kernel.solve(k.data(),r[1],b,pi,c.data())<0)++failed;
      std::array<long double,3> totals={0,0,0};
      for(int j=0;j<9;++j){if(!std::isfinite(c[j]) || c[j]<0)++failed;
        if(r[6+2*j]>1e-10*largest && std::fabs(c[j]/r[6+2*j]-1)>worst){worst=std::fabs(c[j]/r[6+2*j]-1);worstRow=count;worstMode=mode;worstSpecies=j;}
        totals[0]+=f[j].nPb*static_cast<long double>(c[j]);totals[1]+=f[j].nBi*static_cast<long double>(c[j]);totals[2]+=f[j].nI*static_cast<long double>(c[j]);}
      for(int a=0;a<3;++a)if(b[a]>0)balance=std::max(balance,double(std::fabs(totals[a]/b[a]-1)));else if(totals[a]!=0)++failed;
    }++count;
  }
  std::cout<<count<<" points, both frozen/cold and local/warm polished: species error "<<worst<<"; element error "<<balance<<"; failures "<<failed<<'\n';
  std::cout<<"worst row "<<worstRow<<" mode "<<worstMode<<" species "<<worstSpecies<<"\n";
  return count!=936 || failed || worst>1e-9 || balance>1e-14;
}
