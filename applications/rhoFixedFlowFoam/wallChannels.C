#include "wallChannels.H"
#include "phaseChangeKinetics.H"
#include "calculatedFvPatchFields.H"
#include "boundBox.H"
#include <fstream>
#include <sstream>
#include <cstdint>
using namespace Foam;
namespace {
scalar global(const LESTO::compensatedSum& sum) {
  scalar v=sum.value();reduce(v,sumOp<scalar>());return v;
}
const char* atoms[]={"Pb","Bi","I"};
scalar atomsOf(const LESTO::chemicalFormula& f,int e) {
  auto a=f.find(atoms[e]);return a==f.end()?0:a->second;
}
void keys(const dictionary& cfg,const wordList& allowed) {
  for(const entry& key:cfg)if(!allowed.found(key.keyword()))
    FatalIOErrorInFunction(cfg)<<"Unknown wall channel/source key "<<key.keyword()<<exit(FatalIOError);
}
}
LESTO::wallChannels::wallChannels(const fvMesh& mesh,const fluidThermo& thermo,const volScalarField& rho,
 PtrList<volScalarField>& fields,const wordList& names,const List<word>& state,
 const List<scalar>& masses,const dictionary& dict,interfaceExchange& exchange)
: mesh_(mesh),thermo_(thermo),rho_(rho),fields_(fields),exchange_(exchange),masses_(masses),
 formulas_(names.size()),flows_(names.size()),release_(names.size()),pending_(names.size()),su_(names.size()),sp_(names.size()) {
  if(dict.get<word>("model")!="HKS")FatalIOErrorInFunction(dict)<<"Wall channels require model HKS"<<exit(FatalIOError);
  IOdictionary props(IOobject("speciesTransportProperties",mesh.time().constant(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE,false));
  const dictionary& all=dict.subDict("wallChannels");
  if(all.empty())FatalIOErrorInFunction(all)<<"Empty wallChannels"<<exit(FatalIOError);
  OStringStream config;config.precision(17);config<<"ledgerVersion 2; "<<all;
  if(dict.isDict("gasSources"))config<<dict.subDict("gasSources");
  std::ostringstream signature;signature<<config.str().c_str();
  elementsOfCell_.resize(mesh.nCells());
  for(label e=0;e<exchange.nWallElements();++e)elementsOfCell_[exchange.wallCells()[e]].push_back(e);
  forAll(names,si) {
    formulas_[si]=readFormula(props.subDict(names[si]),masses[si]);
    signature<<names[si]<<' '<<phaseChangeLedger::exactDecimal(masses[si]);
    for(const auto& atom:formulas_[si])signature<<atom.first<<' '<<atom.second;
    flows_[si].resize(mesh.nCells(),0);release_[si].resize(mesh.nCells(),0);pending_[si].resize(mesh.nCells(),0);
    su_.set(si,new volScalarField::Internal(IOobject("wallSu_"+names[si],mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),mesh,dimensionedScalar(dimDensity/dimTime,Zero)));
    sp_.set(si,new volScalarField::Internal(IOobject("wallSp_"+names[si],mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE,false),mesh,dimensionedScalar(dimDensity/dimTime,Zero)));
  }
  auto species=[&](const word& name,const word& kind) {
    const label si=names.find(name);
    if(si<0 || state[si]!=kind || formulas_[si].empty())
      FatalIOErrorInFunction(all)<<"Channel species "<<name<<" needs state "<<kind<<", formula and molarMass"<<exit(FatalIOError);
    return si;
  };
  auto readTable=[&](const dictionary& cfg,const wordList& columns,const word& units) {
    keys(cfg,wordList({"file","column","phases","units","interpolation","outOfRange"}));
    fileName file(cfg.get<fileName>("file"));file.replaceAll("<constant>",mesh.time().globalPath()/mesh.time().constant());file.expand();
    std::ifstream input(file.c_str());
    if(!input)FatalIOErrorInFunction(cfg)<<"Cannot read "<<file<<exit(FatalIOError);
    signature<<input.rdbuf();std::vector<std::string> col;for(const word& c:columns)col.push_back(c);
    try {
      return new vapourPressureTable(vapourPressureTable::read(file,col,vapourPressureTable::units(units),
        vapourPressureTable::interpolation(cfg.getOrDefault<word>("interpolation","logInverseT")),
        vapourPressureTable::outOfRange(cfg.getOrDefault<word>("outOfRange","fatal"))));
    } catch(const std::exception& e) {
      FatalIOErrorInFunction(cfg)<<e.what()<<exit(FatalIOError);return static_cast<vapourPressureTable*>(nullptr);
    }
  };
  for(const Foam::entry& entry:all) {
    const dictionary& cfg=entry.dict();
    keys(cfg,wordList({"reactant","condensed","stoichiometry","vapourPressure","equilibrium","reversible","HKS"}));
    const label r=species(cfg.get<word>("reactant"),"gas"),s=species(cfg.get<word>("condensed"),"solid");
    if(exchange.pairOfGas(r)>=0)FatalIOErrorInFunction(cfg)<<"Channel reactant already belongs to a legacy pair"<<exit(FatalIOError);
    for(const Foam::entry& pair:dict.subDict("pairs"))if(pair.dict().get<word>("condensed")==names[s])
      FatalIOErrorInFunction(cfg)<<"A channel condensate cannot also belong to a legacy pair"<<exit(FatalIOError);
    label storei=-1;forAll(stores_,i)if(stores_[i].species==s)storei=i;
    if(storei<0) {
      storei=stores_.size();stores_.resize(storei+1);stores_.set(storei,new store);
      auto& d=stores_[storei];d.species=s;d.mass=masses[s];d.formula=formulas_[s];d.deposit.resize(exchange.nWallElements(),0);
    }
    for(const channel& c:channels_)if(c.reactant==r && c.condensate==storei)
      FatalIOErrorInFunction(cfg)<<"Duplicate reactant/condensate channel"<<exit(FatalIOError);
    const label ci=channels_.size();channels_.resize(ci+1);channels_.set(ci,new channel);
    auto& c=channels_[ci];c.name=entry.keyword();c.reactant=r;c.condensate=storei;
    if(cfg.isDict("stoichiometry")) {
      c.reactive=true;c.irreversible=!cfg.getOrDefault<bool>("reversible",false);
      if(!c.irreversible)FatalIOErrorInFunction(cfg)<<"Reactive wall channels require reversible no"<<exit(FatalIOError);
      const dictionary& nu=cfg.subDict("stoichiometry");const scalar nr=nu.get<scalar>(names[r]);
      if(!std::isfinite(nr) || nr>=0)FatalIOErrorInFunction(nu)<<"Reactant stoichiometry must be finite and negative"<<exit(FatalIOError);
      c.reactantNu=-nr;c.nuCondensate=nu.get<scalar>(names[s])/c.reactantNu;
      chemicalFormula residual;
      for(const Foam::entry& term:nu) {
        const scalar coefficient=nu.get<scalar>(term.keyword());const label si=names.find(term.keyword());
        if(si<0 || formulas_[si].empty() || !std::isfinite(coefficient) || coefficient==0 || (si!=r && coefficient<0))
          FatalIOErrorInFunction(nu)<<"Invalid channel stoichiometry"<<exit(FatalIOError);
        if(si!=r && si!=s) {
          if(state[si]!="gas")FatalIOErrorInFunction(nu)<<"Co-products must be gases"<<exit(FatalIOError);
          c.products.emplace_back(si,coefficient/c.reactantNu);
        }
        for(const auto& atom:formulas_[si])residual[atom.first]+=coefficient*atom.second;
      }
      for(const auto& atom:residual)if(std::fabs(atom.second)>1e-12)
        FatalIOErrorInFunction(nu)<<"Channel does not conserve "<<atom.first<<exit(FatalIOError);
      const dictionary& eq=cfg.subDict("equilibrium");
      if(eq.found("units") || eq.found("phases"))FatalIOErrorInFunction(eq)<<"Reaction table requires dimensionless log10 K (p0=1 bar), with column, not units/phases"<<exit(FatalIOError);
      c.table.reset(readTable(eq,wordList({eq.getOrDefault<word>("column","logK")}),"log10Pa"));
    } else {
      c.irreversible=!cfg.getOrDefault<bool>("reversible",true);c.nuCondensate=cfg.getOrDefault<scalar>("stoichiometry",1);
      for(const auto& atom:formulas_[r])if(!formulas_[s].count(atom.first) || std::fabs(atom.second-c.nuCondensate*formulas_[s].at(atom.first))>1e-12)
        FatalIOErrorInFunction(cfg)<<"Channel formulas do not match stoichiometry"<<exit(FatalIOError);
      if(formulas_[r].size()!=formulas_[s].size())FatalIOErrorInFunction(cfg)<<"Channel formulas differ"<<exit(FatalIOError);
      const dictionary& vp=cfg.subDict("vapourPressure");c.table.reset(readTable(vp,vp.get<wordList>("phases"),vp.getOrDefault<word>("units","log10Pa")));
    }
    if(!std::isfinite(c.nuCondensate) || c.nuCondensate<=0)FatalIOErrorInFunction(cfg)<<"Condensate stoichiometry must be positive"<<exit(FatalIOError);
    const dictionary hks=cfg.subOrEmptyDict("HKS");keys(hks,wordList({"accommodation","Ce","kineticScale"}));
    c.accommodation=hks.getOrDefault<scalar>("accommodation",1);c.Ce=hks.getOrDefault<scalar>("Ce",1);c.scale=hks.getOrDefault<scalar>("kineticScale",1);
    if(!std::isfinite(c.accommodation) || c.accommodation<=0 || c.accommodation>1 || !std::isfinite(c.Ce) || c.Ce<0 || !std::isfinite(c.scale) || c.scale<0)
      FatalIOErrorInFunction(hks)<<"Invalid channel HKS coefficients"<<exit(FatalIOError);
    const label n=exchange.nWallElements();c.a.resize(n);c.g.resize(n);c.upper.resize(n);c.flow.resize(n);c.pEq.resize(n);c.regime.resize(n,label(NONE));
    mesh.setFluxRequired(fields[r].name());
  }
  if(dict.isDict("gasSources"))for(const Foam::entry& entry:dict.subDict("gasSources")) {
    const dictionary& cfg=entry.dict();keys(cfg,wordList({"amounts","box","startTime","duration"}));
    source src;src.start=cfg.get<scalar>("startTime");src.duration=cfg.get<scalar>("duration");src.length=releaseWindowLength(src.start,src.duration);
    if(!std::isfinite(src.start) || !std::isfinite(src.duration) || src.duration<=0 || !std::isfinite(src.length) || src.length<=0)
      FatalIOErrorInFunction(cfg)<<"Invalid gas source window"<<exit(FatalIOError);
    const List<vector> bounds(cfg.get<List<vector>>("box"));
    if(bounds.size()!=2)FatalIOErrorInFunction(cfg)<<"Source box needs two corners"<<exit(FatalIOError);
    boundBox box(bounds[0],bounds[1]);DynamicList<label> cells;compensatedSum volume;
    forAll(mesh.C(),cell)if(box.contains(mesh.C()[cell])) {cells.append(cell);volume.add(mesh.V()[cell]);}
    src.cells.transfer(cells);src.volume=global(volume);
    if(!(src.volume>0))FatalIOErrorInFunction(cfg)<<"Empty gas source box"<<exit(FatalIOError);
    for(const Foam::entry& amount:cfg.subDict("amounts")) {
      const label si=species(amount.keyword(),"gas");const scalar n=cfg.subDict("amounts").get<scalar>(amount.keyword());
      if(!std::isfinite(n) || n<0 || exchange.pairOfGas(si)>=0)
        FatalIOErrorInFunction(cfg)<<"Q2 source amounts must be finite, nonnegative and unpaired"<<exit(FatalIOError);
      src.amounts.emplace_back(si,n);
    }
    sources_.push_back(src);
  }
  std::uint64_t hash=14695981039346656037ULL;
  for(unsigned char b:signature.str()) {hash^=b;hash*=1099511628211ULL;}
  const string key(std::to_string(hash));dictionary& layout=exchange.ledger().layout();const bool restarted=exchange.ledger().restarted();
  if(restarted && layout.getOrDefault<string>("wallChannelsHash",string::null)!=key)
    FatalIOErrorInFunction(all)<<"Wall channel configuration/data changed across restart"<<exit(FatalIOError);
  layout.set("wallChannelsHash",key);
  // Per-rank geometry controls exact restore; mw_ supports changed decomposition.
  hash=14695981039346656037ULL;std::ostringstream geometry;
  for(label e=0;e<exchange.nWallElements();++e) {
    const label cell=exchange.wallCells()[e];geometry<<cell;
    for(int a=0;a<3;++a)geometry<<' '<<phaseChangeLedger::exactDecimal(mesh.C()[cell][a]);
    geometry<<' '<<phaseChangeLedger::exactDecimal(exchange.wallAreas()[e])<<' '<<exchange.wallPatches()[e]<<' '<<exchange.wallFaces()[e];
  }
  for(unsigned char b:geometry.str()) {hash^=b;hash*=1099511628211ULL;}
  fingerprint_=scalar(hash&((1ULL<<52)-1));
  deposits_.reset(new binaryIOList<scalar>(IOobject("phaseChangeWallChannels",mesh.time().timeName(),"uniform",mesh.time(),IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE)));
  accounts_.reset(new binaryGlobalIOList<scalar>(IOobject("phaseChangeChannelCondensates",mesh.time().timeName(),"uniform",mesh.time(),IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE)));
  if(restarted && accounts_->size()!=10*stores_.size())FatalIOErrorInFunction(all)<<"Missing channel condensate accounts"<<exit(FatalIOError);
  if(!restarted) {
    if(!accounts_->empty() || !deposits_->empty())FatalIOErrorInFunction(all)<<"Orphan channel restart state"<<exit(FatalIOError);
    accounts_->resize(10*stores_.size(),0);
  }
  for(scalar value:accounts_())if(!std::isfinite(value))
    FatalIOErrorInFunction(all)<<"Nonfinite channel condensate account"<<exit(FatalIOError);
  forAll(stores_,i)if((*accounts_)[10*i]<0 || (*accounts_)[10*i+1]<0 || (*accounts_)[10*i+8]+(*accounts_)[10*i+9]<0)
    FatalIOErrorInFunction(all)<<"Negative channel inventory/exchange reference"<<exit(FatalIOError);
  const label nWall=exchange.nWallElements();const bool exact=deposits_->size()==1+nWall*stores_.size() && (*deposits_)[0]==fingerprint_;
  forAll(stores_,i) {
    auto& d=stores_[i];IOobject io("mw_"+names[d.species],mesh.time().timeName(),mesh,IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE);
    const bool present=io.typeHeaderOk<volScalarField>(true,false);
    if(restarted && !exact && !present)FatalIOErrorInFunction(all)<<"Restart needs exact channel state or "<<io.name()<<exit(FatalIOError);
    if(present)d.mw.reset(new volScalarField(io,mesh));
    else {io.readOpt(IOobject::NO_READ);d.mw.reset(new volScalarField(io,mesh,dimensionedScalar(dimMass/dimArea,Zero),calculatedFvPatchScalarField::typeName));}
    if(d.mw().dimensions()!=dimMass/dimArea)FatalIOErrorInFunction(all)<<"Channel deposits require kg/m2"<<exit(FatalIOError);
    d.mDep.reset(new volScalarField(IOobject("mDep_"+names[d.species],mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,dimensionedScalar(dimMass/dimArea,Zero),calculatedFvPatchScalarField::typeName));
    for(label e=0;e<nWall;++e) {
      d.deposit[e]=restarted && exact?(*deposits_)[1+i*nWall+e]:exchange.distanceLayer()?d.mw()[exchange.wallCells()[e]]:d.mw().boundaryField()[exchange.wallPatches()[e]][exchange.wallFaces()[e]];
      if(!std::isfinite(d.deposit[e]) || d.deposit[e]<0)FatalIOErrorInFunction(all)<<"Invalid channel deposit"<<exit(FatalIOError);
    }
    const scalar held=storeMass(i);
    if(!restarted)(*accounts_)[10*i]=held;
    else compensatedAdd((*accounts_)[10*i+2],(*accounts_)[10*i+3],held-(*accounts_)[10*i+1]);
    (*accounts_)[10*i+1]=held;
    for(const auto& atom:d.formula)exchange.includeElement(atom.first);
  }
  if(Pstream::master()) {
    const fileName dir(mesh.time().globalPath()/"postProcessing/phaseChangeChannelCondensates"/mesh.time().timeName());mkDir(dir);files_.resize(stores_.size());
    forAll(stores_,i) {
      files_.set(i,new OFstream(dir/(names[stores_[i].species]+".dat")));files_[i].precision(17);
      files_[i]<<"# columns time initial exchanged clamped restart held supplied closure reference"<<nl;
    }
  }
  finish();
  Info<<"Wall channels: "<<channels_.size()<<" channels, "<<stores_.size()<<" shared condensates; co-products after transport; bounded species order"<<nl;
}
bool LESTO::wallChannels::owns(label si)const {
  for(const auto& c:channels_)if(c.reactant==si)return true;
  return false;
}
bool LESTO::wallChannels::uses(label si)const {
  if(owns(si) || hasSource(si))return true;
  for(const auto& c:channels_)for(const auto& p:c.products)if(p.first==si)return true;
  return false;
}
bool LESTO::wallChannels::hasSource(label si)const {
  for(const auto& s:sources_)for(const auto& n:s.amounts)if(n.first==si)return true;
  return false;
}
scalar LESTO::wallChannels::storeMass(label i)const {
  compensatedSum sum;forAll(stores_[i].deposit,e)sum.add(stores_[i].deposit[e]*exchange_.wallAreas()[e]);return global(sum);
}
void LESTO::wallChannels::conversion(label si,scalar mass) {
  for(int e=0;e<3;++e)conversion_[e].add(atomsOf(formulas_[si],e)/masses_[si]*mass);
}
void LESTO::wallChannels::beginStep() {
  for(auto& c:conversion_)c=compensatedSum();
  forAll(fields_,si) {flows_[si]=0;release_[si]=0;pending_[si]=0;}
  const scalar dt=mesh_.time().deltaTValue(),t0=exchange_.ledger().time();
  for(const auto& s:sources_) {
    const scalar overlap=releaseOverlap(t0,t0+dt,s.start,s.duration);
    for(const auto& n:s.amounts) {
      const scalar rate=n.second*masses_[n.first]*overlap/(s.length*dt*s.volume);
      for(label cell:s.cells)release_[n.first][cell]+=rate*mesh_.V()[cell];
    }
  }
  for(auto& c:channels_) {
    const scalarField resistance=exchange_.wallDiffusionResistance(c.reactant);
    forAll(c.a,e) {
    const label cell=exchange_.wallCells()[e];const scalar T=exchange_.wallTemperatures()[e];
    double pressure=0;
    try {pressure=c.table()(T);}catch(const std::exception& error) {FatalErrorInFunction<<c.name<<": "<<error.what()<<exit(FatalError);}
    if(!(pressure>0) || !std::isfinite(pressure))FatalErrorInFunction<<c.name<<": invalid equilibrium table value"<<exit(FatalError);
    if(c.reactive) {
      long double lg=-std::log(pressure);
      for(const auto& p:c.products) {
        const double partial=std::max(double(rho_[cell]*fields_[p.first][cell]/masses_[p.first]*8.314462618*thermo_.T()[cell]),0.0);
        if(partial==0) {lg=-HUGE_VAL;break;}
        lg+=p.second*c.reactantNu*std::log(partial/1e5);
      }
      pressure=1e5*std::exp(double(lg/c.reactantNu));
    }
    c.pEq[e]=pressure;
    scalar G=hksConductance(masses_[c.reactant],T,c.accommodation,c.Ce,c.scale);
    const scalar beta=rho_[cell]*8.314462618*thermo_.T()[cell]/masses_[c.reactant];
    G=resistance[e]>=0?halfCellConductance(G,beta,resistance[e]):0;
    c.a[e]=exchange_.wallAreas()[e]*G*pressure;
    c.g[e]=exchange_.wallAreas()[e]*G*rho_[cell]*8.314462618*thermo_.T()[cell]/masses_[c.reactant];
  }
  }
}
void LESTO::wallChannels::prepare(label si) {
  for(auto& c:channels_)if(c.reactant==si)forAll(c.upper,e)
    c.upper[e]=c.irreversible?0:stores_[c.condensate].deposit[e]*exchange_.wallAreas()[e]/mesh_.time().deltaTValue()*masses_[si]/(c.nuCondensate*stores_[c.condensate].mass);
}
void LESTO::wallChannels::assemble(label si) {
  su_[si].field()=0;sp_[si].field()=0;
  forAll(mesh_.V(),cell)su_[si][cell]=release_[si][cell]/mesh_.V()[cell];
  for(const auto& c:channels_)if(c.reactant==si)forAll(c.regime,e) {
    const label cell=exchange_.wallCells()[e];
    if(c.regime[e]==IMPLICIT) {su_[si][cell]+=c.a[e]/mesh_.V()[cell];sp_[si][cell]+=c.g[e]/mesh_.V()[cell];}
    else if(c.regime[e]==CAPPED)su_[si][cell]+=c.upper[e]/mesh_.V()[cell];
  }
}
void LESTO::wallChannels::linearise(label si,const fvScalarMatrix& eqn) {
  const tmp<volScalarField> H=eqn.H(),A=eqn.A();
  forAll(mesh_.V(),cell) {
    if(elementsOfCell_[cell].empty())continue;
    std::vector<double>a,g,lower,upper;std::vector<int>regime,previous;std::vector<std::pair<label,label>>events;
    forAll(channels_,i) {
      auto& c=channels_[i];if(c.reactant!=si)continue;
      for(label e:elementsOfCell_[cell]) {
        a.push_back(c.a[e]);g.push_back(c.g[e]);lower.push_back(-VGREAT);upper.push_back(c.upper[e]);previous.push_back(c.regime[e]);events.emplace_back(i,e);
      }
    }
    if(events.empty())continue;
    regime.resize(events.size());
    cellPredictor(A()[cell]*mesh_.V()[cell],H()[cell]*mesh_.V()[cell]+release_[si][cell],int(events.size()),a.data(),g.data(),lower.data(),upper.data(),regime.data(),previous.data());
    for(std::size_t i=0;i<events.size();++i)channels_[events[i].first].regime[events[i].second]=regime[i];
  }
  assemble(si);
}
bool LESTO::wallChannels::guard(label si) {
  bool changed=false;
  for(auto& c:channels_)if(c.reactant==si)forAll(c.regime,e) {
    const scalar q=c.a[e]-c.g[e]*fields_[si][exchange_.wallCells()[e]];const int r=boundedRegime(q,-VGREAT,c.upper[e]);
    if(r!=c.regime[e]) {c.regime[e]=r;changed=true;}
  }
  assemble(si);return returnReduce(changed,orOp<bool>());
}
void LESTO::wallChannels::realise(label si) {
  const scalar dt=mesh_.time().deltaTValue();label capped=0;compensatedSum clamped;
  std::vector<compensatedSum> storeClamped(stores_.size()),transfers(stores_.size()),variation(stores_.size());flows_[si]=0;
  for(auto& c:channels_)if(c.reactant==si)forAll(c.flow,e) {
    const label cell=exchange_.wallCells()[e];
    const scalar original=elementFlow(c.regime[e],c.a[e],c.g[e],c.upper[e],fields_[si][cell]);
    const scalar q=admissibleFlow(original,-VGREAT,c.upper[e]);
    if(q!=original)fields_[si][cell]+=(q-original)*dt/(rho_[cell]*mesh_.V()[cell]);
    c.flow[e]=q;flows_[si][cell]+=q;
    auto& d=stores_[c.condensate];const scalar held=d.deposit[e]*exchange_.wallAreas()[e];
    const scalar transfer=q*dt*c.nuCondensate*d.mass/masses_[si];
    if(c.regime[e]==CAPPED) {d.deposit[e]=0;++capped;clamped.add(held-transfer);}
    else d.deposit[e]-=transfer/exchange_.wallAreas()[e];
    if(d.deposit[e]<0 || !std::isfinite(d.deposit[e]))FatalErrorInFunction<<"Channel produced invalid condensate"<<exit(FatalError);
    const scalar residue=(held-d.deposit[e]*exchange_.wallAreas()[e])-transfer;
    storeClamped[c.condensate].add(residue);transfers[c.condensate].add(transfer);variation[c.condensate].add(mag(transfer));
    conversion(d.species,d.deposit[e]*exchange_.wallAreas()[e]-held+residue);conversion(si,q*dt);
    for(const auto& p:c.products)pending_[p.first][cell]+=-q*dt*p.second*masses_[p.first]/masses_[si];
  }
  forAll(stores_,i) {
    compensatedAdd((*accounts_)[10*i+4],(*accounts_)[10*i+5],global(storeClamped[i]));
    compensatedAdd((*accounts_)[10*i+6],(*accounts_)[10*i+7],global(transfers[i]));
    compensatedAdd((*accounts_)[10*i+8],(*accounts_)[10*i+9],global(variation[i]));
  }
  fields_[si].correctBoundaryConditions();reduce(capped,sumOp<label>());
  Info<<fields_[si].name()<<" wall channels: exchange "<<totalRate(si)*dt<<" kg; CAPPED "<<capped<<"; clamped "<<global(clamped)<<" kg"<<nl;
}
void LESTO::wallChannels::products(scalarList& amounts) {
  amounts.resize(fields_.size(),0);
  forAll(fields_,si) {
    compensatedSum transferred;
    forAll(pending_[si],cell)if(pending_[si][cell]!=0) {
      const scalar before=rho_[cell]*fields_[si][cell]*mesh_.V()[cell];
      fields_[si][cell]+=pending_[si][cell]/(rho_[cell]*mesh_.V()[cell]);
      const scalar actual=rho_[cell]*fields_[si][cell]*mesh_.V()[cell]-before;
      transferred.add(actual);conversion(si,actual);
    }
    amounts[si]=global(transferred);if(amounts[si]!=0)fields_[si].correctBoundaryConditions();
  }
}
scalar LESTO::wallChannels::totalRate(label si)const {
  compensatedSum sum;for(scalar q:flows_[si])sum.add(q);return global(sum);
}
scalar LESTO::wallChannels::totalRelease(label si)const {
  compensatedSum sum;for(scalar q:release_[si])sum.add(q);return global(sum)*mesh_.time().deltaTValue();
}
scalar LESTO::wallChannels::exchangeReference(label si)const {
  compensatedSum sum;
  for(const auto& c:channels_)if(c.reactant==si) {
    const label i=c.condensate;
    sum.add(((*accounts_)[10*i+8]+(*accounts_)[10*i+9])*masses_[si]/(c.nuCondensate*stores_[i].mass));
  }
  return sum.value();
}
scalar LESTO::wallChannels::elementConversion(label e)const {return global(conversion_[e]);}
scalar LESTO::wallChannels::elementStore(label e,int entry)const {
  compensatedSum sum;
  forAll(stores_,i) {
    const scalar mass=entry==0?(*accounts_)[10*i]:entry==1?(*accounts_)[10*i+1]:entry==2?(*accounts_)[10*i+2]+(*accounts_)[10*i+3]:(*accounts_)[10*i+4]+(*accounts_)[10*i+5];
    sum.add(atomsOf(stores_[i].formula,e)/stores_[i].mass*mass);
  }
  return sum.value();
}
void LESTO::wallChannels::finish() {
  const label n=exchange_.nWallElements();deposits_->resize(1+stores_.size()*n);(*deposits_)[0]=fingerprint_;
  forAll(stores_,i) {
    auto& d=stores_[i];(*accounts_)[10*i+1]=storeMass(i);auto& Ys=fields_[d.species];Ys.primitiveFieldRef()=0;
    forAll(d.deposit,e) {
      (*deposits_)[1+i*n+e]=d.deposit[e];const label cell=exchange_.wallCells()[e];
      Ys[cell]+=d.deposit[e]*exchange_.wallAreas()[e]/(rho_[cell]*mesh_.V()[cell]);
    }
    Ys.correctBoundaryConditions();exchange_.writeWallStore(d.mw(),d.mDep(),d.deposit);
    if(Pstream::master() && !files_.empty()) {
      const scalar initial=(*accounts_)[10*i],held=(*accounts_)[10*i+1],restart=(*accounts_)[10*i+2]+(*accounts_)[10*i+3],clamped=(*accounts_)[10*i+4]+(*accounts_)[10*i+5],exchanged=(*accounts_)[10*i+6]+(*accounts_)[10*i+7];
      compensatedSum supplied;supplied.add(initial);supplied.add(-exchanged);supplied.add(-clamped);supplied.add(restart);
      files_[i]<<exchange_.ledger().time()<<' '<<initial<<' '<<exchanged<<' '<<clamped<<' '<<restart<<' '<<held<<' '<<supplied.value()<<' '<<held-supplied.value()<<' '<<max(max(mag(initial)+(*accounts_)[10*i+8]+(*accounts_)[10*i+9],held),VSMALL)<<nl;files_[i].flush();
    }
  }
}
void LESTO::wallChannels::addProfiles(const word& element,scalarField& wall,scalar& inventory)const {
  forAll(stores_,i) {
    const auto& d=stores_[i];const auto atom=d.formula.find(element);if(atom==d.formula.end())continue;const scalar k=atom->second/d.mass;
    forAll(wall,e)wall[e]+=k*d.deposit[e]*exchange_.wallAreas()[e];
    inventory+=k*(*accounts_)[10*i+1];
  }
}
void LESTO::wallChannels::writeProfiles(const word& time,scalar clock,const string& reason) {
  forAll(stores_,i) {
    const auto& d=stores_[i];scalarField wall(d.deposit.size());
    forAll(wall,e)wall[e]=d.deposit[e]*exchange_.wallAreas()[e];
    exchange_.writeChannelProfile(time,clock,reason,d.mw().name().substr(3),d.mass,wall,(*accounts_)[10*i+1]);
  }
}
