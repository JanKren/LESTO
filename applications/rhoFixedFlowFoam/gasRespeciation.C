#include "gasRespeciation.H"
#include "vapourPressureTable.H"
#include "matrixDefect.H"
#include "fvc.H"
#include "fvm.H"
#include "emptyFvPatch.H"
#include "OSspecific.H"
#include <fstream>
#include <sstream>
#include <cstdint>
#include <map>
#include <chrono>
using namespace Foam;
namespace {
scalar globalAmount(const LESTO::compensatedSum& local) {
  List<FixedList<scalar,2>> parts(Pstream::nProcs());
  parts[Pstream::myProcNo()] = {local.sum(), local.compensation()};
  Pstream::gatherList(parts);
  LESTO::compensatedSum sum;
  if (Pstream::master()) for (const auto& p : parts) { sum.add(p[0]); sum.add(p[1]); }
  scalar result = sum.value(); Pstream::broadcast(result); return result;
}
}
scalar LESTO::gasRespeciation::value(label i, entry e) const {
  const label k = offset(i,e), n = ids_.size()*NENTRY;
  return (*values_)[k]+(*values_)[n+k];
}
void LESTO::gasRespeciation::add(label i, entry e, scalar x) {
  const label k = offset(i,e), n = ids_.size()*NENTRY;
  compensatedAdd((*values_)[k], (*values_)[n+k], x);
}
void LESTO::gasRespeciation::set(label i, entry e, scalar x) {
  const label k = offset(i,e), n = ids_.size()*NENTRY;
  (*values_)[k]=x; (*values_)[n+k]=0;
}
scalar LESTO::gasRespeciation::inventory(label i) const {
  compensatedSum sum;
  const auto& Y = fields_[ids_[i]].primitiveField();
  forAll(Y,c) sum.add(rho_[c]*Y[c]*mesh_.V()[c]);
  return globalAmount(sum);
}
LESTO::gasRespeciation::gasRespeciation(
  const fvMesh& mesh, const fluidThermo& thermo, const volScalarField& rho,
  PtrList<volScalarField>& fields, const wordList& names, const List<word>& state,
  const List<scalar>& masses, interfaceExchange& exchange)
: mesh_(mesh), thermo_(thermo), rho_(rho), fields_(fields), exchange_(exchange),
  index_(names.size(),-1)
{
  const IOdictionary dict(IOobject("thermochemistryProperties", mesh.time().constant(),
    mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE, false));
  const dictionary* rd = dict.findDict("respeciation");
  const word engine = rd ? rd->getOrDefault<word>("engine","none") : word("none");
  if (engine == "none") {
    if(dict.found("wallChannels") || dict.found("gasSources"))FatalIOErrorInFunction(dict)<<"Wall channels/Q2 sources require respeciation accounts (engine kernel or GEMS)"<<exit(FatalIOError);
    if(!exchange.mock() && exchange.ledger().layout().found("wallChannelsHash"))FatalIOErrorInFunction(dict)<<"Wall channels removed across restart"<<exit(FatalIOError);
    if (!exchange.mock() && exchange.ledger().layout().found("respeciation"))
      FatalErrorInFunction << "Re-speciation engine changed across restart" << exit(FatalError);
    return;
  }
  const wordList known({"engine","species","formationConstants","minTemperature","updateInterval","GEMS"});
  for(const Foam::entry& e:*rd)if(!known.found(e.keyword()))FatalIOErrorInFunction(*rd)<<"Unknown re-speciation key "<<e.keyword()<<exit(FatalIOError);
  if ((engine != "kernel" && engine != "GEMS") || exchange.mock())
    FatalIOErrorInFunction(dict) << "Re-speciation requires engine kernel or GEMS and a phase-change model"
      << exit(FatalIOError);
  if (dict.isDict("HKSCoeffs") && dict.subDict("HKSCoeffs").getOrDefault<word>("speciation","none") == "lagged")
    FatalIOErrorInFunction(dict) << "Re-speciation cannot be combined with HKSCoeffs speciation lagged"
      << exit(FatalIOError);
  minT_ = rd->getOrDefault<scalar>("minTemperature",350);
  interval_ = rd->getOrDefault<label>("updateInterval",1);
  if (!std::isfinite(minT_) || minT_<=0 || interval_<=0)
    FatalIOErrorInFunction(*rd) << "Require positive finite minTemperature and updateInterval" << exit(FatalIOError);
  const wordList selected(rd->get<wordList>("species"));
  if (selected.empty()) FatalIOErrorInFunction(*rd) << "Empty re-speciation species" << exit(FatalIOError);
  IOdictionary props(IOobject("speciesTransportProperties",mesh.time().constant(),mesh,
    IOobject::MUST_READ,IOobject::NO_WRITE,false));
  ids_.resize(selected.size()); masses_.resize(selected.size());
  dictionary record; record.add("engine",engine); record.add("ledgerVersion",1); record.add("species",selected);
  record.add("minTemperature",string(phaseChangeLedger::exactDecimal(minT_))); record.add("updateInterval",interval_);
  dictionary formulaRecord, massRecord;
  std::array<bool,3> used = {false,false,false}, mono = {false,false,false};
  forAll(selected,i) {
    const label si = names.find(selected[i]);
    if (si<0 || state[si]!="gas" || index_[si]>=0)
      FatalIOErrorInFunction(*rd) << "Invalid or duplicate re-speciation species " << selected[i] << exit(FatalIOError);
    ids_[i]=si; index_[si]=i; masses_[i]=masses[si];
    massRecord.add(selected[i],string(phaseChangeLedger::exactDecimal(masses_[i])));
    const chemicalFormula f = readFormula(props.subDict(selected[i]),masses_[i]);
    if (f.empty()) FatalIOErrorInFunction(*rd) << "Formula required for " << selected[i] << exit(FatalIOError);
    gasSpecies g{0,0,0}; int* counts[] = {&g.nPb,&g.nBi,&g.nI};
    const char* atoms[] = {"Pb","Bi","I"};
    dictionary fr;
    for (const auto& a : f) {
      int e=0; while (e<3 && a.first!=atoms[e]) ++e;
      if (e==3 || a.second!=std::floor(a.second))
        FatalIOErrorInFunction(*rd) << "Kernel requires integer Pb/Bi/I formulas" << exit(FatalIOError);
      *counts[e]=int(a.second); used[e]=true; fr.add(a.first,a.second);
    }
    if (g.nPb && g.nBi) FatalIOErrorInFunction(*rd) << "Kernel does not support mixed Pb-Bi species" << exit(FatalIOError);
    for (int e=0;e<3;++e) if (*counts[e]==1 && g.nPb+g.nBi+g.nI==1) mono[e]=true;
    formulas_.push_back(g); formulaRecord.add(selected[i],fr);
    mesh.setFluxRequired(fields[si].name());
    if(exchange.pairOfGas(si)<0)checkTransport(si);
  }
  for (int e=0;e<3;++e) if (used[e]&&!mono[e])
    FatalIOErrorInFunction(*rd) << "Include the monatomic gas of every represented element" << exit(FatalIOError);
  // A paired gas carrying these elements must participate, otherwise the
  // observable and its element ledger would silently omit a reactive store.
  forAll(names,si) if (exchange.pairOfGas(si)>=0 && index_[si]<0)
    FatalIOErrorInFunction(*rd) << "Include every paired gas in respeciation species: " << names[si] << exit(FatalIOError);
  for(int e=0;e<3;++e) if(used[e]) exchange.includeElement(word(e==0?"Pb":e==1?"Bi":"I"));
  record.add("formulas",formulaRecord); record.add("molarMasses",massRecord);
  if(engine=="GEMS") {
    dictionary cfg(dict.subOrEmptyDict("GEMSCoeffs"));
    if(rd->isDict("GEMS"))cfg.merge(rd->subDict("GEMS"));
    gems_.reset(new gemsGasSpeciation(mesh,cfg,selected,formulas_,masses_));
    record.add("backendHash",gems_->signature());
  }
  fileName file(rd->getOrDefault<fileName>("formationConstants",fileName::null));
  std::uint64_t hash=14695981039346656037ULL;
  if(!file.empty()) {
    file.replaceAll("<constant>",mesh.time().globalPath()/mesh.time().constant());file.expand();
    std::ifstream raw(file.c_str());if(!raw)FatalIOErrorInFunction(*rd)<<"Cannot read "<<file<<exit(FatalIOError);
    char ch;while(raw.get(ch)){hash^=static_cast<unsigned char>(ch);hash*=1099511628211ULL;}
  } else if(engine=="kernel")FatalIOErrorInFunction(*rd)<<"Kernel requires formationConstants"<<exit(FatalIOError);
  record.add("formationHash",string(std::to_string(hash)));
  dictionary& layout = exchange.ledger().layout();
  if (exchange.ledger().restarted()) {
    const dictionary* previous=layout.findDict("respeciation");
    bool equal = previous && previous->get<word>("engine")==engine
      && previous->get<label>("ledgerVersion")==1
      && previous->get<wordList>("species")==selected
      && previous->get<string>("minTemperature")==string(phaseChangeLedger::exactDecimal(minT_))
      && previous->get<label>("updateInterval")==interval_
      && previous->get<string>("formationHash")==string(std::to_string(hash));
    if(equal && gems_.valid())equal=previous->getOrDefault<string>("backendHash",string::null)==gems_->signature();
    if(equal) {
      equal=previous->isDict("molarMasses");
      if(equal)for(const word& n:selected)equal=equal && previous->subDict("molarMasses").get<string>(n)==massRecord.get<string>(n);
    }
    if(equal) {
      const dictionary& old=previous->subDict("formulas");
      for(const word& n:selected) {
        if(!old.isDict(n) || old.subDict(n).toc()!=formulaRecord.subDict(n).toc()) {equal=false;break;}
        for(const Foam::entry& e:formulaRecord.subDict(n))
          if(old.subDict(n).get<scalar>(e.keyword())!=formulaRecord.subDict(n).get<scalar>(e.keyword()))equal=false;
      }
    }
    if (!equal)
      FatalIOErrorInFunction(*rd) << "Re-speciation configuration changed across restart" << exit(FatalIOError);
  }
  layout.set("respeciation",record);
  std::vector<std::string> columns; for (const auto& n:selected) columns.push_back(n);
  try {
    autoPtr<vapourPressureTable> table;
    if(!file.empty())table.reset(new vapourPressureTable(vapourPressureTable::read(file,columns,vapourPressureTable::LOG10_PA,
      vapourPressureTable::LOG_INVERSE_T,vapourPressureTable::FATAL)));
    constants_.resize(mesh.nCells());
    std::map<std::pair<scalar,scalar>,std::vector<double>> cache;
    const auto began=std::chrono::steady_clock::now();
    label misses=0;
    forAll(thermo.T(),c) if (thermo.T()[c]>=minT_) {
      constants_[c].resize(ids_.size());
      if(gems_.valid()) {
        const auto key=std::make_pair(thermo.T()[c],thermo.p()[c]);
        auto found=cache.find(key);
        if(found==cache.end()) {
          std::vector<double> data;
          if(!gems_->constants(key.first,key.second,data)) {
            ++misses;
            if(!table.valid())throw std::runtime_error("Frozen T/P outside GEMS grid: provide formationConstants for the kernel fallback");
            data.resize(ids_.size());
            forAll(ids_,i)data[i]=std::log(table->phase(i,key.first));
          }
          found=cache.emplace(key,std::move(data)).first;
        }
        constants_[c]=found->second;
      } else forAll(ids_,i)constants_[c][i]=std::log(table->phase(i,thermo.T()[c]));
      forAll(ids_,i) if (formulas_[i].nPb+formulas_[i].nBi+formulas_[i].nI==1 && std::fabs(constants_[c][i])>1e-12)
        throw std::runtime_error("monatomic formation constants must be zero");
    }
    if(gems_.valid())Info<<"Gas GEMS standard states: "<<cache.size()<<" unique T/P states; grid fallbacks "<<misses
      <<"; startup seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<nl;
  } catch (const std::exception& e) { FatalIOErrorInFunction(*rd) << e.what() << exit(FatalIOError); }
  kernel_.reset(new gasSpeciation(formulas_));
  values_.reset(new binaryGlobalIOList<scalar>(IOobject("phaseChangeSpecies",mesh.time().timeName(),
    "uniform",mesh.time(),IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE)));
  elements_.reset(new binaryGlobalIOList<scalar>(IOobject("phaseChangeElements",mesh.time().timeName(),
    "uniform",mesh.time(),IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE)));
  condensates_.reset(new binaryGlobalIOList<scalar>(IOobject("phaseChangeCondensates",mesh.time().timeName(),
    "uniform",mesh.time(),IOobject::READ_IF_PRESENT,IOobject::AUTO_WRITE)));
  label nCondensates=0; for(label si:ids_) if(exchange.pairOfGas(si)>=0)++nCondensates;
  const bool restart=exchange.ledger().restarted();
  if(!restart && (!values_->empty() || !elements_->empty() || !condensates_->empty()))
    FatalIOErrorInFunction(*rd)<<"Re-speciation accounts exist without the pair restart ledger"<<exit(FatalIOError);
  if(restart && condensates_->size()!=8*nCondensates)
    FatalIOErrorInFunction(*rd)<<"Missing or invalid condensate restart ledger"<<exit(FatalIOError);
  if(!restart)condensates_->resize(8*nCondensates,0);
  label ci=0;
  for(label si:ids_) {
    const label p=exchange.pairOfGas(si);if(p<0)continue;
    const scalar wall=exchange.wallMass(p),sample=exchange.sampleMass(p);
    if(!restart){(*condensates_)[8*ci]=wall;(*condensates_)[8*ci+1]=sample;}
    else (*condensates_)[8*ci+7]+=(wall+sample)-((*condensates_)[8*ci+2]+(*condensates_)[8*ci+3]);
    (*condensates_)[8*ci+2]=wall; (*condensates_)[8*ci+3]=sample;
    ++ci;
  }
  if (restart && (values_->size()!=2*NENTRY*ids_.size() || elements_->size()!=6))
    FatalIOErrorInFunction(*rd) << "Missing or invalid re-speciation restart ledger" << exit(FatalIOError);
  if (!restart) { values_->resize(2*NENTRY*ids_.size(),0); elements_->resize(6,0); }
  flux_.resize(ids_.size(),0); defect_.resize(ids_.size(),0);
  forAll(ids_,i) {
    const scalar held=inventory(i);
    if (!restart) { set(i,INITIAL,held); set(i,HELD,held); }
    else { add(i,RESTART,held-value(i,HELD)); set(i,HELD,held); }
    const label pi=exchange.pairOfGas(ids_[i]);
    if (pi>=0) exchange.ledger().setReaction(pi,value(i,REACTION)+value(i,COPRODUCT));
  }
  if(dict.isDict("wallChannels")) {
    walls_.reset(new wallChannels(mesh,thermo,rho,fields,names,state,masses,dict,exchange));
    forAll(names,si)if(walls_->uses(si) && index_[si]<0)
      FatalIOErrorInFunction(*rd)<<"Include every wall channel/source gas in re-speciation species"<<exit(FatalIOError);
  } else if(dict.found("gasSources") || exchange.ledger().layout().found("wallChannelsHash"))
    FatalIOErrorInFunction(dict)<<"Missing wallChannels for sources/restart"<<exit(FatalIOError);
  active_=true; exchange.attachRespeciation(*this);
  if (Pstream::master()) {
    const fileName base(mesh.time().globalPath()/"postProcessing");
    const fileName dir(base/"phaseChangeSpecies"/mesh.time().timeName()); mkDir(dir);
    files_.resize(ids_.size());
    forAll(ids_,i) { files_.set(i,new OFstream(dir/(selected[i]+".dat"))); files_[i].precision(17);
      files_[i]<<"# columns time initial released transport exchanged solverDefect restart reaction coproduct held reactionVariation supplied closure reference"<<nl; }
    const fileName ed(base/"phaseChangeElements"/mesh.time().timeName()); mkDir(ed);
    elementFiles_.resize(3); const char* atoms[]={"Pb","Bi","I"};
    for(int e=0;e<3;++e) { elementFiles_.set(e,new OFstream(ed/(word(atoms[e])+".dat"))); elementFiles_[e].precision(17);
      elementFiles_[e]<<"# amounts [mol of atoms]; CONVERSION is the measured residue of gas transfers"<<nl
        <<"# columns time initial released transport removed clamped solverDefect restart gas wall sample held supplied closure pairClosure reference conversion reaction"<<nl; }
  }
  Info<<"Re-speciation: engine "<<engine<<"; species "<<selected<<"; minTemperature "<<minT_<<" K; updateInterval "<<interval_<<nl;
}
void LESTO::gasRespeciation::recordSolve(label si,const fvScalarMatrix& eqn) {
  const label i=index_[si]; const auto flux=eqn.flux(); compensatedSum boundary,defect;
  for (const label p:exchange_.ledger().patchIds()) for (scalar f:flux().boundaryField()[p]) boundary.add(f);
  const auto& Y=fields_[si]; scalarField mass(Y.size()),oldMass(Y.size());
  forAll(Y,c) { mass[c]=rho_[c]*Y[c]*mesh_.V()[c]; oldMass[c]=rho_[c]*Y.oldTime()[c]*mesh_.V()[c]; }
  const tmp<fvScalarMatrix> ddt=fvm::ddt(rho_,Y);
  List<long double> residual; scalarField remainder;
  cellDefect(eqn,ddt().source(),oldMass,mass,mesh_.time().deltaTValue(),Y,residual,remainder);
  forAll(residual,cell) {
    if(walls_.valid())residual[cell]+=static_cast<long double>(walls_->rate(si)[cell])+static_cast<long double>(walls_->release(si)[cell]);
    defect.add(scalar(residual[cell]));
  }
  flux_[i]=globalAmount(boundary); defect_[i]=globalAmount(defect);
}
void LESTO::gasRespeciation::applyWallProducts() {
  if(!walls_.valid())return;
  scalarList amounts;walls_->products(amounts);
  forAll(ids_,i) {
    add(i,COPRODUCT,amounts[ids_[i]]);
    const label p=exchange_.pairOfGas(ids_[i]);if(p>=0)exchange_.ledger().setReaction(p,value(i,REACTION)+value(i,COPRODUCT));
  }
  for(int e=0;e<3;++e)compensatedAdd((*elements_)[e],(*elements_)[3+e],walls_->elementConversion(e));
}
void LESTO::gasRespeciation::apply() {
  if (!active_ || (mesh_.time().timeIndex()-1)%interval_!=0) return;
  std::vector<compensatedSum> reaction(ids_.size()),variation(ids_.size()); std::array<compensatedSum,3> conversion;
  std::vector<double> c(ids_.size()),before(ids_.size());
  label calls=0,iterations=0,skipped=0,ipmCalls=0,ipmIterations=0,failed=0,guarded=0; scalar largest=0;
  scalar warmSeconds=0;
  forAll(thermo_.T(),cell) {
    if (constants_[cell].empty()) continue;
    std::array<long double,3> exact={0,0,0};
    forAll(ids_,i) {
      before[i]=rho_[cell]*fields_[ids_[i]][cell]*mesh_.V()[cell];
      const long double amount=static_cast<long double>(rho_[cell])*fields_[ids_[i]][cell]/masses_[i];
      const auto& f=formulas_[i]; exact[0]+=f.nPb*amount;exact[1]+=f.nBi*amount;exact[2]+=f.nI*amount;
    }
    bool invalid=false; for (auto b:exact) invalid=invalid || b<0 || !std::isfinite(b);
    if (invalid) { ++skipped;continue; }
    std::array<double,3> b={double(exact[0]),double(exact[1]),double(exact[2])}, pi={NAN,NAN,NAN};
    const auto cellBegan=std::chrono::steady_clock::now();
    bool polished=false,charged=false;
    if(gems_.valid() && gems_->local()) {
      scalar trace=0;forAll(ids_,i)trace+=max(rho_[cell]*fields_[ids_[i]][cell]/masses_[i],scalar(0));
      int count=0;bool guard=false;
      const bool ok=gems_->warm(cell,thermo_.T()[cell],thermo_.p()[cell],b,trace,pi,count,guard);
      polished=ok;charged=!guard;
      if(guard)++guarded;else {++ipmCalls;ipmIterations+=count;}
      if(!ok){if(!guard)++failed;pi={NAN,NAN,NAN};}
    }
    int it=kernel_->solve(constants_[cell].data(),1e5/(8.314462618*thermo_.T()[cell]),b,pi,c.data());
    if(it<0 && gems_.valid() && gems_->local()){++failed;pi={NAN,NAN,NAN};it=kernel_->solve(constants_[cell].data(),1e5/(8.314462618*thermo_.T()[cell]),b,pi,c.data());}
    // Canonical final refinement prevents an unsaved IPM warm state from
    // changing nearly exhausted minor gases in a restart. It uses the same
    // physical bulk and G0 as the warm polish, with no composition floors.
    if(polished && it>=0){pi={NAN,NAN,NAN};const int canonical=kernel_->solve(constants_[cell].data(),1e5/(8.314462618*thermo_.T()[cell]),b,pi,c.data());it=canonical<0?canonical:it+canonical;}
    if(charged)warmSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-cellBegan).count();
    if (it<0) FatalErrorInFunction<<"Gas speciation failed in cell "<<cell
      <<"; bulk Pb/Bi/I "<<b[0]<<" / "<<b[1]<<" / "<<b[2]
      <<", T "<<thermo_.T()[cell]<<exit(FatalError);
    std::array<long double,3> after={0,0,0};
    forAll(ids_,i) { const auto& f=formulas_[i];
      if (!std::isfinite(c[i]) || c[i]<0) FatalErrorInFunction<<"Nonfinite gas equilibrium"<<exit(FatalError);
      after[0]+=f.nPb*static_cast<long double>(c[i]);after[1]+=f.nBi*static_cast<long double>(c[i]);after[2]+=f.nI*static_cast<long double>(c[i]); }
    for(int e=0;e<3;++e) if (b[e]>0) {
      const scalar error=scalar(std::fabs(after[e]/b[e]-1));
      if(error>1e-14)FatalErrorInFunction<<"Gas equilibrium element residue "<<error<<" in cell "<<cell
        <<"; element "<<e<<", before "<<b[e]<<", after "<<double(after[e])
        <<", T "<<thermo_.T()[cell]<<exit(FatalError);
      largest=max(largest,error);
    }
    ++calls; iterations+=it;
    forAll(ids_,i) {
      fields_[ids_[i]].primitiveFieldRef()[cell]=c[i]*masses_[i]/rho_[cell];
      const scalar mass=rho_[cell]*fields_[ids_[i]][cell]*mesh_.V()[cell];
      const scalar delta=mass-before[i]; reaction[i].add(delta); variation[i].add(mag(delta));
      const auto& f=formulas_[i]; const int nu[]={f.nPb,f.nBi,f.nI};
      for(int e=0;e<3;++e) conversion[e].add(nu[e]*delta/masses_[i]);
    }
  }
  forAll(ids_,i) { add(i,REACTION,globalAmount(reaction[i])); add(i,REACTION_VARIATION,globalAmount(variation[i])); fields_[ids_[i]].correctBoundaryConditions();
    const label p=exchange_.pairOfGas(ids_[i]); if(p>=0) exchange_.ledger().setReaction(p,value(i,REACTION)+value(i,COPRODUCT)); }
  for(int e=0;e<3;++e) compensatedAdd((*elements_)[e],(*elements_)[3+e],globalAmount(conversion[e]));
  reduce(calls,sumOp<label>());reduce(iterations,sumOp<label>());reduce(skipped,sumOp<label>());reduce(largest,maxOp<scalar>());
  if(gems_.valid() && gems_->local()){reduce(ipmCalls,sumOp<label>());reduce(ipmIterations,sumOp<label>());reduce(failed,sumOp<label>());reduce(guarded,sumOp<label>());reduce(warmSeconds,sumOp<scalar>());
    Info<<"Gas GEMS local: calls "<<ipmCalls<<"; iterations "<<ipmIterations<<"; failures "<<failed<<"; guards "<<guarded<<"; mean call us "<<(ipmCalls?1e6*warmSeconds/ipmCalls:0)<<nl;}
  scalarList maximum(ids_.size(),0),negative(ids_.size(),0);
  forAll(ids_,i)forAll(thermo_.T(),cell) {
    const scalar Y=fields_[ids_[i]][cell];
    maximum[i]=max(maximum[i],rho_[cell]*Y/masses_[i]*8.314462618*thermo_.T()[cell]/thermo_.p()[cell]);
    if(Y<0)negative[i]+=-rho_[cell]*Y*mesh_.V()[cell];
  }
  Pstream::listCombineGather(maximum,maxEqOp<scalar>());Pstream::broadcast(maximum);
  Pstream::listCombineGather(negative,plusEqOp<scalar>());Pstream::broadcast(negative);
  forAll(ids_,i)Info<<fields_[ids_[i]].name()<<" re-speciation monitors: negative mass "<<negative[i]
    <<" kg; max mole fraction "<<maximum[i]<<nl;
  Info<<"Re-speciation: cells "<<calls<<"; calls "<<calls<<"; iterations "<<iterations
      <<"; failures 0; invalid-element cells "<<skipped<<"; max element residue "<<largest<<nl;
}
void LESTO::gasRespeciation::book() {
  if(!active_)return;
  const scalar dt=mesh_.time().deltaTValue();
  if(walls_.valid())walls_->finish();
  label ci=0;
  forAll(ids_,i) {
    const label p=exchange_.pairOfGas(ids_[i]);
    if(p>=0) {
      compensatedSum flux;for(scalar f:exchange_.stepTransport(p))flux.add(f);
      add(i,RELEASED,exchange_.stepReleased(p));add(i,TRANSPORT,flux.value()*dt);
      add(i,DEFECT,exchange_.stepDefect(p)*dt);add(i,EXCHANGED,exchange_.stepExchange(p)*dt);
      (*condensates_)[8*ci+2]=exchange_.cachedWallMass(p);(*condensates_)[8*ci+3]=exchange_.cachedSampleMass(p);
      (*condensates_)[8*ci+4]=-value(i,EXCHANGED);
      (*condensates_)[8*ci+5]=exchange_.ledger()(p,phaseChangeLedger::REMOVED);
      (*condensates_)[8*ci+6]=exchange_.ledger()(p,phaseChangeLedger::CLAMPED);++ci;
    } else {
      add(i,TRANSPORT,flux_[i]*dt); add(i,DEFECT,defect_[i]*dt);
      if(walls_.valid()){add(i,EXCHANGED,walls_->totalRate(ids_[i])*dt);add(i,RELEASED,walls_->totalRelease(ids_[i]));}
    }
    const scalar held=inventory(i);set(i,HELD,held);
    compensatedSum supplied;for(entry e:{INITIAL,RELEASED,EXCHANGED,RESTART,REACTION,COPRODUCT})supplied.add(value(i,e));
    supplied.add(-value(i,TRANSPORT));supplied.add(-value(i,DEFECT));
    const scalar ref=max(mag(value(i,INITIAL))+mag(value(i,RELEASED))+value(i,REACTION_VARIATION)+mag(value(i,COPRODUCT))+(walls_.valid()?walls_->exchangeReference(ids_[i]):0),VSMALL);
    if(Pstream::master()) { auto& os=files_[i];os<<exchange_.ledger().time();
      for(int e=0;e<NENTRY;++e)os<<' '<<value(i,entry(e));
      os<<' '<<supplied.value()<<' '<<held-supplied.value()<<' '<<ref<<nl;os.flush(); }
  }
}
void LESTO::gasRespeciation::writeElements() {
  if(!active_ || !Pstream::master())return;
  auto& l=exchange_.ledger();
  for(int e=0;e<3;++e) {
    std::array<compensatedSum,16> t;
    forAll(ids_,i) {
      const auto& f=formulas_[i];const int nu[]={f.nPb,f.nBi,f.nI};const scalar k=nu[e]/masses_[i];if(!k)continue;
      const label p=exchange_.pairOfGas(ids_[i]);
      // Wall/sample inventories have already been reduced by balance().
      if(p>=0) {
        for(auto pair:std::vector<std::pair<int,phaseChangeLedger::entry>>{{0,phaseChangeLedger::INITIAL},{1,phaseChangeLedger::RELEASED},{3,phaseChangeLedger::REMOVED},{4,phaseChangeLedger::CLAMPED},{5,phaseChangeLedger::SOLVER_DEFECT},{6,phaseChangeLedger::RESTART}})
          t[pair.first].add(k*l(p,pair.second));
        t[2].add(k*l.totalTransport(p));
        // Total held minus gas gives the condensed stores without collectives.
        t[8].add(k*exchange_.cachedWallMass(p));t[9].add(k*exchange_.cachedSampleMass(p));
        t[13].add(k*(l(p,phaseChangeLedger::HELD)-l.supplied(p)));t[14].add(k*l.reference(p));
      } else {
        t[0].add(k*value(i,INITIAL));t[1].add(k*value(i,RELEASED));t[2].add(k*value(i,TRANSPORT));
        t[5].add(k*value(i,DEFECT));t[6].add(k*value(i,RESTART));
        t[14].add(k*(mag(value(i,INITIAL))+mag(value(i,RELEASED))));
      }
      t[7].add(k*value(i,HELD));t[15].add(k*value(i,REACTION)+value(i,COPRODUCT));
    }
    if(walls_.valid()) {
      t[0].add(walls_->elementStore(e,0));t[8].add(walls_->elementStore(e,1));t[6].add(walls_->elementStore(e,2));t[4].add(walls_->elementStore(e,3));
      t[14].add(mag(walls_->elementStore(e,0)));
    }
    compensatedSum held,supplied;for(int j:{7,8,9})held.add(t[j].value());
    supplied.add(t[0].value());supplied.add(t[1].value());for(int j:{2,3,4,5})supplied.add(-t[j].value());supplied.add(t[6].value());
    const scalar conv=(*elements_)[e]+(*elements_)[3+e];supplied.add(conv);
    auto& os=elementFiles_[e];os<<l.time();for(int j=0;j<10;++j)os<<' '<<t[j].value();
    os<<' '<<held.value()<<' '<<supplied.value()<<' '<<held.value()-supplied.value()<<' '<<t[13].value()<<' '<<max(t[14].value(),mag(held.value()))<<' '<<conv<<' '<<t[15].value()<<nl;os.flush();
  }
}
void LESTO::gasRespeciation::addGasProfiles(const word& element,scalarField& gas,scalar& inventory)const {
  const int e=element=="Pb"?0:element=="Bi"?1:2;
  forAll(ids_,i) if(exchange_.pairOfGas(ids_[i])<0) {
    const auto& f=formulas_[i];const int nu[]={f.nPb,f.nBi,f.nI};const scalar k=nu[e]/masses_[i];
    forAll(gas,c)gas[c]+=k*rho_[c]*fields_[ids_[i]][c]*mesh_.V()[c];
    inventory+=k*value(i,HELD);
  }
}
void LESTO::gasRespeciation::checkTransport(label si) const {
  const word& name=fields_[si].name();
  ITstream& ddt=mesh_.ddtScheme("ddt("+rho_.name()+','+name+')');ddt.rewind();const word scheme(ddt);ddt.rewind();
  ITstream& div=mesh_.divScheme("div(phi,"+name+')');div.rewind();const word convection(div);const string text(div.toString());div.rewind();
  if(scheme!="Euler" || convection!="Gauss" || text.find("linearUpwind")!=string::npos || text.find("LUST")!=string::npos)
    FatalErrorInFunction<<"Ledgered gas "<<name<<" needs Euler and conservative Gauss convection without explicit source corrections"<<exit(FatalError);
  const word finalName(name+"Final");
  const bool haveFinal=mesh_.solution().solversDict().found(finalName);
  if(mesh_.solution().solverDict(haveFinal?finalName:name).getOrDefault<scalar>("relTol",0)!=0)
    FatalErrorInFunction<<"Ledgered gas "<<name<<" needs relTol 0 for its final solve"<<exit(FatalError);
}
