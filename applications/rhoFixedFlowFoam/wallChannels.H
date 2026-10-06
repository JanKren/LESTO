/* M10d: bounded gas channels sharing physical wall condensates. */
#ifndef LESTO_wallChannels_H
#define LESTO_wallChannels_H
#include "interfaceExchange.H"
#include "chemicalFormula.H"
#include "binaryIOList.H"
#include "compensatedSum.H"
#include <vector>
#include <array>
namespace LESTO {
class wallChannels {
  struct store {
    Foam::label species;
    Foam::scalar mass;
    chemicalFormula formula;
    Foam::scalarField deposit;
    Foam::autoPtr<Foam::volScalarField> mw,mDep;
  };
  struct channel {
    Foam::word name;
    Foam::label reactant,condensate;
    Foam::scalar nuCondensate=1,accommodation=1,Ce=1,scale=1;
    bool irreversible=false,reactive=false;
    std::vector<std::pair<Foam::label,Foam::scalar>> products;
    Foam::scalar reactantNu=1;
    Foam::autoPtr<vapourPressureTable> table;
    Foam::scalarField a,g,upper,flow,pEq;
    Foam::labelList regime;
  };
  struct source {
    Foam::scalar start,duration,length,volume;
    Foam::labelList cells;
    std::vector<std::pair<Foam::label,Foam::scalar>> amounts;
  };
  const Foam::fvMesh& mesh_;
  const Foam::fluidThermo& thermo_;
  const Foam::volScalarField& rho_;
  Foam::PtrList<Foam::volScalarField>& fields_;
  interfaceExchange& exchange_;
  Foam::scalarList masses_;
  std::vector<chemicalFormula> formulas_;
  Foam::PtrList<store> stores_;
  Foam::PtrList<channel> channels_;
  std::vector<source> sources_;
  std::vector<std::vector<Foam::label>> elementsOfCell_;
  Foam::List<Foam::scalarField> flows_,release_,pending_;
  Foam::PtrList<Foam::volScalarField::Internal> su_,sp_;
  Foam::autoPtr<binaryIOList<Foam::scalar>> deposits_;
  Foam::autoPtr<binaryGlobalIOList<Foam::scalar>> accounts_;
  Foam::PtrList<Foam::OFstream> files_;
  Foam::scalar fingerprint_=0;
  std::array<compensatedSum,3> conversion_;
  void assemble(Foam::label);
  Foam::scalar storeMass(Foam::label) const;
  void conversion(Foam::label,Foam::scalar);
public:
  wallChannels(const Foam::fvMesh&,const Foam::fluidThermo&,const Foam::volScalarField&,
    Foam::PtrList<Foam::volScalarField>&,const Foam::wordList&,const Foam::List<Foam::word>&,
    const Foam::List<Foam::scalar>&,const Foam::dictionary&,interfaceExchange&);
  bool owns(Foam::label si) const;
  bool uses(Foam::label si) const;
  bool hasSource(Foam::label si) const;
  void beginStep();
  void prepare(Foam::label);
  void linearise(Foam::label,const Foam::fvScalarMatrix&);
  bool guard(Foam::label);
  void realise(Foam::label);
  void products(Foam::scalarList&);
  void finish();
  const Foam::volScalarField::Internal& Su(Foam::label si)const{return su_[si];}
  const Foam::volScalarField::Internal& Sp(Foam::label si)const{return sp_[si];}
  const Foam::scalarField& rate(Foam::label si)const{return flows_[si];}
  const Foam::scalarField& release(Foam::label si)const{return release_[si];}
  Foam::scalar totalRate(Foam::label)const;
  Foam::scalar totalRelease(Foam::label)const;
  Foam::scalar exchangeReference(Foam::label)const;
  Foam::scalar elementStore(Foam::label e,int entry)const;
  Foam::scalar elementConversion(Foam::label e)const;
  void writeProfiles(const Foam::word&,Foam::scalar,const Foam::string&);
  void addProfiles(const Foam::word&,Foam::scalarField&,Foam::scalar&)const;
};
}
#endif
