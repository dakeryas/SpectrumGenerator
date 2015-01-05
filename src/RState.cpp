#include "RState.hpp"

RState::RState(const FourMomentum& P, const int PDG):State(P),PDG(PDG){
  
}

RState::RState(const FourMomentum& P, const int PDG, const vector<State*>& children):State(P, children), PDG(PDG){
  
  IncorporateRatios(vector<double>(children.size(), 1/children.size()));
  
}

RState::RState(const FourMomentum& P, const int PDG, const vector<State*>& children, const vector<double>& r):State(P, children, r),PDG(PDG){
  
}

RState::RState(const double m0, const int PDG):RState(FourMomentum(m0), PDG){

}

RState::RState(const double m0, const int PDG, const vector<State*>& children):RState(FourMomentum(m0), PDG, children){
 
  IncorporateRatios(vector<double>(children.size(), 1/children.size()));
  
}

RState::RState(const double m0, const int PDG, const vector<State*>& children, const vector<double>& r):RState(FourMomentum(m0), PDG, children, r){
  
}

void RState::SetChildrenP(){
  
  for(State* child : children) child->SetP(P);
  
}

bool RState::IsStable(){
  
  return rsum.empty();
  
}

int RState::GetPDG() const {
  
  return PDG;
  
}

void RState::Process(TRandom3&){
  
  SetChildrenP();
  
}