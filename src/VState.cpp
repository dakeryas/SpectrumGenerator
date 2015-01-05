#include "VState.hpp"

double VState::TrFct(const double x, const double y, const double z){
  
  double t = pow(x,2)+pow(y,2)+pow(z,2)-2*x*y-2*x*z-2*y*z;
  if(t<1e-9) return 0; //we don't want negative values due to computation errors either 
  else return t;
  
}

VState::VState(const FourMomentum& P, const vector<State*>& children):State(P, children){
  
  InitialiseStdRatios();
  
}

VState::VState(const FourMomentum& P, const vector<State*>& children, const vector<double>& r):State(P, children, r){

}


VState::VState(const double m0, const vector<State*>& children):State(FourMomentum(m0), children){
  
  InitialiseStdRatios();
  
}


VState::VState(const double m0, const vector<State*>& children, const vector<double>& r):State(FourMomentum(m0), children, r){
  
}

void VState::InitialiseStdRatios(){
  
  vector<double> r(children.size(), 1);//you can go anywhere
  for(unsigned k = 0; k<r.size(); ++k) if(children[k]->IsStable()) r[k] = 0;//unless some children are already stable
  IncorporateRatios(r);
  
}

State* VState::GetChildOfPDG(const int pdg){
  
  return *find_if(children.begin(),children.end(),[&](State* child){return child->GetPDG() == pdg;});
  
}

State* VState::GetHeaviestChild(){
  
  return *max_element(children.begin(),children.end(),[](State* child1, State* child2){return child1->GetMass() < child2->GetMass();});
  
}

unsigned VState::GetNbody() const{
  
  return children.size();
  
}

vector<State*> VState::GetStableChildren(){
  
  vector<State*> stable;
  for(State* child : children) if(child->IsStable()) stable.push_back(child);
  return stable;
  
}

bool VState::CanDecay(){

  double totalChildrenMass = 0;
  for(State* child : children) totalChildrenMass+=child->GetMass();
  if(GetMass() > totalChildrenMass) return true;
  else return false;
  
}

void VState::TwobodyDecay(TRandom3& ran){
  
  double m0 = GetMass();
  double m1 = children.front()->GetMass();
  double m2 = children.back()->GetMass();
  
  double p = sqrt(TrFct(pow(m0,2),pow(m1,2),pow(m2,2)))/2/m0;
  FourMomentum P1(sqrt(pow(m1,2)+pow(p,2)), p, ran);
  FourMomentum P2(sqrt(pow(m2,2)+pow(p,2)), -P1.Getpx(),-P1.Getpy(),-P1.Getpz());

  if(P.GetKinetic()/m0>1e-9){

    P1.Boost(P.GetReverse()); //take the opposite of the spatial components to go back to the lab frame
    P2.Boost(P.GetReverse());
   
  }
  
  children.front()->SetP(P1);
  children.back()->SetP(P2);

}

void VState::Process(TRandom3& ran){

  unsigned n = GetNbody();
  if(n==2) TwobodyDecay(ran);
  else if(n>2){

    InvMassGenerator MGen(this, ran);
    vector<double> M = MGen.GetMassSet();

    vector<State*> vs(n-1);//you have to decay P first, so there is one more virtual state than invariant masses
    vs[n-2] = new VState(M[n-3], {children[n-2], children[n-1]});//children[children.size()-2] acesses the element before the last
    for(unsigned k = 0; k<n-3; ++k) vs[n-3-k] = new VState(M[n-4-k], {vs[n-2-k], children[n-3-k]});  
    vs[0] = new VState(P, {vs[1], children[0]});

    for(State* st : vs) st->Process(ran);
    
    for(State* st : vs) delete st;
    
  }
    
}
