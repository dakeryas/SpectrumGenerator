#include "State.hpp"

ostream& operator<<(ostream& output, const State* s){//for input masses in MeV sets the file into GeV

  output<<"PDG = "<<setw(10)<<left<<s->GetPDG()<<"Mass = "<<setw(10)<<left<<s->GetMass()<<"m0 = "<<setw(10)<<left<<s->GetAverageMass()<<endl<<s->GetP();
  return output;
  
}

vector<double> TrSum(const vector<double> br, const unsigned end){//does a triangle sum up to the desired end(usually the size of br)
  
  vector<double> br_sum;
  
  double s = 0;
  br_sum.push_back(s);
  
  for(unsigned k = 0; k<end; ++k){
    
    s = s+br[k];
    br_sum.push_back(s);

  }
  
  if(br_sum.back() != 0) for(double& b : br_sum) b = b/br_sum.back(); //normalise to one
  
  return br_sum;
  
}

State::State():P(0), m0(P.GetInvariant()),children(),rsum(){
  
}

State::State(const FourMomentum& P):P(P), m0(P.GetInvariant()),children(),rsum(){
  
}

State::State(const FourMomentum& P, const vector<State*>& children):P(P), m0(P.GetInvariant()),children(children),rsum(vector<double>(children.size(), 0)){
  
}

State::State(const FourMomentum& P, const vector<State*>& children, const vector<double>& r):P(P), m0(P.GetInvariant()),children(children),rsum(){

  IncorporateRatios(r);
  
}

void State::IncorporateRatios(vector<double> r){
  
  unsigned end = children.size();
  
  if(r.size()<end) for(unsigned k = 0; k<end-r.size(); ++k) r.push_back(1./end); //if for some reason there were more children than ratios given, add the most simple ones
  rsum = TrSum(r, end);
  
}

State* State::Next(TRandom3& ran){
    
    if(children.size() != 0){
  
      double r = ran.Uniform();
      unsigned selected_branch = 0;
      
      for(unsigned k = 0; k<rsum.size()-1; ++k){
	
	if(rsum[k] <r && r <rsum[k+1]) selected_branch = k;
	
      }

      return children[selected_branch];
      
    }
    else return this;
    
}

const FourMomentum& State::GetP() const{
  
  return P;
  
}

double State::GetKinetic() const{
  
  return P.GetKinetic();
  
}

double State::GetAverageMass() const{
  
  return m0;

}

double State::GetMass() const{
  
  return P.GetInvariant();
  
}

const State* State::GetChild(unsigned k) const{
  
  if(k<children.size()) return children[k];
  else return this;
  
}

const vector<State*>& State::GetChildren() const{
  
  return children;
  
}

void State::SetP(const FourMomentum& Pnew){
  
  P = Pnew;
  
}