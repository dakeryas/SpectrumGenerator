#include "InvMassGenerator.hpp"
#include "VState.hpp"

double TwobodyPhaseSpc(double m0, double m1, double m2){

  return sqrt(VState::TrFct(pow(m0, 2), pow(m1, 2), pow(m2, 2)))/(8*3.14159*m0);
  
}

InvMassGenerator::InvMassGenerator():n(0),m(),ran(TRandom3 (0)), Mmin(), M(), y(0), dmax(1){
  
}

InvMassGenerator::InvMassGenerator(const VState* s, TRandom3& ran):n(s->GetNbody()),m(),ran(ran), Mmin(), M(), y(0), dmax(1){

  if(n>2){
    
    InitialiseMasses(s);
    ComputeMaxDistribution();
    
  }
  
}

void InvMassGenerator::InitialiseMasses(const VState* s){
  
  m.resize(n+1);
  m[0] = s->GetMass();
  for(unsigned k = 0; k<n; ++k) m[k+1] = s->GetChild(k)->GetMass();
  
  Mmin.resize(n);
  Mmin[0] = m[0];//the first mass shouldn't change, it's the sqrt(P_squared) that decays
  Mmin[n-1] = m[n];
  for(unsigned k = n-2; k>0; --k) Mmin[k] = Mmin[k+1] + m[k+1];
  
  M.resize(n);
  M[n-1] = m[n];//easier to compute f(M) (distribution) this way
  
}

void InvMassGenerator::ComputeMaxDistribution(){

  double massSum = 0;
  
  for(unsigned k = 1; k<n; ++k){
    
    dmax *= TwobodyPhaseSpc(m[0] - massSum, Mmin[k], m[k]);
    massSum += m[k];
    
  }

}

double InvMassGenerator::GetDistribution(){

  double space = 1/dmax; //so the probablity density is normalised
  for(unsigned k = 0; k<n-1; ++k) space *= TwobodyPhaseSpc(M[k], M[k+1], m[k+1]);//the father and n-2 virutal decays account for n-1 decays and hence for n-1 phase space products
  return space;//we don't need the squares of the masses in the TwobodyPhaseSpc because the recurrence formula gives f(mqsr) dmqsr which translates into f(m) dm
  
}

void InvMassGenerator::Generate(){
  
  double Mmax = m[0];
  
  for(unsigned k = 0; k<n-1; ++k){

    M[k] = ran.Uniform(Mmin[k], Mmax);
    Mmax = M[k] - m[k+1];
    
  }
  
  y = ran.Uniform();

  while(y > GetDistribution()){
    
    Mmax = m[0];
  
    for(unsigned k = 0; k<n-1; ++k){//rejection method if y > f(x), then pick another set(x, y)

      M[k] = ran.Uniform(Mmin[k], Mmax);
      Mmax = M[k] - m[k+1];
      
    }
    
    y = ran.Uniform();
  
  }
  
}

vector<double> InvMassGenerator::GetMassSet(){
  
  Generate();
  return vector<double>(M.begin()+1, M.end()-1);//remove the first element and the last element (M.end()-1 is the last one, that we don't want); M size is at least 3
  
}