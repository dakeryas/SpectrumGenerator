#ifndef STATE
#define STATE

#include <vector>
#include <algorithm>
#include "FourMomentum.hpp"

class State { //variables used in the decay from a particle/state of mass m1 and width w1 to particles of masses m2 and m3 and corresponding widths

friend ostream& operator<<(ostream& output, const State* s);
  
protected:
  FourMomentum P;
  double m0; //average mass of the state /!\ doesn't not correspond to the current mass (P.GetInvariant()) for RStates
  vector<State*> children;
  vector<double> rsum; //the triangle sum of the branching ratios associated to the decay paths (the size shouldn't exceed that of decays)
  void IncorporateRatios(vector<double> r); //turns the branching ratios into a useful form 'rsum', a shrunk copy of r according to decays.size() is made first
  
public:
  State();
  State(const FourMomentum& P);
  State(const FourMomentum& P, const vector<State*>& children);//particularly meant for VStates
  State(const FourMomentum& P, const vector<State*>& children, const vector<double>& r);
  virtual ~State(){}
  State* Next(TRandom3& ran); //select a decay path for the state
  const FourMomentum& GetP() const;
  double GetKinetic() const;
  double GetMass() const;//return P.GetInvariant()
  const State* GetChild(unsigned k) const;
  const vector<State*>& GetChildren() const;
  void SetP(const FourMomentum& Pnew);
  virtual void Process(TRandom3&)=0;//decay for virtual states, propagate P to the children for RStates
  //for VStates
  virtual vector<State*> GetStableChildren(){return vector<State*>(0);}
  virtual bool CanDecay(){return true;}
  //for RStates
  virtual bool IsStable(){return false;}
  virtual int GetPDG() const {return 0;}
  virtual void AdaptMass(TRandom3&){}//sets P to a random rest value according to w for RStates /!\ the resulting may be negative
 
};

#endif