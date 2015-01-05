#ifndef FOURMOMENTUM_H
#define FOURMOMENTUM_H

#include <cmath>
#include <iostream>
#include <iomanip>  
#include <TRandom3.h>

using namespace::std;

class FourMomentum {
  
friend ostream& operator<<(ostream& output, const FourMomentum& P);
  
  double E,px,py,pz;

public:
  FourMomentum(); 
  FourMomentum(const double m); //generate the fourmomentum to its rest mass m
  FourMomentum(const double E, const double px, const double py, const double pz); //create a fourmomentum with the given values
  FourMomentum(const double E, const double p, TRandom3& ran);//initialised with a random direction
  FourMomentum operator-();//returns a FourMomentum with all opposite components
  void SetRest(const double m);
  void SetDir(const FourMomentum& Q); //set(px, py, pz) according to(qx, qy, qz)
  void SetE(const double Earg);
  void SetRandomDir(const double p, TRandom3& ran); //set a random direction for(px, py, pz) of norm p
  FourMomentum GetReverse();//returns a FourMomentum with oposite spatial parts
  const double& GetE() const;
  const double& Getpx() const;
  const double& Getpy() const;
  const double& Getpz() const;
  double Getp() const;
  double GetKinetic() const; 
  double GetSquare() const;
  double GetInvariant() const;//returns sqrt(square())
  void Boost(const FourMomentum& P);//apply a  boost according to P to a fourmomentum
  
};

#endif