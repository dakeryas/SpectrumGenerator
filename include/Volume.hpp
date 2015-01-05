#ifndef VOLUME_H
#define VOLUME_H

#include <iostream>
#include <iomanip>  
#include "Point.hpp"

using namespace std;

class Volume {

protected:
  const double pi = 3.14159;
  Point center;
  double density;
  
public:
  
friend ostream& operator<<(ostream& output, const Volume& volume);

  Volume();
  Volume(const Point& center, const double density);
  virtual ~Volume(){}
  void SetCenter(Point center);
  void SetDensity(const double density);
  virtual double GetCapacity()=0;
  double GetMassOfMaterial(){return density*GetCapacity();}

};

#endif