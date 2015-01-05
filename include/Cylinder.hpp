#ifndef CYLINDER_H
#define CYLINDER_H

#include "Volume.hpp"

class Cylinder: public Volume {
  
public:
  
friend ostream& operator<<(ostream& output, const Cylinder& volume);

  double r,z;
  Cylinder();
  Cylinder(const Point& center, const double r, const double z, const double density = 1);
  ~Cylinder(){}
  double GetCapacity();

};

#endif