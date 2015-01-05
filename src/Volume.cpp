#include "Volume.hpp"

ostream& operator<<(ostream& output, const Volume& vol){
  output<<"Center at "<<setw(10)<<left<<vol.center<<"Density = "<<setw(10)<<left<<vol.density;
  return output;
  
}

Volume::Volume():center(0,0,0),density(1){
  
}

Volume::Volume(const Point& center, const double density):center(center),density(density){

}

void Volume::SetCenter(Point center){

  this->center = center;
  
}


void Volume::SetDensity(const double density){
  
  this->density = density;

}

