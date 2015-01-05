#include "Cylinder.hpp"

ostream& operator<<(ostream& output, const Cylinder& vol){
  output<<vol<<"r = "<<setw(10)<<left<<vol.r<<"z = "<<setw(10)<<left<<vol.z;
  return output;
  
}

Cylinder::Cylinder():Volume(),r(0),z(0){
  
}

Cylinder::Cylinder(const Point& center, const double r, const double z, const double density):Volume(center, density),r(r),z(z){

}

double Cylinder::GetCapacity(){
  
  return pi*pow(r,2)*z;
  
}




