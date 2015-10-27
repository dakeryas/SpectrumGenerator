#include "FourMomentum.hpp"

std::ostream& operator<<(std::ostream& output, const FourMomentum& P){

  output<<"E = "<<std::setw(12)<<std::left<<P.GetE()<<"px = "<<std::setw(12)<<std::left<<P.Getpx()<<"py = "<<std::setw(12)<<std::left<<P.Getpy()<<"pz = "<<std::setw(12)<<std::left<<P.Getpz();
  return output;
  
}

FourMomentum::FourMomentum():E(0),px(0),py(0),pz(0){
  
}

FourMomentum::FourMomentum(double m):E(m),px(0),py(0),pz(0){
  
}

FourMomentum::FourMomentum(double E, double px, double py, double pz):E(E),px(px),py(py),pz(pz){
  
}

FourMomentum::FourMomentum(double E, double p, TRandom3& ran):E(E){
  
  SetRandomDir(p, ran);
  
}

FourMomentum FourMomentum::operator-(){

  return FourMomentum(-E, -px, -py, -pz);
  
}

FourMomentum FourMomentum::GetReverse(){
  
  return FourMomentum (E, -px, -py, -pz);
  
}

void FourMomentum::SetDir(const FourMomentum& Q){ 
  
  px = Q.Getpx();
  py = Q.Getpy();
  pz = Q.Getpz();

}

void FourMomentum::SetE(double Earg){ 
  
  E = Earg;

}

void FourMomentum::SetRandomDir(double p, TRandom3& ran){
  
  double cs = ran.Uniform(-1,1);
  double phi = ran.Uniform(0,2*3.14159);
  px = p*sqrt(1-pow(cs,2))*cos(phi);
  py = p*sqrt(1-pow(cs,2))*sin(phi);
  pz = p*cs;
  
}

double FourMomentum::GetSquare() const{
  
  return pow(E,2)-pow(px,2)-pow(py,2)-pow(pz,2);
  
}

double FourMomentum::GetInvariant() const{
  
  if(abs(E-Getp())<1e-6) return 0;//if for the machine, we have zero, do return 0
  else return sqrt(GetSquare());
  
} 

double FourMomentum::Getp() const{
  
  return sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
  
}

double FourMomentum::GetKinetic() const{
  
  return E-GetInvariant();
  
}

const double& FourMomentum::GetE() const{
  
  return E;
  
}

const double& FourMomentum::Getpx() const{
  
  return px;
  
}

const double& FourMomentum::Getpy() const{
  
  return py;
  
}

const double& FourMomentum::Getpz() const{
  
  return pz;
  
}


void FourMomentum::SetRest(double m){ 
  
  E = m;
  px = 0;
  py = 0;
  pz = 0;
  
}

void FourMomentum::Boost(const FourMomentum& P){
  
  double gamma = P.E/P.GetInvariant();
  double b = sqrt(1-1/pow(gamma,2));
  double bx = b*P.Getpx()/P.Getp();
  double by = b*P.Getpy()/P.Getp();
  double bz = b*P.Getpz()/P.Getp();
  
  double E_boosted = gamma*(E-bx*px-by*py-bz*pz);
  double px_boosted = -gamma*bx*E+(1+(gamma-1)*bx*bx/pow(b,2))*px+((gamma-1)*by*bx/pow(b,2))*py+((gamma-1)*bz*bx/pow(b,2))*pz;
  double py_boosted = -gamma*by*E+((gamma-1)*bx*by/pow(b,2))*px+(1+(gamma-1)*by*by/pow(b,2))*py+((gamma-1)*bz*by/pow(b,2))*pz;
  double pz_boosted = -gamma*bz*E+((gamma-1)*bx*bz/pow(b,2))*px+((gamma-1)*by*bz/pow(b,2))*py+(1+(gamma-1)*bz*bz/pow(b,2))*pz;
  
  E = E_boosted;
  px = px_boosted;
  py = py_boosted;
  pz = pz_boosted;
  
}