#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace constants{

  namespace dchooz{
    
    const double R_target = 1150; //target cylinder's radius in mm
    const double h_target = 2458; //target cylinder's total height in mm
    const double R_GC = 1708;//gamma catcher radius in mm
    const double h_GC = 3572; //gamma catcher total height in mm
    
    const double target_mass = 8287; //NT mass in kg
    const double target_C_wt_density = 0.86;
    const double target_C_mass = target_mass * target_C_wt_density;
    const double GC_actual_capacity = 22542;//actual volume
    const double GC_density = 0.804;
    const double GC_mass = GC_actual_capacity*GC_density;
    const double GC_C_wt_density = 0.853;
    const double GC_C_mass = GC_mass * GC_C_wt_density;
    const double probCosmogenicInNT = target_C_mass/(target_C_mass+GC_C_mass);
    
  }

  namespace PDG{
    
    const int _e = 11; 
    const int _alpha = 9802004; //He4 IDHEP code
    const int _n = 2112; //neutron code
    const int _gamma = 22;//for photons
    const int _3H = 9801003;
    const int _7Li = 9803007;
    const int _5He = 9802005;
    const int _8Be = 9804008;
    const int _8Li = 9803008;
    const int _9Be = 9804009;
    const int _9Li = 9803009;
    const int _8He = 9802008;
    
  }

  namespace mass{
    
    const double _e = 0.511; //electron mass in MeV
    const double _atomic = 931.4941;
    const double _n = 939.565;
    const double _gamma = 0;
    const double _alpha = 3728.401; //atomic mass of He4
    
    const double _9Li = 24.9549+9*_atomic;
    const double _9Be = 11.3484+9*_atomic; //atomic masses, be careful
    const double _5He = 11.2312+5*_atomic;
    const double _8Be = 4.9416+8*_atomic;
  //    const double _9Be = 9.0121822*_atomic;//atomic masses, be careful
  //    const double _5He = 5.012220*_atomic;
  //    const double _8Be = 8.00530510*_atomic;
    
    const double _8He = 31.6096+8*_atomic;
    const double _8Li = 20.9457+8*_atomic;
    const double _3H = 14.9498+3*_atomic;
    const double _7Li = 14.9070+7*_atomic;

    const double _9Be243 = _9Be+2.4294;
    const double _9Be278 = _9Be+2.78;
    const double _9Be794 = _9Be+7.94;
    const double _9Be1128 = _9Be+11.28;
    const double _9Be1181 = _9Be+11.81;
    const double _5He374 = _5He+3.74-2.467;
    const double _8Be303 = _8Be+3.03-1.6654;
    const double _8Be1135 = _8Be+11.35-1.6654;
    
    const double _8Li321 = _8Li+3.21;
    const double _8Li54 = _8Li+5.4;
    const double _8Li967 = _8Li+9.67;
    const double _7Li048 = _7Li+0.477612;
    
  }
  
  namespace width{
    
    const double _9Be243 = 0.78*1e-3;// 0.78 kev only
    const double _9Be278 = 1.10;
    const double _9Be794 = 1;
    const double _9Be1128 = 0.575;
    const double _9Be1181 = 0.4;
    const double _5He = 0.648;
    const double _5He374 = _5He;
    const double _8Be = 5.57*1e-6; //5.57 eV
    const double _8Be303 = 1.513;
    const double _8Be1135 = 3.5;
    
    const double _8Li321 = 1;
    const double _8Li54 = 0.65;
    const double _8Li967 = 1;
    
  }
  
}

#endif