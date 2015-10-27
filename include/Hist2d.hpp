#ifndef HIST2D_H
#define HIST2D_H

#include <cmath>
#include <TH2F.h>
#include "Hist.hpp"

class Hist2d { 
    
  std::string name;
  std::vector<double> xedge;
  std::vector<double> yedge;
  std::vector<double> value;
  void buildConstantBiningX(unsigned numberOfBinsX, double lowerEdgeX, double upperEdgeX);
  void buildConstantBiningY(unsigned numberOfBinsY, double lowerEdgeY, double upperEdgeY);
  const std::vector<double>& getXedge() const;
  const std::vector<double>& getYedge() const;
  
public:
  Hist2d();
  Hist2d(const std::string& name, const std::vector<double>& xedge, const std::vector<double>& yedge);
  Hist2d(const std::string& name, const std::vector<double>& xedge, const std::vector<double>& yedge, const std::vector<double>& value);
  Hist2d(const std::string& name, unsigned numberOfBinsX, double lowerEdgeX, double upperEdgeX, unsigned numberOfBinsY, double lowerEdgeY, double upperEdgeY);
  Hist2d(const TH2F* hROOT);//copy the contents of the ROOT histogram
  const std::string& getName() const;
  const double& getBinXLowEdge(unsigned i) const;
  double getBinXCenter(unsigned i) const;
  const double& getBinXUpEdge(unsigned i) const;
  double getBinXWidth(unsigned i) const; 
  const double& getBinYLowEdge(unsigned i) const;
  double getBinYCenter(unsigned i) const;
  const double& getBinYUpEdge(unsigned i) const;
  double getBinYWidth(unsigned i) const; 
  unsigned getNumberOfBinsX() const;
  unsigned getNumberOfBinsY() const;
  unsigned getNumberOfBins() const;
  const double& getBinContent(unsigned i, unsigned j) const;
  unsigned getIndexOfBinXContaining(double x) const;
  double getLowEdgeOfBinXContaining(double x) const;
  double getUpEdgeOfBinXContaining(double x) const;
  unsigned getIndexOfBinYContaining(double y) const;
  double getLowEdgeOfBinYContaining(double y) const;
  double getUpEdgeOfBinYContaining(double y) const;
  double getValueAt(double x, double y) const;
  double getIntegral() const;
  std::pair<unsigned, unsigned> getIndicesOfMaximum() const;
  double getMaximumValue() const;
  double getXRangeWidth() const;
  double getYRangeWidth() const;
  void getInterpolatedSliceAtX(double x, Hist& interpolated) const; //the result is stored in interpolated and assumes that interpolated.getNumberOfBins() = getNumberOfBinsY()
  void getInterpolatedSliceAtY(double y, Hist& interpolated) const;
  double getRandomFromInterpolatedSliceAtX(TRandom3& ran, double x) const;
  double getRandomFromInterpolatedSliceAtY(TRandom3& ran, double y) const;
  void setBinContent(double z, unsigned i, unsigned j);

};

std::ostream& operator<<(std::ostream& output, const Hist2d& h); 

#endif