#ifndef HIST_H
#define HIST_H

#include <vector>
#include <string>
#include <algorithm>
#include <ostream>
#include <iomanip>
#include <TRandom3.h>

class Hist { 
  
  std::string name;
  std::vector<double> edge;
  std::vector<double> value;
  enum randomType{rejection, transform};
  void buildConstantBining(unsigned numberOfBins, double lowerEdge, double upperEdge);
  double getRandomFromRejection(TRandom3& ran) const;
  double getRandomFromTransform(TRandom3& ran) const;
  std::vector<double> getCumulativeIntegral() const;
  
public:
  Hist(const std::string& name, const std::vector<double>& edge);
  Hist(const std::string& name, const std::vector<double>& edge, const std::vector<double>& value);
  Hist(const std::string& name, unsigned numberOfBins, double lowerEdge, double upperEdge);
  const std::string& getName() const;
  const double& getBinLowEdge(unsigned i) const;
  const double& getBinUpEdge(unsigned i) const;
  double getBinWidth(unsigned i) const; 
  unsigned getNumberOfBins() const;
  const double& getBinContent(unsigned i) const;
  unsigned getIndexOfBinContaining(double x) const;
  double getLowEdgeOfBinContaining(double x) const;
  double getUpEdgeOfBinContaining(double x) const;
  double getValueAt(double x) const;
  double getIntegral() const;
  unsigned getIndexOfMaximum() const;
  double getMaximumValue() const;
  double getRangeWidth() const;
  double getRandom(TRandom3& ran, randomType choice = transform) const;
  void setBinContent(double y, unsigned i);

};

std::ostream& operator<<(std::ostream& output, const Hist& h); 

#endif