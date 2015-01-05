#include "Detector.hpp"

ostream& operator<<(ostream& output, const Detector& detector){
  
  for(const Volume* vol: detector.volumes) output<<vol;
  return output;
  
}

Detector::Detector(){
  
}
