#include <thread>
#include "TFile.h"
#include "Constants.hpp"
#include "Event.hpp"

using namespace constants;

struct hepGeneratorLi{
  
  const char* outfilename;
  const unsigned nEvents;
  const vector<double> rb; 
  const vector<double> rs;
  const vector<Hist2d> betaMap;
  
  hepGeneratorLi(const char* outfilename, const unsigned nEvents, const vector<double>& rb, const vector<double>& rs, const vector<Hist2d>& betaMap):outfilename(outfilename),nEvents(nEvents),rb(rb),rs(rs), betaMap(betaMap) {}
  void operator()(){
    
  //Li9 tree 
    State* electron = new RState(mass::_e, PDG::_e);
    State* neutron  = new RState(mass::_n, PDG::_n);
    State* alpha  = new RState(mass::_alpha, PDG::_alpha);
    State* alpha_bis = new RState(mass::_alpha, PDG::_alpha);
    
    State* sa_n = new VState(mass::_5He, {alpha_bis, neutron}); 
    State* sa_n_374 = new VState(mass::_5He374, {alpha_bis, neutron});
    State* sa_a = new VState(mass::_8Be, {alpha, alpha_bis});
    State* sa_a_303 = new VState(mass::_8Be303, {alpha, alpha_bis});
    State* sa_a_1135 = new VState(mass::_8Be1135, {alpha, alpha_bis});
    
    State* s5He = new GState(mass::_5He, PDG::_5He, {sa_n}, {1}, width::_5He);
    State* s5He374 = new GState(mass::_5He374, PDG::_5He, {sa_n_374}, {1}, width::_5He374);
    State* s8Be = new GState(mass::_8Be, PDG::_8Be, {sa_a}, {1}, width::_8Be);
    State* s8Be303 = new GState(mass::_8Be303, PDG::_8Be, {sa_a_303}, {1}, width::_8Be303);
    State* s8Be1135 = new GState(mass::_8Be1135, PDG::_8Be, {sa_a_1135}, {1}, width::_8Be1135);
    
    State* saan_243 = new VState(mass::_9Be243, {alpha, alpha_bis, neutron});
    State* s5He_alpha_243 = new VState(mass::_9Be243, {s5He, alpha});
    State* s8Be_n_243 = new VState(mass::_9Be243, {s8Be, neutron});
    
    State* saan_278 = new VState(mass::_9Be278, {alpha, alpha_bis, neutron});
    State* s5He_alpha_278 = new VState(mass::_9Be278, {s5He, alpha});
    State* s8Be_n_278 = new VState(mass::_9Be278, {s8Be, neutron});
    State* s8Be303_n_278 = new VState(mass::_9Be278, {s8Be303, neutron});
    
    State* saan_794 = new VState(mass::_9Be794, {alpha, alpha_bis, neutron});
    State* s5He_alpha_794 = new VState(mass::_9Be794, {s5He, alpha});
    State* s5He374_alpha_794 = new VState(mass::_9Be794, {s5He374, alpha});
    State* s8Be_n_794 = new VState(mass::_9Be794, {s8Be, neutron});
    State* s8Be303_n_794 = new VState(mass::_9Be794, {s8Be303, neutron});

    State* saan_1128 = new VState(mass::_9Be1128, {alpha, alpha_bis, neutron});
    State* s5He_alpha_1128 = new VState(mass::_9Be1128, {s5He, alpha});
    State* s5He374_alpha_1128 = new VState(mass::_9Be1128, {s5He374, alpha});
    State* s8Be_n_1128 = new VState(mass::_9Be1128, {s8Be, neutron});
    State* s8Be303_n_1128 = new VState(mass::_9Be1128, {s8Be303, neutron});
    State* s8Be1135_n_1128 = new VState(mass::_9Be1128, {s8Be1135, neutron});
    
    State* saan_1181 = new VState(mass::_9Be1181, {alpha, alpha_bis, neutron});
    State* s5He_alpha_1181 = new VState(mass::_9Be1181, {s5He, alpha});
    State* s5He374_alpha_1181 = new VState(mass::_9Be1181, {s5He374, alpha});
    State* s8Be_n_1181 = new VState(mass::_9Be1181, {s8Be, neutron});
    State* s8Be303_n_1181 = new VState(mass::_9Be1181, {s8Be303, neutron});
    State* s8Be1135_n_1181 = new VState(mass::_9Be1181, {s8Be1135, neutron});

    State* s9Be243 = new GState(mass::_9Be243, PDG::_9Be, {saan_243, s8Be_n_243, s5He_alpha_243}, rs, width::_9Be243);
    State* s9Be278 = new GState(mass::_9Be278, PDG::_9Be, {saan_278, s8Be_n_278, s5He_alpha_278, s8Be303_n_278}, rs, width::_9Be278);
    State* s9Be794 = new GState(mass::_9Be794, PDG::_9Be, {saan_794, s8Be_n_794, s5He_alpha_794, s8Be303_n_794, s5He374_alpha_794}, rs, width::_9Be794);
    State* s9Be1128 = new GState(mass::_9Be1128, PDG::_9Be, {saan_1128, s8Be_n_1128, s5He_alpha_1128, s8Be303_n_1128, s5He374_alpha_1128, s8Be1135_n_1128}, rs, width::_9Be1128);
    State* s9Be1181 = new GState(mass::_9Be1181, PDG::_9Be, {saan_1181, s8Be_n_1181, s5He_alpha_1181, s8Be303_n_1181, s5He374_alpha_1181, s8Be1135_n_1181}, rs, width::_9Be1181);
    
    State* s9Be243_e = new BState(mass::_9Li, {s9Be243, electron}, betaMap[0]);
    State* s9Be278_e = new BState(mass::_9Li, {s9Be278, electron}, betaMap[1]);
    State* s9Be794_e = new BState(mass::_9Li, {s9Be794, electron}, betaMap[2]);
    State* s9Be1128_e = new BState(mass::_9Li, {s9Be1128, electron}, betaMap[3]);
    State* s9Be1181_e = new BState(mass::_9Li, {s9Be1181, electron}, betaMap[4]);

    State* s9Li = new RState(mass::_9Li, PDG::_9Li, {s9Be243_e, s9Be278_e, s9Be794_e, s9Be1128_e, s9Be1181_e}, rb);
    
    ofstream out(outfilename); //stream to write the hepevevent text file
    TRandom3 ran(0);
    
    Event lithium(s9Li, dchooz::R_target, dchooz::h_target);

    for(unsigned i = 0; i<nEvents; ++i){

      if(ran.Rndm()<dchooz::probCosmogenicInNT) lithium.SetInCylinder(dchooz::R_target, dchooz::h_target);
      else lithium.SetInHollowCylinder(dchooz::R_target, dchooz::R_GC, dchooz::h_target, dchooz::h_GC);
      lithium.Process();
      out<<lithium;

    }

    out.close();

    //manually release the 9Li tree
    delete electron;delete neutron;delete alpha;delete alpha_bis;
    delete sa_n;delete sa_n_374;delete sa_a;delete sa_a_303;delete sa_a_1135;
    delete s5He;delete s5He374;delete s8Be;delete s8Be303;delete s8Be1135;
    delete saan_243;delete s5He_alpha_243;delete s8Be_n_243;
    delete saan_278;delete s5He_alpha_278;delete s8Be_n_278;delete s8Be303_n_278;
    delete saan_794;delete s5He_alpha_794;delete s5He374_alpha_794;delete s8Be_n_794;delete s8Be303_n_794;
    delete saan_1128;delete s5He_alpha_1128;delete s5He374_alpha_1128;delete s8Be_n_1128;delete s8Be303_n_1128;delete s8Be1135_n_1128;
    delete saan_1181;delete s5He_alpha_1181;delete s5He374_alpha_1181;delete s8Be_n_1181;delete s8Be303_n_1181;delete s8Be1135_n_1181;
    delete s9Be243;delete s9Be278;delete s9Be794;delete s9Be1128;delete s9Be1181;
    delete s9Be243_e;delete s9Be278_e;delete s9Be794_e;delete s9Be1128_e;delete s9Be1181_e;
    delete s9Li;

  };
  
};

struct hepGeneratorHe{
  
  const char* outfilename;
  const unsigned nEvents;
  const vector<double> rb; 
  const vector<double> rs;
  vector<Hist2d> betaMap;
  
  hepGeneratorHe(const char* outfilename, const unsigned nEvents, const vector<double>& rb, const vector<double>& rs, const vector<Hist2d>& betaMap):outfilename(outfilename),nEvents(nEvents),rb(rb),rs(rs), betaMap(betaMap) {}
  void operator()(){
    
  //Li9 tree 
    State* electron = new RState(mass::_e, PDG::_e);
    State* neutron  = new RState(mass::_n, PDG::_n);
    State* alpha  = new RState(mass::_alpha, PDG::_alpha);
    State* alpha_bis = new RState(mass::_alpha, PDG::_alpha);
    
    State* sa_n = new VState(mass::_5He, {alpha_bis, neutron}); 
    
    State* s5He = new GState(mass::_5He, PDG::_5He, {sa_n}, {1}, width::_5He);

  //He8 Tree 
    State* gamma = new RState(mass::_gamma, PDG::_gamma);
    State* triton = new RState(mass::_3H, PDG::_3H);
    State* s7Li = new RState(mass::_7Li, PDG::_7Li);

    State* s7Li_gamma_0748 = new VState(mass::_7Li048, {s7Li, gamma});

    State* s7Li0748 = new RState(mass::_7Li048, PDG::_7Li, {s7Li_gamma_0748});

    State* s7Li_n_321 = new VState(mass::_8Li321, {s7Li, neutron});
    State* s7Li0748_n_321 = new VState(mass::_8Li321, {s7Li0748, neutron});
    
    State* s7Li_n_54 = new VState(mass::_8Li54, {s7Li, neutron});
    State* s7Li0748_n_54 = new VState(mass::_8Li54, {s7Li0748, neutron});
    State* sa_t_n_54 = new VState(mass::_8Li54, {alpha, triton, neutron});
    
    State* s7Li_n_967 = new VState(mass::_8Li967, {s7Li, neutron});
    State* s7Li0748_n_967 = new VState(mass::_8Li967, {s7Li0748, neutron});
    State* sa_t_n_967 = new VState(mass::_8Li967, {alpha, triton, neutron});
    State* s5He_t_967 = new VState(mass::_8Li967, {s5He, triton});
    
    State* s8Li321 = new GState(mass::_8Li321, PDG::_8Li, {s7Li_n_321, s7Li0748_n_321}, rs, width::_8Li321);
    State* s8Li54 = new GState(mass::_8Li54, PDG::_8Li, {s7Li_n_54, s7Li0748_n_54, sa_t_n_54}, rs, width::_8Li54);
    State* s8Li967 = new GState(mass::_8Li967, PDG::_8Li, {s7Li_n_967, s7Li0748_n_967, sa_t_n_967, s5He_t_967}, rs, width::_8Li967);
    
    State* s8Li321_e = new BState(mass::_9Li, {s8Li321, electron}, betaMap[0]);
    State* s8Li54_e = new BState(mass::_9Li, {s8Li54, electron}, betaMap[1]);
    State* s8Li967_e = new BState(mass::_9Li, {s8Li967, electron}, betaMap[2]);
    
    State* s8He = new RState(mass::_8He, PDG::_8He, {s8Li321_e, s8Li54_e, s8Li967_e}, rb);
    
    ofstream out(outfilename); //stream to write the hepevevent text file
    TRandom3 ran(0);

    Event helium(s8He, dchooz::R_target, dchooz::h_target);
      
    for(unsigned i = 0; i<nEvents; ++i){
      
      if(ran.Rndm()<dchooz::probCosmogenicInNT) helium.SetInCylinder(dchooz::R_target, dchooz::h_target); //set the position member to a random point in a cylinder of radius R and height h
      else helium.SetInHollowCylinder(dchooz::R_target, dchooz::R_GC, dchooz::h_target, dchooz::h_GC);
      helium.Process();
      out<<helium;

    }
    
    out.close();
    
    //manually release the 8He tree
    delete electron; delete neutron; delete alpha; delete alpha_bis;
    delete sa_n;
    delete s5He;
    delete gamma;delete triton;delete s7Li;
    delete s7Li_gamma_0748;
    delete s7Li0748;
    delete s7Li_n_321;delete s7Li0748_n_321;
    delete s7Li_n_54;delete s7Li0748_n_54;delete sa_t_n_54;
    delete s7Li_n_967;delete s7Li0748_n_967;delete sa_t_n_967;delete s5He_t_967;
    delete s8Li321;delete s8Li54;delete s8Li967;
    delete s8Li321_e;delete s8Li54_e;delete s8Li967_e;
    delete s8He;

  };
  
};

template <class T> void launchThreads(vector<T>& workers){

  unsigned nThreads = thread::hardware_concurrency();//get the number of working cores
  if(nThreads == 0) nThreads = 1;
  vector<thread> threads(nThreads);
  
  auto itWk = workers.begin();
  while(itWk != workers.end()){
    
    auto itTh = threads.begin(); 
    while(itTh != threads.end() && itWk != workers.end()){

      *itTh = thread(*itWk);//start each thread
      ++itTh;
      ++itWk;
      
    }
    
    for(thread& th : threads) if(th.joinable()) th.join();//join them all to the current thread
    
  }
  
}

void generateLiBranches(const unsigned nEvents = 1e4){
    
  vector<hepGeneratorLi> workers;//speed up the branch generations with a thread per branch
    
  vector<Hist2d> betaMapLi (5);
  TFile LiFile("spectres_beta/cartes_li9.root");//TFile's are not thread safe
  betaMapLi[0] = Hist2d(static_cast<TH2F*>(LiFile.Get("carte_0")));
  betaMapLi[1] = Hist2d(static_cast<TH2F*>(LiFile.Get("carte_1")));
  betaMapLi[2] = Hist2d(static_cast<TH2F*>(LiFile.Get("carte_2")));
  betaMapLi[3] = Hist2d(static_cast<TH2F*>(LiFile.Get("carte_3")));
  betaMapLi[4] = Hist2d(static_cast<TH2F*>(LiFile.Get("carte_4")));
  LiFile.Close();
  
  vector<double> rbLi {1, 0, 0, 0, 0};
  vector<double> rsLi {1, 0, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_243_aan.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 1, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_243_8Be.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 1, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_243_5He.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  
  rbLi = {0, 1, 0, 0, 0};
  rsLi = {1, 0, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_278_aan.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 1, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_278_8Be.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 1, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_278_5He.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 0, 1, 0, 0};
  workers.push_back(hepGeneratorLi("li_278_8Be303.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  
  rbLi = {0, 0, 1, 0, 0};
  rsLi = {1, 0, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_794_aan.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 1, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_794_8Be.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 1, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_794_5He.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 0, 1, 0, 0};
  workers.push_back(hepGeneratorLi("li_794_8Be303.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 0, 0, 1, 0};
  workers.push_back(hepGeneratorLi("li_794_5He374.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  
  rbLi = {0, 0, 0, 1, 0};
  rsLi = {1, 0, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_1128_aan.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 1, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_1128_8Be.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 1, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_1128_5He.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 0, 1, 0, 0};
  workers.push_back(hepGeneratorLi("li_1128_8Be303.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 0, 0, 1, 0};
  workers.push_back(hepGeneratorLi("li_1128_5He374.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 0, 0, 0, 1};
  workers.push_back(hepGeneratorLi("li_1128_8Be1135.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  
  rbLi = {0, 0, 0, 0, 1};
  rsLi = {1, 0, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_1181_aan.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 1, 0, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_1181_8Be.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 1, 0, 0, 0};
  workers.push_back(hepGeneratorLi("li_1181_5He.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 0, 1, 0, 0};
  workers.push_back(hepGeneratorLi("li_1181_8Be303.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 0, 0, 1, 0};
  workers.push_back(hepGeneratorLi("li_1181_5He374.txt", nEvents, rbLi, rsLi,betaMapLi)); 
  rsLi = {0, 0, 0, 0, 0, 1};
  workers.push_back(hepGeneratorLi("li_1181_8Be1135.txt", nEvents, rbLi, rsLi,betaMapLi)); 

  launchThreads(workers);
  
}

void generateHeBranches(const unsigned nEvents = 1e4){
    
  vector<hepGeneratorHe> workers;//speed up the branch generations with a thread per branch
    
  vector<Hist2d> betaMapHe (3);
  TFile HeFile("spectres_beta/cartes_he8.root");//TFile's are not thread safe
  betaMapHe[0] = Hist2d(static_cast<TH2F*>(HeFile.Get("carte_0")));
  betaMapHe[1] = Hist2d(static_cast<TH2F*>(HeFile.Get("carte_1")));
  betaMapHe[2] = Hist2d(static_cast<TH2F*>(HeFile.Get("carte_2")));
  HeFile.Close();
  
  vector<double> rbHe {1, 0, 0};
  vector<double> rsHe {1, 0, 0, 0};
  workers.push_back(hepGeneratorHe("he_321_7Li.txt", nEvents, rbHe, rsHe,betaMapHe)); 
  rsHe = {0, 1, 0, 0};
  workers.push_back(hepGeneratorHe("he_321_7Li048.txt", nEvents, rbHe, rsHe,betaMapHe)); 
  
  rbHe = {0, 1, 0};
  rsHe = {1, 0, 0, 0};
  workers.push_back(hepGeneratorHe("he_54_7Li.txt", nEvents, rbHe, rsHe,betaMapHe)); 
  rsHe = {0, 1, 0, 0};
  workers.push_back(hepGeneratorHe("he_54_7Li048.txt", nEvents, rbHe, rsHe,betaMapHe)); 
  rsHe = {0, 0, 1, 0};
  workers.push_back(hepGeneratorHe("he_54_atn.txt", nEvents, rbHe, rsHe,betaMapHe)); 
  
  rbHe = {0, 0, 1};
  rsHe = {1, 0, 0, 0};
  workers.push_back(hepGeneratorHe("he_967_7Li.txt", nEvents, rbHe, rsHe,betaMapHe)); 
  rsHe = {0, 1, 0, 0};
  workers.push_back(hepGeneratorHe("he_967_7Li048.txt", nEvents, rbHe, rsHe,betaMapHe)); 
  rsHe = {0, 0, 1, 0};
  workers.push_back(hepGeneratorHe("he_967_atn.txt", nEvents, rbHe, rsHe,betaMapHe)); 
  rsHe = {0, 0, 0, 1};
  workers.push_back(hepGeneratorHe("he_967_5He.txt", nEvents, rbHe, rsHe,betaMapHe)); 

  launchThreads(workers);
  
}

int main(int argc, char* argv[]){
  
  if(argc == 2) generateLiBranches(stod(argv[1]));
  else if(argc == 3 && string(argv[2]) == "helium") generateHeBranches(stod(argv[1]));
  else generateLiBranches();
  return 0;

}
