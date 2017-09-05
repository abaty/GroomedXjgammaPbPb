#ifndef MIXING
#define MIXING

#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <iostream>

class Mixing{
  public:

    inline void setTriggerName(std::string name);
    inline void setVzMatchingWindow(float window);   
    inline void setEvtPlaneMatchingWindow(float window);   
    inline void setMaxAttempts(int n);
 
    //scans until it finds an event with matching parameters;
    void getEvent(bool reset = true);
    void getEvent(int _hiBin, bool reset = true);
    void getEvent(int _hiBin, float _vz, bool reset = true); 
    void getEvent(int _hiBin, float _vz, float _evtPlane, bool reset = true);
    int getHibin();
    int getCurrentEvtIndx();


    Mixing(std::vector< std::string > fList);
    ~Mixing();
 
  private:
    std::string mixTriggerName = "HLT_HIL1MinimumBiasHF2AND_part1_v1";
    float vzMatchingWindow = 0.5;
    float evtPlaneMatchingWindow = 0.2;

    unsigned int currentFileIndx = 0;
    int currentEvtIndx = -1;
    int nAttempts = 0;
    int maximumAttempts = 10000000;
   
    TFile * f;
    void openFile(int indx);
    void closeFile();
    void getNextFile();

    TTree * hlt, * evt, * skim;
    void setTrees();
    void setBranches();
    
    //branch variables 
    int trigger, HBHE, pVtx, hfcoinc3, pCluster;
    int hiBin;
    float evtPlane[29];
    float vz;


    std::vector< std::string > mixFileList;
};

inline void Mixing::setTriggerName(std::string name){
  mixTriggerName = name;
}

inline void Mixing::setVzMatchingWindow(float window){
  vzMatchingWindow = window;
}

inline void Mixing::setEvtPlaneMatchingWindow(float window){
  evtPlaneMatchingWindow = window;
}   

inline void Mixing::setMaxAttempts(int n){
  maximumAttempts = n;
}

void Mixing::openFile(int indx){ 
  f = TFile::Open(mixFileList.at(indx).c_str(),"read");
}

void Mixing::closeFile(){ 
  f->Close();
}

void Mixing::getNextFile(){
  closeFile();
  currentFileIndx++;
  if(currentFileIndx == mixFileList.size()){
    currentFileIndx = 0;//wrap around
  }

  currentEvtIndx = 0;//reset evt indx
  openFile(currentFileIndx);
  setTrees();
  setBranches(); 
}

void Mixing::setTrees(){
  evt = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
  hlt = (TTree*)f->Get("hltanalysis/HltTree");
  skim = (TTree*)f->Get("skimanalysis/HltTree");
}

void Mixing::setBranches(){
  hlt->SetBranchAddress(mixTriggerName.c_str(),&trigger);
  evt->SetBranchAddress("hiBin",&hiBin);
  evt->SetBranchAddress("vz",&vz);
  evt->SetBranchAddress("hiEvtPlanes",&evtPlane);
  skim->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHE);
  skim->SetBranchAddress("phfCoincFilter3",&hfcoinc3);
  skim->SetBranchAddress("pclusterCompatibilityFilter",&pCluster);
  skim->SetBranchAddress("pprimaryVertexFilter",&pVtx);
}

void Mixing::getEvent(bool reset){
  if(reset) nAttempts = 0;
  while(true){
    nAttempts++;
    currentEvtIndx++;
    if(nAttempts > maximumAttempts) break;
    
    if(currentEvtIndx >= evt->GetEntries()){
      getNextFile();     
    }
    evt->GetEntry(currentEvtIndx);
    skim->GetEntry(currentEvtIndx); 
    if(!(trigger && HBHE && hfcoinc3 && pCluster && pVtx)) continue;
    evt->GetEntry(currentEvtIndx); 
  
    break; 
  }
}

void Mixing::getEvent(int _hiBin, bool reset){
  if(reset) nAttempts = 0;
  while(true){
    getEvent(false);
    if(nAttempts > maximumAttempts) break;
    if(hiBin != _hiBin) continue;
    break;
  }
}

void Mixing::getEvent(int _hiBin, float _vz, bool reset){
  if(reset) nAttempts = 0;
  while(true){
    getEvent(_hiBin, false);
    if(nAttempts > maximumAttempts) break;
    if(TMath::Abs(vz-_vz) > vzMatchingWindow) continue;
    break;
  }
}

void Mixing::getEvent(int _hiBin, float _vz, float _evtPlane, bool reset){
  if(reset) nAttempts = 0;
  while(true){
    getEvent(_hiBin, _vz, false);
    if(nAttempts > maximumAttempts) break;
    if(TMath::ACos(TMath::Cos(evtPlane[8]-_evtPlane)) > evtPlaneMatchingWindow) continue;
    break;
  }
}

int Mixing::getHibin(){
  return hiBin;
}

int Mixing::getCurrentEvtIndx(){
  return currentEvtIndx;
}

Mixing::Mixing(std::vector< std::string >  fList){
  std::cout << "Initializing Mixing object!" << std::endl; 
  mixFileList = fList;
  if(mixFileList.size() > 0){
    openFile(0);
    setTrees();
    setBranches();
  }
  else{
    std::cout << "empty mixing file list!" << std::endl;
  }
}

Mixing::~Mixing(){
  closeFile();
}
#endif
