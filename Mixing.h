#ifndef MIXING
#define MIXING

#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <iostream>
#include "Tools.h" // needed for dR function

class Mixing{
  public:

    inline void setTriggerNames(std::vector< std::string >  names);
    inline void setJetCollection(std::string name);
    inline void setVzMatchingWindow(float window);   
    inline void setEvtPlaneMatchingWindow(float window);   
    inline void setMaxAttempts(int n);
    inline void setSubjetCollections(std::vector<std::string> vec); 

    //scans until it finds an event with matching parameters;
    void getEvent(bool reset = true);
    void getEvent(int _hiBin, bool reset = true);
    void getEvent(int _hiBin, float _vz, bool reset = true); 
    void getEvent(int _hiBin, float _vz, float _evtPlane, bool reset = true);
    int getHibin();
    int getCurrentEvtIndx();

    void getBack2BackJets(std::vector< float > &pts, float etaCut, float phoPhi, float dPhiCut, float ptCut);
    void getSubjets(std::vector< float>  &mixdR12, int nSubjetTree, float matchingCut, float etaCut, float phoPhi, float dPhiCut, float ptCut); 

    Mixing(std::vector< std::string > fList, int startingIndx = 0);
    ~Mixing();
 
  private:
    std::vector< std::string > mixTriggerNames;
    std::string jetCollection = "akCs4PFJetAnalyzer";
    std::vector< std::string > subjetCollections;
    float vzMatchingWindow = 5;
    float evtPlaneMatchingWindow = 0.2;

    unsigned int currentFileIndx = 0;
    int currentEvtIndx = -1;
    int nAttempts = 0;
    int maximumAttempts = 10000000;
    int nSubjectTrees = 0;
   
    TFile * f;
    void openFile(int indx);
    void closeFile();
    void getNextFile();

    TTree * hlt, * evt, * skim, * jet, * subjet;
    void setTrees();
    void setBranches();
    
    //branch variables
    int triggerList[10] = {0}; 
    int trigger, HBHE, pVtx, hfcoinc3, pCluster;
    int hiBin;
    float evtPlane[29];
    float vz;
    int nref;
    float jetPt[1000];
    float jetEta[1000];
    float jetPhi[1000];

    int groomnref;
    float groomjtphi[1000];
    float groomjteta[1000];
    std::vector< std::vector< float >> * subjetEta = 0;
    std::vector< std::vector< float >> * subjetPhi = 0;
    std::vector< std::vector< float >> * subjetPt = 0;

    std::vector< std::string > mixFileList;
};

inline void Mixing::setTriggerNames(std::vector< std::string >  names){
  if(names.size()>10) std::cout << "Warning, more than 10 triggers specified, increase size of array in Mixing class!" << std::endl;
  mixTriggerNames = names;
}

inline void Mixing::setJetCollection(std::string name){
  jetCollection = name;
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

inline void Mixing::setSubjetCollections(std::vector<std::string> vec){
  subjetCollections = vec;
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
  jet = (TTree*)f->Get(Form("%s/t",jetCollection.c_str())); 
}

void Mixing::setBranches(){
  for(unsigned int i = 0; i<mixTriggerNames.size(); i++) hlt->SetBranchAddress(mixTriggerNames.at(i).c_str(),&(triggerList[i]));
  evt->SetBranchAddress("hiBin",&hiBin);
  evt->SetBranchAddress("vz",&vz);
  evt->SetBranchAddress("hiEvtPlanes",&evtPlane);
  skim->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHE);
  skim->SetBranchAddress("phfCoincFilter3",&hfcoinc3);
  skim->SetBranchAddress("pclusterCompatibilityFilter",&pCluster);
  skim->SetBranchAddress("pprimaryVertexFilter",&pVtx);

  jet->SetBranchAddress("nref",&nref);
  jet->SetBranchAddress("jtpt",&jetPt);
  jet->SetBranchAddress("jteta",&jetEta);
  jet->SetBranchAddress("jtphi",&jetPhi);
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
    hlt->GetEntry(currentEvtIndx);
    trigger = 0;
    for(unsigned int i = 0; i<mixTriggerNames.size(); i++) trigger = trigger || triggerList[i];
    skim->GetEntry(currentEvtIndx); 
    if(!(trigger && HBHE && hfcoinc3 && pCluster && pVtx)) continue;
    evt->GetEntry(currentEvtIndx); 
    if(TMath::Abs(vz)>15) continue;
  
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

void Mixing::getBack2BackJets(std::vector< float > &pts, float etaCut, float phoPhi, float dPhiCut, float ptCut){
  jet->GetEntry(currentEvtIndx);
  for(int i = 0; i<nref; i++){
    if(jetPt[i]<ptCut) continue;
    if(TMath::Abs(jetEta[i])>etaCut) continue;
    if(TMath::ACos(TMath::Cos(jetPhi[i]-phoPhi)) < dPhiCut) continue;
    pts.push_back(jetPt[i]);
  }
}

void Mixing::getSubjets(std::vector< float>&  mixdR12, int  nSubjetTree, float matchingCut, float etaCut, float phoPhi, float dPhiCut, float ptCut){
  jet->GetEntry(currentEvtIndx); 
  subjet = (TTree*) f->Get(Form("%s/t",subjetCollections.at(nSubjetTree).c_str()));
  subjet->SetBranchAddress("nref",&groomnref);
  subjet->SetBranchAddress("jtphi",&groomjtphi);
  subjet->SetBranchAddress("jteta",&groomjteta);
  subjet->SetBranchAddress("jtSubJetEta",&subjetEta);
  subjet->SetBranchAddress("jtSubJetPhi",&subjetPhi);
  subjet->SetBranchAddress("jtSubJetPt",&subjetPt);
  subjet->GetEntry(currentEvtIndx); 
 
  for(int k = 0; k<nref; k++){
    if(jetPt[k]<ptCut) continue;
    if(TMath::Abs(jetEta[k])>etaCut) continue;
    if(TMath::ACos(TMath::Cos(jetPhi[k]-phoPhi)) < dPhiCut) continue;
    bool isMatched = false;

    for(int m = 0; m<groomnref; m++){
      if( dR(jetEta[k], jetPhi[k], groomjteta[m], groomjtphi[m]) > matchingCut) continue;
      isMatched = true;
      if(subjetEta->at(m).size()==2){
        mixdR12.push_back(dR(subjetEta->at(m).at(0),subjetPhi->at(m).at(0),subjetEta->at(m).at(1),subjetPhi->at(m).at(1)));
      }else{
        mixdR12.push_back(-999); 
      }
      break;      
    }
    if(!isMatched){//fix by lowering pt cut on forest
      mixdR12.push_back(-999); 
    }
  }
} 

Mixing::Mixing(std::vector< std::string >  fList, int startingIndx){
  std::cout << "Initializing Mixing object!" << std::endl; 
  mixFileList = fList;
  if(mixFileList.size() > 0){
    openFile(startingIndx);
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
