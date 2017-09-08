#ifndef SKIMSETTINGS
#define SKIMSETTINGS

#include <string>
#include <iostream>
#include "TH1F.h"
#include "TMutex.h"

class SkimSettings{
  public:

  float phoEtaCut = 1.44;
  float phoPtCut = 40;

  std::string jetTree = "akCs4PFJetAnalyzer";
  float jetEtaCut = 1.6;
  float jetPtCut = 15;
  float dPhiCut = 2.74889357; //7pi/8

  float vzCut = 15;
  std::string trigger = "HLT_HISinglePhoton40_Eta1p5_v1";

  float groomedJetMatchingCut = 0.35;
  static const int nSubJetTrees = 2;
  std::string subJetTreeNames[nSubJetTrees] = {"akCsSoftDrop4PFJetAnalyzer","akCsSoftDropZ05B154PFJetAnalyzer"};
  std::string subJetoutBranchNames[nSubJetTrees] = {"Z0B0p1","Z0p5B1p5"};

  int nEvts = -1;
  int nMixEvts = 1;

  SkimSettings();
  
  private:
  
};
  
SkimSettings::SkimSettings(){
  std::cout << "Initializing Settings!" << std::endl; 
}
#endif
