#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "skimSettings.h"
#include "Tools.h"

void gammaJetSkim(std::vector< std::string > inputFiles, int job){
  SkimSettings s = SkimSettings();

  int trig = 0;
  int hfCoinc = 0; 
  int HBHE = 0;
  int pVtx = 0;
  int clusterCompat = 0;
  
  unsigned int run = 0;
  ULong64_t evtNum = 0;
  int hiBin = -1;
  float hiEvtPlanes[29];
  float vz = -99;
  float evtPlane;

  std::vector< float > * seedTime = 0;
  std::vector< float > * swissCrx = 0;
  std::vector< float > * sigmaIEtaIEta = 0;
  std::vector< float > * ecalR4 = 0;
  std::vector< float > * hcalR4 = 0;
  std::vector< float > * trkR4 = 0;
  std::vector< float > * hOverE = 0;
  std::vector< float > * phoEt = 0;
  std::vector< float > * phoPhi = 0;
  std::vector< float > * phoEta = 0;
  float gammaPt, gammaEta, gammaPhi, gammaSigIEIE;

  int nref;
  float jtpt[1000];
  float jtphi[1000];
  float jteta[1000];
  std::vector<float> jetPt;
  std::vector<float> jetEta;
  std::vector<float> jetPhi;

  int groomnref;
  float groomjtphi[1000];
  float groomjteta[1000];
  std::vector< std::vector< float >> * subjetEta = 0;
  std::vector< std::vector< float >> * subjetPhi = 0;
  std::vector< std::vector< float >> * subjetPt = 0;
  std::vector< float > subjet1Pt[s.nSubJetTrees];//output vectors
  std::vector< float > subjet1Eta[s.nSubJetTrees];
  std::vector< float > subjet1Phi[s.nSubJetTrees];
  std::vector< float > subjet2Pt[s.nSubJetTrees];
  std::vector< float > subjet2Eta[s.nSubJetTrees];
  std::vector< float > subjet2Phi[s.nSubJetTrees];
  std::vector< float > dR12[s.nSubJetTrees];

  TFile * output = TFile::Open(Form("output_%d.root",job),"recreate");
  TTree * outTree = new TTree("skim","Photon+GroomedJetSubStructure skim");
  outTree->Branch("run",&run);
  outTree->Branch("evt",&evtNum);
  outTree->Branch("vz",&vz);
  outTree->Branch("hiBin",&hiBin);
  outTree->Branch("evtPlane",&evtPlane);
  outTree->Branch("phoPt",&gammaPt);
  outTree->Branch("phoEta",&gammaEta);
  outTree->Branch("phoPhi",&gammaPhi);
  outTree->Branch("phoSigIEIE",&gammaSigIEIE);
  outTree->Branch("jetPt",&jetPt);
  outTree->Branch("jetEta",&jetEta);
  outTree->Branch("jetPhi",&jetPhi);
  for(int i = 0; i<s.nSubJetTrees; i++){
    outTree->Branch(Form("subjet1Pt_%s",s.subJetoutBranchNames[i].c_str()),&(subjet1Pt[i]));
    outTree->Branch(Form("subjet1Eta_%s",s.subJetoutBranchNames[i].c_str()),&(subjet1Eta[i]));
    outTree->Branch(Form("subjet1Phi_%s",s.subJetoutBranchNames[i].c_str()),&(subjet1Phi[i]));
    outTree->Branch(Form("subjet2Pt_%s",s.subJetoutBranchNames[i].c_str()),&(subjet2Pt[i]));
    outTree->Branch(Form("subjet2Eta_%s",s.subJetoutBranchNames[i].c_str()),&(subjet2Eta[i]));
    outTree->Branch(Form("subjet2Phi_%s",s.subJetoutBranchNames[i].c_str()),&(subjet2Phi[i]));
    outTree->Branch(Form("dR12_%s",s.subJetoutBranchNames[i].c_str()),&(dR12[i]));
  }

  for(unsigned int file = 0; file<inputFiles.size(); file++){
    TFile * f = TFile::Open(inputFiles.at(file).c_str(),"read");
    std::cout << "File: " << file << "/" << inputFiles.size() << std::endl;
    
    TTree * hlt = (TTree*)f->Get("hltanalysis/HltTree"); 
    hlt->SetBranchAddress(s.trigger.c_str(),&trig);

    TTree * skim = (TTree*)f->Get("skimanalysis/HltTree");
    skim->SetBranchAddress("phfCoincFilter3",&hfCoinc);
    skim->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHE);
    skim->SetBranchAddress("pprimaryVertexFilter",&pVtx); 
    skim->SetBranchAddress("pclusterCompatibilityFilter",&clusterCompat); 
    
    TTree * evt = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
    evt->SetBranchAddress("run",&run);
    evt->SetBranchAddress("evt",&evtNum);
    evt->SetBranchAddress("vz",&vz);
    evt->SetBranchAddress("hiBin",&hiBin);
    evt->SetBranchAddress("hiEvtPlanes",&hiEvtPlanes);//use the [8] element (9th in the list, index 8)
    hlt->AddFriend(skim);

    TTree * photons = (TTree*)f->Get("ggHiNtuplizer/EventTree");
    photons->SetBranchAddress("pho_seedTime",&seedTime);
    photons->SetBranchAddress("pho_swissCrx",&swissCrx);
    photons->SetBranchAddress("phoSigmaIEtaIEta",&sigmaIEtaIEta);
    photons->SetBranchAddress("pho_ecalClusterIsoR4",&ecalR4);
    photons->SetBranchAddress("pho_hcalRechitIsoR4",&hcalR4);
    photons->SetBranchAddress("pho_trackIsoR4PtCut20",&trkR4);
    photons->SetBranchAddress("phoHoverE",&hOverE);
    photons->SetBranchAddress("phoEt",&phoEt);
    photons->SetBranchAddress("phoPhi",&phoPhi);
    photons->SetBranchAddress("phoEta",&phoEta);

    TTree * jets = (TTree*)f->Get(Form("%s/t",s.jetTree.c_str()));
    jets->SetBranchAddress("nref",&nref);
    jets->SetBranchAddress("jtpt",&jtpt);
    jets->SetBranchAddress("jtphi",&jtphi);
    jets->SetBranchAddress("jteta",&jteta);

    TTree * subjets[s.nSubJetTrees];
    for(int k = 0; k<s.nSubJetTrees; k++) subjets[k] = (TTree*)f->Get(Form("%s/t",s.subJetTreeNames[k].c_str()));

    //event loop
    for(int i = 0; i<((s.nEvts<0)?hlt->GetEntries():(TMath::Min(s.nEvts,(int)hlt->GetEntries()))); i++){
      hlt->GetEntry(i);
      if(!trig) continue;
      if(!(hfCoinc && HBHE && pVtx && clusterCompat)) continue;

      evt->GetEntry(i);
      if(!(TMath::Abs(vz)<s.vzCut)) continue;
      evtPlane = hiEvtPlanes[8];

      photons->GetEntry(i);
      int leadPhoIndx = -1;
      for(unsigned int j = 0; j<phoEt->size(); j++){
        if( TMath::Abs(phoEta->at(j)) > s.phoEtaCut ) continue;
        if(!(TMath::Abs(seedTime->at(j))<3 && swissCrx->at(j)<0.9 && sigmaIEtaIEta->at(j)>0.002)) continue;//spike rejection    
        if(leadPhoIndx == -1) leadPhoIndx = j;
        else if( phoEt->at(j) > phoEt->at(leadPhoIndx) ) leadPhoIndx = j;
      }     
      if(leadPhoIndx==-1) continue;//no leading photon that was good
      if(phoEt->at(leadPhoIndx) < s.phoPtCut) continue; 
      if(!((ecalR4->at(leadPhoIndx)+hcalR4->at(leadPhoIndx)+trkR4->at(leadPhoIndx)) < 1 && (hOverE->at(leadPhoIndx)<0.1))) continue;//isolation of leading
      //leading not in the signal or sideband region 
      if(!(sigmaIEtaIEta->at(leadPhoIndx)<0.01 || (sigmaIEtaIEta->at(leadPhoIndx)>0.011 && sigmaIEtaIEta->at(leadPhoIndx)<0.017 ))) continue;
      gammaPt = phoEt->at(leadPhoIndx);
      gammaEta = phoEta->at(leadPhoIndx);
      gammaPhi = phoPhi->at(leadPhoIndx);
      gammaSigIEIE = sigmaIEtaIEta->at(leadPhoIndx);

      jets->GetEntry(i);
      jetPt.clear();
      jetEta.clear();
      jetPhi.clear();
      for(int j = 0; j<nref; j++){
        if(jtpt[j] < s.jetPtCut) break;//pt cut
        if(TMath::Abs(jteta[j]) > s.jetEtaCut) continue; //eta cut
        if(dPhi( jtphi[j] , phoPhi->at(leadPhoIndx)) < s.dPhiCut) continue; //dphi cut
        jetPt.push_back(jtpt[j]);
        jetEta.push_back(jteta[j]);
        jetPhi.push_back(jtphi[j]);
      }
      if(jetPt.size()==0) continue;

      for(int j = 0; j<s.nSubJetTrees; j++){
        subjets[j]->SetBranchAddress("nref",&groomnref);
        subjets[j]->SetBranchAddress("jtphi",&groomjtphi);
        subjets[j]->SetBranchAddress("jteta",&groomjteta);
        subjets[j]->SetBranchAddress("jtSubJetEta",&subjetEta);
        subjets[j]->SetBranchAddress("jtSubJetPhi",&subjetPhi);
        subjets[j]->SetBranchAddress("jtSubJetPt",&subjetPt);
        subjets[j]->GetEntry(i);
          
        subjet1Pt[j].clear();
        subjet1Eta[j].clear();
        subjet1Phi[j].clear();
        subjet2Pt[j].clear();
        subjet2Eta[j].clear();
        subjet2Phi[j].clear();
        dR12[j].clear();

        for(unsigned int k = 0; k<jetEta.size(); k++){
          bool isMatched = false;
          for(int m = 0; m<groomnref; m++){
            if( dR(jetEta.at(k), jetPhi.at(k), groomjteta[m], groomjtphi[m]) > s.groomedJetMatchingCut) continue;
            isMatched = true;
            subjet1Pt[j].push_back(subjetPt->at(m).at(0));
            subjet1Eta[j].push_back(subjetEta->at(m).at(0));
            subjet1Phi[j].push_back(subjetPhi->at(m).at(0));
            if(subjetEta->at(m).size()==2){
              subjet2Pt[j].push_back(subjetPt->at(m).at(1));
              subjet2Eta[j].push_back(subjetEta->at(m).at(1));
              subjet2Phi[j].push_back(subjetPhi->at(m).at(1));
              dR12[j].push_back(dR(subjetEta->at(m).at(0),subjetPhi->at(m).at(0),subjetEta->at(m).at(1),subjetPhi->at(m).at(1)));
            }else{
              subjet2Pt[j].push_back(-999);
              subjet2Eta[j].push_back(-999);
              subjet2Phi[j].push_back(-999);
              dR12[j].push_back(-999); 
            }
            break;      
          }
          if(!isMatched){//fix by lowering pt cut on forest
            std::cout << "Error: unmatched jet!" << std::endl;
            std::cout << jetPt.at(k) <<  " " << jetEta.at(k) << " " << jetPhi.at(k)<< " "  << std::endl;
            for(int m = 0; m<groomnref; m++){
              std::cout << "     " <<  groomjteta[m] << " " << groomjtphi[m] << std::endl;
            }
            subjet1Pt[j].push_back(-999);
            subjet1Eta[j].push_back(-999);
            subjet1Phi[j].push_back(-999);
            subjet2Pt[j].push_back(-999);
            subjet2Eta[j].push_back(-999);
            subjet2Phi[j].push_back(-999);
            dR12[j].push_back(-999); 
          }
        }
      }
      outTree->Fill();
   }//end event loop
    f->Close();
  }//end file loop
  output->Write();
  output->Close();  
}

//**************************************************************************
int main(int argc, const char* argv[])
{
  if(argc != 5)
  {
    std::cout << "Usage: <job> <totalJobs> <fileList> <dummyVar>" << std::endl;
    return 1;
  }  

  int job = std::atoi(argv[1]);
  int totalJobs = std::atoi(argv[2]);
  std::string fList = argv[3];
  int dummy = std::atoi(argv[4]);
  std::string buffer;
  std::vector<std::string> listOfFiles;
  std::ifstream inFile(fList.data());

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    int line = 0;
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      if(line%totalJobs==job) listOfFiles.push_back(buffer);
      line++;
    }
  }

  dummy = dummy+1;

  gammaJetSkim(listOfFiles,job);

  return 0; 
}
