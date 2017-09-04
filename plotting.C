#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TColor.h"
#include "TAttAxis.h"
#include "TAttLine.h"
#include "TLegend.h"
#include "TAttMarker.h"
#include <vector>
#include <iostream>
#include <string>

void sidebandSubtr(TH1D * signalClone, TH1D * sideband, float purity){
  TH1D * temp = new TH1D();
  temp = (TH1D*)sideband->Clone("temp");
  temp->Scale(1-purity);
  signalClone->Add(temp,-1);
  signalClone->Scale(1./purity);
  delete temp;
}

void plotting(){
  TH1::SetDefaultSumw2();

  //setup
  TFile * in = TFile::Open("skim_sept1.root","read");
  TTree * t = (TTree*)in->Get("skim");

  /*
  int hiBin;
  std::vector<float> * jetPt = 0;
  std::vector<float> * phoPt = 0;
  std::vector<float> * sigIEIE = 0;
  std::vector<float> * DR12_Z0B0p1 = 0;
  std::vector<float> * DR12_Z0p5B1p5 = 0;

  t->SetBranchAddress("jetPt",&jetPt);
  t->SetBranchAddress("phoPt",&phoPt);
  t->SetBranchAddress("hiBin",&hiBin);
  t->SetBranchAddress("phoSigIEIE",&sigIEIE);
  t->SetBranchAddress("dR12_Z0B0p1",&DR12_Z0B0p1);
  t->SetBranchAddress("dR12_Z0p5B1p5",&DR12_Z0p5B1p5);
  */

  TCanvas * c = new TCanvas("c","c",800,600);

  const int nCentBins = 7;
  const int centBinsLow[nCentBins] = {0,0,30,0,10,30,50};
  const int centBinsHigh[nCentBins] = {100,30,100,10,30,50,100};
  const int nPhoPtBins = 8;
  const int phoBinLow[nPhoPtBins] = {40,60,40,50,60,80,80,100};
  const int phoBinHigh[nPhoPtBins] = {999,999,50,60,80,999,100,999};

  //magic numbers from Table 10 of CMS AN-2016/054 v12 (gamma-jet xjg analysis)
  const float purity[nPhoPtBins][nCentBins] = {{0.704402, 0.695052, 0.745314, 0.670393, 0.712452, 0.747344, 0.737428},{0.725267, 0.708896, 0.785853, 0.681689, 0.730339, 0.775276, 0.825922},{0.692733, 0.684606, 0.730031, 0.654991, 0.700615, 0.733907, 0.714753},{0.707091,0.699483, 0.748273, 0.672784, 0.722749, 0.754953, 0.726901},{0.719437, 0.702283, 0.785353, 0.676977, 0.721539, 0.77048, 0.835029},{0.737224, 0.721948, 0.789327, 0.685621, 0.751769, 0.787902, 0.816627},{0.722065, 0.702349, 0.787935, 0.67868, 0.725783, 0.771699, 0.850336},{0.758906, 0.749149, 0.799366, 0.700252, 0.792531, 0.821031, 0.797929}};


  //plots
  TFile * out = TFile::Open("outputPlots.root","recreate");
  TH1D * nPhotons_signal[nCentBins][nPhoPtBins]; 
  TH1D * nPhotons_sideband[nCentBins][nPhoPtBins]; 
  TH1D * xjg_DR12_LT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH1D * xjg_DR12_GT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH1D * xjg_DR12_LT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * xjg_DR12_GT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * side_xjg_DR12_LT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH1D * side_xjg_DR12_GT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH1D * side_xjg_DR12_LT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * side_xjg_DR12_GT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * subtr_xjg_DR12_LT_Z0B0p1[3*nCentBins*nPhoPtBins]; 
  TH1D * subtr_xjg_DR12_GT_Z0B0p1[3*nCentBins*nPhoPtBins]; 
  TH1D * subtr_xjg_DR12_LT_Z0p5B1p5[3*nCentBins*nPhoPtBins]; 
  TH1D * subtr_xjg_DR12_GT_Z0p5B1p5[3*nCentBins*nPhoPtBins];
  for(int p = 0; p<nPhoPtBins; p++){//photon loop
    for(int c = 0; c<nCentBins; c++){//cent loop
      nPhotons_signal[c][p] = new TH1D(Form("nPhotons_signal_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("nPhotons_signal_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),1,0,2);
      t->Draw(Form("1>>nPhotons_signal_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("phoSigIEIE<0.01 && hiBin>=%d && hiBin<%d && phoPt>%f && phoPt<%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));
      nPhotons_sideband[c][p] = new TH1D(Form("nPhotons_sideband_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("nPhotons_sideband_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),1,0,2);
      t->Draw(Form("1>>nPhotons_sideband_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("phoSigIEIE>0.011 && phoSigIEIE<0.017 && hiBin>=%d && hiBin<%d && phoPt>%f && phoPt<%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));

      for(int i = 0; i<3; i++){//dR loop
        std::cout << p << " " << c << " " << i << std::endl;
        xjg_DR12_LT_Z0B0p1[i][c][p] = new TH1D(Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),16,0,2);
        xjg_DR12_GT_Z0B0p1[i][c][p] = new TH1D(Form("xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),16,0,2);
        xjg_DR12_LT_Z0p5B1p5[i][c][p] = new TH1D(Form("xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),16,0,2);
        xjg_DR12_GT_Z0p5B1p5[i][c][p] = new TH1D(Form("xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),16,0,2);
      
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && dR12_Z0B0p1<%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && dR12_Z0B0p1>%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && dR12_Z0p5B1p5<%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && dR12_Z0p5B1p5>%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        
        //sideband
        side_xjg_DR12_LT_Z0B0p1[i][c][p] = new TH1D(Form("side_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),16,0,2);
        side_xjg_DR12_GT_Z0B0p1[i][c][p] = new TH1D(Form("side_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),16,0,2);
        side_xjg_DR12_LT_Z0p5B1p5[i][c][p] = new TH1D(Form("side_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),16,0,2);
        side_xjg_DR12_GT_Z0p5B1p5[i][c][p] = new TH1D(Form("side_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),16,0,2);
      
        t->Draw(Form("jetPt/phoPt>>side_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && dR12_Z0B0p1<%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>side_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && dR12_Z0B0p1>%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>side_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && dR12_Z0p5B1p5<%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>side_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017  && phoPt>%f && phoPt<%f && dR12_Z0p5B1p5>%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        
        //normalizing by number of photons
        xjg_DR12_LT_Z0B0p1[i][c][p]->Scale(1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        xjg_DR12_GT_Z0B0p1[i][c][p]->Scale(1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        xjg_DR12_LT_Z0p5B1p5[i][c][p]->Scale(1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        xjg_DR12_GT_Z0p5B1p5[i][c][p]->Scale(1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        side_xjg_DR12_LT_Z0B0p1[i][c][p]->Scale(1./(float)nPhotons_sideband[c][p]->GetBinContent(1)); 
        side_xjg_DR12_GT_Z0B0p1[i][c][p]->Scale(1./(float)nPhotons_sideband[c][p]->GetBinContent(1));          
        side_xjg_DR12_LT_Z0p5B1p5[i][c][p]->Scale(1./(float)nPhotons_sideband[c][p]->GetBinContent(1)); 
        side_xjg_DR12_GT_Z0p5B1p5[i][c][p]->Scale(1./(float)nPhotons_sideband[c][p]->GetBinContent(1));

        //weird indexing is due to a crash
        subtr_xjg_DR12_LT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p] = (TH1D*) xjg_DR12_LT_Z0B0p1[i][c][p]->Clone(Form("subtr_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        sidebandSubtr(subtr_xjg_DR12_LT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p], side_xjg_DR12_LT_Z0B0p1[i][c][p], purity[p][c]);
        subtr_xjg_DR12_GT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p] = (TH1D*) xjg_DR12_GT_Z0B0p1[i][c][p]->Clone(Form("subtr_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        sidebandSubtr(subtr_xjg_DR12_GT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p], side_xjg_DR12_GT_Z0B0p1[i][c][p], purity[p][c]);
        subtr_xjg_DR12_LT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p] = (TH1D*) xjg_DR12_LT_Z0p5B1p5[i][c][p]->Clone(Form("subtr_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        sidebandSubtr(subtr_xjg_DR12_LT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p], side_xjg_DR12_LT_Z0B0p1[i][c][p], purity[p][c]);
        subtr_xjg_DR12_GT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p] = (TH1D*) xjg_DR12_GT_Z0p5B1p5[i][c][p]->Clone(Form("subtr_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        sidebandSubtr(subtr_xjg_DR12_GT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p], side_xjg_DR12_GT_Z0B0p1[i][c][p], purity[p][c]);
      }
    }
  }
  std::cout << "writing!" << std::endl;
  out->Write(); 
  std::cout << "done!" << std::endl;
}
