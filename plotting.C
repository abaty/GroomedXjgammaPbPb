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

  const int nCentBins = 9;
  const float centBinsLow[nCentBins] = {0,0,0,10,10,20,30,30,50};
  const float centBinsHigh[nCentBins] = {10,20,30,30,50,50,50,100,100};
  const int nPhoPtBins = 9;
  const float phoBinLow[nPhoPtBins] = {40,60,80,100,40,60,80,60,40};
  const float phoBinLow[nPhoPtBins] = {60,80,100,999,80,100,999,999,999};

  //plots
  TFile * out = new TFile::Open("plots.root","recreate");
  TH1D * xjg_DR12_LT_Z0B0p1[3][nCentBins][nPho]; 
  TH1D * xjg_DR12_GT_Z0B0p1[3][nCentBins][nPho]; 
  TH1D * xjg_DR12_LT_Z0p5Z1p5[3][nCentBins][nPho]; 
  TH1D * xjg_DR12_GT_Z0p5B1p5[3][nCentBins][npho]; 
  for(int p = 0; p<nPho; p++){
    for(int c = 0; c<nCentBins; c++){
      for(int i = 0; i<3; i++){
        std::cout << p << " " << c << " " << i << std::endl;
        xjg_DR12_LT_Z0B0p1[i][c][p] = new TH1D(Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),20,0,2);
        xjg_DR12_GT_Z0B0p1[i][c][p] = new TH1D(Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),20,0,2);
        xjg_DR12_LT_Z0p1B1p5[i][c][p] = new TH1D(Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),20,0,2);
        xjg_DR12_GT_Z0p1B1p5[i][c][p] = new TH1D(Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),20,0,2);
      
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%d && phoPt<%d && dR12_Z0B0p1<%f",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p],i/(float)10));
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%d && phoPt<%d && dR12_Z0B0p1>%f",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p],i/(float)10));
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_LT_Z0p5Z1p5_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%d && phoPt<%d && dR12_Z0p5Z1p5<%f",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p],i/(float)10));
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_GT_Z0p5Z1p5_R0p%d_%d_%d_%d_%d",i,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%d && phoPt<%d && dR12_Z0p5Z1p5>%f",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p],i/(float)10));
      }
    }
  }
  
  



}
