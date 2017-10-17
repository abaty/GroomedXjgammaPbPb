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

void formatPlot(TH1D * h, int color = 0, int lineStyle = 1, int markerStyle = 8){
  h->SetTitle("");
  h->GetXaxis()->SetTitle("x_{j#gamma}");
  h->GetYaxis()->SetTitle("#frac{1}{N_{#gamma}} #frac{dN_{j#gamma}}{dx_{j#gamma}}");
  h->SetMarkerStyle(markerStyle);
  h->GetYaxis()->SetRangeUser(-0.2,1.8);
  if(color==0){
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
  }
  if(color==1){
    h->SetMarkerColor(kBlue);
    h->SetLineColor(kBlue);
  }
  if(color==2){
    h->SetMarkerColor(kRed);
    h->SetLineColor(kRed);
  }
  h->SetLineStyle(lineStyle);
  h->SetStats(0);
}

void plotting(){
  TH1::SetDefaultSumw2();

  //setup
  TFile * in = TFile::Open("skim_sept1.root","read");
  TTree * t = (TTree*)in->Get("skim");

  const int nXJGBins = 16;
  const float maxXJG = 2;

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
  TH2D * xjg_Z0B0p1[nCentBins][nPhoPtBins]; 
  TH1D * xjg_DR12_LT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * xjg_DR12_GT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH2D * xjg_Z0p5B1p5[nCentBins][nPhoPtBins]; 
  TH1D * mixxjg_DR12_LT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH1D * mixxjg_DR12_GT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH2D * mixxjg_Z0B0p1[nCentBins][nPhoPtBins]; 
  TH1D * mixxjg_DR12_LT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * mixxjg_DR12_GT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH2D * mixxjg_Z0p5B1p5[nCentBins][nPhoPtBins]; 
  TH1D * mixsubxjg_DR12_LT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH1D * mixsubxjg_DR12_GT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH2D * mixsubxjg_Z0B0p1[nCentBins][nPhoPtBins]; 
  TH1D * mixsubxjg_DR12_LT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * mixsubxjg_DR12_GT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH2D * mixsubxjg_Z0p5B1p5[nCentBins][nPhoPtBins]; 
  TH1D * side_xjg_DR12_LT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH1D * side_xjg_DR12_GT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH2D * side_xjg_Z0B0p1[nCentBins][nPhoPtBins]; 
  TH1D * side_xjg_DR12_LT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * side_xjg_DR12_GT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH2D * side_xjg_Z0p5B1p5[nCentBins][nPhoPtBins]; 
  TH1D * mixside_xjg_DR12_LT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH1D * mixside_xjg_DR12_GT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH2D * mixside_xjg_Z0B0p1[nCentBins][nPhoPtBins]; 
  TH1D * mixside_xjg_DR12_LT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * mixside_xjg_DR12_GT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH2D * mixside_xjg_Z0p5B1p5[nCentBins][nPhoPtBins]; 
  TH1D * mixsubside_xjg_DR12_LT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH1D * mixsubside_xjg_DR12_GT_Z0B0p1[3][nCentBins][nPhoPtBins]; 
  TH2D * mixsubside_xjg_Z0B0p1[nCentBins][nPhoPtBins]; 
  TH1D * mixsubside_xjg_DR12_LT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH1D * mixsubside_xjg_DR12_GT_Z0p5B1p5[3][nCentBins][nPhoPtBins]; 
  TH2D * mixsubside_xjg_Z0p5B1p5[nCentBins][nPhoPtBins]; 
  TH1D * subtr_xjg_DR12_LT_Z0B0p1[3*nCentBins*nPhoPtBins]; 
  TH1D * subtr_xjg_DR12_GT_Z0B0p1[3*nCentBins*nPhoPtBins]; 
  TH2D * subtr_xjg_Z0B0p1[nCentBins][nPhoPtBins]; 
  TH1D * subtr_xjg_DR12_LT_Z0p5B1p5[3*nCentBins*nPhoPtBins]; 
  TH1D * subtr_xjg_DR12_GT_Z0p5B1p5[3*nCentBins*nPhoPtBins];
  TH2D * subtr_xjg_Z0p5B1p5[nCentBins][nPhoPtBins];

  for(int p = 0; p<nPhoPtBins; p++){//photon loop
    for(int c = 0; c<nCentBins; c++){//cent loop
      nPhotons_signal[c][p] = new TH1D(Form("nPhotons_signal_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("nPhotons_signal_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),1,0,2);
      t->Draw(Form("1>>nPhotons_signal_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("phoSigIEIE<0.01 && hiBin>=%d && hiBin<%d && phoPt>%f && phoPt<%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));
      nPhotons_sideband[c][p] = new TH1D(Form("nPhotons_sideband_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("nPhotons_sideband_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),1,0,2);
      t->Draw(Form("1>>nPhotons_sideband_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("phoSigIEIE>0.011 && phoSigIEIE<0.017 && hiBin>=%d && hiBin<%d && phoPt>%f && phoPt<%f",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));

      //2D correlations
      xjg_Z0B0p1[c][p] = new TH2D(Form("xjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG,16,0,0.8); 
      mixxjg_Z0B0p1[c][p] = new TH2D(Form("mixxjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixxjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG,16,0,0.8); 
      side_xjg_Z0B0p1[c][p] = new TH2D(Form("side_xjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG,16,0,0.8); 
      mixside_xjg_Z0B0p1[c][p] = new TH2D(Form("mixside_xjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixside_xjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG,16,0,0.8); 
      t->Draw(Form("dR12_Z0B0p1:jetPt/phoPt>>xjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));
      t->Draw(Form("mixdR12_Z0B0p1:mixJetPt/phoPt>>mixxjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f  && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));
      t->Draw(Form("dR12_Z0B0p1:jetPt/phoPt>>side_xjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));
      t->Draw(Form("mixdR12_Z0B0p1:mixJetPt/phoPt>>mixside_xjg_Z0B0p1_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));

      xjg_Z0p5B1p5[c][p] = new TH2D(Form("xjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG,16,0,0.8); 
      mixxjg_Z0p5B1p5[c][p] = new TH2D(Form("mixxjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixxjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG,16,0,0.8); 
      side_xjg_Z0p5B1p5[c][p] = new TH2D(Form("side_xjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG,16,0,0.8); 
      mixside_xjg_Z0p5B1p5[c][p] = new TH2D(Form("mixside_xjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixside_xjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG,16,0,0.8); 
      t->Draw(Form("dR12_Z0p5B1p5:jetPt/phoPt>>xjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));
      t->Draw(Form("mixdR12_Z0p5B1p5:mixJetPt/phoPt>>mixxjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f  && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));
      t->Draw(Form("dR12_Z0p5B1p5:jetPt/phoPt>>side_xjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));
      t->Draw(Form("mixdR12_Z0p5B1p5:mixJetPt/phoPt>>mixside_xjg_Z0p5B1p5_%d_%d_%d_%d",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p]));


      for(int i = 0; i<3; i++){//dR loop
        std::cout << p << " " << c << " " << i << std::endl;
        xjg_DR12_LT_Z0B0p1[i][c][p] = new TH1D(Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        xjg_DR12_GT_Z0B0p1[i][c][p] = new TH1D(Form("xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        xjg_DR12_LT_Z0p5B1p5[i][c][p] = new TH1D(Form("xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        xjg_DR12_GT_Z0p5B1p5[i][c][p] = new TH1D(Form("xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        mixxjg_DR12_LT_Z0B0p1[i][c][p] = new TH1D(Form("mixxjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixxjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        mixxjg_DR12_GT_Z0B0p1[i][c][p] = new TH1D(Form("mixxjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixxjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        mixxjg_DR12_LT_Z0p5B1p5[i][c][p] = new TH1D(Form("mixxjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixxjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        mixxjg_DR12_GT_Z0p5B1p5[i][c][p] = new TH1D(Form("mixxjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixxjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
      
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && dR12_Z0B0p1<%f && dR12_Z0B0p1>0 && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && dR12_Z0B0p1>%f && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && dR12_Z0p5B1p5<%f && dR12_Z0p5B1p5>0 && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && dR12_Z0p5B1p5>%f && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        
        t->Draw(Form("mixJetPt/phoPt>>mixxjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && mixdR12_Z0B0p1<%f && mixdR12_Z0B0p1>0 && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("mixJetPt/phoPt>>mixxjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && mixdR12_Z0B0p1>%f && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("mixJetPt/phoPt>>mixxjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && mixdR12_Z0p5B1p5<%f && mixdR12_Z0p5B1p5>0 && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("mixJetPt/phoPt>>mixxjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE<0.01 && phoPt>%f && phoPt<%f && mixdR12_Z0p5B1p5>%f && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
      
        mixsubxjg_DR12_LT_Z0B0p1[i][c][p] = (TH1D*)  xjg_DR12_LT_Z0B0p1[i][c][p]->Clone(Form("mixsubxjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        mixsubxjg_DR12_GT_Z0B0p1[i][c][p] = (TH1D*)  xjg_DR12_GT_Z0B0p1[i][c][p]->Clone(Form("mixsubxjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        mixsubxjg_DR12_LT_Z0p5B1p5[i][c][p] = (TH1D*)xjg_DR12_LT_Z0p5B1p5[i][c][p]->Clone(Form("mixsubxjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        mixsubxjg_DR12_GT_Z0p5B1p5[i][c][p] = (TH1D*)xjg_DR12_GT_Z0p5B1p5[i][c][p]->Clone(Form("mixsubxjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        mixsubxjg_DR12_LT_Z0B0p1[i][c][p]->Add(mixxjg_DR12_LT_Z0B0p1[i][c][p],-1);
        mixsubxjg_DR12_GT_Z0B0p1[i][c][p]->Add(mixxjg_DR12_GT_Z0B0p1[i][c][p],-1); 
        mixsubxjg_DR12_LT_Z0p5B1p5[i][c][p]->Add(mixxjg_DR12_LT_Z0p5B1p5[i][c][p],-1);
        mixsubxjg_DR12_GT_Z0p5B1p5[i][c][p]->Add(mixxjg_DR12_GT_Z0p5B1p5[i][c][p],-1);
         
 
        //sideband
        side_xjg_DR12_LT_Z0B0p1[i][c][p] = new TH1D(Form("side_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        side_xjg_DR12_GT_Z0B0p1[i][c][p] = new TH1D(Form("side_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        side_xjg_DR12_LT_Z0p5B1p5[i][c][p] = new TH1D(Form("side_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        side_xjg_DR12_GT_Z0p5B1p5[i][c][p] = new TH1D(Form("side_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("side_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        mixside_xjg_DR12_LT_Z0B0p1[i][c][p] = new TH1D(Form("mixside_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixside_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        mixside_xjg_DR12_GT_Z0B0p1[i][c][p] = new TH1D(Form("mixside_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixside_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        mixside_xjg_DR12_LT_Z0p5B1p5[i][c][p] = new TH1D(Form("mixside_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixside_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
        mixside_xjg_DR12_GT_Z0p5B1p5[i][c][p] = new TH1D(Form("mixside_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixside_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),nXJGBins,0,maxXJG);
      
        t->Draw(Form("jetPt/phoPt>>side_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && dR12_Z0B0p1<%f && dR12_Z0B0p1>0 && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>side_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && dR12_Z0B0p1>%f && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>side_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && dR12_Z0p5B1p5<%f && dR12_Z0p5B1p5>0 && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("jetPt/phoPt>>side_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017  && phoPt>%f && phoPt<%f && dR12_Z0p5B1p5>%f && jetPt>30",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("mixJetPt/phoPt>>mixside_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && mixdR12_Z0B0p1<%f && mixdR12_Z0B0p1>0 && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("mixJetPt/phoPt>>mixside_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && mixdR12_Z0B0p1>%f && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("mixJetPt/phoPt>>mixside_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017 && phoPt>%f && phoPt<%f && mixdR12_Z0p5B1p5<%f && mixdR12_Z0p5B1p5>0 && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        t->Draw(Form("mixJetPt/phoPt>>mixside_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]),Form("mixJetWeight*(hiBin>=%d && hiBin<%d && phoSigIEIE>0.011 && phoSigIEIE<0.017  && phoPt>%f && phoPt<%f && mixdR12_Z0p5B1p5>%f && mixJetPt>30)",centBinsLow[c],centBinsHigh[c],(float)phoBinLow[p],(float)phoBinHigh[p],i/(float)10+0.1));
        
        mixsubside_xjg_DR12_LT_Z0B0p1[i][c][p] = (TH1D*)  side_xjg_DR12_LT_Z0B0p1[i][c][p]->Clone(Form("mixsubside_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        mixsubside_xjg_DR12_GT_Z0B0p1[i][c][p] = (TH1D*)  side_xjg_DR12_GT_Z0B0p1[i][c][p]->Clone(Form("mixsubside_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        mixsubside_xjg_DR12_LT_Z0p5B1p5[i][c][p] = (TH1D*)side_xjg_DR12_LT_Z0p5B1p5[i][c][p]->Clone(Form("mixsubside_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        mixsubside_xjg_DR12_GT_Z0p5B1p5[i][c][p] = (TH1D*)side_xjg_DR12_GT_Z0p5B1p5[i][c][p]->Clone(Form("mixsubside_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        mixsubside_xjg_DR12_LT_Z0B0p1[i][c][p]->Add(mixside_xjg_DR12_LT_Z0B0p1[i][c][p],-1);
        mixsubside_xjg_DR12_GT_Z0B0p1[i][c][p]->Add(mixside_xjg_DR12_GT_Z0B0p1[i][c][p],-1); 
        mixsubside_xjg_DR12_LT_Z0p5B1p5[i][c][p]->Add(mixside_xjg_DR12_LT_Z0p5B1p5[i][c][p],-1);
        mixsubside_xjg_DR12_GT_Z0p5B1p5[i][c][p]->Add(mixside_xjg_DR12_GT_Z0p5B1p5[i][c][p],-1);
        
        //normalizing by number of photons
        xjg_DR12_LT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        xjg_DR12_GT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        xjg_DR12_LT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        xjg_DR12_GT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        mixxjg_DR12_LT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        mixxjg_DR12_GT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        mixxjg_DR12_LT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        mixxjg_DR12_GT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        mixsubxjg_DR12_LT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        mixsubxjg_DR12_GT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        mixsubxjg_DR12_LT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        mixsubxjg_DR12_GT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_signal[c][p]->GetBinContent(1));
        side_xjg_DR12_LT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1)); 
        side_xjg_DR12_GT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1));          
        side_xjg_DR12_LT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1)); 
        side_xjg_DR12_GT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1));
        mixside_xjg_DR12_LT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1)); 
        mixside_xjg_DR12_GT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1));          
        mixside_xjg_DR12_LT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1)); 
        mixside_xjg_DR12_GT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1));
        mixsubside_xjg_DR12_LT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1)); 
        mixsubside_xjg_DR12_GT_Z0B0p1[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1));          
        mixsubside_xjg_DR12_LT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1)); 
        mixsubside_xjg_DR12_GT_Z0p5B1p5[i][c][p]->Scale((float)nXJGBins/(float)maxXJG*1./(float)nPhotons_sideband[c][p]->GetBinContent(1));

        //weird indexing is due to a crash
        subtr_xjg_DR12_LT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p] = (TH1D*) mixsubxjg_DR12_LT_Z0B0p1[i][c][p]->Clone(Form("subtr_xjg_DR12_LT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        sidebandSubtr(subtr_xjg_DR12_LT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p], mixsubside_xjg_DR12_LT_Z0B0p1[i][c][p], purity[p][c]);
        subtr_xjg_DR12_GT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p] = (TH1D*) mixsubxjg_DR12_GT_Z0B0p1[i][c][p]->Clone(Form("subtr_xjg_DR12_GT_Z0B0p1_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        sidebandSubtr(subtr_xjg_DR12_GT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p], mixsubside_xjg_DR12_GT_Z0B0p1[i][c][p], purity[p][c]);
        subtr_xjg_DR12_LT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p] = (TH1D*) mixsubxjg_DR12_LT_Z0p5B1p5[i][c][p]->Clone(Form("subtr_xjg_DR12_LT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        sidebandSubtr(subtr_xjg_DR12_LT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p], mixsubside_xjg_DR12_LT_Z0B0p1[i][c][p], purity[p][c]);
        subtr_xjg_DR12_GT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p] = (TH1D*) mixsubxjg_DR12_GT_Z0p5B1p5[i][c][p]->Clone(Form("subtr_xjg_DR12_GT_Z0p5B1p5_R0p%d_%d_%d_%d_%d",i+1,centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p]));
        sidebandSubtr(subtr_xjg_DR12_GT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p], mixsubside_xjg_DR12_GT_Z0B0p1[i][c][p], purity[p][c]);
        
        formatPlot(subtr_xjg_DR12_LT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p],2,1,20);
        formatPlot(subtr_xjg_DR12_LT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p],1,1,20);
        formatPlot(subtr_xjg_DR12_GT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p],2,1,25);
        formatPlot(subtr_xjg_DR12_GT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p],1,1,25);
      }
    }
  }
  std::cout << "writing!" << std::endl;
  out->Write(); 
  std::cout << "done!" << std::endl;
  std::cout << "Proceeding to Plotting..." << std::endl;

  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  for(int i = 0; i<3; i++){
    for(int j = i; j<3; j++){
      for(int p = 0; p<nPhoPtBins; p++){//photon loop
        for(int c = 0; c<nCentBins; c++){//cent loop
          subtr_xjg_DR12_LT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p]->Draw("p");
          subtr_xjg_DR12_LT_Z0p5B1p5[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p]->Draw("p same");
          subtr_xjg_DR12_GT_Z0B0p1[j*nCentBins*nPhoPtBins+c*nPhoPtBins+p]->Draw("p same");
          subtr_xjg_DR12_GT_Z0p5B1p5[j*nCentBins*nPhoPtBins+c*nPhoPtBins+p]->Draw("p same");  
  
          TLegend * l = new TLegend(0.6,0.5,0.9,0.9);
          l->AddEntry((TObject*)0,Form("%d%%-%d%%",centBinsLow[c],centBinsHigh[c]),"");
          l->AddEntry((TObject*)0,Form("%d<p_{T}^{#gamma}<%d",phoBinLow[p],phoBinHigh[p]),"");
          l->AddEntry(subtr_xjg_DR12_LT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p],Form("z=0, #beta=0.1, #DeltaR<0.%d",i+1),"p");
          l->AddEntry(subtr_xjg_DR12_GT_Z0B0p1[i*nCentBins*nPhoPtBins+c*nPhoPtBins+p],Form("z=0, #beta=0.1, #DeltaR>0.%d",j+1),"p");
          l->AddEntry(subtr_xjg_DR12_LT_Z0p5B1p5[j*nCentBins*nPhoPtBins+c*nPhoPtBins+p],Form("z=0.5, #beta=1.5, #DeltaR<0.%d",i+1),"p");
          l->AddEntry(subtr_xjg_DR12_GT_Z0p5B1p5[j*nCentBins*nPhoPtBins+c*nPhoPtBins+p],Form("z=0.5, #beta=1.5, #DeltaR>0.%d",j+1),"p");
          
          l->Draw("same");
          c1->SaveAs(Form("img/%d_%d_%d_%d_LT%d_GT%d.png",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p],i+1,j+1));
          c1->SaveAs(Form("img/%d_%d_%d_%d_LT%d_GT%d.pdf",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p],i+1,j+1));
          c1->SaveAs(Form("img/%d_%d_%d_%d_LT%d_GT%d.C",centBinsLow[c],centBinsHigh[c],phoBinLow[p],phoBinHigh[p],i+1,j+1));
          delete l;
        }
      }
    }
  }
}
