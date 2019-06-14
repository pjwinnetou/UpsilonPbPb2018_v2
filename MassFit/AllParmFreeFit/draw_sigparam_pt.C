#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <string>
#include <math.h>

#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TMath.h>
#include <math.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>

#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>

#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooCategory.h>
#include <RooGenericPdf.h>
#include <RooFFTConvPdf.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooConstVar.h>

#include "/home/deathold/work/CMS/analysis/Upsilon_v2/upsilonV2/commonUtility.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/upsilonV2/cutsAndBinUpsilonV2.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/upsilonV2/Style_jaebeom.h"

using namespace std;
using namespace RooFit;

void draw_sigparam_pt(int fver=1)
{
  /////////////////////////////////////////////////////////
  //// set style
  /////////////////////////////////////////////////////////
  
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameLineColor(kBlack);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);

  gStyle->SetTextSize(0.04);
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.048,"xyz");

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.12) ; 
  
  /////////////////////////////////////////////////////////
  //// binning setting
  /////////////////////////////////////////////////////////
  
  int tmpBinNum;
  double tmpArr1[4] = {0.0, 3.0, 6.0, 50.0};
  double tmpArr2[8] = {0.0, 2.0, 4.0, 6.0, 8.0, 12.0, 20.0, 50.0};

  if(fver == 1) tmpBinNum = 3;
  if(fver == 2) tmpBinNum = 7;

  const int nBin = tmpBinNum;
  
  double binArr[nBin+1];

  cout << "nBin = " << nBin << endl;
  for (int ib =0; ib < nBin+1; ib ++ ) {
    if(fver==1) binArr[ib] = tmpArr1[ib]; 
    if(fver==2) binArr[ib] = tmpArr2[ib]; 
    cout << ib <<"th bin = " << binArr[ib] << endl;
  }

  /////////////////////////////////////////////////////////
  //// Open RooDataFile
  /////////////////////////////////////////////////////////
 
  //file and ws
  TFile *fileIn[nBin];
  RooWorkspace* ws[nBin];

  // parameters 
  double alpha1s[nBin];
  double alpha1sErr[nBin];
  double n1s[nBin];
  double n1sErr[nBin];
  double sigma1s_1[nBin];
  double sigma1s_1Err[nBin];
  double f1s[nBin];
  double f1sErr[nBin];
  double x1s[nBin];
  double x1sErr[nBin];
  
  for (int ib =0; ib < nBin; ib ++ ) {
    //// read files
    if(fver==1) fileIn[ib]= new TFile(Form("fitRes/Folder1/fitresults_upsilon_DoubleCB_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality0-200.root",binArr[ib],binArr[ib+1]),"read");
    if(fver==2) fileIn[ib]= new TFile(Form("fitRes/Folder2/fitresults_upsilon_DoubleCB_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality0-200.root",binArr[ib],binArr[ib+1]),"read");
    cout << ib << "th file = " << fileIn[ib]->GetName() << endl;
    if (fileIn[ib]->IsZombie()) { cout << "CANNOT open data root file\n"; }
    fileIn[ib]->cd();
    ws[ib]= (RooWorkspace*)fileIn[ib]->Get("workspace");
    ws[ib]->Print();

    // get parameters
    alpha1s[ib]=ws[ib]->var("alpha1s_1")->getVal();
    alpha1sErr[ib]=ws[ib]->var("alpha1s_1")->getError();
    n1s[ib]=ws[ib]->var("n1s_1")->getVal();
    n1sErr[ib]=ws[ib]->var("n1s_1")->getError();
    sigma1s_1[ib]=ws[ib]->var("sigma1s_1")->getVal();
    sigma1s_1Err[ib]=ws[ib]->var("sigma1s_1")->getError();
    f1s[ib]=ws[ib]->var("f1s")->getVal();
    f1sErr[ib]=ws[ib]->var("f1s")->getError();
    x1s[ib]=ws[ib]->var("x1s")->getVal();
    x1sErr[ib]=ws[ib]->var("x1s")->getError();
  }

  //Parameter limit
  TFile *lfile;
  if(fver==1) lfile = new TFile("fitRes/Folder1/parameter_limit_Folder1.root","read");
  if(fver==2) lfile = new TFile("fitRes/Folder2/parameter_limit_Folder2.root","read");

  TH1D* h_alphaUp = (TH1D*) lfile -> Get("hpt_alphaUp");
  TH1D* h_alphaLow = (TH1D*) lfile -> Get("hpt_alphaLow");
  TH1D* h_nUp = (TH1D*) lfile -> Get("hpt_nUp");
  TH1D* h_nLow = (TH1D*) lfile -> Get("hpt_nLow");
  TH1D* h_fUp = (TH1D*) lfile -> Get("hpt_fUp");
  TH1D* h_fLow = (TH1D*) lfile -> Get("hpt_fLow");
  TH1D* h_xUp = (TH1D*) lfile -> Get("hpt_xUp");
  TH1D* h_xLow = (TH1D*) lfile -> Get("hpt_xLow");
  TH1D* h_sigmaUp = (TH1D*) lfile -> Get("hpt_sigmaUp");
  TH1D* h_sigmaLow = (TH1D*) lfile -> Get("hpt_sigmaLow");

  //// histogram
  TH1D* h1_alpha1s = new TH1D("h1_alpha1s","h1_alpha1s;p_{T} (GeV/c);alpha1s",nBin,binArr); 
  TH1D* h1_n1s = new TH1D("h1_n1s","h1_n1s;p_{T} (GeV/c);n1s",nBin,binArr); 
  TH1D* h1_sigma1s_1 = new TH1D("h1_sigma1s_1","h1_sigma1s_1;p_{T} (GeV/c);sigma1s",nBin,binArr); 
  TH1D* h1_f1s = new TH1D("h1_f1s","h1_f1s;p_{T} (GeV/c);f1s",nBin,binArr); 
  TH1D* h1_x1s = new TH1D("h1_x1s","h1_x1s;p_{T} (GeV/c);x1s",nBin,binArr); 
  
  for (int ib =0; ib < nBin; ib ++ ) {
    h1_alpha1s->SetBinContent(ib+1,alpha1s[ib]);   
    h1_alpha1s->SetBinError(ib+1,alpha1sErr[ib]);   
    h1_n1s->SetBinContent(ib+1,n1s[ib]);   
    h1_n1s->SetBinError(ib+1,n1sErr[ib]);   
    h1_sigma1s_1->SetBinContent(ib+1,sigma1s_1[ib]);   
    h1_sigma1s_1->SetBinError(ib+1,sigma1s_1Err[ib]);   
    h1_f1s->SetBinContent(ib+1,f1s[ib]);   
    h1_f1s->SetBinError(ib+1,f1sErr[ib]);   
    h1_x1s->SetBinContent(ib+1,x1s[ib]);   
    h1_x1s->SetBinError(ib+1,x1sErr[ib]);   
  }

  //// actual draw
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.04);
  
  TCanvas* c_alpha1s = new TCanvas("c_alpha1s","c_alpha1s",600,600);
  c_alpha1s->cd();
  h1_alpha1s->GetXaxis()->CenterTitle(1);
  h1_alpha1s->GetYaxis()->CenterTitle(1);
  h1_alpha1s->GetYaxis()->SetRangeUser(0,5.5);
  h1_alpha1s->SetMarkerColor(kRed);
  h1_alpha1s->SetLineColor(kRed);
  h1_alpha1s->SetMarkerStyle(kFullCircle);
  h1_alpha1s->Draw("pe");
  h_alphaUp->SetLineColor(kBlack);
  h_alphaUp->SetLineStyle(kDashed);
  h_alphaUp->Draw("hist same");
  h_alphaLow->SetLineColor(kBlack);
  h_alphaLow->SetLineStyle(kDashed);
  h_alphaLow->Draw("hist same");
  c_alpha1s->SaveAs(Form("FitInfoPlot/pt_alpha1s_folder%d.pdf",fver));
  
  TCanvas* c_n1s = new TCanvas("c_n1s","c_n1s",600,600);
  c_n1s->cd();
  h1_n1s->GetXaxis()->CenterTitle(1);
  h1_n1s->GetYaxis()->CenterTitle(1);
  h1_n1s->GetYaxis()->SetRangeUser(0,10);
  h1_n1s->SetMarkerColor(kRed);
  h1_n1s->SetLineColor(kRed);
  h1_n1s->SetMarkerStyle(kFullCircle);
  h1_n1s->Draw("pe");
  h_nUp->SetLineColor(kBlack);
  h_nUp->SetLineStyle(kDashed);
  h_nUp->Draw("hist same");
  h_nLow->SetLineColor(kBlack);
  h_nLow->SetLineStyle(kDashed);
  h_nLow->Draw("hist same");
  c_n1s->SaveAs(Form("FitInfoPlot/pt_n1s_folder%d.pdf",fver));
 
  TCanvas* c_sigma1s_1 = new TCanvas("c_sigma1s_1","c_sigma1s_1",600,600);
  c_sigma1s_1->cd();
  h1_sigma1s_1->GetXaxis()->CenterTitle(1);
  h1_sigma1s_1->GetYaxis()->CenterTitle(1);
  h1_sigma1s_1->GetYaxis()->SetRangeUser(0,0.31);
  h1_sigma1s_1->SetMarkerColor(kRed);
  h1_sigma1s_1->SetLineColor(kRed);
  h1_sigma1s_1->SetMarkerStyle(kFullCircle);
  h1_sigma1s_1->Draw("pe");
  h_sigmaUp->SetLineColor(kBlack);
  h_sigmaUp->SetLineStyle(kDashed);
  h_sigmaUp->Draw("hist same");
  h_sigmaLow->SetLineColor(kBlack);
  h_sigmaLow->SetLineStyle(kDashed);
  h_sigmaLow->Draw("hist same");
  c_sigma1s_1->SaveAs(Form("FitInfoPlot/pt_sigma1s_folder%d.pdf",fver));
  
  TCanvas* c_f1s = new TCanvas("c_f1s","c_f1s",600,600);
  c_f1s->cd();
  h1_f1s->GetXaxis()->CenterTitle(1);
  h1_f1s->GetYaxis()->CenterTitle(1);
  h1_f1s->GetYaxis()->SetRangeUser(0,1.0);
  h1_f1s->SetMarkerColor(kRed);
  h1_f1s->SetLineColor(kRed);
  h1_f1s->SetMarkerStyle(kFullCircle);
  h1_f1s->Draw("pe");
  h_fUp->SetLineColor(kBlack);
  h_fUp->SetLineStyle(kDashed);
  h_fUp->Draw("hist same");
  h_fLow->SetLineColor(kBlack);
  h_fLow->SetLineStyle(kDashed);
  h_fLow->Draw("hist same");
//  latex->DrawLatex(0.45,0.83,Form("%s, binning for #Upsilon(%dS)",szAA.Data(),states));
  c_f1s->SaveAs(Form("FitInfoPlot/pt_f1s_folder%d.pdf",fver));
  
  TCanvas* c_x1s = new TCanvas("c_x1s","c_x1s",600,600);
  c_x1s->cd();
  h1_x1s->GetXaxis()->CenterTitle(1);
  h1_x1s->GetYaxis()->CenterTitle(1);
  h1_x1s->GetYaxis()->SetRangeUser(0,1);
  h1_x1s->SetMarkerColor(kRed);
  h1_x1s->SetLineColor(kRed);
  h1_x1s->SetMarkerStyle(kFullCircle);
  h1_x1s->Draw("pe");
  h_xUp->SetLineColor(kBlack);
  h_xUp->SetLineStyle(kDashed);
  h_xUp->Draw("hist same");
  h_xLow->SetLineColor(kBlack);
  h_xLow->SetLineStyle(kDashed);
  h_xLow->Draw("hist same");
  c_x1s->SaveAs(Form("FitInfoPlot/pt_x1s_folder%d.pdf",fver));

  TFile *wf = new TFile(Form("FitInfoPlot/ParmInfo_Pt_folder%d.root",fver),"recreate");
  wf->cd();
  h1_alpha1s->Write();
  h1_n1s->Write();
  h1_sigma1s_1->Write();
  h1_f1s->Write();
  h1_x1s->Write();

  double avg_alpha = 0;
  double avg_n = 0;
  double avg_f = 0;
  double avg_sigma = 0;
  double avg_x = 0;

  for(int ib = 1; ib<=nBin; ib++){
    avg_alpha += h1_alpha1s->GetBinContent(ib);
    avg_n     += h1_n1s->GetBinContent(ib);
    avg_sigma += h1_sigma1s_1->GetBinContent(ib);
    avg_f     += h1_f1s->GetBinContent(ib);
    avg_x     += h1_x1s->GetBinContent(ib);
  }
  
  avg_alpha = avg_alpha/h1_alpha1s->GetNbinsX();
  avg_n     = avg_n/h1_n1s->GetNbinsX();
  avg_sigma = avg_sigma/h1_sigma1s_1->GetNbinsX();
  avg_f     = avg_f/h1_f1s->GetNbinsX();
  avg_x     = avg_x/h1_x1s->GetNbinsX();
  
  double Err_alpha = 0;
  double Err_n = 0;
  double Err_sigma = 0;
  double Err_f = 0;
  double Err_x = 0;
  
  for(int ib = 1; ib<=nBin; ib++){
    Err_alpha += (avg_alpha - h1_alpha1s   -> GetBinContent(ib))*(avg_alpha - h1_alpha1s   -> GetBinContent(ib));
    Err_n     += (avg_n     - h1_n1s       -> GetBinContent(ib))*(avg_n     - h1_n1s       -> GetBinContent(ib));
    Err_sigma += (avg_sigma - h1_sigma1s_1 -> GetBinContent(ib))*(avg_sigma - h1_sigma1s_1 -> GetBinContent(ib));
    Err_f     += (avg_f     - h1_f1s       -> GetBinContent(ib))*(avg_f     - h1_f1s       -> GetBinContent(ib));
    Err_x     += (avg_x     - h1_x1s       -> GetBinContent(ib))*(avg_x     - h1_x1s       -> GetBinContent(ib));
  }

  Err_alpha = TMath::Sqrt( Err_alpha / (nBin*(nBin-1)) );
  Err_n = TMath::Sqrt( Err_n / (nBin*(nBin-1)) );
  Err_sigma = TMath::Sqrt( Err_sigma / (nBin*(nBin-1)) );
  Err_f = TMath::Sqrt( Err_f / (nBin*(nBin-1)) );
  Err_x = TMath::Sqrt( Err_x / (nBin*(nBin-1)) );
  
  TH1D* hParm = new TH1D("hParm",";;;",5,0,5);
  hParm->SetBinContent(1, avg_alpha);
  hParm->SetBinContent(2, avg_n);
  hParm->SetBinContent(3, avg_sigma);
  hParm->SetBinContent(4, avg_f);
  hParm->SetBinContent(5, avg_x);
  hParm->SetBinError(1, Err_alpha);
  hParm->SetBinError(2, Err_n);
  hParm->SetBinError(3, Err_sigma);
  hParm->SetBinError(4, Err_f);
  hParm->SetBinError(5, Err_x);

  hParm->Write();  
}
