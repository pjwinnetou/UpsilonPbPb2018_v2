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

void draw_SigPDF_cent(int fver =5)
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
  double tmpArr4[5] = {0, 10, 30, 50, 100};
  double tmpArr5[10] = {0, 5, 10, 20, 30, 40, 50, 60, 80, 100};

  if(fver == 4) tmpBinNum = 4;
  if(fver == 5) tmpBinNum = 9;

  const int nBin = tmpBinNum;

  
  //String Set
  const char *str_name4[4] = {"010", "1030", "3050", "50100"};
  const char *str_name_leg4[4] = {"0-10", "10-30", "30-50", "50-100"};
  int str_color4[4] = {kBlack, kRed, kViolet, kBlue};
  
  const char *str_name5[9] = {"05", "510", "1020" ,"2030", "3040", "4050", "5060", "6080", "80100"};
  const char *str_name_leg5[9] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-80", "80-100"};
  int str_color5[9] = {kBlack, kGray, kGreen+2, kSpring, kYellow+1, kOrange+7, kRed, kViolet, kBlue};
  
  const char *str_name[nBin];
  const char *str_name_leg[nBin];
  int str_color[nBin]; 
  double binArr[nBin+1];
  for(int ib=0; ib<nBin+1;ib++){
    if(fver==4){ 
      binArr[ib] = tmpArr4[ib]; 
      str_name[ib] = str_name4[ib];
      str_name_leg[ib] = str_name_leg4[ib];
      str_color[ib] = str_color4[ib];
    }
    if(fver==5){
      binArr[ib] = tmpArr5[ib]; 
      str_name[ib] = str_name5[ib];
      str_name_leg[ib] = str_name_leg5[ib];
      str_color[ib] = str_color5[ib];
    }
  }


  /////////////////////////////////////////////////////////
  //// Open RooDataFile
  /////////////////////////////////////////////////////////
 
  //file and ws
  TFile *fileIn[nBin];
  RooWorkspace* ws[nBin];

  //PDF
  RooAbsPdf *cb1s[nBin];

  TString perc="%";

  for (int ib =0; ib < nBin; ib ++ ) {
    //// read files
    fileIn[ib]= new TFile(Form("fitRes/Folder%d/fitresults_upsilon_DoubleCB_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality%.f-%.f.root",fver,binArr[ib]*2,binArr[ib+1]*2),"read");
    cout << ib << "th file = " << fileIn[ib]->GetName() << endl;
    if (fileIn[ib]->IsZombie()) { cout << "CANNOT open data root file\n"; }
    fileIn[ib]->cd();
    ws[ib]= (RooWorkspace*)fileIn[ib]->Get("workspace");
    ws[ib]->Print();

    //PDF
    cb1s[ib] = ws[ib]->pdf("cb1s");

  }
 
  RooRealVar *mass = new RooRealVar("mass","dimuon mass", 8,14);
  const int nMassBin = 60;

  //// RooPlot
  RooPlot* myPlot = mass->frame(nMassBin);

  for (int ib =0; ib < nBin; ib ++ ) {
    cb1s[ib]->plotOn(myPlot,Name(Form("Cent%s",str_name[ib])),LineColor(str_color[ib]),LineStyle(1),LineWidth(2));
  }
/*
  cb1s[0]->plotOn(myPlot,Name("Cent010"),LineColor(kBlack),LineStyle(1),LineWidth(2));
  cb1s[1]->plotOn(myPlot,Name("Cent1030"),LineColor(kBlue),LineStyle(1),LineWidth(2));
  cb1s[2]->plotOn(myPlot,Name("Cent3050"),LineColor(kRed),LineStyle(1),LineWidth(2));
  cb1s[3]->plotOn(myPlot,Name("Cent50100"),LineColor(kViolet),LineStyle(1),LineWidth(2));
*/
  //// actual draw
  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.04);


  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  gPad->SetLeftMargin(0.1);
  myPlot->Draw();
  myPlot->GetYaxis()->SetTitleOffset(1.5);
  myPlot->GetYaxis()->SetTitle("");
  myPlot->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  myPlot->GetXaxis()->CenterTitle();
  myPlot->GetXaxis()->SetTitleSize(0.048);
  myPlot->GetXaxis()->SetTitleOffset(0.9);
 
  double xpos1 = 0.52;
  double xpos2 = 0.88;
  double ypos1 = 0.45;
  double ypos2 = 0.9;
  if(fver==5) ypos1 = 0.2;
  TLegend *leg = new TLegend(xpos1,ypos1,xpos2,ypos2);
  leg->SetHeader("Signal PDF");
  leg->SetTextSize(15);
  leg->SetTextFont(43);
  leg->SetBorderSize(0);
  for(int ib=0; ib<nBin; ib++){
    leg->AddEntry(myPlot->findObject(Form("Cent%s",str_name[ib])),Form("cent. %s%s",str_name_leg[ib],perc.Data()),"l");
  } 
  leg->Draw("same");
  c1->SaveAs(Form("FitInfoPlot/SigShape_cent_f%d.pdf",fver));
 
  
}

