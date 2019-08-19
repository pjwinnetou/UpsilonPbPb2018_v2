#ifndef COMMONOPT_Style_jaebeom_H
#define COMMONOPT_Style_jaebeom_H

#include <Riostream.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TGraphAsymmErrors.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLatex.h>
#include "TPaletteAxis.h"
#include <TCut.h>
#include <TString.h>
#include <TMath.h>
#include <math.h>

#include <TArrow.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TLorentzRotation.h>
#include <TVector3.h>
#include <TRandom.h>

////// calculation with error propagation

void DivideValue(Double_t num, Double_t numErr, Double_t den, Double_t denErr, Double_t* frac, Double_t* fracErr){
  *frac = num/den;
  *fracErr = (*frac) * TMath::Sqrt( TMath::Power(numErr/num,2) + TMath::Power(denErr/den,2) );
}

void MultiplyValue(Double_t a, Double_t aErr, Double_t b, Double_t bErr, Double_t* res, Double_t* resErr){
  *res = a*b;
  *resErr = (*res) * TMath::Sqrt( TMath::Power(aErr/a,2) + TMath::Power(bErr/b,2) );
}

void AddValue(Double_t a, Double_t aErr, Double_t b, Double_t bErr, Double_t* res, Double_t* resErr){
  *res = a+b;
  *resErr = 1 * TMath::Sqrt( TMath::Power(aErr,2) + TMath::Power(bErr,2) );
}

void SubtractValue(Double_t a, Double_t aErr, Double_t b, Double_t bErr, Double_t* res, Double_t* resErr){
  *res = a-b;
  *resErr = 1 * TMath::Sqrt( TMath::Power(aErr,2) + TMath::Power(bErr,2) );
}

////// draw lines

void dashedLine(Double_t x1=0,Double_t y1=0,Double_t x2=1,Double_t y2=1,Int_t color=1, Double_t width=1)
{
   TLine* t1 = new TLine(x1,y1,x2,y2);
   t1->SetLineWidth(width);
   t1->SetLineStyle(7);
   t1->SetLineColor(color);
   t1->Draw();
}

void solidLine(Double_t x1=0,Double_t y1=0,Double_t x2=1,Double_t y2=1,Int_t color=1, Double_t width=1)
{
  TLine* t1 = new TLine(x1,y1,x2,y2);
  t1->SetLineWidth(width);
  t1->SetLineStyle(1);
  t1->SetLineColor(color);
  t1->Draw();
}

////// SetStyle

void SetHistStyle(TH1* h, Int_t c, Int_t m) {
  Int_t colorArr[] = { kGray+3, kRed+2, kBlue+1, kOrange+7, kGreen+3, kAzure+9, kViolet-1, kGreen+1,kBlack };
  Int_t markerFullArr[] = {kFullCircle, kFullTriangleUp, kFullTriangleDown, kFullSquare, kFullStar, kFullDiamond};
  Int_t markerOpenArr[] = {kOpenCircle, kOpenTriangleUp, kOpenTriangleDown, kOpenSquare, kOpenStar, kOpenDiamond};
  h-> SetMarkerColor(colorArr[c]);
  if(m<10) h-> SetMarkerStyle(markerFullArr[m]);
  else h-> SetMarkerStyle(markerOpenArr[m-10]);
  h-> SetMarkerSize(1.0);
  h-> SetLineColor(colorArr[c]);
  h-> SetLineWidth(1.);
  h-> GetXaxis()->CenterTitle();
  h-> GetYaxis()->CenterTitle();
}

void SetHistStyleSmall(TH1* h, Int_t c, Int_t m) {
  Int_t colorArr[] = { kGray+3, kRed+1, kBlue+1, kOrange+7, kGreen+3, kAzure+9, kViolet-1, kGreen+1,kBlack };
  Int_t markerFullArr[] = {kFullCircle, kFullTriangleUp, kFullTriangleDown, kFullSquare, kFullStar, kFullDiamond};
  Int_t markerOpenArr[] = {kOpenCircle, kOpenTriangleUp, kOpenTriangleDown, kOpenSquare, kOpenStar, kOpenDiamond};
  h-> SetMarkerColor(colorArr[c]);
  if(m<10) h-> SetMarkerStyle(markerFullArr[m]);
  else h-> SetMarkerStyle(markerOpenArr[m-10]);
  h-> SetMarkerSize(0.3);
  h-> SetLineColor(colorArr[c]);
  h-> SetLineWidth(1.);
  h-> GetXaxis()->CenterTitle();
  h-> GetYaxis()->CenterTitle();
}

void SetHistStyle2D(TH2* h, Int_t c, Int_t m) {
  Int_t colorArr[] = { kGray+3, kRed+2, kBlue+1, kOrange+7, kGreen+3, kAzure+9, kViolet-1, kGreen+1,kBlack };
  Int_t markerFullArr[] = {kFullCircle, kFullTriangleUp, kFullTriangleDown, kFullSquare, kFullStar, kFullDiamond};
  Int_t markerOpenArr[] = {kOpenCircle, kOpenTriangleUp, kOpenTriangleDown, kOpenSquare, kOpenStar, kOpenDiamond};
/*  h-> SetMarkerColor(colorArr[c]);
  if(m<10) h-> SetMarkerStyle(markerFullArr[m]);
  else h-> SetMarkerStyle(markerOpenArr[m-10]);
  h-> SetMarkerSize(1.2);
  h-> SetLineColor(colorArr[c]);
  h-> SetLineWidth(1.);
*/  h-> GetXaxis()->CenterTitle();
  h-> GetYaxis()->CenterTitle();
}

void SetHistStyle2(TH1* h, Int_t c, Int_t m) {
  Int_t colorArr[] = { kGray+3, kRed+2, kBlue+1, kOrange+7, kAzure+9, kViolet-4, kGreen+1, kPink+9, kBlack };
  Int_t markerFullArr[] = {kFullCircle, kFullTriangleUp, kFullTriangleDown, kFullSquare, kFullStar, kFullDiamond};
  //Int_t markerOpenArr[] = {kOpenCircle, kOpenTriangleUp, kOpenTriangleDown, kOpenSquare, kOpenStar, kOpenDiamond};
  Int_t markerOpenArr[] = {kOpenCircle, kOpenTriangleUp, kOpenTriangleDown, kOpenSquare, kOpenCircle, kOpenTriangleUp, kOpenTriangleDown, kOpenSquare};

  h-> SetMarkerColor(colorArr[c]);
  if(m<10) h-> SetMarkerStyle(markerFullArr[m]);
  else h-> SetMarkerStyle(markerOpenArr[m-10]);
  h-> SetMarkerSize(1.2);
  h-> SetLineColor(colorArr[c]);
  h-> SetLineWidth(1);
}

void SetGraphStyle(TGraph* gr, Int_t c, Int_t m) {
  //Int_t colorArr[] = { kPink-6, kBlue-2, kGreen+3, kViolet-6, kBlack };
  Int_t colorArr[] = { kPink-6, kBlue-3, kGreen+2, kViolet-6, kBlack };
  Int_t markerFullArr[] = {kFullCircle, kFullSquare, kFullDiamond, kFullCross};
  Int_t markerOpenArr[] = {kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross};
  gr-> SetMarkerColor(colorArr[c]);
  if(m<4) gr-> SetMarkerStyle(markerFullArr[m]);
  else gr-> SetMarkerStyle(markerOpenArr[m-4]);
  if (m==2 || m==3) gr-> SetMarkerSize(2.6);
  else gr-> SetMarkerSize(1.5);
  gr-> SetLineColor(colorArr[c]);
  gr-> SetLineWidth(1);
}

void SetGraphStyle2(TGraph* gr, Int_t c, Int_t m) {
  Int_t colorArr[] = { kGray+3, kRed+2, kBlue+1, kOrange+7, kGreen+3, kAzure+9, kViolet-1, kGreen+1,kBlack };
  Int_t markerFullArr[] = {kFullCircle, kFullTriangleUp, kFullTriangleDown, kFullSquare, kFullStar, kFullDiamond};
  Int_t markerOpenArr[] = {kOpenCircle, kOpenTriangleUp, kOpenTriangleDown, kOpenSquare, kOpenStar, kOpenDiamond};
  gr-> SetMarkerColor(colorArr[c]);
  if(m<4) gr-> SetMarkerStyle(markerFullArr[m]);
  else gr-> SetMarkerStyle(markerOpenArr[m-4]);
  if (m==2 || m==3) gr-> SetMarkerSize(1.5);
  else gr-> SetMarkerSize(1.5);
  gr-> SetLineColor(colorArr[c]);
  gr-> SetLineWidth(1.);
  gr-> GetXaxis()->CenterTitle();
  gr-> GetYaxis()->CenterTitle();
}

void SetGraphStyleSmall(TGraph* gr, Int_t c, Int_t m) {
  Int_t colorArr[] = { kGray+3, kRed+2, kBlue+1, kOrange+7, kGreen+3, kAzure+9, kViolet-1, kGreen+1,kBlack };
  Int_t markerFullArr[] = {kFullCircle, kFullTriangleUp, kFullTriangleDown, kFullSquare, kFullStar, kFullDiamond};
  Int_t markerOpenArr[] = {kOpenCircle, kOpenTriangleUp, kOpenTriangleDown, kOpenSquare, kOpenStar, kOpenDiamond};
  gr-> SetMarkerColor(colorArr[c]);
  if(m<4) gr-> SetMarkerStyle(markerFullArr[m]);
  else gr-> SetMarkerStyle(markerOpenArr[m-4]);
  if (m==2 || m==3) gr-> SetMarkerSize(2.6);
  else gr-> SetMarkerSize(0.8);
  gr-> SetLineColor(colorArr[c]);
  gr-> SetLineWidth(1.);
  gr-> GetXaxis()->CenterTitle();
  gr-> GetYaxis()->CenterTitle();
}

void SetGraphStyleOpen(TGraph* gr, Int_t c, Int_t m, Int_t s) {
  //Int_t colorArr[] = { kPink-6, kBlue-2, kGreen+3, kViolet-6, kBlack };
  Int_t colorArr[] = { kPink-6, kBlue-3, kGreen+2, kViolet-6, kBlack };
  Int_t markerFullArr[] = {kFullCircle, kFullSquare, kFullDiamond, kFullCross};
  Int_t markerOpenArr[] = {kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross};
  gr-> SetMarkerColor(colorArr[c]);
  gr-> SetMarkerStyle(markerOpenArr[s]);
  if (m==2 || m==3) gr-> SetMarkerSize(2.6);
  else gr-> SetMarkerSize(1.5);
  gr-> SetLineColor(colorArr[c]);
  gr-> SetLineWidth(1);
}

void SetGraphStyleSys(TGraph* gr, Int_t c) {
  //Int_t colorArr[] = { kPink-6, kBlue-2, kGreen+3, kViolet-6, kBlack };
  Int_t colorArr[] = { kPink-6, kBlue-3, kGreen+2, kViolet-6, kBlack };
  Int_t colorArrSys[] = { kRed-10, kBlue-10, kGreen-10, kMagenta-10, kGray};
  gr-> SetLineColor(colorArr[c]);
  gr-> SetFillColorAlpha(colorArrSys[c],0.5);
  gr-> SetLineWidth(1);
}


void SetLegendStyle(TLegend* l) {
  l->SetFillColor(0);
  l->SetFillStyle(4000);
  l->SetBorderSize(0);
  l->SetMargin(0.2);
  l->SetTextSize(0.029);
  l->SetTextFont(42);
}

void SetTextStyle(TPaveText* t) {
  t->SetFillColor(0);
  t->SetFillStyle(4000);
  t->SetBorderSize(0);
  t->SetMargin(0.2);
  t->SetTextSize(0.035);
}


#endif
