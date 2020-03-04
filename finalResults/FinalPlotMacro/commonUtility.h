#ifndef COMMONUTILITY_Yongsun_H
#define COMMONUTILITY_Yongsun_H
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <TString.h>

#include <ctime>

#include <TROOT.h>
#include "TSystem.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>
#include "TRandom3.h"
#include "TClonesArray.h"
#include <TAxis.h>
#include <cmath>
#include <TLorentzRotation.h>
#include <TMath.h>
#include <math.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include "TRandom3.h"
#include <TCanvas.h>
#include <TPaveStats.h>

#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TDatime.h>
#include <TStyle.h>
#include <TPaveStats.h>


struct kinem {
  double pt;
  double eta;
  double phi;
};

int debugIt=1;

void changeLine() {
  cout << "#################################################" << endl<< endl << endl;
}

void debug() {
  cout << " debugging..line" << debugIt << endl;
  debugIt++;
}


bool acceptance(double muPt, double muEta) {
  return ( 	( (TMath::Abs(muEta) <= 1.2) && (muPt >=3.5) ) || 
		( (TMath::Abs(muEta) > 1.2)  && (TMath::Abs(muEta) <= 2.1) && (muPt >= 5.47-1.89*(TMath::Abs(muEta))) ) || 
		( (TMath::Abs(muEta) > 2.1)  && (TMath::Abs(muEta) <= 2.4) && (muPt >= 1.5) )  
		) ; 
}

bool CaloMatchingCut(float trkPt, float trkEta, float pfEcal, float pfHcal)
{
  return ( trkPt<20 || ((pfEcal+pfHcal)/TMath::CosH(trkEta)>0.2*trkPt && (pfEcal+pfHcal)/TMath::CosH(trkEta)>trkPt-80) );
}




void cleverCanvasSaving(TCanvas* c, TString s,TString format="gif") {
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
}

double getDPHI( double phi1, double phi2) {
  double dphi = phi1 - phi2;
 
  if ( dphi > 3.141592653589 )
    dphi = dphi - 2. * 3.141592653589;
  if ( dphi <= -3.141592653589 ) 
    dphi = dphi + 2. * 3.141592653589;
  
  if ( fabs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << endl;
  }
  
  return dphi;
}

double getAbsDphi( double phi1, double phi2) {
  return fabs( getDPHI(phi1, phi2) ) ;
}


double getDR( double eta1, double phi1, double eta2, double phi2){ 
  double theDphi = getDPHI( phi1, phi2);
  double theDeta = eta1 - eta2;
  return sqrt ( theDphi*theDphi + theDeta*theDeta);
}

void divideWOerr( TH1* h1, TH1* h2) {  //by Yongsun Jan 26 2012                                                                              
  if ( h1->GetNbinsX() != h2->GetNbinsX() ) {
    cout << " different bin numbers!!" << endl;
    return;
  }

  for ( int i=1 ; i<=h1->GetNbinsX() ; i++) {
    if ( (h2->GetBinContent(i) == 0 ) ) {
      h1->SetBinContent(i, 0);
      h1->SetBinError (i, 0);
    }
    else {
      float newV = h1->GetBinContent(i)/ h2->GetBinContent(i);
      float newE = h1->GetBinError(i)  / h2->GetBinContent(i);
      h1->SetBinContent(i, newV);
      h1->SetBinError  (i, newE);
    }
  }
}



void AddBinError( TH1* h=0, int binNumber=0 , float val=0){
  float valO = h->GetBinError(binNumber);
  float newVal = sqrt( valO*valO + val*val);
  h->SetBinError(binNumber,newVal);
}

void drawSys(TH1 *h,double *sys, int theColor= kYellow, int fillStyle = -1, int lineStyle = -1)
{
   for (int i=1;i<=h->GetNbinsX();i++)
      {
	 double val = h->GetBinContent(i);
	 double err = val * sys[i-1];
	 TBox *b = new TBox(h->GetBinLowEdge(i),val-err,h->GetBinLowEdge(i+1),val+err);
	 //      b->SetFillStyle(3001);                                                                                                                                                                              
	 b->SetLineColor(theColor);
	 b->SetFillColor(theColor);
	 if ( fillStyle > -1 ) b->SetFillStyle(fillStyle);
	 if ( lineStyle > -1 ) b->SetLineStyle(lineStyle);
	 
	 b->Draw();
      }
}

void drawSys(TGraph *h, double *sys, double width=5, int theColor= kYellow, int fillStyle = -1, int lineStyle = -1)
{
  for (int i=0;i<h->GetN();i++)
    {
      double val;
      double theX;
      h->GetPoint(i,theX,val);
      double err = val * sys[i];
      TBox *b = new TBox(theX-width,val-err,theX+width,val+err);
      
      b->SetLineColor(theColor);
      b->SetFillColor(theColor);
      if ( fillStyle > -1 ) b->SetFillStyle(fillStyle);
      if ( lineStyle > -1 ) b->SetLineStyle(lineStyle);
      
      b->Draw();
    }
}



void drawSysAbs(TH1 *h,TH1 *sys, int theColor= kYellow, int fillStyle = -1, int lineStyle = -1)
{
   for (int i=1;i<=h->GetNbinsX();i++)
      {
         double val = h->GetBinContent(i);
         double err = fabs(sys->GetBinContent(i));
	 if (err == 0  ) continue;
	 TBox *b = new TBox(h->GetBinLowEdge(i),val-err,h->GetBinLowEdge(i+1),val+err);
         b->SetLineColor(theColor);
         b->SetFillColor(theColor);
         if ( fillStyle > -1 ) b->SetFillStyle(fillStyle);
         if ( lineStyle > -1 ) b->SetLineStyle(lineStyle);

         b->Draw();
      }
}

void integerizeTH1(TH1* h1) {
   for ( int j = 0 ; j <= h1->GetNbinsX()+1 ; j++) {
      float vTemp = h1->GetBinContent(j);
      h1->SetBinContent( j , (int)vTemp ) ;
   }
}

void multiplyBonA(TH1* h1, TH1* h2) {
   if ( h1->GetNbinsX() != h2->GetNbinsX() ) { 
      cout << " different bin numbers.." << endl ;
      return ;
   }
   
   for ( int j=1 ; j<= h1->GetNbinsX(); j++) {
      float v1 = h1->GetBinContent(j);
      float e1 = h1->GetBinError(j);
      float v2 = h2->GetBinContent(j);
      //      float e2 = h2->GetBinError(j);
      
      float v3 = v1 * v2;
      float e3 = e1 * v2;
      h1->SetBinContent(j, v3);
      h1->SetBinError(j, e3);
   }
}

void drawPatch(float x1, float y1, float x2, float y2, int color =0, int style=1001, TString ops=""){
  TLegend *t1=new TLegend(x1,y1,x2,y2,NULL,ops.Data());
   if ( color ==0) t1->SetFillColor(kWhite);
   else t1->SetFillColor(color);
   t1->SetBorderSize(0);
   t1->SetFillStyle(style);
   t1->Draw("");
}

void drawErrorBox(float x1,float y1, float x2, float y2, int theColor=kSpring+8)
{
   TBox* tt = new TBox(x1,y1,x2,y2);
   tt->SetFillColor(theColor);
   tt->SetFillStyle(3001);
   tt->Draw();
}

void drawErrorBand(TH1* h, double* err, int theColor=kSpring+8)
{
   for ( int j=1 ; j<= h->GetNbinsX()-1; j++) {
      double theCont = h->GetBinContent(j);
      double theErr  = err[j] * h->GetBinContent(j);
      TBox* tt = new TBox(h->GetBinLowEdge(j), theCont-theErr, h->GetBinLowEdge(j)+ h->GetBinWidth(j), theCont+theErr);
      tt->SetFillColor(theColor);
      tt->SetFillStyle(3001);
      tt->Draw();
   }
}



void drawText(const char *text, float xp, float yp, int textColor=kBlack, int textSize=18){
   TLatex *tex = new TLatex(xp,yp,text);
   tex->SetTextFont(43);
   //   if(bold)tex->SetTextFont(43);
   tex->SetTextSize(textSize);
   tex->SetTextColor(textColor);
   tex->SetLineWidth(1);
   tex->SetNDC();
   tex->Draw();
}

void drawText2(const char *text, float xp, float yp, int textSize=18){
   TLatex *tex = new TLatex(xp,yp,text);
   tex->SetTextFont(43);
   tex->SetTextSize(textSize);
   tex->SetTextColor(kBlack);
   tex->SetLineWidth(1);
   tex->Draw();
}

void jumSun(double x1=0,double y1=0,double x2=1,double y2=1,int color=1, double width=1)
{
   TLine* t1 = new TLine(x1,y1,x2,y2);
   t1->SetLineWidth(width);
   t1->SetLineStyle(7);
   t1->SetLineColor(color);
   t1->Draw();
}


void onSun(double x1=0,double y1=0,double x2=1,double y2=1,int color=1, double width=1)
{
  TLine* t1 = new TLine(x1,y1,x2,y2);
  t1->SetLineWidth(width);
  t1->SetLineStyle(1);
  t1->SetLineColor(color);
  t1->Draw();
}
void regSun(double x1=0,double y1=0,double x2=1,double y2=1,int color=1, double width=1)
{
   TLine* t1 = new TLine(x1,y1,x2,y2);
   t1->SetLineWidth(width);
   t1->SetLineStyle(1);
   t1->SetLineColor(color);
   t1->Draw();
}

void mcStyle1(TH1* h=0) {
   h->SetLineColor(kRed);
   h->SetFillColor(kRed-9);
   h->SetFillStyle(3004);
}
void mcStyle2(TH1* h=0) {
   h->SetLineColor(kBlue);
   h->SetFillColor(kAzure-8);
   h->SetFillStyle(3005);
}

void mcStyle3(TH1* h=0) {
  h->SetLineColor(kBlue);
  h->SetFillColor(kAzure-8);
  h->SetFillStyle(3005);
  h->SetMarkerSize(0);
}

void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns,
                          const Int_t rows, const Float_t leftOffset=0.,
                          const Float_t bottomOffset=0.,
                          const Float_t leftMargin=0.2,
                          const Float_t bottomMargin=0.2,
                          const Float_t edge=0.05) {
  if (canv==0) {
    cout << "make MultiPanelCanvas    Got null canvas." << endl<<endl<<endl<<endl;
    return;
  }
   canv->Clear();

   TPad* pad[columns][rows];
   
   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth =
      (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
			   (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
      (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			  (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
   
   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
   
   for(Int_t i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
   Int_t ct = 0;
   for(Int_t i=rows-2;i>0;i--) {
      Ylow[i] = Yup[rows-1] + ct*PadHeight;
      Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
      ct++;
   }

   TString padName;
   for(Int_t i=0;i<columns;i++) {
      for(Int_t j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
			      Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);
	 
         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);
	 
         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}

void makeEfficiencyCanvas(TCanvas*& canv, const Int_t columns,
                          const Float_t leftOffset=0.05,
                          const Float_t bottomOffset=0.01,
                          const Float_t leftMargin=0.25,
                          const Float_t bottomMargin=0.25,
                          const Float_t edge=0.01) {

  
  const Int_t rows = 2;
  
  if (canv==0) {
    cout <<"makeMultiPanelCanvas              Got null canvas." <<endl<<endl<<endl;
    return;
  }
  canv->Clear();
  
   TPad* pad[columns][rows];
   
   Float_t Xlow[columns];
   Float_t Xup[columns];
   Float_t Ylow[rows];
   Float_t Yup[rows];
   Float_t PadWidth =
      (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
			   (1.0/(1.0-edge))+(Float_t)columns-2.0);
   Float_t PadHeight =
      (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			  (1.0/(1.0-edge))+(Float_t)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
   
   Yup[0] = 1;
   Ylow[rows-1] = bottomOffset;
  
   Ylow[0] = 1.0-PadHeight/(1.0-edge) - 0.2;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin) - 0.2;
   
   for(Int_t i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
 

   TString padName;
   for(Int_t i=0;i<columns;i++) {
      for(Int_t j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
			      Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);
	 
         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);
	 
         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}


void twikiSave(TCanvas* c=0, TString  name="",int w=0,int h=0)
{
   if ( w==0) w = c->GetWindowWidth();
   if ( h==0) h = c->GetWindowHeight();
   
   c->SaveAs(name.Data());
   cout << Form(" <br/>   <img src=\"%%ATTACHURLPATH%%/%s\" alt=\"%s\" width='%d' height='%d'/>",name.Data(),name.Data(),w,h)<< endl;
}

void centralityBinning(float *b=0)
{
  b[0]=      0;
  b[1]=  5.045;
  b[2]=  7.145;
  b[3]=  8.755;
  b[4]= 10.105;
  b[5]=  11.294;
  b[6]= 12.373;
  b[7]=  13.359;
  b[8]= 14.283;
  b[9]=  15.202;
  b[10] = 100;
 }
void handsomeTH2( TH2 *a=0)
{
   a->GetYaxis()->SetTitleOffset(1.25);
   a->GetXaxis()->CenterTitle();
   a->GetYaxis()->CenterTitle();
}

void handsomeTG1( TGraphErrors *a=0, int col =1, float size=1, int markerstyle=20)
{
  a->SetMarkerColor(col);
  a->SetMarkerSize(size);
  a->SetMarkerStyle(markerstyle);
  a->SetLineColor(col);
  a->GetYaxis()->SetTitleOffset(1.25);
  a->GetXaxis()->CenterTitle();
  a->GetYaxis()->CenterTitle();
}

void handsomeTH1( TH1 *a=0, int col =1, float size=1, int markerstyle=20)
{
  a->SetMarkerColor(col);
  a->SetMarkerSize(size);
  a->SetMarkerStyle(markerstyle);
  a->SetLineColor(col);
  a->GetYaxis()->SetTitleOffset(1.25);
  a->GetXaxis()->CenterTitle();
  a->GetYaxis()->CenterTitle();
}
void fixedFontAxis(TGaxis * ax)
{
   ax->SetLabelFont(43);
   ax->SetLabelOffset(0.01);
   ax->SetLabelSize(22);
   ax->SetTitleFont(43);
   ax->SetTitleSize(14);
   ax->SetTitleOffset(2);
}

void fixedFontHist(TH1 * h, Float_t xoffset=1.3, Float_t yoffset=1.2)
{
   h->SetLabelFont(43,"X");
   h->SetLabelFont(43,"Y");
   //h->SetLabelOffset(0.01);
   h->SetLabelSize(22);
   h->SetTitleFont(43);
   h->SetTitleSize(16);
   h->SetLabelSize(18,"Y");
   h->SetLabelSize(18,"X");
   h->SetTitleFont(43,"Y");
   h->SetTitleSize(20,"Y");
   h->SetTitleSize(20,"X");
   h->SetTitleOffset(xoffset,"X");
   h->SetTitleOffset(yoffset,"Y");
   h->GetXaxis()->CenterTitle();
   h->GetYaxis()->CenterTitle();
}

void handsomeTH1Fill( TH1 *a=0, int nFill=1) {
   handsomeTH1(a,nFill);
   a->SetFillColor(nFill);
}
void handsomeTGraph(TGraphAsymmErrors* a, int col=1)
{
   a->SetLineColor(col);
   a->SetMarkerColor(col);
   a->SetMarkerSize(1);
   a->SetMarkerStyle(20);
}

void TH1ScaleByWidth(TH1* h=0) 
{
   int nBins = h->GetNbinsX();
   // cout << "Start scaling by width" << endl;
   for ( int j=1; j<=nBins ;j++)
      {
         double theWidth = h->GetBinWidth(j);
         //      cout << "width = " << theWidth << ",   " ;                                                                                  
	 double cont = h->GetBinContent(j);
         double err =  h->GetBinError(j);
         h->SetBinContent(j, cont/theWidth);
         h->SetBinError  (j, err/theWidth);
      }
   //   cout << endl;
}

void scaleInt( TH1 *a=0, double normFac=1., double minX=-999.21231, double maxX=-999.21231)
{
  float fac=0;
  int lowBin=1; 
  int highBin=a->GetNbinsX();
  if ( minX != -999.21231)
    lowBin = a->FindBin(minX);
  if ( maxX != -999.21231)
    highBin=a->FindBin(maxX);

  fac =  a->Integral( lowBin, highBin);
  if ( fac>0) a->Scale(normFac/fac);
}

void scaleIntWidth( TH1 *a=0, double normFac=1., double minX=-999.21231, double maxX=-999.21231)
{
  float fac=0;
  int lowBin=1; 
  int highBin=a->GetNbinsX();
  if ( minX != -999.21231)
    lowBin = a->FindBin(minX);
  if ( maxX != -999.21231)
    highBin=a->FindBin(maxX);

  fac =  a->Integral( lowBin, highBin);
  if ( fac>0) a->Scale(normFac/fac);
  
  int nBins = a->GetNbinsX();
  for(int j=1;j<=nBins;j++){
    double theWidth = a->GetBinWidth(j);
    double cont = a->GetBinContent(j);
    double err = a->GetBinError(j);
    a->SetBinContent(j, cont/theWidth);
    a->SetBinError(j, err/theWidth);
  }

}



double goodIntegral( TH1 *a=0, int lower=-123, int upper=-123)
{
   int nBins = a->GetNbinsX();
   
   if ( (lower==-123) || (upper==-123)) {
      lower = 1;
      upper = nBins;
   }
   double tempInt=0;
   for ( int j=lower; j<=upper; j++) {
      tempInt = tempInt + a->GetBinContent(j) * a->GetBinWidth(j);
   }
   return tempInt;
}

double goodIntegralError( TH1 *a=0, int lower=-123, int upper=-123)
{
   int nBins = a->GetNbinsX();
   if ( (lower==-123) || (upper==-123)) {
      lower = 1;
      upper = nBins;
   }
   
   double tempInt=0;
   for ( int j=lower; j<=upper; j++) {
      tempInt = tempInt + a->GetBinError(j) * a->GetBinError(j) * a->GetBinWidth(j) *  a->GetBinWidth(j);
   }
   return sqrt(tempInt);
}




void handsomeTH1Sumw2( TH1 *a=0, int col =1, float size=1, int markerstyle=20)
{
   handsomeTH1(a,col,size,markerstyle);
   a->Sumw2();
}



void handsomeTH1N( TH1 *a=0, int col =1)
{
   handsomeTH1(a,col);
   a->Scale(1./a->GetEntries());
}


void handsomeTH1OnlyColor( TH1 *a=0, int col =1)
{
   a->SetMarkerColor(col);
   a->SetLineColor(col);
   a->GetYaxis()->SetTitleOffset(1.25);
}


void easyLeg( TLegend *a=0 , TString head="")
{
  a->SetBorderSize(0);
  a->SetHeader(head);
  a->SetTextFont(42);
  //  a->SetTextSize(17);
  a->SetLineColor(1);
  a->SetLineStyle(1);
  a->SetLineWidth(1);
  a->SetFillColor(0);
  a->SetFillStyle(0);

}



double cleverRange(TH1* h,float fac=1.2, float minY=1.e-3)
{
   float maxY =  fac * h->GetBinContent(h->GetMaximumBin());
   //   cout <<" range will be set as " << minY << " ~ " << maxY << endl; 
   h->SetAxisRange(minY,maxY,"Y");
   return maxY;
}


double getCleverRange(TH1* h)
{
  double maxY = -1000000;
  for ( int ibin = 1 ; ibin <= h->GetNbinsX() ; ibin++) {
    if (maxY < h->GetBinContent(ibin) ) 
      maxY = h->GetBinContent(ibin);
  }
  return maxY;
}

double cleverRange(TH1* h,TH1* h2, float fac=1.2, float minY=1.e-3)
{
  float maxY1 =  fac * h->GetBinContent(h->GetMaximumBin());
  float maxY2 =  fac * h2->GetBinContent(h2->GetMaximumBin());
  
  //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;                                                                    
  h->SetAxisRange(minY,max(maxY1,maxY2),"Y");
  h2->SetAxisRange(minY,max(maxY1,maxY2),"Y");
  return max(maxY1,maxY2);
}

void cleverRangeLog(TH1* h,float fac=1.2, float theOrder=1.e-4)
{
  float maxY =  fac * h->GetBinContent(h->GetMaximumBin());
  //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;                                               
  h->SetAxisRange(maxY * theOrder,maxY,"Y");
}



TF1* cleverGaus(TH1* h, TString title="h", float c = 2.5, bool quietMode=true)
{
   if ( h->GetEntries() == 0 )
      {
	TF1 *fit0 = new TF1(title.Data(),"gaus",-1,1);
	 fit0->SetParameters(0,0,0);
	 return fit0;
      }
   
   int peakBin  = h->GetMaximumBin();
   double peak =  h->GetBinCenter(peakBin);
   double sigma = h->GetRMS();
  
   TF1 *fit1 = new TF1(title.Data(),"gaus",peak-c*sigma,peak+c*sigma);
   if (quietMode) h->Fit(fit1,"LL M O Q R");
   else    h->Fit(fit1,"LL M O Q R");
   return fit1;
}




void drawCMS(float px, float py, float nLumi) {
   TLatex *cms = new TLatex(px,py,"CMS Preliminary");
   cms->SetTextFont(63);
   cms->SetTextSize(15);
   cms->SetNDC();
   cms->Draw();
   TLatex *lumi = new TLatex(px+0.35,py,Form("#intL dt = %.1f #mub^{-1}",nLumi));
   lumi->SetTextFont(63);
   lumi->SetTextSize(13);
   lumi->SetNDC();
   lumi->Draw();
}

void drawCMSpp(float px, float py, float nLumi) {
  TLatex *cms = new TLatex(px,py,"CMS Preliminary");
  cms->SetTextFont(63);
  cms->SetTextSize(15);
  cms->SetNDC();
  cms->Draw();
  TLatex *lumi = new TLatex(px+0.35,py,Form("#intL dt = %.1f #nb^{-1}",nLumi));
  lumi->SetTextFont(63);
  lumi->SetTextSize(13);
  lumi->SetNDC();
  lumi->Draw();
}



void drawCMS2(float px, float py, float nLumi, int textSize=15) {
   TLatex *cms = new TLatex(px,py,"CMS Preliminary");
   cms->SetTextFont(63);
   cms->SetTextSize(textSize);
   cms->SetNDC();
   cms->Draw();
   TLatex *lumi = new TLatex(px,py-0.08,Form("#intL dt = %.1f #mub^{-1}",nLumi));
   lumi->SetTextFont(63);
   lumi->SetTextSize(textSize-2);
   lumi->SetNDC();
   lumi->Draw();
}
void drawCMS2011(float px, float py, float nLumi=3.8, int textSize=15) {
   TLatex *cms = new TLatex(px,py,"CMS Preliminary");
   cms->SetTextFont(63);
   cms->SetTextSize(textSize);
   cms->SetNDC();
   cms->Draw();
   TLatex *lumi = new TLatex(px,py-0.08,Form("L_{Int} dt = %.0f #mub^{-1}",nLumi));
   lumi->SetTextFont(63);
   lumi->SetTextSize(textSize-2);
   lumi->SetNDC();
   lumi->Draw();
}

void drawCMS3(float px, float py, float nLumi, int textSize=15) {
   TLatex *cms = new TLatex(px,py,"CMS");
   cms->SetTextFont(63);
   cms->SetTextSize(textSize);
   cms->SetNDC();
   cms->Draw();
   TLatex *lumi = new TLatex(px,py-0.08,Form("PbPb  #sqrt{s_{NN}}=2.76TeV  #intL dt = %.0f #mub^{-1}",nLumi));
   lumi->SetTextFont(63);
   lumi->SetTextSize(textSize-2);
   lumi->SetNDC();
   lumi->Draw();
}
void drawCMS4(float px, float py, float nLumi, int textSize=15) {
   TLatex *cms = new TLatex(px,py,"CMS Preliminary");
   cms->SetTextFont(63);
   cms->SetTextSize(textSize);
   cms->SetNDC();
   cms->Draw();
   
   TLatex *pbpb = new TLatex(px,py-0.08,"PbPb  #sqrt{s_{NN}}=2.76TeV");
   pbpb->SetTextFont(63);
   pbpb->SetTextSize(textSize-2);
   pbpb->SetNDC();
   pbpb->Draw();
   
   TLatex *lumi = new TLatex(px,py-0.16,Form("#intL dt = %.1f #mub^{-1}",nLumi));
   lumi->SetTextFont(63);
   lumi->SetTextSize(textSize-2);
   lumi->SetNDC();
   lumi->Draw();
   
   
}



void getNiceBins( TH1* h=0, int nDiv=4) {
   int nBins = h->GetNbinsX();
   double allInt = h->Integral(1,nBins);
   cout<< " All integral = " << allInt<< endl;

   TH1F* hacc = (TH1F*)h->Clone(Form("%s_accu",h->GetName()));
   hacc->Reset();
   double acc =0;
   

   int j=0;
   for ( int i=1 ; i<= nBins ; i++ ) {
      acc = acc + h->GetBinContent(i);
      hacc->SetBinContent( i , acc ) ;
       if ( ( hacc->GetBinContent(i) > j*allInt/nDiv ) && ( hacc->GetBinContent(i-1) <= j*allInt/nDiv)) {
	 j++;
	 cout << j << "th bin = " << i <<  "    value  = " << h->GetBinCenter(i) << endl;
       }
   }     
   cout << " acc = " << acc << endl;
   TCanvas* c1 = new TCanvas(Form("%s_c1Acc",h->GetName()),"",400,400);
   c1->Draw();
   handsomeTH1(hacc,1);
   handsomeTH1(h,2);
   hacc->DrawCopy();
   h->DrawCopy("same");
  
   
}

void stripErr(TH1* theHist=0) {
   for ( int ibin = 1 ; ibin <= theHist->GetNbinsX() ; ibin++)
      theHist->SetBinError(ibin,0);
}


double getPolyArea(TH1* h1, TH1* h2, double minX, double maxX) { 
  if ( (h1->GetNbinsX() != h2->GetNbinsX()) || ( h1->FindBin(minX) != h2->FindBin(minX)) || (h1->FindBin(maxX) != h2->FindBin(maxX)) ) 
    {
      cout <<" binnings are not matched!!! " << endl;
      return -1;
    }
  
  int minBin = h1->FindBin(minX);
  int maxBin = h1->FindBin(maxX);
  cout << " getPolyArea : from " << h1->GetBinLowEdge(minBin) << " to " << h1->GetBinLowEdge(maxBin+1) << endl;
  double area = 0;
  for ( int ibin = minBin ; ibin <= maxBin ; ibin++) {
    area = area + ( h2->GetBinContent(ibin) - h1->GetBinContent(ibin) ) * h1->GetBinWidth(ibin) ;
  }
  
  return area;
}

double getPolyAreaErr(TH1* h1, TH1* h2, double minX, double maxX) {
  if ( (h1->GetNbinsX() != h2->GetNbinsX()) || ( h1->FindBin(minX) != h2->FindBin(minX)) || (h1->FindBin(maxX) != h2->FindBin(maxX)) )
    {
      cout <<" binnings are not matched!!! " << endl;
      return -1;
    }

  int minBin = h1->FindBin(minX);
  int maxBin = h1->FindBin(maxX);
  double err2 = 0;
  for ( int ibin = minBin ; ibin <= maxBin ; ibin++) {
    double addErr2 = (h2->GetBinError(ibin)*h2->GetBinError(ibin) + h1->GetBinError(ibin)*h1->GetBinError(ibin)) * (h1->GetBinWidth(ibin)*h1->GetBinWidth(ibin)) ;
    err2 = err2 + addErr2;
  }
  return sqrt(err2);
}

TH1D* getShiftedTH1D(TH1D* h1, double shift) {
  double newBinning[200];
  for ( int ibin = 0 ; ibin <= h1->GetNbinsX() ; ibin++) {
    newBinning[ibin] = h1->GetBinLowEdge(ibin+1) + shift ;
  }
  TH1D* res = new TH1D(Form("%s_shiftedby%f",h1->GetName(),(float)shift), "", h1->GetNbinsX(), newBinning);
  for ( int ibin = 0 ; ibin <= h1->GetNbinsX() ; ibin++) {
    res->SetBinContent( ibin, h1->GetBinContent(ibin));
    res->SetBinError(ibin, h1->GetBinError(ibin));
  }
  return res;
}



TString getDateAndTime()
{
  time_t currentTime;
  struct tm *localTime;

  time( &currentTime );                   // Get the current time
  localTime = localtime( &currentTime );  // Convert the current time to the local time

  int Day    = localTime->tm_mday;
  int Month  = localTime->tm_mon + 1;
  int Year   = localTime->tm_year + 1900;
  int Hour   = localTime->tm_hour;
  int Min    = localTime->tm_min;
  //  int Sec    = localTime->tm_sec;
  return Form("%d%d%d%d%d",Year,Month,Day,Hour,Min);
}




float getNcollFrom40Bin(int cBin) { 
  if (cBin == 0) return  1747.86 ;
  if (cBin == 1) return  1567.53 ;
  if (cBin == 2) return  1388.39 ;
  if (cBin == 3) return  1231.77 ;
  if (cBin == 4) return  1098.2 ;
  if (cBin == 5) return  980.439 ;
  if (cBin == 6) return  861.609 ;
  if (cBin == 7) return  766.042 ;
  if (cBin == 8) return  676.515 ;
  if (cBin == 9) return  593.473 ;
  if (cBin == 10) return  521.912 ;
  if (cBin == 11) return  456.542 ;
  if (cBin == 12) return  398.546 ;
  if (cBin == 13) return  346.647 ;
  if (cBin == 14) return  299.305 ;
  if (cBin == 15) return  258.344 ;
  if (cBin == 16) return  221.216 ;
  if (cBin == 17) return  188.677 ;
  if (cBin == 18) return  158.986 ;
  if (cBin == 19) return  134.7 ;
  if (cBin == 20) return  112.547 ;
  if (cBin == 21) return  93.4537 ;
  if (cBin == 22) return  77.9314 ;
  if (cBin == 23) return  63.5031 ;
  if (cBin == 24) return  52.0469 ;
  if (cBin == 25) return  42.3542 ;
  if (cBin == 26) return  33.9204 ;
  if (cBin == 27) return  27.3163 ;
  if (cBin == 28) return  21.8028 ;
  if (cBin == 29) return  17.2037 ;
  if (cBin == 30) return  13.5881 ;
  if (cBin == 31) return  10.6538 ;
  if (cBin == 32) return  8.35553 ;
  if (cBin == 33) return  6.40891 ;
  if (cBin == 34) return  5.13343 ;
  if (cBin == 35) return  3.73215 ;
  if (cBin == 36) return  3.06627 ;
  if (cBin == 37) return  2.41926 ;
  if (cBin == 38) return  2.11898 ;
  if (cBin == 39) return  1.76953 ;
  return -100000;
}


float  getNpart(int ibin) {
  if (ibin ==0) return  393.633;
  if (ibin ==1) return  368.819;
  if (ibin ==2) return  343.073;
  if (ibin ==3) return  317.625;
  if (ibin ==4) return  292.932;
  if (ibin ==5) return  271.917;
  if (ibin ==6) return  249.851;
  if (ibin ==7) return  230.72;
  if (ibin ==8) return  212.465;
  if (ibin ==9) return  194.752;
  if (ibin ==10) return  178.571;
  if (ibin ==11) return  163.23;
  if (ibin ==12) return  149.187;
  if (ibin ==13) return  136.011;
  if (ibin ==14) return  123.414;
  if (ibin ==15) return  111.7;
  if (ibin ==16) return  100.831;
  if (ibin ==17) return  90.7831;
  if (ibin ==18) return  80.9823;
  if (ibin ==19) return  72.6236;
  if (ibin ==20) return  64.1508;
  if (ibin ==21) return  56.6284;
  if (ibin ==22) return  49.9984;
  if (ibin ==23) return  43.3034;
  if (ibin ==24) return  37.8437;
  if (ibin ==25) return  32.6659;
  if (ibin ==26) return  27.83;
  if (ibin ==27) return  23.7892;
  if (ibin ==28) return  20.1745;
  if (ibin ==29) return  16.8453;
  if (ibin ==30) return  14.0322;
  if (ibin ==31) return  11.602;
  if (ibin ==32) return  9.52528;
  if (ibin ==33) return  7.6984;
  if (ibin ==34) return  6.446;
  if (ibin ==35) return  4.96683;
  if (ibin ==36) return  4.23649;
  if (ibin ==37) return  3.50147;
  if (ibin ==38) return  3.16107;
  if (ibin ==39) return  2.7877;
  return -100000;
}





#endif 
