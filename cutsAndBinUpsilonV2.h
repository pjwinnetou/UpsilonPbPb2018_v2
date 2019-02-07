#ifndef CutAndBinCollection_C
#define CutAndBinCollection_C

#include <TF1.h>
#include <TCut.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <iostream>
#include <TLine.h>
#include <TMath.h>
#include <TTree.h>
#include <math.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>


float glbMuPtCut = 4; // for acceptance
const int nPtBins1s  = 6;   double ptBin1s[nPtBins1s+1] = {0,2,4,6,9,12,30};
//const int nPtBins1sMC  = 6;   double ptBin1sMC[nPtBins1sMC+1] = {0,2,4,6,9,12,30};
//const int nPtBins1s  = 3;   double ptBin1s[nPtBins1s+1] = {0,4,9,12,30};
//const int nPtBins1sMC  = 3;   double ptBin1sMC[nPtBins1sMC+1] = {0,4,9,12,30};
const int nPtBins1sMC  = 60;  double ptBin1sMC[nPtBins1sMC+1] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30};  

const int nPtBins2s  = 3;   double ptBin2s[nPtBins2s+1] = {0,4,9,30};
//const int nPtBins2s  = 4;   double ptBin2s[nPtBins2s+1] = {0,4,9,12,30};
//const int nPtBins2sMC  = 3;   double ptBin2sMC[nPtBins2sMC+1] = {0,4,9,30};
//const int nPtBins2sMC  = 4;   double ptBin2sMC[nPtBins2sMC+1] = {0,4,9,12,30};
//const int nPtBins2sMC  = 6;   double ptBin2sMC[nPtBins2sMC+1] = {0,2,4,6,9,12,30};
const int nPtBins2sMC  = 60;  double ptBin2sMC[nPtBins2sMC+1] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30};  

const int nPtBins3s  = 2;   double ptBin3s[nPtBins3s+1] = {0,6,30};
//const int nPtBins3s  = 2;   double ptBin3s[nPtBins3s+1] = {0,15,30}; 

const int nYBins1S  = 6;   double yBin1S[nYBins1S+1] ={0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
const int nYBins2S  = 3;   double yBin2S[nYBins2S+1] ={0, 0.8, 1.6, 2.4};
const int nYBins3S  = 2;   double yBin3S[nYBins3S+1] ={0, 1.2, 2.4};

const int nYBins  = 2;   double yBin[nYBins+1] ={0, 1.2, 2.4}; // for event reweighting

const int nCentBins1s  = 9;   double centBin1s[nCentBins1s+1] = {0,10,20,40,60,80,100,120,140,200}; 
const int nCentBins2s  = 9;   double centBin2s[nCentBins2s+1] = {0,10,20,40,60,80,100,120,140,200};
const int nCentBins3s  = 2;   double centBin3s[nCentBins3s+1] = {0,60,200};
//const int nCentBins3s  = 2;   double centBin3s[nCentBins3s+1] = {0,60,200}; 

// Glauber variables https://twiki.cern.ch/twiki/pub/CMS/HiCentrality2016/AN-15-080_temp_20161206.pdf

double nPart1s[nCentBins1s]   = {8.3, 30.6,53.9, 87.0, 131.4, 189.2, 264.2, 333.3, 384.3}; // HIN-16-008 paper
//double nPart1s[nCentBins1s]   = {15.47,30.59,53.85,86.95,131.4,189.2,264.3,333.4,384.4}; // HIN-16-008 paper
double nPart2s[nCentBins1s]   = {8.3, 30.6,53.9, 87.0, 131.4, 189.2, 264.2, 333.3, 384.3}; // HIN-16-008 paper
double nPart3s[nCentBins3s]   = {46.81, 270.7};
//double nPart3s[nCentBins3s]   = {46.81, 270.7};
double nColl1s[nCentBins1s]   = {1819,1432,1005,606,349,186,90.7,40.1,7.67}; 
double nColl2s[nCentBins2s]   = {1819,1432,1005,606,349,186,90.7,40.1,7.67}; 
double nColl3s[nCentBins3s]   = {1079, 98.36};
//double nColl3s[nCentBins3s]   = {1079, 98.36};  

/*
##Weighting function##
For 1S
TF1* fWgtPP1 = new TF1("fWgtPP1","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )", 0, 30);
Fot 2S, 3S
TF1* fWgtPP2 = new TF1("fWgtPP2","( [0] + [1]*x )", 0, 30);

##Parameters##
  ##Nominal##
fWgtPP1->SetParameters( 200.759, 7.09569, 25.3727, -4.51979 );
fWgtPP2->SetParameters( 0.569212, 0.0637386 );
fWgtAA1->SetParameters( 255.074, -93.4016, 44.2256, -4.81048 );
fWgtAA2->SetParameters( 0.778896, 0.0209981 );
  ##Variation1(Param+error)
fWgtPP1->SetParameters( 254.905, 17.5645, 26.6928, -4.04741 );
fWgtPP2->SetParameters( 0.592215, 0.066485 );
fWgtAA1->SetParameters( 320.738, -77.2557, 48.3587, -4.21744 );
fWgtAA2->SetParameters( 1.08026, 0.0468831 );
  ##Variation2(Param-error)
fWgtPP1->SetParameters( 146.613, -3.37308, 24.0526, -4.99218 );
fWgtPP2->SetParameters( 0.546209, 0.0609922 );
fWgtAA1->SetParameters( 189.409, -109.547, 40.0926, -5.40352 );
fWgtAA2->SetParameters( 0.477534, -0.00488679 );
*/

//Upperlimit value
//double lower68_pt1  = 0.00E+00   ;
//double lower68_pt2  = 0.00E+00   ;
//double lower68_y1   = 0.00E+00   ;
//double lower68_y2   = 0.00E+00   ;
//double lower68_c2   = 0.00E+00   ;
//double lower68_c1   = 0.03618991 ;
//double lower68_cint = 0.00E+00   ;
//
//double upper68_pt1  = 0.081095397  ;
//double upper68_pt2  = 0.051848485  ;
//double upper68_y1   = 0.043153559  ;
//double upper68_y2   = 0.085515537  ;
//double upper68_c2   = 0.037200857  ;
//double upper68_c1   = 0.173370454  ;
//double upper68_cint = 0.04200117   ;
//
//double lower95_pt1 = 0;
//double lower95_pt2 = 0;
//double lower95_y1 = 0;
//double lower95_y2 = 0;
//double lower95_c2 = 0;
//double lower95_c1 = 0;
//double lower95_cint = 0;
//
//double upper95_pt1 =  0.148847232 ;
//double upper95_pt2 =  0.095778193 ;
//double upper95_y1 =   0.08403307  ;
//double upper95_y2 =   0.121629254 ;
//double upper95_c2 =   0.075640382 ;
//double upper95_c1 =   0.242656483 ;
//double upper95_cint = 0.075095649 ;

//Upperlimit value 2018_0219
//double lower68_pt1  = 0.00E+00     ;
//double lower68_pt2  = 0.00E+00     ;
//double lower68_y1   = 0.00E+00     ;
//double lower68_y2   = 0.00E+00     ;
//double lower68_c2   = 0.00E+00     ;
//double lower68_c1   = 0.036170831  ;
//double lower68_cint = 0.00E+00     ;
//
//double upper68_pt1  = 0.075375655  ;
//double upper68_pt2  = 0.051848485  ;
//double upper68_y1   = 0.043153559  ;
//double upper68_y2   = 0.067868394  ;
//double upper68_c2   = 0.037200857  ;
//double upper68_c1   = 0.173381445  ;
//double upper68_cint = 0.048198279  ;
//
//double lower95_pt1 = 0;
//double lower95_pt2 = 0;
//double lower95_y1 = 0;
//double lower95_y2 = 0;
//double lower95_c2 = 0;
//double lower95_c1 = 0;
//double lower95_cint = 0;
//
//double upper95_pt1 =  0.137622193  ;
//double upper95_pt2 =  0.095778193  ;
//double upper95_y1 =   0.08403307   ;
//double upper95_y2 =   0.12882883   ;
//double upper95_c2 =   0.075640382  ;
//double upper95_c1 =   0.242656483  ;
//double upper95_cint = 0.090624869  ;

//Upperlimit value 2018_0420
double lower68_pt1  = 0.00E+00     ;
double lower68_pt2  = 0.00E+00     ;
double lower68_y1   = 0.00E+00     ;
double lower68_y2   = 0.00E+00     ;
double lower68_c2   = 0.00E+00     ;
double lower68_c1   = 0.048408272   ;
double lower68_cint = 0.00E+00     ;

double upper68_pt1  = 0.091511678 ;
double upper68_pt2  = 0.051286058 ;
double upper68_y1   = 0.057943188 ;
double upper68_y2   = 0.06871151  ;
double upper68_c2   = 0.038991035 ;
double upper68_c1   = 0.18472074  ;
double upper68_cint = 0.058262811 ;

double lower95_pt1 = 0;
double lower95_pt2 = 0;
double lower95_y1 = 0;
double lower95_y2 = 0;
double lower95_c2 = 0;
double lower95_c1 = 0;
double lower95_cint = 0;

double upper95_pt1 =  0.155279103;
double upper95_pt2 =  0.094318138;
double upper95_y1 =   0.099952628;
double upper95_y2 =   0.125994869;
double upper95_c2 =   0.076087933;
double upper95_c1 =   0.253908129;
double upper95_cint = 0.094618599;
//Upperlimit CrossSection value
//double lower68XS_pt1 = 0;
//double lower68XS_pt2 = 0;
//double lower68XS_y1 = 0;
//double lower68XS_y2 = 0;
//
//double upper68XS_pt1 = 0.001083787  ;
//double upper68XS_pt2 = 0.000233267  ;
//double upper68XS_y1  = 0.008072324  ;
//double upper68XS_y2  = 0.009555016  ;
//
//double lower95XS_pt1 = 0;
//double lower95XS_pt2 = 0;
//double lower95XS_y1 = 0;
//double lower95XS_y2 = 0;
//
//double upper95XS_pt1 = 0.001989228   ;
//double upper95XS_pt2 = 0.000437013   ;
//double upper95XS_y1 =  0.015721179   ;
//double upper95XS_y2 =  0.018107819   ;

//Upperlimit CrossSection value 2018_0219
//double lower68XS_pt1 = 0;
//double lower68XS_pt2 = 0;
//double lower68XS_y1 = 0;
//double lower68XS_y2 = 0;
//
//double upper68XS_pt1 = 0.0010069    ;
//double upper68XS_pt2 = 0.000235107  ;
//double upper68XS_y1  = 0.009949407  ;
//double upper68XS_y2  = 0.010055214    ;
//
//double lower95XS_pt1 = 0;
//double lower95XS_pt2 = 0;
//double lower95XS_y1 = 0;
//double lower95XS_y2 = 0;
//
//double upper95XS_pt1 = 0.001838352  ;
//double upper95XS_pt2 = 0.000440664  ;
//double upper95XS_y1 =  0.019399529  ;
//double upper95XS_y2 =  0.01917852    ;

//Upperlimit CrossSection value 2018_0219
//double lower68XS_pt1 = 0;
//double lower68XS_pt2 = 0;
//double lower68XS_y1 = 0;
//double lower68XS_y2 = 0;
//
//double upper68XS_pt1 = 0.001347207 ;
//double upper68XS_pt2 = 0.000181412 ;
//double upper68XS_y1  = 0.010154671 ;
//double upper68XS_y2  = 0.00991101  ;
//
//double lower95XS_pt1 = 0;
//double lower95XS_pt2 = 0;
//double lower95XS_y1 = 0;
//double lower95XS_y2 = 0;
//
//double upper95XS_pt1 = 0.001961598 ;
//double upper95XS_pt2 = 0.000336498 ;
//double upper95XS_y1 =  0.017516703 ;
//double upper95XS_y2 =  0.018173723 ;

//Upperlimit CrossSection value 2018_0507
double lower68XS_pt1 = 0;
double lower68XS_pt2 = 0;
double lower68XS_y1 = 0;
double lower68XS_y2 = 0;

double upper68XS_pt1 = 0.001347207  ;
double upper68XS_pt2 = 0.000181412  ;
double upper68XS_y1  = 0.010154671  ;
double upper68XS_y2  = 0.009911025  ;

double lower95XS_pt1 = 0;
double lower95XS_pt2 = 0;
double lower95XS_y1 = 0;
double lower95XS_y2 = 0;

double upper95XS_pt1 = 0.002285294  ;
double upper95XS_pt2 = 0.000336498  ;
double upper95XS_y1 =  0.017516703  ;
double upper95XS_y2 =  0.018173723  ;
// TAA Value
double TAA1s[nCentBins1s+1] = {25.98, 20.46, 14.35, 8.66, 4.978, 2.66, 1.296, 0.5729, 0.1095, 5.607};
double TAA2s[nCentBins2s+1] = {25.98, 20.46, 14.35, 8.66, 4.978, 2.66, 1.296, 0.5729, 0.1095, 5.607};
double TAA3s[nCentBins3s+1] = {15.41, 1.405, 5.607};

// TAA Unc 
double TAA_unc1s[nCentBins1s+1] = {0.017, 0.017, 0.02, 0.028, 0.04, 0.058, 0.081, 0.11, 0.18, 0.089};
double TAA_unc2s[nCentBins2s+1] = {0.017, 0.017, 0.02, 0.028, 0.04, 0.058, 0.081, 0.11, 0.18, 0.089};
double TAA_unc3s[nCentBins3s+1] = {0.022, 0.12, 0.089};

double TAA_unc1sLo[nCentBins1s+1] = {0.029, 0.03, 0.032, 0.038, 0.049, 0.066, 0.088, 0.112, 0.102, 0.034};
double TAA_unc2sLo[nCentBins2s+1] = {0.029, 0.03, 0.032, 0.038, 0.049, 0.066, 0.088, 0.112, 0.102, 0.034};
double TAA_unc3sLo[nCentBins2s+1] = {0.031, 0.044, 0.034};

double TAA_unc1sHi[nCentBins1s+1] = {0.018, 0.019, 0.023, 0.033, 0.047, 0.068, 0.092, 0.124, 0.164, 0.028};
double TAA_unc2sHi[nCentBins2s+1] = {0.018, 0.019, 0.023, 0.033, 0.047, 0.068, 0.092, 0.124, 0.164, 0.028};
double TAA_unc3sHi[nCentBins2s+1] = {0.021, 0.067, 0.028};

const double inel_cross_PbPb = 6740;
//const double NumberOfMBColl = 2366003000;
const double NumberOfMBColl =  2454000000;
//const double NumberOfMBColl = 2484303150;
const double NumberOfMBColl1 = 3092000000;
//const double NumberOfMBColl1 = 3284093053;
//const double inel_cross_PbPb = 7716;

// lumi Unc 
double lumi_unc_pp = 0.023;
double nMB_unc = TMath::Sqrt(0.02*0.02+0.01*0.01);

struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z, PiPlus, KaPlus; };
ParticleMass pdgMass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188, 0.139570, 0.49367 };

struct valErr { float val, err; } ; 

int kPPDATA = 0 ;
int kPADATA = 1 ;
int kAADATA = 2 ; // L1 doubleMu 0
int kPPMC = 3 ;
int kPAMC = 4 ;
int kAAMC = 5 ;
int kAADATAPeri = 6 ;
int kAADATACentL3 = 7 ;
int kPPMCUps1S = 8 ;
int kPPMCUps2S = 9 ;
int kPPMCUps3S = 10 ;
int kAAMCUps1S = 11 ;
int kAAMCUps2S = 12 ;
int kAAMCUps3S = 13 ;
int kPPAADATASIMUL = 20 ; // 2 and 0 simultaneous fit
int kPPAADATAPeriSIMUL = 60 ; // 6 and 0 simultaneous fit

TString getCollID( int collid ) {
  if ( collid == kPPDATA ) return "PP_DATA";
  else if ( collid == kPADATA ) return "PA_DATA";
  else if ( collid == kAADATA ) return "AA_DATA";
  else if ( collid == kPPMC ) return "PP_MC";
  else if ( collid == kPAMC ) return "PA_MC";
  else if ( collid == kAAMC ) return "AA_MC";
//  else if ( collid == kAADATAPeri ) return "AA_DATA_PeriL1";
  else if ( collid == kAADATAPeri ) return "AA_DATA";
  else if ( collid == kAADATACentL3 ) return "AA_DATA_CentL3";
  else if ( collid == kPPMCUps1S ) return "PP_MC_Ups1S";
  else if ( collid == kPPMCUps2S ) return "PP_MC_Ups2S";
  else if ( collid == kPPMCUps3S ) return "PP_MC_Ups3S";
  else if ( collid == kAAMCUps1S ) return "AA_MC_Ups1S";
  else if ( collid == kAAMCUps2S ) return "AA_MC_Ups2S";
  else if ( collid == kAAMCUps3S ) return "AA_MC_Ups3S";
  else if ( collid == kPPAADATASIMUL ) return "PP_AA_DATA_SIMUL";
  else if ( collid == kPPAADATAPeriSIMUL ) return "PP_AA_DATA_PeriL1_SIMUL";

  else return "none";
}

int kEPl2HF = 0;
int kEPOppositeHF = 1;
int kEPSameSideHF = 2;


TString getEPSel( int eventPln) {
  if ( eventPln == kEPl2HF)  return "BothHFs";
  else if ( eventPln == kEPOppositeHF ) return "OppositeHF" ;
  else if ( eventPln == kEPSameSideHF ) return "SameSideHF" ;
  else return "none";
}


int kSoftMuCut = 0;
int kHighPtMuCut = 0;




class DiMuon {
 public:
 DiMuon() :
  run(0),   lumi(0), event(0), cBin(0), ep2(0), dphiEp2(0),
    vz(-99),  mass(-1), pt(-1), y(999), phi(999), eta(999),
    pt1(-1), eta1(-1), phi1(-1),        
    pt2(-1), eta2(-1), phi2(-1), weight0(0), weight(0),       
    oniaIndex(-1), softFlag(0), highPtFlag(0),
    qxa(0), qya(0),
    qxb(0), qyb(0),
    qxc(0), qyc(0),
    qxdimu(0), qydimu(0)

    {}
  
  int run;
  int lumi;
  int event;
  int cBin;
  float ep2;
  float dphiEp2;
  float vz;
  float mass;
  float pt;
  float y;
  float phi;    
  float eta;
  float pt1; 
  float eta1;
  float phi1;
  float pt2;
  float eta2;
  float phi2;    
  float weight0;
  float weight;
  int oniaIndex;
  int softFlag;
  int highPtFlag;
  float qxa;
  float qya;
  float qxb;
  float qyb;
  float qxc;
  float qyc;
  float qxdimu;
  float qydimu;

  void clear() {
    run = -99;  lumi=-99; event=-99; cBin=-99; ep2=-99, dphiEp2=-99; 
    vz=-99;     mass = -99; pt=-99; y=-99; phi=-99; eta=-99;      
    pt1=-99; eta1=-99; phi1=-99; pt2=-99; eta2=-99; phi2=-99; weight0=-99, weight=-99;
    oniaIndex=-1; softFlag=-1; highPtFlag=-1; 
    qxa=0;  qya=0; 
    qxb=0;  qyb=0; 
    qxc=0;  qyc=0; 
    qxdimu=0; qydimu=0;
  }

};
TString branchString = "run/I:lumi:event:cBin:ep2/F:dphiEp2:vz:mass:pt:y:phi:eta:pt1:eta1:phi1:pt2:eta2:phi2:weight0:weight:oniaIndex/I:softFlag:highPtFlag:qxa/F:qya:qxb:qyb:qxc:qyc:qxdimu:qydimu";



// Upsilon nominal bins
const int nPtBinsUps = 2;   double ptBinUps[nPtBinsUps+1] = {0, 5,     100};
const int  nYBinsUps = 2;   double yBinUps[nYBinsUps+1] =   {0, 1.2,   2.4};
const int nPBinsUps  = 3;   double pBinUps[nPBinsUps+1] =   {0, 0.167, 0.333,  0.5};


TString getKineLabel(int collId, float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut_, int cLow, int cHigh, float dphiEp2Low, float dphiEp2High) {
  TString kineLabel = Form("%s_pt%.1f-%.1f_y%.1f-%.1f_muPt%.1f",getCollID(collId).Data(), ptLow,ptHigh, yLow, yHigh, (float)muPtCut_) ;
  if ( (collId == kAADATA) || (collId == kPADATA) || (collId == kAAMC) || (collId == kPAMC) || (collId == kAADATAPeri ) || ( collId == kAADATACentL3) || (collId == kAAMCUps1S) || ( collId==kAAMCUps2S) || (collId == kAAMCUps3S) || (collId == kPPAADATASIMUL) || (collId == kPPAADATAPeriSIMUL))
    kineLabel = kineLabel+ Form("_centrality%d-%d_dphiEp_%.2fPI_%.2fPI",(int)cLow, (int)cHigh, (float)dphiEp2Low, (float)dphiEp2High ) ;
  return kineLabel;
}

#endif
