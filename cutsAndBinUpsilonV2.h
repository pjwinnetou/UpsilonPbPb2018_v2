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

int kMuGlb = 1;
int kMuGlbTrk = 2;
  
int selAll = 0;
int selSig = 1;
int selBkg = 2;
int selGlb = 3;
int selGlbTrk = 4;
int selGlbNTrk = 5;
int selGlbSig = 6;
int selGlbBkg = 7;
int selGlbTrkSig = 8;
int selGlbTrkBkg = 9;
int selGlbNTrkSig = 10;
int selGlbNTrkBkg = 11;
int selGlbOTrk = 12;
int selGlbSoftID = 13;
int selGlbTrkSoftID = 14;
int selGlbNTrkSoftID = 15;
int selGlbOTrkSoftID = 16;

const int Ntrig = 4;
int kTrigJpsi = 12;
int kTrigUps = 13;
int kTrigL1DBOS40100 = 0;
int kTrigL1DB50100 = 1;
TString fTrigName[Ntrig] = {"Jpsi", "Ups", "L1DoubleMuOpenOS40100", "L1DoubleMuOpen50100"};


Double_t findNcoll(int hiBin) {
  const int nbins = 200;
  const Double_t Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  return Ncoll[hiBin];
};

bool IsAcceptanceQQ(double pt, double eta)
{
  return ( (fabs(eta)<1.2 && pt>=3.5) ||
      (1.2<=fabs(eta) && fabs(eta)<2.1 && pt>=5.47-1.89*fabs(eta)) ||
      (2.1<=fabs(eta) && fabs(eta)<2.4 && pt>=1.5) );
}

bool IsAcceptanceNoTrig(double pt, double eta)
{
  return ( (fabs(eta)<0.3 && pt>=3.4) ||
      (0.3<=fabs(eta) && fabs(eta)<1.1 && pt>=3.3) ||
      (1.1<=fabs(eta) && fabs(eta)<1.4 && pt>=7.7-4.0*fabs(eta)) ||
      (1.4<=fabs(eta) && fabs(eta)<1.55 && pt>=2.1) ||
      (1.55<=fabs(eta) && fabs(eta)<2.2 && pt>=4.25-1.39*fabs(eta)) ||
      (2.2<=fabs(eta) && fabs(eta)<2.4 && pt>=1.2) );
}



const double inel_cross_PbPb = 6740;
const double NumberOfMBColl =  2454000000;
const double NumberOfMBColl1 = 3092000000;

// lumi Unc 
double lumi_unc_pp = 0.023;
double nMB_unc = TMath::Sqrt(0.02*0.02+0.01*0.01);

struct ParticleMass { double JPsi, Psi2S, Y1S, Y2S, Y3S, Z, PiPlus, KaPlus; };
ParticleMass pdgMass = {3.096, 3.686, 9.460, 10.023, 10.355, 91.188, 0.139570, 0.49367 };

int fpol1 = 1;
int fpol2 = 2;
int fpol3 = 3;


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

TString getKineLabelJpsi(float ptLow, float ptHigh, float yLow, float yHigh, int cLow, int cHigh) {
  TString kineLabel = Form("pt%.1f-%.1f_y%.1f-%.1f",ptLow,ptHigh, yLow, yHigh) ;
    kineLabel = kineLabel+ Form("_centrality%d-%d",(int)cLow, (int)cHigh) ;
  return kineLabel;
}
TString getKineLabel(float ptLow, float ptHigh, float yLow, float yHigh, float muPtCut_, int cLow, int cHigh) {
  TString kineLabel = Form("pt%.1f-%.1f_y%.1f-%.1f_muPt%.1f",ptLow,ptHigh, yLow, yHigh, (float)muPtCut_) ;
    kineLabel = kineLabel+ Form("_centrality%d-%d",(int)cLow, (int)cHigh) ;
  return kineLabel;
}

#endif
