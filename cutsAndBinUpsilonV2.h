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
int kTrigL1DB50100 = 2;
TString fTrigName[Ntrig] = {"Jpsi", "Ups", "L1DoubleMuOpenOS40100", "L1DoubleMuOpen50100"};

const Int_t nBins = 200; // table of bin edges
const Double_t binTable[nBins+1] = {0, 10.5072, 11.2099, 11.8364, 12.478, 13.1194, 13.7623, 14.4081, 15.0709, 15.7532, 16.4673, 17.1881, 17.923, 18.673, 19.4865, 20.3033, 21.1536, 22.0086, 22.9046, 23.8196, 24.7924, 25.8082, 26.8714, 27.9481, 29.0828, 30.2757, 31.5043, 32.8044, 34.1572, 35.6142, 37.1211, 38.6798, 40.3116, 42.0398, 43.8572, 45.6977, 47.6312, 49.6899, 51.815, 54.028, 56.3037, 58.7091, 61.2024, 63.8353, 66.5926, 69.3617, 72.2068, 75.2459, 78.3873, 81.5916, 84.9419, 88.498, 92.1789, 95.9582, 99.8431, 103.739, 107.78, 111.97, 116.312, 120.806, 125.46, 130.269, 135.247, 140.389, 145.713, 151.212, 156.871, 162.729, 168.762, 174.998, 181.424, 188.063, 194.907, 201.942, 209.19, 216.683, 224.37, 232.291, 240.43, 248.807, 257.416, 266.256, 275.348, 284.668, 294.216, 304.053, 314.142, 324.488, 335.101, 345.974, 357.116, 368.547, 380.283, 392.29, 404.564, 417.122, 429.968, 443.116, 456.577, 470.357, 484.422, 498.78, 513.473, 528.479, 543.813, 559.445, 575.411, 591.724, 608.352, 625.344, 642.686, 660.361, 678.371, 696.749, 715.485, 734.608, 754.068, 773.846, 794.046, 814.649, 835.608, 856.972, 878.719, 900.887, 923.409, 946.374, 969.674, 993.435, 1017.62, 1042.21, 1067.28, 1092.72, 1118.64, 1144.96, 1171.71, 1198.98, 1226.67, 1254.82, 1283.46, 1312.65, 1342.21, 1372.27, 1402.85, 1433.93, 1465.49, 1497.62, 1530.29, 1563.49, 1597.22, 1631.49, 1666.37, 1701.8, 1737.75, 1774.35, 1811.51, 1849.29, 1887.75, 1926.79, 1966.6, 2006.97, 2047.99, 2089.71, 2132.1, 2175.23, 2219.17, 2263.72, 2309.2, 2355.43, 2402.47, 2450.33, 2499.05, 2548.66, 2599.16, 2650.59, 2703.03, 2756.32, 2810.75, 2866.27, 2922.91, 2980.54, 3039.47, 3099.53, 3160.98, 3223.66, 3287.71, 3353.18, 3420.34, 3489.13, 3559.72, 3632.06, 3706.18, 3782.42, 3860.78, 3941.42, 4024.52, 4110.27, 4199.4, 4292.8, 4394.49, 4519.52, 5199.95};

const Double_t binTable_up[nBins+1] = {0, 10.5075, 11.2107, 11.838, 12.4797, 13.1213, 13.7641, 14.4124, 15.0745, 15.7577, 16.473, 17.1939, 17.9297, 18.6812, 19.4958, 20.3143, 21.1648, 22.0218, 22.9159, 23.8328, 24.8059, 25.8204, 26.89, 27.9702, 29.1042, 30.3022, 31.528, 32.8347, 34.1896, 35.6439, 37.1542, 38.7172, 40.3518, 42.091, 43.9053, 45.7415, 47.6853, 49.7457, 51.8755, 54.0983, 56.3594, 58.7848, 61.2861, 63.9228, 66.6825, 69.4421, 72.297, 75.3547, 78.4967, 81.6977, 85.0755, 88.6211, 92.3058, 96.1071, 99.9975, 104.065, 108.272, 112.512, 116.906, 121.601, 126.465, 131.482, 136.866, 142.229, 147.786, 153.546, 159.571, 165.586, 171.902, 178.419, 185.063, 191.856, 199.055, 206.261, 213.999, 221.719, 229.671, 237.84, 246.088, 254.828, 263.883, 272.907, 282.236, 291.925, 301.519, 311.477, 321.691, 332.153, 342.892, 353.878, 365.161, 376.742, 388.577, 400.684, 413.075, 425.746, 438.711, 451.989, 465.556, 479.45, 493.608, 508.077, 522.891, 538.003, 553.415, 569.151, 585.216, 601.601, 618.354, 635.422, 652.84, 670.599, 688.699, 707.161, 726.014, 745.185, 764.687, 784.557, 804.838, 825.489, 846.537, 867.951, 889.752, 911.955, 934.588, 957.52, 980.912, 1004.73, 1028.94, 1053.57, 1078.67, 1104.17, 1130.07, 1156.39, 1183.2, 1210.47, 1238.17, 1266.38, 1295.02, 1324.16, 1353.71, 1383.77, 1414.35, 1445.41, 1477, 1509.09, 1541.74, 1574.88, 1608.59, 1642.83, 1677.66, 1713.07, 1748.98, 1785.47, 1822.63, 1860.33, 1898.72, 1937.73, 1977.42, 2017.71, 2058.62, 2100.25, 2142.57, 2185.56, 2229.38, 2273.91, 2319.2, 2365.33, 2412.22, 2459.94, 2508.52, 2557.98, 2608.35, 2659.61, 2711.86, 2765, 2819.23, 2874.58, 2930.97, 2988.46, 3047.12, 3106.95, 3168.15, 3230.6, 3294.37, 3359.58, 3426.47, 3494.95, 3565.21, 3637.21, 3711.03, 3786.91, 3864.85, 3945.11, 4027.8, 4113.06, 4201.73, 4294.72, 4395.9, 4520.5, 5199.95};

const Double_t binTable_down[nBins+1] =  {0, 10.5071, 11.2094, 11.8357, 12.4763, 13.117, 13.7597, 14.4049, 15.0671, 15.7491, 16.4622, 17.1812, 17.9144, 18.6674, 19.4797, 20.2963, 21.1435, 21.9974, 22.8928, 23.8068, 24.7805, 25.7931, 26.8556, 27.9308, 29.0638, 30.2582, 31.4795, 32.7816, 34.1349, 35.5834, 37.0941, 38.6474, 40.2782, 42.0035, 43.8112, 45.6576, 47.5758, 49.6381, 51.6667, 53.7353, 55.8903, 58.1259, 60.4528, 62.8712, 65.3859, 67.9968, 70.7065, 73.5231, 76.4519, 79.4922, 82.6461, 85.9264, 89.3269, 92.8562, 96.5212, 100.322, 104.262, 108.344, 112.585, 116.971, 121.521, 126.225, 131.09, 136.127, 141.328, 146.721, 152.284, 158.014, 163.935, 170.054, 176.372, 182.878, 189.602, 196.532, 203.653, 211.017, 218.599, 226.387, 234.418, 242.667, 251.16, 259.886, 268.852, 278.071, 287.498, 297.2, 307.184, 317.409, 327.894, 338.66, 349.686, 360.996, 372.607, 384.508, 396.669, 409.133, 421.86, 434.906, 448.258, 461.916, 475.906, 490.16, 504.74, 519.663, 534.911, 550.453, 566.322, 582.525, 599.08, 615.968, 633.211, 650.805, 668.76, 687.048, 705.707, 724.774, 744.163, 763.9, 783.999, 804.528, 825.432, 846.746, 868.429, 890.523, 913.007, 935.952, 959.211, 982.919, 1007.08, 1031.63, 1056.62, 1082.08, 1107.96, 1134.24, 1160.99, 1188.22, 1215.91, 1244.06, 1272.69, 1301.85, 1331.45, 1361.51, 1392.07, 1423.18, 1454.77, 1486.93, 1519.57, 1552.81, 1586.55, 1620.87, 1655.79, 1691.26, 1727.27, 1763.93, 1801.12, 1838.97, 1877.47, 1916.61, 1956.45, 1996.89, 2038.04, 2079.84, 2122.35, 2165.52, 2209.53, 2254.24, 2299.83, 2346.19, 2393.31, 2441.28, 2490.16, 2539.86, 2590.57, 2642.16, 2694.74, 2748.23, 2802.81, 2858.47, 2915.33, 2973.2, 3032.28, 3092.56, 3154.24, 3217.19, 3281.45, 3347.18, 3414.6, 3483.65, 3554.56, 3627.2, 3701.66, 3778.25, 3856.97, 3937.98, 4021.48, 4107.62, 4197.21, 4291.05, 4393.19, 4518.6, 5199.95};


Int_t getHiBinFromhiHF(const Double_t hiHF)
{
  Int_t binPos = -1;
  for(int i = 0; i < nBins; ++i){
    if(hiHF >= binTable[i] && hiHF < binTable[i+1]){
      binPos = i;
      break;
    }
  }

  binPos = nBins - 1 - binPos;

  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins));
}

Int_t getHiBinFromhiHF_Up(const Double_t hiHF)
{
  Int_t binPos = -1;
  for(int i = 0; i < nBins; ++i){
    if(hiHF >= binTable_up[i] && hiHF < binTable_up[i+1]){
      binPos = i;
      break;
    }
  }

  binPos = nBins - 1 - binPos;

  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins));
}

Int_t getHiBinFromhiHF_Down(const Double_t hiHF)
{
  Int_t binPos = -1;
  for(int i = 0; i < nBins; ++i){
    if(hiHF >= binTable_down[i] && hiHF < binTable_down[i+1]){
      binPos = i;
      break;
    }
  }

  binPos = nBins - 1 - binPos;
  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBins));
}


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
