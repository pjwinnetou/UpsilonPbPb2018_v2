#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../Style_jaebeom.h"
using namespace std;


//// do NOT use "hadded" ttrees!! (e.g.6-90 GeV) 

TLegend *leg = new TLegend(0.55,0.2, 0.85,0.4,NULL,"brNDC");
void mergeEightInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D*h5=0, TH1D*h6=0, TH1D*h7=0, TH1D*h8=0, const char* str = " ");
void mergeSevenInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D*h5=0, TH1D*h6=0, TH1D*h7=0);
void mergeSixInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D*h5=0, TH1D*h6=0);
void mergeFiveInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0, TH1D*h5=0);
void mergeFourInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0, TH1D* h4=0);
void mergeThreeInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0, TH1D* h3=0);
void mergeTwoInQuad( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void mergeTwoInQuadCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);
void subtractTwo( TH1D* h0=0, TH1D* h1=0, TH1D* h2=0);
void subtractTwoCent( TH1D* h0=0, TH1D* hAA=0, TH1D* hPP=0);

void mergeSys() { 
  
  TH1::SetDefaultSumw2();

 
  TH1D* hpt_Cent1090[10];
  TH1D* hpt_Cent010[10];
  TH1D* hpt_Cent1050[10];
  TH1D* hpt_Cent5090[10];

  TH1D* hCent_Pt050[10];
  TH1D* hCent_Pt03[10];
  TH1D* hCent_Pt36[10];
  TH1D* hCent_Pt610[10];
  TH1D* hCent_Pt1050[10];
  
  TH1D* hInt[10];

  TH1D* hpt_Cent560_y024[10];
  TH1D* hpt_Cent560_y012[10];
  TH1D* hpt_Cent560_y1224[10];


  // 1 : efficiency
  TFile* f1 = new TFile("Efficiency/Efficiency_sys.root");
  hpt_Cent1090[1]       = (TH1D*) f1->Get("hPtSys_int");
  hpt_Cent010[1]        = (TH1D*) f1->Get("hPtSys_cent010");
  hpt_Cent1050[1]       = (TH1D*) f1->Get("hPtSys_cent1050");
  hpt_Cent5090[1]       = (TH1D*) f1->Get("hPtSys_cent5090");

  hCent_Pt050[1]        = (TH1D*) f1->Get("hCentSys_int");
  hCent_Pt03[1]         = (TH1D*) f1->Get("hCentSys_pt03");
  hCent_Pt36[1]         = (TH1D*) f1->Get("hCentSys_pt36");
  hCent_Pt610[1]        = (TH1D*) f1->Get("hCentSys_pt610");
  hCent_Pt1050[1]        = (TH1D*) f1->Get("hCentSys_pt1050");
  hInt[1]               = (TH1D*) f1->Get("hSys_int");
  
  hpt_Cent560_y024[1]   = (TH1D*) f1->Get("hPtSys_rap024");
  hpt_Cent560_y012[1]   = (TH1D*) f1->Get("hPtSys_rap012");
  hpt_Cent560_y1224[1]  = (TH1D*) f1->Get("hPtSys_rap1224");


  // 2 : acceptance
  TFile* f2 = new TFile("Acceptance/Acceptance_sys.root");
  hpt_Cent1090[2]   = (TH1D*) f2->Get("hPtSys_int");
  hpt_Cent010[2]    = (TH1D*) f2->Get("hPtSys_cent010");
  hpt_Cent1050[2]   = (TH1D*) f2->Get("hPtSys_cent1050");
  hpt_Cent5090[2]  = (TH1D*) f2->Get("hPtSys_cent5090");

  hCent_Pt050[2]   = (TH1D*) f2->Get("hCentSys_int");
  hCent_Pt03[2]    = (TH1D*) f2->Get("hCentSys_pt03");
  hCent_Pt36[2]    = (TH1D*) f2->Get("hCentSys_pt36");
  hCent_Pt610[2]        = (TH1D*) f2->Get("hCentSys_pt610");
  hCent_Pt1050[2]        = (TH1D*) f2->Get("hCentSys_pt1050");
  hInt[2]   = (TH1D*) f2->Get("hSys_int");
  
  hpt_Cent560_y024[2]   = (TH1D*) f2->Get("hPtSys_rap024");
  hpt_Cent560_y012[2]   = (TH1D*) f2->Get("hPtSys_rap012");
  hpt_Cent560_y1224[2]  = (TH1D*) f2->Get("hPtSys_rap1224");


  // 3 : signal Par
  TFile* f3 = new TFile("SignalParVariation/SignalParVar_sys.root");
  hpt_Cent1090[3]   = (TH1D*) f3->Get("hPtSys_int");
  hpt_Cent010[3]    = (TH1D*) f3->Get("hPtSys_cent010");
  hpt_Cent1050[3]   = (TH1D*) f3->Get("hPtSys_cent1050");
  hpt_Cent5090[3]  = (TH1D*) f3->Get("hPtSys_cent5090");

  hCent_Pt050[3]   = (TH1D*) f3->Get("hCentSys_int");
  hCent_Pt03[3]    = (TH1D*) f3->Get("hCentSys_pt03");
  hCent_Pt36[3]    = (TH1D*) f3->Get("hCentSys_pt36");
  hCent_Pt610[3]        = (TH1D*) f3->Get("hCentSys_pt610");
  hCent_Pt1050[3]        = (TH1D*) f3->Get("hCentSys_pt1050");
  hInt[3]   = (TH1D*) f3->Get("hSys_int");

  hpt_Cent560_y024[3]   = (TH1D*) f3->Get("hPtSys_rap024");
  hpt_Cent560_y012[3]   = (TH1D*) f3->Get("hPtSys_rap012");
  hpt_Cent560_y1224[3]  = (TH1D*) f3->Get("hPtSys_rap1224");


  // 4 : signal PDF var
  TFile* f4 = new TFile("SignalPDFVariation/AllParFreeFit/SignalPDF_sys.root");
  hpt_Cent1090[4]   = (TH1D*) f4->Get("hPtSys_int");
  hpt_Cent010[4]    = (TH1D*) f4->Get("hPtSys_cent010");
  hpt_Cent1050[4]   = (TH1D*) f4->Get("hPtSys_cent1050");
  hpt_Cent5090[4]  = (TH1D*) f4->Get("hPtSys_cent5090");

  hCent_Pt050[4]   = (TH1D*) f4->Get("hCentSys_int");
  hCent_Pt03[4]    = (TH1D*) f4->Get("hCentSys_pt03");
  hCent_Pt36[4]    = (TH1D*) f4->Get("hCentSys_pt36");
  hCent_Pt610[4]        = (TH1D*) f4->Get("hCentSys_pt610");
  hCent_Pt1050[4]        = (TH1D*) f4->Get("hCentSys_pt1050");
  hInt[4]   = (TH1D*) f4->Get("hSys_int");
  
  hpt_Cent560_y024[4]   = (TH1D*) f4->Get("hPtSys_rap024");
  hpt_Cent560_y012[4]   = (TH1D*) f4->Get("hPtSys_rap012");
  hpt_Cent560_y1224[4]  = (TH1D*) f4->Get("hPtSys_rap1224");


  // 5 : background PDF var
  TFile* f5 = new TFile("BackgroundPDFVariation/AllParmFreeFit/BackgroundPDF_sys.root");
  hpt_Cent1090[5]   = (TH1D*) f5->Get("hPtSys_int");
  hpt_Cent010[5]    = (TH1D*) f5->Get("hPtSys_cent010");
  hpt_Cent1050[5]   = (TH1D*) f5->Get("hPtSys_cent1050");
  hpt_Cent5090[5]  = (TH1D*) f5->Get("hPtSys_cent5090");

  hCent_Pt050[5]   = (TH1D*) f5->Get("hCentSys_int");
  hCent_Pt03[5]    = (TH1D*) f5->Get("hCentSys_pt03");
  hCent_Pt36[5]    = (TH1D*) f5->Get("hCentSys_pt36");
  hCent_Pt610[5]        = (TH1D*) f5->Get("hCentSys_pt610");
  hCent_Pt1050[5]        = (TH1D*) f5->Get("hCentSys_pt1050");
  hInt[5]   = (TH1D*) f5->Get("hSys_int");
  
  hpt_Cent560_y024[5]   = (TH1D*) f5->Get("hPtSys_rap024");
  hpt_Cent560_y012[5]   = (TH1D*) f5->Get("hPtSys_rap012");
  hpt_Cent560_y1224[5]  = (TH1D*) f5->Get("hPtSys_rap1224");



  // 6 : v2 background func. var
  TFile* f6 = new TFile("Background_v2Variation/SigAllFreeFitFix/V2Bkg_sys_max.root");
  hpt_Cent1090[6]   = (TH1D*) f6->Get("hPtSys_int");
  hpt_Cent010[6]    = (TH1D*) f6->Get("hPtSys_cent010");
  hpt_Cent1050[6]   = (TH1D*) f6->Get("hPtSys_cent1050");
  hpt_Cent5090[6]  = (TH1D*) f6->Get("hPtSys_cent5090");

  hCent_Pt050[6]   = (TH1D*) f6->Get("hCentSys_int");
  hCent_Pt03[6]    = (TH1D*) f6->Get("hCentSys_pt03");
  hCent_Pt36[6]    = (TH1D*) f6->Get("hCentSys_pt36");
  hCent_Pt610[6]        = (TH1D*) f6->Get("hCentSys_pt610");
  hCent_Pt1050[6]        = (TH1D*) f6->Get("hCentSys_pt1050");
  hInt[6]   = (TH1D*) f6->Get("hSys_int");
  
  hpt_Cent560_y024[6]   = (TH1D*) f6->Get("hPtSys_rap024");
  hpt_Cent560_y012[6]   = (TH1D*) f5->Get("hPtSys_rap012");  // don't use. Take from bkg pdf var
  hpt_Cent560_y1224[6]  = (TH1D*) f5->Get("hPtSys_rap1224"); // don't use. Take from bkg pdf var


  // 7 : Event Selection
  TFile* f7 = new TFile("EventSelection/EventSelection_sys.root");
  hpt_Cent1090[7]   = (TH1D*) f7->Get("hPtSys_int");
  hpt_Cent010[7]    = (TH1D*) f7->Get("hPtSys_cent010");
  hpt_Cent1050[7]   = (TH1D*) f7->Get("hPtSys_cent1050");
  hpt_Cent5090[7]  = (TH1D*) f7->Get("hPtSys_cent5090");

  hCent_Pt050[7]   = (TH1D*) f7->Get("hCentSys_int");
  hCent_Pt03[7]    = (TH1D*) f7->Get("hCentSys_pt03");
  hCent_Pt36[7]    = (TH1D*) f7->Get("hCentSys_pt36");
  hCent_Pt610[7]        = (TH1D*) f7->Get("hCentSys_pt610");
  hCent_Pt1050[7]        = (TH1D*) f7->Get("hCentSys_pt1050");
  hInt[7]          = (TH1D*) f7->Get("hSys_int");
  
  hpt_Cent560_y024[7]   = (TH1D*) f7->Get("hPtSys_rap024");
  hpt_Cent560_y012[7]   = (TH1D*) f7->Get("hPtSys_rap012");
  hpt_Cent560_y1224[7]  = (TH1D*) f7->Get("hPtSys_rap1224");



  // 8 : TnP Correction
  TFile* f8 = new TFile("TnP/TnP_sys.root");
  hpt_Cent1090[8]   = (TH1D*) f8->Get("hPtSys_int");
  hpt_Cent010[8]    = (TH1D*) f8->Get("hPtSys_cent010");
  hpt_Cent1050[8]   = (TH1D*) f8->Get("hPtSys_cent1050");
  hpt_Cent5090[8]   = (TH1D*) f8->Get("hPtSys_cent5090");

  hCent_Pt050[8]   = (TH1D*) f8->Get("hCentSys_int");
  hCent_Pt03[8]    = (TH1D*) f8->Get("hCentSys_pt03");
  hCent_Pt36[8]    = (TH1D*) f8->Get("hCentSys_pt36");
  hCent_Pt610[8]        = (TH1D*) f8->Get("hCentSys_pt610");
  hCent_Pt1050[8]        = (TH1D*) f8->Get("hCentSys_pt1050");
  hInt[8]          = (TH1D*) f8->Get("hSys_int");
  
  hpt_Cent560_y024[8]   = (TH1D*) f8->Get("hPtSys_rap024");
  hpt_Cent560_y012[8]   = (TH1D*) f8->Get("hPtSys_rap012");
  hpt_Cent560_y1224[8]  = (TH1D*) f8->Get("hPtSys_rap1224");


  hpt_Cent1090[0]  = (TH1D*) hpt_Cent1090[5]  -> Clone("hpt_Cent1090_merged"); hpt_Cent1090[0]->Reset();
  hpt_Cent010[0]   = (TH1D*) hpt_Cent010[5]   -> Clone("hpt_Cent010_merged"); hpt_Cent010[0]->Reset();
  hpt_Cent1050[0]  = (TH1D*) hpt_Cent1050[5]  -> Clone("hpt_Cent1050_merged"); hpt_Cent1050[0]->Reset();
  hpt_Cent5090[0] = (TH1D*) hpt_Cent5090[5] -> Clone("hpt_Cent5090_merged"); hpt_Cent5090[0]->Reset();

  hCent_Pt050[0]  = (TH1D*) hCent_Pt050[5]  -> Clone("hCent_Pt050_merged"); hCent_Pt050[0]->Reset(); 
  hCent_Pt03[0]   = (TH1D*) hCent_Pt03[5]   -> Clone("hCent_Pt03_merged");  hCent_Pt03[0] ->Reset(); 
  hCent_Pt36[0]   = (TH1D*) hCent_Pt36[5]   -> Clone("hCent_Pt36_merged");  hCent_Pt36[0] ->Reset(); 
  hCent_Pt610[0]  = (TH1D*) hCent_Pt610[5]  -> Clone("hCent_Pt610_merged"); hCent_Pt610[0]->Reset(); 
  hCent_Pt1050[0]  = (TH1D*) hCent_Pt1050[5]  -> Clone("hCent_Pt1050_merged"); hCent_Pt1050[0]->Reset(); 
  hInt[0]         = (TH1D*) hInt[5]  -> Clone("hSys_int_merged"); hInt[0]->Reset();
  
  hpt_Cent560_y024[0]   = (TH1D*) hpt_Cent560_y024[5] -> Clone("hpt_Cent560_y024_merged"); hpt_Cent560_y024[0]->Reset();
  hpt_Cent560_y012[0]   = (TH1D*) hpt_Cent560_y012[5] -> Clone("hpt_Cent560_y012_merged"); hpt_Cent560_y012[0]->Reset();
  hpt_Cent560_y1224[0]  = (TH1D*) hpt_Cent560_y1224[5] -> Clone("hpt_Cent560_y1224_merged"); hpt_Cent560_y1224[0]->Reset();


  // Merge uncertainties for cross-section 
/*
  mergeSixInQuad(hpt_Cent090[0], hpt_Cent090[1], hpt_Cent090[2], hpt_Cent090[3], hpt_Cent090[4], hpt_Cent090[5], hpt_Cent090[6]);  
  mergeSixInQuad(hpt_Cent010[0], hpt_Cent010[1], hpt_Cent010[2], hpt_Cent010[3], hpt_Cent010[4], hpt_Cent010[5], hpt_Cent010[6]);  
  mergeSixInQuad(hpt_Cent1030[0], hpt_Cent1030[1], hpt_Cent1030[2], hpt_Cent1030[3], hpt_Cent1030[4], hpt_Cent1030[5], hpt_Cent1030[6]);  
  mergeSixInQuad(hpt_Cent3050[0], hpt_Cent3050[1], hpt_Cent3050[2], hpt_Cent3050[3], hpt_Cent3050[4], hpt_Cent3050[5], hpt_Cent3050[6]);  
  mergeSixInQuad(hpt_Cent5090[0], hpt_Cent5090[1], hpt_Cent5090[2], hpt_Cent5090[3], hpt_Cent5090[4], hpt_Cent5090[5], hpt_Cent5090[6]);  

  mergeSixInQuad(hCent_Pt050[0], hCent_Pt050[1], hCent_Pt050[2], hCent_Pt050[3], hCent_Pt050[4], hCent_Pt050[5], hCent_Pt050[6]);  
  mergeSixInQuad(hCent_Pt03[0], hCent_Pt03[1], hCent_Pt03[2], hCent_Pt03[3], hCent_Pt03[4], hCent_Pt03[5], hCent_Pt03[6]);  
  mergeSixInQuad(hCent_Pt36[0], hCent_Pt36[1], hCent_Pt36[2], hCent_Pt36[3], hCent_Pt36[4], hCent_Pt36[5], hCent_Pt36[6]);  
  mergeSixInQuad(hCent_Pt650[0], hCent_Pt650[1], hCent_Pt650[2], hCent_Pt650[3], hCent_Pt650[4], hCent_Pt650[5], hCent_Pt650[6]);  
  
  mergeSixInQuad(hInt[0], hInt[1], hInt[2], hInt[3], hInt[4], hInt[5], hInt[6]);  
  */
  
  const char* str_cent[3] = {"Cent. 0-10%", "Cent. 10-50%", "Cent. 50-90%"};
  const char* str_pt[4] = {"p_{T} 0-3 GeV", "p_{T} 3-6 GeV", "p_{T} 6-10 GeV", "p_{T} 10-50 GeV"};
  const char* str_centc = "Cent. 5-60%";
  const char* str_cent_int = "Cent. 10-90%";
  const char* str_pt_int = "p_{T} 0-50 GeV";
  const char* str_int = "p_{T} 0-50 GeV, Cent. 10-90%";

  mergeEightInQuad(hpt_Cent1090[0],  hpt_Cent1090[1],  hpt_Cent1090[2],  hpt_Cent1090[3],  hpt_Cent1090[4],  hpt_Cent1090[5],  hpt_Cent1090[6],  hpt_Cent1090[7], hpt_Cent1090[8], str_cent_int);  
  mergeEightInQuad(hpt_Cent010[0] ,  hpt_Cent010[1] ,  hpt_Cent010[2] ,  hpt_Cent010[3] ,  hpt_Cent010[4] ,  hpt_Cent010[5] ,  hpt_Cent010[6] ,  hpt_Cent010[7] , hpt_Cent010[8] , str_cent[0]);  
  mergeEightInQuad(hpt_Cent1050[0],  hpt_Cent1050[1],  hpt_Cent1050[2],  hpt_Cent1050[3],  hpt_Cent1050[4],  hpt_Cent1050[5],  hpt_Cent1050[6],  hpt_Cent1050[7], hpt_Cent1050[8], str_cent[1]);  
  mergeEightInQuad(hpt_Cent5090[0],  hpt_Cent5090[1],  hpt_Cent5090[2],  hpt_Cent5090[3],  hpt_Cent5090[4],  hpt_Cent5090[5],  hpt_Cent5090[6],  hpt_Cent5090[7], hpt_Cent5090[8], str_cent[2]);  
                                                                                                                                                         
  mergeEightInQuad(hCent_Pt050[0] ,  hCent_Pt050[1] ,   hCent_Pt050[2] ,   hCent_Pt050[3] ,   hCent_Pt050[4] ,   hCent_Pt050[5] ,   hCent_Pt050[6] ,   hCent_Pt050[7] ,  hCent_Pt050[8] , str_pt_int);  
  mergeEightInQuad(hCent_Pt03[0]  ,  hCent_Pt03[1]  ,   hCent_Pt03[2]  ,   hCent_Pt03[3]  ,   hCent_Pt03[4]  ,   hCent_Pt03[5]  ,   hCent_Pt03[6]  ,   hCent_Pt03[7]  ,  hCent_Pt03[8]  , str_pt[0]);  
  mergeEightInQuad(hCent_Pt36[0]  ,  hCent_Pt36[1]  ,   hCent_Pt36[2]  ,   hCent_Pt36[3]  ,   hCent_Pt36[4]  ,   hCent_Pt36[5]  ,   hCent_Pt36[6]  ,   hCent_Pt36[7]  ,  hCent_Pt36[8]  , str_pt[1]);  
  mergeEightInQuad(hCent_Pt610[0] ,  hCent_Pt610[1] ,   hCent_Pt610[2] ,   hCent_Pt610[3] ,   hCent_Pt610[4] ,   hCent_Pt610[5] ,   hCent_Pt610[6] ,   hCent_Pt610[7] ,  hCent_Pt610[8] , str_pt[2]);  
  mergeEightInQuad(hCent_Pt1050[0],  hCent_Pt1050[1],   hCent_Pt1050[2],   hCent_Pt1050[3],   hCent_Pt1050[4],   hCent_Pt1050[5],   hCent_Pt1050[6],   hCent_Pt1050[7],  hCent_Pt1050[8], str_pt[3]);  
  
  mergeEightInQuad(hInt[0], hInt[1], hInt[2], hInt[3], hInt[4], hInt[5], hInt[6], hInt[7], hInt[8], str_int);  
  
  mergeEightInQuad(hpt_Cent560_y024[0], hpt_Cent560_y024[1], hpt_Cent560_y024[2], hpt_Cent560_y024[3], hpt_Cent560_y024[4], hpt_Cent560_y024[5], hpt_Cent560_y024[6], hpt_Cent560_y024[7], hpt_Cent560_y024[8], str_centc);  
//  mergeEightInQuad(hpt_Cent560_y012[0], hpt_Cent560_y012[1], hpt_Cent560_y012[2], hpt_Cent560_y012[3], hpt_Cent560_y012[4], hpt_Cent560_y012[5], hpt_Cent560_y012[6], hpt_Cent560_y012[7], hpt_Cent560_y012[8]);  
//  mergeEightInQuad(hpt_Cent560_y1224[0], hpt_Cent560_y1224[1], hpt_Cent560_y1224[2], hpt_Cent560_y1224[3], hpt_Cent560_y1224[4], hpt_Cent560_y1224[5], hpt_Cent560_y1224[6], hpt_Cent560_y1224[7], hpt_Cent560_y1224[8]);  


  TFile* fout = new TFile("merged_sys.root","recreate");
  hpt_Cent1090[0]->Write();
  hpt_Cent010[0]->Write();
  hpt_Cent1050[0]->Write();
  hpt_Cent5090[0]->Write();

  hCent_Pt050[0]->Write();
  hCent_Pt03[0]->Write();
  hCent_Pt36[0]->Write();
  hCent_Pt610[0]->Write();
  hCent_Pt1050[0]->Write();
  hInt[0]->Write();

  hpt_Cent560_y024[0]->Write();
  hpt_Cent560_y012[0]->Write();
  hpt_Cent560_y1224[0]->Write();


  fout->Close();

}


void mergeEightInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7, TH1D* h8, const char* str) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a6 = h6->GetBinContent(i);
    float a7 = h7->GetBinContent(i);
    float a8 = h8->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5 + a6*a6 + a7*a7 + a8*a8);
    h0->SetBinContent( i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  h0->SetAxisRange(0,0.033,"Y");
  h0->SetYTitle("Difference Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3,        4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4,        5); h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h5,        6); h5->SetLineWidth(2); h5->DrawCopy("hist same");
  handsomeTH1(h6,        7); h6->SetLineWidth(2); h6->DrawCopy("hist same");
  handsomeTH1(h7,        14); h7->SetLineWidth(2); h7->DrawCopy("hist same");
  handsomeTH1(h8,        20); h8->SetLineWidth(2); h8->DrawCopy("hist same");
  

  TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
  easyLeg(leg1,Form("%s, #varUpsilon(1S)",h0->GetName()));
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"efficiency","l");
  leg1->AddEntry(h2,"Acceptance","l");
  leg1->AddEntry(h3,"Signal Parameter","l");
  leg1->AddEntry(h4,"Background PDF","l");
  leg1->AddEntry(h5,"Signal PDF","l");
  leg1->AddEntry(h6,"v_{2} background PDF","l");
  leg1->AddEntry(h7,"Event Selection","l");
  leg1->AddEntry(h8,"Tag-And-Probe","l");
  leg1->Draw();
  
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.035);
  globtex->DrawLatex(0.15, 0.87, Form("%s",str));
  c0->SaveAs(Form("pdfFiles/Unc_%s.pdf", h0->GetName() ));
  


  // 7 : Event Selection & Contamination
  // 6 : v2 bkg pdf
  // 5 : CB+Gaus PDF  
  // 4 : background PDF
  // 3 : signal Par
  // 2 : acceptance
  // 1 : efficiency

  
}

void mergeSevenInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6, TH1D* h7) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a6 = h6->GetBinContent(i);
    float a7 = h7->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5 + a6*a6 + a7*a7);
    h0->SetBinContent( i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  h0->SetAxisRange(0,0.033,"Y");
  h0->SetYTitle("Difference Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3,        4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4,        5); h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h5,        6); h5->SetLineWidth(2); h5->DrawCopy("hist same");
  handsomeTH1(h6,        7); h6->SetLineWidth(2); h6->DrawCopy("hist same");
  handsomeTH1(h7,        13); h7->SetLineWidth(2); h7->DrawCopy("hist same");
  

  TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
  easyLeg(leg1,h0->GetName());
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"Tag-And-Probe","l");
  //leg1->AddEntry(h1,"efficiency","l");
  leg1->AddEntry(h2,"Acceptance","l");
  leg1->AddEntry(h3,"Signal Parameter","l");
  leg1->AddEntry(h4,"Background PDF","l");
  leg1->AddEntry(h5,"Signal PDF","l");
  leg1->AddEntry(h6,"v_{2} background PDF","l");
  leg1->AddEntry(h7,"Event Selection","l");
  leg1->Draw();
  c0->SaveAs(Form("pdfFiles/Unc_%s.pdf", h0->GetName() ));
  // 7 : Event Selection & Contamination
  // 6 : v2 bkg pdf
  // 5 : CB+Gaus PDF  
  // 4 : background PDF
  // 3 : signal Par
  // 2 : acceptance
  // 1 : efficiency

  
}

void mergeSixInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6) {
  gStyle->SetOptStat(0);
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a6 = h6->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5 + a6*a6);
    h0->SetBinContent( i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);

  h0->SetAxisRange(0,0.03,"Y");
  h0->SetYTitle("Difference Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3,        4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4,        5); h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h5,        6); h5->SetLineWidth(2); h5->DrawCopy("hist same");
  handsomeTH1(h6,        7); h6->SetLineWidth(2); h6->DrawCopy("hist same");
  

  TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
//  easyLeg(leg1,title.Data());
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"efficiency","l");
  leg1->AddEntry(h2,"Acceptance","l");
  leg1->AddEntry(h3,"Signal Parameter","l");
  leg1->AddEntry(h4,"Background PDF","l");
  leg1->AddEntry(h5,"Signal PDF","l");
  leg1->AddEntry(h6,"v_{2} background PDF","l");
  leg1->Draw();
  c0->SaveAs(Form("pdfFiles/Unc_%s.pdf", h0->GetName() ));
  // 6 : v2 bkg pdf
  // 5 : CB+Gaus PDF  
  // 4 : background PDF
  // 3 : signal Par
  // 2 : acceptance
  // 1 : efficiency

  
}

void mergeFiveInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D *h3, TH1D* h4, TH1D* h5) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a5 = h5->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 + a5*a5);
    h0->SetBinContent( i, a0);
  } 

  TCanvas* c0 = new TCanvas("c_mergedSys","",400,400);
  

  h0->SetAxisRange(-0.5,1.1,"Y");
  h0->SetYTitle("Relative Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2);   h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3,4);         h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4,6);         h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h5,11);         h5->SetLineWidth(2); h5->DrawCopy("hist same");
  
  TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
  leg1->AddEntry(h0,"Total","l");
  leg1->AddEntry(h1,"efficiency","l");
  leg1->AddEntry(h2,"Acceptance","l");
  leg1->AddEntry(h3,"Background PDF","l");
  leg1->AddEntry(h4,"Signal PDF","l");
  leg1->AddEntry(h5,"TAA Uncertainty","l");
  leg1->Draw();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.11);
  c0->SaveAs(Form("pdfFiles/postCWR_free_%s_.pdf", h0->GetName() ));
  
  // 5 : TAA uncertainty
  // 4 : CB+Gaus PDF  
  // 3 : background PDF
  // 2 : acceptance
  // 1 : efficiency



}

void mergeFourInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a4 = h4->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 + a4*a4 );
    h0->SetBinContent( i, a0);
  } 
}

void mergeThreeInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a3 = h3->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2 + a3*a3 );
    h0->SetBinContent( i, a0);
  } 
}

void mergeTwoInQuad( TH1D* h0, TH1D* h1, TH1D* h2) {
  for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
    float a1 = h1->GetBinContent(i);
    float a2 = h2->GetBinContent(i);
    float a0 = sqrt( a1*a1 + a2*a2);
    h0->SetBinContent( i, a0);
  } 
}

void subtractTwo( TH1D* h0, TH1D* h1, TH1D* h2) {
  if ( ( h0->GetNbinsX() != h1->GetNbinsX() ) ) {
    cout << "Inconsistent bin numbers!! ERROR" << endl;
  }
  else {
    for ( int i=1 ; i<= h0->GetNbinsX() ;i++){ 
      float a1 = h1->GetBinContent(i);
      float a2 = h2->GetBinContent(i);
      float a0 = (1. + a1) / ( 1. + a2) - 1; 
      h0->SetBinContent( i, a0);
    } 
  }
}
void mergeTwoInQuadCent( TH1D* h0, TH1D* hAA, TH1D* hPP) {
  if ( (hPP->GetNbinsX() != 1 ) )  {
    cout << "Number of hPP bins are not 1!! ERROR" << endl;
  }
  else if ( ( h0->GetNbinsX() != hAA->GetNbinsX() ) ) {
    cout << "Inconsistent bin numbers!! ERROR" << endl;
  }
  else  {
    for ( int i=1 ; i<= h0->GetNbinsX() ;i++){
      float a1 = hAA->GetBinContent(i);
      float a2 = hPP->GetBinContent(1);
      float a0 = sqrt( a1*a1 + a2*a2); 
      h0->SetBinContent( i, a0);
    }
  }
}

void subtractTwoCent( TH1D* h0, TH1D* hAA, TH1D* hPP) {
  if ( (hPP->GetNbinsX() != 1 ) )  {
    cout << "Number of hPP bins are not 1!! ERROR" << endl;
  }
  else if ( ( h0->GetNbinsX() != hAA->GetNbinsX() ) ) {
    cout << "Inconsistent bin numbers!! ERROR" << endl;
  }
  else  {
    for ( int i=1 ; i<= h0->GetNbinsX() ;i++){
      float a1 = hAA->GetBinContent(i);
      float a2 = hPP->GetBinContent(1);
      float a0 = (1. + a1) / ( 1. + a2) - 1;
      h0->SetBinContent( i, a0);
    }
  }
}
