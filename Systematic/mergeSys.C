#include "../commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../Style_jaebeom.h"
using namespace std;


//// do NOT use "hadded" ttrees!! (e.g.6-100 GeV) 

TLegend *leg = new TLegend(0.55,0.2, 0.85,0.4,NULL,"brNDC");
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

 
  TH1D* hpt_Cent0100[10];
  TH1D* hpt_Cent010[10];
  TH1D* hpt_Cent1030[10];
  TH1D* hpt_Cent3050[10];
  TH1D* hpt_Cent50100[10];

  TH1D* hCent_Pt050[10];
  TH1D* hCent_Pt03[10];
  TH1D* hCent_Pt36[10];
  TH1D* hCent_Pt650[10];
  
  TH1D* hInt[10];

/*
  // 1 : efficiency
  TFile* f1 = new TFile("../Efficiency/sys.root");
  hpt_Cent0100[1]   = (TH1D*) f1->Get("hPtSys_int");
  hpt_Cent010[1]    = (TH1D*) f1->Get("hPtSys_cent010");
  hpt_Cent1030[1]   = (TH1D*) f1->Get("hPtSys_cent1030");
  hpt_Cent3050[1]   = (TH1D*) f1->Get("hPtSys_cent3050");
  hpt_Cent50100[1]  = (TH1D*) f1->Get("hPtSys_cent50100");

  hCent_Pt050[1]   = (TH1D*) f1->Get("hCentSys_int");
  hCent_Pt03[1]    = (TH1D*) f1->Get("hCentSys_pt03");
  hCent_Pt36[1]    = (TH1D*) f1->Get("hCentSys_pt36");
  hCent_Pt650[1]   = (TH1D*) f1->Get("hCentSys_pt650");
  hInt[1]   = (TH1D*) f1->Get("hSys_int");


  // 2 : acceptance
  TFile* f2 = new TFile("../Acceptance/sys.root");
  hpt_Cent0100[2]   = (TH1D*) f2->Get("hPtSys_int");
  hpt_Cent010[2]    = (TH1D*) f2->Get("hPtSys_cent010");
  hpt_Cent1030[2]   = (TH1D*) f2->Get("hPtSys_cent1030");
  hpt_Cent3050[2]   = (TH1D*) f2->Get("hPtSys_cent3050");
  hpt_Cent50100[2]  = (TH1D*) f2->Get("hPtSys_cent50100");

  hCent_Pt050[2]   = (TH1D*) f2->Get("hCentSys_int");
  hCent_Pt03[2]    = (TH1D*) f2->Get("hCentSys_pt03");
  hCent_Pt36[2]    = (TH1D*) f2->Get("hCentSys_pt36");
  hCent_Pt650[2]   = (TH1D*) f2->Get("hCentSys_pt650");
  hInt[2]   = (TH1D*) f2->Get("hSys_int");


  // 3 : signal Par
  TFile* f3 = new TFile(Form("SignalVariation/sys.root");
  hpt_Cent0100[3]   = (TH1D*) f3->Get("hPtSys_int");
  hpt_Cent010[3]    = (TH1D*) f3->Get("hPtSys_cent010");
  hpt_Cent1030[3]   = (TH1D*) f3->Get("hPtSys_cent1030");
  hpt_Cent3050[3]   = (TH1D*) f3->Get("hPtSys_cent3050");
  hpt_Cent50100[3]  = (TH1D*) f3->Get("hPtSys_cent50100");

  hCent_Pt050[3]   = (TH1D*) f3->Get("hCentSys_int");
  hCent_Pt03[3]    = (TH1D*) f3->Get("hCentSys_pt03");
  hCent_Pt36[3]    = (TH1D*) f3->Get("hCentSys_pt36");
  hCent_Pt650[3]   = (TH1D*) f3->Get("hCentSys_pt650");
  hInt[3]   = (TH1D*) f3->Get("hSys_int");
*/

  // 4 : signal PDF var
  TFile* f4 = new TFile("SignalPDFVariation/AllParFreeFit/SignalPDFVar_sys.root");
  hpt_Cent0100[4]   = (TH1D*) f4->Get("hPtSys_int");
  hpt_Cent010[4]    = (TH1D*) f4->Get("hPtSys_cent010");
  hpt_Cent1030[4]   = (TH1D*) f4->Get("hPtSys_cent1030");
  hpt_Cent3050[4]   = (TH1D*) f4->Get("hPtSys_cent3050");
  hpt_Cent50100[4]  = (TH1D*) f4->Get("hPtSys_cent50100");

  hCent_Pt050[4]   = (TH1D*) f4->Get("hCentSys_int");
  hCent_Pt03[4]    = (TH1D*) f4->Get("hCentSys_pt03");
  hCent_Pt36[4]    = (TH1D*) f4->Get("hCentSys_pt36");
  hCent_Pt650[4]   = (TH1D*) f4->Get("hCentSys_pt650");
  hInt[4]   = (TH1D*) f4->Get("hSys_int");

  cout << "asdas" << endl;

  // 5 : background PDF var
  TFile* f5 = new TFile("BackgroundPDFVariation/AllParmFreeFit/BackgroundPDFVar_sys.root");
  hpt_Cent0100[5]   = (TH1D*) f5->Get("hPtSys_int");
  hpt_Cent010[5]    = (TH1D*) f5->Get("hPtSys_cent010");
  hpt_Cent1030[5]   = (TH1D*) f5->Get("hPtSys_cent1030");
  hpt_Cent3050[5]   = (TH1D*) f5->Get("hPtSys_cent3050");
  hpt_Cent50100[5]  = (TH1D*) f5->Get("hPtSys_cent50100");

  hCent_Pt050[5]   = (TH1D*) f5->Get("hCentSys_int");
  hCent_Pt03[5]    = (TH1D*) f5->Get("hCentSys_pt03");
  hCent_Pt36[5]    = (TH1D*) f5->Get("hCentSys_pt36");
  hCent_Pt650[5]   = (TH1D*) f5->Get("hCentSys_pt650");
  hInt[5]   = (TH1D*) f5->Get("hSys_int");
  cout << "asdas" << endl;


  // 6 : v2 background PDF var
  TFile* f6 = new TFile("Background_v2Variation/SigAllFreeFitFix/V2BkgFuncVar_sys.root");
  hpt_Cent0100[6]   = (TH1D*) f6->Get("hPtSys_int");
  hpt_Cent010[6]    = (TH1D*) f6->Get("hPtSys_cent010");
  hpt_Cent1030[6]   = (TH1D*) f6->Get("hPtSys_cent1030");
  hpt_Cent3050[6]   = (TH1D*) f6->Get("hPtSys_cent3050");
  hpt_Cent50100[6]  = (TH1D*) f6->Get("hPtSys_cent50100");

  hCent_Pt050[6]   = (TH1D*) f6->Get("hCentSys_int");
  hCent_Pt03[6]    = (TH1D*) f6->Get("hCentSys_pt03");
  hCent_Pt36[6]    = (TH1D*) f6->Get("hCentSys_pt36");
  hCent_Pt650[6]   = (TH1D*) f6->Get("hCentSys_pt650");
  hInt[6]   = (TH1D*) f6->Get("hSys_int");
  cout << "asdas" << endl;



  hpt_Cent0100[0]  = (TH1D*) hpt_Cent0100[5]  -> Clone("hpt_Cent0100_merged"); hpt_Cent0100[0]->Reset();
  hpt_Cent010[0]   = (TH1D*) hpt_Cent010[5]   -> Clone("hpt_Cent010_merged"); hpt_Cent010[0]->Reset();
  hpt_Cent1030[0]  = (TH1D*) hpt_Cent1030[5]  -> Clone("hpt_Cent1030_merged"); hpt_Cent1030[0]->Reset();
  hpt_Cent3050[0]  = (TH1D*) hpt_Cent3050[5]  -> Clone("hpt_Cent3050_merged"); hpt_Cent3050[0]->Reset();
  hpt_Cent50100[0] = (TH1D*) hpt_Cent50100[5] -> Clone("hpt_Cent50100_merged"); hpt_Cent50100[0]->Reset();

  hCent_Pt050[0]  = (TH1D*) hCent_Pt050[5]  -> Clone("hCent_Pt050_merged"); hCent_Pt050[0]->Reset(); 
  hCent_Pt03[0]   = (TH1D*) hCent_Pt03[5]   -> Clone("hCent_Pt03_merged");  hCent_Pt03[0] ->Reset(); 
  hCent_Pt36[0]   = (TH1D*) hCent_Pt36[5]   -> Clone("hCent_Pt36_merged");  hCent_Pt36[0] ->Reset(); 
  hCent_Pt650[0]  = (TH1D*) hCent_Pt650[5]  -> Clone("hCent_Pt650_merged"); hCent_Pt650[0]->Reset(); 
  hCent_Pt650[0]  = (TH1D*) hCent_Pt650[5]  -> Clone("hCent_Pt650_merged"); hCent_Pt650[0]->Reset(); 
  hInt[0]  = (TH1D*) hInt[5]  -> Clone("hSys_int_merged"); hInt[0]->Reset();
  cout << "asdas" << endl;


  // Merge uncertainties for cross-section 
/*
  mergeSixInQuad(hpt_Cent0100[0], hpt_Cent0100[1], hpt_Cent0100[2], hpt_Cent0100[3], hpt_Cent0100[4], hpt_Cent0100[5], hpt_Cent0100[6]);  
  mergeSixInQuad(hpt_Cent010[0], hpt_Cent010[1], hpt_Cent010[2], hpt_Cent010[3], hpt_Cent010[4], hpt_Cent010[5], hpt_Cent010[6]);  
  mergeSixInQuad(hpt_Cent1030[0], hpt_Cent1030[1], hpt_Cent1030[2], hpt_Cent1030[3], hpt_Cent1030[4], hpt_Cent1030[5], hpt_Cent1030[6]);  
  mergeSixInQuad(hpt_Cent3050[0], hpt_Cent3050[1], hpt_Cent3050[2], hpt_Cent3050[3], hpt_Cent3050[4], hpt_Cent3050[5], hpt_Cent3050[6]);  
  mergeSixInQuad(hpt_Cent50100[0], hpt_Cent50100[1], hpt_Cent50100[2], hpt_Cent50100[3], hpt_Cent50100[4], hpt_Cent50100[5], hpt_Cent50100[6]);  

  mergeSixInQuad(hCent_Pt050[0], hCent_Pt050[1], hCent_Pt050[2], hCent_Pt050[3], hCent_Pt050[4], hCent_Pt050[5], hCent_Pt050[6]);  
  mergeSixInQuad(hCent_Pt03[0], hCent_Pt03[1], hCent_Pt03[2], hCent_Pt03[3], hCent_Pt03[4], hCent_Pt03[5], hCent_Pt03[6]);  
  mergeSixInQuad(hCent_Pt36[0], hCent_Pt36[1], hCent_Pt36[2], hCent_Pt36[3], hCent_Pt36[4], hCent_Pt36[5], hCent_Pt36[6]);  
  mergeSixInQuad(hCent_Pt650[0], hCent_Pt650[1], hCent_Pt650[2], hCent_Pt650[3], hCent_Pt650[4], hCent_Pt650[5], hCent_Pt650[6]);  
  */
  mergeThreeInQuad(hpt_Cent0100[0],   hpt_Cent0100[4], hpt_Cent0100[5],  hpt_Cent0100[6]);
  mergeThreeInQuad(hpt_Cent010[0],     hpt_Cent010[4],  hpt_Cent010[5],   hpt_Cent010[6]);  
  mergeThreeInQuad(hpt_Cent1030[0],   hpt_Cent1030[4], hpt_Cent1030[5],  hpt_Cent1030[6]);  
  mergeThreeInQuad(hpt_Cent3050[0],   hpt_Cent3050[4], hpt_Cent3050[5],  hpt_Cent3050[6]);  
  mergeThreeInQuad(hpt_Cent50100[0], hpt_Cent50100[4],hpt_Cent50100[5], hpt_Cent50100[6]);  

  mergeThreeInQuad(hCent_Pt050[0],     hCent_Pt050[4],  hCent_Pt050[5],   hCent_Pt050[6]);  
  mergeThreeInQuad(hCent_Pt03[0],       hCent_Pt03[4],   hCent_Pt03[5],    hCent_Pt03[6]);  
  mergeThreeInQuad(hCent_Pt36[0],       hCent_Pt36[4],   hCent_Pt36[5],    hCent_Pt36[6]);  
  mergeThreeInQuad(hCent_Pt650[0],     hCent_Pt650[4],  hCent_Pt650[5],   hCent_Pt650[6]);  
  mergeThreeInQuad(hInt[0], hInt[4], hInt[5], hInt[6]);  


  TFile* fout = new TFile("merged_sys.root","recreate");
  hpt_Cent0100[0]->Write();
  hpt_Cent010[0]->Write();
  hpt_Cent1030[0]->Write();
  hpt_Cent3050[0]->Write();
  hpt_Cent50100[0]->Write();

  hCent_Pt050[0]->Write();
  hCent_Pt03[0]->Write();
  hCent_Pt36[0]->Write();
  hCent_Pt650[0]->Write();
  hInt[0]->Write();

  fout->Close();


}


void mergeSixInQuad( TH1D* h0, TH1D* h1, TH1D* h2, TH1D* h3, TH1D* h4, TH1D* h5, TH1D* h6) {
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

  h0->SetAxisRange(0,0.4,"Y");
  h0->SetAxisRange(0,0.7,"Y");
  h0->SetYTitle("Relative Uncertainty");
  handsomeTH1(h0,        1); h0->SetLineWidth(2); h0->DrawCopy("hist");
  handsomeTH1(h1,        2); h1->SetLineWidth(2); h1->DrawCopy("hist same");
  handsomeTH1(h2,        3); h2->SetLineWidth(2); h2->DrawCopy("hist same");
  handsomeTH1(h3,        4); h3->SetLineWidth(2); h3->DrawCopy("hist same");
  handsomeTH1(h4,        5); h4->SetLineWidth(2); h4->DrawCopy("hist same");
  handsomeTH1(h5,        6); h5->SetLineWidth(2); h5->DrawCopy("hist same");
  handsomeTH1(h6,        7); h6->SetLineWidth(2); h6->DrawCopy("hist same");
  

  TLegend *leg1 = new TLegend(0.55,0.6, 0.85,0.9,NULL,"brNDC");
  //easyLeg(leg1,title.Data());
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
