#include "../../../commonUtility.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TFile.h"
#include "TROOT.h"
#include "../../../cutsAndBinUpsilonV2.h"
#include <TGraphErrors.h>

using namespace std;

double getSys(TString fileName="Dir1", TString fileName2="Dir2");

void compute_sys_diff()
{
  const int nPtBin = 3;
  double PtBinArr[nPtBin+1] = {0, 3, 6, 50};
  const int nCentBin = 4;
  double CentBinArr[nCentBin+1] = {0, 10, 30, 50, 100};

  TH1D* hPtSys_cent[nCentBin];
  TH1D* hCentSys_pt[nPtBin];
  for(int ib=0; ib<nCentBin; ib++){
    hPtSys_cent[ib] = new TH1D(Form("hPtSys_cent%.f%.f",CentBinArr[ib],CentBinArr[ib+1]),";p_{T};Uncertainty",nPtBin,PtBinArr);
  }
  for(int ib=0; ib<nPtBin; ib++){
    hCentSys_pt[ib] = new TH1D(Form("hCentSys_pt%.f%.f",PtBinArr[ib],PtBinArr[ib+1]),";Centrality (%);Uncertainty",nCentBin,CentBinArr);
  }
  TH1D* hPtSys_int = new TH1D("hPtSys_int",";p_{T};Uncertainty",nPtBin,PtBinArr);
  TH1D* hCentSys_int = new TH1D("hCentSys_int",";Centrality (%);Uncertainty",nCentBin,CentBinArr);

  TH1D* hSys_int = new TH1D("hSys_int",";Integrated;Uncertainty",1,0,1);
  
  TString dirNom = "/home/deathold/work/CMS/analysis/Upsilon_v2/upsilonV2/SimFit/SigAllFreeFitFix/FitResult";
  TString dirAlt = "/home/deathold/work/CMS/analysis/Upsilon_v2/upsilonV2/Systematic/SignalPDFVariation/AllParFreeFit/FixSignalFromAllFreeFit/FitResult";

  double v2_sys = 0;

  cout << "Bin" << endl;
  for(int ib=0; ib<nCentBin; ib++){
    cout << "Cent : " << CentBinArr[ib] << "-" << CentBinArr[ib+1] << endl; 
    for(int ipt =0; ipt<nPtBin; ipt++){
      cout << "Pt : " << PtBinArr[ipt] << "-" << PtBinArr[ipt+1] << "GeV/c" << endl;
      TString fileNom = dirNom+Form("/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",PtBinArr[ipt],PtBinArr[ipt+1],CentBinArr[ib]*2,CentBinArr[ib+1]*2);
      TString fileAlt = dirAlt+Form("/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",PtBinArr[ipt],PtBinArr[ipt+1],CentBinArr[ib]*2,CentBinArr[ib+1]*2);
      v2_sys = getSys(fileNom,fileAlt);
      hPtSys_cent[ib]  -> SetBinContent(ipt+1, v2_sys); 
      hCentSys_pt[ipt] -> SetBinContent(ib+1, v2_sys); 
      v2_sys=0;
    }
  }

  //integrated
  for(int ib=0; ib<nCentBin; ib++){
    TString fileNom = dirNom+Form("/SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",CentBinArr[ib]*2,CentBinArr[ib+1]*2);
    TString fileAlt = dirAlt+Form("/SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",CentBinArr[ib]*2,CentBinArr[ib+1]*2);
    v2_sys = getSys(fileNom,fileAlt);
    hCentSys_int->SetBinContent(ib+1, v2_sys);
    v2_sys=0;
  }

  for(int ib=0; ib<nPtBin; ib++){
    TString fileNom = dirNom+Form("/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality0-200_m8-14_OS.root",PtBinArr[ib],PtBinArr[ib+1]);
    TString fileAlt = dirAlt+Form("/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality0-200_m8-14_OS.root",PtBinArr[ib],PtBinArr[ib+1]);
    v2_sys = getSys(fileNom,fileAlt);
    hPtSys_int->SetBinContent(ib+1, v2_sys);
    v2_sys=0;
  }
  
  TString fileIntNom = dirNom + "/SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality0-200_m8-14_OS.root";
  TString fileIntAlt = dirAlt + "/SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality0-200_m8-14_OS.root";
  double v2_sys_int = getSys(fileIntNom, fileIntAlt);
  hSys_int->SetBinContent(1,v2_sys_int);
  

  TFile *wf = new TFile("SignalPDFVar_sys.root","recreate");
  wf->cd();
  for(int ib=0; ib<nPtBin; ib++){
    hCentSys_pt[ib]->Write();
  }
  for(int ib=0; ib<nCentBin; ib++){
    hPtSys_cent[ib]->Write();
  }
  hCentSys_int->Write();
  hPtSys_int->Write();
  hSys_int->Write();

}


double getSys(TString fileName, TString fileName2)
{ 
  double v2_nom = -999;
  double v2_alt = -999;
  double x = -999;
  double sys = -999; 
  TString fin_nom;
  fin_nom = Form("%s",fileName.Data());
  TFile* f1 = new TFile(fin_nom.Data() );
  TGraphErrors *g1 = (TGraphErrors*) f1->Get("v2vspt");
  g1->GetPoint(0,x,v2_nom);
  if(x==-999 || v2_nom == -999){cout << "ERROR!!" << endl; return 1;}

  x = -999;
  TString fin_alt;
  fin_alt = Form("%s",fileName2.Data());
  TFile* f2 = new TFile(fin_alt.Data() );
  TGraphErrors *g2 = (TGraphErrors*) f2->Get("v2vspt");
  g2->GetPoint(0,x,v2_alt);
  if(x==-999 || v2_nom == -999 || v2_alt == -999){cout << "ERROR!!" << endl; return 1;}
  
  sys = fabs(v2_alt-v2_nom);
  cout << "Unc : " << Form("%.3f",100*double(sys)) << "%%" << endl;
  return sys;
}


