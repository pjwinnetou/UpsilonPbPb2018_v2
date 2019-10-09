#include "../../commonUtility.h"
#include "TText.h"
#include "TLine.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TFile.h"
#include "TROOT.h"
#include "../../cutsAndBinUpsilonV2.h"
#include <TGraphErrors.h>

using namespace std;

double getSys(TString fileName="Dir1", TString fileName2="Dir2");

void compute_sys_Acc_diff()
{
  const int nPtBin = 4;
  const int nPtBinY = 3;
  double PtBinArr[nPtBin+1] = {0, 3, 6, 10, 50};
  const int nCentBin = 3;
  double CentBinArr[nCentBin+1] = {0, 10, 50, 90};
  const int nCentBinInt = 4;
  double CentBinArrInt[nCentBinInt+1] = {0, 10, 30, 50, 90};
  
  const int nYBin = 2;
  double YBinArr[nYBin+1] = {0, 1.2, 2.4};
  double PtBinArrC[nPtBin+1] = {0, 3, 6, 15};

  TH1D* hPtSys_cent[nCentBin];
  TH1D* hCentSys_pt[nPtBin];
  for(int ib=0; ib<nCentBin; ib++){
    hPtSys_cent[ib] = new TH1D(Form("hPtSys_cent%.f%.f",CentBinArr[ib],CentBinArr[ib+1]),";p_{T};Uncertainty",nPtBin,PtBinArr);
  }
  for(int ib=0; ib<nPtBin; ib++){
    hCentSys_pt[ib] = new TH1D(Form("hCentSys_pt%.f%.f",PtBinArr[ib],PtBinArr[ib+1]),";Centrality (%);Uncertainty",nCentBin,CentBinArr);
  }
  TH1D* hPtSys_int = new TH1D("hPtSys_int",";p_{T};Uncertainty",nPtBin,PtBinArr);
  TH1D* hCentSys_int = new TH1D("hCentSys_int",";Centrality (%);Uncertainty",nCentBinInt,CentBinArrInt);

  TH1D* hSys_int = new TH1D("hSys_int",";Integrated;Uncertainty",1,0,1);


  TString strRap[nYBin] = {"012","1224"};
  TH1D* hPtSys_rap[nYBin];
  for(int irap=0; irap<nYBin; irap++){
    hPtSys_rap[irap] = new TH1D(Form("hPtSys_rap%s",strRap[irap].Data()),";p_{T};Uncertainty",nPtBinY,PtBinArrC);
  }
  TH1D* hPtSys_rapInt = new TH1D("hPtSys_rap024",";p_{T};Uncertainty",nPtBinY,PtBinArrC);
  
  TString dirNom = "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/SimFit/SigAllFreeFitFix/FitResult";
  TString dirAlt = "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Systematic/Acceptance/FitResult";

  double v2_sys = 0;

  cout << "Bin" << endl;
  for(int ib=0; ib<nCentBin; ib++){
    cout << "Cent : " << CentBinArr[ib] << "-" << CentBinArr[ib+1] << endl; 
    for(int ipt =0; ipt<nPtBin; ipt++){
      cout << "Pt : " << PtBinArr[ipt] << "-" << PtBinArr[ipt+1] << "GeV/c" << endl;
      TString fileNom = dirNom+Form("/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",PtBinArr[ipt],PtBinArr[ipt+1],CentBinArr[ib]*2,CentBinArr[ib+1]*2);
      TString fileAlt = dirAlt+Form("/SimFitResult_SysAcc_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",PtBinArr[ipt],PtBinArr[ipt+1],CentBinArr[ib]*2,CentBinArr[ib+1]*2);
      v2_sys = getSys(fileNom,fileAlt);
      hPtSys_cent[ib]  -> SetBinContent(ipt+1, v2_sys); 
      hCentSys_pt[ipt] -> SetBinContent(ib+1, v2_sys); 
      v2_sys=0;
    }
  }

  //integrated
  cout << endl;
  cout << "Integrated Pt 0-50 GeV/c" << endl;
  for(int ib=0; ib<nCentBinInt; ib++){
    cout << "Cent : " << CentBinArrInt[ib] << "-" << CentBinArrInt[ib+1] << endl;
    TString fileNom = dirNom+Form("/SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",CentBinArrInt[ib]*2,CentBinArrInt[ib+1]*2);
    TString fileAlt = dirAlt+Form("/SimFitResult_SysAcc_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",CentBinArrInt[ib]*2,CentBinArrInt[ib+1]*2);
    v2_sys = getSys(fileNom,fileAlt);
    hCentSys_int->SetBinContent(ib+1, v2_sys);
    v2_sys=0;
  }

  cout << endl;
  cout << "Integrated Cent 10-90%" << endl;
  for(int ib=0; ib<nPtBin; ib++){
    cout << "Cent : " << PtBinArr[ib] << "-" << PtBinArr[ib+1] << " GeV/c" << endl;
    TString fileNom = dirNom+Form("/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality20-180_m8-14_OS.root",PtBinArr[ib],PtBinArr[ib+1]);
    TString fileAlt = dirAlt+Form("/SimFitResult_SysAcc_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality20-180_m8-14_OS.root",PtBinArr[ib],PtBinArr[ib+1]);
    v2_sys = getSys(fileNom,fileAlt);
    hPtSys_int->SetBinContent(ib+1, v2_sys);
    v2_sys=0;
  }
  
  TString fileIntNom = dirNom + "/SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality20-180_m8-14_OS.root";
  TString fileIntAlt = dirAlt + "/SimFitResult_SysAcc_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality20-180_m8-14_OS.root";
  double v2_sys_int = getSys(fileIntNom, fileIntAlt);
  hSys_int->SetBinContent(1,v2_sys_int);
  
  
  //Comp ALICE bin
  v2_sys=0;
  cout << endl;
  cout << endl;
  cout << "Bins for comparison with ALICE" << endl;
  cout << endl;
  for(int ib=0; ib<nYBin; ib++){
    cout << "Rapidity : " << YBinArr[ib] << "-" << YBinArr[ib+1] << endl;
    for(int ipt=0; ipt<nPtBinY; ipt++){
      TString fileNom = dirNom+Form("/SimFitResult_pt%.1f-%.1f_y%.1f-%.1f_muPt3.5_centrality10-120_m8-14_OS.root",PtBinArrC[ipt],PtBinArrC[ipt+1],YBinArr[ib],YBinArr[ib+1]);
      TString fileAlt = dirAlt+Form("/SimFitResult_SysAcc_pt%.1f-%.1f_y%.1f-%.1f_muPt3.5_centrality10-120_m8-14_OS.root",PtBinArrC[ipt],PtBinArrC[ipt+1],YBinArr[ib],YBinArr[ib+1]);
      v2_sys = getSys(fileNom,fileAlt);
      hPtSys_rap[ib]->SetBinContent(ipt+1, v2_sys);
      v2_sys=0;
      if(ib==0){
        cout << "Full Rap" << endl;
        fileNom = dirNom+Form("/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality10-120_m8-14_OS.root",PtBinArrC[ipt],PtBinArrC[ipt+1]);
        fileAlt = dirAlt+Form("/SimFitResult_SysAcc_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality10-120_m8-14_OS.root",PtBinArrC[ipt],PtBinArrC[ipt+1]);
        v2_sys = getSys(fileNom,fileAlt);
        hPtSys_rapInt->SetBinContent(ipt+1, v2_sys);
        v2_sys=0;
        cout << endl;
      }
    }
  }


  TFile *wf = new TFile("Acceptance_sys.root","recreate");
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


  for(int irap=0; irap<nYBin; irap++){
    hPtSys_rap[irap]->Write();
  }
  hPtSys_rapInt->Write();

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
  cout << "Unc : " << Form("%.5f",double(sys)) << endl;
  return sys;
}

