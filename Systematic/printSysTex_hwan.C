#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>
#include <Riostream.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TMath.h>
#include <TChain.h>
#include <fstream>
#include <iostream>

void printSysTex_hwan()
{
	using namespace std;

	/*
	TString Acc = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/Acceptance/Acceptance_sys.root";
	TString bkgPdf = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/BackgroundPDFVariation/AllParmFreeFit/BackgroundPDFVar_sys.root";
	TString bkgv2 = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/Background_v2Variation/SigAllFreeFitFix/V2BkgFuncVar_sys.root";
	TString EvntSel = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/EventSelection/EventSelection_sys.root";
	TString SigPdf = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/SignalPDFVariation/AllParFreeFit/SignalPDFVar_sys.root";
	TString SigPar = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/SignalParVariation/SignalParVar_sys.root";
	TString TnP = "/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/TnP/TnP_sys.root";
	*/

	const int nSyst = 9;
	const char* fname[nSyst] = {
	"./Acceptance/Acceptance_sys.root",
	"./Efficiency/Efficiency_sys.root",
	"./TnP/TnP_sys.root",
	"./EventSelection/EventSelection_sys.root",
	"./SignalParVariation/SignalParVar_sys.root",
	"./SignalPDFVariation/AllParFreeFit/SignalPDF_sys.root",
	"./Background_v2Variation/SigAllFreeFitFix/V2Bkg_sys_max.root",
	"./BackgroundPDFVariation/AllParmFreeFit/BackgroundPDF_sys.root",
	"./merged_sys.root",
	};
	const int nCent = 3;
  double centBin[nCent+1] = {10, 30, 50, 90};
	const char* CentName[5] = {
		"hCentSys_pt03",
  		"hCentSys_pt36",
		"hCentSys_pt610",
  		"hCentSys_pt1050",
  		"hCentSys_int",
  		//"hSys_int"
	};
	const char* CentNameMerged[5] = {
		"hCent_Pt03_merged",
  		"hCent_Pt36_merged",
  		"hCent_Pt610_merged",
  		"hCent_Pt1050_merged",
  		"hCent_Pt050_merged",
  		//"hSys_int"
	};
	const int nPt = 4;
  const char* ptBinName[5]= {"0 -- 3", "3 -- 6", "6 -- 10", "6 -- 50", "0 -- 50"};
	const char* ptName[4] = {
  		"hPtSys_cent1030",
  		"hPtSys_cent3050",
  		"hPtSys_cent5090",
		"hPtSys_int"
	};
	const char* ptNameMerged[4] = {
  		"hpt_Cent1030_merged",
  		"hpt_Cent3050_merged",
  		"hpt_Cent5090_merged",
  		"hpt_Cent1090_merged"
	};

	double unCent[nCent];
	double unPt[nPt];
	double a[nCent];
	double b[nPt];
	TString SystName[nSyst] = { "Acc", "BkgPDF", "BkgV2", "EventSel", "SigPDF", "SigPar", "TnP", "Merged" };

	TFile *file[nSyst];
	TH1D *hCent[nCent][nSyst];
	TH1D *hpt[nPt][nSyst];
  TH1D *hInt[nSyst];
  for(int i=0; i<nSyst; i++){
    file[i] = new TFile(fname[i]);
    if(i<nSyst-1){
      for(int icent = 0; icent<nCent; icent++){
        hCent[icent][i] = (TH1D*) file[i]->Get(CentName[icent]);
      }
      for(int ipt=0; ipt<nPt;ipt++){
        hpt[ipt][i] = (TH1D*) file[i]->Get(ptName[ipt]);
      }
      hInt[i] = (TH1D*) file[i]->Get("hSys_int");
    }
    else if(i == nSyst-1){
      for(int icent = 0; icent<nCent; icent++){
        hCent[icent][i] = (TH1D*) file[i]->Get(CentNameMerged[icent]);
      }
      for(int ipt=0; ipt<nPt;ipt++){
        hpt[ipt][i] = (TH1D*) file[i]->Get(ptNameMerged[ipt]);
      }
      hInt[i] = (TH1D*) file[i]->Get("hSys_int_merged");
    }
  }

  for(int icent = 0; icent<nCent+1;icent++){
    //for(int ipt = 0; ipt < nPt-1; ipt++){
    for(int ipt = 0; ipt < nPt; ipt++){
      for(int i=0; i<nSyst; i++){
        if(ipt==0 && icent<nCent && i==0) cout << Form("\\multirow{3}{*}{%.f-%.f\\%%}     &  %s  & ",centBin[icent],centBin[icent+1],ptBinName[ipt]);
        else if(ipt==0 && icent==nCent && i==0) cout << Form("\\multirow{3}{*}{10-90\\%}      &  %s  & ",ptBinName[ipt]);
        else if(ipt!=0 && i==0) cout << Form("                            &  %s  & ",ptBinName[ipt]);
        if(i<nSyst-1 && icent<nCent){
          if(ipt<nPt-1) cout << Form("%.8f",hpt[icent][i]->GetBinContent(ipt+1)) << " & ";
          if(ipt==nPt-1) cout << Form("%.8f",hCent[nCent-1][i]->GetBinContent(icent+1)) << " & ";
        }
        else if(i==nSyst-1 && icent<nCent){
          if(ipt<nPt-1) cout << Form("%.8f",hpt[icent][i]->GetBinContent(ipt+1)) << "  \\\\ " << endl;
          if(ipt==nPt-1) cout << Form("%.8f",hCent[nCent-1][i]->GetBinContent(1+icent)) << "  \\\\ \\hline " << endl;
        }
        if(icent==nCent){
          if(i<nSyst-1 && ipt<nPt-1) cout << Form("%.8f",hpt[nCent][i]->GetBinContent(ipt+1)) << " & ";
          if(i<nSyst-1 && ipt==nPt-1) cout << Form("%.8f",hInt[i]->GetBinContent(1)) << " & ";
          if(i==nSyst-1 && ipt<nPt-1) cout << Form("%.8f",hpt[nCent][i]->GetBinContent(ipt+1)) << "  \\\\ " << endl;
          if(i==nSyst-1 && ipt==nPt-1) cout << Form("%.8f",hInt[i]->GetBinContent(1)) << "  \\\\ \\hline " << endl;
        }
      } 
    }
  }
  
    







}
