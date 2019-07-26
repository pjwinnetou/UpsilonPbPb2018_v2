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

void totalSyst()
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

	const int nSyst = 7;
	const char* fname[nSyst] = {
	"/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/Acceptance/Acceptance_sys.root",
	"/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/BackgroundPDFVariation/AllParmFreeFit/BackgroundPDFVar_sys.root",
	"/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/Background_v2Variation/SigAllFreeFitFix/V2BkgFuncVar_sys.root",
	"/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/EventSelection/EventSelection_sys.root",
	"/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/SignalPDFVariation/AllParFreeFit/SignalPDFVar_sys.root",
	"/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/SignalParVariation/SignalParVar_sys.root",
	"/Users/hwan/tools/2019/CMS/UpsilonPbPb2018_v2/Systematic/TnP/TnP_sys.root",
	};
	const int nCent = 4;
	const char* CentName[nCent] = {
		"hCentSys_pt03",
  		"hCentSys_pt36",
  		"hCentSys_pt650",
  		"hCentSys_int",
  		//"hSys_int"
	};
	const int nPt = 5;
	const char* hPt[nPt] = {
  		"hPtSys_cent010",
  		"hPtSys_cent1030",
  		"hPtSys_cent3050",
  		"hPtSys_cent5090",
  		"hPtSys_int"
	};

	double unCent[nCent];
	double unPt[nPt];
	double a[nCent];
	double b[nPt];
	TString SystName[nSyst] = { "Acc", "BkgPDF", "BkgV2", "EventSel", "SigPDF", "SigPar", "TnP" };

	TFile *file[nSyst];
	TH1D *hCent[nCent];

	for (int i = 0; i < nSyst; i++)
	{
		file[i] = new TFile(fname[i]);
		cout << fname[i] << endl;
		for (int k = 0; k < nCent; k++)
		{
			unCent[k] = 0;
			hCent[k] = (TH1D*)file[i]->Get(CentName[k]);
			for(int c=1; c<5; c++){ unCent[k] = hCent[k]->GetBinContent(c); cout << CentName[k] << "    " << "Cent Bin " << k << "  " << unCent[k] << endl;}
		}
	}
	

		


}
