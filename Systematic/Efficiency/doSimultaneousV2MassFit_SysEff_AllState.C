//Headers{{{
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <RooFit.h>
#include <RooGlobalFunc.h>
#include <RooCategory.h>
#include <RooGenericPdf.h>
#include <RooFFTConvPdf.h>
#include <RooWorkspace.h>
#include <RooBinning.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooHist.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooConstVar.h>

#include <Fit/Fitter.h>
#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Math/WrappedMultiTF1.h>
#include <HFitInterface.h>

#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/commonUtility.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/cutsAndBinUpsilonV2.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/HiEvtPlaneList.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Style_jaebeom.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/tdrstyle.C"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/CMS_lumi_v2mass.C"

//}}}

using namespace std;

const int nParmM = 13;
const int nParmV = 20;
Int_t iparmass[nParmM] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
//Int_t iparvn[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
Int_t iparvn[nParmV] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};


struct GlobalChi2_width
{
	GlobalChi2_width(ROOT::Math::IMultiGenFunction & f1,
						ROOT::Math::IMultiGenFunction & f2):
	fChi2_1(&f1), fChi2_2(&f2) {}

	Double_t operator() (const double *par) const
	{
		Double_t p1[nParmM];
		for(Int_t i = 0; i < nParmM; i++) p1[i] = par[iparmass[i]];
		Double_t p2[nParmV];
		for(Int_t i = 0; i < nParmV; i++) p2[i] = par[iparvn[i]];
		return (*fChi2_1)(p1) + (*fChi2_2)(p2);
	}
	const ROOT::Math::IMultiGenFunction * fChi2_1;
	const ROOT::Math::IMultiGenFunction * fChi2_2;
};

//totalYield{{{
Double_t TotalYield(Double_t* x, Double_t* par)
{
	Double_t N1 = par[0];
	Double_t N2 = par[1];
	Double_t N3 = par[2];
	Double_t Nbkg = par[3];
	Double_t mean = par[4];
	Double_t sigma = par[5];
	Double_t alpha = par[6];
	Double_t n = par[7];
	Double_t ratio = par[8];
	Double_t frac = par[9];
	Double_t Bkgmean = par[10];
	Double_t Bkgsigma = par[11];
	Double_t Bkgp0 = par[12];
	Double_t mean2 = mean*pdgMass.Y2S/pdgMass.Y1S;
	Double_t mean3 = mean*pdgMass.Y3S/pdgMass.Y1S;
	Double_t sigma1_2 = sigma*ratio;
	Double_t sigma2_1 = sigma*pdgMass.Y2S/pdgMass.Y1S;
	Double_t sigma2_2 = sigma*ratio*pdgMass.Y2S/pdgMass.Y1S;
	Double_t sigma3_1 = sigma*pdgMass.Y3S/pdgMass.Y1S;
	Double_t sigma3_2 = sigma*ratio*pdgMass.Y3S/pdgMass.Y1S;

	//t2 > t1
	Double_t U1S_t1 = (x[0]-mean)/sigma;
	Double_t U1S_t2 = (x[0]-mean)/sigma1_2;
	Double_t U2S_t1 = (x[0]-mean2)/sigma2_1;
	Double_t U2S_t2 = (x[0]-mean2)/sigma2_2;
	Double_t U3S_t1 = (x[0]-mean3)/sigma3_1;
	Double_t U3S_t2 = (x[0]-mean3)/sigma3_2;
	if (alpha < 0)
	{

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
	}

	Double_t absAlpha = TMath::Abs(alpha);
	Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
	Double_t b = n/absAlpha - absAlpha;


  Double_t Ups1S_1 = -1;
  Double_t Ups1S_2 = -1;
  Double_t Ups2S_1 = -1;
  Double_t Ups2S_2 = -1;
  Double_t Ups3S_1 = -1;
  Double_t Ups3S_2 = -1;


  if(U1S_t1 > -alpha){
    Ups1S_1 = exp(-U1S_t1*U1S_t1/2.);
  }
  if(U1S_t1 <= -alpha){
    Ups1S_1 = a*TMath::Power((b-U1S_t1),-n);
  }

  if(U1S_t2 > -alpha){
    Ups1S_2 = exp(-U1S_t2*U1S_t2/2.);
  }
  if(U1S_t2 <= -alpha){
    Ups1S_2 = a*TMath::Power((b-U1S_t2),-n);
  }

  if(U2S_t1 > -alpha){
    Ups2S_1 = exp(-U2S_t1*U2S_t1/2.);
  }
  if(U2S_t1 <= -alpha){
    Ups2S_1 = a*TMath::Power((b-U2S_t1),-n);
  }

  if(U2S_t2 > -alpha){
    Ups2S_2 = exp(-U2S_t2*U2S_t2/2.);
  }
  if(U2S_t2 <= -alpha){
    Ups2S_2 = a*TMath::Power((b-U2S_t2),-n);
  }

  if(U3S_t1 > -alpha){
    Ups3S_1 = exp(-U3S_t1*U3S_t1/2.);
  }
  if(U3S_t1 <= -alpha){
    Ups3S_1 = a*TMath::Power((b-U3S_t1),-n);
  }

  if(U3S_t2 > -alpha){
    Ups3S_2 = exp(-U3S_t2*U3S_t2/2.);
  }
  if(U3S_t2 <= -alpha){
    Ups3S_2 = a*TMath::Power((b-U3S_t2),-n);
  }
  
  Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
  Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
  Double_t fN1s_1 = 1./(sigma*(fC+fD));
  Double_t fN1s_2 = 1./(sigma1_2*(fC+fD));
  Double_t fN2s_1 = 1./(sigma2_1*(fC+fD));
  Double_t fN2s_2 = 1./(sigma2_2*(fC+fD));
  Double_t fN3s_1 = 1./(sigma3_1*(fC+fD));
  Double_t fN3s_2 = 1./(sigma3_2*(fC+fD));

  Double_t SigM = N1*(fN1s_1*frac*Ups1S_1 + fN1s_2*(1-frac)*Ups1S_2) + N2*(fN2s_1*frac*Ups2S_1 + fN2s_2*(1-frac)*Ups2S_2) + N3*(fN3s_1*frac*Ups3S_1 + fN3s_2*(1-frac)*Ups3S_2);
  Double_t BkgM = Nbkg*(TMath::Exp(-x[0]/Bkgp0)*(TMath::Erf((x[0]-Bkgmean)/(TMath::Sqrt(2)*Bkgsigma))+1)/2.);
  
  return SigM + BkgM;
}
//}}}


//totalYieldSig{{{
Double_t TotalYieldSig(Double_t* x, Double_t* par)
{
	Double_t N1 = par[0];
	Double_t N2 = par[1];
	Double_t N3 = par[2];
	Double_t mean = par[3];
	Double_t sigma = par[4];
	Double_t alpha = par[5];
	Double_t n = par[6];
	Double_t ratio = par[7];
	Double_t frac = par[8];
	Double_t mean2 = mean*pdgMass.Y2S/pdgMass.Y1S;
	Double_t mean3 = mean*pdgMass.Y3S/pdgMass.Y1S;
	Double_t sigma1_2 = sigma*ratio;
	Double_t sigma2_1 = sigma*pdgMass.Y2S/pdgMass.Y1S;
	Double_t sigma2_2 = sigma*ratio*pdgMass.Y2S/pdgMass.Y1S;
	Double_t sigma3_1 = sigma*pdgMass.Y3S/pdgMass.Y1S;
	Double_t sigma3_2 = sigma*ratio*pdgMass.Y3S/pdgMass.Y1S;

	//t2 > t1
	Double_t U1S_t1 = (x[0]-mean)/sigma;
	Double_t U1S_t2 = (x[0]-mean)/sigma1_2;
	Double_t U2S_t1 = (x[0]-mean2)/sigma2_1;
	Double_t U2S_t2 = (x[0]-mean2)/sigma2_2;
	Double_t U3S_t1 = (x[0]-mean3)/sigma3_1;
	Double_t U3S_t2 = (x[0]-mean3)/sigma3_2;
	if (alpha < 0)
	{

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
	}

	Double_t absAlpha = fabs((Double_t)alpha);
	Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
	Double_t b = n/absAlpha - absAlpha;


  Double_t Ups1S_1 = -1;
  Double_t Ups1S_2 = -1;
  Double_t Ups2S_1 = -1;
  Double_t Ups2S_2 = -1;
  Double_t Ups3S_1 = -1;
  Double_t Ups3S_2 = -1;


  if(U1S_t1 > -alpha){
    Ups1S_1 = exp(-U1S_t1*U1S_t1/2.);
  }
  else if(U1S_t1 <= -alpha){
    Ups1S_1 = a*TMath::Power((b-U1S_t1),-n);
  }

  if(U1S_t2 > -alpha){
    Ups1S_2 = exp(-U1S_t2*U1S_t2/2.);
  }
  else if(U1S_t2 <= -alpha){
    Ups1S_2 = a*TMath::Power((b-U1S_t2),-n);
  }

  if(U2S_t1 > -alpha){
    Ups2S_1 = exp(-U2S_t1*U2S_t1/2.);
  }
  else if(U2S_t1 <= -alpha){
    Ups2S_1 = a*TMath::Power((b-U2S_t1),-n);
  }

  if(U2S_t2 > -alpha){
    Ups2S_2 = exp(-U2S_t2*U2S_t2/2.);
  }
  else if(U2S_t2 <= -alpha){
    Ups2S_2 = a*TMath::Power((b-U2S_t2),-n);
  }

  if(U3S_t1 > -alpha){
    Ups3S_1 = exp(-U3S_t1*U3S_t1/2.);
  }
  else if(U3S_t1 <= -alpha){
    Ups3S_1 = a*TMath::Power((b-U3S_t1),-n);
  }

  if(U3S_t2 > -alpha){
    Ups3S_2 = exp(-U3S_t2*U3S_t2/2.);
  }
  else if(U3S_t2 <= -alpha){
    Ups3S_2 = a*TMath::Power((b-U3S_t2),-n);
  }
  
  Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
  Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
  Double_t fN1s_1 = 1./(sigma*(fC+fD));
  Double_t fN1s_2 = 1./(sigma1_2*(fC+fD));
  Double_t fN2s_1 = 1./(sigma2_1*(fC+fD));
  Double_t fN2s_2 = 1./(sigma2_2*(fC+fD));
  Double_t fN3s_1 = 1./(sigma3_1*(fC+fD));
  Double_t fN3s_2 = 1./(sigma3_2*(fC+fD));

  Double_t SigM = N1*(fN1s_1*frac*Ups1S_1 + fN1s_2*(1-frac)*Ups1S_2) + N2*(fN2s_1*frac*Ups2S_1 + fN2s_2*(1-frac)*Ups2S_2) + N3*(fN3s_1*frac*Ups3S_1 + fN3s_2*(1-frac)*Ups3S_2);
  
  return SigM;
}
//}}}


//totalvn pol2 bkg Upsilon 1S&2S&3S{{{
Double_t Totalvnpol2U1S2S3S(Double_t* x, Double_t* par)
{
	Double_t N1 = par[0];
	Double_t N2 = par[1];
	Double_t N3 = par[2];
	Double_t Nbkg = par[3];
	Double_t mean = par[4];
	Double_t sigma = par[5];
	Double_t alpha = par[6];
	Double_t n = par[7];
	Double_t ratio = par[8];
	Double_t frac = par[9];
	Double_t Bkgmean = par[10];
	Double_t Bkgsigma = par[11];
	Double_t Bkgp0 = par[12];
	Double_t c = par[13];
	Double_t c_2s = par[14];
	Double_t c_3s = par[15];
	Double_t c1 = par[16];
	Double_t c2 = par[17];
	Double_t c3 = par[18];
	Double_t mean2 = mean*pdgMass.Y2S/pdgMass.Y1S;
	Double_t mean3 = mean*pdgMass.Y3S/pdgMass.Y1S;
	Double_t sigma1S_2 = sigma*ratio;
	Double_t sigma2S_1 = sigma*pdgMass.Y2S/pdgMass.Y1S;
	Double_t sigma2S_2 = sigma*ratio*pdgMass.Y2S/pdgMass.Y1S;
	Double_t sigma3S_1 = sigma*pdgMass.Y3S/pdgMass.Y1S;
	Double_t sigma3S_2 = sigma*ratio*pdgMass.Y3S/pdgMass.Y1S;
	
  //t2 > t1
	Double_t U1S_t1 = (x[0]-mean)/sigma;
	Double_t U1S_t2 = (x[0]-mean)/sigma1S_2;
	Double_t U2S_t1 = (x[0]-mean2)/sigma2S_1;
	Double_t U2S_t2 = (x[0]-mean2)/sigma2S_2;
	Double_t U3S_t1 = (x[0]-mean3)/sigma3S_1;
	Double_t U3S_t2 = (x[0]-mean3)/sigma3S_2;
	if (alpha < 0)
	{

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
	}

	Double_t absAlpha = fabs((Double_t)alpha);
	Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
	Double_t b = n/absAlpha - absAlpha;

  Double_t Ups1S_1 = -1;
  Double_t Ups1S_2 = -1;
  Double_t Ups2S_1 = -1;
  Double_t Ups2S_2 = -1;
  Double_t Ups3S_1 = -1;
  Double_t Ups3S_2 = -1;


  if(U1S_t1 > -alpha){
    Ups1S_1 = exp(-U1S_t1*U1S_t1/2.);
  }
  else if(U1S_t1 <= -alpha){
    Ups1S_1 = a*TMath::Power((b-U1S_t1),-n);
  }

  if(U1S_t2 > -alpha){
    Ups1S_2 = exp(-U1S_t2*U1S_t2/2.);
  }
  else if(U1S_t2 <= -alpha){
    Ups1S_2 = a*TMath::Power((b-U1S_t2),-n);
  }

  if(U2S_t1 > -alpha){
    Ups2S_1 = exp(-U2S_t1*U2S_t1/2.);
  }
  else if(U2S_t1 <= -alpha){
    Ups2S_1 = a*TMath::Power((b-U2S_t1),-n);
  }

  if(U2S_t2 > -alpha){
    Ups2S_2 = exp(-U2S_t2*U2S_t2/2.);
  }
  else if(U2S_t2 <= -alpha){
    Ups2S_2 = a*TMath::Power((b-U2S_t2),-n);
  }

  if(U3S_t1 > -alpha){
    Ups3S_1 = exp(-U3S_t1*U3S_t1/2.);
  }
  else if(U3S_t1 <= -alpha){
    Ups3S_1 = a*TMath::Power((b-U3S_t1),-n);
  }

  if(U3S_t2 > -alpha){
    Ups3S_2 = exp(-U3S_t2*U3S_t2/2.);
  }
  else if(U3S_t2 <= -alpha){
    Ups3S_2 = a*TMath::Power((b-U3S_t2),-n);
  }


/*  if( Ups1S_1 == -1 || Ups1S_2 == -1 || Ups2S_1 == -1 || Ups2S_2 == -1 || Ups3S_1 == -1 || Ups3S_2 == -1){
    cout << "ERROR ::: no match for UpsNS value " << endl;
    return -1;
  }
*/
  
  Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
  Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
  Double_t fN1s_1 = 1./(sigma*(fC+fD));
  Double_t fN1s_2 = 1./(sigma1S_2*(fC+fD));
  Double_t fN2s_1 = 1./(sigma2S_1*(fC+fD));
  Double_t fN2s_2 = 1./(sigma2S_2*(fC+fD));
  Double_t fN3s_1 = 1./(sigma3S_1*(fC+fD));
  Double_t fN3s_2 = 1./(sigma3S_2*(fC+fD));

  Double_t SigM = N1*(fN1s_1*frac*Ups1S_1 + fN1s_2*(1-frac)*Ups1S_2) + N2*(fN2s_1*frac*Ups2S_1 + fN2s_2*(1-frac)*Ups2S_2) + N3*(fN3s_1*frac*Ups3S_1 + fN3s_2*(1-frac)*Ups3S_2);
  Double_t SigM1s = N1*(fN1s_1*frac*Ups1S_1 + fN1s_2*(1-frac)*Ups1S_2);
  Double_t SigM2s = N2*(fN2s_1*frac*Ups2S_1 + fN2s_2*(1-frac)*Ups2S_2);
  Double_t SigM3s = N3*(fN3s_1*frac*Ups3S_1 + fN3s_2*(1-frac)*Ups3S_2);
  Double_t BkgM = Nbkg*(TMath::Exp(-x[0]/Bkgp0)*(TMath::Erf((x[0]-Bkgmean)/(TMath::Sqrt(2)*Bkgsigma))+1)/2.);
  
  return c*(SigM1s/(SigM+BkgM)) + c_2s*(SigM2s/(SigM+BkgM)) + c_3s*(SigM3s/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c3 + c2*x[0] + c1*x[0]*x[0]);
  
}
//}}}

//totalvn pol3 bkg Upsilon 1S&2S&3S{{{
Double_t Totalvnpol3U1S2S3S(Double_t* x, Double_t* par)
{
	Double_t N1 = par[0];
	Double_t N2 = par[1];
	Double_t N3 = par[2];
	Double_t Nbkg = par[3];
	Double_t mean = par[4];
	Double_t sigma = par[5];
	Double_t alpha = par[6];
	Double_t n = par[7];
	Double_t ratio = par[8];
	Double_t frac = par[9];
	Double_t Bkgmean = par[10];
	Double_t Bkgsigma = par[11];
	Double_t Bkgp0 = par[12];
	Double_t c = par[13];
	Double_t c_2s = par[14];
	Double_t c_3s = par[15];
	Double_t c1 = par[16];
	Double_t c2 = par[17];
	Double_t c3 = par[18];
	Double_t c4 = par[19];
	Double_t mean2 = mean*pdgMass.Y2S/pdgMass.Y1S;
	Double_t mean3 = mean*pdgMass.Y3S/pdgMass.Y1S;
	Double_t sigma1S_2 = sigma*ratio;
	Double_t sigma2S_1 = sigma*pdgMass.Y2S/pdgMass.Y1S;
	Double_t sigma2S_2 = sigma*ratio*pdgMass.Y2S/pdgMass.Y1S;
	Double_t sigma3S_1 = sigma*pdgMass.Y3S/pdgMass.Y1S;
	Double_t sigma3S_2 = sigma*ratio*pdgMass.Y3S/pdgMass.Y1S;
	
  //t2 > t1
	Double_t U1S_t1 = (x[0]-mean)/sigma;
	Double_t U1S_t2 = (x[0]-mean)/sigma1S_2;
	Double_t U2S_t1 = (x[0]-mean2)/sigma2S_1;
	Double_t U2S_t2 = (x[0]-mean2)/sigma2S_2;
	Double_t U3S_t1 = (x[0]-mean3)/sigma3S_1;
	Double_t U3S_t2 = (x[0]-mean3)/sigma3S_2;
	if (alpha < 0)
	{

    cout << "ERROR ::: alpha variable negative!!!! " << endl;
    return -1;
	}

	Double_t absAlpha = fabs((Double_t)alpha);
	Double_t a = TMath::Power(n/absAlpha,n)*exp(-absAlpha*absAlpha/2.);
	Double_t b = n/absAlpha - absAlpha;

  Double_t Ups1S_1 = -1;
  Double_t Ups1S_2 = -1;
  Double_t Ups2S_1 = -1;
  Double_t Ups2S_2 = -1;
  Double_t Ups3S_1 = -1;
  Double_t Ups3S_2 = -1;


  if(U1S_t1 > -alpha){
    Ups1S_1 = exp(-U1S_t1*U1S_t1/2.);
  }
  else if(U1S_t1 <= -alpha){
    Ups1S_1 = a*TMath::Power((b-U1S_t1),-n);
  }

  if(U1S_t2 > -alpha){
    Ups1S_2 = exp(-U1S_t2*U1S_t2/2.);
  }
  else if(U1S_t2 <= -alpha){
    Ups1S_2 = a*TMath::Power((b-U1S_t2),-n);
  }

  if(U2S_t1 > -alpha){
    Ups2S_1 = exp(-U2S_t1*U2S_t1/2.);
  }
  else if(U2S_t1 <= -alpha){
    Ups2S_1 = a*TMath::Power((b-U2S_t1),-n);
  }

  if(U2S_t2 > -alpha){
    Ups2S_2 = exp(-U2S_t2*U2S_t2/2.);
  }
  else if(U2S_t2 <= -alpha){
    Ups2S_2 = a*TMath::Power((b-U2S_t2),-n);
  }

  if(U3S_t1 > -alpha){
    Ups3S_1 = exp(-U3S_t1*U3S_t1/2.);
  }
  else if(U3S_t1 <= -alpha){
    Ups3S_1 = a*TMath::Power((b-U3S_t1),-n);
  }

  if(U3S_t2 > -alpha){
    Ups3S_2 = exp(-U3S_t2*U3S_t2/2.);
  }
  else if(U3S_t2 <= -alpha){
    Ups3S_2 = a*TMath::Power((b-U3S_t2),-n);
  }


/*  if( Ups1S_1 == -1 || Ups1S_2 == -1 || Ups2S_1 == -1 || Ups2S_2 == -1 || Ups3S_1 == -1 || Ups3S_2 == -1){
    cout << "ERROR ::: no match for UpsNS value " << endl;
    return -1;
  }
*/
  
  Double_t fC = n/absAlpha*1/(n-1)*exp(-absAlpha*absAlpha/2);
  Double_t fD = TMath::Sqrt(TMath::Pi()/2)*(1+TMath::Erf(absAlpha/TMath::Sqrt(2)));
  Double_t fN1s_1 = 1./(sigma*(fC+fD));
  Double_t fN1s_2 = 1./(sigma1S_2*(fC+fD));
  Double_t fN2s_1 = 1./(sigma2S_1*(fC+fD));
  Double_t fN2s_2 = 1./(sigma2S_2*(fC+fD));
  Double_t fN3s_1 = 1./(sigma3S_1*(fC+fD));
  Double_t fN3s_2 = 1./(sigma3S_2*(fC+fD));

  Double_t SigM = N1*(fN1s_1*frac*Ups1S_1 + fN1s_2*(1-frac)*Ups1S_2) + N2*(fN2s_1*frac*Ups2S_1 + fN2s_2*(1-frac)*Ups2S_2) + N3*(fN3s_1*frac*Ups3S_1 + fN3s_2*(1-frac)*Ups3S_2);
  Double_t SigM1s = N1*(fN1s_1*frac*Ups1S_1 + fN1s_2*(1-frac)*Ups1S_2);
  Double_t SigM2s = N2*(fN2s_1*frac*Ups2S_1 + fN2s_2*(1-frac)*Ups2S_2);
  Double_t SigM3s = N3*(fN3s_1*frac*Ups3S_1 + fN3s_2*(1-frac)*Ups3S_2);
  Double_t BkgM = Nbkg*(TMath::Exp(-x[0]/Bkgp0)*(TMath::Erf((x[0]-Bkgmean)/(TMath::Sqrt(2)*Bkgsigma))+1)/2.);
  
  return c*(SigM1s/(SigM+BkgM)) + c_2s*(SigM2s/(SigM+BkgM)) + c_3s*(SigM3s/(SigM+BkgM)) + (1 - SigM/(SigM+BkgM))*(c4 + c3*x[0] + c2*x[0]*x[0] + c1*x[0]*x[0]*x[0]);
  
}
//}}}


//pol2 bkg{{{
Double_t pol2bkg(Double_t* x, Double_t* par)
{
	Double_t c1 = par[0];
	Double_t c2 = par[1];
	Double_t c3 = par[2];

	return c1*x[0]*x[0]+c2*x[0]+c3;
}
//}}}

//pol2 bkg{{{
Double_t pol3bkg(Double_t* x, Double_t* par)
{
	Double_t c1 = par[0];
	Double_t c2 = par[1];
	Double_t c3 = par[2];
	Double_t c4 = par[3];

	return c1*x[0]*x[0]*x[0]+c2*x[0]*x[0]+c3*x[0] + c4;
}
//}}}


void doSimultaneousV2MassFit_SysEff_AllState(int cLow = 20, int cHigh = 180,
                     float ptLow = 0, float ptHigh = 50,
                     float yLow = 0, float yHigh=2.4,
                     float SiMuPtCut = 3.5, float massLow = 8, float massHigh =14, bool dimusign=true, int ibkg_vn_sel = fpol2)
{
  setTDRStyle();
	gStyle->SetOptFit(0000);
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;

//Make directory{{{
	TString mainDIR = gSystem->ExpandPathName(gSystem->pwd());
	TString fDIR = mainDIR + "/FitResult";

	void * dirf = gSystem->OpenDirectory(fDIR.Data());
	if(dirf) gSystem->FreeDirectory(dirf);
	else gSystem->mkdir(fDIR.Data(), kTRUE);

//}}}
  
  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh) ;
  TString dimusignString;
  if(dimusign) dimusignString = "OS";
  else if(!dimusign) dimusignString = "SS";
  
  kineLabel = kineLabel + Form("_m%.0f-%.0f",massLow,massHigh) + "_" + dimusignString;

  TF1* fyieldtot;
  TF1* fvntot;
  Double_t v2;
  Double_t v2e;
  Double_t v2_2s;
  Double_t v2e_2s;
  Double_t v2_3s;
  Double_t v2e_3s;
  Double_t v2_bkg;
  
  int nParmV_;
  int nParBkg;
  if(ibkg_vn_sel == fpol2) {nParmV_ = 19; nParBkg = 3;} 
  else if(ibkg_vn_sel == fpol3) {nParmV_ = 20; nParBkg = 4;}
  else{
    cout << "ERROR ::: No Selection for v2 background function!!!!" << endl;
    return;
  }
  
  //Get yield distribution{{{
  TFile* rf = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/plots/MassV2HistEffSys/Ups_%s_AccW1_EffW1.root",kineLabel.Data()),"read");
  TH1D* h_v2_SplusB = (TH1D*) rf->Get("h_v2_SplusB");  
  TGraphAsymmErrors* g_mass = (TGraphAsymmErrors*) rf->Get("g_mass");  

  TFile *wf = new TFile(Form("%s/SimFitResult_SysEff_%s.root",fDIR.Data(),kineLabel.Data()),"recreate");

  //define function for simultaneous fitting{{{
  TF1* fyield_simul = new TF1("fyield_simul", TotalYield, massLow, massHigh, nParmM);
  TF1* fvn_simul;
  if(ibkg_vn_sel == fpol2) fvn_simul = new TF1("fvn_simul", Totalvnpol2U1S2S3S, massLow, massHigh, nParmV_);
  else if(ibkg_vn_sel == fpol3) fvn_simul = new TF1("fvn_simul", Totalvnpol3U1S2S3S, massLow, massHigh, nParmV_);
  //}}}

  //combine functions{{{
  fyield_simul->SetLineColor(2);
  fyield_simul->SetLineWidth(1);

  fvn_simul->SetLineColor(2);
  fvn_simul->SetLineWidth(1);

  ROOT::Math::WrappedMultiTF1 wmass(*fyield_simul, 1);
  ROOT::Math::WrappedMultiTF1 wvn(*fvn_simul, 1);

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange massrange;

  massrange.SetRange(massLow,massHigh);
  ROOT::Fit::BinData datamass(opt, massrange);
  ROOT::Fit::FillData(datamass, g_mass);

  ROOT::Fit::DataRange vnrange;
  vnrange.SetRange(massLow,massHigh);
  ROOT::Fit::BinData datavn(opt, vnrange);
  ROOT::Fit::FillData(datavn, h_v2_SplusB);

  ROOT::Fit::Chi2Function mass_chi2(datamass, wmass);
  ROOT::Fit::Chi2Function vn_chi2(datavn, wvn);

  GlobalChi2_width globalChi2(mass_chi2, vn_chi2);

  ROOT::Fit::Fitter fitter;
  //}}}


  //Read Signal File
  TFile *sfile = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/MassFit/AllParmFreeFit/fitRes/fitresults_upsilon_DoubleCB_pt%.1f-%.1f_y%.1f-%.1f_muPt3.5_centrality0-180.root",ptLow,ptHigh,yLow,yHigh),"read");
  sfile->cd();
  RooWorkspace *ws = (RooWorkspace*) sfile -> Get("workspace");


  //Get fitting parameter{{{
  Double_t N1_ = ws->var("nSig1s")->getVal();
  Double_t N2_ = ws->var("nSig2s")->getVal();
  Double_t N3_ = ws->var("nSig3s")->getVal();
  Double_t Nbkg_ = ws->var("nBkg")->getVal();
  Double_t mean_ = pdgMass.Y1S;
  Double_t sigma_ = ws->var("sigma1s_1")->getVal();
  Double_t alpha_ = ws->var("alpha1s_1")->getVal();
  Double_t n_ = ws->var("n1s_1")->getVal();
  Double_t ratio_ = ws->var("x1s")->getVal();
  Double_t frac_ = ws->var("f1s")->getVal();
  Double_t Bkgmean_ = ws->var("#mu")->getVal();
  Double_t Bkgsigma_ = ws->var("#sigma")->getVal();
  Double_t Bkgp0_ = ws->var("#lambda")->getVal();
  Double_t c_ = 0.031;
  Double_t c_2s_ = 0;
  Double_t c_3s_ = 0;
  Double_t c1_ = -0.13142884;
  Double_t c2_ = -0.0114097;
  Double_t c3_ = 0.021;
  Double_t c4_ = -0.0011964;
  //}}}

  cout << "sigma_ : " << sigma_ << endl;
  
  Double_t par0[nParmV];

  par0[0] = N1_;
  par0[1] = N2_;
  par0[2] = N3_;
  par0[3] = Nbkg_;
  par0[4] = mean_;
  par0[5] = sigma_;
  par0[6] = alpha_;
  par0[7] = n_;
  par0[8] = ratio_;
  par0[9] = frac_;
  par0[10] = Bkgmean_;
  par0[11] = Bkgsigma_;
  par0[12] = Bkgp0_;
  par0[13] = c_;
  par0[14] = c_2s_;//c_2s_;
  par0[15] = c_3s_;//c_2s_;
  par0[16] = c1_;
  par0[17] = c2_;
  par0[18] = c3_;
  par0[19] = c4_;
  //}}}

  //combined function condition{{{
  Double_t parLimitLow[nParmV]  = {     0,     -2,     -2,        0, mean_ -0.02, 0.00,  0.0,  0.0, 0, 0, par0[10]*0.5,  par0[11]*0.5, par0[12]*0.5, -0.3, -0.4, -0.4, -3, -3, -3, -2};
  Double_t parLimitHigh[nParmV] = {N1_*10, N2_*10, N3_*10, Nbkg_*10, mean_ +0.02,  0.5, 10.0, 10.0, 1, 1, par0[10]*2.5,  par0[11]*2.5, par0[12]*2.5,  0.3,  0.4,  0.4,  3,  3,  3,  2};

  fitter.Config().SetParamsSettings(nParmV_, par0);
  for(int ipar = 0; ipar<nParmV_; ipar++){
    fitter.Config().ParSettings(ipar).SetLimits(parLimitLow[ipar],parLimitHigh[ipar]);
  }

  //Fix Signal Parameter
  fitter.Config().ParSettings(5).Fix();
  fitter.Config().ParSettings(6).Fix();
  fitter.Config().ParSettings(7).Fix();
  fitter.Config().ParSettings(8).Fix();
  fitter.Config().ParSettings(9).Fix();
//  fitter.Config().ParSettings(14).Fix();
//  fitter.Config().ParSettings(15).Fix();
//  fitter.Config().ParSettings(14).Fix();
//  fitter.Config().ParSettings(10).Fix();
//  fitter.Config().ParSettings(11).Fix();
//  fitter.Config().ParSettings(12).Fix();

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2","Migrad");
  //}}}

  fitter.FitFCN(nParmV_, globalChi2, 0, datamass.Size()+datavn.Size(), true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  //Yield fitting result{{{
  fyield_simul->SetFitResult(result, iparmass);
  fyield_simul->SetRange(massrange().first, massrange().second);
  g_mass->GetListOfFunctions()->Add(fyield_simul);

  TF1* fyield_bkg = new TF1("fyield_bkg", "[0]*( TMath::Exp(-x/[3])*(TMath::Erf((x-[1])/(TMath::Sqrt(2)*[2]))+1)/2. )", massLow, massHigh);
  fyield_bkg->FixParameter(0, fyield_simul->GetParameter(3));
  fyield_bkg->FixParameter(1, fyield_simul->GetParameter(10));
  fyield_bkg->FixParameter(2, fyield_simul->GetParameter(11));
  fyield_bkg->FixParameter(3, fyield_simul->GetParameter(12));

  g_mass->GetListOfFunctions()->Add(fyield_bkg);
  //}}}

  //vn fitting result{{{
  fvn_simul->SetFitResult(result, iparvn);
  fvn_simul->SetRange(vnrange().first, vnrange().second);

  TF1* fyield_sig = new TF1("fyield_sig",TotalYieldSig, massLow, massHigh, 9);
  fyield_sig->FixParameter(0, fvn_simul->GetParameter(0));
  fyield_sig->FixParameter(1, fvn_simul->GetParameter(1));
  fyield_sig->FixParameter(2, fvn_simul->GetParameter(2));
  fyield_sig->FixParameter(3, fvn_simul->GetParameter(4));
  fyield_sig->FixParameter(4, fvn_simul->GetParameter(5));
  fyield_sig->FixParameter(5, fvn_simul->GetParameter(6));
  fyield_sig->FixParameter(6, fvn_simul->GetParameter(7));
  fyield_sig->FixParameter(7, fvn_simul->GetParameter(8));
  fyield_sig->FixParameter(8, fvn_simul->GetParameter(9));

  TF1* fvn_bkg;
  if(ibkg_vn_sel == fpol2){
    fvn_bkg = new TF1("fvn_bkg",pol2bkg, massLow, massHigh, nParBkg);
    fvn_bkg->FixParameter(0, fvn_simul->GetParameter(16));
    fvn_bkg->FixParameter(1, fvn_simul->GetParameter(17));
    fvn_bkg->FixParameter(2, fvn_simul->GetParameter(18));
  }
  else if(ibkg_vn_sel == fpol3){
    fvn_bkg = new TF1("fvn_bkg",pol3bkg, massLow, massHigh, nParBkg);
    fvn_bkg->FixParameter(0, fvn_simul->GetParameter(16));
    fvn_bkg->FixParameter(1, fvn_simul->GetParameter(17));
    fvn_bkg->FixParameter(2, fvn_simul->GetParameter(18));
    fvn_bkg->FixParameter(3, fvn_simul->GetParameter(19));
  }

  unsigned int nfpxl = 1100;
//  fvn_bkg->SetNpx(nfpxl);
  fvn_simul->SetNpx(nfpxl);
//  fyield_bkg->SetNpx(nfpxl);
  fyield_simul->SetNpx(nfpxl);

  h_v2_SplusB->GetListOfFunctions()->Add(fvn_simul);
  h_v2_SplusB->GetListOfFunctions()->Add(fvn_bkg);

  v2 = fvn_simul->GetParameter(13);
  v2e = fvn_simul->GetParError(13);
  v2_2s = fvn_simul->GetParameter(14);
  v2e_2s = fvn_simul->GetParError(14);
  v2_3s = fvn_simul->GetParameter(15);
  v2e_3s = fvn_simul->GetParError(15);

  // Drawing 
  fyield_bkg->SetLineColor(kBlue);
  fyield_bkg->SetLineStyle(kDashed);
  fyield_simul->SetLineColor(kBlue);
  fyield_simul->SetLineWidth(2);

  fvn_simul->SetLineColor(kPink-2);
  fvn_simul->SetLineWidth(2);
  fvn_bkg->SetLineColor(kMagenta-2);
  fvn_bkg->SetLineWidth(1);
  fvn_bkg->SetLineStyle(kDashed);

  h_v2_SplusB->GetYaxis()->SetRangeUser(-0.17,0.17);
  h_v2_SplusB->GetYaxis()->SetTitle("v_{2}^{S+B}");
  h_v2_SplusB->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
  h_v2_SplusB->GetYaxis()->SetLabelSize(0.055);
  h_v2_SplusB->GetXaxis()->SetLabelSize(0.055);
  h_v2_SplusB->GetXaxis()->SetTitleSize(0.07);
  h_v2_SplusB->GetYaxis()->SetTitleSize(0.07);
  h_v2_SplusB->GetYaxis()->SetTitleOffset(1.2);
  h_v2_SplusB->GetXaxis()->SetTitleOffset(1.0);
  h_v2_SplusB->GetXaxis()->SetLabelOffset(0.011);
  h_v2_SplusB->GetXaxis()->SetRangeUser(massLow,massHigh);
  h_v2_SplusB->GetXaxis()->SetLimits(massLow,massHigh);
  SetHistStyle(h_v2_SplusB,0,0);
  SetGraphStyle2(g_mass,0,0);
  g_mass->SetMinimum(0);
  g_mass->GetXaxis()->SetLimits(massLow,massHigh);
  g_mass->GetXaxis()->SetRangeUser(massLow,massHigh);
  g_mass->GetYaxis()->SetTitleOffset(1.7);
  g_mass->GetYaxis()->SetTitle("Counts");
  g_mass->GetYaxis()->SetLabelSize(0.055);
  g_mass->GetXaxis()->SetLabelSize(0.055);
  g_mass->GetXaxis()->SetTitleSize(0.07);
  g_mass->GetYaxis()->SetTitleSize(0.07);
  g_mass->GetYaxis()->SetTitleOffset(1.2);
  g_mass->GetXaxis()->SetNdivisions(510);
  
  double sizeTick = 12;

  float pos_x = 0.23;
  float pos_x_mass = 0.53;
  float pos_y = 0.82;
  float pos_y_diff = 0.071;
  int text_color = 1;
  float text_size = 16;
  TString perc = "%";

  TCanvas* c_mass_v2 = new TCanvas("c_mass_v2","",590,750);
  TPad* pad1 = new TPad("pad1","pad1",0,0.5,1.0,1.0);
  TPad* pad2 = new TPad("pad2","pad2",0,0.0,1.0,0.5);
  c_mass_v2->cd();
  pad1->SetTicks(1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.19);
  pad1->SetTopMargin(0.08);
  pad1->Draw();
  pad1->cd();
  double pad1W = pad1->GetWw()*pad1->GetAbsWNDC();
  double pad1H = pad1->GetWh()*pad1->GetAbsHNDC();
  double tickScaleX = (pad1->GetUxmax() - pad1->GetUxmin())/(pad1->GetX2()-pad1->GetX1())*pad1H;
  g_mass->GetXaxis()->SetTickLength(sizeTick/tickScaleX);   
  g_mass->GetXaxis()->SetTitleSize(0);
  g_mass->GetXaxis()->SetLabelSize(0);
  g_mass->GetXaxis()->SetNdivisions(510);
  g_mass->Draw("AP");
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  else if(ptLow!=0) drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow, ptHigh ),pos_x_mass,pos_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow, yHigh ), pos_x_mass,pos_y-pos_y_diff,text_color,text_size);
  drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_x_mass,pos_y-pos_y_diff*2,text_color,text_size);
  drawText("|#eta^{#mu}| < 2.4", pos_x_mass,pos_y-pos_y_diff*3,text_color,text_size);
  drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_x_mass,pos_y-pos_y_diff*4,text_color,text_size);

  pad1->Update();

  c_mass_v2->cd();
  pad2->SetTicks(1,1);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.17);
  pad2->SetLeftMargin(0.19);
  pad2->Draw();
  pad2->cd();
  double pad2W = pad2->GetWw()*pad2->GetAbsWNDC();
  double pad2H = pad2->GetWh()*pad2->GetAbsHNDC();
/*
  TGraphErrors* g_v2_SplusB = new TGraphErrors(h_v2_SplusB);
  g_v2_SplusB->GetYaxis()->SetRangeUser(-0.17,0.17);
  g_v2_SplusB->GetYaxis()->SetTitle("v_{2}^{S+B}");
  g_v2_SplusB->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
  g_v2_SplusB->GetYaxis()->SetLabelSize(0.055);
  g_v2_SplusB->GetXaxis()->SetLabelSize(0.055);
  g_v2_SplusB->GetXaxis()->SetTitleSize(0.07);
  g_v2_SplusB->GetYaxis()->SetTitleSize(0.07);
  g_v2_SplusB->GetYaxis()->SetTitleOffset(1.2);
  g_v2_SplusB->GetXaxis()->SetTitleOffset(1.0);
  g_v2_SplusB->GetXaxis()->SetLabelOffset(0.011);
  g_v2_SplusB->GetXaxis()->SetLimits(massLow,massHigh);
  g_v2_SplusB->GetXaxis()->SetRangeUser(massLow,massHigh);
  */
  h_v2_SplusB->GetXaxis()->SetNdivisions(510);
  tickScaleX = (pad2->GetUxmax() - pad2->GetUxmin())/(pad2->GetX2()-pad2->GetX1())*pad2H;
  h_v2_SplusB->Draw("P");
  pad1->Update();
  pad2->Update();
  h_v2_SplusB->GetXaxis()->SetTickLength(sizeTick/tickScaleX);
  h_v2_SplusB->Draw("P");
  jumSun(massLow,0,massHigh,0,1,1);
  drawText(Form("v_{2}(#Upsilon(1S)) = %.3f #pm %.3f",v2,v2e),pos_x_mass,pos_y,text_color,text_size);
  drawText(Form("v_{2}(#Upsilon(2S)) = %.3f #pm %.3f",v2_2s,v2e_2s),pos_x_mass,pos_y-pos_y/10,text_color,text_size);
  drawText(Form("v_{2}(#Upsilon(3S)) = %.3f #pm %.3f",v2_3s,v2e_3s),pos_x_mass,pos_y-pos_y/10*2,text_color,text_size);
  CMS_lumi_v2mass(pad1,iPeriod,iPos, 1);  
  pad1->Update();
  pad2->Update();
  c_mass_v2->cd();
  pad1->Draw();
  pad2->Draw();
  c_mass_v2->SaveAs(Form("%s/v2Mass_SysEff_%s.pdf",fDIR.Data(),kineLabel.Data()));
/*
  Double_t xmass[200];
  Double_t pullmass[200];

  Float_t Chi2Yield = 0;
  Int_t NdfYield = g_mass->GetN() - fyield_simul->GetNumberFreeParameters();

  for(Int_t ibin = 1; ibin <= g_mass->GetN(); ibin++)
  {
    xmass[ibin] = h_mass->GetBinCenter(ibin);
    pullmass[ibin] = (h_mass->GetBinContent(ibin) - fyield_simul->Eval(xmass[ibin]))/h_mass->GetBinError(ibin);
    Chi2Yield += pullmass[ibin]*pullmass[ibin];
  }

  TLatex* lt1 = new TLatex();
  lt1->SetNDC();
  lt1->DrawLatex(0.49, 0.5, Form("Chi2/ndf = %.f/%d", Chi2Yield, NdfYield));

  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.17);
  pad2->SetLeftMargin(0.19);
  pad2->SetTicks();
  pad2->cd();
  jumSun(massLow,0,massHigh,0,1,1);
  drawText(Form("v_{2}(S) = %.3f #pm %.3f",v2,v2e),pos_x_mass,pos_y,text_color,text_size);

  h_v2_SplusB->Draw("P");
  
  Double_t xv2[200];
  Double_t pullv2[200];
  Double_t v2y[200];
  
  Float_t Chi2v2 = 0;
  Int_t Ndfv2 = h_v2_SplusB->GetNbinsX()-fvn_bkg->GetNumberFreeParameters();
  TGraphErrors *gvn = new TGraphErrors(h_v2_SplusB);
  for(Int_t ibin = 0; ibin < gvn->GetN(); ibin++)
  {
    gvn->GetPoint(ibin, xv2[ibin], v2y[ibin]);
    pullv2[ibin] = (v2y[ibin] - fvn_simul->Eval(xv2[ibin]))/gvn->GetErrorY(ibin);
    if(fabs(pullv2[ibin]) < 100)
    {
      Chi2v2 += pullv2[ibin]*pullv2[ibin];
    }
  }

  lt1->DrawLatex(0.5, 0.38, Form("Chi2/ndf = %.f/%d", Chi2v2, Ndfv2));
  //}}}
*/
/*
  //Pull Distribution 
  TCanvas* c_pull = new TCanvas("c_pull","",590,750);
  TPad* pad_pull1 = new TPad("pad_pull1","pad_pull1",0,0.5,1.0,1.0);
  TPad* pad_pull2 = new TPad("pad_pull2","pad_pull2",0,0.0,1.0,0.5);
  c_pull->cd();
  pad_pull1->SetTicks(1,1);
  pad_pull1->SetBottomMargin(0);
  pad_pull1->SetLeftMargin(0.19);
  pad_pull1->SetTopMargin(0.08);
  pad_pull1->Draw();
  pad_pull1->cd();

  TH1D* h_mass_pull = (TH1D*) h_mass->Clone("h_mass_pull"); h_mass->Reset();
  for(Int_t ibin = 1; ibin <= h_mass->GetNbinsX(); ibin++)
  {
    xmass[ibin] = h_mass->GetBinCenter(ibin);
    pullmass[ibin] = (h_mass->GetBinContent(ibin) - fyield_simul->Eval(xmass[ibin]))/h_mass->GetBinError(ibin);
    Chi2Yield += pullmass[ibin]*pullmass[ibin];
  }

  TLatex* lt1 = new TLatex();
  lt1->SetNDC();
  lt1->DrawLatex(0.49, 0.5, Form("Chi2/ndf = %.f/%d", Chi2Yield, NdfYield));

  
  Double_t xv2[200];
  Double_t pullv2[200];
  Double_t v2y[200];
  
  Float_t Chi2v2 = 0;
  Int_t Ndfv2 = h_v2_SplusB->GetNbinsX()-fvn_bkg->GetNumberFreeParameters();
  TGraphErrors *gvn = new TGraphErrors(h_v2_SplusB);
  for(Int_t ibin = 0; ibin < gvn->GetN(); ibin++)
  {
    gvn->GetPoint(ibin, xv2[ibin], v2y[ibin]);
    pullv2[ibin] = (v2y[ibin] - fvn_simul->Eval(xv2[ibin]))/gvn->GetErrorY(ibin);
    if(fabs(pullv2[ibin]) < 100)
    {
      Chi2v2 += pullv2[ibin]*pullv2[ibin];
    }
  }

  lt1->DrawLatex(0.5, 0.38, Form("Chi2/ndf = %.f/%d", Chi2v2, Ndfv2));
  //}}}
*/
  wf->cd();
  //store individual function{{{
  fyieldtot = (TF1*) fyield_simul->Clone();
  fyieldtot->SetName("massfit");
  fyieldtot->Write();

  fvntot = (TF1*) fvn_simul->Clone();
  fvntot->SetName("vnfit");
  fvntot->Write();

  fyield_bkg->Write();
  fyield_sig->Write();
  fvn_bkg->Write();
  //}}}

  //v2_bkg[ipt] = fvn_simul->GetParameter(14) + fvn_simul->GetParameter(14) * U1S_mass;
  //get Chi2{{{
  Double_t ptar = (ptHigh+ptLow)/2;
  TGraphErrors* v2plot = new TGraphErrors();
  v2plot->SetPoint(0,ptar,v2);
  v2plot->SetPointError(0,0,v2e);
  v2plot->SetPoint(1,ptar*2,v2_2s);
  v2plot->SetPointError(1,0,v2e_2s);
  v2plot->SetPoint(2,ptar*3,v2_3s);
  v2plot->SetPointError(2,0,v2e_3s);
  v2plot->SetTitle("");
  v2plot->SetName("v2vspt");
  v2plot->Write();
}
