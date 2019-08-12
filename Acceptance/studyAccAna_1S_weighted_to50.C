//#include "commonUtility.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
//#include "cutsAndBin.h"
using namespace std;
// Editor : Dongho Moon (dmoon@cern.ch)

// Define weighting functions for pp and AA (Data/MC)
TF1* fWgtPP1 = new TF1("fWgtPP1","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2])))");
TF1* fWgtPP2 = new TF1("fWgtPP2","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2])))");
TF1* fWgtAA1 = new TF1("fWgtAA1","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*9.460))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*9.460))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(9.460*9.460+x*x)-9.460)/([2]*[3])),-[2])))");
TF1* fWgtAA2 = new TF1("fWgtAA2","(([0]-1)*([0]-2)*([2]*[3]*([2]*[3]+([2]-2)*10.023))*TMath::Power((1+(TMath::Sqrt(10.023*10.023+x*x)-10.023)/([0]*[1])),-[0])/(([0]*[1]*([0]*[1]+([0]-2)*10.023))*(([2]-1)*([2]-2))*TMath::Power((1+(TMath::Sqrt(10.032*10.023+x*x)-10.023)/([2]*[3])),-[2])))");


// Define weighting functions for systematics
TF1 *fWgtHigh = new TF1("fWgtHigh","[0]",0.0,50.0);
TF1 *fWgtLow = new TF1("fWgtLow","[0]",0.0,50.0);

// Define up and down 20 % with pt
TF1 *fWgtPtp = new TF1("fWgtPtp","0.8 + 0.0133*x",0.0,50.0); // 20 % up
TF1 *fWgtPtm = new TF1("fWgtPtm","1.2 - 0.0133*x",0.0,50.0); // 20 % down
TF1 *fWgtYp = new TF1("fWgtYp","0.8 + 0.167*x",0.0,2.4); // 20 % up
TF1 *fWgtYm = new TF1("fWgtYm","1.2 - 0.167*x",0.0,2.4); // 20 % down

TF1* fdNdpTWgt = new TF1("fdNdpTWgt","([0]+[1]*x+[2]*x*x)/((x-[3])*(x-[3])*(x-[3]))");


TH1F *hPadPt = new TH1F("hPadPt",";p_{T} GeV/c;Acceptance",10,0.0,50.0);
TH1F *hPadPt1 = new TH1F("hPadPt1",";p_{T} GeV/c;Acceptance",10,0.0,50.0);
TH1F *hPadY = new TH1F("hPadY",";|y|;Acceptance",10,0.0,2.4);
TH1F *hPadY1 = new TH1F("hPadY1",";|y|;Acceptance",10,0.0,2.4);

TLatex *lt1 = new TLatex();

void studyAccAna_1S_weighted_to50(int Cat_ = 0, int Wgtopt=1) { // Cat_ == 0 (1S),  Cat_ == 1 (2S), Cat_ == 2 (3S)
//  gROOT->Macro("~/rootlogon.C");
  gStyle->SetOptStat(0);

  cout<<" Calculating acceptance starts !!!"<<endl;
  if(Cat_ == 0) cout<<" %%%%% Upsilon 1S %%%%% "<<endl;
  if(Cat_ == 1) cout<<" %%%%% Upsilon 2S %%%%% "<<endl;
  if(Cat_ == 2) cout<<" %%%%% Upsilon 3S %%%%% "<<endl;

  hPadPt->GetYaxis()->CenterTitle();
  hPadPt->GetYaxis()->SetTitleOffset(1.3);
  hPadPt->GetXaxis()->CenterTitle();
  hPadY->GetYaxis()->CenterTitle();
  hPadY->GetYaxis()->SetTitleOffset(1.3);
  hPadY->GetXaxis()->CenterTitle();

  hPadPt1->GetYaxis()->CenterTitle();
  hPadPt1->GetYaxis()->SetTitleOffset(1.3);
  hPadPt1->GetXaxis()->CenterTitle();
  hPadY1->GetYaxis()->CenterTitle();
  hPadY1->GetYaxis()->SetTitleOffset(1.3);
  hPadY1->GetXaxis()->CenterTitle();


  lt1->SetNDC();

  fWgtPP1->SetParameters( 0.988141, 3.0971, 1.81891, 10.0239);
  fWgtPP2->SetParameters( 11.518, 7.53196, 2.38444, 2.68481);
  fWgtAA1->SetParameters( 1.0001, 5.1, 2.0024, 12.4243);
  fWgtAA2->SetParameters( 3.46994, 11.8612, 2.10006, 3.25859);

  TH1::SetDefaultSumw2();
  //// modify by hand according to the pt range of the sample
  int nPtBins = 1;  double* ptBin = NULL;  int nYBins = 1;  double *yBin = NULL;
  const int nPtBins1 = 1;   double ptBin1[nPtBins1+1] = {0,50};             const int nYBins1 = 1;   double yBin1[nYBins1+1] = {0,2.4}; // integrated
  const int nPtBins2 = 6;   double ptBin2[nPtBins2+1] = {0,2,4,6,9,12,50};  const int nYBins2 = 1;   double yBin2[nYBins2+1] = {0,2.4}; // 1S
  const int nPtBins3 = 1;   double ptBin3[nPtBins3+1] = {0,50};             const int nYBins3 = 6;   double yBin3[nYBins3+1] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4}; // 1S
  const int nPtBins4 = 3;   double ptBin4[nPtBins4+1] = {0,4,9,50};         const int nYBins4 = 1;   double yBin4[nYBins4+1] = {0,2.4}; // 2S
  const int nPtBins5 = 1;   double ptBin5[nPtBins5+1] = {0,50};             const int nYBins5 = 3;   double yBin5[nYBins5+1] = {0,0.8,1.6,2.4}; // 2S
  const int nPtBins6 = 2;   double ptBin6[nPtBins6+1] = {0,6,50};           const int nYBins6 = 1;   double yBin6[nYBins6+1] = {0,2.4}; // 3S
  const int nPtBins7 = 1;   double ptBin7[nPtBins7+1] = {0,50};             const int nYBins7 = 2;   double yBin7[nYBins7+1] = {0,1.2,2.4}; // 3S

  // bin 1 : integrated
  // bin 2, 4 : pT 
  // bin 3, 5 : y

  TFile *in; // input skimed files
//  if(Cat_ == 0) in = new TFile("skimedForAcc_MC_Ups1S_20170808.root");//0606
  if(Cat_ == 0) in = new TFile("/afs/cern.ch/user/g/goni/public/Acceptance/skimedForAcc_MC_Ups1S_20170808.root");//0606
  if(Cat_ == 1) in = new TFile("skimedForAcc_MC_Ups2S_20170606.root");
  if(Cat_ == 2) in = new TFile("skimedForAcc_MC_Ups3S_20170606.root");
  
  TFile *out; // define output file
//  if(Cat_ == 0) out = new TFile("acceptance_wgt_1S_20180610.root","RECREATE");//20170111
//  if(Cat_ == 0) out = new TFile(Form("acceptance_wgt_1S_20190810_Wgtopt%d_0_to_50.root",Wgtopt),"RECREATE");//20170111
	if(Cat_ == 0) {
		if (wgtopt==0) out = new TFile("acceptance_wgt_1S_pt0_50_20190810_nonweighted.root","RECREATE");
	   if (wgtopt==1) out = new TFile("acceptance_wgt_1S_pt0_50_20190810_dNdptWeighted.root","RECREATE");
	}
  if(Cat_ == 1) out = new TFile("acceptance_wgt_2S_20170111.root","RECREATE");
  if(Cat_ == 2) out = new TFile("acceptance_wgt_3S_20170111.root","RECREATE");

  TFile *fweight;
  if(Cat_ == 0) fweight = new TFile("Func_dNdpT.root");

  TF1* fitRatio = (TF1*)fweight->Get("fitRatio");
  double fdNdpTWgt_p0 = fitRatio->GetParameter(0);
  double fdNdpTWgt_p1 = fitRatio->GetParameter(1);
  double fdNdpTWgt_p2 = fitRatio->GetParameter(2);
  double fdNdpTWgt_p3 = fitRatio->GetParameter(3);

	std::cout << fdNdpTWgt_p0 << std::endl;
	std::cout << fdNdpTWgt_p1 << std::endl;
	std::cout << fdNdpTWgt_p2 << std::endl;
	std::cout << fdNdpTWgt_p3 << std::endl;




  fdNdpTWgt->SetParameters( fdNdpTWgt_p0, fdNdpTWgt_p1, fdNdpTWgt_p2, fdNdpTWgt_p3);	
//TF1* fdNdpTWgt = new TF1("fdNdpTWgt","([0]+[1]*x+[2]*x*x)/((x-[3])*(x-[3])*(x-[3]))");

	std::cout << fdNdpTWgt->GetParameter(0) << std::endl;
	std::cout << fdNdpTWgt->GetParameter(1) << std::endl;
	std::cout << fdNdpTWgt->GetParameter(2) << std::endl;
	std::cout << fdNdpTWgt->GetParameter(3) << std::endl;





  // 1 : denominator, 2 : numerator 
  TH1F *hIntAccNoW = new TH1F("hIntAccNoW",";;",nYBins1,yBin1);
  TH1F *hIntAccNoW1 = new TH1F("hIntAccNoW1",";;",nYBins1,yBin1);
  TH1F *hIntAccNoW2 = new TH1F("hIntAccNoW2",";;",nYBins1,yBin1);
  TH1F *hIntAccPP = new TH1F("hIntAccPP",";;",nYBins1,yBin1);
  TH1F *hIntAccPP1 = new TH1F("hIntAccPP1",";;",nYBins1,yBin1);
  TH1F *hIntAccPP2 = new TH1F("hIntAccPP2",";;",nYBins1,yBin1);
  TH1F *hIntAccAA = new TH1F("hIntAccAA",";;",nYBins1,yBin1);
  TH1F *hIntAccAA1 = new TH1F("hIntAccAA1",";;",nYBins1,yBin1);
  TH1F *hIntAccAA2 = new TH1F("hIntAccAA2",";;",nYBins1,yBin1);

  hIntAccNoW->Sumw2();
  hIntAccNoW1->Sumw2();
  hIntAccNoW2->Sumw2();
  hIntAccPP->Sumw2();
  hIntAccPP1->Sumw2();
  hIntAccPP2->Sumw2();
  hIntAccAA->Sumw2();
  hIntAccAA1->Sumw2();
  hIntAccAA2->Sumw2();


  if(Cat_ == 0) {
    nPtBins = nPtBins2; ptBin = ptBin2;
    nYBins = nYBins3; yBin = yBin3;
  }else if(Cat_ == 1){
    nPtBins = nPtBins4; ptBin = ptBin4;
    nYBins = nYBins5; yBin = yBin5;
  }else if(Cat_ == 2){
    nPtBins = nPtBins6; ptBin = ptBin6;
    nYBins = nYBins7; yBin = yBin7;
  }
	
/*
  TH1F *hptAccNoW = new TH1F("hptAccNoW",";;",nPtBins,ptBin);
  TH1F *hptAccNoW1 = new TH1F("hptAccNoW1",";;",nPtBins,ptBin);
  TH1F *hptAccNoW2 = new TH1F("hptAccNoW2",";;",nPtBins,ptBin);
  TH1F *hptAccPP = new TH1F("hptAccPP",";;",nPtBins,ptBin);
  TH1F *hptAccPP1 = new TH1F("hptAccPP1",";;",nPtBins,ptBin);
  TH1F *hptAccPP2 = new TH1F("hptAccPP2",";;",nPtBins,ptBin);
  TH1F *hptAccAA = new TH1F("hptAccAA",";;",nPtBins,ptBin);
  TH1F *hptAccAA1 = new TH1F("hptAccAA1",";;",nPtBins,ptBin);
  TH1F *hptAccAA2 = new TH1F("hptAccAA2",";;",nPtBins,ptBin);
*/
  TH1F *hptAccNoW = new TH1F("hptAccNoW",";;",50,0.0,50.0);
  TH1F *hptAccNoW1 = new TH1F("hptAccNoW1",";;",50,0.0,50.0);
  TH1F *hptAccNoW2 = new TH1F("hptAccNoW2",";;",50,0.0,50.0);
  TH1F *hptAccPP = new TH1F("hptAccPP",";;",50,0.0,50.0);
  TH1F *hptAccPP1 = new TH1F("hptAccPP1",";;",50,0.0,50.0);
  TH1F *hptAccPP2 = new TH1F("hptAccPP2",";;",50,0.0,50.0);
  TH1F *hptAccAA = new TH1F("hptAccAA",";;",50,0.0,50.0);
  TH1F *hptAccAA1 = new TH1F("hptAccAA1",";;",50,0.0,50.0);
  TH1F *hptAccAA2 = new TH1F("hptAccAA2",";;",50,0.0,50.0);


  hptAccNoW->Sumw2();
  hptAccNoW1->Sumw2();
  hptAccNoW2->Sumw2();
  hptAccPP->Sumw2();
  hptAccPP1->Sumw2();
  hptAccPP2->Sumw2();
  hptAccAA->Sumw2();
  hptAccAA1->Sumw2();
  hptAccAA2->Sumw2();

  TH1F *hptAccPPpTp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPPpTp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPPpTm1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPPpTm2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPPyp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPPyp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPPym1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccPPym2 = (TH1F*)hptAccNoW->Clone();

  TH1F *hptAccAApTp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccAApTp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccAApTm1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccAApTm2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccAAyp1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccAAyp2 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccAAym1 = (TH1F*)hptAccNoW->Clone();
  TH1F *hptAccAAym2 = (TH1F*)hptAccNoW->Clone();


  TH1F *hrapAccNoW = new TH1F("hrapAccNoW",";;",nYBins,yBin);
  TH1F *hrapAccNoW1 = new TH1F("hrapAccNoW1",";;",nYBins,yBin);
  TH1F *hrapAccNoW2 = new TH1F("hrapAccNoW2",";;",nYBins,yBin);
  TH1F *hrapAccPP = new TH1F("hrapAccPP",";;",nYBins,yBin);
  TH1F *hrapAccPP1 = new TH1F("hrapAccPP1",";;",nYBins,yBin);
  TH1F *hrapAccPP2 = new TH1F("hrapAccPP2",";;",nYBins,yBin);
  TH1F *hrapAccAA = new TH1F("hrapAccAA",";;",nYBins,yBin);
  TH1F *hrapAccAA1 = new TH1F("hrapAccAA1",";;",nYBins,yBin);
  TH1F *hrapAccAA2 = new TH1F("hrapAccAA2",";;",nYBins,yBin);

  hrapAccNoW->Sumw2();
  hrapAccNoW1->Sumw2();
  hrapAccNoW2->Sumw2();
  hrapAccPP->Sumw2();
  hrapAccPP1->Sumw2();
  hrapAccPP2->Sumw2();
  hrapAccAA->Sumw2();
  hrapAccAA1->Sumw2();
  hrapAccAA2->Sumw2();

  TH1F *hrapAccPPpTp1 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccPPpTp2 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccPPpTm1 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccPPpTm2 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccPPyp1 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccPPyp2 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccPPym1 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccPPym2 = (TH1F*)hrapAccNoW->Clone();

  TH1F *hrapAccAApTp1 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccAApTp2 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccAApTm1 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccAApTm2 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccAAyp1 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccAAyp2 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccAAym1 = (TH1F*)hrapAccNoW->Clone();
  TH1F *hrapAccAAym2 = (TH1F*)hrapAccNoW->Clone();


  //getting trees from skimed files
  float mass = 0.0, pt = 0.0, phi = 0.0, y = 0.0, eta = 0.0;
  float pt1 = 0.0, phi1 = 0.0, eta1 = 0.0;
  float pt2 = 0.0, phi2 = 0.0, eta2 = 0.0;

  TBranch *b_mass;
  TBranch *b_y;
  TBranch *b_pt;
  TBranch *b_phi;
  TBranch *b_eta;
  TBranch *b_pt1;
  TBranch *b_phi1;
  TBranch *b_eta1;
  TBranch *b_pt2;
  TBranch *b_phi2;
  TBranch *b_eta2;

  TTree *mm = (TTree*)in->Get("mmGen");
  mm->SetBranchAddress("mass", &mass, &b_mass);
  mm->SetBranchAddress("pt", &pt, &b_pt);
  mm->SetBranchAddress("phi", &phi, &b_phi);
  mm->SetBranchAddress("y", &y, &b_y);
  mm->SetBranchAddress("eta", &eta, &b_eta);
  mm->SetBranchAddress("pt1", &pt1, &b_pt1);
  mm->SetBranchAddress("eta1", &eta1, &b_eta1);
  mm->SetBranchAddress("phi1", &phi1, &b_phi1);
  mm->SetBranchAddress("pt2", &pt2, &b_pt2);
  mm->SetBranchAddress("eta2", &eta2, &b_eta2);
  mm->SetBranchAddress("phi2", &phi2, &b_phi2);



  int nEntries = mm->GetEntries();
  //nEntries = 100;

  int cnt1 = 0; // just for counting for denominator
  int cnt2 = 0; // just for counting for numerator
  for(int i = 0; i < nEntries; i++){
	 
    mm->GetEntry(i);
    if (i%200000==0) cout << ">>>>> EVENT " << i << " / " << mm->GetEntries() <<  endl;
    //cout << ">>>>> EVENT " << i << " / " << mm->GetEntries() <<  endl;
    if(fabs(y) < 2.4&&pt<50.0){
      double ppwgt = fWgtPP1->Eval(pt); // apply weighting factor from functions above-defined for pp
      double aawgt = fWgtAA1->Eval(pt); // apply weighting factor from functions above-defined for AA
      double pTwgtp = fWgtPtp->Eval(pt);
      double pTwgtm = fWgtPtm->Eval(pt);
      double ywgtp = fWgtYp->Eval(fabs(y));
      double ywgtm = fWgtYm->Eval(fabs(y));
      double shigh = fWgtHigh->Eval(pt);
      double slow = fWgtLow->Eval(pt);

		double Wgt_temp;

		if (Wgtopt==0) Wgt_temp=1.0;
      else if (Wgtopt==1) {Wgt_temp=fdNdpTWgt->Eval(pt);}//std::cout << "new weight : " << Wgt_temp << std::endl;}
//  fitRatio->SetParameters( fdNdpTWgt_p0, fdNdpTWgt_p1, fdNdpTWgt_p2, fdNdpTWgt_p3);	
//TF1* fdNdpTWgt = new TF1("fdNdpTWgt","([0]+[1]*x+[2]*x*x)/((x-[3])*(x-[3])*(x-[3]))");

		ppwgt=ppwgt*Wgt_temp;
		aawgt=aawgt*Wgt_temp;



      hIntAccNoW1->Fill(fabs(y),Wgt_temp);
      hptAccNoW1->Fill(pt,Wgt_temp);
      hrapAccNoW1->Fill(fabs(y),Wgt_temp);
      hIntAccPP1->Fill(fabs(y),ppwgt);
      hptAccPP1->Fill(pt,ppwgt);
      hrapAccPP1->Fill(fabs(y),ppwgt);
      hIntAccAA1->Fill(fabs(y),aawgt);
      hptAccAA1->Fill(pt,aawgt);
      hrapAccAA1->Fill(fabs(y),aawgt);

      // for systematics with 20 % up and down
      hptAccPPpTp1->Fill(pt,ppwgt*pTwgtp);
      hptAccPPpTm1->Fill(pt,ppwgt*pTwgtm);
      hptAccPPyp1->Fill(pt,ppwgt*ywgtp);
      hptAccPPym1->Fill(pt,ppwgt*ywgtm);
      hrapAccPPpTp1->Fill(fabs(y),ppwgt*pTwgtp);
      hrapAccPPpTm1->Fill(fabs(y),ppwgt*pTwgtm);
      hrapAccPPyp1->Fill(fabs(y),ppwgt*ywgtp);
      hrapAccPPym1->Fill(fabs(y),ppwgt*ywgtm);

      hptAccAApTp1->Fill(pt,aawgt*pTwgtp);
      hptAccAApTm1->Fill(pt,aawgt*pTwgtm);
      hptAccAAyp1->Fill(pt,aawgt*ywgtp);
      hptAccAAym1->Fill(pt,aawgt*ywgtm);
      hrapAccAApTp1->Fill(fabs(y),aawgt*pTwgtp);
      hrapAccAApTm1->Fill(fabs(y),aawgt*pTwgtm);
      hrapAccAAyp1->Fill(fabs(y),aawgt*ywgtp);
      hrapAccAAym1->Fill(fabs(y),aawgt*ywgtm);
      cnt1++;
      if(pt1>3.5 && pt2>3.5 && fabs(eta1) < 2.4 && fabs(eta2) < 2.4){
        hIntAccNoW2->Fill(fabs(y),Wgt_temp);
        hptAccNoW2->Fill(pt,Wgt_temp);
        hrapAccNoW2->Fill(fabs(y),Wgt_temp);
        hIntAccPP2->Fill(fabs(y),ppwgt);
        hptAccPP2->Fill(pt,ppwgt);
        hrapAccPP2->Fill(fabs(y),ppwgt);
        hIntAccAA2->Fill(fabs(y),aawgt);
        hptAccAA2->Fill(pt,aawgt);
        hrapAccAA2->Fill(fabs(y),aawgt);

        // for systematics with 20 % up and down
        hptAccPPpTp2->Fill(pt,ppwgt*pTwgtp);
        hptAccPPpTm2->Fill(pt,ppwgt*pTwgtm);
        hptAccPPyp2->Fill(pt,ppwgt*ywgtp);
        hptAccPPym2->Fill(pt,ppwgt*ywgtm);
        hrapAccPPpTp2->Fill(fabs(y),ppwgt*pTwgtp);
        hrapAccPPpTm2->Fill(fabs(y),ppwgt*pTwgtm);
        hrapAccPPyp2->Fill(fabs(y),ppwgt*ywgtp);
        hrapAccPPym2->Fill(fabs(y),ppwgt*ywgtm);

        hptAccAApTp2->Fill(pt,aawgt*pTwgtp);
        hptAccAApTm2->Fill(pt,aawgt*pTwgtm);
        hptAccAAyp2->Fill(pt,aawgt*ywgtp);
        hptAccAAym2->Fill(pt,aawgt*ywgtm);
        hrapAccAApTp2->Fill(fabs(y),aawgt*pTwgtp);
        hrapAccAApTm2->Fill(fabs(y),aawgt*pTwgtm);
        hrapAccAAyp2->Fill(fabs(y),aawgt*ywgtp);
        hrapAccAAym2->Fill(fabs(y),aawgt*ywgtm);

        cnt2++;
      }
    }
  }
  cout<<"cnt1 : "<<cnt1<<", cnt2 : "<<cnt2<<endl;
  // Cloning for dividing
  TH1F * hIntAccNoW2Cp = (TH1F*)hIntAccNoW2->Clone();
  TH1F * hptAccNoW2Cp = (TH1F*)hptAccNoW2->Clone();
  TH1F * hrapAccNoW2Cp = (TH1F*)hrapAccNoW2->Clone();
  TH1F * hIntAccPP2Cp = (TH1F*)hIntAccPP2->Clone();
  TH1F * hptAccPP2Cp = (TH1F*)hptAccPP2->Clone();
  TH1F * hrapAccPP2Cp = (TH1F*)hrapAccPP2->Clone();
  TH1F * hIntAccAA2Cp = (TH1F*)hIntAccAA2->Clone();
  TH1F * hptAccAA2Cp = (TH1F*)hptAccAA2->Clone();
  TH1F * hrapAccAA2Cp = (TH1F*)hrapAccAA2->Clone();
  hIntAccNoW2Cp->Divide(hIntAccNoW1);
  hptAccNoW2Cp->Divide(hptAccNoW1);
  hrapAccNoW2Cp->Divide(hrapAccNoW1);
  hIntAccPP2Cp->Divide(hIntAccPP1);
  hptAccPP2Cp->Divide(hptAccPP1);
  hrapAccPP2Cp->Divide(hrapAccPP1);
  hIntAccAA2Cp->Divide(hIntAccAA1);
  hptAccAA2Cp->Divide(hptAccAA1);
  hrapAccAA2Cp->Divide(hrapAccAA1);

  // for systematics with 20 % up and down
  hptAccPPpTp2->Divide(hptAccPPpTp1);
  hptAccPPpTm2->Divide(hptAccPPpTm1);
  hptAccPPyp2->Divide(hptAccPPyp1);
  hptAccPPym2->Divide(hptAccPPym1);
  hrapAccPPpTp2->Divide(hrapAccPPpTp1);
  hrapAccPPpTm2->Divide(hrapAccPPpTm1);
  hrapAccPPyp2->Divide(hrapAccPPyp1);
  hrapAccPPym2->Divide(hrapAccPPym1);

  hptAccAApTp2->Divide(hptAccAApTp1);
  hptAccAApTm2->Divide(hptAccAApTm1);
  hptAccAAyp2->Divide(hptAccAAyp1);
  hptAccAAym2->Divide(hptAccAAym1);
  hrapAccAApTp2->Divide(hrapAccAApTp1);
  hrapAccAApTm2->Divide(hrapAccAApTm1);
  hrapAccAAyp2->Divide(hrapAccAAyp1);
  hrapAccAAym2->Divide(hrapAccAAym1);


  hptAccPPpTp2->SetName("hptAccPPSys1");
  hptAccPPpTm2->SetName("hptAccPPSys2");
  hptAccPPyp2->SetName("hptAccPPSys3");
  hptAccPPym2->SetName("hptAccPPSys4");
  hrapAccPPpTp2->SetName("hrapAccPPSys1");
  hrapAccPPpTm2->SetName("hrapAccPPSys2");
  hrapAccPPyp2->SetName("hrapAccPPSys3");
  hrapAccPPym2->SetName("hrapAccPPSys4");


  hptAccAApTp2->SetName("hptAccAASys1");
  hptAccAApTm2->SetName("hptAccAASys2");
  hptAccAAyp2->SetName("hptAccAASys3");
  hptAccAAym2->SetName("hptAccAASys4");
  hrapAccAApTp2->SetName("hrapAccAASys1");
  hrapAccAApTm2->SetName("hrapAccAASys2");
  hrapAccAAyp2->SetName("hrapAccAASys3");
  hrapAccAAym2->SetName("hrapAccAASys4");


  TCanvas *c1 = new TCanvas("c1","",1200,400);
  c1->Divide(3,1);
  c1->cd(1);
  hPadPt->Draw();
  hptAccNoW2Cp->Draw("same E");
  c1->cd(2);
  hPadPt->Draw();
  hptAccPP2Cp->Draw("same E");
  c1->cd(3);
  hPadPt->Draw();
  hptAccAA2Cp->Draw("same E");

  if(Cat_ == 0) {c1->SaveAs("plot_acc_1S_pt.png");c1->SaveAs("plot_acc_1S_pt.pdf");}
  if(Cat_ == 1) {c1->SaveAs("plot_acc_2S_pt.png");c1->SaveAs("plot_acc_2S_pt.pdf");}
  if(Cat_ == 2) {c1->SaveAs("plot_acc_3S_pt.png");c1->SaveAs("plot_acc_3S_pt.pdf");}

  c1->cd(1);
  hPadY->Draw();
  hrapAccNoW2Cp->Draw("same E");
  c1->cd(2);
  hPadY->Draw();
  hrapAccPP2Cp->Draw("same E");
  c1->cd(3);
  hPadY->Draw();
  hrapAccAA2Cp->Draw("same E");

  if(Cat_ == 0) {c1->SaveAs("plot_acc_1S_rap.png");c1->SaveAs("plot_acc_1S_rap.pdf");}
  if(Cat_ == 1) {c1->SaveAs("plot_acc_2S_rap.png");c1->SaveAs("plot_acc_2S_rap.pdf");}
  if(Cat_ == 2) {c1->SaveAs("plot_acc_3S_rap.png");c1->SaveAs("plot_acc_3S_rap.pdf");}
  
  TF1 *fline = new TF1("fline","1",0.0,50);
  fline->SetLineColor(kRed+1);
  fline->SetLineStyle(2);
  fline->SetLineWidth(2);

  // Ratio
  TH1F *hptAccNoWRat = (TH1F*)hptAccNoW2Cp->Clone();
  TH1F *hptAccPPRat = (TH1F*)hptAccPP2Cp->Clone();
  TH1F *hptAccAARat = (TH1F*)hptAccAA2Cp->Clone();
  hptAccPPRat->Divide(hptAccNoWRat);
  hptAccAARat->Divide(hptAccNoWRat);

  TH1F *hrapAccNoWRat = (TH1F*)hrapAccNoW2Cp->Clone();
  TH1F *hrapAccPPRat = (TH1F*)hrapAccPP2Cp->Clone();
  TH1F *hrapAccAARat = (TH1F*)hrapAccAA2Cp->Clone();
  hrapAccPPRat->Divide(hrapAccNoWRat);
  hrapAccAARat->Divide(hrapAccNoWRat);


  TCanvas *c2 = new TCanvas("c2","",800,400);
  c2->Divide(2,1);

  hPadPt1->GetYaxis()->SetRangeUser(0.8,1.2);
  hPadY1->GetYaxis()->SetRangeUser(0.8,1.2);

  c2->cd(1);
  hPadPt1->Draw();
  fline->Draw("same");
  hptAccPPRat->Draw("same");
  c2->cd(2);
  hPadPt1->Draw();
  fline->Draw("same");
  hptAccAARat->Draw("same");

  if(Cat_ == 0) {c2->SaveAs("plot_acc_ratio_1S_pt.png");c2->SaveAs("plot_acc_ratio_1S_pt.pdf");}
  if(Cat_ == 1) {c2->SaveAs("plot_acc_ratio_2S_pt.png");c2->SaveAs("plot_acc_ratio_2S_pt.pdf");}
  if(Cat_ == 2) {c2->SaveAs("plot_acc_ratio_3S_pt.png");c2->SaveAs("plot_acc_ratio_3S_pt.pdf");}

  c2->cd(1);
  hPadY1->Draw();
  fline->Draw("same");
  hrapAccPPRat->Draw("same");
  c2->cd(2);
  hPadY1->Draw();
  fline->Draw("same");
  hrapAccAARat->Draw("same");

  if(Cat_ == 0) {c2->SaveAs("plot_acc_ratio_1S_rap.png");c2->SaveAs("plot_acc_ratio_1S_rap.pdf");}
  if(Cat_ == 1) {c2->SaveAs("plot_acc_ratio_2S_rap.png");c2->SaveAs("plot_acc_ratio_2S_rap.pdf");}
  if(Cat_ == 2) {c2->SaveAs("plot_acc_ratio_3S_rap.png");c2->SaveAs("plot_acc_ratio_3S_rap.pdf");}

  if(Cat_ == 0) cout<<" %%%%% Upsilon 1S %%%%% "<<endl;
  if(Cat_ == 1) cout<<" %%%%% Upsilon 2S %%%%% "<<endl;
  if(Cat_ == 2) cout<<" %%%%% Upsilon 3S %%%%% "<<endl;

  // Acceptance values and systematics
  for(int i = 0; i < hIntAccNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Acceptance :: Integrated %%%"<<endl;
    cout<<"No Weighted "<<Form("%0.3f",hIntAccNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hIntAccNoW2Cp->GetBinError(i+1))<<endl;
    cout<<"pp "<<Form("%0.3f",hIntAccPP2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hIntAccPP2Cp->GetBinError(i+1))<<endl;
    cout<<"AA "<<Form("%0.3f",hIntAccAA2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hIntAccAA2Cp->GetBinError(i+1))<<endl;
  }

  cout<<""<<endl;
  for(int i = 0; i < hIntAccNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Systematics :: Integrated %%%"<<endl;
    cout<<"pp "<<Form("%0.3f",fabs(hIntAccPP2Cp->GetBinContent(i+1)-hIntAccNoW2Cp->GetBinContent(i+1))/hIntAccPP2Cp->GetBinContent(i+1)*100)<<endl;
    cout<<"AA "<<Form("%0.3f",fabs(hIntAccAA2Cp->GetBinContent(i+1)-hIntAccNoW2Cp->GetBinContent(i+1))/hIntAccAA2Cp->GetBinContent(i+1)*100)<<endl;
  }

  cout<<""<<endl;
  cout<<"%%% Pt  %%%"<<endl;
  for(int i = 0; i < hptAccNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"No Weighted :: "<<i<<" :: "<<Form("%0.3f",hptAccNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccNoW2Cp->GetBinError(i+1))<<endl;
    cout<<"pp :: "<<i<<" :: "<<Form("%0.3f",hptAccPP2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccPP2Cp->GetBinError(i+1))<<endl;
    cout<<"AA :: "<<i<<" :: "<<Form("%0.3f",hptAccAA2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hptAccAA2Cp->GetBinError(i+1))<<endl;
  }

  cout<<""<<endl;
  for(int i = 0; i < hptAccNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Systematics :: Pt %%%"<<endl;
    cout<<"pp "<<Form("%0.3f",fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccNoW2Cp->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100)<<endl;
    cout<<"AA "<<Form("%0.3f",fabs(hptAccAA2Cp->GetBinContent(i+1)-hptAccNoW2Cp->GetBinContent(i+1))/hptAccAA2Cp->GetBinContent(i+1)*100)<<endl;
  }

  cout<<""<<endl;
  cout<<"%%% Rapidity  %%%"<<endl;
  for(int i = 0; i < hrapAccNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"No Weighted :: "<<i<<" "<<Form("%0.3f",hrapAccNoW2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hrapAccNoW2Cp->GetBinError(i+1))<<endl;
    cout<<"pp :: "<<i<<" :: "<<Form("%0.3f",hrapAccPP2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hrapAccPP2Cp->GetBinError(i+1))<<endl;
    cout<<"AA :: "<<i<<" :: "<<Form("%0.3f",hrapAccAA2Cp->GetBinContent(i+1))<<" +- "<<Form("%0.4f",hrapAccAA2Cp->GetBinError(i+1))<<endl;
  }

  cout<<""<<endl;
  for(int i = 0; i < hrapAccNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Systematics :: Rapidity %%%"<<endl;
    cout<<"pp "<<Form("%0.3f",fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccNoW2Cp->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100)<<endl;
    cout<<"AA "<<Form("%0.3f",fabs(hrapAccAA2Cp->GetBinContent(i+1)-hrapAccNoW2Cp->GetBinContent(i+1))/hrapAccAA2Cp->GetBinContent(i+1)*100)<<endl;
  }

  // Calculation of Systematics for 20 % up and down

  cout<<""<<endl;
  cout<<"%%% Calculation of Systematics for 20 % up and down %%% "<<endl;
  cout<<""<<endl;
  cout<<"%%% Systematics :: Pt Bins %%%"<<endl;
  for(int i = 0; i < hptAccNoW2Cp->GetXaxis()->GetNbins(); i++){
    //cout<<"%%% Systematics :: Pt Bins %%%"<<endl;
    double sys1 = max(fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPPpTp2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100,fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPPpTm2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100);
    double sys2 = max(fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPPyp2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100,fabs(hptAccPP2Cp->GetBinContent(i+1)-hptAccPPym2->GetBinContent(i+1))/hptAccPP2Cp->GetBinContent(i+1)*100);
    double sys3 = max(fabs(hptAccAA2Cp->GetBinContent(i+1)-hptAccAApTp2->GetBinContent(i+1))/hptAccAA2Cp->GetBinContent(i+1)*100,fabs(hptAccAA2Cp->GetBinContent(i+1)-hptAccAApTm2->GetBinContent(i+1))/hptAccAA2Cp->GetBinContent(i+1)*100);
    double sys4 = max(fabs(hptAccAA2Cp->GetBinContent(i+1)-hptAccAAyp2->GetBinContent(i+1))/hptAccAA2Cp->GetBinContent(i+1)*100,fabs(hptAccAA2Cp->GetBinContent(i+1)-hptAccAAym2->GetBinContent(i+1))/hptAccAA2Cp->GetBinContent(i+1)*100);
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    cout<<"PP 1 : "<<Form("%0.3f",sys1)<<", PP 2 : "<<Form("%0.3f",sys2)<<", PP total : "<<Form("%0.3f",sqrt(sys2*sys2+sys1*sys1))<<endl;
    cout<<"AA 1 : "<<Form("%0.3f",sys3)<<", AA 2 : "<<Form("%0.3f",sys4)<<", AA total : "<<Form("%0.3f",sqrt(sys3*sys3+sys4*sys4))<<endl;
  }

  cout<<""<<endl;
  cout<<"%%% Systematics :: Rapidity Bins %%%"<<endl;
  for(int i = 0; i < hrapAccNoW2Cp->GetXaxis()->GetNbins(); i++){
    cout<<"%%% Bin "<<i+1<<" %%%"<<endl;
    double sys1 = max(fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPPpTp2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100,fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPPpTm2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100);
    double sys2 = max(fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPPyp2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100,fabs(hrapAccPP2Cp->GetBinContent(i+1)-hrapAccPPym2->GetBinContent(i+1))/hrapAccPP2Cp->GetBinContent(i+1)*100);
    double sys3 = max(fabs(hrapAccAA2Cp->GetBinContent(i+1)-hrapAccAApTp2->GetBinContent(i+1))/hrapAccAA2Cp->GetBinContent(i+1)*100,fabs(hrapAccAA2Cp->GetBinContent(i+1)-hrapAccAApTm2->GetBinContent(i+1))/hrapAccAA2Cp->GetBinContent(i+1)*100);
    double sys4 = max(fabs(hrapAccAA2Cp->GetBinContent(i+1)-hrapAccAAyp2->GetBinContent(i+1))/hrapAccAA2Cp->GetBinContent(i+1)*100,fabs(hrapAccAA2Cp->GetBinContent(i+1)-hrapAccAAym2->GetBinContent(i+1))/hrapAccAA2Cp->GetBinContent(i+1)*100);
    cout<<"PP 1 : "<<Form("%0.3f",sys1)<<", PP 2 : "<<Form("%0.3f",sys2)<<", PP total : "<<Form("%0.3f",sqrt(sys2*sys2+sys1*sys1))<<endl;
    cout<<"AA 1 : "<<Form("%0.3f",sys3)<<", AA 2 : "<<Form("%0.3f",sys4)<<", AA total : "<<Form("%0.3f",sqrt(sys3*sys3+sys4*sys4))<<endl;
  }

  out->cd();
/*
  hIntAccNoW1->Write();
  hptAccNoW1->Write();
  hrapAccNoW1->Write();
  hIntAccPP1->Write();
  hptAccPP1->Write();
  hrapAccPP1->Write();
  hIntAccAA1->Write();
  hptAccAA1->Write();
  hrapAccAA1->Write();

  hIntAccNoW2->Write();
  hptAccNoW2->Write();
  hrapAccNoW2->Write();
  hIntAccPP2->Write();
  hptAccPP2->Write();
  hrapAccPP2->Write();
  hIntAccAA2->Write();
  hptAccAA2->Write();
  hrapAccAA2->Write();
*/


  if(Cat_ == 0) {
    hIntAccNoW2Cp->SetName("hIntAccNoW1S");
    hptAccNoW2Cp->SetName("hptAccNoW1S");
    hrapAccNoW2Cp->SetName("hrapAccNoW1S");
    hIntAccPP2Cp->SetName("hIntAccPP1S");
    hptAccPP2Cp->SetName("hptAccPP1S");
    hrapAccPP2Cp->SetName("hrapAccPP1S");
    hIntAccAA2Cp->SetName("hIntAccAA1S");
    hptAccAA2Cp->SetName("hptAccAA1S");
    hrapAccAA2Cp->SetName("hrapAccAA1S");
  }

  if(Cat_ == 1) {
    hIntAccNoW2Cp->SetName("hIntAccNoW2S");
    hptAccNoW2Cp->SetName("hptAccNoW2S");
    hrapAccNoW2Cp->SetName("hrapAccNoW2S");
    hIntAccPP2Cp->SetName("hIntAccPP2S");
    hptAccPP2Cp->SetName("hptAccPP2S");
    hrapAccPP2Cp->SetName("hrapAccPP2S");
    hIntAccAA2Cp->SetName("hIntAccAA2S");
    hptAccAA2Cp->SetName("hptAccAA2S");
    hrapAccAA2Cp->SetName("hrapAccAA2S");
  }

  if(Cat_ == 2) {
    hIntAccNoW2Cp->SetName("hIntAccNoW3S");
    hptAccNoW2Cp->SetName("hptAccNoW3S");
    hrapAccNoW2Cp->SetName("hrapAccNoW3S");
    hIntAccPP2Cp->SetName("hIntAccPP3S");
    hptAccPP2Cp->SetName("hptAccPP3S");
    hrapAccPP2Cp->SetName("hrapAccPP3S");
    hIntAccAA2Cp->SetName("hIntAccAA3S");
    hptAccAA2Cp->SetName("hptAccAA3S");
    hrapAccAA2Cp->SetName("hrapAccAA3S");
  }
  hptAccNoW2Cp->Write();
  std::cout << "hptAccNoW2Cp->GetEntries() : " << hptAccNoW2Cp->GetEntries() << std::endl;
  std::cout << "hptAccNoW2Cp->Integral() : " << hptAccNoW2Cp->Integral() << std::endl;


/*
  hIntAccNoW2Cp->Write();
  hptAccNoW2Cp->Write();
  hrapAccNoW2Cp->Write();
  hIntAccPP2Cp->Write();
  hptAccPP2Cp->Write();
  hrapAccPP2Cp->Write();
  hIntAccAA2Cp->Write();
  hptAccAA2Cp->Write();
  hrapAccAA2Cp->Write();

  hptAccNoWRat->Write();
  hptAccPPRat->Write();
  hptAccAARat->Write();
  hrapAccNoWRat->Write();
  hrapAccPPRat->Write();
  hrapAccAARat->Write();

  hptAccPPpTp2->Write();
  hptAccPPpTm2->Write();
  hptAccPPyp2->Write();
  hptAccPPym2->Write();
  hrapAccPPpTp2->Write();
  hrapAccPPpTm2->Write();
  hrapAccPPyp2->Write();
  hrapAccPPym2->Write();

  hptAccAApTp2->Write();
  hptAccAApTm2->Write();
  hptAccAAyp2->Write();
  hptAccAAym2->Write();
  hrapAccAApTp2->Write();
  hrapAccAApTm2->Write();
  hrapAccAAyp2->Write();
  hrapAccAAym2->Write();
*/
//  out->Write();

}
