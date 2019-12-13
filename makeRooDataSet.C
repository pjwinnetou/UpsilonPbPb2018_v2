#include <iostream>
#include "commonUtility.h"
#include "cutsAndBinUpsilonV2.h"
#include "HiEvtPlaneList.h"
#include "Style_jaebeom.h"
#include "tdrstyle.C"
#include "CMS_lumi_v2mass.C"
#include "rootFitHeaders.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"

using namespace std;
using namespace RooFit;

using namespace hi;

double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);
void GetHistSqrt(TH1D* h1 =0, TH1D* h2=0);

void makeRooDataSet(bool isMC = true, bool fAccW = false, bool fEffW = false, int state=2)
{
  //Basic Setting
  gStyle->SetOptStat(0);
  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  TString dimusignString;

  //READ Input Skimmed File
  TFile *rf;
  if(isMC){
    if(state==1) rf = new TFile("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/skimmedFiles/OniaFlowSkim_UpsTrig_DB_isMC1_HFNom_Y1S_190801.root","read");
    else if(state==2) rf = new TFile("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/skimmedFiles/OniaFlowSkim_UpsTrig_DB_isMC1_HFNom_Y2S_190801.root","read");
  }
  else if(!isMC) rf = new TFile("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/skimmedFiles/OniaFlowSkim_UpsTrig_DBPD_isMC0_190710.root","read");
  TTree *tree = (TTree*) rf -> Get("mmepevt");

  
  //Get Correction histograms
  bool isTnP = true;
  TFile *fEff = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Efficiency/mc_eff_vs_pt_TnP%d_PtW1_OfficialMC_Y%dS_muPtCut3.5.root",isTnP,state),"read");
  TH1D* hEffPt[3];
  hEffPt[0] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_Cent010",isTnP)); 
  hEffPt[1] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_Cent1050",isTnP)); 
  hEffPt[2] = (TH1D*) fEff -> Get(Form("mc_eff_vs_pt_TnP%d_Cent5090",isTnP)); 
  TFile *fAcc = new TFile(Form("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Acceptance/acceptance_wgt_%dS_pt0_50_20190813_dNdptWeighted.root",state),"read");
  TH1D* hAccPt = (TH1D*) fAcc -> Get(Form("hptAccNoW%dS",state)); 


  //SetBranchAddress
  const int nMaxDimu = 1000;
  float mass[nMaxDimu];
  float qxa[nMaxDimu];
  float qya[nMaxDimu];
  float qxb[nMaxDimu]; 
  float qyb[nMaxDimu];
  float qxc[nMaxDimu];
  float qyc[nMaxDimu];
  float qxdimu[nMaxDimu];
  float qydimu[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu]; 
  float pt1[nMaxDimu];
  float pt2[nMaxDimu];
  float eta[nMaxDimu];
  float eta1[nMaxDimu];
  float eta2[nMaxDimu];
  Int_t cBin;
  Int_t event; 
  Int_t nDimu; 
  float vz;
  int recoQQsign[nMaxDimu];
  double weight;

  TBranch *b_event;
  TBranch *b_cBin;
  TBranch *b_nDimu;
  TBranch *b_vz;
  TBranch *b_mass;
  TBranch *b_recoQQsign;
  TBranch *b_pt;
  TBranch *b_y;
  TBranch *b_eta;
  TBranch *b_eta1;
  TBranch *b_eta2;
  TBranch *b_pt1;
  TBranch *b_pt2;
  TBranch *b_qxa;
  TBranch *b_qxb;
  TBranch *b_qxc;
  TBranch *b_qxdimu;
  TBranch *b_qya;
  TBranch *b_qyb;
  TBranch *b_qyc;
  TBranch *b_qydimu;
  TBranch *b_weight;
  

  tree -> SetBranchAddress("event", &event, &b_event);
  tree -> SetBranchAddress("cBin", &cBin, &b_cBin);
  tree -> SetBranchAddress("nDimu", &nDimu, &b_nDimu);
  tree -> SetBranchAddress("vz", &vz, &b_vz);
  tree -> SetBranchAddress("recoQQsign", recoQQsign, &b_recoQQsign);
  tree -> SetBranchAddress("mass", mass, &b_mass);
  tree -> SetBranchAddress("y", y, &b_y);
  tree -> SetBranchAddress("pt", pt, &b_pt);
  tree -> SetBranchAddress("pt1", pt1, &b_pt1);
  tree -> SetBranchAddress("pt2", pt2, &b_pt2);
  tree -> SetBranchAddress("eta", eta, &b_eta);
  tree -> SetBranchAddress("eta1", eta1, &b_eta1);
  tree -> SetBranchAddress("eta2", eta2, &b_eta2);
  tree -> SetBranchAddress("qxa", qxa, &b_qxa);
  tree -> SetBranchAddress("qxb", qxb, &b_qxb);
  tree -> SetBranchAddress("qxc", qxc, &b_qxc);
  tree -> SetBranchAddress("qxdimu", qxdimu, &b_qxdimu);
  tree -> SetBranchAddress("qya", qya, &b_qya);
  tree -> SetBranchAddress("qyb", qyb, &b_qyb);
  tree -> SetBranchAddress("qyc", qyc, &b_qyc);
  tree -> SetBranchAddress("qydimu", qydimu, &b_qydimu);
  tree -> SetBranchAddress("weight", &weight, &b_weight);
  
  ////////////////////////////////////////////////////////////////////////
  //////////////////  RooDataSet 
  ////////////////////////////////////////////////////////////////////////
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* evtWeight = new RooRealVar("weight","corr weight", 0, 10000,"");
  RooRealVar* recoQQ = new RooRealVar("recoQQsign","qq sign",-1,3,"");
//  RooRealVar* ctau3D = new RooRealVar("ctau3D","ctau3D dimuon",-20,20,"cm");
  RooRealVar* NumDimu = new RooRealVar("NumDimu","number of dimuon",0,100,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var,*evtWeight);
  argSet->add(*cBinVar); argSet->add(*recoQQ); argSet->add(*NumDimu); //argSet->add(*ctau3D);
  
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);


  int nDimuPass=0;
  int nDimu_one=0;
  int nDimu_more=0;
  int nDimu_all=0;

  double weight_acc = 1;
  double weight_eff = 1;
  

  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;
  
  double SiMuPtCut = 3.5;
  //Begin Loop
  for(int i=0; i<nEvt; i++){
    tree->GetEntry(i);
    nDimuPass=0;
    
    if(fabs(vz)>=15) continue;

    //Remove double candidate
    for(int j=0; j<nDimu; j++){
      if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && (double)pt1[j]>SiMuPtCut&&(double)pt2[j]>SiMuPtCut&&abs((double)eta1[j])<2.4&&abs((double)eta2[j])<2.4)) continue;
      nDimuPass++;
    }

    nDimu_all++;
    if(nDimuPass>1) {nDimu_more++; continue;}
    if(nDimuPass==1) nDimu_one++;

    // Fill Dimuon Loop
    for(int j=0; j<nDimu; j++){
        if(! ((double)pt[j]<50 && abs((double)y[j])<2.4 && (double)pt1[j]>SiMuPtCut&&(double)pt2[j]>SiMuPtCut&&abs((double)eta1[j])<2.4&&abs((double)eta2[j])<2.4)) continue;
        weight_acc = 1;
        weight_eff = 1;
        if(fAccW){weight_acc = getAccWeight(hAccPt, pt[j]);} 
        if(fEffW){ 
          if(cBin<20) weight_eff = getEffWeight(hEffPt[0], pt[j]);
          if(cBin>=20 && cBin<100) weight_eff = getEffWeight(hEffPt[1], pt[j]);
          if(cBin>=100 && cBin<180) weight_eff = getEffWeight(hEffPt[2], pt[j]);
        }
        double weight_ = weight * weight_eff * weight_acc;
        recoQQ->setVal((int)recoQQsign[j]);     
        massVar->setVal( (double)mass[j] ) ;
        ptVar->setVal(   (double)pt[j]   ) ;
        yVar->setVal(    (double)y[j]    ) ;
        pt1Var->setVal(  (double)pt1[j]  ) ;
        eta1Var->setVal( (double)eta1[j] ) ;
        pt2Var->setVal(  (double)pt2[j]  ) ;
        eta2Var->setVal( (double)eta2[j] ) ;
        cBinVar->setVal( (double)cBin ) ;
        evtWeight->setVal( (double)weight_ ) ;
        NumDimu->setVal((int)nDimu);
        dataSet->add( *argSet);
    }
  }

  cout << "All : " << nDimu_all << endl;
  cout << "more than one dimuon : " << nDimu_more << endl;
  cout << "one dimuon : " << nDimu_one << endl;
  
  TFile *wf = new TFile(Form("skimmedFiles/OniaRooDataSet_isMC%d_Y%dSW.root",isMC,state),"recreate");
  wf->cd();
 dataSet->Write();
}
    

void GetHistSqrt(TH1D* h1, TH1D* h2){
  if(h1->GetNbinsX() != h2->GetNbinsX()){ cout << "Inconsistent # of bins b/w histograms !! " << endl;}
  double content;
  double err;
  for(int i=1; i<=h1->GetNbinsX(); i++){
    content=0;err=0;
    content = h1->GetBinContent(i);
    err = h1->GetBinError(i);
    err = 0.5*err*TMath::Power(content,-0.5);
    h2->SetBinContent(i,TMath::Sqrt(content));
    h2->SetBinError(i,err);
  }
} 

double getAccWeight(TH1D* h, double pt){
  double binN = h->FindBin(pt);
  double weight_ = 1./(h->GetBinContent(binN));
  return weight_;
} 

double getEffWeight(TH1D *h, double pt){
  double binN = h->FindBin(pt);
  double weight_ = 1./(h->GetBinContent(binN));
  return weight_;
} 
