#include <iostream>

#include <TLorentzVector.h>
#include "../commonUtility.h"
#include "../HiEvtPlaneList.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../Style_jaebeom.h"
#include "tnp_weight_lowptPbPb.h"

using namespace std;

void getEfficiency_syst(
  float ptLow = 0.0, float ptHigh = 50.0,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 0, int cHigh = 180, bool isTnP = true, bool isPtWeight = true, int fMuId = -1, int fInnTrk=0, int fTrig = 0, int state=1
  ) {

  gStyle->SetOptStat(0);
  int kTrigSel = 13;

  float muPtCut = 3.5;
  float muEtaCut = 2.4;

  float massLow = 8.0;
  float massHigh = 10.0;

  if(state==2){massLow = 8.5; massHigh = 11;}

  double min = 0;
  double max = ptHigh;
  double binwidth = 1;
  const int numBins = 31;//(max-min)/binwidth;

  //input files
  TString inputMC = "/eos/cms/store/group/phys_heavyions/dileptons/MC2018/PbPb502TeV/TTrees/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC*20190801.root";
  if(state==2) inputMC = "/eos/cms/store/group/phys_heavyions/dileptons/MC2018/PbPb502TeV/TTrees/Upsi2S_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8_v1.root";
  TChain* mytree = new TChain("myTree"); 
  mytree->Add(inputMC.Data());

  //pT reweighting function
  TFile *fPtW = new TFile(Form("../Reweight/WeightedFunc/Func_dNdpT_%dS.root",state),"read");
  TF1* f1 = (TF1*) fPtW->Get("fitRatio");

  
  double ptBin[numBins+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,27,30,34,38,42,46,50};

  TH1D* hpt_reco = new TH1D("hpt_reco","hpt_reco",numBins,ptBin);
  TH1D* hptLowC_reco = new TH1D("hptLowC_reco","hptLowC_reco",numBins,ptBin);
  TH1D* hptMidC_reco = new TH1D("hptMidC_reco","hptMidC_reco",numBins,ptBin);
  TH1D* hptHighC_reco = new TH1D("hptHighC_reco","hptHighC_reco",numBins,ptBin);

  TH1D* hpt_gen = new TH1D("hpt_gen","hpt_gen",numBins,ptBin);
  TH1D* hptLowC_gen = new TH1D("hptLowC_gen","hptLowC_gen",numBins,ptBin);
  TH1D* hptMidC_gen = new TH1D("hptMidC_gen","hptMidC_gen",numBins,ptBin);
  TH1D* hptHighC_gen = new TH1D("hptHighC_gen","hptHighC_gen",numBins,ptBin);

  TH1D* hpt_reco_NoTrig = new TH1D("hpt_reco_NoTrig","hpt_reco_NoTrig",numBins,ptBin);
  TH1D* hptLowC_reco_NoTrig = new TH1D("hptLowC_reco_NoTrig","hptLowC_reco_NoTrig",numBins,ptBin);
  TH1D* hptMidC_reco_NoTrig = new TH1D("hptMidC_reco_NoTrig","hptMidC_reco_NoTrig",numBins,ptBin);
  TH1D* hptHighC_reco_NoTrig = new TH1D("hptHighC_reco_NoTrig","hptHighC_reco_NoTrig",numBins,ptBin);

  hpt_reco->Sumw2();
  hptLowC_reco->Sumw2();
  hptMidC_reco->Sumw2();
  hptHighC_reco->Sumw2();

  hpt_gen->Sumw2();
  hptLowC_gen->Sumw2();
  hptMidC_gen->Sumw2();
  hptHighC_gen->Sumw2();

  hpt_reco_NoTrig->Sumw2();
  hptLowC_reco_NoTrig->Sumw2();
  hptMidC_reco_NoTrig->Sumw2();
  hptHighC_reco_NoTrig->Sumw2();

  hpt_reco->SetTitle("Reco: Full Centrality");
  hptLowC_reco->SetTitle("Reco: Centrality 0-10%");
  hptMidC_reco->SetTitle("Reco: Centrality 10-50%");
  hptHighC_reco->SetTitle("Reco: Centrality 50-90%");

  hpt_reco->GetXaxis()->SetTitle("pt");
  hptLowC_reco->GetXaxis()->SetTitle("pt");
  hptMidC_reco->GetXaxis()->SetTitle("pt");
  hptHighC_reco->GetXaxis()->SetTitle("pt");

  hpt_gen->SetTitle("Gen: Full Centrality");
  hptLowC_gen->SetTitle("Gen: Centrality 0-10%");
  hptMidC_gen->SetTitle("Gen: Centrality 10-50%");
  hptHighC_gen->SetTitle("Gen: Centrality 50-90%");

  hpt_gen->GetXaxis()->SetTitle("pt");
  hptLowC_gen->GetXaxis()->SetTitle("pt");
  hptMidC_gen->GetXaxis()->SetTitle("pt");
  hptHighC_gen->GetXaxis()->SetTitle("pt");

  const int maxBranchSize = 1000;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Gen_QQ_size;
  Int_t           Gen_mu_size;
  TClonesArray    *Gen_QQ_4mom;
  TClonesArray    *Gen_mu_4mom;
  ULong64_t       Gen_QQ_trig[maxBranchSize];   //[Gen_QQ_size]
  Float_t         Gen_QQ_VtxProb[maxBranchSize];   //[Gen_QQ_size]
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Gen_QQ_size;   //!
  TBranch        *b_Gen_mu_size;   //!
  TBranch        *b_Gen_QQ_4mom;   //!
  TBranch        *b_Gen_mu_4mom;   //!
  TBranch        *b_Gen_QQ_trig;   //!
  TBranch        *b_Gen_QQ_VtxProb;   //!

  Gen_QQ_4mom = 0; Gen_mu_4mom = 0;
  mytree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
  mytree->SetBranchAddress("Gen_mu_size", &Gen_mu_size, &b_Gen_mu_size);
  mytree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
  mytree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  //  muon id 
  Int_t           Gen_QQ_mupl_idx[maxBranchSize];
  Int_t           Gen_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Gen_QQ_mupl_idx;
  TBranch        *b_Gen_QQ_mumi_idx;
  mytree->SetBranchAddress("Gen_QQ_mupl_idx",Gen_QQ_mupl_idx,&b_Gen_QQ_mupl_idx);
  mytree->SetBranchAddress("Gen_QQ_mumi_idx",Gen_QQ_mumi_idx,&b_Gen_QQ_mumi_idx);

  Int_t           Gen_mu_charge[maxBranchSize];
  TBranch        *b_Gen_mu_charge;   //!
  mytree->SetBranchAddress("Gen_mu_charge", Gen_mu_charge, &b_Gen_mu_charge);


  Int_t           Reco_QQ_size;
  Int_t           Reco_mu_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  ULong64_t       Reco_mu_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_mu_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Reco_QQ_4mom = 0; Reco_mu_4mom = 0;
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_mu_trig", Reco_mu_trig, &b_Reco_mu_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id 
  Int_t           Reco_QQ_mupl_idx[maxBranchSize];
  Int_t           Reco_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Reco_QQ_mupl_idx;
  TBranch        *b_Reco_QQ_mumi_idx;
  mytree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);

  Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dxy;   //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dz;   //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Int_t Reco_mu_whichGen[maxBranchSize];
  TBranch *b_Reco_mu_whichGen;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  mytree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
  mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);

  
  TLorentzVector* JP_Gen= new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;

  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  double weight = 1;
  double tnp_weight = 1;
  double tnp_trig_weight_mupl = -1;
  double tnp_trig_weight_mumi = -1;
  double pt_weight = 1;

  int kL2filter = 38;
  int kL3filter = 39;

  int count =0;
  int counttnp =0;
  const int nevt = mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    if(Centrality > cHigh || Centrality < cLow) continue;
    weight = findNcoll(Centrality) * Gen_weight;
    
    for(int igen = 0; igen<Gen_QQ_size; igen++){
      JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(igen);
      mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
      mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

      if(!( fabs(JP_Gen->Rapidity())<2.4 && (mupl_Gen->Pt()>muPtCut && fabs(mupl_Gen->Eta())<2.4) && (mumi_Gen->Pt()>muPtCut && fabs(mumi_Gen->Eta())<2.4) )) continue;
      if(Gen_mu_charge[Gen_QQ_mupl_idx[igen]]*Gen_mu_charge[Gen_QQ_mumi_idx[igen]]>0) continue;

      pt_weight = 1;
      if(isPtWeight) pt_weight = f1->Eval(JP_Gen->Pt()); 

      hpt_gen->Fill(JP_Gen->Pt(),weight*pt_weight);
      if(Centrality < 20) hptLowC_gen -> Fill(JP_Gen->Pt(), weight*pt_weight);
      else if(Centrality > 20 && Centrality < 100) hptMidC_gen -> Fill(JP_Gen->Pt(), weight*pt_weight);
      else if(Centrality > 100 && Centrality < 180) hptHighC_gen -> Fill(JP_Gen->Pt(), weight*pt_weight);
    }

  
    if(!((HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))) ) continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      
      if( !((Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel))) ) continue;

      if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
      if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;

      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)));

      bool muplSoft = (  //(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) &&
          passMuonTypePl        //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      bool mumiSoft = ( //(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  && 
          passMuonTypeMi       //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      if ( !(muplSoft && mumiSoft) ) continue;   
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) continue;
      if(Reco_QQ_sign[irqq]!=0) continue;  
      
      if(!( fabs(JP_Reco->Rapidity())<2.4 && (mupl_Reco->Pt()>muPtCut && fabs(mupl_Reco->Eta())<2.4) && (mumi_Reco->Pt()>muPtCut && fabs(mumi_Reco->Eta())<2.4) && fabs(JP_Reco->Rapidity())<2.4 && JP_Reco->Pt()<50  && JP_Reco->M()>massLow && JP_Reco->M()<massHigh)) continue;
     
      count++;
      if(isTnP){
       tnp_weight = 1;
       tnp_trig_weight_mupl = -1;
       tnp_trig_weight_mumi = -1;
       tnp_weight = tnp_weight * tnp_weight_muid_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), fMuId) * tnp_weight_muid_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), fMuId); //mu id
       tnp_weight = tnp_weight * tnp_weight_trk_pbpb(mupl_Reco->Eta(), fInnTrk) * tnp_weight_trk_pbpb(mumi_Reco->Eta(), fInnTrk); //inner tracker

       //Trigger part
       if(!((Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) && (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ){
//         cout << "irqq : " << irqq << " - iev : " << iev << endl;
//         cout << "TnP ERROR !!!! ::: No matched L2 filter1 " << endl;
         continue;
       }
       bool mupl_L2Filter = ( (Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ? true : false ;
       bool mupl_L3Filter = ( (Reco_mu_trig[Reco_QQ_mupl_idx[irqq]]&((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)) ) ? true : false ;
       bool mumi_L2Filter = ( (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL2filter))) == ((ULong64_t)pow(2, kL2filter)) ) ? true : false ;
       bool mumi_L3Filter = ( (Reco_mu_trig[Reco_QQ_mumi_idx[irqq]]&((ULong64_t)pow(2, kL3filter))) == ((ULong64_t)pow(2, kL3filter)) ) ? true : false ;
       if(mupl_L2Filter == false || mumi_L2Filter == false){ cout << "TnP ERROR !!!! ::: No matched L2 filter2 " << endl; cout << endl;} 

       bool mupl_isL2 = (mupl_L2Filter && !mupl_L3Filter) ? true : false;
       bool mupl_isL3 = (mupl_L2Filter && mupl_L3Filter) ? true : false;
       bool mumi_isL2 = (mumi_L2Filter && !mumi_L3Filter) ? true : false;
       bool mumi_isL3 = (mumi_L2Filter && mumi_L3Filter) ? true : false;
       bool SelDone = false;

       if( mupl_isL2 && mumi_isL3){
         tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, fTrig);
         tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, fTrig);
         SelDone = true;
       }
       else if( mupl_isL3 && mumi_isL2){
         tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, fTrig);
         tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, fTrig);
         SelDone = true;
       }
       else if( mupl_isL3 && mumi_isL3){
         int t[2] = {-1,1}; // mupl, mumi
         int l = rand() % (2); 
         //pick up what will be L2
         if(t[l]==-1){
           tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 2, fTrig);
           tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 3, fTrig);
         }
         else if(t[l]==1){
           tnp_trig_weight_mupl = tnp_weight_trg_pbpb(mupl_Reco->Pt(), mupl_Reco->Eta(), 3, fTrig);
           tnp_trig_weight_mumi = tnp_weight_trg_pbpb(mumi_Reco->Pt(), mumi_Reco->Eta(), 2, fTrig);
         }
         else {cout << "ERROR :: No random selection done !!!!" << endl; continue;}
         SelDone = true;
       }    

       if(SelDone == false || (tnp_trig_weight_mupl == -1 || tnp_trig_weight_mumi == -1)){continue;}
       tnp_weight = tnp_weight * tnp_trig_weight_mupl * tnp_trig_weight_mumi;
       counttnp++;
      }

      pt_weight = 1;
      if(isPtWeight) pt_weight = f1->Eval(JP_Reco->Pt());

      hpt_reco->Fill(JP_Reco->Pt(),weight * tnp_weight * pt_weight);
      if(Centrality < 20) hptLowC_reco -> Fill(JP_Reco->Pt(), weight * tnp_weight * pt_weight);
      else if(Centrality > 20 && Centrality < 100) hptMidC_reco -> Fill(JP_Reco->Pt(), weight * tnp_weight * pt_weight);
      else if(Centrality > 100 && Centrality < 180) hptHighC_reco -> Fill(JP_Reco->Pt(), weight * tnp_weight * pt_weight);
    }
  }

  cout << "count " << count << endl;
  cout << "counttnp " << counttnp << endl;
  


  //Draw
  TCanvas * cpt_reco = new TCanvas("cpt_reco","cpt_reco",0,0,400,400);
  cpt_reco->cd();
  hpt_reco->Draw();

  TCanvas * cptLowC_reco = new TCanvas("cptLowC_reco","cptLowC_reco",0,0,400,400);
  cptLowC_reco->cd();
  hptLowC_reco->Draw();

  TCanvas * cptMidC_reco = new TCanvas("cptMidC_reco","cptMidC_reco",400,0,400,400);
  cptMidC_reco->cd(); 
  hptMidC_reco->Draw();

  TCanvas * cptHighC_reco = new TCanvas("cptHighC_reco","cptHighC_reco",800,0,400,400);
  cptHighC_reco->cd();
  hptHighC_reco->Draw();


  //Gen
  TCanvas * cpt_gen = new TCanvas("cpt_gen","cpt_gen",0,400,400,400);
  cpt_gen->cd();
  hpt_gen->Draw();

  TCanvas * cptLowC_gen = new TCanvas("cptLowC_gen","cptLowC_gen",0,400,400,400);
  cptLowC_gen->cd();
  hptLowC_gen->Draw();

  TCanvas * cptMidC_gen = new TCanvas("cptMidC_gen","cptMidC_gen",0,400,400,400);
  cptMidC_gen->cd();
  hptMidC_gen->Draw();

  TCanvas * cptHighC_gen = new TCanvas("cptHighC_gen","cptHighC_gen",0,400,400,400);
  cptHighC_gen->cd();
  hptHighC_gen->Draw();

  cpt_reco->Update();
  cptLowC_reco->Update();
  cptMidC_reco->Update();
  cptHighC_reco->Update();

  //Divide
  TH1D* hpt_eff;
  TH1D* hptLowC_eff;
  TH1D* hptMidC_eff;
  TH1D* hptHighC_eff;

  TH1D* hpt_eff_NoTrig;
  TH1D* hptLowC_eff_NoTrig;
  TH1D* hptMidC_eff_NoTrig;
  TH1D* hptHighC_eff_NoTrig;

  TH1D* hpt_eff_Trig;
  TH1D* hptLowC_eff_Trig;
  TH1D* hptMidC_eff_Trig;
  TH1D* hptHighC_eff_Trig;

  hpt_eff = (TH1D*)hpt_reco->Clone("hpt_eff");
  hptLowC_eff = (TH1D*)hptLowC_reco->Clone("hptLowC_eff");
  hptMidC_eff = (TH1D*)hptMidC_reco->Clone("hptMidC_eff");
  hptHighC_eff = (TH1D*)hptHighC_reco->Clone("hptHighC_eff");

  hpt_eff->Divide(hpt_eff, hpt_gen, 1, 1, "B");
  hptLowC_eff->Divide(hptLowC_eff, hptLowC_gen, 1, 1, "B");
  hptMidC_eff->Divide(hptMidC_eff, hptMidC_gen, 1, 1, "B");
  hptHighC_eff->Divide(hptHighC_eff, hptHighC_gen, 1, 1, "B");

  hpt_eff->SetTitle("Eff: Full Centrality");
  hptLowC_eff->SetTitle("Eff: Centrality 0-10%");
  hptMidC_eff->SetTitle("Eff: Centrality 10-50%");
  hptHighC_eff->SetTitle("Eff: Centrality 50-90%");

  TCanvas * cpt_eff = new TCanvas("cpt_eff","cpt_eff",0,400,400,400);
  cpt_eff->cd();
  hpt_eff->Draw();

  TCanvas * cptLowC_eff = new TCanvas("cptLowC_eff","cptLowC_eff",0,400,400,400);
  cptLowC_eff->cd();
  hptLowC_eff->Draw();

  TCanvas * cptMidC_eff = new TCanvas("cptMidC_eff","cptMidC_eff",400,400,400,400);
  cptMidC_eff->cd();
  hptMidC_eff->Draw();

  TCanvas * cptHighC_eff = new TCanvas("cptHighC_eff","cptHighC_eff",800,400,400,400);
  cptHighC_eff->cd();
  hptHighC_eff->Draw();

  //Save efficiency files for later use.
  hpt_eff->SetName(Form("mc_eff_vs_pt_TnP%d_Cent090",isTnP));
  hptLowC_eff->SetName(Form("mc_eff_vs_pt_TnP%d_Cent010",isTnP));
  hptMidC_eff->SetName(Form("mc_eff_vs_pt_TnP%d_Cent1050",isTnP));
  hptHighC_eff->SetName(Form("mc_eff_vs_pt_TnP%d_Cent5090",isTnP));


  TString fMuIdSysString = "SysNom";
  TString fMuIdStatString = "StatNom";
  if(fMuId==-1) fMuIdSysString = "SysUp";
  else if(fMuId==-2) fMuIdSysString = "SysDo";
  else if(fMuId==1) fMuIdStatString = "StatUp";
  else if(fMuId==2) fMuIdStatString = "StatDo";

  TString fInnTrkSysString = "SysNom";
  TString fInnTrkStatString = "StatNom";
  if(fInnTrk==-1) fInnTrkSysString = "SysUp";
  else if(fInnTrk==-2) fInnTrkSysString = "SysDo";
  else if(fInnTrk==1) fInnTrkStatString = "StatUp";
  else if(fInnTrk==2) fInnTrkStatString = "StatDo";
  
  TString fTrigSysString = "SysNom";
  TString fTrigStatString = "StatNom";
  if(fTrig==-1) fTrigSysString = "SysUp";
  else if(fTrig==-2) fTrigSysString = "SysDo";
  else if(fTrig==1) fTrigStatString = "StatUp";
  else if(fTrig==2) fTrigStatString = "StatDo";


  TString outFileName = Form("mc_eff_Sys_Y%dS_muPtCut%.1f_MuId%s%s_InnTrk%s%s_Trig%s%s.root",state,muPtCut,fMuIdSysString.Data(),fMuIdStatString.Data(),fInnTrkSysString.Data(),fInnTrkStatString.Data(),fTrigSysString.Data(),fTrigStatString.Data());
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hpt_eff->Write();
  hptLowC_eff->Write();
  hptMidC_eff->Write();
  hptHighC_eff->Write();
  outFile->Close();

}
