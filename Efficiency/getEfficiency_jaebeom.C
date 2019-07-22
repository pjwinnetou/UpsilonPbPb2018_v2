#include <iostream>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"
#include "Style_jaebeom.h"

using namespace std;

void getEfficiency_jaebeom(
  float ptLow = 0.0, float ptHigh = 50.0,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 0, int cHigh = 180
  ) {

  gStyle->SetOptStat(0);
  int kTrigSel = 13;

  float muPtCut = 3.5;
  float muEtaCut = 2.4;

  float massLow = 8.0;
  float massHigh = 10.0;

  double min = 0;
  double max = ptHigh;
  double binwidth = 1;
  const int numBins = (max-min)/binwidth;

  //input files
  TString inputMC1 = "/pnfs/knu.ac.kr/data/cms/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC-v1/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_v1.root";
  TString inputMC2 = "/pnfs/knu.ac.kr/data/cms/store/group/phys_heavyions/hidilepton/2018PbPbMC/Ups1S_TuneCP5_HydjetDrumMB_5p02TeV_officialPythia8MC-v1/Upsi1S_TuneCP5_HydjetDrumMB_officialPythia8MC_ext-v1.root";
  TChain* mytree = new TChain("myTree"); 
  mytree->Add(inputMC1.Data());
  mytree->Add(inputMC2.Data());

  TH1D* hpt_reco = new TH1D("hpt_reco","hpt_reco",numBins,min,max);
  TH1D* hptLowestC_reco = new TH1D("hptLowestC_reco","hptLowestC_reco",numBins,min,max);
  TH1D* hptLowC_reco = new TH1D("hptLowC_reco","hptLowC_reco",numBins,min,max);
  TH1D* hptMidC_reco = new TH1D("hptMidC_reco","hptMidC_reco",numBins,min,max);
  TH1D* hptHighC_reco = new TH1D("hptHighC_reco","hptHighC_reco",numBins,min,max);

  TH1D* hpt_gen = new TH1D("hpt_gen","hpt_gen",numBins,min,max);
  TH1D* hptLowestC_gen = new TH1D("hptLowestC_gen","hptLowestC_gen",numBins,min,max);
  TH1D* hptLowC_gen = new TH1D("hptLowC_gen","hptLowC_gen",numBins,min,max);
  TH1D* hptMidC_gen = new TH1D("hptMidC_gen","hptMidC_gen",numBins,min,max);
  TH1D* hptHighC_gen = new TH1D("hptHighC_gen","hptHighC_gen",numBins,min,max);

  hpt_reco->Sumw2();
  hptLowestC_reco->Sumw2();
  hptLowC_reco->Sumw2();
  hptMidC_reco->Sumw2();
  hptHighC_reco->Sumw2();

  hpt_gen->Sumw2();
  hptLowestC_gen->Sumw2();
  hptLowC_gen->Sumw2();
  hptMidC_gen->Sumw2();
  hptHighC_gen->Sumw2();

  hpt_reco->SetTitle("Reco: Full Centrality");
  hptLowestC_reco->SetTitle("Reco: Centrality 0-10%");
  hptLowC_reco->SetTitle("Reco: Centrality 10-30%");
  hptMidC_reco->SetTitle("Reco: Centrality 30-50%");
  hptHighC_reco->SetTitle("Reco: Centrality 50-90%");

  hpt_reco->GetXaxis()->SetTitle("pt");
  hptLowestC_reco->GetXaxis()->SetTitle("pt");
  hptLowC_reco->GetXaxis()->SetTitle("pt");
  hptMidC_reco->GetXaxis()->SetTitle("pt");
  hptHighC_reco->GetXaxis()->SetTitle("pt");

  hpt_gen->SetTitle("Gen: Full Centrality");
  hptLowestC_gen->SetTitle("Gen: Centrality 0-10%");
  hptLowC_gen->SetTitle("Gen: Centrality 10-30%");
  hptMidC_gen->SetTitle("Gen: Centrality 30-50%");
  hptHighC_gen->SetTitle("Gen: Centrality 50-90%");

  hpt_gen->GetXaxis()->SetTitle("pt");
  hptLowestC_gen->GetXaxis()->SetTitle("pt");
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
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Reco_QQ_4mom = 0; Reco_mu_4mom = 0;
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
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
//  mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);

  
  TLorentzVector* JP_Gen= new TLorentzVector;
  TLorentzVector* mupl_Gen = new TLorentzVector;
  TLorentzVector* mumi_Gen = new TLorentzVector;

  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  double weight = 1;
  const int nevt = mytree->GetEntries();
  cout << "Total Events : " << nevt << endl;
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    weight = findNcoll(Centrality) * Gen_weight;
    
    for(int igen = 0; igen<Gen_QQ_size; igen++){
      JP_Gen = (TLorentzVector*) Gen_QQ_4mom->At(igen);
      mupl_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[igen]);
      mumi_Gen = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[igen]);

      if(!( (mupl_Gen->Pt()>3.5 && fabs(mupl_Gen->Eta())<2.4) && (mumi_Gen->Pt()>3.5 && fabs(mumi_Gen->Eta())<2.4) )) continue;
      if(Gen_mu_charge[Gen_QQ_mupl_idx[igen]]*Gen_mu_charge[Gen_QQ_mumi_idx[igen]]>0) continue;

      hpt_gen->Fill(JP_Gen->Pt(),weight);
      if(Centrality < 20) hptLowestC_gen -> Fill(JP_Gen->Pt(), weight);
      else if(Centrality > 20 && Centrality < 60) hptLowC_gen -> Fill(JP_Gen->Pt(), weight);
      else if(Centrality > 60 && Centrality < 100) hptMidC_gen -> Fill(JP_Gen->Pt(), weight);
      else if(Centrality > 100 && Centrality < 180) hptHighC_gen -> Fill(JP_Gen->Pt(), weight);
    }

  
    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      
      if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

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
      
      if(!( (mupl_Reco->Pt()>3.5 && fabs(mupl_Reco->Eta())<2.4) && (mumi_Reco->Pt()>3.5 && fabs(mumi_Reco->Eta())<2.4) && fabs(JP_Reco->Rapidity())<2.4 && JP_Reco->Pt()<50)) continue;
      
      hpt_reco->Fill(JP_Reco->Pt(),weight);
      if(Centrality < 20) hptLowestC_reco -> Fill(JP_Reco->Pt(), weight);
      else if(Centrality > 20 && Centrality < 60) hptLowC_reco -> Fill(JP_Reco->Pt(), weight);
      else if(Centrality > 60 && Centrality < 100) hptMidC_reco -> Fill(JP_Reco->Pt(), weight);
      else if(Centrality > 100 && Centrality < 180) hptHighC_reco -> Fill(JP_Reco->Pt(), weight);
    }
  }

  


  //Draw
  TCanvas * cpt_reco = new TCanvas("cpt_reco","cpt_reco",0,0,400,400);
  cpt_reco->cd();
  hpt_reco->Draw();

  TCanvas * cptLowestC_reco = new TCanvas("cptLowestC_reco","cptLowestC_reco",0,0,400,400);
  cptLowestC_reco->cd();
  hptLowestC_reco->Draw();

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

  TCanvas * cptLowestC_gen = new TCanvas("cptLowestC_gen","cptLowestC_gen",0,400,400,400);
  cptLowestC_gen->cd();
  hptLowestC_gen->Draw();

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
  cptLowestC_reco->Update();
  cptLowC_reco->Update();
  cptMidC_reco->Update();
  cptHighC_reco->Update();

  //Divide
  TH1D* hpt_eff;
  TH1D* hptLowestC_eff;
  TH1D* hptLowC_eff;
  TH1D* hptMidC_eff;
  TH1D* hptHighC_eff;

  hpt_eff = (TH1D*)hpt_reco->Clone("hpt_eff");
  hptLowestC_eff = (TH1D*)hptLowestC_reco->Clone("hptLowestC_eff");
  hptLowC_eff = (TH1D*)hptLowC_reco->Clone("hptLowC_eff");
  hptMidC_eff = (TH1D*)hptMidC_reco->Clone("hptMidC_eff");
  hptHighC_eff = (TH1D*)hptHighC_reco->Clone("hptHighC_eff");

  hpt_eff->Divide(hpt_eff, hpt_gen, 1, 1, "B");
  hptLowestC_eff->Divide(hptLowestC_eff, hptLowestC_gen, 1, 1, "B");
  hptLowC_eff->Divide(hptLowC_eff, hptLowC_gen, 1, 1, "B");
  hptMidC_eff->Divide(hptMidC_eff, hptMidC_gen, 1, 1, "B");
  hptHighC_eff->Divide(hptHighC_eff, hptHighC_gen, 1, 1, "B");

  hpt_eff->SetTitle("Eff: Full Centrality");
  hptLowestC_eff->SetTitle("Eff: Centrality 0-10%");
  hptLowC_eff->SetTitle("Eff: Centrality 10-30%");
  hptMidC_eff->SetTitle("Eff: Centrality 30-50%");
  hptHighC_eff->SetTitle("Eff: Centrality 50-90%");

  TCanvas * cpt_eff = new TCanvas("cpt_eff","cpt_eff",0,400,400,400);
  cpt_eff->cd();
  hpt_eff->Draw();

  TCanvas * cptLowestC_eff = new TCanvas("cptLowestC_eff","cptLowestC_eff",0,400,400,400);
  cptLowestC_eff->cd();
  hptLowestC_eff->Draw();

  TCanvas * cptLowC_eff = new TCanvas("cptLowC_eff","cptLowC_eff",0,400,400,400);
  cptLowC_eff->cd();
  hptLowC_eff->Draw();

  TCanvas * cptMidC_eff = new TCanvas("cptMidC_eff","cptMidC_eff",400,400,400,400);
  cptMidC_eff->cd();
  hptMidC_eff->Draw();

  TCanvas * cptHighC_eff = new TCanvas("cptHighC_eff","cptHighC_eff",800,400,400,400);
  cptHighC_eff->cd();
  hptHighC_eff->Draw();

  //Save canvases
  cpt_reco->SaveAs("Plots/cpt_reco.png");
  cptLowestC_reco->SaveAs("Plots/cptLowestC_reco.png");
  cptLowC_reco->SaveAs("Plots/cptLowC_reco.png");
  cptMidC_reco->SaveAs("Plots/cptMidC_reco.png");
  cptHighC_reco->SaveAs("Plots/cptHighC_reco.png");
  cpt_reco->SaveAs("Plots/cpt_reco.pdf");
  cptLowestC_reco->SaveAs("Plots/cptLowestC_reco.pdf");
  cptLowC_reco->SaveAs("Plots/cptLowC_reco.pdf");
  cptMidC_reco->SaveAs("Plots/cptMidC_reco.pdf");
  cptHighC_reco->SaveAs("Plots/cptHighC_reco.pdf");

  cpt_gen->SaveAs("Plots/cpt_gen.png");
  cptLowestC_gen->SaveAs("Plots/cptLowestC_gen.png");
  cptLowC_gen->SaveAs("Plots/cptLowC_gen.png");
  cptMidC_gen->SaveAs("Plots/cptMidC_gen.png");
  cptHighC_gen->SaveAs("Plots/cptHighC_gen.png");
  cpt_gen->SaveAs("Plots/cpt_gen.pdf");
  cptLowestC_gen->SaveAs("Plots/cptLowestC_gen.pdf");
  cptLowC_gen->SaveAs("Plots/cptLowC_gen.pdf");
  cptMidC_gen->SaveAs("Plots/cptMidC_gen.pdf");
  cptHighC_gen->SaveAs("Plots/cptHighC_gen.pdf");

  cpt_eff->SaveAs("Plots/cpt_eff.png");
  cptLowestC_eff->SaveAs("Plots/cptLowestC_eff.png");
  cptLowC_eff->SaveAs("Plots/cptLowC_eff.png");
  cptMidC_eff->SaveAs("Plots/cptMidC_eff.png");
  cptHighC_eff->SaveAs("Plots/cptHighC_eff.png");
  cpt_eff->SaveAs("Plots/cpt_eff.pdf");
  cptLowestC_eff->SaveAs("Plots/cptLowestC_eff.pdf");
  cptLowC_eff->SaveAs("Plots/cptLowC_eff.pdf");
  cptMidC_eff->SaveAs("Plots/cptMidC_eff.pdf");
  cptHighC_eff->SaveAs("Plots/cptHighC_eff.pdf");

  //Save efficiency files for later use.
  hpt_eff->SetName("mc_eff_vs_pt_noTnP_Cent0100");
  hptLowestC_eff->SetName("mc_eff_vs_pt_noTnP_Cent010");
  hptLowC_eff->SetName("mc_eff_vs_pt_noTnP_Cent1030");
  hptMidC_eff->SetName("mc_eff_vs_pt_noTnP_Cent3050");
  hptHighC_eff->SetName("mc_eff_vs_pt_noTnP_Cent5090");
  TString outFileName = "mc_eff_vs_pt_noTnP_OfficialMC.root";
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hpt_eff->Write();
  hptLowestC_eff->Write();
  hptLowC_eff->Write();
  hptMidC_eff->Write();
  hptHighC_eff->Write();
  outFile->Close();

}
