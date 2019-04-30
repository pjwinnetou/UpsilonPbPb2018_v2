#include <ctime>
#include <iostream>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"
#include "SkimTree_MuonID.h"

static const long MAXTREESIZE = 1000000000000;

void SkimTree_MuonID(int nevt=-1, bool isMC = false, bool isJPsiTrig = true) 
{
  
  TString fMCstr = (isMC) ? "MC" : "DATA" ;

  using namespace std;
  using namespace hi;


  TString fnameData1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part*.root";
  TString fnameData2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
  TString fnameMC = "/eos/cms/store/group/phys_heavyions/dileptons/MC2018/PbPb502TeV/TTrees/Oniatree_MC_PromptJpsi_5TeV_PbPbEmbd_20190326.root";
  
  TChain *mytree;
  if(!isMC){
    mytree = new TChain("myTree");
    mytree->Add(fnameData1.Data());
    mytree->Add(fnameData2.Data());
  }
  else if(isMC){
    mytree = new TChain("hionia/myTree");
    mytree->Add(fnameMC.Data());
  }
  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
  
  mytree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  mytree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  mytree->SetBranchAddress("Reco_mu_normChi2_inner", Reco_mu_normChi2_inner, &b_Reco_mu_normChi2_inner);
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  mytree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  mytree->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);
  mytree->SetBranchAddress("Reco_mu_ptErr_inner", Reco_mu_ptErr_inner, &b_Reco_mu_ptErr_inner);
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);

  Int_t Reco_mu_whichGen[maxBranchSize];
  TBranch *b_Reco_mu_whichGen;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  if(isMC){
    mytree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
    mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);
  }

  float px_, py_, pz_, x_, y_, z_, t_, e_;

  TFile* newfile;
  newfile = new TFile(Form("Onia_MuId_SkimHist_%s_isJPsiTrig%d.root",fMCstr.Data(),isJPsiTrig),"recreate");

  Int_t nMuCand;

  TLorentzVector *MuMu_Reco = new TLorentzVector;
  TLorentzVector *MuPl_Reco = new TLorentzVector;
  TLorentzVector *MuMi_Reco = new TLorentzVector;

  map<TString, TH1D*> hNormChi2Global, hmudxy, hmudxyErr, hmudz, hmudzErr, hnTrkHits, hnMuValHits, hnTrkWMea, hmuTMOneStaTight, hmunPixWMea, hmuStationsMatched, hmunPixValHits, hmuptErrglobal, hmupt, hmumupt, hmuphi, hmueta, hmumuy, hVtxProb, hmudxy_den, hmudz_den, hmunPixWMea_den, hnTrkWMea_den;

  const int nHistType  = 12;
  const char* histType[nHistType] = {"All", "Sig", "Bkg", "Glb", "GlbTrk", "GlbNTrk", "GlbSig", "GlbBkg", "GlbTrkSig", "GlbTrkBkg", "GlbNTrkSig", "GlbNTrkBkg"};
  TObjArray* fOBj[nHistType];

  for(int ihist = 0; ihist<nHistType; ihist++){
    hNormChi2Global[Form("%s",histType[ihist])]   = new TH1D(Form("hNormChi2Global_%s",histType[ihist]),";Reco_mu_normChi2_global;",1000,-10,100);
    hnTrkHits[Form("%s",histType[ihist])]         = new TH1D(Form("hnTrkHits_%s",histType[ihist]),";Reco_mu_nTrkHits;",40,0,40);
    hnMuValHits[Form("%s",histType[ihist])]       = new TH1D(Form("hnMuValHits_%s",histType[ihist]),";Reco_mu_nMuValHits;",55,0,55);
    hmudxy[Form("%s",histType[ihist])]            = new TH1D(Form("hmudxy_%s",histType[ihist]),";Reco_mu_dxy;",50,0,0.3);
    hmudxyErr[Form("%s",histType[ihist])]         = new TH1D(Form("hmudxyErr_%s",histType[ihist]),";Reco_mu_dxyErr;",500,0,0.1);
    hmudz[Form("%s",histType[ihist])]             = new TH1D(Form("hmudz_%s",histType[ihist]),";Reco_mu_dz;",200,0,5);
    hmudzErr[Form("%s",histType[ihist])]          = new TH1D(Form("hmudzErr_%s",histType[ihist]),";Reco_mu_dzErr;",500,0,0.3);
    hnTrkWMea[Form("%s",histType[ihist])]         = new TH1D(Form("hnTrkWMea_%s",histType[ihist]),";Reco_mu_nTrkWMea;",20,0,20);
    hmuTMOneStaTight[Form("%s",histType[ihist])]  = new TH1D(Form("hmuTMOneStaTight_%s",histType[ihist]),";Reco_mu_TMOneStaTight;",2,0,2);
    hmunPixWMea[Form("%s",histType[ihist])]       = new TH1D(Form("hmunPixWMea_%s",histType[ihist]),";Reco_mu_nPixWMea;",7,0,7);
    hmunPixValHits[Form("%s",histType[ihist])]    = new TH1D(Form("hmunPixValHits_%s",histType[ihist]),";Reco_mu_nPixValHits;",12,0,12);
    hmuptErrglobal[Form("%s",histType[ihist])]    = new TH1D(Form("hmuptErrglobal_%s",histType[ihist]),";Reco_mu_ptErr_global;",1000,0,1);
    hmupt[Form("%s",histType[ihist])]             = new TH1D(Form("hmupt_%s",histType[ihist]),";p_{T}^{#mu};",250,0,50);
    hmumupt[Form("%s",histType[ihist])]           = new TH1D(Form("hmumupt_%s",histType[ihist]),";p_{T}^{#mu#mu};",250,0,50);
    hmuphi[Form("%s",histType[ihist])]            = new TH1D(Form("hmuphi_%s",histType[ihist]),";#phi^{#mu};",300,-4,4);
    hmueta[Form("%s",histType[ihist])]            = new TH1D(Form("hmueta_%s",histType[ihist]),";#eta^{#mu};",300,-4,4);
    hmumuy[Form("%s",histType[ihist])]            = new TH1D(Form("hmumuy_%s",histType[ihist]),";y^{#mu#mu};",300,-4,4);
    hVtxProb[Form("%s",histType[ihist])]          = new TH1D(Form("hVtxProb_%s",histType[ihist]),";Reco_QQ_VtxProb;",200,0,1);
    hmudz_den[Form("%s",histType[ihist])]         = new TH1D(Form("hmudz_%s_den",histType[ihist]),";Reco_mu_dz;",200,0,5);
    hmudxy_den[Form("%s",histType[ihist])]        = new TH1D(Form("hmuxy_%s_den",histType[ihist]),";Reco_mu_dz;",200,0,5);
    hmunPixWMea_den[Form("%s",histType[ihist])]   = new TH1D(Form("hmunPixWMea_%s_den",histType[ihist]),";Reco_mu_nPixWMea;",7,0,7);
    hnTrkWMea_den[Form("%s",histType[ihist])]     = new TH1D(Form("hnTrkWMea_%s_den",histType[ihist]),";Reco_mu_nTrkWMea;",20,0,20);
    fOBj[ihist] = new TObjArray();
  }
  
  int kTrigSel = 12;
  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();
  cout << "Total events = " << mytree->GetEntries() << endl;
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);

    if(Reco_mu_size==0) continue;

    nMuCand = 0;
    if(isJPsiTrig && (!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) )) continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      MuMu_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      MuPl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      MuMi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      
      if(Reco_QQ_sign[irqq]!=0 || MuMu_Reco->M() < 2.6 || MuMu_Reco->M()>3.5) continue;
      
      if(isJPsiTrig){
        if(!IsAcceptanceQQ(MuPl_Reco->Pt(), MuPl_Reco->Eta()) || !IsAcceptanceQQ(MuMi_Reco->Pt(), MuMi_Reco->Eta())) continue;
        if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;
      }
      else if(!isJPsiTrig){
        if(!IsAcceptanceNoTrig(MuPl_Reco->Pt(), MuPl_Reco->Eta()) || !IsAcceptanceNoTrig(MuMi_Reco->Pt(), MuMi_Reco->Eta())) continue;
      }

      bool passGlbMuPl = ((Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)))>0) ? true : false;
      bool passTrkMuPl = ((Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)))>0) ? true : false;
      bool passGlbMuMi = ((Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)))>0) ? true : false;
      bool passTrkMuMi = ((Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)))>0) ? true : false;
      
      bool passMuonTypeGlbMuPl = true;
      bool passMuonTypeGlbMuMi = true;
      bool passMuonTypeGlbTrkMuPl = true;
      bool passMuonTypeGlbTrkMuMi = true;
      bool passMuonTypeGlbNTrkMuPl = true;
      bool passMuonTypeGlbNTrkMuMi = true;
      
      passMuonTypeGlbMuPl = passMuonTypeGlbMuPl && passGlbMuPl;
      passMuonTypeGlbMuMi = passMuonTypeGlbMuMi && passGlbMuMi;
      passMuonTypeGlbTrkMuPl = passMuonTypeGlbTrkMuPl && passGlbMuPl && passTrkMuPl;
      passMuonTypeGlbTrkMuMi = passMuonTypeGlbTrkMuMi && passGlbMuMi && passTrkMuMi;
      passMuonTypeGlbNTrkMuPl = passMuonTypeGlbNTrkMuPl && passGlbMuPl && (!passTrkMuPl);
      passMuonTypeGlbNTrkMuMi = passMuonTypeGlbNTrkMuMi && passGlbMuMi && (!passTrkMuMi);

      px_ = MuMu_Reco->Px();
      py_ = MuMu_Reco->Py();
      pz_ = MuMu_Reco->Pz();
      e_ = MuMu_Reco->E();
      x_ = MuMu_Reco->X();
      y_ = MuMu_Reco->Y();
      z_ = MuMu_Reco->Z();
      t_ = MuMu_Reco->T();

      new(Reco_QQ_4mom_Lorentz[irqq]) TLorentzVector(x_,y_,z_,t_);
      
      if(isMC){
        if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      }

      bool isSignal = false;
      bool isBkg = false;
      if( MuMu_Reco->M() > 3.0 && MuMu_Reco->M() < 3.2) isSignal = true;
      else if( (MuMu_Reco->M() > 2.6 && MuMu_Reco->M() < 2.8) || (MuMu_Reco->M() > 3.3 && MuMu_Reco->M() < 3.5) ) isBkg = true;

      double weight = 1;
      if(isMC) weight = findNcoll(Centrality) * Gen_weight;

      bool isHistPass[nHistType] = {true, (isSignal), (isBkg), (passMuonTypeGlbMuMi && passMuonTypeGlbMuPl), (passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl), (passMuonTypeGlbNTrkMuMi && passMuonTypeGlbNTrkMuPl), (isSignal && passMuonTypeGlbMuMi && passMuonTypeGlbMuPl), (isBkg && passMuonTypeGlbMuMi && passMuonTypeGlbMuPl), (isSignal && passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl), (isBkg && passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl), (isSignal && passMuonTypeGlbNTrkMuMi && passMuonTypeGlbNTrkMuPl), (isBkg && passMuonTypeGlbNTrkMuMi && passMuonTypeGlbNTrkMuPl)};

      for(int ihist =0; ihist<nHistType; ihist++){
        if(isHistPass[ihist]){
          hNormChi2Global[histType[ihist]]   -> Fill(Reco_mu_normChi2_global[Reco_QQ_mumi_idx[irqq]],weight);
          hNormChi2Global[histType[ihist]]   -> Fill(Reco_mu_normChi2_global[Reco_QQ_mupl_idx[irqq]],weight);
          hnTrkHits[histType[ihist]]         -> Fill(Reco_mu_nTrkHits[Reco_QQ_mumi_idx[irqq]],weight);
          hnTrkHits[histType[ihist]]         -> Fill(Reco_mu_nTrkHits[Reco_QQ_mupl_idx[irqq]],weight);
          hnMuValHits[histType[ihist]]       -> Fill(Reco_mu_nMuValHits[Reco_QQ_mumi_idx[irqq]],weight);
          hnMuValHits[histType[ihist]]       -> Fill(Reco_mu_nMuValHits[Reco_QQ_mupl_idx[irqq]],weight);
          hmudxyErr[histType[ihist]]         -> Fill(Reco_mu_dxyErr[Reco_QQ_mumi_idx[irqq]],weight);
          hmudxyErr[histType[ihist]]         -> Fill(Reco_mu_dxyErr[Reco_QQ_mupl_idx[irqq]],weight);
          hmudzErr[histType[ihist]]          -> Fill(Reco_mu_dzErr[Reco_QQ_mumi_idx[irqq]],weight);
          hmudzErr[histType[ihist]]          -> Fill(Reco_mu_dzErr[Reco_QQ_mupl_idx[irqq]],weight);
          hmunPixValHits[histType[ihist]]    -> Fill(Reco_mu_nPixValHits[Reco_QQ_mumi_idx[irqq]],weight);
          hmunPixValHits[histType[ihist]]    -> Fill(Reco_mu_nPixValHits[Reco_QQ_mupl_idx[irqq]],weight);
          hmuptErrglobal[histType[ihist]]    -> Fill(Reco_mu_ptErr_global[Reco_QQ_mumi_idx[irqq]],weight);
          hmuptErrglobal[histType[ihist]]    -> Fill(Reco_mu_ptErr_global[Reco_QQ_mupl_idx[irqq]],weight);
          hmupt[histType[ihist]]             -> Fill(MuMi_Reco->Pt(),weight);
          hmupt[histType[ihist]]             -> Fill(MuPl_Reco->Pt(),weight);
          hmumupt[histType[ihist]]           -> Fill(MuMu_Reco->Pt(),weight);
          hmuphi[histType[ihist]]            -> Fill(MuMi_Reco->Phi(),weight);
          hmuphi[histType[ihist]]            -> Fill(MuPl_Reco->Phi(),weight);
          hmueta[histType[ihist]]            -> Fill(MuMi_Reco->Eta(),weight);
          hmueta[histType[ihist]]            -> Fill(MuPl_Reco->Eta(),weight);
          hmumuy[histType[ihist]]            -> Fill(MuMu_Reco->Rapidity(),weight);
          hVtxProb[histType[ihist]]          -> Fill(Reco_QQ_VtxProb[irqq],weight);
          hmuTMOneStaTight[histType[ihist]]  -> Fill(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]],weight);
          hmuTMOneStaTight[histType[ihist]]  -> Fill(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]],weight);
          hnTrkWMea[histType[ihist]]         -> Fill(Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]],weight);
          hnTrkWMea[histType[ihist]]         -> Fill(Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]],weight);
          hmunPixWMea[histType[ihist]]       -> Fill(Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]],weight);
          hmunPixWMea[histType[ihist]]       -> Fill(Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]],weight);
          hmudxy[histType[ihist]]            -> Fill(fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]), weight);
          hmudxy[histType[ihist]]            -> Fill(fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]), weight);
          hmudz[histType[ihist]]             -> Fill(fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]), weight);
          hmudz[histType[ihist]]             -> Fill(fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]), weight);
          
          if(Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]]>5 && Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]]>0 && fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3 && fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20 && Reco_QQ_VtxProb[irqq]>0.01 && MuPl_Reco->Pt()>3.5){
            hmudxy_den[histType[ihist]] -> Fill(fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]),weight);
            hmudz_den[histType[ihist]]  -> Fill(fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]),weight);
            hnTrkWMea_den[histType[ihist]] -> Fill(Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]],weight);
            hmunPixWMea_den[histType[ihist]] -> Fill(Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]],weight);
          }
          if(Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]]>5 && Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]]>0 && fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3 && fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20 && Reco_QQ_VtxProb[irqq]>0.01 && MuMi_Reco->Pt()>3.5){
            hmudxy_den[histType[ihist]] -> Fill(fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]),weight);
            hmudz_den[histType[ihist]]  -> Fill(fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]),weight);
            hnTrkWMea_den[histType[ihist]] -> Fill(Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]],weight);
            hmunPixWMea_den[histType[ihist]] -> Fill(Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]],weight);
          }
        }
      }
      nMuCand++; 

    } // end of dimuon loop

  } //end of event loop
 
  newfile->cd();

  for(int ihist = 0; ihist<nHistType; ihist++){
    fOBj[ihist]  ->  Add(hNormChi2Global[histType[ihist]]); 
    fOBj[ihist]  ->  Add(hnTrkHits[histType[ihist]]);
    fOBj[ihist]  ->  Add(hnMuValHits[histType[ihist]]);    
    fOBj[ihist]  ->  Add(hmudxyErr[histType[ihist]]);       
    fOBj[ihist]  ->  Add(hmudzErr[histType[ihist]]);        
    fOBj[ihist]  ->  Add(hmunPixValHits[histType[ihist]]);  
    fOBj[ihist]  ->  Add(hmuptErrglobal[histType[ihist]]);  
    fOBj[ihist]  ->  Add(hmupt[histType[ihist]]);           
    fOBj[ihist]  ->  Add(hmumupt[histType[ihist]]);         
    fOBj[ihist]  ->  Add(hmuphi[histType[ihist]]);          
    fOBj[ihist]  ->  Add(hmueta[histType[ihist]]);          
    fOBj[ihist]  ->  Add(hmumuy[histType[ihist]]);          
    fOBj[ihist]  ->  Add(hVtxProb[histType[ihist]]);        
    fOBj[ihist]  ->  Add(hmuTMOneStaTight[histType[ihist]]);
    fOBj[ihist]  ->  Add(hnTrkWMea[histType[ihist]]);       
    fOBj[ihist]  ->  Add(hmunPixWMea[histType[ihist]]);     
    fOBj[ihist]  ->  Add(hmudxy[histType[ihist]]);          
    fOBj[ihist]  ->  Add(hmudz[histType[ihist]]);

    fOBj[ihist]  ->  Add(hnTrkWMea_den[histType[ihist]]);       
    fOBj[ihist]  ->  Add(hmunPixWMea_den[histType[ihist]]);     
    fOBj[ihist]  ->  Add(hmudxy_den[histType[ihist]]);          
    fOBj[ihist]  ->  Add(hmudz_den[histType[ihist]]);

    fOBj[ihist] -> Write(Form("f%s",histType[ihist]),TObject::kOverwrite | TObject::kSingleKey);
  }    

  newfile->Close();
}

