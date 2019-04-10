#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"

static const long MAXTREESIZE = 1000000000000;

void SkimTree_MuonID(int nevt=-1, int kMuId = kMuGlb, bool isMC = false, bool isSignal = true) 
{
  TString muIdStr;
  if(kMuId == kMuGlb) muIdStr = "Glb";
  if(kMuId == kMuGlbTrk) muIdStr = "GlbTrk";
  
  TString fMCstr = (isMC) ? "MC" : "DATA" ;
  TString fSigstr = (isSignal) ? "Sig" : "Bkg" ;

  using namespace std;
  using namespace hi;


  TString fnameData1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part*.root";
  TString fnameData2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
  TString fnameMC = "/eos/cms/store/group/phys_heavyions/dileptons/MC2018/PbPb502TeV/TTrees/Oniatree_MC_PromptJpsi_5TeV_PbPbEmbd_20190326.root";
  
  TChain *mytree = new TChain("myTree");
  if(!isMC){
    mytree->Add(fnameData1.Data());
    mytree->Add(fnameData2.Data());
  }
  else if(isMC){mytree->Add(fnameMC.Data());}
  
  const int maxBranchSize = 1000;

  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Reco_QQ_size;
  Int_t           Reco_mu_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_mu_4mom;
  TClonesArray    *Reco_QQ_4mom_ = new TClonesArray("TLorentzVector");
  TClonesArray    &Reco_QQ_4mom_Lorentz = *Reco_QQ_4mom_;
  TClonesArray    *Reco_mu_4mom_ = new TClonesArray("TLorentzVector");
  TClonesArray    &Reco_mu_4mom_Lorentz = *Reco_mu_4mom_;

  ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_mu_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_mu_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Bool_t          Reco_mu_highPurity[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_highPurity;   //!
  mytree->SetBranchAddress("Reco_mu_highPurity", Reco_mu_highPurity, &b_Reco_mu_highPurity);

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
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

  //  muon id 
  Int_t           Reco_QQ_mupl_idx[maxBranchSize];
  Int_t           Reco_QQ_mumi_idx[maxBranchSize];
  TBranch        *b_Reco_QQ_mupl_idx;
  TBranch        *b_Reco_QQ_mumi_idx;
  mytree->SetBranchAddress("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx,&b_Reco_QQ_mupl_idx);
  mytree->SetBranchAddress("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx,&b_Reco_QQ_mumi_idx);

  Int_t           Reco_mu_nTrkHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_mu_normChi2_global[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_mu_nMuValHits[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_mu_StationsMatched[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_TMOneStaTight;   //!

  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
//  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_mu_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_mu_nPixValHits", Reco_mu_nPixValHits, &b_Reco_mu_nPixValHits);
  Float_t         Reco_mu_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
  TBranch        *b_Reco_mu_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_mu_ptErr_global", Reco_mu_ptErr_global, &b_Reco_mu_ptErr_global);

  Int_t           Reco_mu_SelectionType[maxBranchSize];
  TBranch        *b_Reco_mu_SelectionType;
  mytree->SetBranchAddress("Reco_mu_SelectionType", Reco_mu_SelectionType, &b_Reco_mu_SelectionType);

  Int_t Reco_mu_whichGen[maxBranchSize];
  if(isMC){
    TBranch *b_Reco_mu_whichGen;
    mytree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
  }


  TFile* newfile;
  newfile = new TFile(Form("Onia_MuId_SkimHist_%s_%s_%s.root",muIdStr.Data(),fMCstr.Data(),fSigstr.Data()),"recreate");

  Int_t nMuCand;
  TTree* newMytree = new TTree("newTree","skimmedTree");
  newMytree -> Branch("eventNb",&eventNb,"eventNb/i");
  newMytree -> Branch("zVtx",&zVtx,"zVtx/F");
  newMytree -> Branch("Centrality", &Centrality, "Centrality/I");
  newMytree -> Branch("HLTriggers", &HLTriggers, "HLTriggers/l");
  newMytree -> Branch("Reco_QQ_size", &Reco_QQ_size, "Reco_QQ_size/I");
  newMytree -> Branch("nMuCand",&nMuCand,"nMuCand/I");
  newMytree -> Branch("Reco_mu_size", &Reco_mu_size, "Reco_mu_size/I");
  newMytree -> Branch("Reco_QQ_4mom", "TClonesArray", &Reco_QQ_4mom_, 32000, 0);
//  newMytree -> Branch("Reco_mu_4mom", "TClonesArray", &Reco_mu_4mom_, 32000, 0);
  newMytree -> Branch("Reco_QQ_trig", Reco_QQ_trig, "Reco_QQ_trig[Reco_QQ_size]/l");
  newMytree -> Branch("Reco_QQ_VtxProb", Reco_QQ_VtxProb, "Reco_QQ_VtxProb[Reco_QQ_size]/F");
  newMytree -> Branch("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx, "Reco_QQ_mupl_idx[Reco_QQ_size]/I");
  newMytree -> Branch("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx, "Reco_QQ_mumi_idx[Reco_QQ_size]/I");
  newMytree -> Branch("Reco_mu_nTrkHits", Reco_mu_nTrkHits, "Reco_mu_nTrkHits[Reco_mu_size]/I");
  newMytree -> Branch("Reco_mu_normChi2_global", Reco_mu_normChi2_global, "Reco_mu_normChi2_global[Reco_mu_size]/F");
  newMytree -> Branch("Reco_mu_nMuValHits", Reco_mu_nMuValHits, "Reco_mu_nMuValHits[Reco_mu_size]/I");
  newMytree -> Branch("Reco_mu_StationsMatched", Reco_mu_StationsMatched, "Reco_mu_StationsMatched[Reco_mu_size]/I");
  newMytree -> Branch("Reco_mu_dxy", Reco_mu_dxy, "Reco_mu_dxy[Reco_mu_size]/F");
  newMytree -> Branch("Reco_mu_dxyErr", Reco_mu_dxyErr, "Reco_mu_dxyErr[Reco_mu_size]/F");
  newMytree -> Branch("Reco_mu_dz", Reco_mu_dz, "Reco_mu_dz[Reco_mu_size]/F");
  newMytree -> Branch("Reco_mu_dzErr", Reco_mu_dzErr, "Reco_mu_dzErr[Reco_mu_size]/F");
  newMytree -> Branch("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, "Reco_mu_nTrkWMea[Reco_mu_size]/I");
  newMytree -> Branch("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, "Reco_mu_TMOneStaTight[Reco_mu_size]/O");
  newMytree -> Branch("Reco_mu_nPixWMea", Reco_mu_nPixWMea, "Reco_mu_nPixWMea[Reco_mu_size]/I");
  newMytree -> Branch("Reco_QQ_sign", Reco_QQ_sign, "Reco_QQ_sign[Reco_QQ_size]/I");
  newMytree -> Branch("Reco_mu_nPixValHits", Reco_mu_nPixValHits, "Reco_mu_nPixValHits[Reco_mu_size]/I");
  newMytree -> Branch("Reco_mu_ptErr_global", Reco_mu_ptErr_global, "Reco_mu_ptErr_global[Reco_mu_size]/F");
  newMytree -> Branch("Reco_mu_SelectionType", Reco_mu_SelectionType, "Reco_mu_SelectionType[Reco_mu_size]/I");


  float px_, py_, pz_, e_, x_, y_, z_, m_, t_;

  TLorentzVector *MuMu_Reco = new TLorentzVector;
  TLorentzVector *MuPl_Reco = new TLorentzVector;
  TLorentzVector *MuMi_Reco = new TLorentzVector;

  TH1D* hCent = new TH1D("hCent",";Centrality;",200,0,200);
  TH1D* hNormChi2 = new TH1D("hNormChi2",";Reco_mu_normChi2_global;",1000,-10,100);
  TH1D* hnTrkHits = new TH1D("hnTrkHits",";Reco_mu_nTrkHits;",40,0,40);
  TH1D* hnMuValHits = new TH1D("hnMuValHits",";Reco_mu_nMuValHits;",55,0,55);
  TH1D* hmudxy = new TH1D("hmudxy",";Reco_mu_dxy;",100,-0.3,0.3);
  TH1D* hmudxyErr = new TH1D("hmudxyErr",";Reco_mu_dxyErr;",500,0,0.1);
  TH1D* hmudz = new TH1D("hmudz",";Reco_mu_dz;",1000,-10,10);
  TH1D* hmudzErr = new TH1D("hmudzErr",";Reco_mu_dzErr;",500,0,0.3);
  TH1D* hnTrkWMea = new TH1D("hnTrkWMea",";Reco_mu_nTrkWMea;",20,0,20);
  TH1D* hmuTMOneStaTight = new TH1D("hmuTMOneStaTight",";Reco_mu_TMOneStaTight;",2,0,2);
  TH1D* hmunPixWMea = new TH1D("hmunPixWMea",";Reco_mu_nPixWMea;",7,0,7);
  TH1D* hmunPixValHits = new TH1D("hmunPixValHits",";Reco_mu_nPixValHits;",12,0,12);
  TH1D* hmuptErrglobal = new TH1D("hmuptErrglobal",";Reco_mu_ptErr_global;",1000,0,1);
  TH1D* hmupt = new TH1D("hmupt",";p_{T}^{#mu};",250,0,50);
  TH1D* hmumupt = new TH1D("hmumupt",";p_{T}^{#mu#mu};",250,0,50);
  TH1D* hmuphi = new TH1D("hmuphi",";#phi^{#mu};",300,-4,4);
  TH1D* hmueta = new TH1D("hmueta",";#eta^{#mu};",300,-4,4);
  TH1D* hmumuy = new TH1D("hmumuy",";y^{#mu#mu};",300,-4,4);
  TH1D* hVtxProb = new TH1D("hVtxProb",";Reco_QQ_VtxProb;",200,0,1);
  
  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();
  cout << "Total events = " << mytree->GetEntries() << endl;

  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);

    if(Reco_mu_size==0) continue;

    nMuCand = 0;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      MuMu_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      MuPl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      MuMi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      bool passGlbMuPl = ((Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)))>0) ? true : false;
      bool passTrkMuPl = ((Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)))>0) ? true : false;
      bool passGlbMuMi = ((Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)))>0) ? true : false;
      bool passTrkMuMi = ((Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)))>0) ? true : false;

      bool passMuonTypeMuPl = true;
      bool passMuonTypeMuMi = true;

      if(kMuId == kMuGlbTrk){
        passMuonTypeMuPl = passMuonTypeMuPl && passGlbMuPl && passTrkMuPl;
        passMuonTypeMuMi = passMuonTypeMuMi && passGlbMuMi && passTrkMuMi;
      }
      if(kMuId == kMuGlb){
        passMuonTypeMuPl = passMuonTypeMuPl && passGlbMuPl;
        passMuonTypeMuMi = passMuonTypeMuMi && passGlbMuMi;
      }
      if((kMuId!=kMuGlbTrk) && (kMuId!=kMuGlb)){
        cout << "ERROR!!!! no passing muon type selected!!" << endl;
        break;
      }

      bool muPlSoft = ( 
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) &&
          passMuonTypeMuPl      
          ) ; 
      bool muMiSoft = ( 
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.) &&
          passMuonTypeMuMi      
          ) ; 
      if ( !(muPlSoft && muMiSoft) ) 
        continue;   
      
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

      
      bool massPassFlag = false;
      if(isSignal){
        if( MuMu_Reco->M() > 2.9 && MuMu_Reco->M() < 3.2) massPassFlag = true;
      }
      else if(!isSignal){
        if( (MuMu_Reco->M() > 2.6 && MuMu_Reco->M() < 2.8) || (MuMu_Reco->M() > 3.3 && MuMu_Reco->M() < 3.5) ) massPassFlag = true;
      }

      if(!massPassFlag) continue;

      nMuCand++; 
      
      hNormChi2          -> Fill(Reco_mu_normChi2_global[Reco_QQ_mumi_idx[irqq]]);
      hNormChi2          -> Fill(Reco_mu_normChi2_global[Reco_QQ_mupl_idx[irqq]]);
      hnTrkHits          -> Fill(Reco_mu_nTrkHits[Reco_QQ_mumi_idx[irqq]]);
      hnTrkHits          -> Fill(Reco_mu_nTrkHits[Reco_QQ_mupl_idx[irqq]]);
      hnMuValHits        -> Fill(Reco_mu_nMuValHits[Reco_QQ_mumi_idx[irqq]]);
      hnMuValHits        -> Fill(Reco_mu_nMuValHits[Reco_QQ_mupl_idx[irqq]]);
      hmudxy             -> Fill(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]]);
      hmudxy             -> Fill(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]]);
      hmudxyErr          -> Fill(Reco_mu_dxyErr[Reco_QQ_mumi_idx[irqq]]);
      hmudxyErr          -> Fill(Reco_mu_dxyErr[Reco_QQ_mupl_idx[irqq]]);
      hmudz              -> Fill(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]]);
      hmudz              -> Fill(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]]);
      hmudzErr           -> Fill(Reco_mu_dzErr[Reco_QQ_mumi_idx[irqq]]);
      hmudzErr           -> Fill(Reco_mu_dzErr[Reco_QQ_mupl_idx[irqq]]);
      hnTrkWMea          -> Fill(Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]]);
      hnTrkWMea          -> Fill(Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]]);
      hmuTMOneStaTight   -> Fill(Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]);
      hmuTMOneStaTight   -> Fill(Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]);
      hmunPixWMea        -> Fill(Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]]);
      hmunPixWMea        -> Fill(Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]]);
      hmunPixValHits     -> Fill(Reco_mu_nPixValHits[Reco_QQ_mumi_idx[irqq]]);
      hmunPixValHits     -> Fill(Reco_mu_nPixValHits[Reco_QQ_mupl_idx[irqq]]);
      hmuptErrglobal     -> Fill(Reco_mu_ptErr_global[Reco_QQ_mumi_idx[irqq]]);
      hmuptErrglobal     -> Fill(Reco_mu_ptErr_global[Reco_QQ_mupl_idx[irqq]]);

      hVtxProb           -> Fill(Reco_QQ_VtxProb[irqq]);

      hmupt              -> Fill(MuPl_Reco->Pt());
      hmupt              -> Fill(MuMi_Reco->Pt());
      hmuphi             -> Fill(MuPl_Reco->Phi());
      hmuphi             -> Fill(MuMi_Reco->Phi());
      hmueta             -> Fill(MuPl_Reco->Eta());
      hmueta             -> Fill(MuMi_Reco->Eta());
      
      hmumupt            -> Fill(MuMu_Reco->Pt());
      hmumuy             -> Fill(MuMu_Reco->Rapidity());
    
      hCent              -> Fill(Centrality);
      
      
    } // end of dimuon loop

    if(nMuCand!=0) newMytree->Fill(); 
  } //end of event loop

  newfile->cd();
  newMytree->Write();
  hNormChi2       -> Write(); 
  hNormChi2       -> Write(); 
  hnTrkHits       -> Write(); 
  hnTrkHits       -> Write(); 
  hnMuValHits     -> Write(); 
  hnMuValHits     -> Write(); 
  hmudxy          -> Write(); 
  hmudxy          -> Write(); 
  hmudxyErr       -> Write(); 
  hmudxyErr       -> Write(); 
  hmudz           -> Write(); 
  hmudz           -> Write(); 
  hmudzErr        -> Write(); 
  hmudzErr        -> Write(); 
  hnTrkWMea       -> Write(); 
  hnTrkWMea       -> Write(); 
  hmuTMOneStaTight-> Write(); 
  hmuTMOneStaTight-> Write(); 
  hmunPixWMea     -> Write(); 
  hmunPixWMea     -> Write(); 
  hmunPixValHits  -> Write(); 
  hmunPixValHits  -> Write(); 
  hmuptErrglobal  -> Write(); 
  hmuptErrglobal  -> Write(); 
  hVtxProb        -> Write(); 
  hmupt           -> Write(); 
  hmupt           -> Write(); 
  hmuphi          -> Write(); 
  hmuphi          -> Write(); 
  hmueta          -> Write(); 
  hmueta          -> Write(); 

  hmumupt         -> Write(); 
  hmumuy          -> Write(); 
  hCent           -> Write(); 
  newMytree->Close();
}
