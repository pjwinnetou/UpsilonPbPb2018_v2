#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"
#include "SkimTree_MuonID.h"

static const long MAXTREESIZE = 1000000000000;

void SkimTree_MuonID(int nevt=-1, bool isMC = true, bool isJPsiTrig = true) 
{
  
  TString fMCstr = (isMC) ? "MC" : "DATA" ;

  using namespace std;
  using namespace hi;
  bool isSignal;


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
  
  const int nHist = 18;
  const char* histName = {"NormChi2global","NormChi2Inner","nTrkHits","nMuValHits","mudxy","mudxyErr","mudz","mudzErr","nTrkWMea","muStationsMatched","muTMOneStaTight","munPixWMea","munPixValHits","muptErrglobal","muptErrInner","mupt","muphi","mueta"};
  const int nbins[nhist] = {1000, 1000, 40, 55, 100, 500, 1000, 500, 20, 30, 2, 7, 12, 1000, 1000, 250, 300, 300};

  const double lowbinval[nhist]  = {0,0,0,0,-0.3,0,-5,0,0,0,0,0,0,0,0,0,-4,-4};
  const double highbinval[nhist] = {100,100,40,55,0.3,0.1,5,0.3,20,2,7,12,1,1,50,4,4};
  
  const int nHistdimu = 3;
  const char* histNamedimu = {"mumuPt","mumuY","VtxProb"};
  const int nbinsdimu[nhistdimu] = {250, 300, 1000};

  const double lowbinvaldimu[nhist]  = {0,-4,0};
  const double highbinvaldimu[nhist] = {50,4,1};
  
  const int nHistCut = 6;
  const char* histNameCut[nHistCut] = {"mudxy","mudz","nTrkWMea","muTMOneStaTight","munPixWMea","VtxProb"};
  const int nbinsCut[nHistCut] = {100,1000,20,2,7,200};

  const double lowbinvalCut[nHistCut]  = {-0.3,-5,0,0,0,0};
  const double highbinvalCut[nHistCut] = {0.3,5,20,2,7,1};

  map<TString, TH1D*> hAll, hSig, hBkg, hSigSCut, hBkgSCut, hGlbCut, hGlbTrkCut, hSigGlbCut, hSigGlbTrkCut, hBkgGlbCut, hBkgGlbTrkCut, hSigGlbSCut, hSigGlbTrkSCut, hBkgGlbSCut, hBkgGlbTrkSCut;
  
  TObjArray* fAll = new TObjArray();
  TObjArray* fSig = new TObjArray();
  TObjArray* fBkg = new TObjArray();
  TObjArray* fSigSCut = new TObjArray();
  TObjArray* fBkgSCut = new TObjArray();
  TObjArray* fGlbCut = new TObjArray();
  TObjArray* fGlbTrkCut = new TObjArray();
  TObjArray* fSigGlbCut = new TObjArray();
  TObjArray* fSigGlbTrkCut = new TObjArray();
  TObjArray* fBkgGlbCut = new TObjArray();
  TObjArray* fBkgGlbTrkCut = new TObjArray(); 
  
  TObjArray* fSigGlbSCut = new TObjArray();
  TObjArray* fSigGlbTrkSCut = new TObjArray();
  TObjArray* fBkgGlbSCut = new TObjArray();
  TObjArray* fBkgGlbTrkSCut = new TObjArray(); 
  
  for(int i=0; i<nHist; i++){
    hAll[Form("%s",histName[i])] = new TH1D(Form("h%s_all",histName[i]),Form(";%s;",histName[i]), nbins[i],lowBinVal[i],highBinVal[i]);
    hSig[Form("%s",histName[i])] = new TH1D(Form("h%s_Sig",histName[i]),Form(";%s;",histName[i]), nbins[i],lowBinVal[i],highBinVal[i]);
    hBkg[Form("%s",histName[i])] = new TH1D(Form("h%s_Bkg",histName[i]),Form(";%s;",histName[i]), nbins[i],lowBinVal[i],highBinVal[i]);
    hSigGlbCut[Form("%s",histName[i])] = new TH1D(Form("h%s_SigGlbCut",histName[i]),Form(";%s;",histName[i]), nbins[i],lowBinVal[i],highBinVal[i]);
    hBkgGlbCut[Form("%s",histName[i])] = new TH1D(Form("h%s_BkgGlbCut",histName[i]),Form(";%s;",histName[i]), nbins[i],lowBinVal[i],highBinVal[i]);
    hSigGlbTrkCut[Form("%s",histName[i])] = new TH1D(Form("h%s_SigGlbTrkCut",histName[i]),Form(";%s;",histName[i]), nbins[i],lowBinVal[i],highBinVal[i]);
    hBkgGlbTrkCut[Form("%s",histName[i])] = new TH1D(Form("h%s_BkgGlbTrkCut",histName[i]),Form(";%s;",histName[i]), nbins[i],lowBinVal[i],highBinVal[i]);
  }

  for(int i=0; i<nHistdimu; i++){
    hAll[Form("%s",histNamedimu[i])] = new TH1D(Form("h%s_all",histNamedimu[i]),Form(";%s;",histNamedimu[i]), nbinsdimu[i],lowBinValdimu[i],highBinValdimu[i]);
    hSig[Form("%s",histNamedimu[i])] = new TH1D(Form("h%s_Sig",histNamedimu[i]),Form(";%s;",histNamedimu[i]), nbinsdimu[i],lowBinValdimu[i],highBinValdimu[i]);
    hBkg[Form("%s",histNamedimu[i])] = new TH1D(Form("h%s_Bkg",histNamedimu[i]),Form(";%s;",histNamedimu[i]), nbinsdimu[i],lowBinValdimu[i],highBinValdimu[i]);
    hSigGlbCut[Form("%s",histNamedimu[i])] = new TH1D(Form("h%s_SigGlbCut",histNamedimu[i]),Form(";%s;",histNamedimu[i]), nbinsdimu[i],lowBinValdimu[i],highBinValdimu[i]);
    hBkgGlbCut[Form("%s",histNamedimu[i])] = new TH1D(Form("h%s_BkgGlbCut",histNamedimu[i]),Form(";%s;",histNamedimu[i]), nbinsdimu[i],lowBinValdimu[i],highBinValdimu[i]);
    hSigGlbTrkCut[Form("%s",histNamedimu[i])] = new TH1D(Form("h%s_SigGlbTrkCut",histNamedimu[i]),Form(";%s;",histNamedimu[i]), nbinsdimu[i],lowBinValdimu[i],highBinValdimu[i]);
    hBkgGlbTrkCut[Form("%s",histNamedimu[i])] = new TH1D(Form("h%s_BkgGlbTrkCut",histNamedimu[i]),Form(";%s;",histNamedimu[i]), nbinsdimu[i],lowBinValdimu[i],highBinValdimu[i]);
  }

  for(int i=0; i<nHistCut; i++){
    hSigSCut[Form("%s",histNameCut[i])] = new TH1D(Form("h%s_SigSCut",histNameCut[i]),Form(";%s;",histNameCut[i]), nbinsCut[i],lowBinValCut[i],highBinValCut[i]);
    hBkgSCut[Form("%s",histNameCut[i])] = new TH1D(Form("h%s_BkgSCut",histNameCut[i]),Form(";%s;",histNameCut[i]), nbinsCut[i],lowBinValCut[i],highBinValCut[i]);
    hSigGlbSCut[Form("%s",histNameCut[i])] = new TH1D(Form("h%s_SigGlbSCut",histNameCut[i]),Form(";%s;",histNameCut[i]), nbinsCut[i],lowBinValCut[i],highBinValCut[i]);
    hBkgGlbSCut[Form("%s",histNameCut[i])] = new TH1D(Form("h%s_BkgGlbSCut",histNameCut[i]),Form(";%s;",histNameCut[i]), nbinsCut[i],lowBinValCut[i],highBinValCut[i]);
    hSigGlbTrkSCut[Form("%s",histNameCut[i])] = new TH1D(Form("h%s_SigGlbTrkSCut",histNameCut[i]),Form(";%s;",histNameCut[i]), nbinsCut[i],lowBinValCut[i],highBinValCut[i]);
    hBkgGlbTrkSCut[Form("%s",histNameCut[i])] = new TH1D(Form("h%s_BkgGlbTrkSCut",histNameCut[i]),Form(";%s;",histNameCut[i]), nbinsCut[i],lowBinValCut[i],highBinValCut[i]);
  }


  int kTrigSel = 12;
  // event loop start
  if(nevt == -1) nevt = 100000;//mytree->GetEntries();
  cout << "Total events = " << mytree->GetEntries() << endl;

  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);

    if(Reco_mu_size==0) continue;

    nMuCand = 0;
    if(isJpsiTrig && (!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) )) continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      MuMu_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      MuPl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      MuMi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      
      if(Reco_QQ_sign==0 && (MuMu_Reco->M() < 2.6 || MuMu_Reco->M()>3.5)) continue;
      
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
      
      passMuonTypeGlbTrkMuPl = passMuonTypeGlbTrkMuPl && passGlbMuPl && passTrkMuPl;
      passMuonTypeGlbTrkMuMi = passMuonTypeGlbTrkMuMi && passGlbMuMi && passTrkMuMi;
      passMuonTypeGlbMuPl = passMuonTypeGlbMuPl && passGlbMuPl && passTrkMuPl;
      passMuonTypeGlbMuMi = passMuonTypeGlbMuMi && passGlbMuMi && passTrkMuMi;

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

      if( MuMu_Reco->M() > 2.9 && MuMu_Reco->M() < 3.2) isSignal = true;
      else if( (MuMu_Reco->M() > 2.6 && MuMu_Reco->M() < 2.8) || (MuMu_Reco->M() > 3.3 && MuMu_Reco->M() < 3.5) ) isSignal = false;

      double weight = 1;
      if(isMC) weight = findNcoll(Centrality);

      double histFillVarM[nHist] = {Reco_mu_normChi2_global[Reco_QQ_mumi_idx[irqq]], Reco_mu_normChi2_inner[Reco_QQ_mumi_idx[irqq]], Reco_mu_nTrkHits[Reco_QQ_mumi_idx[irqq]], Reco_mu_nMuValHits[Reco_QQ_mumi_idx[irqq]], Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]], Reco_mu_dxyErr[Reco_QQ_mumi_idx[irqq]], Reco_mu_dz[Reco_QQ_mumi_idx[irqq]], Reco_mu_dzErr[Reco_QQ_mumi_idx[irqq]], Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]], Reco_mu_StationsMatched[Reco_QQ_mumi_idx[irqq]], Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]], Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]], Reco_mu_nPixValHits[Reco_QQ_mumi_idx[irqq]], Reco_mu_ptErr_global[Reco_QQ_mumi_idx[irqq]], Reco_mu_ptErr_inner[Reco_QQ_mumi_idx[irqq]], MuMi_Reco->Pt(), MuMi_Reco->Phi(), MuMi_Reco->Eta()};
      double histFillVarP[nHist] = {Reco_mu_normChi2_global[Reco_QQ_mupl_idx[irqq]], Reco_mu_normChi2_inner[Reco_QQ_mupl_idx[irqq]], Reco_mu_nTrkHits[Reco_QQ_mupl_idx[irqq]], Reco_mu_nMuValHits[Reco_QQ_mupl_idx[irqq]], Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]], Reco_mu_dxyErr[Reco_QQ_mupl_idx[irqq]], Reco_mu_dz[Reco_QQ_mupl_idx[irqq]], Reco_mu_dzErr[Reco_QQ_mupl_idx[irqq]], Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]], Reco_mu_StationsMatched[Reco_QQ_mupl_idx[irqq]], Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]], Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]], Reco_mu_nPixValHits[Reco_QQ_mupl_idx[irqq]], Reco_mu_ptErr_global[Reco_QQ_mupl_idx[irqq]], Reco_mu_ptErr_inner[Reco_QQ_mupl_idx[irqq]], muPl_Reco->Pt(), muPl_Reco->Phi(), muPl_Reco->Eta()};
      
      double histFillVardimu[nHistdimu] = {MuMu_Reco->Pt(), MuMu_Reco->Rapidity(), Reco_QQ_VtxProb[irqq]};
      
      double histFillVarCutM[nHistCut] = {Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]], Reco_mu_dz[Reco_QQ_mumi_idx[irqq]], Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]], Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]], Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]], Reco_QQ_VtxProb[irqq]};
      double histFillVarCutP[nHistCut] = {Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]], Reco_mu_dz[Reco_QQ_mupl_idx[irqq]], Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]], Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]], Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]],-999};
      

      for(int ihist = 0; ihist<nHist; ihist++){
        hAll[histName[ihist]]->Fill(histFillVarM[ihist],weight);
        hAll[histName[ihist]]->Fill(histFillVarP[ihist],weight);
        if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl){
          hGlbCut[histName[ihist]]->Fill(histFillVarM[ihist],weight);
          hGlbCut[histName[ihist]]->Fill(histFillVarP[ihist],weight);
        }
        if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl){ 
          hGlbTrkCut[histName[ihist]]->Fill(histFillVarM[ihist],weight);
          hGlbTrkCut[histName[ihist]]->Fill(histFillVarP[ihist],weight);
        }
        if(isSignal){
          hSig[histName[ihist]]->Fill(histFillVarM[ihist],weight);
          hSig[histName[ihist]]->Fill(histFillVarP[ihist],weight);
          if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl){ 
            hSigGlbCut[histName[ihist]]->Fill(histFillVarM[ihist],weight);
            hSigGlbCut[histName[ihist]]->Fill(histFillVarP[ihist],weight);
          }
          if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl){ 
            hSigGlbTrkCut[histName[ihist]]->Fill(histFillVarM[ihist],weight);
            hSigGlbTrkCut[histName[ihist]]->Fill(histFillVarP[ihist],weight);
          }
        }
        else if(!isSignal){
          hBkg[histName[ihist]]->Fill(histFillVarM[ihist],weight);
          hBkg[histName[ihist]]->Fill(histFillVarP[ihist],weight);
          if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl){ 
            hBkgGlbCut[histName[ihist]]->Fill(histFillVarM[ihist],weight);
            hBkgGlbCut[histName[ihist]]->Fill(histFillVarP[ihist],weight);
          }
          if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl){ 
            hBkgGlbTrkCut[histName[ihist]]->Fill(histFillVarM[ihist],weight);
            hBkgGlbTrkCut[histName[ihist]]->Fill(histFillVarP[ihist],weight);
          }
        }
      }
      for(int ihist = 0; ihist<nHistdimu; ihist++){
        hAll[histNamedimu[ihist]]->Fill(histFillVardimu[ihist],weight);
        if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl) hGlbCut[histNamedimu[ihist]]->Fill(histFillVardimu[ihist],weight);
        if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl) hGlbTrkCut[histNamedimu[ihist]]->Fill(histFillVardimu[ihist],weight);
        if(isSignal){
          hSig[histNamedimu[ihist]]->Fill(histFillVardimu[ihist],weight);
          if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl) hSigGlbCut[histNamedimu[ihist]]->Fill(histFillVardimu[ihist],weight);
          if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl) hSigGlbTrkCut[histNamedimu[ihist]]->Fill(histFillVardimu[ihist],weight);
        }
        else if(!isSignal){
          hBkg[histNamedimu[ihist]]->Fill(histFillVardimu[ihist],weight);
          if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl) hBkgGlbCut[histNamedimu[ihist]]->Fill(histFillVardimu[ihist],weight);
          if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl) hBkgGlbTrkCut[histNamedimu[ihist]]->Fill(histFillVardimu[ihist],weight);
        }
      }
      for(int ihist = 0; ihist<nHistCut; ihist++){
        if(isSignal){
          if(VarCut(histNameCut[ihist],irqq) && ihist<5){
            hSigSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
            hSigSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
            if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl){ 
              hSigGlbSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
              hSigGlbSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
            }
            if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl){ 
              hSigGlbTrkSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
              hSigGlbTrkSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
            }
          }
          else if(ihist==5){
            if(Reco_QQ_VtxProb[irqq]>0.01){
              hSigSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
              hSigSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
              if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl){ 
                hSigGlbSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
                hSigGlbSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
              }
              if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl){ 
                hSigGlbTrkSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
                hSigGlbTrkSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
              }
            }
          }
        }
        else if(!isSignal){
          if(VarCut(histNameCut[ihist],irqq) && ihist<5){
            hSigSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
            hSigSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
            if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl){ 
              hSigGlbSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
              hSigGlbSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
            }
            if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl){ 
              hSigGlbTrkSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
              hSigGlbTrkSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
            }
          }
          else if(ihist==5){
            if(Reco_QQ_VtxProb[irqq]>0.01){
              hBkgSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
              hBkgSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
              if(passMuonTypeGlbMuMi && passMuonTypeGlbMuPl){ 
                hBkgGlbSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
                hBkgGlbSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
              }
              if(passMuonTypeGlbTrkMuMi && passMuonTypeGlbTrkMuPl){ 
                hBkgGlbTrkSCut[histNameCut[ihist]]->Fill(histFillVarCutM[ihist],weight);
                hBkgGlbTrkSCut[histNameCut[ihist]]->Fill(histFillVarCutP[ihist],weight);
              }
            }
          }
        }
      }

      nMuCand++; 
      
    } // end of dimuon loop

    if(nMuCand!=0) newMytree->Fill(); 
  } //end of event loop

  for(int ihist = 0; i<nHist; ihist++){
    fAll            -> Add(hAll[histName[ihist]]);  
    fSig            -> Add(hSig[histName[ihist]]); 
    fBkg            -> Add(hBkg[histName[ihist]]);  
    fGlbCut         -> Add(hGlbCut[histName[ihist]]);  
    fGlbTrkCut      -> Add(hGlbTrkCut[histName[ihist]]);  
    fSigGlbCut      -> Add(hSigGlbCut[histName[ihist]]);  
    fSigGlbTrkCut   -> Add(hSigGlbTrkCut[histName[ihist]]);
    fBkgGlbCut      -> Add(hBkgGlbCut[histName[ihist]]);  
    fBkgGlbTrkCut   -> Add(hBkgGlbTrkCut[histName[ihist]]);  
  }
  for(int ihist = 0; i<nHistdimu; ihist++){
    fAll            -> Add(hAll[histNamedimu[ihist]]);  
    fSig            -> Add(hSig[histNamedimu[ihist]]); 
    fBkg            -> Add(hBkg[histNamedimu[ihist]]);  
    fGlbCut         -> Add(hGlbCut[histNamedimu[ihist]]);  
    fGlbTrkCut      -> Add(hGlbTrkCut[histNamedimu[ihist]]);  
    fSigGlbCut      -> Add(hSigGlbCut[histNamedimu[ihist]]);  
    fSigGlbTrkCut   -> Add(hSigGlbTrkCut[histNamedimu[ihist]]);
    fBkgGlbCut      -> Add(hBkgGlbCut[histNamedimu[ihist]]);  
    fBkgGlbTrkCut   -> Add(hBkgGlbTrkCut[histNamedimu[ihist]]);  
  }
  for(int ihist = 0; i<nHistCut; ihist++){
    fSigSCut        -> Add(hSigSCut[histNameCut[ihist]]);  
    fBkgSCut        -> Add(hBlgSCut[histNameCut[ihist]]);
    fSigGlbSCut     -> Add(hSigGlbSCut[histNameCut[ihist]]);  
    fSigGlbTrkSCut  -> Add(hSigGlbTrkSCut[histNameCut[ihist]]); 
    fBkgGlbSCut     -> Add(hBkgGlbSCut[histNameCut[ihist]]); 
    fBkgGlbTrkSCut  -> Add(hBkgGlbTrkSCut[histNameCut[ihist]]); 
  }
  newfile->cd();
  newMytree->Write();
  fAll            -> Write("AllHist", TObject::kOverwrite | TObject::kSingleKey);
  fSig            -> Write("SigHist", TObject::kOverwrite | TObject::kSingleKey);
  fBkg            -> Write("BkgHist", TObject::kOverwrite | TObject::kSingleKey);
  fGlbCut         -> Write("GlbHist", TObject::kOverwrite | TObject::kSingleKey);
  fGlbTrkCut      -> Write("GlbTrkHist", TObject::kOverwrite | TObject::kSingleKey);
  fSigGlbCut      -> Write("SigGlbHist", TObject::kOverwrite | TObject::kSingleKey);
  fSigGlbTrkCut   -> Write("SigGlbTrkHist", TObject::kOverwrite | TObject::kSingleKey);
  fBkgGlbCut      -> Write("BkgGlbHist", TObject::kOverwrite | TObject::kSingleKey);
  fBkgGlbTrkCut   -> Write("BkgGlbTrkHist", TObject::kOverwrite | TObject::kSingleKey);
  fSigSCut        -> Write("SigSHist", TObject::kOverwrite | TObject::kSingleKey);
  fBkgSCut        -> Write("BkgSHist", TObject::kOverwrite | TObject::kSingleKey);
  fSigGlbSCut     -> Write("SigGlbSHist", TObject::kOverwrite | TObject::kSingleKey);
  fSigGlbTrkSCut  -> Write("SigGlbTrkSHist", TObject::kOverwrite | TObject::kSingleKey);
  fBkgGlbSCut     -> Write("BkgGlbSHist", TObject::kOverwrite | TObject::kSingleKey);
  fBkgGlbTrkSCut  -> Write("BkgGlbTrkSHist", TObject::kOverwrite | TObject::kSingleKey);
  newfile->Close();
}


Int_t VarCut(char histName, int irqq){
  if(histName == "mudxy"){
    return (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3 && fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3);
  }
  if(histName == "mudz"){
    return (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20 && fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20);
  }
  if(histName == "nTrkWMea"){
    return (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5 && Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5);
  }
  if(histName == "muTMOneStaTight"){
    return ( (Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) && (Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) );
  }
  if(histName == "munPixWMea"){
    return (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0 && Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0);
  }
  return false;
}


