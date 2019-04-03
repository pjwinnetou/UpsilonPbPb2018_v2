#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"

static const long MAXTREESIZE = 1000000000000;

void histMuonIdComp(bool isMC = false, int muSel = 2) 
{
  using namespace std;
  using namespace hi;

  TString strMuSel;
  if(muSel == 1) strMuSel = "Glb";
  else if(muSel == 2) strMuSel = "GlbTrk";

  TString strMC = (isMC) ? "MC" : "DATA";  

  TFile* f1 = new TFile(Form("Onia_Muon_Skim_%s_%s.root",strMC.Data(),strMuSel.Data()),"read");
  TTree *mytree = (TTree*) f1->Get("newTree");
  
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

  Reco_QQ_4mom = 0;
  Reco_mu_4mom = 0;
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  mytree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
//  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

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


  TFile* newfile = new TFile(Form("hist_%s_%s.root",strMC.Data(),strMuSel.Data()),"recreate");

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
  TH1D* hmuphi = new TH1D("hmuphi",";#phi^{#mu};",300,-4,4);
  TH1D* hmueta = new TH1D("hmueta",";#eta^{#mu};",300,-4,4);

  cout << "OK" << endl;

  TLorentzVector* mu_Reco = new TLorentzVector;

  // event loop start
  int nevt = mytree->GetEntries();

  cout << "Total events = " << mytree->GetEntries() << endl;

  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    
    hCent->Fill(Centrality);

    for (Int_t irqq=0; irqq<Reco_mu_size; ++irqq) 
    {

      mu_Reco = (TLorentzVector*) Reco_mu_4mom->At(irqq); 

      hNormChi2->Fill(Reco_mu_normChi2_global[irqq]);
      hnTrkHits->Fill(Reco_mu_nTrkHits[irqq]);
      hnMuValHits->Fill(Reco_mu_nMuValHits[irqq]);
      hmudxy -> Fill(Reco_mu_dxy[irqq]);
      hmudxyErr -> Fill(Reco_mu_dxyErr[irqq]);
      hmudz -> Fill(Reco_mu_dz[irqq]);
      hmudzErr->Fill(Reco_mu_dzErr[irqq]);
      hnTrkWMea->Fill(Reco_mu_nTrkWMea[irqq]);
      hmuTMOneStaTight->Fill(Reco_mu_TMOneStaTight[irqq]);
      hmunPixWMea->Fill(Reco_mu_nPixWMea[irqq]);
      hmunPixValHits->Fill(Reco_mu_nPixValHits[irqq]);
      hmuptErrglobal->Fill(Reco_mu_ptErr_global[irqq]);

      hmupt->Fill(mu_Reco->Pt());
      hmuphi->Fill(mu_Reco->Phi());
      hmueta->Fill(mu_Reco->Eta());

    } // end of dimuon loop
   
  } //end of event loop
 
  newfile->cd();
  
  hNormChi2        -> Write(); 
  hnTrkHits        -> Write();
  hnMuValHits      -> Write();
  hmudxy           -> Write();
  hmudxyErr        -> Write();
  hmudz            -> Write();
  hmudzErr         -> Write();
  hnTrkWMea        -> Write();
  hmuTMOneStaTight -> Write();
  hmunPixWMea      -> Write();
  hmunPixValHits   -> Write();
  hmuptErrglobal   -> Write();
                  
  hmupt            -> Write();
  hmuphi           -> Write();
  hmueta           -> Write();

  newfile->Close();

} 

