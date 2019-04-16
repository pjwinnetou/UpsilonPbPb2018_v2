#ifndef SkimTree_MuonID_C
#define SkimTree_MuonID_C


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


#endif
