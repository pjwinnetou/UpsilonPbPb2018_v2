#ifndef SkimTree_MuonID_C
#define SkimTree_MuonID_C
#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"


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


Int_t           Reco_QQ_mupl_idx[maxBranchSize];
Int_t           Reco_QQ_mumi_idx[maxBranchSize];
TBranch        *b_Reco_QQ_mupl_idx;
TBranch        *b_Reco_QQ_mumi_idx;

Int_t           Reco_mu_nTrkHits[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_nTrkHits;   //!
Float_t         Reco_mu_normChi2_global[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_normChi2_global;   //!
Float_t         Reco_mu_normChi2_inner[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_normChi2_inner;   //!
Int_t           Reco_mu_nMuValHits[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_nMuValHits;   //!
Int_t           Reco_mu_StationsMatched[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_StationsMatched;   //!
Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
Float_t         Reco_mu_dxyErr[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_dxy;   //!
TBranch        *b_Reco_mu_dxyErr;   //!
Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
Float_t         Reco_mu_dzErr[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_dz;   //!
TBranch        *b_Reco_mu_dzErr;   //!
Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_nTrkWMea;   //!
Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_TMOneStaTight;   //!

Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_mu_nPixWMea;   //!
Int_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_QQ_sign;   //!
Float_t         rpAng[29];   //[nEP]
TBranch        *b_rpAng;   //!
//  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

Int_t           Reco_mu_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_mu_nPixValHits;   //!
Float_t         Reco_mu_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_mu_ptErr_global;   //!
Float_t         Reco_mu_ptErr_inner[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_mu_ptErr_inner;   //!

Int_t           Reco_mu_SelectionType[maxBranchSize];
TBranch        *b_Reco_mu_SelectionType;

#endif
