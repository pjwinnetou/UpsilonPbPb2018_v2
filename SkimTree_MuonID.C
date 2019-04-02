#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"

static const long MAXTREESIZE = 1000000000000;

void SkimTree_MuonID(int nevt=-1) 
{
  int kMuGlb = 1;
  int kMuGlbTrk = 2;
  int kMuId = kMuGlb;
  TString muIdStr;
  if(kMuId == kMuGlb) muIdStr = "Glb";
  if(kMuId == kMuGlbTrk) muIdStr = "GlbTrk";

  using namespace std;
  using namespace hi;


  TString fname1_1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part1.root";
  TString fname1_2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part2.root";
  TString fname1_3 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part3.root";
  TString fname1_4 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part4.root";
  TString fname1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part*.root";
  TString fname2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
  
  TChain *mytree = new TChain("myTree");
  mytree->Add(fname1.Data());
  mytree->Add(fname2.Data());

  
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

  TFile* newfile;
  newfile = new TFile(Form("Onia_Muon_Skim_%s_all.root",muIdStr.Data()),"recreate");

  TTree* newMytree = new TTree("newTree","skimmedTree");
  newMytree -> Branch("eventNb",&eventNb,"eventNb/i");
  newMytree -> Branch("zVtx",&zVtx,"zVtx/F");
  newMytree -> Branch("Centrality", &Centrality, "Centrality/I");
  newMytree -> Branch("HLTriggers", &HLTriggers, "HLTriggers/l");
  newMytree -> Branch("Reco_QQ_size", &Reco_QQ_size, "Reco_QQ_size/I");
  newMytree -> Branch("Reco_mu_size", &Reco_mu_size, "Reco_mu_size/I");
  newMytree -> Branch("Reco_QQ_4mom", &Reco_QQ_4mom);
  newMytree -> Branch("Reco_mu_4mom", &Reco_mu_4mom);
  newMytree -> Branch("Reco_QQ_trig", Reco_QQ_trig, "Reco_QQ_trig/l");
  newMytree -> Branch("Reco_QQ_VtxProb", Reco_QQ_VtxProb, "Reco_QQ_VtxProb/l");
  newMytree -> Branch("Reco_QQ_mupl_idx",Reco_QQ_mupl_idx, "Reco_QQ_mupl_idx/I");
  newMytree -> Branch("Reco_QQ_mumi_idx",Reco_QQ_mumi_idx, "Reco_QQ_mumi_idx/I");
  newMytree -> Branch("Reco_mu_nTrkHits", Reco_mu_nTrkHits, "Reco_mu_nTrkHits/I");
  newMytree -> Branch("Reco_mu_normChi2_global", Reco_mu_normChi2_global, "Reco_mu_normChi2_global/F");
  newMytree -> Branch("Reco_mu_nMuValHits", Reco_mu_nMuValHits, "Reco_mu_nMuValHits/I");
  newMytree -> Branch("Reco_mu_StationsMatched", Reco_mu_StationsMatched, "Reco_mu_StationsMatched/I");
  newMytree -> Branch("Reco_mu_dxy", Reco_mu_dxy, "Reco_mu_dxy/F");
  newMytree -> Branch("Reco_mu_dxyErr", Reco_mu_dxyErr, "Reco_mu_dxyErr/F");
  newMytree -> Branch("Reco_mu_dz", Reco_mu_dz, "Reco_mu_dz/F");
  newMytree -> Branch("Reco_mu_dzErr", Reco_mu_dzErr, "Reco_mu_dzErr/F");
  newMytree -> Branch("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, "Reco_mu_nTrkWMea/I");
  newMytree -> Branch("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, "Reco_mu_TMOneStaTight/O");
  newMytree -> Branch("Reco_mu_nPixWMea", Reco_mu_nPixWMea, "Reco_mu_nPixWMea/I");
  newMytree -> Branch("Reco_QQ_sign", Reco_QQ_sign, "Reco_QQ_sign/I");
  newMytree -> Branch("Reco_mu_nPixValHits", Reco_mu_nPixValHits, "Reco_mu_nPixValHits/I");
  newMytree -> Branch("Reco_mu_ptErr_global", Reco_mu_ptErr_global, "Reco_mu_ptErr_global/F");
  newMytree -> Branch("Reco_mu_SelectionType", Reco_mu_SelectionType, "Reco_mu_SelectionType/I");


  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();

  cout << "Total events = " << mytree->GetEntries() << endl;

  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);

    if(Reco_mu_size==0) continue;

    int nMuCand = 0;
    for (Int_t irqq=0; irqq<Reco_mu_size; ++irqq) 
    {
      bool passGlb = ((Reco_mu_SelectionType[irqq]&((int)pow(2,1)))>0) ? true : false;
      bool passTrk = ((Reco_mu_SelectionType[irqq]&((int)pow(2,3)))>0) ? true : false;

      bool passMuonType = true;

      if(kMuId == kMuGlbTrk){
        passMuonType = passMuonType && passGlb && passTrk;
      }
      if(kMuId == kMuGlb){
        passMuonType = passMuonType && passGlb;
      }
      if(kMuId!=kMuGlbTrk && kMuId!=kMuGlb){
        cout << "ERROR!!!! no passing muon type selected!!" << endl;
        break;
      }

      bool muSoft = ( 
          (Reco_mu_nTrkWMea[irqq] > 5) &&
          (Reco_mu_nPixWMea[irqq] > 0) &&
          (fabs(Reco_mu_dxy[irqq])<0.3) &&
          (fabs(Reco_mu_dz[irqq])<20.) &&
          passMuonType        
          ) ; 
      if ( !(muSoft) ) 
        continue;   
      
      nMuCand++; 
    } // end of dimuon loop
   
    if(nMuCand!=0) newMytree->Fill(); 
  } //end of event loop
 
  newfile->cd();
  newMytree->Write();
  newfile->Close();
  
} 

