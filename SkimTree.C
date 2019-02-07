#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"

static const long MAXTREESIZE = 10000000000;

void SkimTree(int nevt=-1) 
{

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  TString fname = "/eos/cms/store/group/phys_heavyions/dileptons/goni/HIDoubleMu_Run2018A_PromptAOD_v1_OniaTree_v2Ana_Cert_Run_326381_327560_48per_190130.root";
  
  TChain *mytree = new TChain("hionia/myTree");
  mytree->Add(fname.Data());


  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Int_t           Reco_QQ_size;
  TClonesArray    *Reco_QQ_4mom;
  TClonesArray    *Reco_QQ_mupl_4mom;
  TClonesArray    *Reco_QQ_mumi_4mom;
  ULong64_t       Reco_QQ_trig[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_VtxProb[200];   //[Reco_QQ_size]
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;
  TBranch        *b_zVtx;   //!
  TBranch        *b_Centrality;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_Reco_QQ_size;   //!
  TBranch        *b_Reco_QQ_4mom;   //!
  TBranch        *b_Reco_QQ_mupl_4mom;   //!
  TBranch        *b_Reco_QQ_mumi_4mom;   //!
  TBranch        *b_Reco_QQ_trig;   //!
  TBranch        *b_Reco_QQ_VtxProb;   //!

  Bool_t          Reco_QQ_mupl_highPurity[200];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_mumi_highPurity[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_highPurity;   //!
  TBranch        *b_Reco_QQ_mumi_highPurity;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_highPurity", Reco_QQ_mupl_highPurity, &b_Reco_QQ_mupl_highPurity);
  mytree->SetBranchAddress("Reco_QQ_mumi_highPurity", Reco_QQ_mumi_highPurity, &b_Reco_QQ_mumi_highPurity);

  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;
  mytree->SetBranchAddress("runNb", &runNb, &b_runNb);
  mytree->SetBranchAddress("LS", &LS, &b_LS);
  mytree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  mytree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  mytree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  mytree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  mytree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  mytree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  //  mytree->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(kFALSE);
  mytree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  //  mytree->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(kFALSE);
  mytree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  mytree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  mytree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);

  //  muon id 
  Int_t           Reco_QQ_mupl_nTrkHits[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nTrkHits[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nTrkHits[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nTrkHits;   //!
  TBranch        *b_Reco_QQ_mumi_nTrkHits;   //!
  TBranch        *b_Reco_mu_nTrkHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
  mytree->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
  mytree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  Float_t         Reco_QQ_mupl_normChi2_global[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_normChi2_global[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_normChi2_global[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_normChi2_global;   //!
  TBranch        *b_Reco_QQ_mumi_normChi2_global;   //!
  TBranch        *b_Reco_mu_normChi2_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
  mytree->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
  mytree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  Int_t           Reco_QQ_mupl_nMuValHits[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nMuValHits[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nMuValHits[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nMuValHits;   //!
  TBranch        *b_Reco_QQ_mumi_nMuValHits;   //!
  TBranch        *b_Reco_mu_nMuValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
  mytree->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
  mytree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  Int_t           Reco_QQ_mupl_StationsMatched[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_StationsMatched[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_StationsMatched[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_StationsMatched;   //!
  TBranch        *b_Reco_QQ_mumi_StationsMatched;   //!
  TBranch        *b_Reco_mu_StationsMatched;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
  mytree->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
  mytree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  Float_t         Reco_QQ_mupl_dxy[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxy[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dxyErr[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dxyErr[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_dxy[200];   //[Reco_mu_size]
  Float_t         Reco_mu_dxyErr[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_dxy;   //!
  TBranch        *b_Reco_QQ_mumi_dxy;   //!
  TBranch        *b_Reco_QQ_mupl_dxyErr;   //!
  TBranch        *b_Reco_QQ_mumi_dxyErr;   //!
  TBranch        *b_Reco_mu_dxy;   //!
  TBranch        *b_Reco_mu_dxyErr;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  mytree->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  mytree->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
  mytree->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
  mytree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  mytree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  Float_t         Reco_QQ_mupl_dz[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dz[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mupl_dzErr[200];   //[Reco_QQ_size]
  Float_t         Reco_QQ_mumi_dzErr[200];   //[Reco_QQ_size]
  Float_t         Reco_mu_dz[200];   //[Reco_mu_size]
  Float_t         Reco_mu_dzErr[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_dz;   //!
  TBranch        *b_Reco_QQ_mumi_dz;   //!
  TBranch        *b_Reco_QQ_mupl_dzErr;   //!
  TBranch        *b_Reco_QQ_mumi_dzErr;   //!
  TBranch        *b_Reco_mu_dz;   //!
  TBranch        *b_Reco_mu_dzErr;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  mytree->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  mytree->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
  mytree->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
  mytree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  mytree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  Int_t           Reco_QQ_mupl_nTrkWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nTrkWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nTrkWMea[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
  TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
  TBranch        *b_Reco_mu_nTrkWMea;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
  mytree->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
  mytree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  Bool_t          Reco_QQ_mupl_TMOneStaTight[200];   //[Reco_QQ_size]
  Bool_t          Reco_QQ_mumi_TMOneStaTight[200];   //[Reco_QQ_size]
  Bool_t          Reco_mu_TMOneStaTight[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_TMOneStaTight;   //!
  TBranch        *b_Reco_QQ_mumi_TMOneStaTight;   //!
  TBranch        *b_Reco_mu_TMOneStaTight;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
  mytree->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);

  mytree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  Int_t           Reco_QQ_mupl_nPixWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_QQ_mumi_nPixWMea[200];   //[Reco_QQ_size]
  Int_t           Reco_mu_nPixWMea[200];   //[Reco_mu_size]
  TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
  TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
  TBranch        *b_Reco_mu_nPixWMea;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
  mytree->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
  mytree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  Int_t           Reco_QQ_sign[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_sign;   //!
  mytree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  Float_t         rpAng[29];   //[nEP]
  TBranch        *b_rpAng;   //!
  mytree->SetBranchAddress("rpAng", rpAng, &b_rpAng);

  Int_t           Reco_QQ_mupl_nPixValHits[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
  Int_t           Reco_QQ_mumi_nPixValHits[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mumi_nPixValHits;   //!
  mytree->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
  Float_t         Reco_QQ_mupl_ptErr_global[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mupl_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
  Float_t         Reco_QQ_mumi_ptErr_global[200];   //[Reco_QQ_size]
  TBranch        *b_Reco_QQ_mumi_ptErr_global;   //!
  mytree->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);

  ULong64_t           Reco_QQ_mupl_SelectionType[200];
  ULong64_t           Reco_QQ_mumi_SelectionType[200];
  TBranch        *b_Reco_QQ_mupl_SelectionType;
  TBranch        *b_Reco_QQ_mumi_SelectionType;
  mytree->SetBranchAddress("Reco_QQ_mupl_SelectionType", Reco_QQ_mupl_SelectionType, &b_Reco_QQ_mupl_SelectionType);
  mytree->SetBranchAddress("Reco_QQ_mumi_SelectionType", Reco_QQ_mumi_SelectionType, &b_Reco_QQ_mumi_SelectionType);



  TChain *eptree = new TChain("vnanalyzer/tree");
  eptree->Add(fname.Data());
  
  
  const int nEP = 29;  // number of event planes in the tree
  double qx[nEP]; 
  double qy[nEP]; 
  TBranch *b_qx;
  TBranch *b_qy;
  eptree->SetBranchAddress("qx",qx, &b_qx);
  eptree->SetBranchAddress("qy",qy, &b_qy);
  
  TFile* newfile;
  newfile = new TFile("OniaFlow_skim.root","recreate");

  DiMuon dm;
  TTree* mmtree = new TTree("mmep","dimuonAndEventPlanes");
  mmtree->SetMaxTreeSize(MAXTREESIZE);
  mmtree->Branch("mm",&dm.run,branchString.Data());


  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();

  cout << "Total events = " << nevt<<endl;

  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    eptree->GetEntry(iev);
    
    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      dm.clear();      // clear the output tree: 
      dm.run = runNb;
      dm.lumi = LS ;
      dm.event = eventNb ;
      dm.vz = zVtx;
      dm.cBin = Centrality ;

      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_QQ_mupl_4mom->At(irqq);
      mumi_Reco = (TLorentzVector*) Reco_QQ_mumi_4mom->At(irqq);
      double mass = JP_Reco->M();
      double phi = JP_Reco->Phi();
      double eta = JP_Reco->Eta();

      bool muplSoft = ( (Reco_QQ_mupl_TMOneStaTight[irqq]==true) &&
          (Reco_QQ_mupl_nTrkWMea[irqq] > 5) &&
          (Reco_QQ_mupl_nPixWMea[irqq] > 0) &&
          (fabs(Reco_QQ_mupl_dxy[irqq])<0.3) &&
          (fabs(Reco_QQ_mupl_dz[irqq])<20.) &&
          ((Reco_QQ_mupl_SelectionType[irqq]&1)>0)        //			 &&  (Reco_QQ_mupl_highPurity[irqq]==true) 
          ) ; 

      bool mumiSoft = ( (Reco_QQ_mumi_TMOneStaTight[irqq]==true) &&
          (Reco_QQ_mumi_nTrkWMea[irqq] > 5) &&
          (Reco_QQ_mumi_nPixWMea[irqq] > 0) &&
          (fabs(Reco_QQ_mumi_dxy[irqq])<0.3) &&
          (fabs(Reco_QQ_mumi_dz[irqq])<20.)  && 
          ((Reco_QQ_mumi_SelectionType[irqq]&1)>0)        //			 &&  (Reco_QQ_mupl_highPurity[irqq]==true) 
          ) ; 

      bool muplHighPtCut = ( (Reco_QQ_mupl_nMuValHits[irqq]>0) &&
          (Reco_QQ_mupl_StationsMatched[irqq]>1) &&
          (Reco_QQ_mupl_ptErr_global[irqq]/mupl_Reco->Pt() < 0.3 ) && 
          (Reco_QQ_mupl_dxy[irqq]<0.2) &&
          (Reco_QQ_mupl_dz[irqq] <0.5) &&
          (Reco_QQ_mupl_nPixValHits[irqq] > 0 ) &&
          (Reco_QQ_mupl_nTrkWMea[irqq] > 5) 
          );
      bool mumiHighPtCut = ( (Reco_QQ_mumi_nMuValHits[irqq]>0) &&
          (Reco_QQ_mumi_StationsMatched[irqq]>1) &&
          (Reco_QQ_mumi_ptErr_global[irqq]/mumi_Reco->Pt() < 0.3 ) && 
          (Reco_QQ_mumi_dxy[irqq]<0.2) &&
          (Reco_QQ_mumi_dz[irqq] <0.5) &&
          (Reco_QQ_mumi_nPixValHits[irqq] > 0 ) &&
          (Reco_QQ_mumi_nTrkWMea[irqq] > 5) 
          );

      if ( !(muplSoft && mumiSoft) ) 
        continue;   
      
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) 
        continue;
        
      if( Reco_QQ_sign[irqq] != 0) continue;

     /* cout << "Event #: " << iev<<endl;
      cout << "   mass = "<<mass <<endl;
      cout << "   phi = "<<phi<<endl;
      cout << "   eta = "<<eta<<endl;
      cout << "HF plus Q-vector" << endl;
      cout << "   (qx, qy) = ("<<qx[HFp2]<<", "<<qy[HFp2]<<")"<<endl;
      cout << "HF minus Q-vector" << endl;
      cout << "   (qx, qy) = ("<<qx[HFm2]<<", "<<qy[HFm2]<<")"<<endl;
      cout << "Track    Q-vector" << endl;
      cout << "   (qx, qy) = ("<<qx[trackmid2]<<", "<<qy[trackmid2]<<")"<<endl;

      cout <<endl;
     */      
      // Fill the output tree
      if ( eta < 0 )  {  
        dm.qxa = qx[HFp2] ;  
        dm.qya = qy[HFp2] ;  
        dm.qxb = qx[HFm2] ;  
        dm.qyb = qy[HFm2] ;  
      }
      else {
        dm.qxa = qx[HFm2] ;  
        dm.qya = qy[HFm2] ;  
        dm.qxb = qx[HFp2] ;  
        dm.qyb = qy[HFp2] ;  
      }
      
      dm.qxc = qx[trackmid2];
      dm.qyc = qy[trackmid2];
      
      dm.mass   = JP_Reco->M();
      dm.pt     = JP_Reco->Pt();
      dm.phi    = JP_Reco->Phi();
      dm.y      = JP_Reco->Rapidity();
      dm.eta      = JP_Reco->Eta();
      dm.pt1  = mupl_Reco->Pt();
      dm.eta1 = mupl_Reco->Eta();
      dm.phi1 = mupl_Reco->Phi();
      dm.pt2  = mumi_Reco->Pt();
      dm.eta2 = mumi_Reco->Eta();
      dm.phi2 = mumi_Reco->Phi();

      dm.qxdimu = TMath::Cos(2*dm.phi);
      dm.qydimu = TMath::Sin(2*dm.phi);

      mmtree->Fill();
      
    } // end of dimuon loop
    
    
  } //end of event loop
  mmtree->Write();  // Don't need to call Write() for trees
  newfile->Close();
  
} 

