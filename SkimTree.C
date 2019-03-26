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

  TString fname1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part*.root";
  TString fname2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
  
  TChain *mytree = new TChain("hionia/myTree");
  mytree->Add(fname1.Data());
  mytree->Add(fname2.Data());

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


  TChain *eptree = new TChain("vnanalyzer/tree");
  eptree->Add(fname1.Data());
  eptree->Add(fname2.Data());
  
  
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
  
  //Make RooDataSet
  RooRealVar* massVar  = new RooRealVar("mass","mass variable",0,200,"GeV/c^{2}");
  RooRealVar* ptVar    = new RooRealVar("pt","pt variable", 0,100,"GeV/c");
  RooRealVar* yVar     = new RooRealVar("y","rapidity of the dimuon pair", -5,5,"");
  RooRealVar* pt1Var   = new RooRealVar("pt1","pt of muon+", 0,500,"GeV/c");
  RooRealVar* eta1Var  = new RooRealVar("eta1","eta of muon+", -4,4,"");
  RooRealVar* pt2Var   = (RooRealVar*)pt1Var->Clone("pt2");
  RooRealVar* eta2Var  = (RooRealVar*)eta1Var->Clone("eta2");
  RooRealVar* cBinVar   = new RooRealVar("cBin","Centrality bin", -100,500,"");
  RooRealVar* ep2Var   = new RooRealVar("ep2","2nd order event plane", -100,100,"");
  RooRealVar* dphiEp2Var   = new RooRealVar("dphiEp2","Delta Phi from 2nd order event plane", -100,100,"");
  RooRealVar* evtWeight = new RooRealVar("weight","pt weight", 0, 10000,"");
  RooArgSet* argSet    = new RooArgSet(*massVar, *ptVar, *yVar, *pt1Var, *pt2Var, *eta1Var, *eta2Var,*evtWeight);
  argSet->add(*cBinVar);argSet->add(*ep2Var); argSet->add(*dphiEp2Var);
  
  RooDataSet* dataSet  = new RooDataSet("dataset", " a dataset", *argSet);



  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;

  
  int kTrigSel = 13; // 13: for Upsilon Trigger
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
      
      if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;
      
      bool passMuonTypePl = true;
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,1)));
      passMuonTypePl = passMuonTypePl && (Reco_mu_SelectionType[Reco_QQ_mupl_idx[irqq]]&((int)pow(2,3)));

      bool passMuonTypeMi = true;
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,1)));
      passMuonTypeMi = passMuonTypeMi && (Reco_mu_SelectionType[Reco_QQ_mumi_idx[irqq]]&((int)pow(2,3)));

      bool muplSoft = ( (Reco_mu_TMOneStaTight[Reco_QQ_mupl_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mupl_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mupl_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mupl_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mupl_idx[irqq]])<20.) &&
          passMuonTypePl        //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 

      bool mumiSoft = ( (Reco_mu_TMOneStaTight[Reco_QQ_mumi_idx[irqq]]==true) &&
          (Reco_mu_nTrkWMea[Reco_QQ_mumi_idx[irqq]] > 5) &&
          (Reco_mu_nPixWMea[Reco_QQ_mumi_idx[irqq]] > 0) &&
          (fabs(Reco_mu_dxy[Reco_QQ_mumi_idx[irqq]])<0.3) &&
          (fabs(Reco_mu_dz[Reco_QQ_mumi_idx[irqq]])<20.)  && 
          passMuonTypeMi       //			 &&  (Reco_mu_highPurity[Reco_QQ_mupl_idx[irqq]]==true) 
          ) ; 


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
      
      
      //Fill RooDataSet
      massVar->setVal( (double)dm.mass ) ;
      ptVar->setVal(   (double)dm.pt   ) ;
      yVar->setVal(    (double)dm.y    ) ;
      pt1Var->setVal(  (double)dm.pt1  ) ;
      dphiEp2Var->setVal(   (double)dm.dphiEp2  ) ;
      eta1Var->setVal( (double)dm.eta1 ) ;
      pt2Var->setVal(  (double)dm.pt2  ) ;
      eta2Var->setVal( (double)dm.eta2 ) ;
      ep2Var->setVal( (double)dm.ep2 ) ;
      cBinVar->setVal( (double)dm.cBin ) ;
      evtWeight->setVal( (double)dm.weight ) ;
      dataSet->add( *argSet);
      
    } // end of dimuon loop
    
    
  } //end of event loop
  mmtree->Write();  // Don't need to call Write() for trees
  dataSet->Write();
  newfile->Close();
  
} 

