#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinUpsilonV2.h"

static const long MAXTREESIZE = 1000000000000;
double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);

const Int_t nBins = 200; // table of bin edges
const Double_t binTable[nBins+1] = {0, 10.5072, 11.2099, 11.8364, 12.478, 13.1194, 13.7623, 14.4081, 15.0709, 15.7532, 16.4673, 17.1881, 17.923, 18.673, 19.4865, 20.3033, 21.1536, 22.0086, 22.9046, 23.8196, 24.7924, 25.8082, 26.8714, 27.9481, 29.0828, 30.2757, 31.5043, 32.8044, 34.1572, 35.6142, 37.1211, 38.6798, 40.3116, 42.0398, 43.8572, 45.6977, 47.6312, 49.6899, 51.815, 54.028, 56.3037, 58.7091, 61.2024, 63.8353, 66.5926, 69.3617, 72.2068, 75.2459, 78.3873, 81.5916, 84.9419, 88.498, 92.1789, 95.9582, 99.8431, 103.739, 107.78, 111.97, 116.312, 120.806, 125.46, 130.269, 135.247, 140.389, 145.713, 151.212, 156.871, 162.729, 168.762, 174.998, 181.424, 188.063, 194.907, 201.942, 209.19, 216.683, 224.37, 232.291, 240.43, 248.807, 257.416, 266.256, 275.348, 284.668, 294.216, 304.053, 314.142, 324.488, 335.101, 345.974, 357.116, 368.547, 380.283, 392.29, 404.564, 417.122, 429.968, 443.116, 456.577, 470.357, 484.422, 498.78, 513.473, 528.479, 543.813, 559.445, 575.411, 591.724, 608.352, 625.344, 642.686, 660.361, 678.371, 696.749, 715.485, 734.608, 754.068, 773.846, 794.046, 814.649, 835.608, 856.972, 878.719, 900.887, 923.409, 946.374, 969.674, 993.435, 1017.62, 1042.21, 1067.28, 1092.72, 1118.64, 1144.96, 1171.71, 1198.98, 1226.67, 1254.82, 1283.46, 1312.65, 1342.21, 1372.27, 1402.85, 1433.93, 1465.49, 1497.62, 1530.29, 1563.49, 1597.22, 1631.49, 1666.37, 1701.8, 1737.75, 1774.35, 1811.51, 1849.29, 1887.75, 1926.79, 1966.6, 2006.97, 2047.99, 2089.71, 2132.1, 2175.23, 2219.17, 2263.72, 2309.2, 2355.43, 2402.47, 2450.33, 2499.05, 2548.66, 2599.16, 2650.59, 2703.03, 2756.32, 2810.75, 2866.27, 2922.91, 2980.54, 3039.47, 3099.53, 3160.98, 3223.66, 3287.71, 3353.18, 3420.34, 3489.13, 3559.72, 3632.06, 3706.18, 3782.42, 3860.78, 3941.42, 4024.52, 4110.27, 4199.4, 4292.8, 4394.49, 4519.52, 5199.95};

Int_t getHiBinFromhiHF(const Double_t hiHF)
{
  Int_t binPos = -1;
  for(int i = 0; i < nBins; ++i){
    if(hiHF >= binTable[i] && hiHF < binTable[i+1]){
      binPos = i;
      break;
    }
  }

  binPos = nBins - 1 - binPos;

  return (Int_t)(200*((Double_t)binPos)/((Double_t)nBin));
}


void SkimTree_Event(int nevt=-1, bool isMC = false, int kTrigSel = 13) 
{

  using namespace std;
  using namespace hi;

  // Example of using event plane namespace 
  cout << " Index of "<< EPNames[HFm2] << " = " << HFm2 << endl;
  cout << " Index of "<< EPNames[HFp2] << " = " << HFp2 << endl;
  cout << " Index of "<< EPNames[trackmid2] << " = " << trackmid2 << endl;

  //TString fname1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/MinBias/HIMinimumBias_Run2018_Upsilon_PromptReco_v1.root";
  TString fnameData1 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v1_Oniatree_addvn_part*.root";
  TString fnameData2 = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/DoubleMuonPD/PromptAOD_v2_Oniatree_addvn_part*.root";
  TString fnameDataReReco = "/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuon/ReReco_Oniatree_addvn_part*.root";
  TString fnameMC = "/eos/cms/store/group/phys_heavyions/gbak/Ups1SMM_MC/Oniatree_Ups1SMM_5p02TeV_TuneCP5_Embd_RECO_MC.root";

  TChain *mytree = new TChain("myTree");
  if(!isMC){
    mytree->Add(fnameDataReReco.Data());
  }
  else if(isMC){
    mytree->Add(fnameMC.Data());
  }

  const int maxBranchSize = 1000;

  UInt_t          runNb;
  UInt_t          eventNb, LS;
  float           zVtx;
  Int_t           Centrality;
  ULong64_t       HLTriggers;
  Float_t         SumET_HF;
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
  TBranch        *b_SumET_HF;   //!
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
  mytree->SetBranchAddress("SumET_HF", &SumET_HF, &b_SumET_HF);
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
  TBranch *b_Reco_mu_whichGen;
  Float_t Gen_weight;
  TBranch *b_Gen_weight;
  if(isMC){
    mytree->SetBranchAddress("Reco_mu_whichGen",Reco_mu_whichGen, &b_Reco_mu_whichGen);
    mytree->SetBranchAddress("Gen_weight",&Gen_weight, &b_Gen_weight);
  }


  TChain *eptree = new TChain("tree");
  if(!isMC){
    eptree->Add(fnameData1.Data());
    eptree->Add(fnameData2.Data());
  }
  else if(isMC){
    eptree->Add(fnameMC.Data());
  }
  
  
  const int nEP = 29;  // number of event planes in the tree
  double qx[nEP]; 
  double qy[nEP]; 
  TBranch *b_qx;
  TBranch *b_qy;
  eptree->SetBranchAddress("qx",qx, &b_qx);
  eptree->SetBranchAddress("qy",qy, &b_qy);
  
  TFile* newfile;
  newfile = new TFile(Form("OniaFlowSkim_UpsTrig_DBPD_isMC%d_190710.root",isMC),"recreate");

  const static int nMaxDimu = 1000;
  int evt;
  int runN;
  int lumi;
  int cBin;
  int nDimu;
  float vz;
  float mass[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu];
  float phi[nMaxDimu];
  float eta[nMaxDimu];
  float eta1[nMaxDimu];
  float eta2[nMaxDimu];
  float phi1[nMaxDimu];
  float phi2[nMaxDimu];
  float pt1[nMaxDimu];
  float pt2[nMaxDimu];
  float weight0[nMaxDimu]; 
  float weight1[nMaxDimu]; 
  float qxa[nMaxDimu];
  float qya[nMaxDimu];
  float qxb[nMaxDimu];
  float qyb[nMaxDimu];
  float qxc[nMaxDimu];
  float qyc[nMaxDimu];
  float qxdimu[nMaxDimu];
  float qydimu[nMaxDimu];
  float qxmupl[nMaxDimu];
  float qxmumi[nMaxDimu];
  float qymupl[nMaxDimu];
  float qymumi[nMaxDimu];
  int recoQQsign[nMaxDimu];
  double weight = 1;

  TTree* mmevttree = new TTree("mmepevt","dimuonAndEventPlanes in event based");
  mmevttree->SetMaxTreeSize(MAXTREESIZE);
  mmevttree->Branch("event",&evt,"event/I");
  mmevttree->Branch("runN",&runN,"runN/I");
  mmevttree->Branch("lumi",&lumi,"lumi/I");
  mmevttree->Branch("cBin",&cBin,"cBin/I");
  mmevttree->Branch("vz",&vz,"vz/F");
  mmevttree->Branch("nDimu",&nDimu,"nDimu/I");
  mmevttree->Branch("mass",mass,"mass[nDimu]/F");
  mmevttree->Branch("y",y,"y[nDimu]/F");
  mmevttree->Branch("pt",pt,"pt[nDimu]/F");
  mmevttree->Branch("pt1",pt1,"pt1[nDimu]/F");
  mmevttree->Branch("pt2",pt2,"pt2[nDimu]/F");
  mmevttree->Branch("eta",eta,"eta[nDimu]/F");
  mmevttree->Branch("eta1",eta1,"eta1[nDimu]/F");
  mmevttree->Branch("eta2",eta2,"eta2[nDimu]/F");
  mmevttree->Branch("qxa",qxa,"qxa[nDimu]/F");
  mmevttree->Branch("qxb",qxb,"qxb[nDimu]/F");
  mmevttree->Branch("qxc",qxc,"qxc[nDimu]/F");
  mmevttree->Branch("qya",qya,"qya[nDimu]/F");
  mmevttree->Branch("qyb",qyb,"qyb[nDimu]/F");
  mmevttree->Branch("qyc",qyc,"qyc[nDimu]/F");
  mmevttree->Branch("qxdimu",qxdimu,"qxdimu[nDimu]/F");
  mmevttree->Branch("qydimu",qydimu,"qydimu[nDimu]/F");
  mmevttree->Branch("qxmupl",qxmupl,"qxmupl[nDimu]/F");
  mmevttree->Branch("qxmumi",qxmumi,"qxmumi[nDimu]/F");
  mmevttree->Branch("qymupl",qymupl,"qymupl[nDimu]/F");
  mmevttree->Branch("qymumi",qymumi,"qymumi[nDimu]/F");
  mmevttree->Branch("recoQQsign",recoQQsign,"recoQQsign[nDimu]/I");
  mmevttree->Branch("weight",&weight,"weight/D");
      


  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies 
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;



  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();

  cout << "Total events = " << nevt << ", : " << eptree->GetEntries() << endl;

  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);
    eptree->GetEntry(iev);
  
    nDimu = 0;
    
    if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;

    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      runN = runNb;
      evt = eventNb;
      lumi = LS;
      cBin = getHiBinFromhiHF(SumET_HF);
      vz = zVtx;


      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[irqq]);
      mumi_Reco = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[irqq]);
      
      weight = 1.;
      if(isMC) weight = findNcoll(Centrality) * Gen_weight;

      if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;
     
      if(isMC){
        if(Reco_mu_whichGen[Reco_QQ_mupl_idx[irqq]] == -1) continue;
        if(Reco_mu_whichGen[Reco_QQ_mumi_idx[irqq]] == -1) continue;
      }
      
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

      if ( !(muplSoft && mumiSoft) ) 
        continue;   
      
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) 
        continue;
   
      recoQQsign[irqq] = Reco_QQ_sign[irqq];     


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
      if ( JP_Reco->Eta() < 0 )  {  
        qxa[nDimu] = qx[HFp2];
        qya[nDimu] = qy[HFp2];
        qxb[nDimu] = qx[HFm2];
        qyb[nDimu] = qy[HFm2];

      }
      else {
        qxa[nDimu] = qx[HFm2];
        qya[nDimu] = qy[HFm2];
        qxb[nDimu] = qx[HFp2];
        qyb[nDimu] = qy[HFp2];
      }
      

      mass[nDimu] = JP_Reco->M();
      phi[nDimu] = JP_Reco->Phi();
      phi1[nDimu] = mupl_Reco->Phi();
      phi2[nDimu] = mumi_Reco->Phi();
      eta[nDimu] = JP_Reco->Eta();
      y[nDimu] = JP_Reco->Rapidity();
      pt[nDimu] = JP_Reco->Pt();
      pt1[nDimu] = mupl_Reco->Pt();
      pt2[nDimu] = mumi_Reco->Pt();
      eta1[nDimu] = mupl_Reco->Eta();
      eta2[nDimu] = mumi_Reco->Eta();
      qxdimu[nDimu] = TMath::Cos(2*phi[nDimu]);
      qydimu[nDimu] = TMath::Sin(2*phi[nDimu]);
      qxmupl[nDimu] = TMath::Cos(2*phi1[nDimu]);
      qxmumi[nDimu] = TMath::Cos(2*phi2[nDimu]);
      qymupl[nDimu] = TMath::Sin(2*phi1[nDimu]);
      qymumi[nDimu] = TMath::Sin(2*phi2[nDimu]);
      nDimu++;

    } // end of dimuon loop
    
    if(nDimu>0) mmevttree->Fill();  
    
  } //end of event loop
//  mmtree->Write();  // Don't need to call Write() for trees
  newfile->cd();
  mmevttree->Write();
  newfile->Close();
  
} 
