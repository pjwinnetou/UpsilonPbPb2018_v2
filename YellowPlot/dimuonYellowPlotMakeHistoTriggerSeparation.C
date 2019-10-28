//Description{{{
//
// dimuonYellowPlotMakeHisto.C
// CMS Heavy-Ion dilepton group
// Modifications: MCBS, Nov 2015
// Macro that will:
// 1) Open a file containing trees with CMS dimuons, onia trees.
// 2) Book a histogram of the dimuon invariant mass
// 3) Write to a file
// The drawing with all the bells and whistles is left to a follow-up macro
// Needs dimuonYellowPlotDrawHisto.C for the code to do the drawing.
// The goal is for this macro to be run only when the data need to be made into
// a histogram, as this is time consuming.  The macro to draw can then be
// invoked as many times as needed to make plots look as desired, without the
// overhead of the call to the tree->Draw command that can take many minutes.
//
//}}} 

//Headers{{{
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TString.h>
#include <iostream>
//}}}

void dimuonYellowPlotMakeHistoTriggerSeparation(TString histoFileName="dimuonMassYellowPlotHistosTriggerSeparation_RERECO_fixed.root") {
  cout << "dimuonYellowPlotMakeHisto: Starting macro dimuonYellowPlotMakeHisto" << endl;
  
 //Input data files with TFile(obsolete){{{ 
  // 2013 pPb data:
  //inFile = new TFile("pPbData_2nd_PromptReco-v1_GR_P_V43D_pileupRej_newAccCut_tot.root","READ");
  
  // 2013 pp data:
  // inFile = new TFile("All_v2.24_Histos_Runs_211739-211831_GlbGlb_woPileUpRej_muLessPV.root","READ");
  
  // 2011 PbPb data:
  //inFile = new TFile("ForYellowPlot_Jpsi_Histos_181611-183013_lowmass_147invmub.root","READ");
 
  // 2015 pp data, test
  //inFile = new TFile("/afs/cern.ch/user/c/chflores/work/public/ForManuel/OniaTree_Run262167.root");
  //inFile = new TFile("/afs/cern.ch/user/a/anstahll/work/public/OniaTree_262163_262277.root");

  // 2015 pp data, Prompt Reco
  //inFile = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/pp502TeV/TTrees/PromptReco/OniaTree_DoubleMu_PromptReco_262081_262273.root");

  //
  // Open the appropriate input file
  //
  //TFile *inFile1, *inFile2, *inFile3, *inFile4, *inFile5;
  //TFile *inFile1 = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262548_262988.root");
  //TFile *inFile1 = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262548_263729_noCUT.root");
  
  //TFile *inFile1 = TFile::Open("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0_HIRun2015-PromptReco-v1_Run_262548_263757_noCUT.root");
//}}}

//Input data files wiht TChain(Current using){{{
  //TChain* theChain = new TChain("hionia/myTree");
  TChain* theChain = new TChain("myTree");
  // // 2015 PbPb data, 
  // PromptAOD L1DoubleMu0 (and B, C, D Physics Datasets), Jan 2016
  
  //theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/Express/v1/000/326/483/HiOniatree_326483.root");
  //theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/Express/v1/000/326/476/HiOniatree_326476.root");
  //theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/Express/v1/000/326/477/HiOniatree_326477.root");
  //theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/Express/v1/000/326/478/HiOniatree_326478.root");
  //theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/Express/v1/000/326/479/HiOniatree_326479.root");
  //theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/Express/v1/000/326/480/HiOniatree_326480.root");
  //theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/Express/v1/000/326/482/HiOniatree_326482.root");
  //theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/Express/v1/000/326/392/HiOniatree_326392.root");
  //theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/PromptAOD/Oniatree_326381_326941_PromptAOD_EvtSel_Nov25th.root");
  //theChain->Add("/eos/cms/store/group/phys_heavyions/gbak/2018PbPbData/OniaVN/Re_Reco/HIDoubleMuon/crab_HIDoubleMu_Run2018A_Promptv1_OniaTreeGF_VN_v3_Cert_326381_327564_Reduced_ReReco_DoubleMuon_part1/190610_114234/Oniatree_addvn_part1_0000.root");
  theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuon/ReReco_Oniatree_addvn_part1.root");
  theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuon/ReReco_Oniatree_addvn_part2.root");
  theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuon/ReReco_Oniatree_addvn_part3.root");
  theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuon/ReReco_Oniatree_addvn_part4.root");
  theChain->Add("/eos/cms/store/group/phys_heavyions/dileptons/Data2018/PbPb502TeV/TTrees/ReReco/AOD/DoubleMuon/ReReco_Oniatree_addvn_part5.root");
  //cout << "Added A file, entries: " << theChain->GetEntries() << endl;

/*
  theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0_HIRun2015-PromptReco-v1_Run_262548_263757_noCUT.root");
  cout << "Added A file, entries: " << theChain->GetEntries() << endl;

  theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0B_HIRun2015-PromptReco-v1_Run_263322_263757_noCUT.root");
  cout << "Added B file, entries: " << theChain->GetEntries() << endl;

  theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0C_HIRun2015-PromptReco-v1_Run_263322_263757_noCUT.root");
  cout << "Added C file, entries: " << theChain->GetEntries() << endl;

  theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/PromptAOD/OniaTree_HIOniaL1DoubleMu0D_HIRun2015-PromptReco-v1_Run_263322_263757_noCUT.root");
  cout << "Added D file, entries: " << theChain->GetEntries() << endl;
*/

  
  
  // Express Streams used in Nov-Dec 2015
  // theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262548_262640.root");
  // theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262656_262694.root");
  // theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262695_262697.root");
  // theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262698_262699.root");
  // theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262701_262726.root");
  // theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262731_262735.root");
  // theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262783_262784.root");
  // theChain->Add("root://eoscms//eos/cms/store/group/phys_heavyions/dileptons/Data2015/PbPb502TeV/TTrees/ExpressStream/OniaTree_262811_262811.root");
//}}}
  
  // The histograms
  // sse variable size bins, multiplying the bin limit of a previous bin by a scale factor.
  // This will make the bins span about the same width in a log-x scale.
  
  cout << "dimuonYellowPlot: creating bin limit array..." << endl;

//Define bins with diferent bin width{{{
  double bins[100000];
  bins[0] = 0.2;
  int nBins = 0;
  double massBinScaleFactor = 1.015; // Typical values: //1.04; //1.015;
  // The for loop should give enough bins to span up to mass of about 200 GeV/c^2.
  for (int i=1; bins[i-1]<200; i++) {
    bins[i] = bins[i-1]*massBinScaleFactor; 
    nBins++;
  }
//}}}

  cout << "dimuonYellowPlot: booking mass histogram with ";
  cout << nBins << " mass bins, mass limits: ";
  cout << bins[0] << " - " << bins[nBins] << endl;

//Trigger bit and names{{{
  // PbPb trigger bits for 2015:
  // Bit  1: HLT_HIL1DoubleMu0_v1
  // Bit  4: HLT_HIL1DoubleMu10_v1                    Prescale 1
  // Bit 14: HLT_HIL2DoubleMu0_Cent30_OS_NHitQ_v1     Prescale 50
  // Bit 16: HLT_HIL3DoubleMu0_Cent30_OS_m2p5to4p5_v1 Prescale 8 
  // Bit 17: HLT_HIL3DoubleMu0_Cent30_OS_m7to14_v1    Prescale 7 
  // Bit 39: HLT_HIL2Mu20_2HF_v1                      Prescale 2
  // Note, bit 39 in long is 274877906944
  
  // PbPb trigger bits for 2018:
  // Bit  1: HLT_HIL1DoubleMuOpen_v1 
  // Bit 13: HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 Prescale ?
  // Bit 14: HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v1      Prescale ?
  // Bit 24: HLT_HIL3Mu12_v1                            Prescale ?

  // pp trigger bits:
  // Bit  1: HLT_HIL1DoubleMu0_v1
  // Bit  2: HLT_HIL1DoubleMu10_v1
  // Bit  4: HLT_HIL3DoubleMu0_OS_m2p5to4p5_v1
  // Bit  5: HLT_HIL3DoubleMu0_OS_m7to14_v1
  // Bit  6: HLT_HIL2Mu3_NHitQ10_v1
  // Bit  7: HLT_HIL3Mu3_NHitQ15_v1
  // Bit  8: HLT_HIL2Mu5_NHitQ10_v1
  // Bit  9: HLT_HIL3Mu5_NHitQ15_v1
  // Bit 10: HLT_HIL2Mu7_NHitQ10_v1
  // Bit 11: HLT_HIL3Mu7_NHitQ15_v1
  // Bit 12: HLT_HIL2Mu15_v1
  // Bit 13: HLT_HIL3Mu15_v1
  // Bit 14: HLT_HIL2Mu20_v1
  // Bit 15: HLT_HIL3Mu20_v1
//}}}

	const Int_t NOT = 4;//Number Of Triggers
	const Int_t tBits[NOT] = {0, 12, 13, 23};//bit-1
	TH1F* massHistoUnlikeSignBit[NOT];
	TH1F* massHistoLikeSignBit[NOT];
	TH1F* centralityBit[NOT];
	TH1F* maxDaughterMuonPtBit[NOT];
	for(Int_t itrg = 0; itrg < NOT; itrg++)
	{
		massHistoUnlikeSignBit[itrg] = new TH1F(Form("massHistoUnlikeSignBit%d", tBits[itrg]), Form("massHistoUnlikeSignBit%d;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", tBits[itrg]), nBins, bins);
		massHistoLikeSignBit[itrg] = new TH1F(Form("massHistoLikeSignBit%d", tBits[itrg]), Form("massHistoLikeSignBit%d;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})", tBits[itrg]), nBins, bins);
		centralityBit[itrg] = new TH1F(Form("centralityBit%d", tBits[itrg]), Form("centralityBit%d;Centrality Bin;Counts", tBits[itrg]),200,0,200);
		maxDaughterMuonPtBit[itrg] = new TH1F(Form("maxDaughterMuonPtBit%d", tBits[itrg]), Form("maxDaughterMuonPtBit%d;max(p_{T}^{#mu+},p_{T}^{#mu-}) GeV/c; Counts", tBits[itrg]),100,0,100);
	}

/*
//sample{{{
  TH1F* massHistoUnlikeSignBit1 = new TH1F("massHistoUnlikeSignBit1","massHistoUnlikeSignBit1;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})",nBins,bins);
  TH1F* massHistoUnlikeSignBit4 = new TH1F("massHistoUnlikeSignBit4","massHistoUnlikeSignBit4;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})",nBins,bins);

  TH1F* massHistoLikeSignBit4 = new TH1F("massHistoLikeSignBit4","massHistoLikeSignBit4;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})",nBins,bins);

  TH1F* centralityBit4 = new TH1F("centralityBit4","centralityBit4;Centrality Bin;Counts",200,0,200);
  
  TH1F* maxDaughterMuonPtBit4 = new TH1F("maxDaughterMuonPtBit4","maxDaughterMuonPtBit4;max(p_{T}^{#mu+},p_{T}^{#mu-}) GeV/c; Counts",100,0,100);

	TH1F* massHistoLikeSign = new TH1F("massHistoLikeSign","massHistoLikeSign;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})",nBins,bins);
	TH1F* massHistoSubtracted = new TH1F("massHistoSubtracted","massHistoSubtracted;m_{#mu#mu} (GeV/c^{2});Entries/(GeV/c^{2})",nBins,bins);
	TH2F* h3 = new TH2F("h3","h3;#eta_{#mu^{+}} - #eta_{#mu^{-}};#mu^{+}#mu^{-} mass",500,-5.0,5.0,nBins,bins);
//}}}
*/

  cout << "dimuonYellowPlot: Draw histograms from tree..." << endl;

  int year = 2018; // choose year: 2011, 2013, 2015, 2018
  int triggerOption = 2; // 1= pp triggers, 2 = PbPb triggers
  if (year == 2011) {
    theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSign","(Reco_QQ_trig&1)==1&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
  }
  else if (year == 2013) {
    // New variable added, Reco_QQ_type, which should be 2.
    // this works for pPb and for pp data from 2013
    theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSign","(Reco_QQ_trig&1)==1&&Reco_QQ_type==2&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      
  }
  else if (year == 2015) {
//{{{
    // Note: for pp, type 2 is ok. For PbPb type == 0 (GlobalMuon && Tracker muon)
    //theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSign","(Reco_QQ_trig&1)==1&&Reco_QQ_type==0&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
    //theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSign","(Reco_QQ_trig&1)==1&&Reco_QQ_type==0&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4&&runNb>262703","");
    // Separation by triggers, pp
    switch (triggerOption) {
    case 1:
      
      break;
    case 2:
      // Separation by triggers, PbPb
      cout << "Drawing mass, unlike sign... bit 1" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSignBit1","(HLTriggers&(1<<(1-1))==1<<(1-1))&&(Reco_QQ_trig&(1<<(1-1)))==(1<<(1-1))&&Reco_QQ_type==0&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      cout << "Drawing mass, unlike sign... bit 4" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSignBit4","(HLTriggers&(1<<(4-1))==1<<(4-1))&&(Reco_QQ_trig&(1<<(4-1)))==(1<<(4-1))&&Reco_QQ_type==0&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      cout << "Drawing mass, unlike sign... bit 14" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSignBit14","(HLTriggers&(1<<(14-1))==1<<(14-1))&&(Reco_QQ_trig&(1<<(14-1)))==(1<<(14-1))&&Reco_QQ_type==0&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      cout << "Drawing mass, unlike sign... bit 16" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSignBit16","(HLTriggers&(1<<(16-1))==1<<(16-1))&&(Reco_QQ_trig&(1<<(16-1)))==(1<<(16-1))&&Reco_QQ_type==0&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      cout << "Drawing mass, unlike sign... bit 17" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSignBit17","(HLTriggers&(1<<(17-1))==1<<(17-1))&&(Reco_QQ_trig&(1<<(17-1)))==(1<<(17-1))&&Reco_QQ_type==0&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      cout << "Drawing mass, unlike sign... bit 39" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSignBit39","(HLTriggers&(274877906944)==274877906944)&&(Reco_QQ_trig&(274877906944))==(274877906944)&&Reco_QQ_type==0&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      
      // Draw Like Sign Histograms
      cout << "Drawing mass, like sign... bit 4" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoLikeSignBit4","(HLTriggers&(1<<(4-1))==1<<(4-1))&&(Reco_QQ_trig&(1<<(4-1)))==(1<<(4-1))&&Reco_QQ_type==0&&Reco_QQ_sign!=0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      cout << "Drawing mass, like sign... bit 14" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoLikeSignBit14","(HLTriggers&(1<<(14-1))==1<<(14-1))&&(Reco_QQ_trig&(1<<(14-1)))==(1<<(14-1))&&Reco_QQ_type==0&&Reco_QQ_sign!=0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      cout << "Drawing mass, like sign... bit 16" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoLikeSignBit16","(HLTriggers&(1<<(16-1))==1<<(16-1))&&(Reco_QQ_trig&(1<<(16-1)))==(1<<(16-1))&&Reco_QQ_type==0&&Reco_QQ_sign!=0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      cout << "Drawing mass, like sign... bit 17" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoLikeSignBit17","(HLTriggers&(1<<(17-1))==1<<(17-1))&&(Reco_QQ_trig&(1<<(17-1)))==(1<<(17-1))&&Reco_QQ_type==0&&Reco_QQ_sign!=0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      cout << "Drawing mass, like sign... bit 39" << endl;
      theChain->Draw("Reco_QQ_4mom.M()>>massHistoLikeSignBit39","(HLTriggers&(274877906944)==274877906944)&&(Reco_QQ_trig&(274877906944))==(274877906944)&&Reco_QQ_type==0&&Reco_QQ_sign!=0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
      
      // Draw Centrality Distributions
      cout << "Drawing centrality, bit 4" << endl;
      theChain->Draw("Centrality>>centralityBit4","(HLTriggers&(1<<(4-1))==1<<(4-1))&&(Reco_QQ_trig&(1<<(4-1)))==(1<<(4-1))");
      cout << "Drawing centrality, bit 14" << endl;
      theChain->Draw("Centrality>>centralityBit14","(HLTriggers&(1<<(14-1))==1<<(14-1))&&(Reco_QQ_trig&(1<<(14-1)))==(1<<(14-1))");
      cout << "Drawing centrality, bit 16" << endl;
      theChain->Draw("Centrality>>centralityBit16","(HLTriggers&(1<<(16-1))==1<<(16-1))&&(Reco_QQ_trig&(1<<(16-1)))==(1<<(16-1))");
      cout << "Drawing centrality, bit 17" << endl;
      theChain->Draw("Centrality>>centralityBit17","(HLTriggers&(1<<(17-1))==1<<(17-1))&&(Reco_QQ_trig&(1<<(17-1)))==(1<<(17-1))");
      cout << "Drawing centrality, bit 39" << endl;
      theChain->Draw("Centrality>>centralityBit39","(HLTriggers&274877906944==274877906944)&&(Reco_QQ_trig&(274877906944))==(274877906944)");
      
      // Draw max daughter pt Distributions
      cout << "Drawing max daughter pt, bit 4" << endl;
      theChain->Draw("max(Reco_QQ_mumi_4mom.Pt(),Reco_QQ_mupl_4mom.Pt())>>maxDaughterMuonPtBit4","(HLTriggers&(1<<(4-1))==1<<(4-1))&&(Reco_QQ_trig&(1<<(4-1)))==(1<<(4-1))");
      cout << "Drawing max daughter pt, bit 14" << endl;
      theChain->Draw("max(Reco_QQ_mumi_4mom.Pt(),Reco_QQ_mupl_4mom.Pt())>>maxDaughterMuonPtBit14","(HLTriggers&(1<<(14-1))==1<<(14-1))&&(Reco_QQ_trig&(1<<(14-1)))==(1<<(14-1))");
      cout << "Drawing max daughter pt, bit 16" << endl;
      theChain->Draw("max(Reco_QQ_mumi_4mom.Pt(),Reco_QQ_mupl_4mom.Pt())>>maxDaughterMuonPtBit16","(HLTriggers&(1<<(16-1))==1<<(16-1))&&(Reco_QQ_trig&(1<<(16-1)))==(1<<(16-1))");
      cout << "Drawing max daughter pt, bit 17" << endl;
      theChain->Draw("max(Reco_QQ_mumi_4mom.Pt(),Reco_QQ_mupl_4mom.Pt())>>maxDaughterMuonPtBit17","(HLTriggers&(1<<(17-1))==1<<(17-1))&&(Reco_QQ_trig&(1<<(17-1)))==(1<<(17-1))");
      cout << "Drawing max daughter pt, bit 39" << endl;
      theChain->Draw("max(Reco_QQ_mumi_4mom.Pt(),Reco_QQ_mupl_4mom.Pt())>>maxDaughterMuonPtBit39","(HLTriggers&274877906944==274877906944)&&(Reco_QQ_trig&(274877906944))==(274877906944)");
      break;
    default:
      cout << "Switch option not recognized. Not drawing anything!" << endl;
      break;
    }
//}}}
  }
	else if(year == 2018)
	{
		//cout << "Drawing mass, unlike sign... bit " << tBits[itrg] << endl;
		//theChain->Draw("Reco_QQ_4mom.M()>>massHistoUnlikeSignBit4", "(HLTriggers&(1<<(4-1))==1<<(4-1))&&(Reco_QQ_trig&(1<<(4-1)))==(1<<(4-1))&&Reco_QQ_type==0&&Reco_QQ_sign==0&&Reco_QQ_mupl_4mom.Pt()>4&&Reco_QQ_mumi_4mom.Pt()>4","");
                Int_t nevent = theChain->GetEntries();

                Int_t           nTrig;
                Int_t           trigPrescale[26];   //[nTrig]
                ULong64_t       HLTriggers;
                Int_t           Reco_mu_size;
                Int_t           Reco_QQ_size;
                ULong64_t       Reco_QQ_trig[66];
                TClonesArray    *Reco_mu_4mom;
                TClonesArray    *Reco_QQ_4mom;
                Int_t           Reco_QQ_mupl_idx[1000];
                Int_t           Reco_QQ_mumi_idx[1000];
                Int_t           Centrality;
                Int_t           Reco_QQ_type[66];
                Int_t           Reco_QQ_sign[66];

                TBranch *b_Reco_mu_size;
                TBranch *b_Reco_QQ_size;
                TBranch *b_Reco_mu_4mom;
                TBranch *b_Reco_QQ_4mom;
                TBranch *b_Reco_QQ_mupl_idx;
                TBranch *b_Reco_QQ_mumi_idx;
                TBranch *b_Reco_QQ_sign;
                TBranch *b_Centrality;
                TBranch *b_nTrig;
                TBranch *b_trigPrescale;
                TBranch *b_HLTriggers;
                TBranch *b_Reco_QQ_trig;
                TBranch *b_Reco_QQ_type;

                Reco_mu_4mom = 0;
                Reco_QQ_4mom = 0;

                theChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
                theChain->SetBranchAddress("trigPrescale", trigPrescale, &b_trigPrescale);
                theChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
                theChain->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
                theChain->SetBranchAddress("Reco_QQ_type", Reco_QQ_type, &b_Reco_QQ_type);
                theChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
                theChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
                theChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
                theChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
                theChain->SetBranchAddress("Reco_QQ_mupl_idx", &Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
                theChain->SetBranchAddress("Reco_QQ_mumi_idx", &Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
                theChain->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign, &b_Reco_QQ_sign);
                theChain->SetBranchAddress("Centrality", &Centrality, &b_Centrality);


                for(int i=0; i<nevent; i++){
                  theChain->GetEvent(i);
                  if(i%561143==0){cout<<">>>>> EVENT "<<i<<" / "<<theChain->GetEntries()<<" ("<<(int)(100.*i/theChain->GetEntries())<<"%)"<<endl;}
                  for(int j=0; j<Reco_QQ_size; j++){
                    TLorentzVector *RecoQQ4mom = (TLorentzVector*) Reco_QQ_4mom->At(j);
                    TLorentzVector *RecoQQmupl = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mupl_idx[j]);
                    TLorentzVector *RecoQQmumi = (TLorentzVector*) Reco_mu_4mom->At(Reco_QQ_mumi_idx[j]);
                    for(Int_t itrg = 0; itrg < NOT; itrg++)
                    {
                      //include Reco_QQ_trig{{{
                      if(
                          ((HLTriggers&     ((ULong64_t)pow(2,tBits[itrg])))==((ULong64_t)pow(2,tBits[itrg])))&&
                          ((Reco_QQ_trig[j]&((ULong64_t)pow(2,tBits[itrg])))==((ULong64_t)pow(2,tBits[itrg])))&&
                          Reco_QQ_type[j]==1&&
                          Reco_QQ_sign[j]==0&&
                          RecoQQmupl->Pt()>4&&RecoQQmumi->Pt()>4)
                      {
                      //cout << "Drawing mass, unlike sign... bit " << tBits[itrg] << endl;
                      massHistoUnlikeSignBit[itrg]->Fill(RecoQQ4mom->M());
                      }
                      //if(Form("((HLTriggers&pow(2,%d))==pow(2,%d))&&((Reco_QQ_trig&pow(2,%d))==pow(2,%d))&&Reco_QQ_type==1&&Reco_QQ_sign!=0&&RecoQQmupl.Pt()>4&&RecoQQmumi.Pt()>4", tBits[itrg], tBits[itrg], tBits[itrg], tBits[itrg])){
                      ////cout << "Drawing mass, like sign... bit " << tBits[itrg] << endl;
                      //massHistoLikeSignBit[itrg]->Fill(RecoQQ4mom->M());
                      //}
                      //if(Form("((HLTriggers&pow(2,%d))==pow(2,%d))&&((Reco_QQ_trig&pow(2,%d))==pow(2,%d))", tBits[itrg], tBits[itrg], tBits[itrg], tBits[itrg])){
                      ////cout << "Drawing centrality, bit " << tBits[itrg] << endl;
                      //centralityBit[itrg]->Fill(Centrality);
                      //}
                      //if(Form("((HLTriggers&pow(2,%d))==pow(2,%d))&&((Reco_QQ_trig&pow(2,%d))==pow(2,%d))", tBits[itrg], tBits[itrg], tBits[itrg], tBits[itrg])){
                      ////cout << "Drawing max daughter pt, bit " << tBits[itrg] << endl;
                      //maxDaughterMuonPtBit[itrg]->Fill(max(RecoQQmumi->Pt(),RecoQQmupl->Pt()));
                      //}
                      //theChain->Draw(Form("Reco_QQ_4mom.M()>>massHistoUnlikeSignBit%d", tBits[itrg]), Form("((HLTriggers&pow(2,%d))==pow(2,%d))&&((Reco_QQ_trig&pow(2,%d))==pow(2,%d))&&Reco_QQ_type==1&&Reco_QQ_sign==0&&RecoQQmupl.Pt()>4&&RecoQQmumi.Pt()>4", tBits[itrg], tBits[itrg], tBits[itrg], tBits[itrg]),"");
//                      cout << "Drawing mass, like sign... bit " << tBits[itrg] << endl;
//                      theChain->Draw(Form("Reco_QQ_4mom.M()>>massHistoLikeSignBit%d", tBits[itrg]), Form("((HLTriggers&pow(2,%d))==pow(2,%d))&&((Reco_QQ_trig&pow(2,%d))==pow(2,%d))&&Reco_QQ_type==1&&Reco_QQ_sign!=0&&RecoQQmupl.Pt()>4&&RecoQQmumi.Pt()>4", tBits[itrg], tBits[itrg], tBits[itrg], tBits[itrg]),"");
//                      cout << "Drawing centrality, bit " << tBits[itrg] << endl;
//                      theChain->Draw(Form("Centrality>>centralityBit%d", tBits[itrg]), Form("((HLTriggers&pow(2,%d))==pow(2,%d))&&((Reco_QQ_trig&pow(2,%d))==pow(2,%d))", tBits[itrg], tBits[itrg], tBits[itrg], tBits[itrg]), "");
//                      cout << "Drawing max daughter pt, bit " << tBits[itrg] << endl;
//                      theChain->Draw(Form("max(RecoQQmumi.Pt(),RecoQQmupl.Pt())>>maxDaughterMuonPtBit%d", tBits[itrg]), Form("((HLTriggers&pow(2,%d))==pow(2,%d))&&((Reco_QQ_trig&pow(2,%d))==pow(2,%d))", tBits[itrg], tBits[itrg], tBits[itrg], tBits[itrg]), "");
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
                      //}}}
                      /*
                      //without Reco_QQ_trig{{{
                      cout << "Drawing mass, unlike sign... bit " << tBits[itrg] << endl;
                      theChain->Draw(Form("Reco_QQ_4mom.M()>>massHistoUnlikeSignBit%d", tBits[itrg]), Form("((HLTriggers&pow(2,%d))==pow(2,%d))&&Reco_QQ_type==1&&Reco_QQ_sign==0&&RecoQQmupl->Pt()>4&&RecoQQmumi->Pt()>4", tBits[itrg], tBits[itrg]), "");
                      cout << "Drawing mass, like sign... bit " << tBits[itrg] << endl;
                      theChain->Draw(Form("Reco_QQ_4mom.M()>>massHistoLikeSignBit%d", tBits[itrg]), Form("((HLTriggers&pow(2,%d))==pow(2,%d))&&Reco_QQ_type==1&&Reco_QQ_sign!=0&&RecoQQmupl->Pt()>4&&RecoQQmumi->Pt()>4", tBits[itrg], tBits[itrg]), "");
                      cout << "Drawing centrality, bit " << tBits[itrg] << endl;
                      theChain->Draw(Form("Centrality>>centralityBit%d", tBits[itrg]), Form("((HLTriggers&pow(2,%d))==pow(2,%d))", tBits[itrg], tBits[itrg]), "");
                      cout << "Drawing max daughter pt, bit " << tBits[itrg] << endl;
                      theChain->Draw(Form("max(RecoQQmumi->Pt(),RecoQQmupl->Pt())>>maxDaughterMuonPtBit%d", tBits[itrg]), Form("((HLTriggers&pow(2,%d))==pow(2,%d))", tBits[itrg], tBits[itrg]), "");
                      //}}}
                      */
                    }
                  }
                }
        }
  else {
    cout << "Unrecognized year. No drawing, exit macro." << endl;
    return;
  }

	// Output file creation and writing of histogram
	cout << "dimuonYellowPlotMakeHisto: Open output file... " << histoFileName << endl;
	TFile* outFile = new TFile(histoFileName,"RECREATE");

	cout << "Writing histograms .. " << endl;
	for(Int_t itrg = 0; itrg < NOT; itrg++)
	{
		cout << massHistoUnlikeSignBit[itrg]->GetName() << ", " << massHistoUnlikeSignBit[itrg]->GetEntries() << " dimuons " << endl;
		massHistoUnlikeSignBit[itrg]->Write();
		massHistoLikeSignBit[itrg]->Write();
		centralityBit[itrg]->Write();
		maxDaughterMuonPtBit[itrg]->Write();
	}
  
	outFile->Close();
	return;
}
