//The MC is here: /eos/cms/store/group/phys_heavyions/gbak/Ups1SMM_MC/Oniatree_Ups1SMM_5p02TeV_TuneCP5_Embd_RECO_MC.root 3.9G

//Efficiency = reco/gen

/*Steps:

  1. Fit MC Gen and Reco in all the analysis bins.
  2. Take the ratio, eff=Reco/Gen.

Centrality bins: 
10-30%
30-50%
50-100%

*/



void getEfficiency(
  float ptLow = 0.0, float ptHigh = 50.0,
  float yLow = 0.0, float yHigh = 2.4,
  int cLow = 0, int cHigh = 200
  ) {

  gStyle->SetOptStat(0);

  float muPtCut = 3.5;
  float muEtaCut = 2.4;

  float massLow = 8.0;
  float massHigh = 10.0;

  double min = 0;
  double max = ptHigh;
  double binwidth = 1;
  const int numBins = (max-min)/binwidth;

  //input files
  //TFile *inputMC   = TFile::Open("../../Oniatree_Ups1SMM_5p02TeV_TuneCP5_Embd_RECO_MC_190610.root","READ");
  TFile *inputMC   = TFile::Open("/eos/cms/store/group/phys_heavyions/gbak/2018PbPbMC/OniatreeMC_TuneCP5_HydjetDrumMB_5p02TeV_Pythia8_July_8th.root","READ");
  TTree* MC_QQ_Tree = (TTree*)inputMC->Get("myTree");

  TH1D* hpt_reco = new TH1D("hpt_reco","hpt_reco",numBins,min,max);
  TH1D* hptLowestC_reco = new TH1D("hptLowestC_reco","hptLowestC_reco",numBins,min,max);
  TH1D* hptLowC_reco = new TH1D("hptLowC_reco","hptLowC_reco",numBins,min,max);
  TH1D* hptMidC_reco = new TH1D("hptMidC_reco","hptMidC_reco",numBins,min,max);
  TH1D* hptHighC_reco = new TH1D("hptHighC_reco","hptHighC_reco",numBins,min,max);

  TH1D* hpt_gen = new TH1D("hpt_gen","hpt_gen",numBins,min,max);
  TH1D* hptLowestC_gen = new TH1D("hptLowestC_gen","hptLowestC_gen",numBins,min,max);
  TH1D* hptLowC_gen = new TH1D("hptLowC_gen","hptLowC_gen",numBins,min,max);
  TH1D* hptMidC_gen = new TH1D("hptMidC_gen","hptMidC_gen",numBins,min,max);
  TH1D* hptHighC_gen = new TH1D("hptHighC_gen","hptHighC_gen",numBins,min,max);

  hpt_reco->Sumw2();
  hptLowestC_reco->Sumw2();
  hptLowC_reco->Sumw2();
  hptMidC_reco->Sumw2();
  hptHighC_reco->Sumw2();

  hpt_gen->Sumw2();
  hptLowestC_gen->Sumw2();
  hptLowC_gen->Sumw2();
  hptMidC_gen->Sumw2();
  hptHighC_gen->Sumw2();

  hpt_reco->SetTitle("Reco: Full Centrality");
  hptLowestC_reco->SetTitle("Reco: Centrality 0-10%");
  hptLowC_reco->SetTitle("Reco: Centrality 10-30%");
  hptMidC_reco->SetTitle("Reco: Centrality 30-50%");
  hptHighC_reco->SetTitle("Reco: Centrality 50-100%");

  hpt_reco->GetXaxis()->SetTitle("pt");
  hptLowestC_reco->GetXaxis()->SetTitle("pt");
  hptLowC_reco->GetXaxis()->SetTitle("pt");
  hptMidC_reco->GetXaxis()->SetTitle("pt");
  hptHighC_reco->GetXaxis()->SetTitle("pt");

  hpt_gen->SetTitle("Gen: Full Centrality");
  hptLowestC_gen->SetTitle("Gen: Centrality 0-10%");
  hptLowC_gen->SetTitle("Gen: Centrality 10-30%");
  hptMidC_gen->SetTitle("Gen: Centrality 30-50%");
  hptHighC_gen->SetTitle("Gen: Centrality 50-100%");

  hpt_gen->GetXaxis()->SetTitle("pt");
  hptLowestC_gen->GetXaxis()->SetTitle("pt");
  hptLowC_gen->GetXaxis()->SetTitle("pt");
  hptMidC_gen->GetXaxis()->SetTitle("pt");
  hptHighC_gen->GetXaxis()->SetTitle("pt");

  //Reco Cuts
  TString HLTCut = "(HLTriggers&(pow(2, 13))) == (pow(2, 13))";//Index 13

  TString muIdCuts = Form("Reco_mu_highPurity==1 && Reco_mu_dxy<0.3 && Reco_mu_dz<20 && Reco_mu_nPixWMea>0 && Reco_mu_nTrkWMea>5 && Reco_QQ_VtxProb>.01");

  TString singleMuCuts = Form("Reco_mu_4mom->Pt()>%f && abs(Reco_mu_4mom->Eta())<%f",muPtCut,muEtaCut);

  TString diMuCuts = Form("Reco_QQ_4mom->Pt()>%f && Reco_QQ_4mom->Pt()<%f && abs(Reco_QQ_4mom->Rapidity())>%f && abs(Reco_QQ_4mom->Rapidity())<%f && Centrality>=%i && Centrality<%i && Reco_QQ_4mom->M()>%f && Reco_QQ_4mom->M()<%f",ptLow,ptHigh,yLow,yHigh,cLow,cHigh,massLow,massHigh);

  TString cuts = HLTCut + " && " + muIdCuts + " && " + singleMuCuts + " && " + diMuCuts;
  TString cutsLowestC = cuts + " && Centrality>=0 && Centrality<20";
  TString cutsLowC = cuts + " && Centrality>=20 && Centrality<60";
  TString cutsMidC = cuts + " && Centrality>=60 && Centrality<100";
  TString cutsHighC = cuts + " && Centrality>=100 && Centrality<200";

  //Gen Cuts
  TString singleMuCuts_gen = Form("Gen_mu_4mom->Pt()>%f && abs(Gen_mu_4mom->Eta())<%f",muPtCut,muEtaCut);

  TString diMuCuts_gen = Form("Gen_QQ_4mom->Pt()>%f && Gen_QQ_4mom->Pt()<%f && abs(Gen_QQ_4mom->Rapidity())>%f && abs(Gen_QQ_4mom->Rapidity())<%f && Centrality>=%i && Centrality<%i && Gen_QQ_4mom->M()>%f && Gen_QQ_4mom->M()<%f",ptLow,ptHigh,yLow,yHigh,cLow,cHigh,massLow,massHigh);

  TString cuts_gen = singleMuCuts_gen + " && " + diMuCuts_gen;
  TString cutsLowestC_gen = cuts_gen + " && Centrality>=0 && Centrality<20";
  TString cutsLowC_gen = cuts_gen + " && Centrality>=20 && Centrality<60";
  TString cutsMidC_gen = cuts_gen + " && Centrality>=60 && Centrality<100";
  TString cutsHighC_gen = cuts_gen + " && Centrality>=100 && Centrality<200";

  //Draw
  TCanvas * cpt_reco = new TCanvas("cpt_reco","cpt_reco",0,0,400,400);
  cpt_reco->cd();
  //MC_QQ_Tree->Draw("Reco_QQ_4mom->Pt()>>hpt_reco","Reco_QQ_4mom->Pt()<30");
  MC_QQ_Tree->Draw("Reco_QQ_4mom->Pt()>>hpt_reco",cuts.Data());

  TCanvas * cptLowestC_reco = new TCanvas("cptLowestC_reco","cptLowestC_reco",0,0,400,400);
  cptLowestC_reco->cd();
  MC_QQ_Tree->Draw("Reco_QQ_4mom->Pt()>>hptLowestC_reco",cutsLowestC.Data());

  TCanvas * cptLowC_reco = new TCanvas("cptLowC_reco","cptLowC_reco",0,0,400,400);
  cptLowC_reco->cd();
  MC_QQ_Tree->Draw("Reco_QQ_4mom->Pt()>>hptLowC_reco",cutsLowC.Data());

  TCanvas * cptMidC_reco = new TCanvas("cptMidC_reco","cptMidC_reco",400,0,400,400);
  cptMidC_reco->cd(); 
  MC_QQ_Tree->Draw("Reco_QQ_4mom->Pt()>>hptMidC_reco",cutsMidC.Data());

  TCanvas * cptHighC_reco = new TCanvas("cptHighC_reco","cptHighC_reco",800,0,400,400);
  cptHighC_reco->cd();
  MC_QQ_Tree->Draw("Reco_QQ_4mom->Pt()>>hptHighC_reco",cutsHighC.Data());

  //Gen
  TCanvas * cpt_gen = new TCanvas("cpt_gen","cpt_gen",0,400,400,400);
  cpt_gen->cd();
  //MC_QQ_Tree->Draw("Gen_QQ_4mom->Pt()>>hpt_gen","Gen_QQ_4mom->Pt()<30");
  MC_QQ_Tree->Draw("Gen_QQ_4mom->Pt()>>hpt_gen",cuts_gen.Data());

  TCanvas * cptLowestC_gen = new TCanvas("cptLowestC_gen","cptLowestC_gen",0,400,400,400);
  cptLowestC_gen->cd();
  MC_QQ_Tree->Draw("Gen_QQ_4mom->Pt()>>hptLowestC_gen",cutsLowestC_gen.Data());

  TCanvas * cptLowC_gen = new TCanvas("cptLowC_gen","cptLowC_gen",0,400,400,400);
  cptLowC_gen->cd();
  MC_QQ_Tree->Draw("Gen_QQ_4mom->Pt()>>hptLowC_gen",cutsLowC_gen.Data());

  TCanvas * cptMidC_gen = new TCanvas("cptMidC_gen","cptMidC_gen",0,400,400,400);
  cptMidC_gen->cd();
  MC_QQ_Tree->Draw("Gen_QQ_4mom->Pt()>>hptMidC_gen",cutsMidC_gen.Data());

  TCanvas * cptHighC_gen = new TCanvas("cptHighC_gen","cptHighC_gen",0,400,400,400);
  cptHighC_gen->cd();
  MC_QQ_Tree->Draw("Gen_QQ_4mom->Pt()>>hptHighC_gen",cutsHighC_gen.Data());

  cpt_reco->Update();
  cptLowestC_reco->Update();
  cptLowC_reco->Update();
  cptMidC_reco->Update();
  cptHighC_reco->Update();

  //Divide
  TH1D* hpt_eff;
  TH1D* hptLowestC_eff;
  TH1D* hptLowC_eff;
  TH1D* hptMidC_eff;
  TH1D* hptHighC_eff;

  hpt_eff = (TH1D*)hpt_reco->Clone("hpt_eff");
  hptLowestC_eff = (TH1D*)hptLowestC_reco->Clone("hptLowestC_eff");
  hptLowC_eff = (TH1D*)hptLowC_reco->Clone("hptLowC_eff");
  hptMidC_eff = (TH1D*)hptMidC_reco->Clone("hptMidC_eff");
  hptHighC_eff = (TH1D*)hptHighC_reco->Clone("hptHighC_eff");

  hpt_eff->Divide(hpt_gen);
  hptLowestC_eff->Divide(hptLowestC_gen);
  hptLowC_eff->Divide(hptLowC_gen);
  hptMidC_eff->Divide(hptMidC_gen);
  hptHighC_eff->Divide(hptHighC_gen);

  hpt_eff->SetTitle("Eff: Full Centrality");
  hptLowestC_eff->SetTitle("Eff: Centrality 0-10%");
  hptLowC_eff->SetTitle("Eff: Centrality 10-30%");
  hptMidC_eff->SetTitle("Eff: Centrality 30-50%");
  hptHighC_eff->SetTitle("Eff: Centrality 50-100%");

  TCanvas * cpt_eff = new TCanvas("cpt_eff","cpt_eff",0,400,400,400);
  cpt_eff->cd();
  hpt_eff->Draw();

  TCanvas * cptLowestC_eff = new TCanvas("cptLowestC_eff","cptLowestC_eff",0,400,400,400);
  cptLowestC_eff->cd();
  hptLowestC_eff->Draw();

  TCanvas * cptLowC_eff = new TCanvas("cptLowC_eff","cptLowC_eff",0,400,400,400);
  cptLowC_eff->cd();
  hptLowC_eff->Draw();

  TCanvas * cptMidC_eff = new TCanvas("cptMidC_eff","cptMidC_eff",400,400,400,400);
  cptMidC_eff->cd();
  hptMidC_eff->Draw();

  TCanvas * cptHighC_eff = new TCanvas("cptHighC_eff","cptHighC_eff",800,400,400,400);
  cptHighC_eff->cd();
  hptHighC_eff->Draw();

  //Save canvases
  cpt_reco->SaveAs("Plots/cpt_reco.png");
  cptLowestC_reco->SaveAs("Plots/cptLowestC_reco.png");
  cptLowC_reco->SaveAs("Plots/cptLowC_reco.png");
  cptMidC_reco->SaveAs("Plots/cptMidC_reco.png");
  cptHighC_reco->SaveAs("Plots/cptHighC_reco.png");
  cpt_reco->SaveAs("Plots/cpt_reco.pdf");
  cptLowestC_reco->SaveAs("Plots/cptLowestC_reco.pdf");
  cptLowC_reco->SaveAs("Plots/cptLowC_reco.pdf");
  cptMidC_reco->SaveAs("Plots/cptMidC_reco.pdf");
  cptHighC_reco->SaveAs("Plots/cptHighC_reco.pdf");

  cpt_gen->SaveAs("Plots/cpt_gen.png");
  cptLowestC_gen->SaveAs("Plots/cptLowestC_gen.png");
  cptLowC_gen->SaveAs("Plots/cptLowC_gen.png");
  cptMidC_gen->SaveAs("Plots/cptMidC_gen.png");
  cptHighC_gen->SaveAs("Plots/cptHighC_gen.png");
  cpt_gen->SaveAs("Plots/cpt_gen.pdf");
  cptLowestC_gen->SaveAs("Plots/cptLowestC_gen.pdf");
  cptLowC_gen->SaveAs("Plots/cptLowC_gen.pdf");
  cptMidC_gen->SaveAs("Plots/cptMidC_gen.pdf");
  cptHighC_gen->SaveAs("Plots/cptHighC_gen.pdf");

  cpt_eff->SaveAs("Plots/cpt_eff.png");
  cptLowestC_eff->SaveAs("Plots/cptLowestC_eff.png");
  cptLowC_eff->SaveAs("Plots/cptLowC_eff.png");
  cptMidC_eff->SaveAs("Plots/cptMidC_eff.png");
  cptHighC_eff->SaveAs("Plots/cptHighC_eff.png");
  cpt_eff->SaveAs("Plots/cpt_eff.pdf");
  cptLowestC_eff->SaveAs("Plots/cptLowestC_eff.pdf");
  cptLowC_eff->SaveAs("Plots/cptLowC_eff.pdf");
  cptMidC_eff->SaveAs("Plots/cptMidC_eff.pdf");
  cptHighC_eff->SaveAs("Plots/cptHighC_eff.pdf");

  //Save efficiency files for later use.
  hpt_eff->SetName("mc_eff_vs_pt_noTnP_Cent0100");
  hptLowestC_eff->SetName("mc_eff_vs_pt_noTnP_Cent010");
  hptLowC_eff->SetName("mc_eff_vs_pt_noTnP_Cent1030");
  hptMidC_eff->SetName("mc_eff_vs_pt_noTnP_Cent3050");
  hptHighC_eff->SetName("mc_eff_vs_pt_noTnP_Cent50100");
  TString outFileName = "mc_eff_vs_pt_noTnP_20190614.root";
  TFile* outFile = new TFile(outFileName,"RECREATE");
  hpt_eff->Write();
  hptLowestC_eff->Write();
  hptLowC_eff->Write();
  hptMidC_eff->Write();
  hptHighC_eff->Write();
  outFile->Close();

}
