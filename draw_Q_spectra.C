#include <iostream>
#include "commonUtility.h"
#include "cutsAndBinUpsilonV2.h"
#include "HiEvtPlaneList.h"
#include "Style_jaebeom.h"

using namespace std;
using namespace hi;

void GetHistSqrt(TH1D* h1 =0, TH1D* h2=0);

void draw_Q_spectra(int cLow = 20, int cHigh = 80,
                float ptLow = 0, float ptHigh = 30, 
                float yLow = 0, float yHigh=2.4,
                float SiMuPtCut = 4)
{
  gStyle->SetOptStat(0);
  TFile *rf = new TFile("TEST.root","read");
  TTree *tree = (TTree*) rf -> Get("mmepevt");

  TH1::SetDefaultSumw2();
  
  const int nMaxDimu = 1000;
  float mass[nMaxDimu];
  float qxa[nMaxDimu];
  float qya[nMaxDimu];
  float qxb[nMaxDimu]; 
  float qyb[nMaxDimu];
  float qxc[nMaxDimu];
  float qyc[nMaxDimu];
  float qxdimu[nMaxDimu];
  float qydimu[nMaxDimu];
  float pt[nMaxDimu];
  float y[nMaxDimu]; 
  float pt1[nMaxDimu];
  float pt2[nMaxDimu];
  float eta[nMaxDimu];
  float eta1[nMaxDimu];
  float eta2[nMaxDimu];
  Int_t cBin;
  Int_t event; 
  Int_t nDimu; 

  TBranch *b_event;
  TBranch *b_cBin;
  TBranch *b_nDimu;
  TBranch *b_mass;
  TBranch *b_pt;
  TBranch *b_y;
  TBranch *b_eta;
  TBranch *b_eta1;
  TBranch *b_eta2;
  TBranch *b_pt1;
  TBranch *b_pt2;
  TBranch *b_qxa;
  TBranch *b_qxb;
  TBranch *b_qxc;
  TBranch *b_qxdimu;
  TBranch *b_qya;
  TBranch *b_qyb;
  TBranch *b_qyc;
  TBranch *b_qydimu;
  

  tree -> SetBranchAddress("event", &event, &b_event);
  tree -> SetBranchAddress("cBin", &cBin, &b_cBin);
  tree -> SetBranchAddress("nDimu", &nDimu, &b_nDimu);
  tree -> SetBranchAddress("mass", mass, &b_mass);
  tree -> SetBranchAddress("y", y, &b_y);
  tree -> SetBranchAddress("pt", pt, &b_pt);
  tree -> SetBranchAddress("pt1", pt1, &b_pt1);
  tree -> SetBranchAddress("pt2", pt2, &b_pt2);
  tree -> SetBranchAddress("eta", eta, &b_eta);
  tree -> SetBranchAddress("eta1", eta1, &b_eta1);
  tree -> SetBranchAddress("eta2", eta2, &b_eta2);
  tree -> SetBranchAddress("qxa", qxa, &b_qxa);
  tree -> SetBranchAddress("qxb", qxb, &b_qxb);
  tree -> SetBranchAddress("qxc", qxc, &b_qxc);
  tree -> SetBranchAddress("qxdimu", qxdimu, &b_qxdimu);
  tree -> SetBranchAddress("qya", qya, &b_qya);
  tree -> SetBranchAddress("qyb", qyb, &b_qyb);
  tree -> SetBranchAddress("qyc", qyc, &b_qyc);
  tree -> SetBranchAddress("qydimu", qydimu, &b_qydimu);
  

  const int nMassBin = 10;
  float massLow = 8;
  float massHigh = 14;

  float massBinDiff = (massHigh-massLow)/nMassBin;
  float massBin[nMassBin+1];
  for(int i=0; i<=nMassBin; i++){
    massBin[i] = massLow + massBinDiff*i;
  }

  int nQ2DBin = 200;
  double q2Dlow = -250;
  double q2Dhigh = 250;

  TH2D* h_qa_2D = new TH2D("h_qa_2D",";Q_{xA};Q_{yA}",nQ2DBin,q2Dlow,q2Dhigh,nQ2DBin,q2Dlow,q2Dhigh);
  TH2D* h_qb_2D = new TH2D("h_qb_2D",";Q_{xB};Q_{yB}",nQ2DBin,q2Dlow,q2Dhigh,nQ2DBin,q2Dlow,q2Dhigh);
  TH2D* h_qc_2D = new TH2D("h_qc_2D",";Q_{xC};Q_{yC}",nQ2DBin,q2Dlow,q2Dhigh,nQ2DBin,q2Dlow,q2Dhigh);
  TH2D* h_q_2D = new TH2D("h_q_2D",";Q_{x};Q_{y}",nQ2DBin,-1.1,1.1,nQ2DBin,-1.1,1.1);
  

  double mass_low_SB1 = 8; 
  double mass_high_SB1 = 9; 
  double mass_low_SB = 9; 
  double mass_high_SB = 10; 
  double mass_low_SB2 = 10; 
  double mass_high_SB2 = 14; 

  int nQBin = 100;
  const int nMass_div = 3;

  TString fSB[nMass_div] = {"SB1 (8<m<9)","SB2 (10m<14)","S (9<m<10)"};

  TH1D* h_qxa[nMass_div]; 
  TH1D* h_qya[nMass_div]; 
  TH1D* h_qxb[nMass_div]; 
  TH1D* h_qyb[nMass_div]; 
  TH1D* h_qxc[nMass_div]; 
  TH1D* h_qyc[nMass_div]; 
  TH1D* h_qxdimu[nMass_div]; 
  TH1D* h_qydimu[nMass_div];
  for(int i=0; i<nMass_div; i++){
   h_qxa[i] = new TH1D(Form("h_qxa_%d",i),";Q_{xA};Counts",nQBin,q2Dlow,q2Dhigh);
   h_qya[i] = new TH1D(Form("h_qya_%d",i),";Q_{yA};Counts",nQBin,q2Dlow,q2Dhigh);
   h_qxb[i] = new TH1D(Form("h_qxb_%d",i),";Q_{xB};Counts",nQBin,q2Dlow,q2Dhigh);
   h_qyb[i] = new TH1D(Form("h_qyb_%d",i),";Q_{yB};Counts",nQBin,q2Dlow,q2Dhigh);
   h_qxc[i] = new TH1D(Form("h_qxc_%d",i),";Q_{xC};Counts",nQBin,q2Dlow,q2Dhigh);
   h_qyc[i] = new TH1D(Form("h_qyc_%d",i),";Q_{yC};Counts",nQBin,q2Dlow,q2Dhigh);
   h_qxdimu[i] = new TH1D(Form("h_qxdimu_%d",i),";Q_{x};Counts",nQBin,-1.1,1.1);
   h_qydimu[i] = new TH1D(Form("h_qydimu_%d",i),";Q_{y};Counts",nQBin,-1.1,1.1);
  }
  

  int nDimuPass=0;
  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;
  for(int i=0; i<nEvt; i++){
    tree->GetEntry(i);
    nDimuPass=0;
    
    if(!(cBin>cLow&&cBin<cHigh)) continue; 

    for(int j=0; j<nDimu; j++){
      if(mass[j]<massLow || mass[j]>massHigh) continue;
      if(pt[j]>ptLow&&pt[j]<ptHigh&&abs(y[j])<yHigh&&abs(y[j])>yLow&&pt1[j]>SiMuPtCut&&pt2[j]>SiMuPtCut&&abs(eta1[j])<2.4&&abs(eta2[j])<2.4){
        nDimuPass++;
      }
    }

    if(nDimuPass>1) continue;
    
    for(int j=0; j<nDimu; j++){
      if(mass[j]<massLow || mass[j]>massHigh) continue;
      if(pt[j]>ptLow&&pt[j]<ptHigh&&abs(y[j])<yHigh&&abs(y[j])>yLow&&pt1[j]>SiMuPtCut&&pt2[j]>SiMuPtCut&&abs(eta1[j])<2.4&&abs(eta2[j])<2.4){
        h_qa_2D -> Fill(qxa[j],qya[j]);
        h_qb_2D -> Fill(qxb[j],qyb[j]);
        h_qc_2D -> Fill(qxc[j],qyc[j]);
        h_q_2D -> Fill(qxdimu[j],qydimu[j]);
        if(mass[j]>=mass_low_SB1 && mass[j]<mass_high_SB1){
          h_qxa[0]->Fill(qxa[j]);
          h_qya[0]->Fill(qya[j]);
          h_qxb[0]->Fill(qxb[j]);
          h_qyb[0]->Fill(qyb[j]);
          h_qxc[0]->Fill(qxc[j]);
          h_qyc[0]->Fill(qyc[j]);
          h_qxdimu[0]->Fill(qxdimu[j]);
          h_qydimu[0]->Fill(qydimu[j]);
        }
        else if(mass[j]>=mass_low_SB2 && mass[j]<mass_high_SB2){
          h_qxa[1]->Fill(qxa[j]);
          h_qya[1]->Fill(qya[j]);
          h_qxb[1]->Fill(qxb[j]);
          h_qyb[1]->Fill(qyb[j]);
          h_qxc[1]->Fill(qxc[j]);
          h_qyc[1]->Fill(qyc[j]);
          h_qxdimu[1]->Fill(qxdimu[j]);
          h_qydimu[1]->Fill(qydimu[j]);
        }
        else if(mass[j]>=mass_low_SB && mass[j]<mass_high_SB){
          h_qxa[2]->Fill(qxa[j]);
          h_qya[2]->Fill(qya[j]);
          h_qxb[2]->Fill(qxb[j]);
          h_qyb[2]->Fill(qyb[j]);
          h_qxc[2]->Fill(qxc[j]);
          h_qyc[2]->Fill(qyc[j]);
          h_qxdimu[2]->Fill(qxdimu[j]);
          h_qydimu[2]->Fill(qydimu[j]);
        }
      }
    }
  }
  
  TCanvas *c_Q2D_A = new TCanvas("c_Q2D_A","",600,600);
  TCanvas *c_Q2D_B = new TCanvas("c_Q2D_B","",600,600);
  TCanvas *c_Q2D_C = new TCanvas("c_Q2D_C","",600,600);
  TCanvas *c_Q2D = new TCanvas("c_Q2D","",600,600);
  
  gPad->SetLeftMargin(0.17);
  gPad->SetLeftMargin(0.17);
  c_Q2D_A->cd();
  h_qa_2D->Draw("colz");

  gPad->SetLeftMargin(0.17);
  gPad->SetLeftMargin(0.17);
  c_Q2D_B->cd();
  h_qb_2D->Draw("colz");
  
  gPad->SetLeftMargin(0.17);
  gPad->SetLeftMargin(0.17);
  c_Q2D_C->cd();
  h_qc_2D->Draw("colz");
  
  gPad->SetLeftMargin(0.17);
  gPad->SetLeftMargin(0.17);
  c_Q2D->cd();
  h_q_2D->Draw("colz");

  c_Q2D_A->SaveAs("plots/c_Q2D_A.pdf");
  c_Q2D_B->SaveAs("plots/c_Q2D_B.pdf");
  c_Q2D_C->SaveAs("plots/c_Q2D_C.pdf");
  c_Q2D->SaveAs("plots/c_Q2D.pdf");


  TCanvas *c_qxA = new TCanvas("plots/c_qxA","",600,600);
  TCanvas *c_qyA = new TCanvas("plots/c_qyA","",600,600);
  TCanvas *c_qxB = new TCanvas("plots/c_qxB","",600,600);
  TCanvas *c_qyB = new TCanvas("plots/c_qyB","",600,600);
  TCanvas *c_qxC = new TCanvas("plots/c_qxC","",600,600);
  TCanvas *c_qyC = new TCanvas("plots/c_qyC","",600,600);
  TCanvas *c_qxdimu = new TCanvas("plots/c_qxdimu","",600,600);
  TCanvas *c_qydimu = new TCanvas("plots/c_qydimu","",600,600);


  TLegend *leg_qxa = new TLegend(0.6,0.6,0.9,0.9);
  TLegend *leg_qxb = new TLegend(0.6,0.6,0.9,0.9);
  TLegend *leg_qxc = new TLegend(0.6,0.6,0.9,0.9);
  TLegend *leg_qxdimu = new TLegend(0.6,0.6,0.9,0.9);
  TLegend *leg_qya = new TLegend(0.6,0.6,0.9,0.9);
  TLegend *leg_qyb = new TLegend(0.6,0.6,0.9,0.9);
  TLegend *leg_qyc = new TLegend(0.6,0.6,0.9,0.9);
  TLegend *leg_qydimu = new TLegend(0.6,0.6,0.9,0.9);
  SetLegendStyle(leg_qxa);
  SetLegendStyle(leg_qxb);
  SetLegendStyle(leg_qxc);
  SetLegendStyle(leg_qxdimu);
  SetLegendStyle(leg_qya);
  SetLegendStyle(leg_qyb);
  SetLegendStyle(leg_qyc);
  SetLegendStyle(leg_qydimu);

  for(int i=0; i<nMass_div; i++){
    SetHistStyle(h_qxa[i],i,i);
    SetHistStyle(h_qxb[i],i,i);
    SetHistStyle(h_qxc[i],i,i);
    SetHistStyle(h_qxdimu[i],i,i);
    SetHistStyle(h_qya[i],i,i);
    SetHistStyle(h_qyb[i],i,i);
    SetHistStyle(h_qyc[i],i,i);
    SetHistStyle(h_qydimu[i],i,i);

    scaleInt(h_qxa[i]);
    scaleInt(h_qxb[i]);
    scaleInt(h_qxc[i]);
    scaleInt(h_qya[i]);
    scaleInt(h_qyb[i]);
    scaleInt(h_qyc[i]);
    scaleInt(h_qxdimu[i]);
    scaleInt(h_qydimu[i]);

    leg_qxa->AddEntry(h_qxa[i],fSB[i].Data(),"lp");
    leg_qxb->AddEntry(h_qxb[i],fSB[i].Data(),"lp");
    leg_qxc->AddEntry(h_qxc[i],fSB[i].Data(),"lp");
    leg_qxdimu->AddEntry(h_qxdimu[i],fSB[i].Data(),"lp");
    leg_qya->AddEntry(h_qya[i],fSB[i].Data(),"lp");
    leg_qyb->AddEntry(h_qyb[i],fSB[i].Data(),"lp");
    leg_qyc->AddEntry(h_qyc[i],fSB[i].Data(),"lp");
    leg_qydimu->AddEntry(h_qydimu[i],fSB[i].Data(),"lp");

    c_qxA->cd();
    h_qxa[i]->Draw("P same");
    leg_qxa->Draw("same");
    c_qxB->cd();
    h_qxb[i]->Draw("P same");
    leg_qxb->Draw("same");
    c_qxC->cd();
    h_qxc[i]->Draw("P same");
    leg_qxc->Draw("same");
    c_qxdimu->cd();
    h_qxdimu[i]->Draw("P same");
    leg_qxdimu->Draw("same");

    c_qyA->cd();
    h_qya[i]->Draw("P same");
    leg_qya->Draw("same");
    c_qyB->cd();
    h_qyb[i]->Draw("P same");
    leg_qyb->Draw("same");
    c_qyC->cd();
    h_qyc[i]->Draw("P same");
    leg_qyc->Draw("same");
    c_qydimu->cd();
    h_qydimu[i]->Draw("P same");
    leg_qydimu->Draw("same");
  }

  c_qxA->SaveAs("plots/c_qxA.pdf");
  c_qxB->SaveAs("plots/c_qxB.pdf");
  c_qxC->SaveAs("plots/c_qxC.pdf");
  c_qxdimu->SaveAs("plots/c_qxdimu.pdf");
  c_qyA->SaveAs("plots/c_qyA.pdf");
  c_qyB->SaveAs("plots/c_qyB.pdf");
  c_qyC->SaveAs("plots/c_qyC.pdf");
  c_qydimu->SaveAs("plots/c_qydimu.pdf");

}
    

void GetHistSqrt(TH1D* h1, TH1D* h2){
  if(h1->GetNbinsX() != h2->GetNbinsX()){ cout << "Inconsistent # of bins b/w histograms !! " << endl;}
  double content;
  double err;
  for(int i=1; i<=h1->GetNbinsX(); i++){
    content=0;err=0;
    content = h1->GetBinContent(i);
    err = h1->GetBinError(i);
    err = 0.5*err*TMath::Power(content,-0.5);
    h2->SetBinContent(i,TMath::Sqrt(content));
    h2->SetBinError(i,err);
  }
} 
