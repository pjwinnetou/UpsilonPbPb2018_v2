#include <iostream>
#include "commonUtility.h"
#include "cutsAndBinUpsilonV2.h"
#include "HiEvtPlaneList.h"
#include "Style_jaebeom.h"

using namespace std;
using namespace hi;

void makeV2Hist(int cLow = 20, int cHigh = 120,
                float ptLow = 0, float ptHigh = 30, 
                float yLow = 0, float yHigh=2.4,
                float SiMuPtCut = 4)
{
  TFile *rf = new TFile("OniaFlow_skim.root","read");
  TTree *tree = (TTree*) rf -> Get("mmep");

  DiMuon dm;
  TBranch *b_dm;
  tree -> SetBranchAddress("mm", &dm, &b_dm);
  
  Double_t mass;
  Double_t qxa, qya, qxb, qyb, qxc, qyc, qxdimu, qydimu;
  Double_t pt, y, pt1, pt2, eta1, eta2;
  Int_t cBin;
  Int_t event; 

  const int nMassBin = 30;
  double massLow = 8;
  double massHigh = 14;

  int bfevt =-1;
  int afevt =-1;
  
  TH1D* h1 = new TH1D("h1",";m_{#mu^{+}#mu^{-}};v_{2}^{S+B}",nMassBin,massLow,massHigh);
  TH1D* h2 = new TH1D("h2",";m_{#mu^{+}#mu^{-}};Counts",60,massLow,massHigh);

  double v2_1[nMassBin];
  double v2_2[nMassBin];
  double v2_3[nMassBin];
  double v2_4[nMassBin];
  double v2_1_avg[nMassBin];
  double v2_2_avg[nMassBin];
  double v2_3_avg[nMassBin];
  double v2_4_avg[nMassBin];
  int count[nMassBin] = {0};
  Int_t nEvt = tree->GetEntries();
  cout << "nEvt : " << nEvt << endl;
  for(int i=0; i<nEvt; i++){
    dm.clear();
    tree->GetEntry(i);
    bfevt = 0; afevt = 0;
    mass    = dm.mass;
    pt      = dm.pt;
    pt1     = dm.pt1;
    pt2     = dm.pt2;
    eta1    = dm.eta1;
    eta2    = dm.eta2;
    y       = dm.y;
    cBin    = dm.cBin;
    qxa     = dm.qxa;
    qxb     = dm.qxb;
    qxc     = dm.qxc;
    qya     = dm.qya;
    qyb     = dm.qyb;
    qyc     = dm.qyc;
    qxdimu  = dm.qxdimu;
    qydimu  = dm.qydimu;
    event   = dm.event;

    for(int j=0; j<nMassBin; j++){
      if(pt>ptLow&&pt<ptHigh&&abs(y)<yHigh&&abs(y)>yLow&&pt1>SiMuPtCut&&pt2>SiMuPtCut&&abs(eta1)<2.4&&abs(eta2)<2.4 && cBin>cLow&&cBin<cHigh)
      {
        bfevt = event;
        tree->GetEntry(i+1);
        if(pt>ptLow&&pt<ptHigh&&abs(y)<yHigh&&abs(y)>yLow&&pt1>SiMuPtCut&&pt2>SiMuPtCut&&abs(eta1)<2.4&&abs(eta2)<2.4 && cBin>cLow&&cBin<cHigh) afevt = dm.event;
        if(bfevt == afevt) continue;

        if(mass >= (massLow+(massHigh-massLow)/nMassBin*j) && mass < (massLow+(massHigh-massLow)/nMassBin*(j+1)))
        {
          v2_1[j] = qxa*qxdimu + qya*qydimu;
          v2_1_avg[j] += v2_1[j];  
          v2_2[j] = qxa*qxb + qya*qyb;
          v2_2_avg[j] += v2_2[j];  
          v2_3[j] = qxa*qxc + qya*qyc;
          v2_3_avg[j] += v2_3[j];  
          v2_4[j] = qxb*qxc + qyb*qyc;
          v2_4_avg[j] += v2_4[j];  
          count[j]++;
        }
        h2->Fill(mass);
      }
    }
  }

  double v2_final[nMassBin]; 
  
  cout << "Number of upsilon candidate in given bin " << endl;
  for(int ibin=0; ibin<nMassBin; ibin++){
    
    cout << ibin << "th Bin : " << count[ibin] << endl;

    v2_1_avg[ibin] = v2_1_avg[ibin]/count[ibin];
    v2_2_avg[ibin] = v2_2_avg[ibin]/count[ibin];
    v2_3_avg[ibin] = v2_3_avg[ibin]/count[ibin];
    v2_4_avg[ibin] = v2_4_avg[ibin]/count[ibin];

    v2_final[ibin] = v2_1_avg[ibin]/TMath::Sqrt(v2_2_avg[ibin]*v2_3_avg[ibin]/v2_4_avg[ibin]);
    cout << "v2_1_avg " << ibin << " : " << v2_1_avg[ibin] << endl;
    cout << "v2_2_avg " << ibin << " : " << v2_2_avg[ibin] << endl;
    cout << "v2_3_avg " << ibin << " : " << v2_3_avg[ibin] << endl;
    cout << "v2_4_avg " << ibin << " : " << v2_4_avg[ibin] << endl;
    h1->SetBinContent(ibin+1,v2_final[ibin]);
    h1->SetBinError(ibin+1,0);
  }

  SetHistStyle(h1,0,0);
  SetHistStyle(h2,0,0);
  h2->GetYaxis()->SetLimits(0,14000);
//  h2->GetYaxis()->SetRangeUser(0,14000);
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","",600,600);
  gPad->SetLeftMargin(0.17);
  TCanvas *c2 = new TCanvas("c2","",600,600);
  gPad->SetLeftMargin(0.17);
  c1->cd();
  h1->Draw("P");
  c2->cd();
  h2->Draw("P");

  c1->SaveAs("v2_SplusB_test.pdf");
  c2->SaveAs("MassDist.pdf");
}
    


