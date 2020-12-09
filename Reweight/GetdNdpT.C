#include <iostream>
#include "../commonUtility.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../HiEvtPlaneList.h"
#include "../Style_jaebeom.h"
#include "../tdrstyle.C"
#include "../rootFitHeaders.h"
#include "../CMS_lumi_v2mass.C"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"

using namespace std;

void GetdNdpT(int state=2)
{
  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile* rf = new TFile(Form("Spectra_hist%dS.root",state));
  TH1D* hptMC = (TH1D*) rf->Get("hptMC"); 
  TH1D* hptDATA = (TH1D*) rf->Get("hptDATA"); 
  TH1D* hCentMC = (TH1D*) rf->Get("hCentMC"); 
  TH1D* hCentDATA = (TH1D*) rf->Get("hCentDATA"); 

  TH1D* hptRatio = (TH1D*) hptDATA->Clone("hptRatio");
  hptRatio->Divide(hptMC);
  TH1D* hCentRatio = (TH1D*) hCentDATA->Clone("hCentRatio");
  hCentRatio->Divide(hCentMC);

  SetHistStyleSmall(hptDATA,0,0);
  SetHistStyleSmall(hptMC,1,1);
  SetHistStyleSmall(hCentDATA,0,0);
  SetHistStyleSmall(hCentMC,1,1);
  hptDATA->SetLineWidth(2);
  hCentDATA->SetLineWidth(2);
  hptMC->SetLineWidth(2);
  hCentMC->SetLineWidth(2);


  double poslegx1 = 0.7; double poslegy1 = 0.5; 
  double poslegx2 = 0.9; double poslegy2 = 0.7;

  TLegend *legPt = new TLegend(poslegx1, poslegy1, poslegx2, poslegy2);
  SetLegendStyle(legPt);
  legPt->AddEntry(hptDATA,"DATA","le");  
  legPt->AddEntry(hptMC,"MC","le");  
  
  TLegend *legCent = new TLegend(poslegx1, poslegy1, poslegx2, poslegy2);
  SetLegendStyle(legCent);
  legCent->AddEntry(hptDATA,"DATA","le");  
  legCent->AddEntry(hptMC,"MC","le");  

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TPad* padPt1 = new TPad("padPt1","padPt1",0,0.31,1.0,1.0);
  TPad* padPt2 = new TPad("padPt2","padPt2",0,0.0,1.0,0.31);
  padPt1->SetBottomMargin(0.016);
  padPt1->SetLeftMargin(0.17);
  padPt1->SetTopMargin(0.08);
  padPt1->Draw();
  padPt1->cd();
  hptDATA->Draw();
  hptMC->Draw("same");
  legPt->Draw("same");

  double ymax = hptDATA->GetBinContent(hptDATA->GetMaximumBin());
  double ymax_ = hptMC->GetBinContent(hptMC->GetMaximumBin());
  if(ymax < ymax_) ymax = ymax_;

  hptDATA->GetYaxis()->SetLimits(0,ymax*1.2);
  hptDATA->GetYaxis()->SetRangeUser(0,ymax*1.5);
  hptDATA->GetXaxis()->SetTitle("");
  hptDATA->GetXaxis()->SetTitleSize(0);
  hptDATA->GetXaxis()->SetLabelSize(0);

  float pos_x = 0.33; float pos_x_mass = 0.53;
  float pos_y = 0.76; float pos_y_diff = 0.071;
  int text_color = 1; float text_size = 16;
  drawText(Form("#varUpsilon(%dS)",state),pos_x,pos_y,text_color,text_size);
  drawText("Centrality 10-90%",pos_x,pos_y - pos_y_diff,text_color,text_size);
  
  padPt2->SetTopMargin(0);
  padPt2->SetBottomMargin(0.33);
  padPt2->SetLeftMargin(0.17);
  padPt2->SetTicks();
  padPt2->cd();
  hptRatio->Draw();
  TF1 *fitRatio_;
  if(state==1){
   fitRatio_ = new TF1("fitRatio_","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,50);
   fitRatio_->SetParameters( 0.55074e+02, -1.34016e+01, 3.42256e+01, -2.81048e+00 );
  }
  else if(state==2){
    fitRatio_ = new TF1("fitRatio_","([0] + [1]*x)",0,50);
    fitRatio_ -> SetParameters(0.213, 10.2);
  }
  unsigned int nfpxl = 1100;
  fitRatio_->SetNpx(nfpxl);
  hptRatio->Fit(fitRatio_,"REIM","",0,50);



  TF1 *fitRatio;
  if(state==1){
    fitRatio = new TF1("fitRatio","( [0] + [1]*x + [2]*x*x ) / (  (x-[3])*(x-[3])*(x-[3])  )",0,50);
    for(int i=0; i<4; i++){
      fitRatio->FixParameter(i, fitRatio_->GetParameter(i));
      cout << "par " << i << " : " << fitRatio_->GetParameter(i) << endl;
    }
  }
  else if(state==2){
    fitRatio = new TF1("fitRatio","( [0] + [1]*x)",0,50);
    for(int i=0; i<2; i++){
      fitRatio->FixParameter(i, fitRatio_->GetParameter(i));
      cout << "par " << i << " : " << fitRatio_->GetParameter(i) << endl;
    }
  }

  hptRatio->GetYaxis()->SetLimits(0.1,3);
  hptRatio->GetYaxis()->SetRangeUser(0.1,2.96);
  hptRatio->GetYaxis()->SetTitle("DATA/MC"); 
  hptRatio->GetYaxis()->SetTitleSize(0.122); 
  hptRatio->GetYaxis()->SetTitleOffset(0.63); 
  hptRatio->GetYaxis()->SetNdivisions(503); 
  hptRatio->GetYaxis()->SetLabelSize(0.12);
  hptRatio->GetYaxis()->SetTickSize(0.02);
  
  hptRatio->GetXaxis()->SetTitleSize(0.12); 
  hptRatio->GetXaxis()->SetTitleOffset(1.09); 
  hptRatio->GetXaxis()->SetLabelSize(0.12);
  hptRatio->SetMarkerSize(0.8);
  hptRatio->SetLineWidth(2);
  


  padPt2->Update();
  jumSun(0,1,50,1,1,1);

  CMS_lumi_v2mass(padPt1,iPeriod,iPos,1);  
  padPt1->Update();
  padPt2->Update();
  c1->cd();
  padPt1->Draw();
  padPt2->Draw();
  c1->SaveAs(Form("dNdpT_DATAvsMC_%dS.pdf",state));

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  TPad* padCent1 = new TPad("padCent1","padCent1",0,0.31,1.0,1.0);
  TPad* padCent2 = new TPad("padCent2","padCent2",0,0.0,1.0,0.31);
  padCent1->SetBottomMargin(0.016);
  padCent1->SetLeftMargin(0.17);
  padCent1->SetTopMargin(0.08);
  padCent1->Draw();
  padCent1->cd();
  hCentDATA->Draw();
  hCentMC->Draw("same");
  legCent->Draw("same");

  hCentDATA->GetYaxis()->SetTitle("Normalized Counts");
  hCentDATA->GetXaxis()->SetTitle("");
  hCentDATA->GetXaxis()->SetLimits(0,180);
  hCentDATA->GetXaxis()->SetTitleSize(0);
  hCentDATA->GetXaxis()->SetLabelSize(0);
  
  ymax = hCentDATA->GetBinContent(hCentDATA->GetMaximumBin());
  ymax_ = hCentMC->GetBinContent(hCentMC->GetMaximumBin());
  if(ymax < ymax_) ymax = ymax_;

  hCentDATA->GetYaxis()->SetLimits(0,ymax*1.2);
  hCentDATA->GetYaxis()->SetRangeUser(0,ymax*1.15);

  drawText("p_{T} < 50 GeV/c",pos_x_mass,pos_y,text_color,text_size);
  
  padCent2->SetTopMargin(0);
  padCent2->SetBottomMargin(0.33);
  padCent2->SetLeftMargin(0.17);
  padCent2->SetTicks();
  padCent2->cd();
  hCentRatio->Draw();
  hCentRatio->GetYaxis()->SetLimits(0.1,3);
  hCentRatio->GetYaxis()->SetRangeUser(0,4.46);
  hCentRatio->GetXaxis()->SetLimits(0,180);
  hCentRatio->GetYaxis()->SetTitle("DATA/MC"); 
  hCentRatio->GetYaxis()->SetTitleSize(0.122); 
  hCentRatio->GetYaxis()->SetTitleOffset(0.63); 
  hCentRatio->GetYaxis()->SetNdivisions(505); 
  hCentRatio->GetYaxis()->SetLabelSize(0.12);
  hCentRatio->GetYaxis()->SetTickSize(0.02);
  
  hCentRatio->GetXaxis()->SetTitleSize(0.12); 
  hCentRatio->GetXaxis()->SetTitleOffset(1.09); 
  hCentRatio->GetXaxis()->SetLabelSize(0.12);
  hCentRatio->SetMarkerSize(0.8);
  hCentRatio->SetLineWidth(2);


  jumSun(20,1,180,1,1,1);

  CMS_lumi_v2mass(padCent1,iPeriod,iPos,1);  
  padCent1->Update();
  padCent2->Update();
  c2->cd();
  padCent1->Draw();
  padCent2->Draw();
  c2->SaveAs(Form("Cent_DATAvsMC_%dS.pdf",state));

  TFile *wf = new TFile(Form("WeightedFunc/Func_dNdpT_%dS.root",state),"recreate");
  wf->cd();
  fitRatio->Write();
  hptRatio->Write();

} 
