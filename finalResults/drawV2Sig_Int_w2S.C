#include <iostream>
#include "../Style_jaebeom.h"
#include "../commonUtility.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

using namespace std;

void drawV2Sig_Int_w2S(){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  
  const int nBin=4;
  double exsys[nBin] =  {1.5, 1.5,22};
  double xmax = 90; 

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nBin];
  TFile* fIn_int[nBin];
  double centBin[nBin+1] = {0, 10, 30, 50, 90};

	TGraphErrors* gv2[nBin];
	TGraphErrors* gv2_int[nBin];
  for (int is=0; is<nBin; is++){
  	fIn_int[is] = new TFile(Form("../SimFit/SigAllFreeFitFix/FitResult/SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",centBin[is]*2,centBin[is+1]*2),"READ");
    gv2_int[is]=(TGraphErrors*)fIn_int[is]->Get("v2vspt");
  }

  TH1D *final_v2_int = new TH1D("Final_v2_int",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,centBin);
  

  //// get bin width and calculate systematic uncertainties
  double pxtmp_int, pytmp_int, extmp_int, eytmp_int;
  double relsys_int;

  for (int is=0; is<nBin; is++){
    cout << is+1 <<"th bin***************" << endl;
    pxtmp_int=0; pytmp_int=0; extmp_int=0; eytmp_int=0; relsys_int=0;
    gv2_int[is] -> GetPoint(0,pxtmp_int,pytmp_int);
    eytmp_int = gv2_int[is]->GetErrorY(0);
    final_v2_int->SetBinContent(is+1, pytmp_int);
    final_v2_int->SetBinError(is+1, eytmp_int);
  }

  final_v2_int->GetXaxis()->CenterTitle();
  final_v2_int->GetYaxis()->CenterTitle();
  final_v2_int->GetYaxis()->SetTitleOffset(1.2);
  final_v2_int->GetXaxis()->SetTitleOffset(1.);
  final_v2_int->SetMinimum(-0.05);
  final_v2_int->SetMaximum(0.2);

  TGraphErrors* final_v2_int_g = new TGraphErrors(final_v2_int);
  TFile* fSys = new TFile("../Systematic/merged_sys.root","read");
  TH1D* hsys_int = (TH1D*) fSys->Get("hCent_Pt050_merged");
  TGraphErrors *gsys_int = new TGraphErrors();
  for(int ib=0;ib<final_v2_int->GetNbinsX();ib++){
    gsys_int->SetPoint(ib, final_v2_int->GetBinCenter(ib+1), final_v2_int->GetBinContent(ib+1));
    gsys_int->SetPointError(ib, final_v2_int->GetBinWidth(ib+1)/2, hsys_int->GetBinContent(ib+1));
    final_v2_int_g->SetPointError(ib,0,final_v2_int->GetBinError(ib+1));
  }
  SetGraphStyle(final_v2_int_g,1,1);
  SetGraphStyleSys(gsys_int,1);
  
  TFile *fFullInt2S = new TFile("../SimFit/SigAllFreeFitFix/FitResult/2S_SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality0-180_m8-14_OS.root","read");
  
  TGraphErrors* gInt2S = (TGraphErrors*) fFullInt2S -> Get("v2vspt");  
  double xInt2S, yInt2S;
  double xInt2SErr, yInt2SErr;
  
  gInt2S->GetPoint(0, xInt2S, yInt2S);
  yInt2SErr= gInt2S->GetErrorY(0);

  TFile* fSys2S = new TFile("../Systematic/merged_sys_2S.root","read");
  TH1D* hInt2S = (TH1D*) fSys2S->Get("hSys_int_merged");
  double yInt2Ssys = hInt2S->GetBinContent(1);

  TGraphErrors* final_v2_int_g_2S = new TGraphErrors();
  TGraphErrors* gsys_int_2S = new TGraphErrors();
  final_v2_int_g_2S->SetPoint(0,xmax/2,yInt2S);
  final_v2_int_g_2S->SetPointError(0,0,yInt2SErr);
  gsys_int_2S->SetPoint(0,xmax/2,yInt2S);
  gsys_int_2S->SetPointError(0,xmax/2,yInt2Ssys);
  gsys_int_2S->GetXaxis()->CenterTitle();
  gsys_int_2S->GetYaxis()->CenterTitle();
  gsys_int_2S->GetYaxis()->SetTitleOffset(1.2);
  gsys_int_2S->GetXaxis()->SetTitleOffset(1.);
  gsys_int_2S->GetYaxis()->SetTitle("v_{2}^{sig}");
  gsys_int_2S->GetXaxis()->SetTitle("Centrality (%)");
  gsys_int_2S->GetXaxis()->SetRangeUser(0,xmax);
  gsys_int_2S->SetMinimum(-0.05);
  gsys_int_2S->SetMaximum(0.15);
  SetGraphStyle(final_v2_int_g_2S,2,2);
  SetGraphStyleSys(gsys_int_2S,2);

 
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  gPad->SetLeftMargin(0.19);
  gPad->SetBottomMargin(0.17);
  gPad->SetTopMargin(0.065);
  c2->cd();
  final_v2_int_g->SetMarkerStyle(kFullCircle);
  gsys_int->Draw("A5");
  gsys_int->GetXaxis()->CenterTitle();
  gsys_int->GetYaxis()->CenterTitle();
  gsys_int->GetYaxis()->SetTitleOffset(1.2);
  gsys_int->GetXaxis()->SetTitleOffset(1.);
  gsys_int->GetYaxis()->SetTitle("v_{2}^{sig}");
  gsys_int->GetXaxis()->SetTitle("Centrality (%)");
  gsys_int->GetXaxis()->SetRangeUser(0,xmax);
  gsys_int->SetMinimum(-0.05);
  gsys_int->SetMaximum(0.15);
  gsys_int_2S->Draw("5");

  final_v2_int_g->Draw("pe same");
  final_v2_int_g_2S->Draw("pe same");

  TLegend *leg_int= new TLegend(0.74, 0.68, 0.925, 0.78);
  SetLegendStyle(leg_int);
  leg_int->SetTextSize(0.038);
  leg_int->AddEntry(final_v2_int_g,"#Upsilon(1S)","pe");
  leg_int->AddEntry(final_v2_int_g_2S,"#Upsilon(2S)","pe");
  leg_int->Draw("same");
  dashedLine(0.,0.,xmax,0.,1,1);
  
  double sz_init = 0.875; double sz_step = 0.0535;
  TLatex* globtex_int = new TLatex();
  globtex_int->SetNDC();
  globtex_int->SetTextAlign(12); //left-center
  globtex_int->SetTextFont(42);
  globtex_int->SetTextSize(0.040);
  globtex_int->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
  globtex_int->DrawLatex(0.23, sz_init-sz_step, "p_{T}^{#Upsilon} < 50 GeV/c");
  globtex_int->DrawLatex(0.23, sz_init-sz_step*2, "|y| < 2.4");
  
  CMS_lumi_square( c2, iPeriod, iPos );
	c2->Update();
  c2->SaveAs("v2Sig_int_w2S.pdf");

	return;
}
