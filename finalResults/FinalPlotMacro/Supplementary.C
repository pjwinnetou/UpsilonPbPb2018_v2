#include <iostream>
#include "Style_jaebeom.h"
#include "commonUtility.h"
#include "cutsAndBinUpsilonV2.h"
#include "CMS_lumi_square.C"
#include "tdrstyle.C"

using namespace std;

void Supplementary(){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  
  const int nBin=3;
  const int nYBin = 2;
  double exsys[nBin] =  {1.5, 1.5,4.5};

  double xmax = 6; double ymin = -0.08; double ymax = 0.23;
  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile *fSpectra = new TFile("Spectra_v2.root");
  double ptBin[nBin+1] = {0, 3, 6, 15};

  TFile *fALICE = new TFile("HEPData-ins1742764-v1-Table_2.root");
  TH1D* hpy = (TH1D*) fALICE->Get("Table 2/Hist1D_y1");
  TH1D* hpy_stat = (TH1D*) fALICE->Get("Table 2/Hist1D_y1_e1");
  TH1D* hpy_sys = (TH1D*) fALICE->Get("Table 2/Hist1D_y1_e2");

	TGraphErrors* gv2_int = (TGraphErrors*) fSpectra->Get("gr_Pt_YIntC");
	TGraphErrors* gv2_int_sys = (TGraphErrors*) fSpectra->Get("gr_sys_Pt_YIntC");

  TGraphErrors* gv2_alice = new TGraphErrors();
  TGraphErrors* gv2_alice_sys = new TGraphErrors();

  double pxtmp, pytmp, extmp, eytmp, esy, relerr;
  double esx = xmax/(3*8);
  double posx_shift = esx;
  for(int ip = 0; ip<gv2_int->GetN(); ip++){
    pxtmp = 0; pytmp = 0; extmp = 0; eytmp = 0; esy=0; relerr=0;
    gv2_int->GetPoint(ip,pxtmp,pytmp);
    eytmp = gv2_int->GetErrorY(ip);
    relerr = gv2_int_sys->GetErrorY(ip);
    gv2_int->SetPoint(ip, xmax/(nBin*2)+xmax/nBin*ip - posx_shift, pytmp);
    gv2_int->SetPointError(ip, 0, eytmp);
    gv2_int_sys->SetPoint(ip, xmax/(nBin*2)+xmax/nBin*ip - posx_shift, pytmp);
    gv2_int_sys->SetPointError(ip, esx/2, relerr);
  }

  for(int ip=0; ip<hpy->GetNbinsX(); ip++){
    gv2_alice->SetPoint(ip,xmax/(nBin*2)+xmax/nBin*ip + posx_shift,hpy->GetBinContent(ip+1));
    gv2_alice_sys->SetPoint(ip, xmax/(nBin*2)+xmax/nBin*ip + posx_shift,hpy->GetBinContent(ip+1));
    gv2_alice->SetPointError(ip,0,hpy_stat->GetBinContent(ip+1));
    gv2_alice_sys->SetPointError(ip, esx/2,hpy_sys->GetBinContent(ip+1));
  }
  
  double ey_stat = gv2_int->GetErrorY(2);
  double ey_sys = gv2_int_sys->GetErrorY(2);
  double ey_tot = TMath::Sqrt(ey_stat*ey_stat + ey_sys*ey_sys);
  double y_p,x_p;
  gv2_int->GetPoint(2,x_p,y_p);
  cout << y_p << ", " << ey_tot << endl;
  
  SetGraphStyle2(gv2_int,2,0);
  SetGraphStyleSys2(gv2_int_sys,2);
  gv2_int_sys-> SetFillColorAlpha(kBlue+1,0.3);
  SetGraphStyleOpen(gv2_alice,4,0,0);
  SetGraphStyleSys2(gv2_alice_sys,0);
  gv2_alice_sys->SetFillColorAlpha(kGray+3,0.2);

  TLatex* globtex_label = new TLatex();
  globtex_label->SetNDC();
  globtex_label->SetTextAlign(12); //left-center
  globtex_label->SetTextFont(42);
  globtex_label->SetTextSize(0.040);
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  double sz_init = 0.867; double sz_step = 0.0535; double ymin_d = 0.007;
  double lab_posx = 0.282; double lab_posy = 0.14; double lab_pos_diff = 0.269;
  TCanvas *c1= new TCanvas("c1","c1",600,600);
  gPad->SetLeftMargin(0.19);
  gPad->SetBottomMargin(0.17);
  gPad->SetTopMargin(0.065);
  c1->cd();
  gv2_int_sys->SetMarkerStyle(kFullCircle);
  gv2_int_sys->GetXaxis()->CenterTitle();
  gv2_int_sys->GetYaxis()->CenterTitle();
  gv2_int_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2_int_sys->GetXaxis()->SetTitleOffset(1.1);
  gv2_int_sys->GetYaxis()->SetTitle("#it{v}_{2}");
  gv2_int_sys->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  gv2_int_sys->GetXaxis()->SetTitleSize(0.055);
  gv2_int_sys->GetXaxis()->SetLabelSize(0);
  gv2_int_sys->GetXaxis()->SetTickLength(0);
  gv2_int_sys->GetXaxis()->SetRangeUser(0,xmax);
  gv2_int_sys->GetXaxis()->SetLimits(0,xmax);
  gv2_int_sys->SetMinimum(ymin);
  gv2_int_sys->SetMaximum(ymax);

  gv2_int_sys->Draw("A5");
  gv2_alice_sys->Draw("5");
  gv2_int->Draw("P");
  gv2_alice->Draw("P");
  dashedLine(0.,0.,xmax,0.,1,1);

  TString perc = "%";
  TLegend* leg= new TLegend(0.225, 0.59, 0.410, 0.77);
  SetLegendStyle(leg);
  leg->SetTextSize(0.042);
  leg->AddEntry(gv2_int,"|y| < 2.4","pe");
  leg->AddEntry(gv2_alice,"2.5 < |y| < 4 (ALICE arXiv:1907.03169) ","pe");
  leg->Draw("same");
  solidLine(xmax/3,ymin,xmax/3,ymin+ymin_d,1,1);
  solidLine(xmax/3*2,ymin,xmax/3*2,ymin+ymin_d,1,1);

  solidLine(xmax/3,ymax,xmax/3,ymax-ymin_d,1,1);
  solidLine(xmax/3*2,ymax,xmax/3*2,ymax-ymin_d,1,1);

  globtex->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
//  globtex->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
  globtex->DrawLatex(0.23, sz_init-sz_step, "Cent. 5-60 %");
  globtex_label->DrawLatex(lab_posx, lab_posy, "0-3");
  globtex_label->DrawLatex(lab_posx+lab_pos_diff,  lab_posy,"3-6");
  globtex_label->DrawLatex(lab_posx+lab_pos_diff*2, lab_posy, "6-15");
  
  double lab_posx_ = 0.812; double lab_posy_ = 0.724; double lab_pos_diff_ = 0.186;
  TLatex* globtex_int = new TLatex();
  globtex_int->SetNDC();
  globtex_int->SetTextAlign(12); //left-center
  globtex_int->SetTextFont(42);
  globtex_int->SetTextSize(0.040);
  globtex_int->DrawLatex(lab_posx_,lab_posy_,"#bf{#varUpsilon(1S)}"); 
  
  CMS_lumi_square( c1, iPeriod, iPos,1);
  c1->Update();
  c1->SaveAs("v2Sig_pt_CompALICE.pdf");


  /*
     TFile *fFullInt1S = new TFile("../SimFit/SigAllFreeFitFix/FitResult/SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality0-180_m8-14_OS.root","read");
     TFile *fFullInt2S = new TFile("../SimFit/SigAllFreeFitFix/FitResult/2S_SimFitResult_pt0.0-50.0_y0.0-2.4_muPt3.5_centrality0-180_m8-14_OS.root","read");

  TGraphErrors* gInt1S = (TGraphErrors*) fFullInt1S -> Get("v2vspt");  
  TGraphErrors* gInt2S = (TGraphErrors*) fFullInt2S -> Get("v2vspt");  
  double xInt1S, yInt1S, xInt2S, yInt2S;
  double xInt1SErr, yInt1SErr, xInt2SErr, yInt2SErr;
  
  gInt1S->GetPoint(0, xInt1S, yInt1S);
  gInt2S->GetPoint(0, xInt2S, yInt2S);
  yInt1SErr= gInt1S->GetErrorY(0);
  yInt2SErr= gInt2S->GetErrorY(0);

  TFile* fSys2S = new TFile("../Systematic/merged_sys_2S.root","read");
  TH1D* hInt1S = (TH1D*) fSys->Get("hSys_int_merged");
  TH1D* hInt2S = (TH1D*) fSys2S->Get("hSys_int_merged");
  
  double yInt1Ssys = hInt1S->GetBinContent(1);
  double yInt2Ssys = hInt2S->GetBinContent(1);

  cout << "Integrated 1S : " << Form("%.3f",yInt1S) << " +/- " << Form("%.3f",yInt1SErr) << " +/- " << Form("%.3f",yInt1Ssys) << endl;
  cout << "Integrated 2S : " << Form("%.3f",yInt2S) << " +/- " << Form("%.3f",yInt2SErr) << " +/- " << Form("%.3f",yInt2Ssys) << endl;

  TFile *out = new TFile("out_v2_vs_pt.root","recreate");
  out->cd();
  gsys_int->SetName("gr_sys_v2_vs_pt_Cent090");
  final_v2_int_g->SetName("gr_point_v2_vs_pt_Cent090");
  gsys_int->Write(); 
  final_v2_int_g->Write();
  for(int i=0;i<nYBin;i++){
    gsys[i]->Write();
    final_v2_g[i]->Write();
  }
*/

	return;
} 
