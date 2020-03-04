#include <iostream>
#include "Style_jaebeom.h"
#include "commonUtility.h"
#include "cutsAndBinUpsilonV2.h"
#include "CMS_lumi_square.C"
#include "tdrstyle.C"

using namespace std;

void Figure3(){

  setTDRStyle();
  writeExtraText= false;
  int iPeriod = 2;
  int iPos = 33;
  
  const int nBin=4;
  const int nCentBin = 3;
  double exsys[nBin] =  {1.5, 1.5,22};

  double xmax = 50; double ymin = -0.10; double ymax = 0.165;
  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile *fSpectra = new TFile("Spectra_v2.root");
  double ptBin[nBin+1] = {0, 3, 6, 10, 50};
  double centBin[nCentBin+1] = {10,30,50,90};

	TGraphErrors* gv2[nCentBin];
	TGraphErrors* gv2_sys[nCentBin];
	TGraphErrors* gv2_int = (TGraphErrors*) fSpectra->Get("gr_Pt_CentInt");
	TGraphErrors* gv2_int_sys = (TGraphErrors*) fSpectra->Get("gr_sys_Pt_CentInt");

	TGraphErrors* gr_v2[nCentBin];
	TGraphErrors* gr_v2_sys[nCentBin];
	TGraphErrors* gr_v2_int = new TGraphErrors();
	TGraphErrors* gr_v2_int_sys = new TGraphErrors();
  
  double pxtmp, pytmp, extmp, eytmp, esy;
  double esx = xmax/(8*8);
  double esx_binned[nBin] = {1.5, 1.5, 2, 20};
  for(int icent = 0; icent<nCentBin; icent++){
    gv2[icent] = (TGraphErrors*) fSpectra->Get(Form("gr_Pt_Cent%.f%.f",centBin[icent],centBin[icent+1])); 
    gv2_sys[icent] = (TGraphErrors*) fSpectra->Get(Form("gr_sys_Pt_Cent%.f%.f",centBin[icent],centBin[icent+1])); 
    for(int ip = 0; ip<gv2[icent]->GetN(); ip++){
      pxtmp = 0; pytmp = 0; extmp = 0; eytmp = 0; esy=0;
      gv2[icent]->GetPoint(ip,pxtmp,pytmp);
      eytmp = gv2[icent]->GetErrorY(ip);
      esy = gv2_sys[icent]->GetErrorY(ip);
      gv2[icent]->SetPoint(ip, (ptBin[ip]+ptBin[ip+1])/2, pytmp);
      gv2[icent]->SetPointError(ip, 0, eytmp);
      gv2_sys[icent]->SetPoint(ip, (ptBin[ip]+ptBin[ip+1])/2, pytmp);
      gv2_sys[icent]->SetPointError(ip, esx_binned[ip], esy);
    }
    SetGraphStyle(gv2[icent],1,1);
    SetGraphStyleSys(gv2_sys[icent],1);
  }
  
  //Set for Integrated Bin
  for(int ip = 0; ip<gv2_int->GetN(); ip++){
    pxtmp = 0; pytmp = 0; extmp = 0; eytmp = 0; esy=0;
    gv2_int->GetPoint(ip,pxtmp,pytmp);
    eytmp = gv2_int->GetErrorY(ip);
    esy = gv2_int_sys->GetErrorY(ip);
    gv2_int->SetPoint(ip, pxtmp, pytmp);
    gv2_int->SetPointError(ip, 0, eytmp);
    gv2_int_sys->SetPoint(ip, pxtmp, pytmp);
    gv2_int_sys->SetPointError(ip, esx_binned[ip], esy);
    cout << "pxtmp : " << pxtmp << endl;
    cout << "pytmp : " << pytmp << endl;
    cout << endl;
  }
  SetGraphStyle(gv2_int,1,1);
  SetGraphStyleSys(gv2_int_sys,1);

  // Theory
  // pNRQCD
  TFile *th_file_pNRQCD = new TFile("Theory_v2_pNRQCD.root","read");
  TGraphErrors *g_v2_d6 = (TGraphErrors*) th_file_pNRQCD->Get("g_v2_d6");
  g_v2_d6->SetLineStyle(8);
  g_v2_d6->SetLineWidth(3);
  g_v2_d6->SetLineColor(kViolet-1);
  g_v2_d6->SetMarkerStyle(kOpenTriangleUp);
  
  // Xiaojun
  TFile *th_file_xiao = new TFile("Theory_v2_Xiaojun.root","read");
  TGraphErrors *g_v2 = (TGraphErrors*) th_file_xiao->Get("gv2_tot");
  
  g_v2->SetLineWidth(2);
  g_v2 -> SetLineColor(kRed);

  // Du Rapp
  TFile *th_file_DuRapp = new TFile("Theory_v2_DuRapp.root","read");
  TGraphErrors *g_v2TotMin_durapp = (TGraphErrors*) th_file_DuRapp->Get("g_v2TotMin");
  TGraphErrors *g_v2TotMax_durapp = (TGraphErrors*) th_file_DuRapp->Get("g_v2TotMax");
  TGraphErrors *gv2_tot_durapp = (TGraphErrors*) th_file_DuRapp->Get("gv2_tot"); 
  
  g_v2TotMin_durapp->SetLineWidth(1);
  g_v2TotMax_durapp->SetLineWidth(1);
  g_v2TotMin_durapp -> SetLineColor(kGreen+1);
  g_v2TotMax_durapp -> SetLineColor(kGreen+1);

  gv2_tot_durapp->SetFillStyle(3004);
  gv2_tot_durapp->SetFillColor(kGreen+2);
  gv2_tot_durapp->SetLineColor(kGreen+2);

  //Strickland
  TFile *th_file_strickland = new TFile("Theory_v2_aHydro.root","read");
  TGraphErrors *gv2_strickland = (TGraphErrors*) th_file_strickland->Get("g_v2_aHydro");
  gv2_strickland->SetLineWidth(3);
  gv2_strickland->SetLineColor(kOrange+5);
  gv2_strickland->SetLineStyle(2);


  TCanvas* c1[nCentBin];
  TLegend *leg[nCentBin];
  TLegend *leg_th[nCentBin];
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
  double lab_posx = 0.269; double lab_posy = 0.14; double lab_pos_diff = 0.186;
  for(int ic=0; ic<nCentBin; ic++){
    c1[ic]= new TCanvas(Form("c%d",ic),Form("c%d",ic),600,600);
    gPad->SetLeftMargin(0.19);
    gPad->SetBottomMargin(0.17);
    gPad->SetTopMargin(0.065);
    c1[ic]->cd();
    gv2[ic]->SetMarkerStyle(kFullCircle);
    gv2_sys[ic]->Draw("A5");
    gv2_sys[ic]->GetXaxis()->CenterTitle();
    gv2_sys[ic]->GetYaxis()->CenterTitle();
    gv2_sys[ic]->GetYaxis()->SetTitleOffset(1.5);
    gv2_sys[ic]->GetXaxis()->SetTitleOffset(1.1);
    gv2_sys[ic]->GetYaxis()->SetTitle("#it{v}_{2}^{#varUpsilon(1S)}");
    gv2_sys[ic]->GetXaxis()->SetTitle("p_{T}^{#varUpsilon(1S)} (GeV/c)");
    gv2_sys[ic]->GetXaxis()->SetTitleSize(0.055);
//    gv2_sys[ic]->GetXaxis()->SetTitleSize(0);
//    gv2_sys[ic]->GetXaxis()->SetLabelSize(0);
//    gv2_sys[ic]->GetXaxis()->SetTickLength(0);
    gv2_sys[ic]->GetXaxis()->SetRangeUser(0,xmax);
    gv2_sys[ic]->GetXaxis()->SetLimits(0,xmax);
    gv2_sys[ic]->SetMinimum(ymin);
    gv2_sys[ic]->SetMaximum(ymax);

    //Draw
    gv2[ic]->Draw("pe same");
    g_v2_d6->Draw("P");
    dashedLine(0.,0.,xmax,0.,1,1);

    //Text & Line
    TString perc = "%";
    leg[ic]= new TLegend(0.74, 0.70, 0.925, 0.77);
    SetLegendStyle(leg[ic]);
    leg[ic]->SetTextSize(0.044);
    leg[ic]->AddEntry(gv2[ic],"#varUpsilon(1S)","pe");
    leg[ic]->Draw("same");

    leg_th[ic] = new TLegend(0.54, 0.70, 0.725, 0.77);
    SetLegendStyle(leg_th[ic]);
    leg_th[ic]->SetTextSize(0.040);
    
    globtex->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
    globtex->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
    globtex->DrawLatex(0.23, sz_init-sz_step*2, Form("Cent. %.f-%.f%s",centBin[ic],centBin[ic+1],perc.Data()));
    lab_pos_diff = lab_pos_diff - 0.0002;
    gv2_sys[ic]->SetName(Form("gr_sys_v2_vs_pt_Cent%.f%.f",centBin[ic],centBin[ic+1]));
    gv2[ic]->SetName(Form("gr_point_v2_vs_pt_Cent%.f%.f",centBin[ic],centBin[ic+1]));
    CMS_lumi_square( c1[ic], iPeriod, iPos,1 );
    c1[ic]->Update();
    //c1[ic]->SaveAs(Form("compare_theory_pNRQCD_pt_Cent%d.pdf",ic));
  }

  TLatex* globtex_int = new TLatex();
  globtex_int->SetNDC();
  globtex_int->SetTextAlign(12); //left-center
  globtex_int->SetTextFont(42);
  globtex_int->SetTextSize(0.040);

  TCanvas* c2 = new TCanvas("cint","cint",600,600);
  gPad->SetLeftMargin(0.19);
  gPad->SetBottomMargin(0.17);
  gPad->SetTopMargin(0.065);
  c2->cd();
  gv2_int->SetMarkerStyle(kFullCircle);
  gv2_int_sys->Draw("A5");
  gv2_int_sys->GetXaxis()->CenterTitle();
  gv2_int_sys->GetYaxis()->CenterTitle();
  gv2_int_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2_int_sys->GetXaxis()->SetTitleOffset(1.1);
  gv2_int_sys->GetYaxis()->SetTitle("#it{v}_{2}^{#varUpsilon(1S)}");
  gv2_int_sys->GetXaxis()->SetTitle("p_{T}^{#varUpsilon(1S)} (GeV/c)");
  gv2_int_sys->GetXaxis()->SetTitleSize(0.055);
  gv2_int_sys->GetXaxis()->SetLabelSize(0.04);
  gv2_int_sys->GetYaxis()->SetNdivisions(506);
  gv2_int_sys->GetYaxis()->SetLabelSize(0.04);
//  gv2_int_sys->GetXaxis()->SetTickLength(0);
  gv2_int_sys->GetXaxis()->SetLimits(0,xmax);
  gv2_int_sys->GetXaxis()->SetRangeUser(0,xmax);
  ymin = -0.05;
  ymax = 0.15;
  gv2_int_sys->SetMinimum(ymin);
  gv2_int_sys->SetMaximum(ymax);
  
  gv2_int->Draw("pe same");

  g_v2_d6->Draw("L"); // draw pNRQCD
  
  g_v2->Draw("l"); // draw Xiaojun
 
  g_v2TotMin_durapp->Draw("l"); // draw Du Rapp
  g_v2TotMax_durapp->Draw("l"); // draw Du Rapp
  gv2_tot_durapp->Draw("f"); //draw Du Rapp
  
  gv2_strickland->Draw("L"); //draw aHydro


  TLegend *leg_int= new TLegend(0.74, 0.70, 0.925, 0.77);
  SetLegendStyle(leg_int);
  leg_int->SetTextSize(0.044);
  leg_int->AddEntry(gv2_int,"#varUpsilon(1S)","pe");
  leg_int->Draw("same");
  dashedLine(0.,0.,xmax,0.,1,1);
  
  globtex_int->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
  globtex_int->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
  globtex_int->DrawLatex(0.23, sz_init-sz_step*2, "Cent. 10-90%");

  TLegend *leg_th_int = new TLegend(0.41,0.55,0.62,0.72);
  SetLegendStyle(leg_th_int);
  leg_th_int->SetTextSize(0.030);
  leg_th_int->AddEntry(g_v2_d6," Hong, Lee (10-50%)","l");
  leg_th_int->AddEntry(g_v2," Yao (10-90%)","l");
  leg_th_int->AddEntry(gv2_tot_durapp," Du, Rapp (20-40%)","f");
  leg_th_int->AddEntry(gv2_strickland," Bhaduri, Borghini, (10-90%)","l");
  leg_th_int->AddEntry((TObject*)0," Jaiswal, Strickland","");
  leg_th_int->Draw("same");

  
  double lab_posx_ = 0.769; double lab_posy_ = 0.54; double lab_pos_diff_ = 0.186;
  TLatex* globtex_label_ = new TLatex();
  globtex_label_->SetNDC();
  globtex_label_->SetTextAlign(12); //left-center
  globtex_label_->SetTextFont(42);
  globtex_label_->SetTextSize(0.040);
  lab_pos_diff_ =lab_pos_diff_ - 0.0002;
  
 // globtex_label_->DrawLatex(lab_posx_,lab_posy_,"#varUpsilon(1S)"); 
  CMS_lumi_square( c2, iPeriod, iPos,1 );
	c2->Update();
  c2->SaveAs("Figure3.pdf");



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
  for(int i=0;i<nCentBin;i++){
    gsys[i]->Write();
    final_v2_g[i]->Write();
  }
*/

	return;
} 
