#include <iostream>
#include "Style_jaebeom.h"
#include "commonUtility.h"
#include "cutsAndBinUpsilonV2.h"
#include "CMS_lumi_square.C"
#include "tdrstyle.C"

using namespace std;

void Figure2rightand4(){

  setTDRStyle();
  writeExtraText= false;
  int iPeriod = 2;
  int iPos = 33;
  
  const int nBin=4;
  const int nCentBin = 3;
  double exsys[nBin] =  {1.5, 1.5,22};
  TString perc = "%";

  double xmax = 6; double ymin = -0.10; double ymax = 0.165;
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
  double esx = xmax/(3*8);
  for(int icent = 0; icent<nCentBin; icent++){
    gv2[icent] = (TGraphErrors*) fSpectra->Get(Form("gr_Pt_Cent%.f%.f",centBin[icent],centBin[icent+1])); 
    gv2_sys[icent] = (TGraphErrors*) fSpectra->Get(Form("gr_sys_Pt_Cent%.f%.f",centBin[icent],centBin[icent+1])); 
    for(int ip = 0; ip<gv2[icent]->GetN(); ip++){
      pxtmp = 0; pytmp = 0; extmp = 0; eytmp = 0; esy=0;
      gv2[icent]->GetPoint(ip,pxtmp,pytmp);
      eytmp = gv2[icent]->GetErrorY(ip);
      esy = gv2_sys[icent]->GetErrorY(ip);
      gv2[icent]->SetPoint(ip, xmax/(nBin*2)+xmax/nBin*ip, pytmp);
      gv2[icent]->SetPointError(ip, 0, eytmp);
      gv2_sys[icent]->SetPoint(ip, xmax/(nBin*2)+xmax/nBin*ip, pytmp);
      gv2_sys[icent]->SetPointError(ip, esx, esy);
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
    gv2_int->SetPoint(ip, xmax/(nBin*2)+xmax/nBin*ip, pytmp);
    gv2_int->SetPointError(ip, 0, eytmp);
    gv2_int_sys->SetPoint(ip, xmax/(nBin*2)+xmax/nBin*ip, pytmp);
    gv2_int_sys->SetPointError(ip, esx, esy);
    cout << "ip : " << ip << ", " << pytmp << " +/- " << eytmp << " +/- " << esy << endl;
  }
  SetGraphStyle(gv2_int,1,1);
  SetGraphStyleSys(gv2_int_sys,1);

  /*double ey_stat = gv2_int->GetErrorY(2);
  double ey_sys = gv2_int_sys->GetErrorY(2);
  double ey_tot = TMath::Sqrt(ey_stat*ey_stat + ey_sys*ey_sys);
  double y_p,x_p;
  gv2_int->GetPoint(2,x_p,y_p);
  cout << y_p << ", " << ey_tot << endl;
  */
  double ey_stat = gv2[1]->GetErrorY(2);
  double ey_sys = gv2_sys[1]->GetErrorY(2);
  double ey_tot = TMath::Sqrt(ey_stat*ey_stat + ey_sys*ey_sys);
  double y_p,x_p;
  gv2[1]->GetPoint(2,x_p,y_p);
  cout << y_p << ", " << ey_tot << endl;

  TPad* pad1[nCentBin];
  TCanvas* c1 = new TCanvas("c1","",1400,600);
  
  TLegend *leg = new TLegend(0.74, 0.69, 1.1, 0.77);
  SetLegendStyle(leg);
  leg->SetTextSize(0.040);
  
  TLatex* globtex_label = new TLatex();
  globtex_label->SetNDC();
  globtex_label->SetTextAlign(12); //left-center
  globtex_label->SetTextFont(42);
  globtex_label->SetTextSize(0.040);
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.038);
  double sz_init = 0.867; double sz_step = 0.0535; double ymin_d = 0.007;
  double lab_posx = 0.269; double lab_posy = 0.14; double lab_pos_diff = 0.196;
  double pad_x_end1 = 1./3 + 0.04;
  double pad_diff = (1-pad_x_end1)/2; 
  double pad_x1; double pad_x2;

  for(int ic=0; ic<nCentBin; ic++){
    if(ic==0){pad_x1 = 0;  pad_x2=pad_x_end1;}
    else if(ic!=0){ pad_x1 = pad_x_end1 + pad_diff*(ic-1); pad_x2=pad_x_end1 + pad_diff*ic;}
    pad1[ic] = new TPad(Form("pad1_%d",ic),"", pad_x1, 0, pad_x2, 1);
    gv2[ic]->SetMarkerStyle(kFullCircle);
    c1->cd();
    pad1[ic]->Draw();
    pad1[ic]->cd();
    gv2_sys[ic]->Draw("A5");
    gv2_sys[ic]->GetXaxis()->CenterTitle();
    gv2_sys[ic]->GetYaxis()->CenterTitle();
    gv2_sys[ic]->GetYaxis()->SetTitleOffset(1.5);
    gv2_sys[ic]->GetYaxis()->SetTitle("#it{v}_{2}^{#varUpsilon(1S)}");
    gv2_sys[ic]->GetXaxis()->SetTitle("p_{T}^{#varUpsilon(1S)} (GeV/c)");
    if(ic==0){gv2_sys[ic]->GetXaxis()->SetTitleSize(0.0448); gv2_sys[ic]->GetXaxis()->SetTitleOffset(1.35);}
    else if(ic!=0){gv2_sys[ic]->GetXaxis()->SetTitleSize(0.055); gv2_sys[ic]->GetXaxis()->SetTitleOffset(1.1);}
//    gv2_sys[ic]->GetXaxis()->SetTitleSize(0);
    gv2_sys[ic]->GetXaxis()->SetLabelSize(0);
    gv2_sys[ic]->GetXaxis()->SetTickLength(0);
    gv2_sys[ic]->GetXaxis()->SetRangeUser(0,xmax);
    gv2_sys[ic]->GetXaxis()->SetLimits(0,xmax);
    gv2_sys[ic]->SetMinimum(ymin);
    gv2_sys[ic]->SetMaximum(ymax);
    
    gPad->SetLeftMargin(0.19);
    gPad->SetBottomMargin(0.17);
    gPad->SetTopMargin(0.065);
    if(ic!=nCentBin-1)gPad->SetRightMargin(0);
    if(ic!=0){
      gv2_sys[ic]->GetYaxis()->SetTitle(" ");
      gv2_sys[ic]->GetYaxis()->SetTitleSize(0);
      globtex->SetTextSize(0.043);
      globtex->DrawLatex(0.04, sz_init, Form("#bf{Cent. %.f-%.f%s}",centBin[ic],centBin[ic+1],perc.Data()));

      gPad->SetLeftMargin(0);
    }
    
    if(ic==0){
      globtex->DrawLatex(0.23, sz_init-sz_step, "p_{T}^{#mu} > 3.5 GeV/c");
      globtex->DrawLatex(0.23, sz_init-sz_step*1.8, "|y| < 2.4");
      globtex->DrawLatex(0.23, sz_init, Form("#bf{Cent. %.f-%.f%s}",centBin[ic],centBin[ic+1],perc.Data()));
      leg->AddEntry(gv2[0],"#varUpsilon(1S)","pe");
      leg->Draw("same");
    }

    gv2[ic]->Draw("pe same");
    dashedLine(0.,0.,xmax,0.,1,1);
    
    solidLine(xmax/4,ymin,xmax/4,ymin+ymin_d,1,1);
    solidLine(xmax/2,ymin,xmax/2,ymin+ymin_d,1,1);
    solidLine(xmax/4*3,ymin,xmax/4*3,ymin+ymin_d,1,1);
    solidLine(xmax/4,ymax,xmax/4,ymax-ymin_d,1,1);
    solidLine(xmax/2,ymax,xmax/2,ymax-ymin_d,1,1);
    solidLine(xmax/4*3,ymax,xmax/4*3,ymax-ymin_d,1,1);

    lab_pos_diff = 0.196;
    if(ic!=0) globtex_label->SetTextSize(0.047); 
    if(ic==1) {lab_posx = 0.1; lab_pos_diff = 0.25;}
    else if(ic==2) {lab_posx = 0.1; lab_pos_diff = 0.24;}
    globtex_label->DrawLatex(lab_posx, lab_posy, "0-3");
    globtex_label->DrawLatex(lab_posx+lab_pos_diff,  lab_posy,"3-6");
    if(ic==2) lab_pos_diff = lab_pos_diff - 0.003;
    globtex_label->DrawLatex(lab_posx+lab_pos_diff*2, lab_posy, "6-10");
    
    lab_pos_diff = lab_pos_diff - 0.0003;
    if(ic==1) lab_pos_diff = lab_pos_diff - 0.005;
    globtex_label->DrawLatex(lab_posx+lab_pos_diff*3, lab_posy, "10-50");
    gv2_sys[ic]->SetName(Form("gr_sys_v2_vs_pt_Cent%.f%.f",centBin[ic],centBin[ic+1]));
    gv2[ic]->SetName(Form("gr_point_v2_vs_pt_Cent%.f%.f",centBin[ic],centBin[ic+1]));
    if(ic==0) CMS_lumi_square( pad1[ic], 0, iPos, 1 );
    if(ic==2) CMS_lumi_square( pad1[ic], iPeriod, iPos, 0 );
    pad1[ic]->Update();
    c1->Update();
  }
  c1->SaveAs("Figure4.pdf");

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
  gv2_int_sys->GetXaxis()->SetLabelSize(0);
  gv2_int_sys->GetYaxis()->SetNdivisions(506);
  gv2_int_sys->GetYaxis()->SetLabelSize(0.04);
  gv2_int_sys->GetXaxis()->SetTickLength(0);
  gv2_int_sys->GetXaxis()->SetRangeUser(0,xmax);
  gv2_int_sys->GetXaxis()->SetLimits(0,xmax);
  ymin = -0.05;
  ymax = 0.15;
  gv2_int_sys->SetMinimum(ymin);
  gv2_int_sys->SetMaximum(ymax);
  
  gv2_int->Draw("pe same");
  
  TLegend *leg_int= new TLegend(0.74, 0.70, 0.925, 0.77);
  SetLegendStyle(leg_int);
  leg_int->SetTextSize(0.044);
  leg_int->AddEntry(gv2_int,"#varUpsilon(1S)","pe");
  leg_int->Draw("same");
  dashedLine(0.,0.,xmax,0.,1,1);
  solidLine(xmax/4,ymin,xmax/4,ymin+ymin_d,1,1);
  solidLine(xmax/2,ymin,xmax/2,ymin+ymin_d,1,1);
  solidLine(xmax/4*3,ymin,xmax/4*3,ymin+ymin_d,1,1);
  solidLine(xmax/4,ymax,xmax/4,ymax-ymin_d,1,1);
  solidLine(xmax/2,ymax,xmax/2,ymax-ymin_d,1,1);
  solidLine(xmax/4*3,ymax,xmax/4*3,ymax-ymin_d,1,1);
  
  globtex_int->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
  globtex_int->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
  globtex_int->DrawLatex(0.23, sz_init-sz_step*2, "Cent. 10-90%");
  
  double lab_posx_ = 0.269; double lab_posy_ = 0.14; double lab_pos_diff_ = 0.186;
  TLatex* globtex_label_ = new TLatex();
  globtex_label_->SetNDC();
  globtex_label_->SetTextAlign(12); //left-center
  globtex_label_->SetTextFont(42);
  globtex_label_->SetTextSize(0.040);
  globtex_label_->DrawLatex(lab_posx_, lab_posy_, "0-3");
  globtex_label_->DrawLatex(lab_posx_+lab_pos_diff_,  lab_posy_,"3-6");
  globtex_label_->DrawLatex(lab_posx_+lab_pos_diff_*2, lab_posy_, "6-10");
  lab_pos_diff_ =lab_pos_diff_ - 0.0002;
  globtex_label_->DrawLatex(lab_posx_+lab_pos_diff_*3, lab_posy_, "10-50");
  
  CMS_lumi_square( c2, iPeriod, iPos,1 );
	c2->Update();
  c2->SaveAs("Figure2right.pdf");



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
