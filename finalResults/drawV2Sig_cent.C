#include <iostream>
#include "../Style_jaebeom.h"
#include "../commonUtility.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

using namespace std;

void drawV2Sig_cent(){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  
  const int nBin=3;
  const int nCentBin = 4;
  double exsys[nBin] =  {1.5, 1.5,22};

  double xmax = 8; double ymin = -0.05; double ymax = 0.15;
  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile *fSpectra = new TFile("Spectra_v2.root");
  double ptBin[nBin+1] = {0, 3, 6, 50};
  double centBin[nCentBin+1] = {0,10,30,50,90};

	TGraphErrors* gv2[nCentBin];
	TGraphErrors* gv2_sys[nCentBin];
	TGraphErrors* gv2_int = (TGraphErrors*) fSpectra->Get("gr_Cent_PtInt");
	TGraphErrors* gv2_int_sys = (TGraphErrors*) fSpectra->Get("gr_sys_Cent_PtInt");

	TGraphErrors* gr_v2[nCentBin];
	TGraphErrors* gr_v2_sys[nCentBin];
	TGraphErrors* gr_v2_int = new TGraphErrors();
	TGraphErrors* gr_v2_int_sys = new TGraphErrors();
  
  double pxtmp, pytmp, extmp, eytmp, esy;
  double esx = xmax/(3*8);
  for(int ipt = 0; ipt<nBin; ipt++){
    gv2[ipt] = (TGraphErrors*) fSpectra->Get(Form("gr_Cent_Pt%.f%.f",ptBin[ipt],ptBin[ipt+1])); 
    gv2_sys[ipt] = (TGraphErrors*) fSpectra->Get(Form("gr_sys_Cent_Pt%.f%.f",ptBin[ipt],ptBin[ipt+1])); 
    for(int ip = 0; ip<gv2[ipt]->GetN(); ip++){
      pxtmp = 0; pytmp = 0; extmp = 0; eytmp = 0; esy=0;
      gv2[ipt]->GetPoint(ip,pxtmp,pytmp);
      eytmp = gv2[ipt]->GetErrorY(ip);
      esy = gv2_sys[ipt]->GetErrorY(ip);
      gv2[ipt]->SetPoint(ip, xmax/(nCentBin*2)+xmax/nCentBin*ip, pytmp);
      gv2[ipt]->SetPointError(ip, 0, eytmp);
      gv2_sys[ipt]->SetPoint(ip, xmax/(nCentBin*2)+xmax/nCentBin*ip, pytmp);
      gv2_sys[ipt]->SetPointError(ip, esx, esy);
    }
    SetGraphStyle(gv2[ipt],1,1);
    SetGraphStyleSys(gv2_sys[ipt],1);
  }
  
  //Set for Integrated Bin
  for(int ip = 0; ip<gv2_int->GetN(); ip++){
    pxtmp = 0; pytmp = 0; extmp = 0; eytmp = 0; esy=0;
    gv2_int->GetPoint(ip,pxtmp,pytmp);
    eytmp = gv2_int->GetErrorY(ip);
    esy = gv2_int_sys->GetErrorY(ip);
    gv2_int->SetPoint(ip, xmax/(nCentBin*2)+xmax/nCentBin*ip, pytmp);
    gv2_int->SetPointError(ip, 0, eytmp);
    gv2_int_sys->SetPoint(ip, xmax/(nCentBin*2)+xmax/nCentBin*ip, pytmp);
    gv2_int_sys->SetPointError(ip, esx, esy);
  }
  SetGraphStyle(gv2_int,1,1);
  SetGraphStyleSys(gv2_int_sys,1);


  TCanvas* c1[nCentBin];
  TLegend *leg[nCentBin];
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
  for(int ic=0; ic<nBin; ic++){
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
    gv2_sys[ic]->GetXaxis()->SetTitle("Centrality (%)");
    gv2_sys[ic]->GetXaxis()->SetTitleSize(0.055);
//    gv2_sys[ic]->GetXaxis()->SetTitleSize(0);
    gv2_sys[ic]->GetXaxis()->SetLabelSize(0);
    gv2_sys[ic]->GetXaxis()->SetTickLength(0);
    gv2_sys[ic]->GetXaxis()->SetRangeUser(0,xmax);
    gv2_sys[ic]->GetXaxis()->SetLimits(0,xmax);
    gv2_sys[ic]->SetMinimum(ymin);
    gv2_sys[ic]->SetMaximum(ymax);

    gv2[ic]->Draw("pe same");
    dashedLine(0.,0.,xmax,0.,1,1);

    TString perc = "%";
    leg[ic]= new TLegend(0.74, 0.70, 0.925, 0.77);
    SetLegendStyle(leg[ic]);
    leg[ic]->SetTextSize(0.044);
    leg[ic]->AddEntry(gv2[ic],"#varUpsilon(1S)","pe");
    leg[ic]->Draw("same");
    solidLine(xmax/3,ymin,xmax/3,ymin+ymin_d,1,1);
    solidLine(xmax/3*2,ymin,xmax/3*2,ymin+ymin_d,1,1);
    
    solidLine(xmax/3,ymax,xmax/3,ymax-ymin_d,1,1);
    solidLine(xmax/3*2,ymax,xmax/3*2,ymax-ymin_d,1,1);


    globtex->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
    globtex->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
    globtex->DrawLatex(0.23, sz_init-sz_step*2, Form("p_{T}^{#varUpsilon}. %.f-%.f%s",centBin[ic],centBin[ic+1],perc.Data()));
    globtex_label->DrawLatex(lab_posx, lab_posy, "0-3");
    globtex_label->DrawLatex(lab_posx+lab_pos_diff,  lab_posy,"3-6");
    globtex_label->DrawLatex(lab_posx+lab_pos_diff*2, lab_posy, "6-50");
    gv2_sys[ic]->SetName(Form("gr_sys_v2_vs_pt_Cent%.f%.f",centBin[ic],centBin[ic+1]));
    gv2[ic]->SetName(Form("gr_point_v2_vs_pt_Cent%.f%.f",centBin[ic],centBin[ic+1]));
    CMS_lumi_square( c1[ic], iPeriod, iPos );
    c1[ic]->Update();
    c1[ic]->SaveAs(Form("v2Sig_Cent_pt%d.pdf",ic));
  }

  TLatex* globtex_int = new TLatex();
  globtex_int->SetNDC();
  globtex_int->SetTextAlign(12); //left-center
  globtex_int->SetTextFont(42);
  globtex_int->SetTextSize(0.040);




  //Integrated Panel
  TGraphErrors* gv2_int1S = (TGraphErrors*) fSpectra->Get("gr_Int1S");
  TGraphErrors* gv2_int2S = (TGraphErrors*) fSpectra->Get("gr_Int2S");
  TGraphErrors* gv2_int1S_sys = (TGraphErrors*) fSpectra->Get("gr_sys_Int1S");
  TGraphErrors* gv2_int2S_sys = (TGraphErrors*) fSpectra->Get("gr_sys_Int2S");
  
  double xmax_int = 1;
  double pos_d = xmax_int/8;
  double pxint=0; double pyint =0; double eyint = 0; double esysint = 0;
  double xlonger = 120;
    
  gv2_int1S->GetPoint(0, pxint, pyint);
  eyint = gv2_int1S->GetErrorY(0);
  esysint = gv2_int1S_sys->GetErrorY(0);
  gv2_int1S->SetPoint(0,xmax_int/2-pos_d,pyint);
  gv2_int1S->SetPointError(0,0,eyint);
  gv2_int1S_sys->SetPoint(0,xmax_int/2-pos_d,pyint);
  gv2_int1S_sys->SetPointError(0,xmax_int/12,esysint);

  pxint=0; pyint =0; eyint = 0; esysint = 0;
  gv2_int2S->GetPoint(0, pxint, pyint);
  eyint = gv2_int2S->GetErrorY(0);
  esysint = gv2_int2S_sys->GetErrorY(0);
  gv2_int2S->SetPoint(0,xmax_int/2+pos_d,pyint);
  gv2_int2S->SetPointError(0,0,eyint);
  gv2_int2S_sys->SetPoint(0,xmax_int/2+pos_d,pyint);
  gv2_int2S_sys->SetPointError(0,xmax_int/12,esysint);

  SetGraphStyle(gv2_int1S,1,1);
  SetGraphStyleSys(gv2_int1S_sys,1);
  SetGraphStyle(gv2_int2S,2,2);
  SetGraphStyleSys(gv2_int2S_sys,2);
  gv2_int1S_sys->GetXaxis()->SetLimits(0,xmax_int);
  gv2_int1S_sys->GetXaxis()->SetRangeUser(0,xmax_int);
  gv2_int1S_sys->SetMinimum(ymin);
  gv2_int1S_sys->SetMaximum(ymax);
  gv2_int1S_sys->GetXaxis()->SetNdivisions(101);
  gv2_int1S_sys->GetXaxis()->SetLabelSize(0);
  gv2_int1S_sys->GetYaxis()->SetTickLength(0.03*600/xlonger);
  gv2_int1S_sys->GetYaxis()->SetLabelSize(0);
  
  TCanvas *c2 = new TCanvas("cint","", 600,600);
  TPad* pad_diff = new TPad("pad_diff", "",0, 0, 600/(600.+xlonger), 1.0); // vs centrality
  pad_diff->SetRightMargin(0);
  pad_diff->SetLeftMargin(0.20);
  pad_diff->SetBottomMargin(0.13);
  pad_diff->SetTopMargin(0.067);
  TPad* pad_int = new TPad("pad_int", "",600/(600.+xlonger), 0, 1.0, 1.0); // centrality-integrated
  pad_int->SetBottomMargin(0.13);
  pad_int->SetTopMargin(0.067);
  pad_int->SetLeftMargin(0);
  pad_int->SetRightMargin(0.032*600/xlonger);

  c2->cd();
  pad_diff->Draw();
  pad_diff->cd();


  gv2_int->SetMarkerStyle(kFullCircle);
  gv2_int_sys->Draw("A5");
  gv2_int_sys->GetXaxis()->CenterTitle();
  gv2_int_sys->GetYaxis()->CenterTitle();
  gv2_int_sys->GetYaxis()->SetTitleOffset(1.5);
  gv2_int_sys->GetXaxis()->SetTitleOffset(1.1);
  gv2_int_sys->GetYaxis()->SetTitle("#it{v}_{2}^{#varUpsilon(1S)}");
  gv2_int_sys->GetXaxis()->SetTitle("Centrality (%)");
  gv2_int_sys->GetXaxis()->SetTitleSize(0.055);
  gv2_int_sys->GetXaxis()->SetLabelSize(0);
  gv2_int_sys->GetXaxis()->SetTickLength(0);
  gv2_int_sys->GetXaxis()->SetRangeUser(0,xmax);
  gv2_int_sys->GetXaxis()->SetLimits(0,xmax);
  gv2_int_sys->SetMinimum(ymin);
  gv2_int_sys->SetMaximum(ymax);
  
  gv2_int->Draw("pe same");
  
  TLegend *leg_int= new TLegend(0.74, 0.70, 0.925, 0.77);
  SetLegendStyle(leg_int);
  leg_int->SetTextSize(0.044);
  leg_int->AddEntry(gv2_int,"#varUpsilon(1S)","pe");
  leg_int->AddEntry(gv2_int2S,"#varUpsilon(2S)","pe");
  leg_int->Draw("same");
  dashedLine(0.,0.,xmax,0.,1,1);
  ymin_d = 0.007;
  solidLine(xmax/4,ymin,xmax/4,ymin+ymin_d,1,1);
  solidLine(xmax/2,ymin,xmax/2,ymin+ymin_d,1,1);
  solidLine(xmax/4*3,ymin,xmax/4*3,ymin+ymin_d,1,1);
  solidLine(xmax/4,ymax,xmax/4,ymax-ymin_d,1,1);
  solidLine(xmax/2,ymax,xmax/2,ymax-ymin_d,1,1);
  solidLine(xmax/4*3,ymax,xmax/4*3,ymax-ymin_d,1,1);
  
  globtex_int->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
  globtex_int->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
  globtex_int->DrawLatex(0.23, sz_init-sz_step*2, "p_{T}^{#varUpsilon} < 50 GeV/c");
  
  double lab_posx_ = 0.269; double lab_posy_ = 0.103; double lab_pos_diff_ = 0.186;
  TLatex* globtex_label_ = new TLatex();
  globtex_label_->SetNDC();
  globtex_label_->SetTextAlign(12); //left-center
  globtex_label_->SetTextFont(42);
  globtex_label_->SetTextSize(0.040);
  globtex_label_->DrawLatex(lab_posx_, lab_posy_, "0-10");
  globtex_label_->DrawLatex(lab_posx_+lab_pos_diff_,  lab_posy_,"10-30");
  lab_pos_diff_ = lab_pos_diff_ + 0.0015;
  globtex_label_->DrawLatex(lab_posx+lab_pos_diff_*2, lab_posy_, "30-50");
  lab_pos_diff_ = lab_pos_diff_ + 0.0018;
  globtex_label_->DrawLatex(lab_posx+lab_pos_diff_*3, lab_posy_, "50-90");
     
  pad_diff->Update();
  CMS_lumi_square( pad_diff, iPeriod, iPos );
  
  pad_diff->Update();

 
  c2->cd(); 
  pad_int->Draw();
  pad_int->cd();
  gv2_int1S_sys->Draw("A5");
  gv2_int2S_sys->Draw("5");
  gv2_int1S->Draw("P");
  gv2_int2S->Draw("P");
  dashedLine(0.,0.,xmax_int,0.,1,1);
  globtex_label_->SetTextSize(0.21);
  globtex_label_->DrawLatex(lab_posx-0.2+lab_pos_diff_, lab_posy_, "0-90");
  
	c2->Update();
  c2->SaveAs("v2Sig_cent_w2S.pdf");



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
