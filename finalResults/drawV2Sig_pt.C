#include <iostream>
#include "../Style_jaebeom.h"
#include "../commonUtility.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

using namespace std;

void drawV2Sig_pt(int fVer=0){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  
  const int nBin=3;
  double exsys[nBin] =  {1.5, 1.5,22};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nBin];
  TFile* fIn_int[nBin];
  double ptBin[nBin+1] = {0, 3, 6, 50};

	TGraphErrors* gv2[nBin];
	TGraphErrors* gv2_int[nBin];
  for (int is=0; is<nBin; is++){
  	if(fVer==0) fIn[is] = new TFile(Form("../SimFit/SigAllFreeFitFix/FitResult/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality60-100_m8-14_OS.root",ptBin[is],ptBin[is+1]),"READ");
  	if(fVer==1) fIn[is] = new TFile(Form("../SimFit/AllParmFree/FitResult/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality60-100_m8-14_OS.root",ptBin[is],ptBin[is+1]),"READ");
  	if(fVer==0) fIn_int[is] = new TFile(Form("../SimFit/SigAllFreeFitFix/FitResult/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality0-200_m8-14_OS.root",ptBin[is],ptBin[is+1]),"READ");
  	if(fVer==1) fIn_int[is] = new TFile(Form("../SimFit/AllParmFree/FitResult/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality0-200_m8-14_OS.root",ptBin[is],ptBin[is+1]),"READ");
    gv2[is]=(TGraphErrors*)fIn[is]->Get("v2vspt");
    gv2_int[is]=(TGraphErrors*)fIn_int[is]->Get("v2vspt");
  }

  TH1D *final_v2 = new TH1D("Final_v2",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,ptBin);
  TH1D *final_v2_int = new TH1D("Final_v2_int",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,ptBin);
  

  //// get bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;

  for (int is=0; is<nBin; is++){
    cout << is+1 <<"th bin***************" << endl;
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
    gv2[is] -> GetPoint(0,pxtmp,pytmp);
    eytmp = gv2[is]->GetErrorY(0);
    final_v2->SetBinContent(is+1, pytmp);
    final_v2->SetBinError(is+1, eytmp);
  }

  final_v2->GetXaxis()->CenterTitle();
  final_v2->GetYaxis()->CenterTitle();
  final_v2->GetYaxis()->SetTitleOffset(1.2);
  final_v2->GetXaxis()->SetTitleOffset(1.);
  final_v2->SetMinimum(-0.05);
  final_v2->SetMaximum(0.2);

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

  TGraphErrors* final_v2_g = new TGraphErrors(final_v2);
  TFile* fSys = new TFile("../Systematic/merged_sys.root","read");
  TH1D* hsys = (TH1D*) fSys->Get("hpt_Cent3050_merged");
  TGraphErrors *gsys = new TGraphErrors();
  for(int ib=0;ib<final_v2->GetNbinsX();ib++){
    gsys->SetPoint(ib, final_v2->GetBinCenter(ib+1), final_v2->GetBinContent(ib+1));
    gsys->SetPointError(ib, final_v2->GetBinWidth(ib+1)/2, final_v2->GetBinContent(ib+1)*hsys->GetBinContent(ib+1));
    final_v2_g->SetPointError(ib,0,final_v2->GetBinError(ib+1));
  }

  TGraphErrors* final_v2_int_g = new TGraphErrors(final_v2_int);
  TH1D* hsys_int = (TH1D*) fSys->Get("hpt_Cent0100_merged");
  TGraphErrors *gsys_int = new TGraphErrors();
  for(int ib=0;ib<final_v2_int->GetNbinsX();ib++){
    gsys_int->SetPoint(ib, final_v2_int->GetBinCenter(ib+1), final_v2_int->GetBinContent(ib+1));
    gsys_int->SetPointError(ib, final_v2_int->GetBinWidth(ib+1)/2, final_v2_int->GetBinContent(ib+1)*hsys_int->GetBinContent(ib+1));
    final_v2_int_g->SetPointError(ib,0,final_v2_int->GetBinError(ib+1));
  }
  SetGraphStyle(final_v2_g,4,4);
  SetGraphStyleSys(gsys,4);

  SetGraphStyle(final_v2_int_g,4,4);
  SetGraphStyleSys(gsys_int,4);

  double xmax = 50; 
 
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.15);
  c1->cd();
  final_v2_g->SetMarkerStyle(kFullCircle);
  gsys->Draw("A5");
  gsys->GetXaxis()->CenterTitle();
  gsys->GetYaxis()->CenterTitle();
  gsys->GetYaxis()->SetTitleOffset(1.2);
  gsys->GetXaxis()->SetTitleOffset(1.);
  gsys->GetYaxis()->SetTitle("v_{2}^{sig}");
  gsys->GetXaxis()->SetTitle("p_{T}^{#Upsilon(1S)} (GeV/c)");
  gsys->GetXaxis()->SetRangeUser(0,xmax);
  gsys->SetMinimum(-0.05);
  gsys->SetMaximum(0.2);

  final_v2_g->Draw("pe same");

  TLegend *leg= new TLegend(0.57, 0.62, 0.785, 0.74);
  SetLegendStyle(leg);
  leg->AddEntry(final_v2,"#Upsilon(1S)","pe");
  leg->Draw("same");
  dashedLine(0.,0.,xmax,0.,1,1);
  
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  double sz_init = 0.875; double sz_step = 0.0535;
  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
  globtex->DrawLatex(0.22, sz_init-sz_step, "|y| < 2.4");
  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 30-50%");
  
  CMS_lumi_square( c1, iPeriod, iPos );
	c1->Update();
  c1->SaveAs(Form("v2Sig_pt_%d.pdf",fVer));

 
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.15);
  c2->cd();
  final_v2_int_g->SetMarkerStyle(kFullCircle);
  gsys_int->Draw("A5");
  gsys_int->GetXaxis()->CenterTitle();
  gsys_int->GetYaxis()->CenterTitle();
  gsys_int->GetYaxis()->SetTitleOffset(1.2);
  gsys_int->GetXaxis()->SetTitleOffset(1.);
  gsys_int->GetYaxis()->SetTitle("v_{2}^{sig}");
  gsys_int->GetXaxis()->SetTitle("p_{T}^{#Upsilon(1S)} (GeV/c)");
  gsys_int->GetXaxis()->SetRangeUser(0,xmax);
  gsys_int->SetMinimum(-0.05);
  gsys_int->SetMaximum(0.2);

  final_v2_int_g->Draw("pe same");

  TLegend *leg_int= new TLegend(0.57, 0.62, 0.785, 0.74);
  SetLegendStyle(leg_int);
  leg_int->AddEntry(final_v2_int,"#Upsilon(1S)","pe");
  leg_int->Draw("same");
  dashedLine(0.,0.,xmax,0.,1,1);
  
  TLatex* globtex_int = new TLatex();
  globtex_int->SetNDC();
  globtex_int->SetTextAlign(12); //left-center
  globtex_int->SetTextFont(42);
  globtex_int->SetTextSize(0.040);
  globtex_int->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
  globtex_int->DrawLatex(0.22, sz_init-sz_step, "|y| < 2.4");
  globtex_int->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");
  
  CMS_lumi_square( c2, iPeriod, iPos );
	c2->Update();
  c2->SaveAs(Form("v2Sig_pt_int_%d.pdf",fVer));

	return;
}
