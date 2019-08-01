#include <iostream>
#include "../Style_jaebeom.h"
#include "../commonUtility.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

using namespace std;

void draw_Figure3_pt2(){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  
  const int nBin=3;
  const int nCentBin = 4;
  double exsys[nBin] =  {1.5, 1.5,22};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn;
  double ptBin[nBin+1] = {0, 3, 6, 50};
  double centBin[nCentBin+1] = {0,10,30,50,90};

	TGraphErrors* gv2;
	TGraphErrors* gsys;
  
  TH1D *final_v2;
  
  fIn = new TFile("./out_v2_vs_pt.root","READ");

  final_v2 = new TH1D("Final_v2",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,0,3);
  gv2 = (TGraphErrors*)fIn->Get("gr_point_v2_vs_pt_Cent1030");

  //// get bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;

  for (int is=0; is<nBin; is++){
	  cout << is+1 <<"th bin***************" << endl;
	  pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
	  gv2 -> GetPoint(is,pxtmp,pytmp);
	  eytmp = gv2->GetErrorY(is);
	  final_v2->SetBinContent(is+1, pytmp);
	  final_v2->SetBinError(is+1, eytmp);
  }
  final_v2->GetXaxis()->CenterTitle();
  final_v2->GetYaxis()->CenterTitle();
  final_v2->GetYaxis()->SetTitleOffset(1.2);
  final_v2->GetXaxis()->SetTitleOffset(1.);
  final_v2->GetXaxis()->SetLabelSize(0);
  final_v2->GetXaxis()->SetNdivisions(210);
  final_v2->SetMinimum(-0.05);
  final_v2->SetMaximum(0.2);

  double sysX, sysY, sysX_Error, sysY_Error;
  TGraphErrors* final_v2_g = new TGraphErrors(final_v2);
  gsys = (TGraphErrors*)fIn->Get("gr_sys_v2_vs_pt_Cent1030");
  TGraphErrors* gsys2 = new TGraphErrors();
  TH1D *hsys= new TH1D("hsys",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,0,3);

  for(int ib=0;ib<final_v2->GetNbinsX();ib++){
	  sysX=0; sysY=0; sysX_Error=0; sysY_Error=0;
	  gsys->GetPoint(ib,sysX,sysY);
	  sysY_Error = gsys->GetErrorY(ib);
	  gsys2->SetPoint(ib,final_v2->GetBinCenter(ib+1),sysY);
	  gsys2->SetPointError(ib,final_v2->GetBinWidth(ib+1)/2,sysY_Error);
	  hsys->SetBinContent(ib+1, final_v2->GetBinContent(ib+1));
	  hsys->SetBinError(ib+1, gsys2->GetErrorY(ib+1));
	  final_v2_g->SetPointError(ib,0,final_v2->GetBinError(ib+1));
  }
  SetGraphStyle(final_v2_g,1,1);
  SetGraphStyleSys(gsys2,1);
  


  double xmax = 3; 
 
  TCanvas* c1;
  TLegend *leg;
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  double sz_init = 0.867; double sz_step = 0.0535;
  c1= new TCanvas("c1","c1",600,600);
  gPad->SetLeftMargin(0.19);
  gPad->SetBottomMargin(0.17);
  gPad->SetTopMargin(0.065);
  c1->cd();
  final_v2_g->SetMarkerStyle(kFullCircle);
  gsys2->Draw("A5");
  gsys2->GetXaxis()->CenterTitle();
  gsys2->GetYaxis()->CenterTitle();
  gsys2->GetYaxis()->SetTitleOffset(1.5);
  gsys2->GetXaxis()->SetTitleOffset(1.1);
  gsys2->GetYaxis()->SetTitle("v_{2}^{#Upsilon(1S)}");
  gsys2->GetXaxis()->SetTitle("p_{T}^{#Upsilon(1S)} (GeV/c)");
  gsys2->GetXaxis()->SetRangeUser(0,xmax);
  gsys2->SetMinimum(-0.05);
  gsys2->SetMaximum(0.15);
  gsys2->GetXaxis()->SetLabelSize(0);
  gsys2->GetXaxis()->SetNdivisions(210);
  gsys2->GetXaxis()->SetTickLength(0);

  final_v2_g->Draw("pe same");

  TF1 *line = new TF1("line","0",0,3);
  line->SetLineWidth(1);
  line->SetLineStyle(7);
  line->SetLineColor(kBlack);

  line->Draw("same");

  TLine *linebin1 = new TLine(1,-0.05,1,-0.045);
  TLine *linebin2 = new TLine(2,-0.05,2,-0.045);
  TLine *linebin3 = new TLine(3,-0.05,3,-0.045);
  linebin1->Draw("same");
  linebin2->Draw("same");
  //linebin3->Draw("same");


  TString perc = "%";
  leg= new TLegend(0.74, 0.70, 0.925, 0.77);
  SetLegendStyle(leg);
  leg->SetTextSize(0.044);
  leg->AddEntry(final_v2_g,"#Upsilon(1S)","pe");
  leg->Draw("same");
  //globtex->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
  //globtex->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
  //globtex->DrawLatex(0.23, sz_init-sz_step*2, Form("Cent. 0-10%s",perc.Data()));
  globtex->DrawLatex(0.23, sz_init, "|y| < 2.4");
  globtex->DrawLatex(0.23, sz_init-sz_step, Form("Cent. 10-30%s",perc.Data()));
  globtex->DrawLatex(0.30, 0.14, "0-3");
  globtex->DrawLatex(0.55, 0.14, "3-6");
  globtex->DrawLatex(0.80, 0.14, "6-50");
  CMS_lumi_square( c1, iPeriod, iPos );
  c1->Update();
  c1->SaveAs("v2Sig_pt_Cent_2.pdf");


  return;
} 
