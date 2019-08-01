#include <iostream>
#include "../Style_jaebeom.h"
#include "../commonUtility.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

using namespace std;

void draw_Figure2_cent(){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  
  const int nBin=4;
  double exsys[nBin] =  {1.5, 1.5,22};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn;
  TFile* fInt;
  double centBin[nBin+1] = {0, 10, 30, 50, 90};
  double centBin_int[nBin+1] = {0,90};

  TGraphErrors* gv2;
  TGraphErrors* gsys;
  TGraphErrors* gv2_int;
  TGraphErrors* gsys_int;
  TGraphErrors* gv2Int1S;
  TGraphErrors* gv2Int2S;
  TGraphErrors* gv2Int1Ssys;
  TGraphErrors* gv2Int2Ssys;

  fIn = new TFile("./out_v2_vs_cent.root","READ");
  fInt = new TFile("./out_v2_Int.root","READ");
  //fIn_int= new TFile("./out_v2_vs_cent.root","READ");
  gv2 = (TGraphErrors*)fIn->Get("gr_point_v2_vs_cent_pt050");

  // Integrated
  gv2Int1S = (TGraphErrors*)fInt->Get("gr_point_v2_Int_1S");
  gv2Int2S = (TGraphErrors*)fInt->Get("gr_point_v2_Int_2S");
  gv2Int1Ssys = (TGraphErrors*)fInt->Get("gr_sys_v2_Int_1S");
  gv2Int2Ssys = (TGraphErrors*)fInt->Get("gr_sys_v2_Int_2S");
  //	gv2_int = (TGraphErrors*)fIn->Get("gr_point_v2_vs_cent_pt050");

  TH1D *final_v2 = new TH1D("Final_v2",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,0,4);
  //TH1D *final_v2 = new TH1D("Final_v2",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,centBin);
  TH1D *final_v2_int = new TH1D("Final_v2_int",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,centBin);


  //// get bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;

  for (int is=0; is<nBin; is++){
    cout << is+1 <<"th bin***************" << endl;
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
    gv2-> GetPoint(is,pxtmp,pytmp);
    eytmp= gv2->GetErrorY(is);
    final_v2->SetBinContent(is+1, pytmp);
    final_v2->SetBinError(is+1, eytmp);
  }

  final_v2->GetXaxis()->CenterTitle();
  final_v2->GetYaxis()->CenterTitle();
  final_v2->GetYaxis()->SetTitleOffset(1.2);
  final_v2->GetXaxis()->SetTitleOffset(1.);
  final_v2->SetMinimum(-0.05);
  final_v2->SetMaximum(0.15);
  
  final_v2_int->SetBinContent(1,0.065);
  final_v2_int->SetBinError(1,0.09);

  final_v2_int->GetXaxis()->CenterTitle();
  final_v2_int->GetYaxis()->CenterTitle();
  final_v2_int->GetYaxis()->SetTitleOffset(1.2);
  final_v2_int->GetXaxis()->SetTitleOffset(1.);
  final_v2_int->SetMinimum(-0.05);
  final_v2_int->SetMaximum(0.15);

  double sysX, sysY, sysX_Error, sysY_Error;
  TGraphErrors* final_v2_g = new TGraphErrors(final_v2);
  gsys = (TGraphErrors*)fIn->Get("gr_sys_v2_vs_cent_pt050");
  TGraphErrors* gsys2 = new TGraphErrors();
  TH1D *hsys= new TH1D("hsys",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,0,4);
  //TH1D *hsys= new TH1D("hsys",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,centBin);
  for(int ib=0;ib<final_v2->GetNbinsX();ib++){
	sysX=0; sysY=0; sysX_Error=0; sysY_Error=0;
	gsys-> GetPoint(ib,sysX,sysY);
	sysY_Error = gsys->GetErrorY(ib);
	hsys->SetBinContent(ib+1, sysY);
	hsys->SetBinError(ib+1, sysY_Error);
    gsys2->SetPoint(ib,final_v2->GetBinCenter(ib+1),sysY);
    gsys2->SetPointError(ib,final_v2->GetBinWidth(ib+1)/2,sysY_Error);
    final_v2_g->SetPointError(ib,0,final_v2->GetBinError(ib+1));
    //final_v2_g->SetPointError(ib,0,final_v2->GetBinError(ib+1));
	cout << "Syst Error : " << hsys->GetBinContent(ib+1) << endl;
	cout << "Nomi Error : " << hsys->GetBinError(ib+1) << endl;
  }
  SetGraphStyle(final_v2_g,1,1);
  SetGraphStyleSys(gsys,1);
  SetGraphStyleSys(gsys2,1);

  SetGraphStyle(gv2Int1S,1,1);
  SetGraphStyleSys(gv2Int1Ssys,2);
  SetGraphStyle(gv2Int2S,2,2);
  SetGraphStyleSys(gv2Int1Ssys,1);
 

  double xmax = 4.05; // 90, 4 
  double xlonger = 120;
  
  // change x position
  double v2s[2] = {0.0, 0.0};
  double v2Err[2] = {0.0, 0.0};
  double v2Sys[2] = {0.0, 0.0};
  double tmp = 0.0;
  gv2Int1S->GetPoint(0, tmp, v2s[0]);
  gv2Int2S->GetPoint(0, tmp, v2s[1]);
  v2Err[0] = gv2Int1S->GetErrorY(0);
  v2Err[1] = gv2Int2S->GetErrorY(0);
  v2Sys[0] = gv2Int1Ssys->GetErrorY(0);
  v2Sys[1] = gv2Int2Ssys->GetErrorY(0);

  double tmpv2s1 = v2s[0];
  double tmpv2s2 = v2s[1];
  double tmpv2sEr1 = v2Err[0];
  double tmpv2sEr2 = v2Err[1];


  TH1F *hPadInt = new TH1F("hPadInt",";;",1,0,20);
  hPadInt->SetMaximum(0.15);
  hPadInt->SetMinimum(-0.05);

  hPadInt->GetXaxis()->SetNdivisions(101);
  hPadInt->GetXaxis()->SetLabelSize(0);
  hPadInt->GetXaxis()->SetRangeUser(0,20);
  hPadInt->GetYaxis()->SetTickLength(0.03*600/xlonger);
  hPadInt->GetYaxis()->SetLabelSize(0);
  hPadInt->GetXaxis()->SetTitle("");

  double xwth = 0.1;
  double swth = 2.0;
  double intbin1[] = {5.0-xwth, 5.0+xwth};
  double intbin2[] = {15.0-xwth, 15.0+xwth};
  double sysbin1[] = {5.0-swth, 5.0+swth};
  double sysbin2[] = {15.0-swth, 15.0+swth};

  TH1F *hv2Int1S2 = new TH1F("hv2Int1S2",";;",1,intbin1);
  TH1F *hv2Int2S2 = new TH1F("hv2Int2S2",";;",1,intbin2);
  TH1F *hv2Int1Ssys2 = new TH1F("hv2Int1Ssys2",";;",1,sysbin1);
  TH1F *hv2Int2Ssys2 = new TH1F("hv2Int2Ssys2",";;",1,sysbin2);
  hv2Int1S2->SetBinContent(1, v2s[0]);
  hv2Int1S2->SetBinError(1, v2Err[0]);
  hv2Int1Ssys2->SetBinContent(1, v2s[0]);
  hv2Int1Ssys2->SetBinError(1, v2Sys[0]);
  hv2Int2S2->SetBinContent(1, v2s[1]);
  hv2Int2S2->SetBinError(1, v2Err[1]);
  hv2Int2Ssys2->SetBinContent(1, v2s[1]);
  hv2Int2Ssys2->SetBinError(1, v2Sys[1]);

  TGraphErrors* gv2Int1S2 = new TGraphErrors(hv2Int1S2);
  TGraphErrors* gv2Int2S2 = new TGraphErrors(hv2Int2S2);
  TGraphErrors* gv2Int1Ssys2 = new TGraphErrors(hv2Int1Ssys2);
  TGraphErrors* gv2Int2Ssys2 = new TGraphErrors(hv2Int2Ssys2);


  SetGraphStyle(gv2Int1S2,1,1);
  SetGraphStyleSys(gv2Int1Ssys2,1);
  SetGraphStyle(gv2Int2S2,2,2);
  SetGraphStyleSys(gv2Int2Ssys2,2);
  gv2Int1S2->SetMarkerStyle(20);


  TF1 *line = new TF1("line","0",0,1000);
  line->SetLineWidth(1);
  line->SetLineStyle(7);
  line->SetLineColor(kBlack);
 
  TCanvas* c1 = new TCanvas("c1","c1",600+xlonger,600);
  TPad* pad_diff = new TPad("pad_diff","",0,0,600/(600.+xlonger),1.0);
  pad_diff->SetRightMargin(0);
  pad_diff->SetBottomMargin(0.17);
  pad_diff->SetTopMargin(0.065);
  TPad* pad_int = new TPad("pad_int","",600/(600.+xlonger),0,1.0,1.0);
  pad_int->SetBottomMargin(0.17);
  pad_int->SetTopMargin(0.065);
  pad_int->SetLeftMargin(0);
  pad_int->SetRightMargin(0.032*600/xlonger);

  
  gsys2->GetXaxis()->CenterTitle();
  gsys2->GetYaxis()->CenterTitle();
  gsys2->GetYaxis()->SetTitleOffset(1.2);
  gsys2->GetXaxis()->SetTitleOffset(1.);
  gsys2->GetYaxis()->SetTitle("v_{2}^{sig}");
  gsys2->GetXaxis()->SetTitle("Centrality (%)");
  gsys2->GetXaxis()->SetRangeUser(0,xmax);
  gsys2->SetMinimum(-0.05);
  gsys2->SetMaximum(0.15);
  gsys2->GetXaxis()->SetLabelSize(0);
  gsys2->GetXaxis()->SetNdivisions(210); // 210
  gsys2->GetXaxis()->SetNdivisions(210); // 210
  gsys2->GetXaxis()->SetTickLength(0);


  c1->cd();
  pad_diff->Draw();
  pad_diff->cd();
 
  final_v2_g->SetMarkerStyle(kFullCircle);
  gsys2->Draw("A5");
  line->Draw("same");
  
  gsys->GetXaxis()->CenterTitle();
  gsys->GetYaxis()->CenterTitle();
  gsys->GetYaxis()->SetTitleOffset(1.2);
  gsys->GetXaxis()->SetTitleOffset(1.);
  gsys->GetYaxis()->SetTitle("v_{2}^{sig}");
  gsys->GetXaxis()->SetTitle("Centrality (%)");
  gsys->GetXaxis()->SetRangeUser(0,xmax);
  gsys->SetMinimum(-0.05);
  gsys->SetMaximum(0.15);
  

  final_v2_g->Draw("pe same");



  TLegend *leg= new TLegend(0.74, 0.65, 0.925, 0.77);
  SetLegendStyle(leg);
  leg->SetTextSize(0.044);
  leg->AddEntry(final_v2_g,"#Upsilon(1S)","pe");
  leg->AddEntry(gv2Int2S,"#Upsilon(2S)","pe");
  leg->Draw("same");

  TLine *linebin1 = new TLine(1,-0.05,1,-0.045);
  TLine *linebin2 = new TLine(2,-0.05,2,-0.045);
  TLine *linebin3 = new TLine(3,-0.05,3,-0.045);
  linebin1->Draw("same");
  linebin2->Draw("same");
  linebin3->Draw("same");

  pad_diff->Update();


  TLatex* globtex= new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  double sz_init = 0.875; double sz_step = 0.0535;
  //globtex->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
  //globtex->DrawLatex(0.23, sz_init-sz_step, "p_{T}^{#Upsilon} < 50 GeV/c");
  //globtex->DrawLatex(0.23, sz_init-sz_step*2, "|y| < 2.4");
  globtex->DrawLatex(0.23, sz_init, "p_{T}^{#Upsilon} < 50 GeV/c");
  globtex->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
  CMS_lumi_square( pad_diff, iPeriod, iPos );
  globtex->DrawLatex(0.23, 0.14, "0-10");
  globtex->DrawLatex(0.43, 0.14, "10-30");
  globtex->DrawLatex(0.64, 0.14, "30-50");
  globtex->DrawLatex(0.85, 0.14, "50-90");

  c1->cd();
  pad_int->Draw();
  pad_int->cd();

  double sz_allign = 0.034;
  TLatex *lt1 = new TLatex();
  lt1->SetNDC(); 
  lt1->SetTextAlign(22);
  lt1->SetTextFont(42);
  lt1->SetTextSize(0.038*600./xlonger);

  double xmin_int = 0.6;
  double xmax_int = 1.8; 

  gv2Int1Ssys->GetXaxis()->SetNdivisions(101);
  gv2Int1Ssys->GetXaxis()->SetLabelSize(0);
  gv2Int1Ssys->GetYaxis()->SetTickLength(0.03*600/xlonger);
  gv2Int1Ssys->GetYaxis()->SetLabelSize(0);
  gv2Int1S->SetMarkerStyle(kFullCircle);
  
  gv2Int1Ssys->GetXaxis()->SetTitle("");
  gv2Int1Ssys->SetMinimum(-0.05);
  gv2Int1Ssys->SetMaximum(0.15);


  gv2Int1Ssys2->GetXaxis()->SetNdivisions(101);
  gv2Int1Ssys2->GetXaxis()->SetLabelSize(0);
  gv2Int1Ssys2->GetXaxis()->SetRangeUser(0,20);
  gv2Int1Ssys2->GetYaxis()->SetTickLength(0.03*600/xlonger);
  gv2Int1Ssys2->GetYaxis()->SetLabelSize(0);
  gv2Int1Ssys2->GetXaxis()->SetTitle("");


  /*
  gv2Int1Ssys->Draw("A5");
  gv2Int2Ssys->Draw("5");
  gv2Int1S->Draw("pe same");
  gv2Int2S->Draw("pe same");
  */

  THStack* ths = new THStack();
  ths->Add(hPadInt);

  ths->Draw("samenostack");

  //hPadInt->Draw();

  gv2Int1Ssys2->GetYaxis()->SetLimits(-0.05,0.15);
  gv2Int1Ssys2->GetYaxis()->SetRangeUser(-0.05,0.15);
  gv2Int1Ssys2->GetXaxis()->SetLimits(0,20);
  gv2Int1Ssys2->GetXaxis()->SetRangeUser(0,20);
  gv2Int1Ssys2->Draw("A5");
  gv2Int2Ssys2->Draw("5");
  gv2Int1S2->Draw("pe same");
  gv2Int2S2->Draw("pe same");


  line->Draw("same");

  double yext = 0.1;
  lt1->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step-sz_allign+yext, "Cent.");
  lt1->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step*2-sz_allign+yext, "0-90%");

  pad_int->Update();

  c1->Update();
  c1->SaveAs("v2Sig_cent.pdf");
 
  return;
}
