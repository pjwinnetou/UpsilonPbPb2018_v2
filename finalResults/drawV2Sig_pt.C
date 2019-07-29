#include <iostream>
#include "../Style_jaebeom.h"
#include "../commonUtility.h"
#include "../cutsAndBinUpsilonV2.h"
#include "../CMS_lumi_square.C"
#include "../tdrstyle.C"

using namespace std;

void drawV2Sig_pt(){

  setTDRStyle();
  writeExtraText= true;
  int iPeriod = 2;
  int iPos = 33;
  
  const int nBin=3;
  const int nCentBin = 4;
  double exsys[nBin] =  {1.5, 1.5,22};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nBin][nCentBin];
  TFile* fIn_int[nBin];
  double ptBin[nBin+1] = {0, 3, 6, 50};
  double centBin[nCentBin+1] = {0,10,30,50,90};

	TGraphErrors* gv2[nBin][nCentBin];
	TGraphErrors* gv2_int[nBin];
  
  TH1D *final_v2[nCentBin];
  
  for (int is=0; is<nBin; is++){
  	fIn_int[is] = new TFile(Form("../SimFit/SigAllFreeFitFix/FitResult/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality0-180_m8-14_OS.root",ptBin[is],ptBin[is+1]),"READ");
    gv2_int[is]=(TGraphErrors*)fIn_int[is]->Get("v2vspt");
    for(int ic = 0; ic<nCentBin; ic++){
      fIn[is][ic] = new TFile(Form("../SimFit/SigAllFreeFitFix/FitResult/SimFitResult_pt%.1f-%.1f_y0.0-2.4_muPt3.5_centrality%.f-%.f_m8-14_OS.root",ptBin[is],ptBin[is+1],centBin[ic]*2,centBin[ic+1]*2),"READ");
      gv2[is][ic]=(TGraphErrors*)fIn[is][ic]->Get("v2vspt");
    }
  }

  for(int ic=0; ic<nCentBin; ic++){
    final_v2[ic] = new TH1D(Form("Final_v2_Cent%d",ic),";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,ptBin);
  }
  TH1D *final_v2_int = new TH1D("Final_v2_int",";p_{T}^{#Upsilon(1S)} (GeV/c);v_{2}^{Sig}",nBin,ptBin);
  

  //// get bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;

  for(int ic=0; ic<nCentBin;ic++){
    cout << "Cent : " << ic << endl;
    for (int is=0; is<nBin; is++){
      cout << is+1 <<"th bin***************" << endl;
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gv2[is][ic] -> GetPoint(0,pxtmp,pytmp);
      eytmp = gv2[is][ic]->GetErrorY(0);
      final_v2[ic]->SetBinContent(is+1, pytmp);
      final_v2[ic]->SetBinError(is+1, eytmp);
    }
    final_v2[ic]->GetXaxis()->CenterTitle();
    final_v2[ic]->GetYaxis()->CenterTitle();
    final_v2[ic]->GetYaxis()->SetTitleOffset(1.2);
    final_v2[ic]->GetXaxis()->SetTitleOffset(1.);
    final_v2[ic]->SetMinimum(-0.05);
    final_v2[ic]->SetMaximum(0.2);
  }


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

  TGraphErrors* final_v2_g[nCentBin];
  TFile *fSys = new TFile("../Systematic/merged_sys.root","read");
  TH1D  *hSys[nCentBin];
  TGraphErrors *gsys[nCentBin];

  for(int ic = 0; ic<nCentBin; ic++){
   final_v2_g[ic] = new TGraphErrors(final_v2[ic]);
   hSys[ic] = (TH1D*) fSys->Get(Form("hpt_Cent%.f%.f_merged",centBin[ic],centBin[ic+1]));
   gsys[ic] = new TGraphErrors();
   for(int ib=0;ib<final_v2[ic]->GetNbinsX();ib++){
     gsys[ic]->SetPoint(ib, final_v2[ic]->GetBinCenter(ib+1), final_v2[ic]->GetBinContent(ib+1));
     gsys[ic]->SetPointError(ib, final_v2[ic]->GetBinWidth(ib+1)/2, fabs(hSys[ic]->GetBinContent(ib+1)));
     final_v2_g[ic]->SetPointError(ib,0,final_v2[ic]->GetBinError(ib+1));
   }
   SetGraphStyle(final_v2_g[ic],1,1);
   SetGraphStyleSys(gsys[ic],1);
  }

  TGraphErrors* final_v2_int_g = new TGraphErrors(final_v2_int);
  TH1D* hsys_int = (TH1D*) fSys->Get("hpt_Cent090_merged");
  TGraphErrors *gsys_int = new TGraphErrors();
  for(int ib=0;ib<final_v2_int->GetNbinsX();ib++){
    gsys_int->SetPoint(ib, final_v2_int->GetBinCenter(ib+1), final_v2_int->GetBinContent(ib+1));
    gsys_int->SetPointError(ib, final_v2_int->GetBinWidth(ib+1)/2, fabs(hsys_int->GetBinContent(ib+1)));
    final_v2_int_g->SetPointError(ib,0,final_v2_int->GetBinError(ib+1));
  }

  SetGraphStyle(final_v2_int_g,1,1);
  SetGraphStyleSys(gsys_int,1);

  double xmax = 50; 
 
  TCanvas* c1[nCentBin];
  TLegend *leg[nCentBin];
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  double sz_init = 0.867; double sz_step = 0.0535;
  for(int ic=0; ic<nCentBin; ic++){
    c1[ic]= new TCanvas(Form("c%d",ic),Form("c%d",ic),600,600);
    gPad->SetLeftMargin(0.19);
    gPad->SetBottomMargin(0.17);
    gPad->SetTopMargin(0.065);
    c1[ic]->cd();
    final_v2_g[ic]->SetMarkerStyle(kFullCircle);
    gsys[ic]->Draw("A5");
    gsys[ic]->GetXaxis()->CenterTitle();
    gsys[ic]->GetYaxis()->CenterTitle();
    gsys[ic]->GetYaxis()->SetTitleOffset(1.5);
    gsys[ic]->GetXaxis()->SetTitleOffset(1.1);
    gsys[ic]->GetYaxis()->SetTitle("v_{2}^{sig}");
    gsys[ic]->GetXaxis()->SetTitle("p_{T}^{#Upsilon(1S)} (GeV/c)");
    gsys[ic]->GetXaxis()->SetRangeUser(0,xmax);
    gsys[ic]->SetMinimum(-0.05);
    gsys[ic]->SetMaximum(0.15);

    final_v2_g[ic]->Draw("pe same");

    TString perc = "%";
    leg[ic]= new TLegend(0.74, 0.70, 0.925, 0.77);
    SetLegendStyle(leg[ic]);
    leg[ic]->SetTextSize(0.044);
    leg[ic]->AddEntry(final_v2_g[ic],"#Upsilon(1S)","pe");
    leg[ic]->Draw("same");
    dashedLine(0.,0.,xmax,0.,1,1);
    globtex->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
    globtex->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
    globtex->DrawLatex(0.23, sz_init-sz_step*2, Form("Cent. %.f-%.f%s",centBin[ic],centBin[ic+1],perc.Data()));
    gsys[ic]->SetName(Form("gr_sys_v2_vs_pt_Cent%.f%.f",centBin[ic],centBin[ic+1]));
    final_v2_g[ic]->SetName(Form("gr_point_v2_vs_pt_Cent%.f%.f",centBin[ic],centBin[ic+1]));
    CMS_lumi_square( c1[ic], iPeriod, iPos );
    c1[ic]->Update();
    c1[ic]->SaveAs(Form("v2Sig_pt_Cent%d.pdf",ic));
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
  final_v2_int_g->SetMarkerStyle(kFullCircle);
  gsys_int->Draw("A5");
  gsys_int->GetXaxis()->CenterTitle();
  gsys_int->GetYaxis()->CenterTitle();
  gsys_int->GetYaxis()->SetTitleOffset(1.5);
  gsys_int->GetXaxis()->SetTitleOffset(1.1);
  gsys_int->GetYaxis()->SetTitle("v_{2}^{sig}");
  gsys_int->GetXaxis()->SetTitle("p_{T}^{#Upsilon(1S)} (GeV/c)");
  gsys_int->GetXaxis()->SetRangeUser(0,xmax);
  gsys_int->SetMinimum(-0.05);
  gsys_int->SetMaximum(0.15);
  
  final_v2_int_g->Draw("pe same");
  
  TLegend *leg_int= new TLegend(0.74, 0.70, 0.925, 0.77);
  SetLegendStyle(leg_int);
  leg_int->SetTextSize(0.044);
  leg_int->AddEntry(final_v2_int_g,"#Upsilon(1S)","pe");
  leg_int->Draw("same");
  dashedLine(0.,0.,xmax,0.,1,1);
  
  globtex_int->DrawLatex(0.23, sz_init, "p_{T}^{#mu} > 3.5 GeV/c");
  globtex_int->DrawLatex(0.23, sz_init-sz_step, "|y| < 2.4");
  globtex_int->DrawLatex(0.23, sz_init-sz_step*2, "Cent. 0-90%");
  
  CMS_lumi_square( c2, iPeriod, iPos );
	c2->Update();
  c2->SaveAs("v2Sig_pt_int.pdf");




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


	return;
} 
