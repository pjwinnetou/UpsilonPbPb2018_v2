//
// dimuonYellowPlotDrawHisto.C
// Draw the dimuon invariant mass from CMS data
// Uses the output from 

//Headers{{{
#include <TROOT.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <iostream>
#include "CMS/tdrstyle.C"
#include "CMS/CMS_lumi.C"
//}}}

void dimuonYellowPlotDrawHistoOverlayPbPbTriggersExpress(){
  cout << "dimuonYellowPlot: Starting macro dimuonYellowPlot" << endl;
  cout << "dimuonYellowPlot: Setting styles..." << endl;

  // Henceforth, this macro will use only CMS TDR style macro.
  // Margins, fonts, can be handled by adjusting this
  // style macro.  Leave some of the parameters that are typically fiddled with below
  //
  // Note: the title offests and sizes need to be changed directly on the histograms
  // if they were created with different values.
  // Changing them here affect newly created histograms only. Histograms written to file
  // need to be changed each time.
  //

  setTDRStyle();
  //tdrStyle->SetPadTopMargin(0.05); // default 0.08
  //tdrStyle->SetPadBottomMargin(0.11); // default 0.12
  //tdrStyle->SetPadLeftMargin(0.10); // default 0.16 
  //tdrStyle->SetPadRightMargin(0.02); // default 0.04
  //tdrStyle->SetTitleYOffset(1.12); // default 1.25

  //
  // Open the file, this should be the file created from the dimuonYellowPlotMakeHisto.C
  // macro.
  // and create a pointer to the histogram.
  // If the pointer is not created properly, some things will work, but the
  // scaling code won't work, the part inside the for-loop will not recognize the
  // name of the histogram as a variable.
  //

  //TString inFileName1 = "dimuonMassYellowPlotHistosTriggerSeparationPbPbExpress_262548_263729_noCUT.root";
	//TString inFileName1 = "dimuonMassYellowPlotHistosTriggerSeparation381-941.root";
	TString inFileName1 = "dimuonMassYellowPlotHistosTriggerSeparation_RERECO_fixed.root";//input file name
	TString aka = "test";
	cout << "dimuonYellowPlot: Opening histogram file " << inFileName1 << endl;

	TFile *inf = new TFile(inFileName1,"READ");

//Trigger bit and names{{{
	// PbPb trigger bits for 2015:
	// Bit  1: HLT_HIL1DoubleMu0_v1
	// Bit  4: HLT_HIL1DoubleMu10_v1                    Prescale 1
	// Bit 14: HLT_HIL2DoubleMu0_Cent30_OS_NHitQ_v1     Prescale 50 - 130
	// Bit 16: HLT_HIL3DoubleMu0_Cent30_OS_m2p5to4p5_v1 Prescale 5 
	// Bit 17: HLT_HIL3DoubleMu0_Cent30_OS_m7to14_v1    Prescale 7 
	// Bit 39: HLT_HIL2Mu20_2HF_v1                      Prescale 2
	// Note, bit 39 in long is 274877906944
 
	// PbPb trigger bits for 2018:
	// Bit  1: HLT_HIL1DoubleMuOpen_v1
	// Bit 13: HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1 Prescale ?
	// Bit 14: HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v1      Prescale ?
	// Bit 24: HLT_HIL3Mu12_v1                            Prescale ?
//}}}

	const Int_t NOT = 3;//Number Of Triggers
	TH1F* massHistoUnlikeSignBit[NOT];
	TH1F* massHistoLikeSignBit[NOT];
	//const Int_t tBit[NOT] = {0, 12, 13, 23};
	const Int_t tBit[NOT] = {0, 12, 13};
	TString triggerName[NOT] = {"HLT_HIL1DoubleMuOpen_v1",
										"HLT_HIL3Mu0NHitQ10_L2Mu0_MAXdR3p5_M1to5_v1",
										"HLT_HIL3Mu2p5NHitQ10_L2Mu2_M7toinf_v1"};
	//									"HLT_HIL3Mu12_v1"};
	TString triggerNameh[NOT] = {"Double muon inclusive",
										"J/#psi region",
										"#varUpsilon + high masses"};
	unsigned int BinsTogether = 2;//nuber of bins to be merged
	for(Int_t itrg = 0; itrg < NOT; itrg++)
	{

		massHistoUnlikeSignBit[itrg] = (TH1F*) inf->Get(Form("massHistoUnlikeSignBit%d", tBit[itrg]));
		massHistoLikeSignBit[itrg] = (TH1F*) inf->Get(Form("massHistoLikeSignBit%d", tBit[itrg]));
		massHistoUnlikeSignBit[itrg]->Rebin(BinsTogether); 
		massHistoLikeSignBit[itrg]->Rebin(BinsTogether); 

/*
		massHistoUnlikeSign1[itrg] = (TH1F*) inf->Get(Form("massHistoUnlikeSignHLT%d", tBit[itrg]));
		massHistoUnlikeSign2[itrg] = (TH1F*) inf->Get(Form("massHistoUnlikeSignReco%d", tBit[itrg]));
		massHistoUnlikeSignBit[itrg] = (TH1F*) inf->Get(Form("massHistoUnlikeSignBoth%d", tBit[itrg]));
		massHistoUnlikeSign1[itrg]->Rebin(BinsTogether); 
		massHistoUnlikeSign2[itrg]->Rebin(BinsTogether); 
		massHistoUnlikeSignBit[itrg]->Rebin(BinsTogether); 
*/
	}

/*
	TH1F* massHistoUnlikeSignBit4  = (TH1F*) inf->Get("massHistoUnlikeSignBit4");
	TH1F* massHistoLikeSignBit4  = (TH1F*) inf->Get("massHistoLikeSignBit4");

	massHistoUnlikeSignBit4->Rebin(BinsTogether); 
	massHistoLikeSignBit4->Rebin(BinsTogether); 
*/
  
  if (massHistoUnlikeSignBit[0]==0x0) { 
    cout << "Histogram not found in the input file. Fatal error, can't continue! Bailing out! Agh!" << endl;
    return;
  }
/*
  cout << "Bit 1,  HLT_HIL1DoubleMu0_v1 : " <<  massHistoUnlikeSignBit1->GetName() << ", entries: " << massHistoUnlikeSignBit1->GetEntries() << endl;
  cout << "Bit 16, HLT_HIL3DoubleMu0_Cent30_OS_m2p5to4p5_v1 : " << massHistoUnlikeSignBit16->GetName() << ", entries: " << massHistoUnlikeSignBit16->GetEntries() << endl;
*/

/*
	const Int_t tFColor[NOT] = {kTeal-5, kYellow, kRed-5, kBlue-5,
										kViolet-5};
	const Int_t tLColor[NOT] = {kGreen+3, kYellow+1, kRed+3, kBlue+3,
										kViolet+3};
*/
	//const Int_t tFColor[NOT] = {393, 208, 215};
	const Int_t tFColor[NOT] = {215, 208, 393};
	//const Int_t tLColor[NOT] = {393+1, 208-1, 215-2};
	const Int_t tLColor[NOT] = {215-2, 208-1, 393+1};

  // Set up the Canvas
	TCanvas* yellowPlot = new TCanvas("yellowPlot","yellowPlot",500,500);
	yellowPlot->cd();
	yellowPlot->SetLogx();
	yellowPlot->SetLogy();


	for(Int_t itrg = 0; itrg < NOT; itrg++)
	{
		massHistoUnlikeSignBit[itrg]->SetFillColor(tFColor[itrg]);
		massHistoUnlikeSignBit[itrg]->SetLineWidth(2);
		massHistoUnlikeSignBit[itrg]->SetLineColor(tLColor[itrg]);
		if(itrg == 0)
		{
			massHistoUnlikeSignBit[itrg]->GetYaxis()->SetTitleOffset(1.33);
			massHistoUnlikeSignBit[itrg]->GetYaxis()->SetTitleSize(0.04);
			massHistoUnlikeSignBit[itrg]->GetYaxis()->SetTitleFont(42);
			massHistoUnlikeSignBit[itrg]->GetYaxis()->CenterTitle();
			massHistoUnlikeSignBit[itrg]->GetYaxis()->SetRangeUser(1e-1,4e5);
			massHistoUnlikeSignBit[itrg]->GetXaxis()->SetTitleOffset(1.);
			massHistoUnlikeSignBit[itrg]->GetXaxis()->SetTitleSize(0.045);
			massHistoUnlikeSignBit[itrg]->GetXaxis()->SetTitleFont(42);
			massHistoUnlikeSignBit[itrg]->GetXaxis()->CenterTitle();
			massHistoUnlikeSignBit[itrg]->GetXaxis()->SetRangeUser(0.5, 2e2);
                        massHistoUnlikeSignBit[itrg]->SetFillStyle(4050);
                        massHistoUnlikeSignBit[itrg]->SetFillColorAlpha(tFColor[itrg], 0.35);
			massHistoUnlikeSignBit[itrg]->Draw("");
		}
		else massHistoUnlikeSignBit[itrg]->Draw("same");
	}

	yellowPlot->SaveAs(Form("YellowPlot_%s.pdf", aka.Data()));
	yellowPlot->SaveAs(Form("YellowPlot_%s.png", aka.Data()));

  // The second histogram can cover the axes if it has a fill color.
  // Add the following line so that the axes and tick marks are visible.
	//gPad->RedrawAxis();

/*
//scale by width{{{
  // Scale by bin width, needs to be done bin by bin because bin widths are different.
  for (int i=1; i<=massHistoUnlikeSignBit1->GetNbinsX(); i++) {
    massHistoUnlikeSignBit4->SetBinContent(i, massHistoUnlikeSignBit4->GetBinContent(i)/massHistoUnlikeSignBit4->GetBinWidth(i));
    
    // Like sign
    massHistoLikeSignBit4->SetBinContent(i, massHistoLikeSignBit4->GetBinContent(i)/massHistoLikeSignBit4->GetBinWidth(i));
  }
//}}}
*/

  //
  // Options to be used with the CMS_lumi macro
  //
  
	writeExtraText = true;       // if extra text
	extraText      = "Preliminary";  // default extra text is "Preliminary"
	lumi_PbPb2011 = "PbPb 166 #mub^{-1}  pp 5.4 pb^{-1}";
  //
  // iPeriod options: 99 for pPb 5.02 TeV, 101 for PbPb 2011, 102 for pp 2013, 104 for pp 5.02 TeV, 105 for PbPb 5.02 TeV
  //
  int iPeriod      = 105; 
  lumiTextOffset   = 0.3; // default 0.28

  //lumi_5TeV = "pp, ~20 pb^{-1} lumi";
  //lumi_5TeV = "PbPb [262548-263729], Online Express Stream";
  lumi_5TeV = "PbPb [326483]";

  // Call the CMS_lumi macro to draw:
  // CMS preliminary, aligned on the right and justified (iPos=33, third argument)
  // integrated luminosity (drawn on top left, out of frame, or use lumiTextOffset)
  // center of mass energy (drawn on top right, out of frame, or use lumiTextOffset)
  
  //CMS_lumi(yellowPlot,iPeriod,33);

  // Note, the code for putting the CMS logo is only used
  // in the style macro if one does not write CMS preliminary.
  // For writing both, one needs to put one of them by hand.
  // This is a copy-paste of the code to insert the CMS logo, changing it to
  // a color logo. Needs to have the CMS-Color.gif file in the same directory.
  /*
  float H = gPad->GetWh();
  float W = gPad->GetWw();
  float l = gPad->GetLeftMargin();
  float t = gPad->GetTopMargin();
  float r = gPad->GetRightMargin();
  float b = gPad->GetBottomMargin();
  float posX_ =   l + 0.025*(1-l-r)*W/H;
  float posY_ = 1-t - 0.025*(1-t-b);
  float xl_0 = posX_;
  float yl_0 = posY_ - 0.1;
  float xl_1 = posX_ + 0.1*H/W;
  float yl_1 = posY_;
  TASImage* CMS_logo = new TASImage("CMS-Color.gif");
  TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
  pad_logo->Draw();
  pad_logo->cd();
  CMS_logo->Draw("X");
  pad_logo->Modified();
  yellowPlot->cd();
  */

  // Booking the TLatex class and set the parameters
  // CMS Preliminary (not used anymore) left as a dummy
  // Each of the peaks is done with a call to DrawLatex, which
  // makes a copy of the TLatex class style, so the style should
  // be set with the first TLatex instance, and that style will 
  // be inherited by each of the copies.

  TLatex* latex = new TLatex(0.62,0.88,"CMS Preliminary");
  latex->SetNDC();
  latex->SetTextColor(kBlack);
  latex->SetTextFont(42);
  //latex->DrawLatex(0.24,0.88,"#rho, #omega");
  //latex->DrawLatex(0.32,0.88,"#phi");
  
//{{{
	cout << "dimuonYellowPlot: Scaling histograms..." << endl;

	TCanvas* yellowPlotScale = new TCanvas("yellowPlotScale","yellowPlot",500,500);
	yellowPlotScale->cd();
	yellowPlotScale->SetLogx();
	yellowPlotScale->SetLogy();

//bin width scale{{{
	for(Int_t itrg = 0; itrg < NOT; itrg++)
	{
		for(Int_t ibin = 1; ibin <= massHistoUnlikeSignBit[0]->GetNbinsX(); ibin++)
		{
			massHistoUnlikeSignBit[itrg]->SetBinContent(ibin, massHistoUnlikeSignBit[itrg]->GetBinContent(ibin)/massHistoUnlikeSignBit[itrg]->GetBinWidth(ibin));
			massHistoLikeSignBit[itrg]->SetBinContent(ibin, massHistoLikeSignBit[itrg]->GetBinContent(ibin)/massHistoLikeSignBit[itrg]->GetBinWidth(ibin));
		}
	}
//}}}

	for(Int_t itrg = 0; itrg < NOT; itrg++)
	{
		if(itrg == 0)
		{
			massHistoUnlikeSignBit[itrg]->GetYaxis()->SetRangeUser(1e-1,8e7);
			massHistoUnlikeSignBit[itrg]->Draw("");
		}
		else massHistoUnlikeSignBit[itrg]->Draw("same");
	}

//Draw legend{{{
	TLegend* leg = new TLegend(0.2,0.73,0.75,0.88);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->SetHeader("Trigger selections");

	for(Int_t itrg = 0; itrg < NOT; itrg++)
	{
		leg->AddEntry(massHistoUnlikeSignBit[itrg], Form("%s", triggerNameh[itrg].Data()), "F");
	}
	leg->Draw();
//}}}

	massHistoUnlikeSignBit[0]->Draw("same");

//Draw text{{{
	latex->SetTextSize(0.04);
	latex->DrawLatex(0.17,0.94,"PbPb 1.6 nb^{-1}");
	latex->DrawLatex(0.70,0.94,"#sqrt{s_{NN}} = 5.02 TeV");
	latex->DrawLatex(0.43,0.68,"J/#psi");
	//latex->DrawLatex(0.32,0.70,"#psi(2S)");
	latex->DrawLatex(0.58,0.68,"#varUpsilon(1,2,3S)");
	latex->DrawLatex(0.84,0.54,"Z");
	latex->DrawLatex(0.25,0.25,"p_{T}^{#mu} > 4 GeV/c");

	Double_t txtweight = yellowPlotScale->GetTopMargin();
	latex->SetTextAlign(11);
	latex->SetTextFont(61);
	latex->SetTextSize(0.75*txtweight);
	latex->DrawLatex(0.73, 0.85, "CMS");
	latex->SetTextFont(52);
	latex->SetTextSize(0.75*0.76*txtweight);
	latex->DrawLatex(0.73, 0.80, "Preliminary");
//}}}

	gPad->RedrawAxis();
	yellowPlotScale->SaveAs(Form("YellowPlot_binwidth_scaled_%s.pdf", aka.Data()));
	yellowPlotScale->SaveAs(Form("YellowPlot_binwidth_scaled_%s.png", aka.Data()));
//}}}

/*
//Draw each trigger separately{{{
  TCanvas* dimuonMassBit4 = new TCanvas("dimuonMassBit4","dimuonMassBit4",500,500);
  massHistoUnlikeSignBit4->SetMaximum(1e3);
  massHistoUnlikeSignBit4->Draw();

  massHistoLikeSignBit4->SetFillColor(kOrange-3);
  massHistoLikeSignBit4->SetLineWidth(2);
  massHistoLikeSignBit4->SetLineColor(kOrange+7);

  massHistoLikeSignBit4->Draw("same");
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->RedrawAxis();

  TLegend* legMassBit4 = new TLegend(0.5,0.7,0.9,0.9);
  legMassBit4->SetBorderSize(0);
  legMassBit4->SetFillColor(0);
  legMassBit4->SetTextSize(0.03);
  legMassBit4->AddEntry(massHistoUnlikeSignBit4,"HIL1DoubleMu10","F");
  legMassBit4->AddEntry(massHistoLikeSignBit4,"Like-sign","F");
  legMassBit4->Draw();

  TCanvas* dimuonMassBit14 = new TCanvas("dimuonMassBit14","dimuonMassBit14",500,500);
  //massHistoUnlikeSignBit14->SetMaximum(3e4);
  massHistoUnlikeSignBit14->Draw();

  massHistoLikeSignBit14->SetFillColor(kRed-4);
  massHistoLikeSignBit14->SetLineWidth(2);
  massHistoLikeSignBit14->SetLineColor(kRed+1);

  massHistoLikeSignBit14->Draw("same");

  gPad->SetLogy();
  gPad->SetLogx();
  gPad->RedrawAxis();

  TLegend* legMassBit14 = new TLegend(0.35,0.7,0.9,0.9);
  legMassBit14->SetBorderSize(0);
  legMassBit14->SetFillColor(0);
  legMassBit14->SetTextSize(0.03);
  legMassBit14->AddEntry(massHistoUnlikeSignBit14,"HIL2DoubleMu0_Cent30_OS_NHitQ","F");
  legMassBit14->AddEntry(massHistoLikeSignBit14,"Like-sign","F");
  legMassBit14->Draw();



  TCanvas* dimuonMassBit16 = new TCanvas("dimuonMassBit16","dimuonMassBit16",500,500);
  massHistoUnlikeSignBit16->SetMaximum(1e5);
  massHistoUnlikeSignBit16->GetXaxis()->SetRangeUser(.1,2e2);
  massHistoUnlikeSignBit16->Draw();

  massHistoLikeSignBit16->SetFillColor(kAzure-4);
  massHistoLikeSignBit16->SetLineWidth(2);
  massHistoLikeSignBit16->SetLineColor(kAzure+3);

  massHistoLikeSignBit16->Draw("same");
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->RedrawAxis();

  TLegend* legMassBit16 = new TLegend(0.15,0.7,0.6,0.9);
  legMassBit16->SetBorderSize(0);
  legMassBit16->SetFillColor(0);
  legMassBit16->SetTextSize(0.03);
  legMassBit16->AddEntry(massHistoUnlikeSignBit16,"HIL3DoubleMu0_Cent30_OS_m2p5to4p5","F");
  legMassBit16->AddEntry(massHistoLikeSignBit16,"Like-sign","F");
  legMassBit16->Draw();


  TCanvas* dimuonMassBit17 = new TCanvas("dimuonMassBit17","dimuonMassBit17",500,500);
  massHistoUnlikeSignBit17->SetMaximum(1e5);
  massHistoUnlikeSignBit17->GetXaxis()->SetRangeUser(.1,2e2);
  massHistoUnlikeSignBit17->Draw();

  massHistoLikeSignBit17->SetFillColor(kSpring);
  massHistoLikeSignBit17->SetLineWidth(2);
  massHistoLikeSignBit17->SetLineColor(kSpring-7);

  massHistoLikeSignBit17->Draw("same");
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->RedrawAxis();

  TLegend* legMassBit17 = new TLegend(0.15,0.7,0.6,0.9);
  legMassBit17->SetBorderSize(0);
  legMassBit17->SetFillColor(0);
  legMassBit17->SetTextSize(0.03);
  legMassBit17->AddEntry(massHistoUnlikeSignBit17,"HIL3DoubleMu0_Cent30_OS_m7to14","F");
  legMassBit17->AddEntry(massHistoLikeSignBit17,"Like-sign","F");
  legMassBit17->Draw();

  TCanvas* dimuonMassBit39 = new TCanvas("dimuonMassBit39","dimuonMassBit39",500,500);
  massHistoUnlikeSignBit39->SetMaximum(1e4);
  massHistoUnlikeSignBit39->Draw();

  massHistoLikeSignBit39->SetFillColor(kMagenta-4);
  massHistoLikeSignBit39->SetLineWidth(2);
  massHistoLikeSignBit39->SetLineColor(kMagenta+3);

  massHistoLikeSignBit39->Draw("same");
  gPad->SetLogy();
  gPad->SetLogx();
  gPad->RedrawAxis();

  TLegend* legMassBit39 = new TLegend(0.15,0.75,0.9,0.9);
  legMassBit39->SetBorderSize(0);
  legMassBit39->SetFillColor(0);
  legMassBit39->SetTextSize(0.03);
  legMassBit39->AddEntry(massHistoUnlikeSignBit39,"HIL2Mu20_2HF","F");
  legMassBit39->AddEntry(massHistoLikeSignBit39,"Like-sign","F");
  legMassBit39->Draw();
//}}}

//Draw centralitiy distribution{{{
  TCanvas* centralityBit4cnv = new TCanvas("centralityBit4cnv","centralityBit4cnv",500,500);
  
  centralityBit4->SetFillColor(kYellow);
  centralityBit4->SetLineWidth(2);
  centralityBit4->SetLineColor(kYellow+1);

  centralityBit14->SetFillColor(kRed-5);
  centralityBit14->SetLineWidth(2);
  centralityBit14->SetLineColor(kRed+3);

  centralityBit16->SetFillColor(kBlue-5);
  centralityBit16->SetLineWidth(2);
  centralityBit16->SetLineColor(kBlue+3);

  centralityBit17->SetFillColor(kTeal-5);
  centralityBit17->SetLineWidth(2);
  centralityBit17->SetLineColor(kGreen+3);

  centralityBit39->SetFillColor(kViolet-5);
  centralityBit39->SetLineWidth(2);
  centralityBit39->SetLineColor(kViolet+3);

  centralityBit4->Draw();
  gPad->SetLogy();

  TLegend* legCentBit4 = new TLegend(0.5,0.7,0.9,0.9);
  legCentBit4->SetBorderSize(0);
  legCentBit4->SetFillColor(0);
  legCentBit4->SetTextSize(0.03);
  legCentBit4->AddEntry(centralityBit4,"HIL1DoubleMu10","F");
  legCentBit4->Draw();
//}}}

//Draw centrality distribution for each trigger separately{{{
  TCanvas* centralityBit14cnv = new TCanvas("centralityBit14cnv","centralityBit14cnv",500,500);
  centralityBit14->Draw();
  gPad->SetLogy();

  TLegend* legCentBit14 = new TLegend(0.2,0.75,0.9,0.9);
  legCentBit14->SetBorderSize(0);
  legCentBit14->SetFillColor(0);
  legCentBit14->SetTextSize(0.03);
  legCentBit14->AddEntry(centralityBit14,"HIL2DoubleMu0_Cent30_OS_NHitQ","F");
  legCentBit14->Draw();

  TCanvas* centralityBit16cnv = new TCanvas("centralityBit16cnv","centralityBit16cnv",500,500);
  centralityBit16->Draw();
  gPad->SetLogy();

  TLegend* legCentBit16 = new TLegend(0.2,0.75,0.9,0.9);
  legCentBit16->SetBorderSize(0);
  legCentBit16->SetFillColor(0);
  legCentBit16->SetTextSize(0.03);
  legCentBit16->AddEntry(centralityBit16,"HIL3DoubleMu0_Cent30_OS_m2p5to4p5","F");
  legCentBit16->Draw();

  TCanvas* centralityBit17cnv = new TCanvas("centralityBit17cnv","centralityBit17cnv",500,500);
  centralityBit17->Draw();
  gPad->SetLogy();

  TLegend* legCentBit17 = new TLegend(0.2,0.75,0.9,0.9);
  legCentBit17->SetBorderSize(0);
  legCentBit17->SetFillColor(0);
  legCentBit17->SetTextSize(0.03);
  legCentBit17->AddEntry(centralityBit17,"HIL3DoubleMu0_Cent30_OS_m7to14","F");
  legCentBit17->Draw();

  TCanvas* centralityBit39cnv = new TCanvas("centralityBit39cnv","centralityBit39cnv",500,500);
  centralityBit39->Draw();
  gPad->SetLogy();
  
  TLegend* legCentBit39 = new TLegend(0.5,0.7,0.9,0.9);
  legCentBit39->SetBorderSize(0);
  legCentBit39->SetFillColor(0);
  legCentBit39->SetTextSize(0.03);
  legCentBit39->AddEntry(centralityBit39,"HIL2Mu20_2HF","F");
  legCentBit39->Draw();
//}}}

//Draw Max daughter muon pt{{{
  maxDaughterMuonPtBit4->SetFillColor(kYellow);
  maxDaughterMuonPtBit4->SetLineWidth(2);
  maxDaughterMuonPtBit4->SetLineColor(kYellow+1);

  maxDaughterMuonPtBit14->SetFillColor(kRed-5);
  maxDaughterMuonPtBit14->SetLineWidth(2);
  maxDaughterMuonPtBit14->SetLineColor(kRed+3);

  maxDaughterMuonPtBit16->SetFillColor(kBlue-5);
  maxDaughterMuonPtBit16->SetLineWidth(2);
  maxDaughterMuonPtBit16->SetLineColor(kBlue+3);

  maxDaughterMuonPtBit17->SetFillColor(kTeal-5);
  maxDaughterMuonPtBit17->SetLineWidth(2);
  maxDaughterMuonPtBit17->SetLineColor(kGreen+3);

  maxDaughterMuonPtBit39->SetFillColor(kViolet-5);
  maxDaughterMuonPtBit39->SetLineWidth(2);
  maxDaughterMuonPtBit39->SetLineColor(kViolet+3);
//}}}

//Draw max dauther muon pt for each trigger separately{{{
  TCanvas* maxDaughterMuonPtBit4cnv = new TCanvas("maxDaughterMuonPtBit4cnv","max daughter muon pt, Bit4",500,500);

  maxDaughterMuonPtBit4->Draw();
  gPad->SetLogy();

  TLegend* legDauMuPtBit4 = new TLegend(0.5,0.7,0.9,0.9);
  legDauMuPtBit4->SetBorderSize(0);
  legDauMuPtBit4->SetFillColor(0);
  legDauMuPtBit4->SetTextSize(0.03);
  legDauMuPtBit4->AddEntry(maxDaughterMuonPtBit4,"HIL1DoubleMu10","F");
  legDauMuPtBit4->Draw();

  TCanvas* maxDaughterMuonPtBit14cnv = new TCanvas("maxDaughterMuonPtBit14cnv","max daughter muon pt, Bit14",500,500);

  maxDaughterMuonPtBit14->Draw();
  gPad->SetLogy();

  TLegend* legDauMuPtBit14 = new TLegend(0.2,0.75,0.9,0.9);
  legDauMuPtBit14->SetBorderSize(0);
  legDauMuPtBit14->SetFillColor(0);
  legDauMuPtBit14->SetTextSize(0.03);
  legDauMuPtBit14->AddEntry(maxDaughterMuonPtBit14,"HIL2DoubleMu0_Cent30_OS_NHitQ","F");
  legDauMuPtBit14->Draw();

  TCanvas* maxDaughterMuonPtBit16cnv = new TCanvas("maxDaughterMuonPtBit16cnv","max daughter muon pt, Bit16",500,500);

  maxDaughterMuonPtBit16->Draw();
  gPad->SetLogy();

  TLegend* legDauMuPtBit16 = new TLegend(0.2,0.75,0.9,0.9);
  legDauMuPtBit16->SetBorderSize(0);
  legDauMuPtBit16->SetFillColor(0);
  legDauMuPtBit16->SetTextSize(0.03);
  legDauMuPtBit16->AddEntry(maxDaughterMuonPtBit16,"HIL3DoubleMu0_Cent30_OS_m2p5to4p5","F");
  legDauMuPtBit16->Draw();

  TCanvas* maxDaughterMuonPtBit17cnv = new TCanvas("maxDaughterMuonPtBit17cnv","max daughter muon pt, Bit17",500,500);

  maxDaughterMuonPtBit17->Draw();
  gPad->SetLogy();

  TLegend* legDauMuPtBit17 = new TLegend(0.2,0.75,0.9,0.9);
  legDauMuPtBit17->SetBorderSize(0);
  legDauMuPtBit17->SetFillColor(0);
  legDauMuPtBit17->SetTextSize(0.03);
  legDauMuPtBit17->AddEntry(maxDaughterMuonPtBit17,"HIL3DoubleMu0_Cent30_OS_m7to14","F");
  legDauMuPtBit17->Draw();

  TCanvas* maxDaughterMuonPtBit39cnv = new TCanvas("maxDaughterMuonPtBit39cnv","max daughter muon pt, Bit39",500,500);

  maxDaughterMuonPtBit39->Draw();
  gPad->SetLogy();
  
  TLegend* legDauMuPtBit39 = new TLegend(0.5,0.7,0.9,0.9);
  legDauMuPtBit39->SetBorderSize(0);
  legDauMuPtBit39->SetFillColor(0);
  legDauMuPtBit39->SetTextSize(0.03);
  legDauMuPtBit39->AddEntry(maxDaughterMuonPtBit39,"HIL2Mu20_2HF","F");
  legDauMuPtBit39->Draw();
//}}}
*/
  return;
}
