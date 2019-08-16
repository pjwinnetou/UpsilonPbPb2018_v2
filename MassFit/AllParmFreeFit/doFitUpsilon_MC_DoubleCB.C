#include <iostream>
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/rootFitHeaders.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/commonUtility.h"
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/cutsAndBinUpsilonV2.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/PsetCollection.h"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/CMS_lumi_massPull.C"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/tdrstyle.C"
#include "/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/Style_jaebeom.h"

using namespace std;
using namespace RooFit;
void doFitUpsilon_MC_DoubleCB(
       float ptLow=0, float ptHigh=50, 
       float yLow=0, float yHigh=2.4,
       int cLow=0, int cHigh=180,
       float SiMuPtCut=3.5,
       bool recoSign=true,
       bool fixParameters=0  )
{
  float dphiEp2Low = 0 ;
  float dphiEp2High = 100 ;
  

  using namespace RooFit;
  gStyle->SetEndErrorSize(0);
 
  TString SignalCB = "Double";

  float massLow = 8.5; 
  float massHigh = 10;

  float massLowForPlot = massLow;    
  float massHighForPlot = massHigh;
  int   nMassBin  = (massHigh-massLow)*100;

  float SiMuEtaCut = 2.4;

  TFile* f1 = new TFile("/home/deathold/work/CMS/analysis/Upsilon_v2/UpsilonPbPb2018_v2/skimmedFiles/OniaTree_Skim_MC_UpsTrig_190801.root");
 
  TString recoSignString;
  if(recoSign==true) recoSignString = " && recoQQsign==0";
  else if(recoSign==false) recoSignString = " && recoQQsign!=0";

  TString kineLabel = getKineLabel (ptLow, ptHigh, yLow, yHigh, SiMuPtCut, cLow, cHigh) ;
  TString kineCut = Form("pt>%.2f && pt<%.2f && abs(y)>%.2f && abs(y)<%.2f",ptLow, ptHigh, yLow, yHigh);
  if (SiMuPtCut>0) kineCut = kineCut + Form(" && (pt1>%.2f) && (pt2>%.2f) && abs(eta1)<%.1f && abs(eta2)<%.1f", (float)SiMuPtCut, (float)SiMuPtCut, SiMuEtaCut, SiMuEtaCut );
  kineCut = kineCut + Form(" && (cBin>=%d && cBin<%d)",cLow, cHigh) + recoSignString;
  
  
  TTree* tree = (TTree*) f1->Get("mm");
  RooDataSet *dataset = (RooDataSet*)f1->Get("dataset");
  RooWorkspace *ws = new RooWorkspace("workspace");
  ws->import(*dataset);
  cout << "####################################" << endl;
  RooDataSet *reducedDS = (RooDataSet*)dataset->reduce(RooArgSet(*(ws->var("mass")), *(ws->var("pt")), *(ws->var("y"))), kineCut.Data() );
  reducedDS->SetName("reducedDS");
  ws->import(*reducedDS);
  ws->data("reducedDS")->Print();
  ws->var("mass")->setRange(massLow, massHigh);
  ws->var("mass")->Print();

  TCanvas* c1 =  new TCanvas("canvas2","My plots",4,45,550,520);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.16, 0.98, 1.0);
  pad1->SetTicks(1,1);
  pad1->Draw(); pad1->cd();
  c1->SetLeftMargin(2.6);

  RooPlot* myPlot = ws->var("mass")->frame(nMassBin); // bins
  ws->data("reducedDS")->plotOn(myPlot,Name("dataHist"));
  RooRealVar mean1s("m_{#Upsilon(1S)}","mean of the signal gaussian mass PDF",pdgMass.Y1S, pdgMass.Y1S -0.1, pdgMass.Y1S + 0.1 ) ;
  RooRealVar mRatio21("mRatio21","mRatio21",pdgMass.Y2S / pdgMass.Y1S );
  RooRealVar mRatio31("mRatio31","mRatio31",pdgMass.Y3S / pdgMass.Y1S );
  RooFormulaVar mean2s("mean2s","m_{#Upsilon(1S)}*mRatio21", RooArgSet(mean1s,mRatio21) );
  RooFormulaVar mean3s("mean3s","m_{#Upsilon(1S)}*mRatio31", RooArgSet(mean1s,mRatio31) );
         
  int collId = kAADATA; 
  PSetUpsAndBkg initPset = getUpsilonPsets( collId, ptLow, ptHigh, yLow, yHigh, cLow, cHigh, SiMuPtCut) ; 
  initPset.SetMCSgl();

  RooRealVar    sigma1s_1("sigma1s_1","width/sigma of the signal gaussian mass DF",0.1, 0.01, 0.4);
  RooFormulaVar sigma2s_1("sigma2s_1","@0*@1",RooArgList(sigma1s_1,mRatio21) );
  RooFormulaVar sigma3s_1("sigma3s_1","@0*@1",RooArgList(sigma1s_1,mRatio31) );

  RooRealVar *x1s = new RooRealVar("x1s","sigma ratio ", 0.5, 0, 1);

  RooFormulaVar sigma1s_2("sigma1s_2","@0*@1",RooArgList(sigma1s_1, *x1s) );
  RooFormulaVar sigma2s_2("sigma2s_2","@0*@1",RooArgList(sigma1s_2,mRatio21) );
  RooFormulaVar sigma3s_2("sigma3s_2","@0*@1",RooArgList(sigma1s_2,mRatio31) );
  
  RooRealVar alpha1s_1("alpha1s_1","tail shift", 2.522 , 0.5, 5.5);
  RooFormulaVar alpha2s_1("alpha2s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_1("alpha3s_1","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha1s_2("alpha1s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha2s_2("alpha2s_2","1.0*@0",RooArgList(alpha1s_1) );
  RooFormulaVar alpha3s_2("alpha3s_2","1.0*@0",RooArgList(alpha1s_1) );

  RooRealVar n1s_1("n1s_1","power order", 1.705 , 0.5, 5.5);
  RooFormulaVar n2s_1("n2s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_1("n3s_1","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n1s_2("n1s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n2s_2("n2s_2","1.0*@0",RooArgList(n1s_1) );
  RooFormulaVar n3s_2("n3s_2","1.0*@0",RooArgList(n1s_1) );

  RooRealVar *f1s = new RooRealVar("f1s","1S CB fraction", 0.5, 0., 1);
  RooFormulaVar f2s("f2s","1.0*@0",RooArgList(*f1s) );
  RooFormulaVar f3s("f3s","1.0*@0",RooArgList(*f1s) );


  
  cout << " n1s_1.getVal() = " << n1s_1.getVal() << endl;
  RooCBShape* cb1s_1 = new RooCBShape("cball1s_1", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_1, alpha1s_1, n1s_1);
  RooCBShape* cb1s_2 = new RooCBShape("cball1s_2", "cystal Ball", *(ws->var("mass")), mean1s, sigma1s_2, alpha1s_2, n1s_2);

  RooAddPdf*  cb1s = new RooAddPdf("cb1s","Signal 1S",RooArgList(*cb1s_1,*cb1s_2), RooArgList(*f1s) );
  
  RooRealVar ch4_k1("ch4_k1","ch4_k1",0.2,-4,4);
  RooRealVar ch4_k2("ch4_k2","ch4_k2",-0.1,-4,4);
  RooChebychev *bkg = new RooChebychev("cPol1Bkg","Background",*(ws->var("mass")),RooArgSet(ch4_k1,ch4_k2));

  RooRealVar *nSig1s= new RooRealVar("nSig1s"," 1S signals",0,10000000);
  RooRealVar *nBkg= new RooRealVar("nBkg"," Bkg signals",0,1000000);
  
  RooAddPdf* model = new RooAddPdf();
  model = new RooAddPdf("model","1S+2S+3S + Bkg",RooArgList(*cb1s, *bkg),RooArgList(*nSig1s,*nBkg));

  ws->import(*model);


  RooPlot* myPlot2 = (RooPlot*)myPlot->Clone();
  ws->data("reducedDS")->plotOn(myPlot2,Name("dataOS_FIT"),MarkerSize(.8));
 
  RooFitResult* fitRes2 = ws->pdf("model")->fitTo(*reducedDS,Save(), Hesse(kTRUE),Range(massLow, massHigh),Timer(kTRUE),Extended(kTRUE));
  ws->pdf("model")->plotOn(myPlot2,Name("modelHist"));
  ws->pdf("model")->plotOn(myPlot2,Name("Sig1S"),Components(RooArgSet(*cb1s)),LineColor(kOrange+7),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Name("Sig1S_CB1"),Components(RooArgSet(*cb1s_1)),LineColor(kGreen-2),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Name("Sig1S_CB2"),Components(RooArgSet(*cb1s_2)),LineColor(kMagenta+1),LineWidth(2),LineStyle(2));
  ws->pdf("model")->plotOn(myPlot2,Name("bkgPDF"),Components(RooArgSet(*bkg)),LineColor(kBlue),LineStyle(kDashed),LineWidth(2));

  myPlot2->SetFillStyle(4000);
  myPlot2->SetAxisRange(massLowForPlot, massHighForPlot,"X");
  myPlot2->GetYaxis()->SetTitleOffset(1.4);
  myPlot2->GetYaxis()->CenterTitle();
  myPlot2->GetYaxis()->SetTitleSize(0.048);
  myPlot2->GetXaxis()->SetLabelSize(0);
  myPlot2->GetXaxis()->SetRangeUser(8,14);
  myPlot2->GetXaxis()->SetTitleSize(0);
  myPlot2->Draw();
  fitRes2->Print("v");
  Double_t theNLL = fitRes2->minNll();
  cout << " *** NLL : " << theNLL << endl;
  TString perc = "%";

  float pos_text_x = 0.31;
  float pos_text_y = 0.78;
  float pos_y_diff = 0.056;
  float text_size = 15;
  int text_color = 1;
  if(ptLow==0) drawText(Form("p_{T}^{#mu#mu} < %.f GeV/c",ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else if(ptLow == 2.5 && ptHigh==5) drawText(Form("%.1f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  else drawText(Form("%.f < p_{T}^{#mu#mu} < %.f GeV/c",ptLow,ptHigh ),pos_text_x,pos_text_y,text_color,text_size);
  if(yLow==0) drawText(Form("|y^{#mu#mu}| < %.1f",yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  else if(yLow!=0) drawText(Form("%.1f < |y^{#mu#mu}| < %.1f",yLow,yHigh ), pos_text_x,pos_text_y-pos_y_diff,text_color,text_size);
  if(collId != kPPDATA && collId != kPPMCUps1S && collId != kPPMCUps2S) 
  {
      drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
      drawText("|#eta^{#mu}| < 2.4 GeV/c", pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
      drawText(Form("Centrality %d-%d%s",cLow/2,cHigh/2,perc.Data()),pos_text_x,pos_text_y-pos_y_diff*4,text_color,text_size);
  }
  else {
    drawText(Form("p_{T}^{#mu} > %.1f GeV/c", SiMuPtCut ), pos_text_x,pos_text_y-pos_y_diff*2,text_color,text_size);
    drawText("|#eta^{#mu}| < 2.4 GeV/c", pos_text_x,pos_text_y-pos_y_diff*3,text_color,text_size);
  }  
//  drawText(Form("Signal Function : %s CB", SignalCB.Data() ), 0.55,0.54,1,14);

  TLegend* fitleg = new TLegend(0.68,0.42,0.88,0.7); fitleg->SetTextSize(15);
  fitleg->SetTextFont(43);
  fitleg->SetBorderSize(0);
  fitleg->AddEntry(myPlot2->findObject("dataOS_FIT"),"Data","pe");
  fitleg->AddEntry(myPlot2->findObject("modelHist"),"Total fit","l");
  fitleg->AddEntry(myPlot2->findObject("Sig1S"),"signal","l");
  fitleg->AddEntry(myPlot2->findObject("bkgPDF"),"background","l");
  fitleg->Draw("same");

  // PULL 

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 0.98, 0.23);
  pad2->SetTopMargin(0); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.63); 
  pad1->SetLeftMargin(0.15);
  pad2->SetLeftMargin(0.15);
  pad2->SetTicks(1,1);
  pad2->cd();
  
  RooHist* hpull = myPlot2->pullHist("dataHist","modelHist");
  hpull->SetMarkerSize(0.8);
  RooPlot* pullFrame = ws->var("mass")->frame(Title("Pull Distribution")) ;
  pullFrame->addPlotable(hpull,"P") ;
  pullFrame->SetTitleSize(0);
  pullFrame->GetYaxis()->SetTitleOffset(0.31) ;
  pullFrame->GetYaxis()->SetTitle("Pull") ;
  pullFrame->GetYaxis()->SetTitleSize(0.17) ;
  pullFrame->GetYaxis()->SetLabelSize(0.13) ;
  pullFrame->GetYaxis()->SetRangeUser(-4.5,4.5) ;
//  pullFrame->GetYaxis()->SetLimits(-6,6) ;
  pullFrame->GetYaxis()->CenterTitle();

  pullFrame->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  pullFrame->GetXaxis()->SetTitleOffset(1.20) ;
  pullFrame->GetXaxis()->SetLabelOffset(0.04) ;
  pullFrame->GetXaxis()->SetLabelSize(0.183) ;
  pullFrame->GetXaxis()->SetTitleSize(0.25) ;
  pullFrame->GetXaxis()->CenterTitle();
 // pullFrame->GetXaxis()->SetTitleFont(43);
 // pullFrame->GetYaxis()->SetTitleFont(43);
  
  pullFrame->GetYaxis()->SetTickSize(0.02);
  pullFrame->GetYaxis()->SetNdivisions(505);
  pullFrame->GetXaxis()->SetTickSize(0.03);
  pullFrame->Draw() ;

  
  double chisq = 0;
  int nFullBinsPull = 0;
  int nBins = nMassBin; 
  double *ypull = hpull->GetY();
  for(int i=0;i<nBins;i++)
  {
    if(ypull[i] == 0) continue;
    chisq += TMath::Power(ypull[i],2);
    nFullBinsPull++;
  }

  int numFitPar = fitRes2->floatParsFinal().getSize();
  int ndf = nFullBinsPull - numFitPar;

  cout << Form("chi^{2}/ndf : %.3f / %d ",chisq,ndf ) << endl;

  TLine *l1 = new TLine(massLow,0,massHigh,0);
  l1->SetLineStyle(9);
  l1->Draw("same");
  pad1->Update();
              
  
  TH1D* outh = new TH1D("fitResults","fit result",20,0,20);

  outh->GetXaxis()->SetBinLabel(1,"Upsilon1S");
  
  float temp1 = ws->var("nSig1s")->getVal();  
  float temp1err = ws->var("nSig1s")->getError();  
  
  outh->SetBinContent(1,  temp1 ) ;
  outh->SetBinError  (1,  temp1err ) ;

  cout << "1S signal    =  " << outh->GetBinContent(1) << " +/- " << outh->GetBinError(1) << endl;

  setTDRStyle();
  writeExtraText = true;
  extraText = "Preliminary";

  TString label;
  label="";
  CMS_lumi_massPull(pad1, 2 ,33);


  pad1->Update();
  pad2->Update();

  c1->cd();
  pad1->Draw();
  pad2->Draw();

  pad1->Update();
  pad2->Update();

  TFile* outf = new TFile(Form("fitResMC/fitresults_upsilon_DoubleCB_%s.root",kineLabel.Data()),"recreate");
  outh->Write();
  c1->SaveAs(Form("fitResMC/fitresults_upsilon_DoubleCB_%s.pdf",kineLabel.Data()));
  c1->Write();
  ws->Write();
  cout << "N, alpha, sigma1s, M0, f, X double CB for data " << endl;
  //  void setSignalParMC(float MCn_, float MCalpha_, float MCsigma1S_, float MCm0_, float MCf_, float MCx_)
  cout << Form(" else if ( binMatched( %.f, %.f, %.f, %.1f, %.1f) ) {setSignalParMC(",SiMuPtCut, ptLow, ptHigh, yLow, yHigh);
  cout <<  ws->var("n1s_1")->getVal() << ", " <<  ws->var("alpha1s_1")->getVal() << ", "<<  ws->var("sigma1s_1")->getVal() << ", " <<  ws->var("m_{#Upsilon(1S)}")->getVal() << ", " <<  ws->var("f1s")->getVal() << ", "<<  ws->var("x1s")->getVal() << " );} " << endl;
//  outf->Close();


  ///  cout parameters :
  /*
  cout << "N, alpha, sigma1s, M0, f, X double CB for data " << endl;
  cout << "if ( (SiMuPtCut==(float)"<< SiMuPtCut<<") &&  ( ptLow == (float)"<< ptLow <<" ) && (ptHigh == (float)"<<ptHigh<<" ) && (yLow == (float)"<<yLow<<" ) && (yHigh == (float)"<<yHigh<<" ) )" << endl;

  //  void setSignalParMC(float MCn_, float MCalpha_, float MCsigma1S_, float MCm0_, float MCf_, float MCx_)
  cout << " {ret.setParMC( " ;
  cout <<  ws->var("n1s_1")->getVal() << ", " <<  ws->var("alpha1s_1")->getVal() << ", "<<  ws->var("sigma1s_1")->getVal() << ", " << endl;
  cout <<  ws->var("m_{#Upsilon(1S)}")->getVal() << ", " <<  ws->var("f1s")->getVal() << ", "<<  ws->var("x1s")->getVal() << " );} " << endl;
  */

} 
