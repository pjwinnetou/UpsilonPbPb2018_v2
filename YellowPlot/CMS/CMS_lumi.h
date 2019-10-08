#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TBox.h"
#include "TASImage.h"
#include <iostream>

//
// Global variables
//

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = true;//false;
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.5;
float lumiTextOffset   = 0.28;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.01;  // in case of ticks

float relPosX    = 0.045;
float relPosY    = 0.045;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "20.1 fb^{-1}";
TString lumi_8TeV  = "19.7 fb^{-1}";
TString lumi_7TeV  = "5.1 fb^{-1}";
TString lumi_5TeV  = "pPb 34.6 nb^{-1}";
TString lumi_PbPb2011  = "PbPb 166 #mub^{-1}";
TString lumi_PbPb2015  = "PbPb ?.? #mub^{-1}";
TString lumi_pp2013  = "pp 5.4 pb^{-1}";
TString lumi_pp2015  = "pp [Express 262163-262277]";

bool drawLogo      = false;

void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10, TString label = lumi_pp2015 );

