#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TASImage.h"

//
// Global variables
//

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = false;
//TString extraText   = "Supplementary";
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
//float lumiTextSize     = 0.6;
float lumiTextSize     = 0.6*1.5; // KYO
float lumiTextOffset   = 0.2;
//float cmsTextSize      = 0.75;
float cmsTextSize      = 0.75*1.5*0.6; // KYO
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_pp502TeV  = "28.0 pb^{-1}";
TString lumi_pPb502TeV  = "34.6 nb^{-1}";
TString lumi_PbPb502TeV  = "1.6 nb^{-1}";
TString lumi_PbPb502TeV_projected  = "10 nb^{-1} projection";
TString lumi_sqrtS = "";

bool drawLogo      = false;

void CMS_lumi_massPull( TPad* pad, int iPeriod=3, int iPosX=10 );

