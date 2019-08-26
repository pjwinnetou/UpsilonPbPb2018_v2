#include <iostream>
#include "TFile.h"
#include "TString.h"
#include "TH2D.h"
#include "TDirectory.h"

static const int NIOV = 10;

unsigned int runs[] = {181611, 181671, 181701, 181801, 182001, 182201, 182401, 182601, 182801, 182951, 183101};


class TrackEfficiency{
  
public:
  TrackEfficiency(TString effOption);
  double getEfficiencies(double pt, double cent, double phi, double eta);
  int setRunNumber(unsigned int runno_);
  
private:
  TFile * f_eff;
  
  TH2D * eff[10];
  int IOVDir;
  
  
  bool DATA;
  
};

//-----------------------------------------------------------------------------
//                                  METHODS
//-----------------------------------------------------------------------------

TrackEfficiency::TrackEfficiency(TString effOption){
        
        f_eff = new TFile("Eff.root");
	setRunNumber(runs[0]);    
 
    
    
}
//-----------------------------------------------------------------------------

int TrackEfficiency::setRunNumber(unsigned int runno_){
    
  IOVDir = -1;
  for(int iIOV = 0; iIOV < NIOV; iIOV++){
    if(runno_ >= runs[iIOV] && runno_ < runs[iIOV+1]) {
      IOVDir = iIOV;
      break;
    }
  }
  if(IOVDir<0) return 0;
  for(int i = 0; i<10; i++) {
    eff[i] = (TH2D *) f_eff->Get(Form("%d/Eff_%d_%d",runs[IOVDir],10*i,10*i+10));
  }
  return 1;
}

//-----------------------------------------------------------------------------

double TrackEfficiency::getEfficiencies(double pt, double cent, double phi, double eta){
    
  int icent = cent/10;
  if(icent>=0 && icent<10) {
    return eff[icent]->GetBinContent( eff[icent]->GetXaxis()->FindBin(phi),eff[icent]->GetYaxis()->FindBin(eta));
    
  }
  
  return 1.;
}

//-----------------------------------------------------------------------------





    
