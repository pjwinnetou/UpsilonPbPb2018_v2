#ifndef tnp_weight_lowptPbPb_h
#define tnp_weight_lowptPbPb_h
#include "TMath.h"

// IN THIS FILE YOU WILL FIND:
// +++++++++++++++++++++++++++++++++++++++
// - Trigger: (tnp_weight_trg_pbpb)
//   * filterId = 0: Jpsi L2 filter
//   * filterId = 1: Jpsi L3 filter
//   * filterId = 2: Upsi L2 filter
//   * filterId = 3: Upsi L3 filter
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// - MuID: (tnp_weight_muid_pbpb)
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// - Inner tracking: (tnp_weight_trk_pbpb)
//   * idx = 0:  nominal
//   * idx = -1: syst variation, +1 sigma
//   * idx = -2: syst variation, -1 sigma
//   * idx = +1: stat variation, +1 sigma
//   * idx = +2: stat variation, -1 sigma
// +++++++++++++++++++++++++++++++++++++++

double tnp_weight_muid_pbpb(double pt, double eta, int idx=0);
double tnp_weight_trk_pbpb(double eta, int idx=0);
double tnp_weight_trg_pbpb(double pt, double eta, int filterId=0,int idx=0);

///////////////////////////////////////////////////
//             M u I D    P b P b                //
///////////////////////////////////////////////////

double tnp_weight_muid_pbpb(double pt, double eta, int idx) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  // SF for 0 < |eta| < 1.2
  if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
    if (x >= 3.5 && x <4) {num = 0.946125; den = 0.95062; statUp = 0.00377662; statDown = 0.00386775;}
    else if (x >= 4 && x <4.5) {num = 0.973366; den = 0.975256; statUp = 0.00284959; statDown = 0.00295956;}
    else if (x >= 4.5 && x <5) {num = 0.984622; den = 0.983749; statUp = 0.00263774; statDown = 0.00277141;}
    else if (x >= 5 && x <5.5) {num = 0.985061; den = 0.98948; statUp = 0.00285151; statDown = 0.00301546;}
    else if (x >= 5.5 && x <6.5) {num = 0.988937; den = 0.99051; statUp = 0.0020037; statDown = 0.00211913;}
    else if (x >= 6.5 && x <8) {num = 0.991257; den = 0.992598; statUp = 0.0017347; statDown = 0.00185756;}
    else if (x >= 8 && x <10.5) {num = 0.986505; den = 0.9885; statUp = 0.00213856; statDown = 0.00228346;}
    else if (x >= 10.5 && x <14) {num = 0.975903; den = 0.967524; statUp = 0.00364064; statDown = 0.0038785;}
    else if (x >= 14 && x <18) {num = 0.979888; den = 0.978318; statUp = 0.00488479; statDown = 0.00538685;}
    else {num = 0.983224; den = 0.984461; statUp = 0.0056352; statDown = 0.00635297;}
  }
  // SF for 1.2 < |eta| < 1.8
  if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
    if (x >= 2.37 && x <3) {num = 0.983353; den = 0.981734; statUp = 0.00668073; statDown = 0.00681507;}
    else if (x >= 3 && x <3.5) {num = 0.985622; den = 0.984888; statUp = 0.00411601; statDown = 0.00425201;}
    else if (x >= 3.5 && x <4) {num = 0.980436; den = 0.989513; statUp = 0.00390367; statDown = 0.00408228;}
    else if (x >= 4 && x <4.5) {num = 0.986109; den = 0.991777; statUp = 0.00358665; statDown = 0.00382294;}
    else if (x >= 4.5 && x <5) {num = 0.992053; den = 0.993039; statUp = 0.00353145; statDown = 0.00382873;}
    else if (x >= 5 && x <6) {num = 0.989106; den = 0.994324; statUp = 0.00288884; statDown = 0.00311446;}
    else if (x >= 6 && x <7.5) {num = 0.995133; den = 0.997452; statUp = 0.00241041; statDown = 0.00267219;}
    else if (x >= 7.5 && x <10) {num = 0.998692; den = 0.996495; statUp = 0.00130754; statDown = 0.00252507;}
    else if (x >= 10 && x <15) {num = 0.979734; den = 0.971008; statUp = 0.0202658; statDown = 0.00568435;}
    else {num = 0.996433; den = 0.987155; statUp = 0.0035675; statDown = 0.00767196;}
  }
  // SF for 1.8 < |eta| < 2.1
  if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
    if (x >= 1.8 && x <2.5) {num = 0.970771; den = 0.955555; statUp = 0.010087; statDown = 0.0102019;}
    else if (x >= 2.5 && x <3) {num = 0.984571; den = 0.983608; statUp = 0.00598521; statDown = 0.00625035;}
    else if (x >= 3 && x <3.5) {num = 0.991532; den = 0.99085; statUp = 0.00485266; statDown = 0;}
    else if (x >= 3.5 && x <4) {num = 0.982929; den = 0.995395; statUp = 0.00477394; statDown = 0.00523748;}
    else if (x >= 4 && x <4.5) {num = 0.997191; den = 0.998022; statUp = 0.00280897; statDown = 0.00374539;}
    else if (x >= 4.5 && x <5.5) {num = 0.996918; den = 0.997783; statUp = 0.00226149; statDown = 0.00267054;}
    else if (x >= 5.5 && x <7) {num = 1; den = 0.999082; statUp = 9.82446e-09; statDown = 0.00120742;}
    else if (x >= 7 && x <9) {num = 0.996645; den = 0.999904; statUp = 0.00294462; statDown = 0.00358784;}
    else if (x >= 9 && x <12) {num = 0.997479; den = 0.997288; statUp = 0.00252069; statDown = 0.00478073;}
    else {num = 0.989358; den = 0.997918; statUp = 0.00514686; statDown = 0.00685444;}
  }
  // SF for 2.1 < |eta| < 2.4
  if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
    if (x >= 1.8 && x <2.2) {num = 0.89788; den = 0.929933; statUp = 0.0167719; statDown = 0.0168054;}
    else if (x >= 2.2 && x <2.7) {num = 0.940431; den = 0.968575; statUp = 0.0107689; statDown = 0.0109474;}
    else if (x >= 2.7 && x <3.2) {num = 0.972522; den = 0.980868; statUp = 0.0084141; statDown = 0.00875246;}
    else if (x >= 3.2 && x <3.7) {num = 0.977761; den = 0.988464; statUp = 0.00677253; statDown = 0.00722554;}
    else if (x >= 3.7 && x <4.7) {num = 0.997206; den = 0.994088; statUp = 0.00279445; statDown = 0.00494402;}
    else if (x >= 4.7 && x <8) {num = 0.996334; den = 0.99794; statUp = 0.0025871; statDown = 0.00292285;}
    else if (x >= 8 && x <11) {num = 1; den = 0.999768; statUp = 6.83212e-09; statDown = 0.00252398;}
    else if (x >= 11 && x <14) {num = 0.978986; den = 0.998297; statUp = 0.0099437; statDown = 0.0135461;}
    else {num = 1; den = 0.996328; statUp = 1.46346e-07; statDown = 0.00540197;}
  }

  if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
    // syst uncertainties
    if (x >= 3.5 && x < 4) syst = 0.000855038;
    if (x >= 4 && x < 4.5) syst = 0.000240768;
    if (x >= 4.5 && x < 5) syst = 0.000899317;
    if (x >= 5 && x < 5.5) syst = 0.00109601;
    if (x >= 5.5 && x < 6.5) syst = 0.000311285;
    if (x >= 6.5 && x < 8) syst = 0.000476892;
    if (x >= 8 && x < 10.5) syst = 0.00122726;
    if (x >= 10.5 && x < 14) syst = 0.00274564;
    if (x >= 14 && x < 18) syst = 0.00133458;
    if (x >= 18 && x < 30) syst = 0.00365302;
  }
  if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
    // syst uncertainties
    if (x >= 2.37 && x < 3) syst = 0.00121011;
    if (x >= 3 && x < 3.5) syst = 0.00177203;
    if (x >= 3.5 && x < 4) syst = 0.00181968;
    if (x >= 4 && x < 4.5) syst = 0.00154391;
    if (x >= 4.5 && x < 5) syst = 0.00183399;
    if (x >= 5 && x < 6) syst = 0.000750531;
    if (x >= 6 && x < 7.5) syst = 0.00111009;
    if (x >= 7.5 && x < 10) syst = 0.00123526;
    if (x >= 10 && x < 15) syst = 0.00156111;
    if (x >= 15 && x < 30) syst = 0.00311647;
  }
  if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
    // syst uncertainties
    if (x >= 1.8 && x < 2.5) syst = 0.00425782;
    if (x >= 2.5 && x < 3) syst = 0.00100318;
    if (x >= 3 && x < 3.5) syst = 0.00476823;
    if (x >= 3.5 && x < 4) syst = 0.00107436;
    if (x >= 4 && x < 4.5) syst = 0.00280821;
    if (x >= 4.5 && x < 5.5) syst = 0.000743448;
    if (x >= 5.5 && x < 7) syst = 9.14067e-07;
    if (x >= 7 && x < 9) syst = 0.00197499;
    if (x >= 9 && x < 12) syst = 0.000842135;
    if (x >= 12 && x < 20) syst = 0.000789339;
  }
  if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
    // syst uncertainties
    if (x >= 1.8 && x < 2.2) syst = 0.0234659;
    if (x >= 2.2 && x < 2.7) syst = 0.00575045;
    if (x >= 2.7 && x < 3.2) syst = 0.00231199;
    if (x >= 3.2 && x < 3.7) syst = 0.00333113;
    if (x >= 3.7 && x < 4.7) syst = 0.00183657;
    if (x >= 4.7 && x < 8) syst = 0.00226061;
    if (x >= 8 && x < 11) syst = 2.45089e-07;
    if (x >= 11 && x < 14) syst = 0.00132322;
    if (x >= 14 && x < 20) syst = 1.46346e-07;
  }
  double syst_factor = 0; double stat_factor = 0;
  if (idx == -1) syst_factor = syst;
  if (idx == -2) syst_factor = -1*syst;
  if (idx == +1) stat_factor = statUp;
  if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}

///////////////////////////////////////////////////
//              T R G     P b P b                //
///////////////////////////////////////////////////

double tnp_weight_trg_pbpb(double pt, double eta, int filterId,int idx) {
  double x = pt;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  if (filterId==0) { //L2 Jpsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.663241; den = 0.612464; statUp = 0.00756585; statDown = 0.00760942;}
      else if (x >= 4 && x <4.5) {num = 0.85917; den = 0.842014; statUp = 0.00553614; statDown = 0.0056399;}
      else if (x >= 4.5 && x <5) {num = 0.916521; den = 0.900241; statUp = 0.00469998; statDown = 0.00484238;}
      else if (x >= 5 && x <5.5) {num = 0.934641; den = 0.919918; statUp = 0.0045738; statDown = 0.00475841;}
      else if (x >= 5.5 && x <6.5) {num = 0.947297; den = 0.940729; statUp = 0.00343701; statDown = 0.00356287;}
      else if (x >= 6.5 && x <8) {num = 0.949896; den = 0.95601; statUp = 0.00340062; statDown = 0.00353618;}
      else if (x >= 8 && x <10.5) {num = 0.950506; den = 0.96243; statUp = 0.00363768; statDown = 0.00379628;}
      else if (x >= 10.5 && x <14) {num = 0.947321; den = 0.964831; statUp = 0.00478944; statDown = 0.00505914;}
      else if (x >= 14 && x <18) {num = 0.940314; den = 0.966093; statUp = 0.00760886; statDown = 0.00820262;}
      else {num = 0.943817; den = 0.957341; statUp = 0.00874496; statDown = 0.00954934;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.698036; den = 0.656794; statUp = 0.0116241; statDown = 0.0117424;}
      else if (x >= 3 && x <3.5) {num = 0.829209; den = 0.793366; statUp = 0.00781434; statDown = 0.00789364;}
      else if (x >= 3.5 && x <4) {num = 0.900878; den = 0.867213; statUp = 0.00646057; statDown = 0;}
      else if (x >= 4 && x <4.5) {num = 0.917489; den = 0.917832; statUp = 0.00630865; statDown = 0.00655031;}
      else if (x >= 4.5 && x <5) {num = 0.935786; den = 0.936854; statUp = 0.00635004; statDown = 0.00668515;}
      else if (x >= 5 && x <6) {num = 0.940635; den = 0.950382; statUp = 0.00519892; statDown = 0.00546059;}
      else if (x >= 6 && x <7.5) {num = 0.942288; den = 0.955739; statUp = 0.00542; statDown = 0.00571559;}
      else if (x >= 7.5 && x <10) {num = 0.953266; den = 0.959009; statUp = 0.00550351; statDown = 0.00587165;}
      else if (x >= 10 && x <15) {num = 0.948259; den = 0.959592; statUp = 0.00686182; statDown = 0.0074032;}
      else {num = 0.913688; den = 0.954037; statUp = 0.0126313; statDown = 0.0139534;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.641279; den = 0.556324; statUp = 0.031358; statDown = 0.0311058;}
      else if (x >= 2 && x <2.5) {num = 0.746463; den = 0.740273; statUp = 0.014056; statDown = 0.0141737;}
      else if (x >= 2.5 && x <3) {num = 0.882257; den = 0.893319; statUp = 0.0101587; statDown = 0.0104338;}
      else if (x >= 3 && x <3.5) {num = 0.908732; den = 0.927004; statUp = 0.00969244; statDown = 0.0101057;}
      else if (x >= 3.5 && x <4) {num = 0.918071; den = 0.943802; statUp = 0.00962387; statDown = 0.010125;}
      else if (x >= 4 && x <4.5) {num = 0.907847; den = 0.946961; statUp = 0.0108044; statDown = 0.0115161;}
      else if (x >= 4.5 && x <5.5) {num = 0.930721; den = 0.954231; statUp = 0.00791393; statDown = 0.00835311;}
      else if (x >= 5.5 && x <6.5) {num = 0.892041; den = 0.952605; statUp = 0.0117125; statDown = 0.0123847;}
      else if (x >= 6.5 && x <8) {num = 0.917612; den = 0.946287; statUp = 0.0198739; statDown = 0;}
      else if (x >= 8 && x <9.5) {num = 0.914016; den = 0.940696; statUp = 0.0135905; statDown = 0.0147506;}
      else if (x >= 9.5 && x <13) {num = 0.915476; den = 0.936928; statUp = 0.0139215; statDown = 0.0151525;}
      else {num = 0.905898; den = 0.924578; statUp = 0.0221738; statDown = 0.0243314;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.789011; den = 0.794511; statUp = 0.0182577; statDown = 0.018358;}
      else if (x >= 2.2 && x <2.7) {num = 0.854198; den = 0.867878; statUp = 0.0140189; statDown = 0.0142627;}
      else if (x >= 2.7 && x <3.2) {num = 0.878207; den = 0.893527; statUp = 0.0129086; statDown = 0.0132743;}
      else if (x >= 3.2 && x <3.7) {num = 0.877411; den = 0.911546; statUp = 0.0121215; statDown = 0.0125934;}
      else if (x >= 3.7 && x <4.7) {num = 0.885503; den = 0.912297; statUp = 0.0108807; statDown = 0.011314;}
      else if (x >= 4.7 && x <6.5) {num = 0.912669; den = 0.928676; statUp = 0.00953789; statDown = 0.010012;}
      else if (x >= 6.5 && x <8.5) {num = 0.900634; den = 0.922926; statUp = 0.014265; statDown = 0.0151886;}
      else if (x >= 8.5 && x <11) {num = 0.869202; den = 0.916807; statUp = 0.0198993; statDown = 0.021416;}
      else {num = 0.891181; den = 0.915654; statUp = 0.0229342; statDown = 0.0252136;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.00136539;
      if (x >= 4 && x < 4.5) syst = 0.00105145;
      if (x >= 4.5 && x < 5) syst = 0.00132265;
      if (x >= 5 && x < 5.5) syst = 0.000788367;
      if (x >= 5.5 && x < 6.5) syst = 0.000562329;
      if (x >= 6.5 && x < 8) syst = 0.000155103;
      if (x >= 8 && x < 10.5) syst = 0.000158629;
      if (x >= 10.5 && x < 14) syst = 0.000280705;
      if (x >= 14 && x < 18) syst = 0.00108834;
      if (x >= 18 && x < 30) syst = 0.00305146;
    }
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00329732;
      if (x >= 3 && x < 3.5) syst = 0.00520027;
      if (x >= 3.5 && x < 4) syst = 0.00240912;
      if (x >= 4 && x < 4.5) syst = 0.00222157;
      if (x >= 4.5 && x < 5) syst = 0.0017746;
      if (x >= 5 && x < 6) syst = 0.00167846;
      if (x >= 6 && x < 7.5) syst = 0.000384188;
      if (x >= 7.5 && x < 10) syst = 0.00184184;
      if (x >= 10 && x < 15) syst = 0.000586472;
      if (x >= 15 && x < 30) syst = 0.00495382;
    }
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.0133628;
      if (x >= 2 && x < 2.5) syst = 0.0347588;
      if (x >= 2.5 && x < 3) syst = 0.00227196;
      if (x >= 3 && x < 3.5) syst = 0.00511758;
      if (x >= 3.5 && x < 4) syst = 0.00466125;
      if (x >= 4 && x < 4.5) syst = 0.0105296;
      if (x >= 4.5 && x < 5.5) syst = 0.00242737;
      if (x >= 5.5 && x < 6.5) syst = 0.00472684;
      if (x >= 6.5 && x < 8) syst = 0.00762573;
      if (x >= 8 && x < 9.5) syst = 0.00627021;
      if (x >= 9.5 && x < 13) syst = 0.00231563;
      if (x >= 13 && x < 20) syst = 0.019257;
    }
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.0690381;
      if (x >= 2.2 && x < 2.7) syst = 0.0073301;
      if (x >= 2.7 && x < 3.2) syst = 0.00969524;
      if (x >= 3.2 && x < 3.7) syst = 0.00401152;
      if (x >= 3.7 && x < 4.7) syst = 0.0131915;
      if (x >= 4.7 && x < 6.5) syst = 0.00669816;
      if (x >= 6.5 && x < 8.5) syst = 0.00311553;
      if (x >= 8.5 && x < 11) syst = 0.00999295;
      if (x >= 11 && x < 20) syst = 0.0151298;
    }
  }
  if (filterId==1) { //L3 Jpsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.0981413; den = 0.071216; statUp = 0.00475984; statDown = 0.00464341;}
      else if (x >= 4 && x <4.5) {num = 0.309366; den = 0.234954; statUp = 0.00730881; statDown = 0.00724846;}
      else if (x >= 4.5 && x <5) {num = 0.49696; den = 0.427459; statUp = 0.00850388; statDown = 0.00850324;}
      else if (x >= 5 && x <5.5) {num = 0.646567; den = 0.569803; statUp = 0.00897182; statDown = 0.00902901;}
      else if (x >= 5.5 && x <6.5) {num = 0.717727; den = 0.665688; statUp = 0.00690496; statDown = 0.00696443;}
      else if (x >= 6.5 && x <8) {num = 0.770814; den = 0.736851; statUp = 0.00662338; statDown = 0.00670205;}
      else if (x >= 8 && x <10.5) {num = 0.791509; den = 0.777511; statUp = 0.0068552; statDown = 0.00695591;}
      else if (x >= 10.5 && x <14) {num = 0.826768; den = 0.814073; statUp = 0.00832433; statDown = 0.00852041;}
      else if (x >= 14 && x <18) {num = 0.799276; den = 0.820805; statUp = 0.0131407; statDown = 0.0135476;}
      else {num = 0.844468; den = 0.83657; statUp = 0.0140064; statDown = 0.0146299;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.331992; den = 0.302908; statUp = 0.011463; statDown = 0.0112805;}
      else if (x >= 3 && x <3.5) {num = 0.429327; den = 0.424718; statUp = 0.00993062; statDown = 0.00987992;}
      else if (x >= 3.5 && x <4) {num = 0.543359; den = 0.527324; statUp = 0.00926375; statDown = 0.00951398;}
      else if (x >= 4 && x <4.5) {num = 0.590264; den = 0.603494; statUp = 0.0113814; statDown = 0.0114274;}
      else if (x >= 4.5 && x <5) {num = 0.638602; den = 0.645253; statUp = 0.0126665; statDown = 0.0127591;}
      else if (x >= 5 && x <6) {num = 0.671782; den = 0.678571; statUp = 0.0104761; statDown = 0.0105625;}
      else if (x >= 6 && x <7.5) {num = 0.683073; den = 0.719277; statUp = 0.0110239; statDown = 0.0111476;}
      else if (x >= 7.5 && x <10) {num = 0.746848; den = 0.753451; statUp = 0.0114305; statDown = 0.0116301;}
      else if (x >= 10 && x <15) {num = 0.772519; den = 0.800301; statUp = 0.0133498; statDown = 0.013682;}
      else {num = 0.768513; den = 0.831634; statUp = 0.0213549; statDown = 0.0221413;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.0302057; den = 0.0291672; statUp = 0.0103319; statDown = 0.00934636;}
      else if (x >= 2 && x <2.5) {num = 0.173305; den = 0.159501; statUp = 0.0110778; statDown = 0.0107941;}
      else if (x >= 2.5 && x <3) {num = 0.470274; den = 0.464604; statUp = 0.0151192; statDown = 0.0150115;}
      else if (x >= 3 && x <3.5) {num = 0.661007; den = 0.654318; statUp = 0.0158638; statDown = 0.0159609;}
      else if (x >= 3.5 && x <4) {num = 0.707769; den = 0.711729; statUp = 0.0158855; statDown = 0.0161;}
      else if (x >= 4 && x <4.5) {num = 0.726282; den = 0.741631; statUp = 0.0174238; statDown = 0.017772;}
      else if (x >= 4.5 && x <5.5) {num = 0.763392; den = 0.796057; statUp = 0.0131827; statDown = 0.013438;}
      else if (x >= 5.5 && x <6.5) {num = 0.767094; den = 0.82346; statUp = 0.016378; statDown = 0.0168236;}
      else if (x >= 6.5 && x <8) {num = 0.810425; den = 0.833902; statUp = 0.0145812; statDown = 0.0150811;}
      else if (x >= 8 && x <9.5) {num = 0.819571; den = 0.840655; statUp = 0.0190366; statDown = 0.01987;}
      else if (x >= 9.5 && x <13) {num = 0.829677; den = 0.856559; statUp = 0.0192692; statDown = 0.0202277;}
      else {num = 0.862606; den = 0.873057; statUp = 0.0264458; statDown = 0.0282459;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.0980329; den = 0.092124; statUp = 0.0109739; statDown = 0.0105493;}
      else if (x >= 2.2 && x <2.7) {num = 0.264092; den = 0.248237; statUp = 0.0154033; statDown = 0.0150975;}
      else if (x >= 2.7 && x <3.2) {num = 0.417324; den = 0.407684; statUp = 0.0182234; statDown = 0.0179965;}
      else if (x >= 3.2 && x <3.7) {num = 0.486753; den = 0.49461; statUp = 0.0184525; statDown = 0.0183476;}
      else if (x >= 3.7 && x <4.7) {num = 0.569753; den = 0.574351; statUp = 0.0169356; statDown = 0.0169599;}
      else if (x >= 4.7 && x <6.5) {num = 0.674727; den = 0.651274; statUp = 0.0161099; statDown = 0.0162599;}
      else if (x >= 6.5 && x <8.5) {num = 0.720956; den = 0.710909; statUp = 0.0216217; statDown = 0.0220983;}
      else if (x >= 8.5 && x <11) {num = 0.714358; den = 0.74928; statUp = 0.0275205; statDown = 0.0283679;}
      else {num = 0.748962; den = 0.778156; statUp = 0.0324675; statDown = 0.0339396;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.000650706;
      if (x >= 4 && x < 4.5) syst = 0.00108325;
      if (x >= 4.5 && x < 5) syst = 0.00298052;
      if (x >= 5 && x < 5.5) syst = 0.00341277;
      if (x >= 5.5 && x < 6.5) syst = 0.000613358;
      if (x >= 6.5 && x < 8) syst = 0.000661917;
      if (x >= 8 && x < 10.5) syst = 0.000756817;
      if (x >= 10.5 && x < 14) syst = 0.000783331;
      if (x >= 14 && x < 18) syst = 0.00218267;
      if (x >= 18 && x < 30) syst = 0.00267133;
    }
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00235318;
      if (x >= 3 && x < 3.5) syst = 0.0041557;
      if (x >= 3.5 && x < 4) syst = 0.00434833;
      if (x >= 4 && x < 4.5) syst = 0.00151395;
      if (x >= 4.5 && x < 5) syst = 0.00325726;
      if (x >= 5 && x < 6) syst = 0.00151752;
      if (x >= 6 && x < 7.5) syst = 0.0019946;
      if (x >= 7.5 && x < 10) syst = 0.00321667;
      if (x >= 10 && x < 15) syst = 0.00350192;
      if (x >= 15 && x < 30) syst = 0.00441023;
    }
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.0108482;
      if (x >= 2 && x < 2.5) syst = 0.00318395;
      if (x >= 2.5 && x < 3) syst = 0.00449944;
      if (x >= 3 && x < 3.5) syst = 0.00272957;
      if (x >= 3.5 && x < 4) syst = 0.00779276;
      if (x >= 4 && x < 4.5) syst = 0.00622458;
      if (x >= 4.5 && x < 5.5) syst = 0.00114197;
      if (x >= 5.5 && x < 6.5) syst = 0.0077565;
      if (x >= 6.5 && x < 8) syst = 0.00157525;
      if (x >= 8 && x < 9.5) syst = 0.00373986;
      if (x >= 9.5 && x < 13) syst = 0.00445251;
      if (x >= 13 && x < 20) syst = 0.0284466;
    }
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.0294959;
      if (x >= 2.2 && x < 2.7) syst = 0.00108294;
      if (x >= 2.7 && x < 3.2) syst = 0.0152933;
      if (x >= 3.2 && x < 3.7) syst = 0.00300734;
      if (x >= 3.7 && x < 4.7) syst = 0.00349077;
      if (x >= 4.7 && x < 6.5) syst = 0.0159082;
      if (x >= 6.5 && x < 8.5) syst = 0.00321438;
      if (x >= 8.5 && x < 11) syst = 0.00306389;
      if (x >= 11 && x < 20) syst = 0.0125999;
    }
  }
  if (filterId==2) { //L2 Upsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.687436; den = 0.628026; statUp = 0.00741958; statDown = 0.00747367;}
      else if (x >= 4 && x <4.5) {num = 0.879668; den = 0.855425; statUp = 0.00516048; statDown = 0.0052745;}
      else if (x >= 4.5 && x <5) {num = 0.936813; den = 0.91093; statUp = 0.00414716; statDown = 0.00429608;}
      else if (x >= 5 && x <5.5) {num = 0.948744; den = 0.930006; statUp = 0.00405143; statDown = 0.00424179;}
      else if (x >= 5.5 && x <6.5) {num = 0.964058; den = 0.95091; statUp = 0.00286448; statDown = 0.00299563;}
      else if (x >= 6.5 && x <8) {num = 0.967607; den = 0.967695; statUp = 0.0027257; statDown = 0.00286633;}
      else if (x >= 8 && x <10.5) {num = 0.977261; den = 0.975624; statUp = 0.0024585; statDown = 0.00262625;}
      else if (x >= 10.5 && x <14) {num = 0.975199; den = 0.979933; statUp = 0.00326164; statDown = 0.00354865;}
      else if (x >= 14 && x <18) {num = 0.978416; den = 0.983164; statUp = 0.00450526; statDown = 0.00515026;}
      else {num = 0.986566; den = 0.984365; statUp = 0.00413127; statDown = 0.00501923;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.707798; den = 0.663284; statUp = 0.0115376; statDown = 0.0116219;}
      else if (x >= 3 && x <3.5) {num = 0.841787; den = 0.79966; statUp = 0.00756481; statDown = 0.00765363;}
      else if (x >= 3.5 && x <4) {num = 0.915148; den = 0.874984; statUp = 0.0058153; statDown = 0;}
      else if (x >= 4 && x <4.5) {num = 0.932092; den = 0.925267; statUp = 0.00578591; statDown = 0.0060332;}
      else if (x >= 4.5 && x <5) {num = 0.95135; den = 0.943497; statUp = 0.005532; statDown = 0.00587852;}
      else if (x >= 5 && x <6) {num = 0.968164; den = 0.958147; statUp = 0.00392457; statDown = 0.00419381;}
      else if (x >= 6 && x <7.5) {num = 0.962352; den = 0.96425; statUp = 0.00435955; statDown = 0.00466282;}
      else if (x >= 7.5 && x <10) {num = 0.979864; den = 0.969778; statUp = 0.00356325; statDown = 0.00395684;}
      else if (x >= 10 && x <15) {num = 0.983282; den = 0.973169; statUp = 0.00384401; statDown = 0.00442388;}
      else {num = 0.985422; den = 0.979848; statUp = 0.00507078; statDown = 0.00689321;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.651185; den = 0.562127; statUp = 0.0311006; statDown = 0.0309056;}
      else if (x >= 2 && x <2.5) {num = 0.75568; den = 0.744834; statUp = 0.0139774; statDown = 0.0139217;}
      else if (x >= 2.5 && x <3) {num = 0.895072; den = 0.898086; statUp = 0.00967712; statDown = 0.00997311;}
      else if (x >= 3 && x <3.5) {num = 0.911965; den = 0.932775; statUp = 0.0093472; statDown = 0.00978825;}
      else if (x >= 3.5 && x <4) {num = 0.933577; den = 0.949087; statUp = 0.008866; statDown = 0.00937714;}
      else if (x >= 4 && x <4.5) {num = 0.92098; den = 0.952348; statUp = 0.0100317; statDown = 0.0107328;}
      else if (x >= 4.5 && x <5.5) {num = 0.945015; den = 0.959924; statUp = 0.00698531; statDown = 0.00744697;}
      else if (x >= 5.5 && x <6.5) {num = 0.915225; den = 0.961517; statUp = 0.0105111; statDown = 0.0112266;}
      else if (x >= 6.5 && x <8) {num = 0.942614; den = 0.958878; statUp = 0.00845552; statDown = 0.00920477;}
      else if (x >= 8 && x <9.5) {num = 0.943992; den = 0.95865; statUp = 0.0111118; statDown = 0.0123342;}
      else if (x >= 9.5 && x <13) {num = 0.953176; den = 0.956512; statUp = 0.0100647; statDown = 0.0114277;}
      else {num = 0.9711; den = 0.959808; statUp = 0.0112001; statDown = 0.0137596;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.80688; den = 0.804533; statUp = 0.017851; statDown = 0.0179625;}
      else if (x >= 2.2 && x <2.7) {num = 0.864384; den = 0.87583; statUp = 0.0134953; statDown = 0.0137959;}
      else if (x >= 2.7 && x <3.2) {num = 0.892179; den = 0.904052; statUp = 0.0121752; statDown = 0.0125694;}
      else if (x >= 3.2 && x <3.7) {num = 0.896585; den = 0.926574; statUp = 0.0112721; statDown = 0.0117765;}
      else if (x >= 3.7 && x <4.7) {num = 0.910998; den = 0.928332; statUp = 0.00974335; statDown = 0.010226;}
      else if (x >= 4.7 && x <6.5) {num = 0.938687; den = 0.949857; statUp = 0.00795089; statDown = 0.00845509;}
      else if (x >= 6.5 && x <8.5) {num = 0.941825; den = 0.954442; statUp = 0.0109965; statDown = 0.0143609;}
      else if (x >= 8.5 && x <11) {num = 0.931535; den = 0.959109; statUp = 0.0136929; statDown = 0.0156881;}
      else {num = 0.969056; den = 0.967084; statUp = 0.00977435; statDown = 0.0126053;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.00119995;
      if (x >= 4 && x < 4.5) syst = 0.000801484;
      if (x >= 4.5 && x < 5) syst = 0.00142786;
      if (x >= 5 && x < 5.5) syst = 0.000859141;
      if (x >= 5.5 && x < 6.5) syst = 0.000855793;
      if (x >= 6.5 && x < 8) syst = 0.000338442;
      if (x >= 8 && x < 10.5) syst = 0.000905661;
      if (x >= 10.5 && x < 14) syst = 0.000193737;
      if (x >= 14 && x < 18) syst = 0.000621028;
      if (x >= 18 && x < 30) syst = 0.0029276;
    }
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00301699;
      if (x >= 3 && x < 3.5) syst = 0.0051637;
      if (x >= 3.5 && x < 4) syst = 0.00271564;
      if (x >= 4 && x < 4.5) syst = 0.00128082;
      if (x >= 4.5 && x < 5) syst = 0.00105614;
      if (x >= 5 && x < 6) syst = 0.00120191;
      if (x >= 6 && x < 7.5) syst = 0.000729975;
      if (x >= 7.5 && x < 10) syst = 0.00139352;
      if (x >= 10 && x < 15) syst = 0.00151879;
      if (x >= 15 && x < 30) syst = 0.00138277;
    }
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.0234362;
      if (x >= 2 && x < 2.5) syst = 0.00781699;
      if (x >= 2.5 && x < 3) syst = 0.0020642;
      if (x >= 3 && x < 3.5) syst = 0.00494294;
      if (x >= 3.5 && x < 4) syst = 0.00372959;
      if (x >= 4 && x < 4.5) syst = 0.0101533;
      if (x >= 4.5 && x < 5.5) syst = 0.00248577;
      if (x >= 5.5 && x < 6.5) syst = 0.00480156;
      if (x >= 6.5 && x < 8) syst = 0.00535204;
      if (x >= 8 && x < 9.5) syst = 0.00407749;
      if (x >= 9.5 && x < 13) syst = 0.000734987;
      if (x >= 13 && x < 20) syst = 0.0152108;
    }
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.0547508;
      if (x >= 2.2 && x < 2.7) syst = 0.0439035;
      if (x >= 2.7 && x < 3.2) syst = 0.0100721;
      if (x >= 3.2 && x < 3.7) syst = 0.00486924;
      if (x >= 3.7 && x < 4.7) syst = 0.0164241;
      if (x >= 4.7 && x < 6.5) syst = 0.0045128;
      if (x >= 6.5 && x < 8.5) syst = 0.00615735;
      if (x >= 8.5 && x < 11) syst = 0.00521994;
      if (x >= 11 && x < 20) syst = 0.00496602;
    }
  }
  if (filterId==3) { //L3 Upsi
    // SF for 0 < |eta| < 1.2
    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      if (x >= 3.5 && x <4) {num = 0.0981413; den = 0.0714076; statUp = 0.00475984; statDown = 0.00464341;}
      else if (x >= 4 && x <4.5) {num = 0.309591; den = 0.234967; statUp = 0.00731017; statDown = 0.00724988;}
      else if (x >= 4.5 && x <5) {num = 0.49696; den = 0.427491; statUp = 0.00850388; statDown = 0.00850324;}
      else if (x >= 5 && x <5.5) {num = 0.646567; den = 0.569805; statUp = 0.00897182; statDown = 0.00902901;}
      else if (x >= 5.5 && x <6.5) {num = 0.717727; den = 0.665698; statUp = 0.00690496; statDown = 0.00696443;}
      else if (x >= 6.5 && x <8) {num = 0.771046; den = 0.736859; statUp = 0.00662127; statDown = 0.00670002;}
      else if (x >= 8 && x <10.5) {num = 0.792067; den = 0.777534; statUp = 0.00684886; statDown = 0.00695048;}
      else if (x >= 10.5 && x <14) {num = 0.826589; den = 0.814236; statUp = 0.00832558; statDown = 0.00852162;}
      else if (x >= 14 && x <18) {num = 0.800339; den = 0.820918; statUp = 0.0131166; statDown = 0.0135246;}
      else {num = 0.846856; den = 0.837225; statUp = 0.0139208; statDown = 0.0145458;}
    }
    // SF for 1.2 < |eta| < 1.8
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      if (x >= 2.37 && x <3) {num = 0.307823; den = 0.284114; statUp = 0.0111767; statDown = 0.0110334;}
      else if (x >= 3 && x <3.5) {num = 0.429139; den = 0.424849; statUp = 0.00992916; statDown = 0.00988034;}
      else if (x >= 3.5 && x <4) {num = 0.54449; den = 0.527662; statUp = 0.45551; statDown = 0.00969908;}
      else if (x >= 4 && x <4.5) {num = 0.591156; den = 0.604174; statUp = 0.0113833; statDown = 0.0114203;}
      else if (x >= 4.5 && x <5) {num = 0.639967; den = 0.645913; statUp = 0.0126515; statDown = 0.0127542;}
      else if (x >= 5 && x <6) {num = 0.673449; den = 0.679156; statUp = 0.0104608; statDown = 0.010556;}
      else if (x >= 6 && x <7.5) {num = 0.685635; den = 0.720417; statUp = 0.0109999; statDown = 0.0111263;}
      else if (x >= 7.5 && x <10) {num = 0.749011; den = 0.754465; statUp = 0.0113948; statDown = 0.0115976;}
      else if (x >= 10 && x <15) {num = 0.773253; den = 0.801728; statUp = 0.0133331; statDown = 0.0136666;}
      else {num = 0.773038; den = 0.833024; statUp = 0.0211603; statDown = 0.0220621;}
    }
    // SF for 1.8 < |eta| < 2.1
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      if (x >= 1.8 && x <2) {num = 0.00155461; den = 0.000463755; statUp = 0.00211758; statDown = 0.00108585;}
      else if (x >= 2 && x <2.5) {num = 0.00555387; den = 0.00858885; statUp = 0.00315843; statDown = 0.00555387;}
      else if (x >= 2.5 && x <3) {num = 0.451703; den = 0.447277; statUp = 0.0150343; statDown = 0.0149314;}
      else if (x >= 3 && x <3.5) {num = 0.655861; den = 0.650893; statUp = 0.0159254; statDown = 0.0160134;}
      else if (x >= 3.5 && x <4) {num = 0.706533; den = 0.710335; statUp = 0.0158931; statDown = 0.0161213;}
      else if (x >= 4 && x <4.5) {num = 0.726282; den = 0.741482; statUp = 0.0174238; statDown = 0.017772;}
      else if (x >= 4.5 && x <5.5) {num = 0.764391; den = 0.796199; statUp = 0.013155; statDown = 0.0134088;}
      else if (x >= 5.5 && x <6.5) {num = 0.769821; den = 0.824468; statUp = 0.016319; statDown = 0.0167672;}
      else if (x >= 6.5 && x <8) {num = 0.811763; den = 0.834174; statUp = 0.0145482; statDown = 0.0150501;}
      else if (x >= 8 && x <9.5) {num = 0.819571; den = 0.841319; statUp = 0.0190366; statDown = 0.01987;}
      else if (x >= 9.5 && x <13) {num = 0.829677; den = 0.857122; statUp = 0.0192692; statDown = 0.0202277;}
      else {num = 0.874981; den = 0.874474; statUp = 0.0258047; statDown = 0.0277093;}
    }
    // SF for 2.1 < |eta| < 2.4
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      if (x >= 1.8 && x <2.2) {num = 0.00366548; den = 0.000413051; statUp = 0.00211949; statDown = 0.00167583;}
      else if (x >= 2.2 && x <2.7) {num = 0.116176; den = 0.109916; statUp = 0.0108872; statDown = 0.0105043;}
      else if (x >= 2.7 && x <3.2) {num = 0.413123; den = 0.401827; statUp = 0.0181733; statDown = 0.017947;}
      else if (x >= 3.2 && x <3.7) {num = 0.482035; den = 0.491334; statUp = 0.0184417; statDown = 0.0183373;}
      else if (x >= 3.7 && x <4.7) {num = 0.568894; den = 0.573223; statUp = 0.0169342; statDown = 0.0169658;}
      else if (x >= 4.7 && x <6.5) {num = 0.675048; den = 0.651616; statUp = 0.0160907; statDown = 0.0162558;}
      else if (x >= 6.5 && x <8.5) {num = 0.722882; den = 0.711847; statUp = 0.0215971; statDown = 0.022054;}
      else if (x >= 8.5 && x <11) {num = 0.714358; den = 0.750096; statUp = 0.0275205; statDown = 0.0283679;}
      else {num = 0.753355; den = 0.779864; statUp = 0.0323455; statDown = 0.0336041;}
    }

    if (fabs(eta) >= 0 && fabs(eta) < 1.2) {
      // syst uncertainties
      if (x >= 3.5 && x < 4) syst = 0.000650706;
      if (x >= 4 && x < 4.5) syst = 0.0010869;
      if (x >= 4.5 && x < 5) syst = 0.00298052;
      if (x >= 5 && x < 5.5) syst = 0.00341277;
      if (x >= 5.5 && x < 6.5) syst = 0.000613358;
      if (x >= 6.5 && x < 8) syst = 0.000658119;
      if (x >= 8 && x < 10.5) syst = 0.000756756;
      if (x >= 10.5 && x < 14) syst = 0.000662617;
      if (x >= 14 && x < 18) syst = 0.00220571;
      if (x >= 18 && x < 30) syst = 0.00215326;
    }
    if (fabs(eta) >= 1.2 && fabs(eta) < 1.8) {
      // syst uncertainties
      if (x >= 2.37 && x < 3) syst = 0.00406324;
      if (x >= 3 && x < 3.5) syst = 0.00422745;
      if (x >= 3.5 && x < 4) syst = 0.00493964;
      if (x >= 4 && x < 4.5) syst = 0.0015019;
      if (x >= 4.5 && x < 5) syst = 0.00349953;
      if (x >= 5 && x < 6) syst = 0.00165421;
      if (x >= 6 && x < 7.5) syst = 0.00195686;
      if (x >= 7.5 && x < 10) syst = 0.00305233;
      if (x >= 10 && x < 15) syst = 0.00341103;
      if (x >= 15 && x < 30) syst = 0.00425449;
    }
    if (fabs(eta) >= 1.8 && fabs(eta) < 2.1) {
      // syst uncertainties
      if (x >= 1.8 && x < 2) syst = 0.00154707;
      if (x >= 2 && x < 2.5) syst = 0.00239348;
      if (x >= 2.5 && x < 3) syst = 0.00578521;
      if (x >= 3 && x < 3.5) syst = 0.0019864;
      if (x >= 3.5 && x < 4) syst = 0.00826595;
      if (x >= 4 && x < 4.5) syst = 0.00622458;
      if (x >= 4.5 && x < 5.5) syst = 0.00155048;
      if (x >= 5.5 && x < 6.5) syst = 0.00738518;
      if (x >= 6.5 && x < 8) syst = 0.00155169;
      if (x >= 8 && x < 9.5) syst = 0.00373986;
      if (x >= 9.5 && x < 13) syst = 0.00445251;
      if (x >= 13 && x < 20) syst = 0.028681;
    }
    if (fabs(eta) >= 2.1 && fabs(eta) < 2.4) {
      // syst uncertainties
      if (x >= 1.8 && x < 2.2) syst = 0.000519911;
      if (x >= 2.2 && x < 2.7) syst = 0.0288676;
      if (x >= 2.7 && x < 3.2) syst = 0.013137;
      if (x >= 3.2 && x < 3.7) syst = 0.00153582;
      if (x >= 3.7 && x < 4.7) syst = 0.00361851;
      if (x >= 4.7 && x < 6.5) syst = 0.0155374;
      if (x >= 6.5 && x < 8.5) syst = 0.00321391;
      if (x >= 8.5 && x < 11) syst = 0.00306389;
      if (x >= 11 && x < 20) syst = 0.0199929;
    }
  }
  double syst_factor = 0; double stat_factor = 0;
  if (idx == -1) syst_factor = syst;
  if (idx == -2) syst_factor = -1*syst;
  if (idx == +1) stat_factor = statUp;
  if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}

///////////////////////////////////////////////////
//              T R K     P b P b                //
///////////////////////////////////////////////////

double tnp_weight_trk_pbpb(double eta, int idx) {
  double x = eta;
  double num=1, den=1, syst=0, statUp=0, statDown=0;
  //SF in eta bins
  if (x >= -2.4 && x < -1.6) {num = 0.994498; den = 0.998413; statUp = 0.00287386; statDown = 0.00290888;}
  if (x >= -1.6 && x < -1.2) {num = 0.973539; den = 0.967322; statUp = 0.00399501; statDown = 0.00407843;}
  if (x >= -1.2 && x < -0.9) {num = 0.964465; den = 0.970816; statUp = 0.00861188; statDown = 0.00869453;}
  if (x >= -0.9 && x < -0.6) {num = 0.96081; den = 0.974407; statUp = 0.0223599; statDown = 0.00682405;}
  if (x >= -0.6 && x < -0.3) {num = 0.964464; den = 0.97802; statUp = 0.00612474; statDown = 0.00621291;}
  if (x >= -0.3 && x < 0.3) {num = 0.963862; den = 0.966583; statUp = 0.00496914; statDown = 0.00503378;}
  if (x >= 0.3 && x < 0.6) {num = 0.956897; den = 0.967248; statUp = 0.00672757; statDown = 0.0068202;}
  if (x >= 0.6 && x < 0.9) {num = 0.964172; den = 0.966882; statUp = 0.00735892; statDown = 0.00754429;}
  if (x >= 0.9 && x < 1.2) {num = 0.961874; den = 0.955987; statUp = 0.00987473; statDown = 0.0099638;}
  if (x >= 1.2 && x < 1.6) {num = 0.964754; den = 0.964653; statUp = 0.0042287; statDown = 0.00430601;}
  if (x >= 1.6 && x < 2.4) {num = 0.999937; den = 0.998771; statUp = 6.33084e-05; statDown = 0.00310832;}

  // syst uncertainties
  if (x >= -2.4 && x < -1.6) syst = 0.00129036;
  if (x >= -1.6 && x < -1.2) syst = 0.00376932;
  if (x >= -1.2 && x < -0.9) syst = 0.00125496;
  if (x >= -0.9 && x < -0.6) syst = 0.00240006;
  if (x >= -0.6 && x < -0.3) syst = 0.00228604;
  if (x >= -0.3 && x < 0.3) syst = 0.00494002;
  if (x >= 0.3 && x < 0.6) syst = 0.00529974;
  if (x >= 0.6 && x < 0.9) syst = 0.0023114;
  if (x >= 0.9 && x < 1.2) syst = 0.012059;
  if (x >= 1.2 && x < 1.6) syst = 0.00278996;
  if (x >= 1.6 && x < 2.4) syst = 0.00104774;

  double syst_factor = 0; double stat_factor = 0;
  if (idx == -1) syst_factor = syst;
  if (idx == -2) syst_factor = -1*syst;
  if (idx == +1) stat_factor = statUp;
  if (idx == +2) stat_factor = -1*statDown;
  return ((num+syst_factor+stat_factor)/den);
}

#endif
