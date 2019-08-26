#ifndef QN
#define QN
#ifndef INFILE
#define INFILE
string infile = "data.root";
#endif
#ifndef DEBUG
#define DEBUG
bool debug = false;
#endif
#ifndef GRAPH
#define GRAPH
bool graph = true;
#endif
#include "HiEvtPlaneList.h"
using namespace hi;
using namespace std;
const  int ncent = 13;
const  int centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,100};
string getNextCent(int mincent,int maxcent, int & ccount){
  int imin = -1;
  int imax = -1;
  for(int i = 0; i<ncent; i++) {
    if(mincent==centbins[i]) imin=i;
    if(maxcent==centbins[i]) imax=i;
  }
  string cbin = to_string(centbins[imin+ccount])+"_"+to_string(centbins[imin+ccount+1]);
  if(imin+ccount+1>imax) cbin="stop";
  ++ccount;
  return cbin;
}
#include "epselection.C"
#include "qAB.C"
TH1D * qn(string anal="N2", string side="+",  int mincent = 15, int maxcent = 20, double etamin = -0.8, double etamax = 0.8, bool threesub = true, bool decor = false, int decorOff=0,int marker = 20, int color = kBlue){
  if(debug) {
    cout<<"Input File: "<<infile<<"   anal: "<<anal<<endl;
    if(decor) {
      cout<<"decor true"<<endl;
    } else {
      cout<<"decor false"<<endl;
    } 
  }
  int A = -1;
  int B = -1;
  int C = -1;
  int xref = -1;
  string harmonicAnal = anal; 
  string harmonicSide = "A";
  if(side.find("-")!=std::string::npos) side = "B";
  epselection(anal,side,etamin,true,decor,A,B,C);
  double AB = qAB(A,B,mincent,maxcent);
  double AC = qAB(A,C,mincent,maxcent);
  double BC = qAB(B,C,mincent,maxcent);
  double rAB = resAB(A,B,mincent,maxcent);
  double rAC = resAB(A,C,mincent,maxcent);
  double rBC = resAB(B,C,mincent,maxcent);
  double res = sqrt(fabs(AB));
  double epres = sqrt(fabs(rAB));
  if(threesub) res = sqrt(fabs(AB)*fabs(AC)/fabs(BC));
  if(threesub) epres = sqrt(fabs(rAB)*fabs(rAC)/fabs(rBC));
  if(anal.find("N112")!=std::string::npos) {
    AB = 1;
    AC = 1;
    BC = 1;
    rAB = 1;
    rAC = 1;
    rBC = 1;
  }
  if(debug) cout<<etamin<<"\t"<<etamax<<"\t"<<EPNames[A]<<"\t"<<EPNames[B]<<"\t"<<EPNames[C]<<"\tEPres: "<<epres<<endl;
  TFile * fin = new TFile(infile.data());
  string tmp = Form("vnanalyzer/Harmonics/%d_%d/%s/q%s",mincent,maxcent,harmonicAnal.data(),harmonicSide.data());
  cout<<tmp<<endl;
  TH2D * qn = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%d_%d/%s/q%s",mincent,maxcent,harmonicAnal.data(),harmonicSide.data()));
  TH2D * wn = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%d_%d/%s/wn%s",mincent,maxcent,harmonicAnal.data(),harmonicSide.data()));
  cout<<infile<<"\t"<<fin<<"\t"<<qn<<"\t"<<wn<<endl;
  int ietamin = qn->GetYaxis()->FindBin(etamin+0.01);
  int ietamax = qn->GetYaxis()->FindBin(etamax-0.01);
  TH1D * qn1d = (TH1D *) qn->ProjectionX(Form("%s_%s_%d_%d_clone",harmonicAnal.data(),harmonicSide.data(),mincent,maxcent), ietamin,ietamax);
  qn1d->Divide((TH1D *) wn->ProjectionX(Form("%s_%s_%d_%d",harmonicAnal.data(),harmonicSide.data(),mincent,maxcent),ietamin,ietamax));
  qn1d->Scale(1/res);
  qn1d->SetDirectory(0);
  qn1d->SetMarkerStyle(marker);
  qn1d->SetMarkerColor(color);
  qn1d->SetLineColor(color);
  fin->Close();
  if(graph) qn1d->Draw();
  return qn1d;
}

#endif
