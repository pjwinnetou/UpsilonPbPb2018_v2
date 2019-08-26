#ifndef QAB
#define QAB
#ifndef INFILE
#define INFILE
string infile = "data.root";
#endif
#ifndef DEBUG
#define DEBUG
bool debug = false;
#endif
#include "HiEvtPlaneList.h"
using namespace hi;
using namespace std;
int fbeg[4] = {HFm1,HFm2,HFm3,HFm4};
double qAB(int Aname=HFm2, int Bname=HFp2,  int mincent=15, int maxcent=20){
  TFile * fin = new TFile(infile.data());
  TH2D *res;
  TH2D *resw;
  if(EPOrder[Aname]==EPOrder[Bname]) {
    res = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%d_%d/res%d",mincent,maxcent,EPOrder[Aname]));
    resw = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%d_%d/resw%d",mincent,maxcent,EPOrder[Aname]));
  }  else if(EPOrder[Aname]==1 && EPOrder[Bname]==2) {
    res = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%d_%d/res12",mincent,maxcent));
    resw = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%d_%d/resw12",mincent,maxcent));
  }
  res->Divide(resw);
  int indx1 = Aname-fbeg[EPOrder[Aname]-1];
  int indx2 = Bname-fbeg[EPOrder[Bname]-1];
  if(indx2<indx1) {
    int hold = indx1;
    indx1 = indx2;
    indx2 = hold;
  }

  double ret = res->GetBinContent(indx1+1,indx2+1);
  cout<<"ret = "<<ret<<endl;
  fin->Close();
  return ret;

}
double resAB(int Aname=HFm2, int Bname=HFp2,  int mincent=15, int maxcent=20){
  TFile * fin = new TFile(infile.data());
  TH2D *res;
  TH2D *resw;
  if(EPOrder[Aname]==EPOrder[Bname]) {
    res = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%d_%d/resep%d",mincent,maxcent,EPOrder[Aname]));
    resw = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%d_%d/rescnt%d",mincent,maxcent,EPOrder[Aname]));
  }  else if(EPOrder[Aname]==1 && EPOrder[Bname]==2) {
    res = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%d_%d/resep12",mincent,maxcent));
    resw = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%d_%d/rescnt12",mincent,maxcent));
  }
  res->Divide(resw);
  int indx1 = Aname-fbeg[EPOrder[Aname]-1];
  int indx2 = Bname-fbeg[EPOrder[Bname]-1];

  if(indx2<indx1) {
    int hold = indx1;
    indx1 = indx2;
    indx2 = hold;
  }
  double ret = res->GetBinContent(indx1+1,indx2+1);
  if(debug) cout<<EPNames[Aname]<<"/"<<EPNames[Bname]<<" EP Resolution "<<mincent<<"\t"<<maxcent<<"\t"<<ret<<endl;
  fin->Close();
  return ret;

}
#endif
