bool epselection(string anal, string side, double etamin, bool threesub, bool decor, int & A, int & B, int & C){
  bool isgood = false;
  int ieta = (int) ( (etamin+2.4)/0.4 + 0.01 );
  A = -1;
  B = -1;
  C = -1;
  cout<<"epselection: "<<anal<<"\t"<<side<<"\t"<<ieta<<endl;

  //+++ N2
  if(anal.find("N2") != std::string::npos) {
    if(side.find("+") != std::string::npos) {
      A = HFp2;
      B = HFm2;
      C = trackmid2;
    } else {
      A = HFm2;
      B = HFp2;
      C = trackmid2;
    }
  }
  //+++ N3
  if(anal.find("N3") != std::string::npos) {
    if(side.find("+") != std::string::npos) {
      A = HFp3;
      B = HFm3;
      C = trackmid3;
    } else {
      A = HFm3;
      B = HFp3;
      C = trackmid3;
    }
  }
  //+++ N4
  if(anal.find("N4") != std::string::npos) {
    if(side.find("+") != std::string::npos) {
      A = HFp4;
      B = HFm4;
      C = trackmid4;
    } else {
      A = HFm4;
      B = HFp4;
      C = trackmid4;
    }
  }
  cout<<"A: "<<A<<endl;
  cout<<"B: "<<B<<endl;
  cout<<"C: "<<C<<endl;

  if(A>0 && B>0 && C>0) isgood = true;
  return isgood;
}
