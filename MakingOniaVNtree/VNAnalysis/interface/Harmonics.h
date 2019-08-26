void Fill_N(int anal,double order, int bin, TH2D *qxtrk_, TH2D * qytrk_, TH2D * qcnt_, double Ax, double Ay, double Bx, double By, double Cx, double Cy, double wA, double wB, double wC){
  if(pow(Ax,2)+pow(Ay,2) < 1e-6) return;
  if(pow(Bx,2)+pow(By,2) < 1e-6) return;
  if(pow(Cx,2)+pow(Cy,2) < 1e-6) return;
  
  qanal[anal].thA[bin][0]->Fill(TMath::ATan2(Ay,Ax)/order);
  qanal[anal].thB[bin][0]->Fill(TMath::ATan2(By,Bx)/order);
  qanal[anal].thC[bin][0]->Fill(TMath::ATan2(Cy,Cx)/order);
  qanal[anal].thA2[bin][0]->Fill(0);
  qanal[anal].thB2[bin][0]->Fill(0);
  qanal[anal].thC2[bin][0]->Fill(0);
  qanal[anal].thn[bin][0]->Fill(TMath::ATan2(qytrk_->GetBinContent(5,9),qxtrk_->GetBinContent(5,9))/order);


  qanal[anal].qA[bin][0]->Add(qxtrk_,Ax);
  qanal[anal].qA[bin][0]->Add(qytrk_,Ay);
  qanal[anal].qB[bin][0]->Add(qxtrk_,Bx);
  qanal[anal].qB[bin][0]->Add(qytrk_,By);
  qanal[anal].wnA[bin][0]->Add(qcnt_,wA);
  qanal[anal].wnB[bin][0]->Add(qcnt_,wB);
  qanal[anal].qBA[bin][0]->Fill(0.,Bx*Ax + By*Ay);
  qanal[anal].qCA[bin][0]->Fill(0.,Cx*Ax + Cy*Ay);
  qanal[anal].qCB[bin][0]->Fill(0.,Cx*Bx + Cy*By);
  qanal[anal].qBAcnt[bin][0]->Fill(0.,wB*wA);
  qanal[anal].qCAcnt[bin][0]->Fill(0.,wC*wA);
  qanal[anal].qCBcnt[bin][0]->Fill(0.,wC*wB);
  
  int j=(int)(ran->Uniform(0,9.999))+1;
  qanal[anal].qA[bin][j]->Add(qxtrk_,Ax);
  qanal[anal].qA[bin][j]->Add(qytrk_,Ay);
  qanal[anal].qB[bin][j]->Add(qxtrk_,Bx);
  qanal[anal].qB[bin][j]->Add(qytrk_,By);
  qanal[anal].wnA[bin][j]->Add(qcnt_,wA);
  qanal[anal].wnB[bin][j]->Add(qcnt_,wB);
  qanal[anal].qBA[bin][j]->Fill(0.,Bx*Ax + By*Ay);
  qanal[anal].qCA[bin][j]->Fill(0.,Cx*Ax + Cy*Ay);
  qanal[anal].qCB[bin][j]->Fill(0.,Cx*Bx + Cy*By);
  qanal[anal].qBAcnt[bin][j]->Fill(0.,wB*wA);
  qanal[anal].qCAcnt[bin][j]->Fill(0.,wC*wA);
  qanal[anal].qCBcnt[bin][j]->Fill(0.,wC*wB);
}

