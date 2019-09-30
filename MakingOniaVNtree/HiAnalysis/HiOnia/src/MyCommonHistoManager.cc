#include "HiAnalysis/HiOnia/interface/MyCommonHistoManager.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

void
MyCommonHistoManager::Add(std::string theAppendix, std::string theName2) {
  if (!hName.empty()) {
    if (!theName2.empty()) {
      theAppendix += ("_" + theName2);
    }
    myHistos.push_back(new MyCommonHistograms(hName + "_" + theAppendix));
    ret = theCategories.insert( pair<std::string,int>(theAppendix, myHistos.size()-1) );
    if (ret.second==false) {
      edm::LogWarning("MyCommonHistograms") << "[MyCommonHistoManager::Add] Histogram with name: " << hName + "_" + theAppendix << " already existed at: " << ret.first->second;
    }
  }
  else {
    edm::LogWarning("MyCommonHistograms") << "[MyCommonHistoManager::Add] Warning: Histograms do not have a name. Please run SetName(std::string) first" << endl;
    return;
  }

  return;
}

void
MyCommonHistoManager::Print() {
  for (idx = theCategories.begin(); idx!=theCategories.end(); ++idx) {
    cout << (*idx).first << " => " << (*idx).second << endl;
  }
  return;
}

void
MyCommonHistoManager::Fill(const reco::Candidate *p, std::string theFullName) {
  idx = theCategories.find(theFullName);
  if (idx!=theCategories.end()) {
    myHistos.at(idx->second)->Fill(p);
  }
  else {
    edm::LogWarning("MyCommonHistograms") << "[MyCommonHistoManager::Fill] Histogram " << hName << " - " << theFullName << " not found." << endl;
  }

  return;
}

void
MyCommonHistoManager::Fill(const reco::Candidate *p1, const reco::Candidate *p2, std::string theFullName, std::string theName2) {
  idx = theCategories.find(theFullName + "_" + theName2);
  if (idx!=theCategories.end()) {
    myHistos.at(idx->second)->Fill(p1, p2, theName2);
  }
  else {
edm::LogWarning("MyCommonHistograms") << "[MyCommonHistoManager::Fill] Histogram " << hName << " - " << theFullName + "_" + theName2 << " not found." << endl;
  }

  return;
}

void
MyCommonHistoManager::Write(TFile *outf) {
  for (unsigned int i=0; i<myHistos.size(); ++i) {
    myHistos.at(i)->Write(outf);
  }

  return;
}

MyCommonHistograms::MyCommonHistograms(std::string theFullName) {
  hName = theFullName;
  hLabel = MakeLabel(hName);

  if (hName.find("Jpsi")!=std::string::npos ||
      hName.find("Upsilon")!=std::string::npos ||
      hName.find("Chi")!=std::string::npos) {
    useRapidity = true;
  }
  else {
    useRapidity = false;
  }

  booked1ParticleHistos = false;
  booked2ParticleHistos = false;

  // Define a default binning
  theMassBinning = new binning(300,  2.0,   5.0);
  theEBinning    = new binning(400,  0.0,  20.0);
  thePtBinning   = new binning(400,  0.0,  20.0);
  theEtaBinning  = new binning(128, -3.2,   3.2);
  thePhiBinning  = new binning(128, -3.2,   3.2);
  theCentBinning = new binning(40,   0.0, 100.0);
  theCtauBinning = new binning(400, -1.0,   3.0);

  the3dEBinning   = new binning(40,  0.0, 20.0);
  the3dPtBinning  = new binning(40,  0.0, 20.0);
  the3dEtaBinning = new binning(16, -3.2,  3.2);
}

void
MyCommonHistograms::BookParticleHistos() {
  edm::LogInfo("MyCommonHistograms")<<"[MyCommonHistograms::BookParticleHistos] " << hName;

  hMass = new TH1F(("h" + hName + "_mass").c_str(),
		   ("Invariant mass of " + hLabel + ";m_{inv} (" + hLabel + ") [GeV/c^{2}];counts").c_str(),
		   theMassBinning->GetNbins(), theMassBinning->GetMinVal(), theMassBinning->GetMaxVal());

  hE = new TH1F(("h" + hName + "_E").c_str(),
		("Energy of " + hLabel + ";E (" + hLabel + ") [GeV];counts").c_str(),
		theEBinning->GetNbins(), theEBinning->GetMinVal(), theEBinning->GetMaxVal());
  hPt = new TH1F(("h" + hName + "_pt").c_str(),
		 ("Transverse momentum of " + hLabel + ";p_{T} (" + hLabel + ") [GeV/c];counts").c_str(),
		 thePtBinning->GetNbins(), thePtBinning->GetMinVal(), thePtBinning->GetMaxVal());
  hPhi = new TH1F(("h" + hName + "_phi").c_str(),
		  ("Aziumuthal angle of " + hLabel + ";#phi (" + hLabel + ") [rad];counts").c_str(),
		  thePhiBinning->GetNbins(), thePhiBinning->GetMinVal(), thePhiBinning->GetMaxVal());
  if (useRapidity) {
    hEta = new TH1F(("h" + hName + "_eta").c_str(),
		    ("Rapidity of " + hLabel + ";y (" + hLabel + ");counts").c_str(),
		    theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal());
    hE_Eta = new TH2F(("h" + hName + "_E_eta").c_str(),
		      ("Energy vs. rapidity of " + hLabel + ";y (" + hLabel + ");E (" + hLabel + ") [GeV];counts").c_str(),
		      theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal(),
		      theEBinning->GetNbins(), theEBinning->GetMinVal(), theEBinning->GetMaxVal());
    hPt_Eta = new TH2F(("h" + hName + "_pt_eta").c_str(),
		       ("Transverse momentum vs. rapidity of " + hLabel + ";y (" + hLabel + ");p_{T} (" + hLabel + ") [GeV/c];counts").c_str(),
		       theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal(),
		       thePtBinning->GetNbins(), thePtBinning->GetMinVal(), thePtBinning->GetMaxVal());
    hPhi_Eta = new TH2F(("h" + hName + "_phi_eta").c_str(),
			("Azimuthal angle vs. rapidity of " + hLabel + ";y (" + hLabel + ");#phi (" + hLabel + ") [rad];counts").c_str(),
			theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal(),
			thePhiBinning->GetNbins(), thePhiBinning->GetMinVal(), thePhiBinning->GetMaxVal());
    hEta_Mass = new TH2F(("h" + hName + "_eta_mass").c_str(),
			 ("Rapidity vs. invariant mass of " + hLabel + ";m_{inv} (" + hLabel + ") [GeV/c^{2}];y (" + hLabel + ");counts").c_str(),
			 theMassBinning->GetNbins(), theMassBinning->GetMinVal(), theMassBinning->GetMaxVal(),
			 theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal());

  }
  else {
    hEta = new TH1F(("h" + hName + "_eta").c_str(),
		    ("Pseudo-rapidity of " + hLabel + ";#eta (" + hLabel + ");counts").c_str(),
		    theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal());
    hE_Eta = new TH2F(("h" + hName + "_E_eta").c_str(),
		      ("Energy vs. pseudo-rapidity of " + hLabel + ";#eta (" + hLabel + ");E (" + hLabel + ") [GeV];counts").c_str(),
		      theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal(),
		      theEBinning->GetNbins(), theEBinning->GetMinVal(), theEBinning->GetMaxVal());
    hPt_Eta = new TH2F(("h" + hName + "_pt_eta").c_str(),
		       ("Transverse momentum vs. pseudo-rapidity of " + hLabel + ";#eta (" + hLabel + ");p_{T} (" + hLabel + ") [GeV/c];counts").c_str(),
		       theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal(),
		       thePtBinning->GetNbins(), thePtBinning->GetMinVal(), thePtBinning->GetMaxVal());
    hPhi_Eta = new TH2F(("h" + hName + "_phi_eta").c_str(),
			("Azimuthal angle vs. pseudo-rapidity of " + hLabel + ";#eta (" + hLabel + ");#phi (" + hLabel + ") [rad];counts").c_str(),
			theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal(),
			thePhiBinning->GetNbins(), thePhiBinning->GetMinVal(), thePhiBinning->GetMaxVal());
    hEta_Mass = new TH2F(("h" + hName + "_eta_mass").c_str(),
			 ("Pseudo-rapidity vs. invariant mass of " + hLabel + ";m_{inv} (" + hLabel + ") [GeV/c^{2}];#eta (" + hLabel + ");counts").c_str(),
			 theMassBinning->GetNbins(), theMassBinning->GetMinVal(), theMassBinning->GetMaxVal(),
			 theEtaBinning->GetNbins(), theEtaBinning->GetMinVal(), theEtaBinning->GetMaxVal());

  }
    
  hE_Phi = new TH2F(("h" + hName + "_E_phi").c_str(),
		    ("Energy vs. azimuthal angle of " + hLabel + ";#phi (" + hLabel + ") [rad];E (" + hLabel + ") [GeV];counts").c_str(),
		    thePhiBinning->GetNbins(), thePhiBinning->GetMinVal(), thePhiBinning->GetMaxVal(),
		    theEBinning->GetNbins(), theEBinning->GetMinVal(), theEBinning->GetMaxVal());
  hPt_Phi = new TH2F(("h" + hName + "_pt_phi").c_str(),
		     ("Transverse momentum vs. azimuthal angle of " + hLabel + ";#phi (" + hLabel + ") [rad];p_{T} (" + hLabel + ") [GeV/c];counts").c_str(),
		     thePhiBinning->GetNbins(), thePhiBinning->GetMinVal(), thePhiBinning->GetMaxVal(),
		     thePtBinning->GetNbins(), thePtBinning->GetMinVal(), thePtBinning->GetMaxVal());

  hPhi_Mass = new TH2F(("h" + hName + "_phi_mass").c_str(),
		       ("Azimuthal angle vs. invariant mass of " + hLabel + ";m_{inv} (" + hLabel + ") [GeV/c^{2}];#phi (" + hLabel + ") [rad];counts").c_str(),
		       thePhiBinning->GetNbins(), thePhiBinning->GetMinVal(), thePhiBinning->GetMaxVal(),
		       theMassBinning->GetNbins(), theMassBinning->GetMinVal(), theMassBinning->GetMaxVal());

  hE_Mass = new TH2F(("h" + hName + "_E_mass").c_str(),
		     ("Energy vs. invariant mass of " + hLabel + ";m_{inv} (" + hLabel + ") [GeV/c^{2}];E (" + hLabel + ") [GeV];counts").c_str(),
		     theMassBinning->GetNbins(), theMassBinning->GetMinVal(), theMassBinning->GetMaxVal(),
		     theEBinning->GetNbins(), theEBinning->GetMinVal(), theEBinning->GetMaxVal());
  hPt_Mass = new TH2F(("h" + hName + "_pt_mass").c_str(),
		      ("Transverse momentum vs. invariant mass of " + hLabel + ";m_{inv} (" + hLabel + ") [GeV/c^{2}];p_{T} (" + hLabel + ") [GeV/c];counts").c_str(),
		      theMassBinning->GetNbins(), theMassBinning->GetMinVal(), theMassBinning->GetMaxVal(),
		      thePtBinning->GetNbins(), thePtBinning->GetMinVal(), thePtBinning->GetMaxVal());

  hCent_Mass = new TH2F(("h" + hName + "_cent_mass").c_str(),
		      ("Centrality vs. invariant mass of " + hLabel + ";m_{inv} (" + hLabel + ") [GeV/c^{2}];Centrality [%];counts").c_str(),
		      theMassBinning->GetNbins(), theMassBinning->GetMinVal(), theMassBinning->GetMaxVal(),
		      theCentBinning->GetNbins(), theCentBinning->GetMinVal(), theCentBinning->GetMaxVal());

  hCtau_Mass = new TH2F(("h" + hName + "_ctau_mass").c_str(),
			("c#tau vs. invariant mass of " + hLabel + ";m_{inv} (" + hLabel + ") [GeV/c^{2}];l_{" + hLabel + "} [mm];counts").c_str(),
		      theMassBinning->GetNbins(), theMassBinning->GetMinVal(), theMassBinning->GetMaxVal(),
		      theCtauBinning->GetNbins(), theCtauBinning->GetMinVal(), theCtauBinning->GetMaxVal());

  hMass->Sumw2();
  hE->Sumw2();
  hPt->Sumw2();
  hPhi->Sumw2();
  hEta->Sumw2();

  hE_Eta->Sumw2();
  hPt_Eta->Sumw2();
  hPhi_Eta->Sumw2();

  hE_Phi->Sumw2();
  hPt_Phi->Sumw2();

  hE_Mass->Sumw2();
  hPt_Mass->Sumw2();
  hEta_Mass->Sumw2();
  hPhi_Mass->Sumw2();
  hCent_Mass->Sumw2();
  hCtau_Mass->Sumw2();

  
  booked1ParticleHistos = true;
  return;
}


void
MyCommonHistograms::BookParticleHistos(std::string hName2) {
  edm::LogInfo("MyCommonHistograms")<<"[MyCommonHistograms::BookParticleHistos] " << hName;
  hLabel2 = MakeLabel(hName2);

  /*
    add stuff like J/psi vs. muon
   */

  booked2ParticleHistos = true;
  return;
}

void
MyCommonHistograms::Fill(const reco::Candidate *p) {
  if (!booked1ParticleHistos) BookParticleHistos();

  hMass->Fill(p->mass());
  hE->Fill(p->energy());
  hPt->Fill(p->pt());
  hPhi->Fill(p->phi());

  if (useRapidity) {
    hEta->Fill(p->rapidity());
    
    hE_Eta->Fill(p->rapidity(), p->energy());
    hPt_Eta->Fill(p->rapidity(), p->pt());
    hPhi_Eta->Fill(p->rapidity(), p->phi());

    hEta_Mass->Fill(p->mass(), p->rapidity());
  }
  else {
    hEta->Fill(p->eta());
    
    hE_Eta->Fill(p->eta(), p->energy());
    hPt_Eta->Fill(p->eta(), p->pt());
    hPhi_Eta->Fill(p->eta(), p->phi());

    hEta_Mass->Fill(p->mass(), p->eta());
  }

  hE_Phi->Fill(p->phi(), p->energy());
  hPt_Phi->Fill(p->phi(), p->pt());
  
  hPhi_Mass->Fill(p->mass(), p->phi());
  hE_Mass->Fill(p->mass(), p->energy());
  hPt_Mass->Fill(p->mass(), p->pt());

  double theCentrality;
  double theCtau;
  //  double theCtauErr;

  if (hName.find("Jpsi")!=std::string::npos) {
    const pat::CompositeCandidate * cand = dynamic_cast<const pat::CompositeCandidate*> (p->clone());

    theCentrality = 2.5*cand->userInt("centBin");
    hCent_Mass->Fill(p->mass(), theCentrality);

    theCtau = 10.0*cand->userFloat("ppdlPV");
    //    theCtauErr = 10.0*cand->userFloat("ppdlErrPV");
  
    
    hCtau_Mass->Fill(p->mass(), theCtau);
  }

  return;
}

void
MyCommonHistograms::Fill(const reco::Candidate *p1, const reco::Candidate *p2, std::string hName2) {
  if (!booked2ParticleHistos) BookParticleHistos(hName2);  

  /*
    add stuff
   */

  return;
}

void
MyCommonHistograms::Write(TFile *outf) {

  if (!outf->IsOpen()) {
    cerr << "No output file defined. Writing output to DummyHistos.root" << endl;
    outf = new TFile("DummyHistos.root","RECREATE");
  }


  if (booked1ParticleHistos) {
    if(outf->FindKey(hName.c_str())==0)
      outf->mkdir(hName.c_str())->cd();

    hMass->Write();
    hE->Write();
    hPt->Write();
    hPhi->Write();
    hEta->Write();

    hE_Eta->Write();
    hPt_Eta->Write();
    hPhi_Eta->Write();

    hE_Phi->Write();
    hPt_Phi->Write();

    hE_Mass->Write();
    hPt_Mass->Write();
    hEta_Mass->Write();
    hPhi_Mass->Write();
    hCent_Mass->Write();
    hCtau_Mass->Write();
  }
  
  if (booked2ParticleHistos) {
    if(outf->FindKey(hName.c_str())==0)
      outf->mkdir(hName.c_str())->cd();
    /*
      add stuff
     */
  }

  outf->cd();
  return;
}


std::string
MyCommonHistograms::MakeLabel(std::string name) {
  std::string label;

  size_t stop;

  stop = name.find("_");
  if(stop==std::string::npos)
    stop = name.length();

  if (name.find("GenPhoton")<stop)
    label = "#gamma_{gen}";
  else if (name.find("GenMuon")<stop)
    label = "#mu_{gen}";
  else if (name.find("GenJpsi")<stop)
    label = "J/#psi_{gen}";
  else if (name.find("GenChic")<stop)
    label = "#chi_{c, gen}";
  else if (name.find("GenUpsilon")<stop)
    label = "#Upsilon_{gen}";
  else if (name.find("GenChib")<stop)
    label = "#chi_{b, gen}}";
  else if (name.find("RecoPhoton")<stop)
    label = "#gamma_{reco}";
  else if (name.find("RecoMuon")<stop)
    label = "#mu_{reco}";
  else if (name.find("GlobalMuon")<stop)
    label = "#mu_{glb}";
  else if (name.find("TrackerMuon")<stop)
    label = "#mu_{trk}";
  else if (name.find("CaloMuon")<stop)
    label = "#mu_{cal}";
  else if (name.find("RecoJpsi")<stop)
    label = "J/#psi_{reco}";
  else if (name.find("GlbGlbJpsi")<stop)
    label = "J/#psi (global+global)";
  else if (name.find("GlbTrkJpsi")<stop)
    label = "J/#psi (global+tracker)";
  else if (name.find("TrkTrkJpsi")<stop)
    label = "J/#psi (tracker+tracker)";
  else if (name.find("GlbCalJpsi")<stop)
    label = "J/#psi (global+calo)";
  else if (name.find("RecoChic")<stop)
    label = "#chi_{c, reco}";
  else if (name.find("GlbGlbChic")<stop)
    label = "#chi_{c} (global+global)";
  else if (name.find("GlbTrkChic")<stop)
    label = "#chi_{c} (global+tracker)";
  else if (name.find("TrkTrkChic")<stop)
    label = "#chi_{c} (tracker+tracker)";
  else if (name.find("GlbCalChic")<stop)
    label = "#chi_{c} (global+calo)";
  else if (name.find("RecoUpsilon")<stop)
    label = "#Upsilon_{reco}";
  else if (name.find("RecoChib")<stop)
    label = "#chi_{b, reco}}";
  else 
    label = "some particle";

  return label;
}
