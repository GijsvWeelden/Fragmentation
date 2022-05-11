Int_t is_charged(TParticle *part) {
  Int_t abs_kf = abs(part->GetPdgCode());
  
  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf != 310 && abs_kf != 3122 && abs_kf != 3112  && abs_kf != 3222  && abs_kf != 3312 && abs_kf != 3322 && abs_kf != 3334 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16 && abs_kf != 3124)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

void analysePythia(Char_t const *indir = "output/events_40_80_50", Char_t const *foutname = "qPythiaJetSpectra.root") {
  gSystem->AddDynamicPath("/usr/lib64");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  //gSystem->Load("libTENDER");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCAL");
  
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  //gSystem->Load("libSISConePlugin");
  //gSystem->Load("libCDFConesPlugin");
  gSystem->Load("libfastjetcontribfragile");
  gSystem->Load("libPWGJEEMCALJetTasks");

  gSystem->Load("libqpythia");
  gSystem->Load("libEGPythia6");
  gSystem->Load("liblhapdf");
  gSystem->Load("libAliPythia6");


  const Float_t eta_jet_max = 1.0;

  const Int_t chargedOnly = 0;  // charged only or full jets
  // Constituents are charged only

  const Float_t pt_cut_rho = 1;

  const Int_t nPtJetBins = 250;
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 250;
  TFile *fout = new TFile(foutname,"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);
 
  TH1F *hPtPion = new TH1F("hPtPion","charged pion pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  TH1F *hPtCharged = new TH1F("hPtCharged","charged particle pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  TH1F *hPtAll = new TH1F("hPtAll","all particle pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  //const Float_t Rvals[]={0.2,0.3,0.4,0.5};
  const Float_t Rvals[]={0.02,0.05,0.1,0.15};
  //const Float_t Rvals[]={0.4,0.6,0.8,1.0};
  //static const Int_t nR = sizeof(*Rvals)/sizeof(Float_t);
  const Int_t nR = 4;
  TH3F *hJetPtEtaPhi[nR] = {0};
  TH3F *hJetPtEtaMass[nR] = {0};
  TH3F *hLeadJetPtEtaPhi[nR] = {0};
  TH3F *hPtJetDeltaEtaDeltaPhiN[nR] = {0};
  TH3F *hPtJetDeltaEtaDeltaPhiPt[nR] = {0};
  TH2F *hPtJetRN[nR] = {0};
  TH2F *hPtJetRPt[nR] = {0};
  TH2F *hPtJetZ[nR] = {0};
  TH2F *hPtJetPt[nR] = {0};
  TH2F *hPtJetXi[nR] = {0};
  Float_t max_delta_eta = 2.5;
  for (Int_t iR = 0; iR < nR; iR++) {
    hJetPtEtaPhi[iR] = new TH3F(Form("hJetPtEtaPhi_R%.0f",10*Rvals[iR]),Form("hJetPtEtaPhi R=%.1f;p_{T,jet} (GeV);#eta;#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,32,0,2*TMath::Pi());
    hJetPtEtaMass[iR] = new TH3F(Form("hJetPtEtaMass_R%.0f",10*Rvals[iR]),Form("hJetPtEtaMass R=%.1f;p_{T,jet} (GeV);#eta;jet mass (GeV/c^{2})",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,200,0,50);
    hLeadJetPtEtaPhi[iR] = new TH3F(Form("hLeadJetPtEtaPhi_R%.0f",10*Rvals[iR]),Form("hLeadJetPtEtaPhi R=%.1f;p_{T,jet} (GeV);#eta;#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,32,0,2*TMath::Pi());
    hPtJetDeltaEtaDeltaPhiN[iR] = new TH3F(Form("hPtJetDeltaEtaDeltaPhiN_R%.0f",10*Rvals[iR]),Form("hPtJetDeltaEtaDeltaPhiN R=%.1f;p_{T,jet} (GeV);#Delta#eta;#Delta#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,100,-max_delta_eta,max_delta_eta,126,-TMath::Pi(),TMath::Pi());
    hPtJetDeltaEtaDeltaPhiPt[iR] = new TH3F(Form("hPtJetDeltaEtaDeltaPhiPt_R%.0f",10*Rvals[iR]),Form("hPtJetDeltaEtaDeltaPhiPt R=%.1f;p_{T,jet} (GeV);#Delta#eta;#Delta#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,100,-max_delta_eta,max_delta_eta,126,-TMath::Pi(),TMath::Pi());
    hPtJetDeltaEtaDeltaPhiPt[iR]->Sumw2();

    hPtJetRN[iR] = new TH2F(Form("hPtJetRN_R%.0f",10*Rvals[iR]),Form("hPtJetRN R=%.1f;p_{T,jet} (GeV);R",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,30,0,1.5);
    hPtJetRPt[iR] = new TH2F(Form("hPtJetRPt_R%.0f",10*Rvals[iR]),Form("hPtJetRN R=%.1f;p_{T,jet} (GeV);R;p_{T}",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,30,0,1.5);
    hPtJetRPt[iR]->Sumw2();

    hPtJetPt[iR] = new TH2F(Form("hPtJetPt_R%.0f",10*Rvals[iR]),Form("hPtJetPt R=%.1f;p_{T,jet} (GeV);p_{T,part}",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
    hPtJetZ[iR] = new TH2F(Form("hPtJetZ_R%.0f",10*Rvals[iR]),Form("hPtJetZ R=%.1f;p_{T,jet} (GeV);z=p_{T,part}/p_{T,jet}",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,50,0,1);
    hPtJetXi[iR] = new TH2F(Form("hPtJetXi_R%.0f",10*Rvals[iR]),Form("hPtJetXi R=%.1f;p_{T,jet} (GeV);#xi=log(p_{T,jet}/p_{T,part})",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,70,0,7);
  }
  cout << "fout.ls() after creating histos" << endl;
  fout->ls();


  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

  AliEmcalJetFinder *jetFinder = new AliEmcalJetFinder();
  jetFinder->SetJetAlgorithm(0);  // 0 = anti-kt  1 = kt
  jetFinder->SetRecombSheme(0);  // 0 = E_scheme


  void *inDirP = gSystem->OpenDirectory(indir);
  Char_t *subdir = 0;
  while (subdir = gSystem->GetDirEntry(inDirP)) {
    if (strcmp(subdir,".")==0 || strcmp(subdir,"..") == 0)
      continue;

    TString curDir(indir);
    curDir += "/";
    curDir += subdir;

    // skip first 50 dirs; temp test
    Int_t num = atoi(subdir);
    /*
    if (num<=50)
      continue;
    */
    cout << "input directory: " << curDir << " num " << num << endl;
     
    AliRunLoader *rl = AliRunLoader::Open(curDir+"/galice.root");  
    if (rl == 0)
      continue;
    
    rl->LoadKinematics();

    cout << "fout.ls() after LoadKinematics" << endl;
    fout->ls();

    TObjArray plistSel(1000);
    
    Int_t iEvent = 0;
    while (rl->GetEvent(iEvent)==0) {

      TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet

      AliStack *stack = rl->Stack();

      plistSel.Clear();

      Int_t n_trk = 0;
      Int_t n_part = stack->GetNtrack();
      for (Int_t i_part = 0; i_part < n_part; i_part++) {
	TParticle *part=(TParticle*) stack->Particle(i_part);

	if (part->GetStatusCode() >= 10)  // Not a final state particle 
	  continue;

	if (chargedOnly) {
	  if (!is_charged(part)) 
	    continue;
	}
	else { // full jets
	  if (part->GetPdgCode()==12 || part->GetPdgCode()==14 || part->GetPdgCode()==16 ) // reject neutrinos
	    continue;
	}

	if (abs(part->GetPdgCode()) == 211) {
	  hPtPion->Fill(part->Pt());
	}
	hPtAll->Fill(part->Pt());
        if (is_charged(part)) {
	  hPtCharged->Fill(part->Pt());
        }
	plistSel.AddLast(part);

	//}
	// Only use charged tracks?

	n_trk++;
      }
      //cout << n_trk << " particles filled" << endl;
      //cout << JetFinderEvent.GetNCalTrkTracks() << " in event" << endl;

      for (Int_t iR = 0; iR < nR; iR++) {

	for (Int_t iPart = 0; iPart < plistSel.GetEntriesFast(); iPart++) {
	  TParticle *part = (TParticle*) plistSel.At(iPart);
	  jetFinder->AddInputVector(part->Px(), part->Py(), part->Pz(), part->Energy());
	}
	jetFinder->SetRadius(Rvals[iR]); 
	//header->SetGhostEtaMax(2);
	//header->SetGhostArea(0.05);

	jetFinder->FindJets();

	Float_t lead_pt = -1;
 
	Float_t lead_phi = -999;
	Float_t lead_eta = -999;
	for (Int_t iJet = 0; iJet < jetFinder->GetNumberOfJets(); iJet++) {
	  AliEmcalJet *jet = jetFinder->GetJet(iJet);
	  //cout << "\t jet " << iJet << " pt " << jet->Pt() << " eta " << jet->Eta() << " phi " << jet->Phi() << endl; 
	  if (TMath::Abs(jet->Eta()) < eta_jet_max) {
	    hJetPtEtaPhi[iR]->Fill(jet->Pt(),jet->Eta(),jet->Phi());
	    hJetPtEtaMass[iR]->Fill(jet->Pt(),jet->Eta(),jet->M());
	    if (jet->Pt() > lead_pt) {
	      lead_pt = jet->Pt();
	      lead_eta = jet->Eta();
	      lead_phi = jet->Phi();
	    }
	  }
	}
	//cout << "lead_pt " << lead_pt << endl;
	if (lead_pt > 0) {
	  hLeadJetPtEtaPhi[iR]->Fill(lead_pt,lead_eta,lead_phi);
	  //cout << "aliplist entries: " << aliplist.GetEntriesFast() << " fast " << aliplist.GetEntriesFast() << " size " << aliplist.GetSize() << endl;
	  //Int_t n_trk2 = 0;
	  //Float_t sum_pt = 0;
	  for (Int_t ipart = 0; ipart < plistSel.GetEntriesFast(); ipart++) {
	    TParticle *part = (TParticle*) plistSel.At(ipart);
	    //cout << "particle " << part << " phi " << part->Pt() << endl;
	    if (!is_charged(part))
	      continue;
	    Float_t deta = part->Eta() - lead_eta;
	    if (TMath::Abs(deta) < max_delta_eta) {
	      Float_t dphi = part->Phi() - lead_phi;
	      if (dphi < -TMath::Pi())
		dphi += 2*TMath::Pi();
	      if (dphi > TMath::Pi())
		dphi -= 2*TMath::Pi();
	      //cout << "fill" << hPtJetDeltaEtaDeltaPhi[iR] << endl;
	      hPtJetDeltaEtaDeltaPhiN[iR]->Fill(lead_pt, deta, dphi);
	      hPtJetDeltaEtaDeltaPhiPt[iR]->Fill(lead_pt, deta, dphi, part->Pt());
	      if (TMath::Abs(dphi) < 1.5 && TMath::Abs(deta) < 1.5) {
		Float_t Rdist =  TMath::Sqrt(dphi*dphi + deta*deta);
		if (part->Pt() > pt_cut_rho) {
		  hPtJetRN[iR]->Fill(lead_pt, Rdist);
		  hPtJetRPt[iR]->Fill(lead_pt, Rdist, part->Pt()/lead_pt); // Divide by pt, like CMS does
		}
	      
		//if (TMath::Abs(dphi) < Rvals[iR] && TMath::Abs(deta) < Rvals[iR] && Rdit < Rvals[iR])
		if (Rdist < Rvals[iR])
		  {
		    hPtJetPt[iR]->Fill(lead_pt, part->Pt());
		    hPtJetZ[iR]->Fill(lead_pt, part->Pt()/lead_pt);
		    hPtJetXi[iR]->Fill(lead_pt, log(lead_pt/part->Pt()));
		    //sum_pt+= part->Pt();
		  }
	      }
	    }
	    //n_trk2++;
	  }
	  //cout << "tracks for deta, dphi " << n_trk2 << endl;
	  //cout << "sum pt in cone " << sum_pt << endl;
	}
      }
      iEvent++;
    }
    delete rl;

    cout << "fout.ls() after deleting runloader" << endl;
    fout->ls();

    // copy event and xsec info
    TFile fin(curDir+"/pythia_xsec.root");
    hNEventIn = (TH1D*) fin.Get("hNEvent");
    hNEvent->Fill(0.5,hNEventIn->GetBinContent(1));  
    hXSecIn = (TProfile*) fin.Get("hXSec");
    hXSec->Fill(0.5,hXSecIn->GetBinContent(1));  

    cout << "fout.ls() after filling xsec" << endl;
    fout->ls();

  }
  fout->cd();
  cout << "fout.ls() after everything" << endl;
  fout->ls();


  fout->Write();
}
