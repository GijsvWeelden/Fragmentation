Int_t is_charged(TParticle *part) {
  Int_t abs_kf = abs(part->GetPdgCode());
  
  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16 && abs_kf != 3124)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

void pythia6JetSpectra(Int_t nEvent = 50, Char_t const *foutname = "qPythiaJetSpectra.root", Float_t pthard_min=10, Float_t pthard_max=100, Float_t qHat = 50., Float_t e_cms = 2760) {
  // NB: qHat is mean qhat in units 0.1 GeV^2/fm

  // example macro for on-the-fly generation of PYTHIA6 events
  // and analysis with new reader interface and FastJet
  // M. van Leeuwen

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  gSystem->Load("libTENDER");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCAL");
  
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  //gSystem->Load("libsiscone");
  //gSystem->Load("libsiscone_spherical");
  //gSystem->Load("libfastjetplugins");
  gSystem->Load("libSISConePlugin");
  gSystem->Load("libCDFConesPlugin");
  //gSystem->Load("libfastjetcontribfragile");
  gSystem->Load("libPWGJEEMCALJetTasks");

  gSystem->Load("libpythia6.so");
  gSystem->Load("libEGPythia6.so");
  gSystem->Load("libAliPythia6.so");


  const Float_t eta_jet_max = 1.0;

  const Int_t chargedOnly = 1;  // charged only or full jets

  const Int_t nPtJetBins = 250;
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 250;
  TFile *fout = new TFile(foutname,"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);
 
  TH1F *hPtPion = new TH1F("hPtPion","charged pion pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  TH1F *hPtCharged = new TH1F("hPtCharged","charged particle pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  const Float_t Rvals[]={0.2,0.3,0.4};
  //static const Int_t nR = sizeof(*Rvals)/sizeof(Float_t);
  const Int_t nR = 3;
  TH3F *hJetPtEtaPhi[nR] = {0};
  TH3F *hJetPtEtaMass[nR] = {0};
  TH3F *hLeadJetPtEtaPhi[nR] = {0};
  TH3F *hPtJetDeltaEtaDeltaPhiN[nR] = {0};
  TH3F *hPtJetDeltaEtaDeltaPhiPt[nR] = {0};
  TH2F *hPtJetRN[nR] = {0};
  TH2F *hPtJetRPt[nR] = {0};
  TH2F *hPtJetZ[nR] = {0};
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

    hPtJetZ[iR] = new TH2F(Form("hPtJetZ_R%.0f",10*Rvals[iR]),Form("hPtJetZ R=%.1f;p_{T,jet} (GeV);z=p_{T,part}/p_{T,jet}",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,50,0,1);
  }

  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

  AliPythia6 *pythia=AliPythia6::Instance();

  pythia->SetCKIN(3,pthard_min);   // minimum hard pt
  pythia->SetCKIN(4,pthard_max);  // maximum hard pt

  pythia->SetMDCY(pythia->Pycomp(111),1,0);  // switch off pi0 decay

  pythia->Initialize("CMS","p","p",e_cms);

  if (qHat != 0) { // Switch on quenching
    // Set qhat : <qhat> ~ k and <qhat>=1.7 GeV^2/fm for k=0.6*10^6 (default AliPythia)
    Float_t myk=qHat*0.6e6/17.;
    pythia->InitQuenching(0,0.1,myk,0);
  }


  AliEmcalJetFinder *jetFinder = new AliEmcalJetFinder();
  jetFinder->SetJetAlgorithm(0);  // 0 = anti-kt  1 = kt
  jetFinder->SetRecombSheme(0);  // 0 = E_scheme

  TClonesArray *plist=new TClonesArray("TParticle",1000);
  TObjArray plistSel(1000);

  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet

    // Switch off hadronization to generate only partonic scattering 
    pythia->SetMSTJ(1, 0);
    //pythia->GenerateEvent(); 
    pythia->Pyevnt();

    if (qHat != 0)
      pythia->Quench();

    // Do hadronization
    pythia->SetMSTJ(1, 1);
    pythia->Pyexec();
    
    //pythia->ImportParticles();
    //TObjArray *plist=pythia->GetListOfParticles();

    hNEvent->Fill(0.5);

    pythia->GetParticles(plist);

    plistSel.Clear();

    Int_t n_trk = 0;
    Int_t n_part = plist->GetEntries();
    for (Int_t i_part = 0; i_part < n_part; i_part++) {
      TParticle *part=(TParticle*) plist->At(i_part);

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
      hPtCharged->Fill(part->Pt());
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
	      hPtJetRN[iR]->Fill(lead_pt, Rdist);
	      hPtJetRPt[iR]->Fill(lead_pt, Rdist, part->Pt());
	      
	      //if (TMath::Abs(dphi) < Rvals[iR] && TMath::Abs(deta) < Rvals[iR] && Rdit < Rvals[iR])
	      if (Rdist < Rvals[iR])
		{
		  hPtJetZ[iR]->Fill(lead_pt, part->Pt()/lead_pt);
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
  }

  hXSec->Fill(0.5,pythia->GetPARI(1));

  fout->Write();
}
