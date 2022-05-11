Int_t is_charged(TParticle *part) {
  Int_t abs_kf = abs(part->GetPdgCode());
  
  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf != 310 && abs_kf != 3122 && abs_kf != 3112  && abs_kf != 3222  && abs_kf != 3312 && abs_kf != 3322 && abs_kf != 3334 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16 && abs_kf != 3124)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

void pythia6JetSpectra(Int_t nEvent = 100, Char_t const *foutname = "PythiaJetSpectra.root", Float_t weight_pow = 6, Float_t e_cms = 2760) {
  //
  // example macro for on-the-fly generation of PYTHIA6 events
  // and analysis with new reader interface and FastJet
  // M. van Leeuwen

  /*
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCAL");
  gSystem->AddDynamicPath("/usr/lib64");
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  //gSystem->Load("libSISConePlugin");
  //gSystem->Load("libCDFConesPlugin");
  gSystem->Load("libfastjetcontribfragile");
  gSystem->Load("libPWGJEEMCALJetTasks");
  */ 

  //gSystem->Load("libqpythia");
  //gSystem->Load("libpythia6_4_28"); // Does not have pyevwt?
  gSystem->Load("libpythia6");
  gSystem->Load("libEGPythia6");
  gSystem->Load("liblhapdf");
  gSystem->Load("libAliPythia6");

  const Float_t max_eta = 1.0; // particles
  const Int_t nEtaBins = 20;

  const Float_t eta_jet_max = 1.0; //0.5;

  const Int_t chargedOnly = 1;  // charged only or full jets

  const Int_t nPtJetBins = 250;
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 250;
  TFile *fout = new TFile(foutname,"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; ;N",1,0,1);
  hNEvent->Sumw2();
  TH1F *hNTrials = new TH1F("hNTrials","number of trials; ;N",1,0,1);
  hNTrials->Sumw2();

  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);
 
  TH1F *hPtHard = new TH1F("hPtHard","PThard;p_{T,hard} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  hPtHard->Sumw2();
  TH1F *hPtHardUnweighted = new TH1F("hPtHardUnweighted","PThard (unweighted);p_{T,hard} (GeV)",nPtJetBins,minPtJet,maxPtJet);

  TH1F *hPtPion = new TH1F("hPtPion","charged pion pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  hPtPion->Sumw2();
  TH2F *hEtaPtCharged = new TH2F("hEtaPtCharged","charged particle pt spectrum;#eta;p_{T} (GeV)",nEtaBins,-max_eta,max_eta,nPtJetBins,minPtJet,maxPtJet);
  hEtaPtCharged->Sumw2();

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  const Float_t Rvals[]={0.2,0.3,0.4};
  //static const Int_t nR = sizeof(*Rvals)/sizeof(Float_t);
  const Int_t nR = 3;
  TH3F *hJetPtEtaPhiUnweighted[nR] = {0};
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
    hJetPtEtaPhiUnweighted[iR] = new TH3F(Form("hJetPtEtaPhiUnweighted_R%.0f",10*Rvals[iR]),Form("hJetPtEtaPhiUnweighted R=%.1f;p_{T,jet} (GeV);#eta;#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,32,0,2*TMath::Pi());
    hJetPtEtaPhi[iR] = new TH3F(Form("hJetPtEtaPhi_R%.0f",10*Rvals[iR]),Form("hJetPtEtaPhi R=%.1f;p_{T,jet} (GeV);#eta;#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,32,0,2*TMath::Pi());
    hJetPtEtaPhi[iR]->Sumw2();
    hJetPtEtaMass[iR] = new TH3F(Form("hJetPtEtaMass_R%.0f",10*Rvals[iR]),Form("hJetPtEtaMass R=%.1f;p_{T,jet} (GeV);#eta;jet mass (GeV/c^{2})",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,200,0,50);
    hJetPtEtaMass[iR]->Sumw2();
    hLeadJetPtEtaPhi[iR] = new TH3F(Form("hLeadJetPtEtaPhi_R%.0f",10*Rvals[iR]),Form("hLeadJetPtEtaPhi R=%.1f;p_{T,jet} (GeV);#eta;#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,32,0,2*TMath::Pi());
    hLeadJetPtEtaPhi[iR]->Sumw2();
    hPtJetDeltaEtaDeltaPhiN[iR] = new TH3F(Form("hPtJetDeltaEtaDeltaPhiN_R%.0f",10*Rvals[iR]),Form("hPtJetDeltaEtaDeltaPhiN R=%.1f;p_{T,jet} (GeV);#Delta#eta;#Delta#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,100,-max_delta_eta,max_delta_eta,126,-TMath::Pi(),TMath::Pi());
    hPtJetDeltaEtaDeltaPhiN[iR]->Sumw2();
    hPtJetDeltaEtaDeltaPhiPt[iR] = new TH3F(Form("hPtJetDeltaEtaDeltaPhiPt_R%.0f",10*Rvals[iR]),Form("hPtJetDeltaEtaDeltaPhiPt R=%.1f;p_{T,jet} (GeV);#Delta#eta;#Delta#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,100,-max_delta_eta,max_delta_eta,126,-TMath::Pi(),TMath::Pi());
    hPtJetDeltaEtaDeltaPhiPt[iR]->Sumw2();

    hPtJetRN[iR] = new TH2F(Form("hPtJetRN_R%.0f",10*Rvals[iR]),Form("hPtJetRN R=%.1f;p_{T,jet} (GeV);R",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,30,0,1.5);
    hPtJetRN[iR]->Sumw2();
    hPtJetRPt[iR] = new TH2F(Form("hPtJetRPt_R%.0f",10*Rvals[iR]),Form("hPtJetRN R=%.1f;p_{T,jet} (GeV);R;p_{T}",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,30,0,1.5);
    hPtJetRPt[iR]->Sumw2();

    hPtJetZ[iR] = new TH2F(Form("hPtJetZ_R%.0f",10*Rvals[iR]),Form("hPtJetZ R=%.1f;p_{T,jet} (GeV);z=p_{T,part}/p_{T,jet}",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,50,0,1);
  }

  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

  AliGenPythia *pythia=new AliGenPythia(1);
  //pythia->SetProcess(kPyZgamma+1); //kPyJetsHardOnly); not sure why enum does not work..
  pythia->SetProcess(kPyJets);
  //pythia->SetPtHard(pthard_min, pthard_max);
  //pythia->SetPtHard(5, 300);
  pythia->SetPtHard(5, 57);
  pythia->SetWeightPower(weight_pow);
  pythia->SetEnergyCMS(e_cms);
  pythia->Init();

  AliEmcalJetFinder *jetFinder = new AliEmcalJetFinder();
  jetFinder->SetJetAlgorithm(0);  // 0 = anti-kt  1 = kt
  jetFinder->SetRecombSheme(1);  // 0 = E_scheme  1 = pt_scheme

  jetFinder->SetJetMaxEta(eta_jet_max);
  jetFinder->SetTrackMaxEta(max_eta);

  AliStack *stack = new AliStack();
  // There is a big mess with the connection between runloader, AliRun and the gAlice pointer. 
  // This order of things seems to work...
  AliRunLoader *dummyrl = new AliRunLoader("dummyEvtFolder");
  dummyrl->MakeHeader();
  dummyrl->SetEventNumber(0);
  gAlice->SetRunLoader(dummyrl);

  pythia->SetStack(stack);

  TObjArray plistSel(1000);

  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    if ((iEvent%1000)==0) {
      cout << "Event " << iEvent << endl;
      ProcInfo_t proc_i;
      gSystem->GetProcInfo(&proc_i);
      cout << " cpu sys " << proc_i.fCpuSys << " user " << proc_i.fCpuUser << endl; 
    }

    stack->Reset();

    TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet

    pythia->Generate();
    AliGenPythiaEventHeader *header = dynamic_cast<AliGenPythiaEventHeader*>(dummyrl->GetHeader()->GenEventHeader());
    //cout << "pthard " << header->GetPtHard() << " Weight " << header->EventWeight() << endl;
    Double_t evt_wght = header->EventWeight();
    
    hNEvent->Fill(0.5,evt_wght); //0.5);
    hNTrials->Fill(0.5,header->Trials());
    hPtHard->Fill(header->GetPtHard(),evt_wght);
    hPtHardUnweighted->Fill(header->GetPtHard());
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
	hPtPion->Fill(part->Pt(),evt_wght);
      }
      if (chargedOnly || is_charged(part))
         hEtaPtCharged->Fill(part->Eta(),part->Pt(),evt_wght);
      if (TMath::Abs(part->Eta()) < max_eta)
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
	  hJetPtEtaPhiUnweighted[iR]->Fill(jet->Pt(),jet->Eta(),jet->Phi());
	  hJetPtEtaPhi[iR]->Fill(jet->Pt(),jet->Eta(),jet->Phi(),evt_wght);
	  hJetPtEtaMass[iR]->Fill(jet->Pt(),jet->Eta(),jet->M(),evt_wght);
          if (jet->Pt() > lead_pt) {
	    lead_pt = jet->Pt();
	    lead_eta = jet->Eta();
	    lead_phi = jet->Phi();
          }
	}
      }
      //cout << "lead_pt " << lead_pt << endl;
      if (lead_pt > 0) {
	hLeadJetPtEtaPhi[iR]->Fill(lead_pt,lead_eta,lead_phi,evt_wght);
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
	    hPtJetDeltaEtaDeltaPhiN[iR]->Fill(lead_pt, deta, dphi, evt_wght);
	    hPtJetDeltaEtaDeltaPhiPt[iR]->Fill(lead_pt, deta, dphi, part->Pt()*evt_wght);
	    if (TMath::Abs(dphi) < 1.5 && TMath::Abs(deta) < 1.5) {
	      Float_t Rdist =  TMath::Sqrt(dphi*dphi + deta*deta);
	      hPtJetRN[iR]->Fill(lead_pt, Rdist, evt_wght);
	      hPtJetRPt[iR]->Fill(lead_pt, Rdist, part->Pt()*evt_wght);
	      
	      //if (TMath::Abs(dphi) < Rvals[iR] && TMath::Abs(deta) < Rvals[iR] && Rdit < Rvals[iR])
	      if (Rdist < Rvals[iR])
		{
		  hPtJetZ[iR]->Fill(lead_pt, part->Pt()/lead_pt, evt_wght);
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

  hXSec->Fill(0.5,pythia->GetXsection());

  cout << hNEvent->GetBinContent(1) << " events, " << hNTrials->GetBinContent(1) << " trials, xsec " << pythia->GetXsection() << endl;
  
  fout->Write();

  ProcInfo_t proc_i;
  gSystem->GetProcInfo(&proc_i);
  cout << "Total: cpu sys " << proc_i.fCpuSys << " user " << proc_i.fCpuUser << endl; 
}
