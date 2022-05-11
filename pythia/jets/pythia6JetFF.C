Int_t is_charged(TParticle *part) {
  Int_t abs_kf = abs(part->GetPdgCode());
  
  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16 && abs_kf != 3124)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

void pythia6JetFF(Int_t nEvent = 50, Char_t const *foutname = "pythiaJetFF.root", Float_t pthard_min=10, Float_t pthard_max=100, Float_t e_cms = 7000) {

  // Macro that runs pythia (in pt hard bins) and fills some
  // basic histos for fragmentation functions.
  //
  // Based on the example macro in $ALICE_ROOT/PWGJE/macros/examples
  //
  // Marco van Leeuwen, Nikhef and UU

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  // Fastjet v2.4
  //gSystem->Load("libsiscone");
  //gSystem->Load("libSISConePlugin");
  // Fastjet v3.0
  gSystem->Load("libsiscone"); 
  gSystem->Load("libsiscone_sperical"); 
  gSystem->Load("libfastjetplugins"); 
  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");

  gSystem->Load("libpythia6.so");
  gSystem->Load("libEGPythia6.so");
  gSystem->Load("libAliPythia6.so");


  const Float_t jetR = 0.4;
  const Int_t nEtaBins = 50;
  const Float_t maxEtaJet = 5.0;

  const Int_t chargedOnly = 0;  // charged only or full jets

  const Int_t nPtJetBins = 100;
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 100;
  TFile *fout = new TFile(foutname,"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);
 
  TH1F *hEtaPtPion = new TH1F("hEtaPtPion","charged pion pt spectrum;#eta;p_{T} (GeV)",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet);
  TH1F *hPtD = new TH1F("hPtD","D^{0} pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);

  TH3F *hJetPtEtaPhi = new TH3F("hJetPtEtaPhi","hJetPtEtaPhi;p_{T,jet} (GeV);#eta;#phi",nPtJetBins,minPtJet,maxPtJet,nEtaBins,-maxEtaJet,maxEtaJet,32,0,2*TMath::Pi());
  TH3F *hEtaPtJetPtCh = new TH3F("hEtaPtJetPtCh","hEtaPtJetPtCh;p_{T,jet} (GeV);p_{T,part} (GeV)",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
  TH3F *hEtaPtJetPtPion = new TH3F("hEtaPtJetPtPion","hEtaPtJetPtPion;p_{T,jet} (GeV);p_{T,part} (GeV)",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
  TH3F *hEtaPtJetZPion = new TH3F("hEtaPtJetZPion","hEtaPtJetZPion;p_{T,jet} (GeV);z",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet,44,0,1.1);
  TH3F *hEtaPtJetPtGamma = new TH3F("hEtaPtJetPtGamma","hEtaPtJetPtGamma;p_{T,jet} (GeV);p_{T,part} (GeV)",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
  TH3F *hEtaPtJetZGamma = new TH3F("hEtaPtJetZGamma","hEtaPtJetZGamma;p_{T,jet} (GeV);z",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet,44,0,1.1);
  TH3F *hEtaPtJetPtKaon = new TH3F("hEtaPtJetPtKaon","hEtaPtJetPtKaon;p_{T,jet} (GeV);p_{T,part} (GeV)",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
  TH3F *hEtaPtJetPtProton = new TH3F("hEtaPtJetPtProton","hEtaPtJetPtProton;p_{T,jet} (GeV);p_{T,part} (GeV)",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
  TH3F *hEtaPtJetPtD = new TH3F("hEtaPtJetPtD","hEtaPtJetPtD;p_{T,jet} (GeV);p_{T,part} (GeV)",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
  TH3F *hEtaPtJetZD = new TH3F("hEtaPtJetZD","hEtaPtJetZD;p_{T,jet} (GeV);z",nEtaBins,-maxEtaJet,maxEtaJet,nPtJetBins,minPtJet,maxPtJet,12,0,1.2);

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

  AliPythia6 *pythia=AliPythia6::Instance();

  pythia->SetCKIN(3,pthard_min);   // minimum hard pt
  pythia->SetCKIN(4,pthard_max);  // maximum hard pt

  pythia->SetMDCY(pythia->Pycomp(111),1,0);  // switch off pi0 decay

  pythia->Initialize("CMS","p","p",e_cms);


  AliFastJetHeaderV1 *header = new AliFastJetHeaderV1;
  header->SetBGMode(0);
  //  header->SetRadius(0.4);
  header->SetRparam(jetR); 
  header->SetGhostEtaMax(maxEtaJet);
  header->SetJetEtaMin(-maxEtaJet);
  header->SetJetEtaMax(maxEtaJet);
  //header->SetGhostArea(0.05);
  header->SetAlgorithm(2); // antikt_algorithm = 2, kt = 0 (see fastjet/fastjet/JetDefinition.hh

  AliFastJetFinder *FastJet = new AliFastJetFinder;
  FastJet->SetJetHeader(header);

  AliAODEvent *aod = new AliAODEvent();
  aod->CreateStdContent();
  FastJet->ConnectAOD(aod);

  AliJetCalTrkEvent JetFinderEvent(0,1);
  TClonesArray *plist = new TClonesArray("TParticle");
  TClonesArray aliplist("AliMCParticle",1000);

  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet

    pythia->GenerateEvent();
    
    hNEvent->Fill(0.5);

    pythia->GetParticles(plist);

    aliplist.Clear();
    JetFinderEvent.Clear();

    Int_t n_part = plist->GetEntries();
    for (Int_t i_part = 0; i_part < n_part; i_part++) {
      part=(TParticle*) plist->At(i_part);

      if (abs(part->GetPdgCode()) == 421) {
	hPtD->Fill(part->Pt());
      }
      if (abs(part->GetPdgCode()) != 421 && part->GetStatusCode() >= 10)  // Not a final state particle (Keep Dzeros, if any)
	continue;

      if (chargedOnly) {
       if (!is_charged(part) && abs(part->GetPdgCode()) != 421) // skip neutrals if 'chargedOnly' is set, but keep D0
        continue;
      }
      else { // full jets
        if (part->GetPdgCode()==12 || part->GetPdgCode()==14 || part->GetPdgCode()==16 ) // reject neutrinos
          continue;
      }

      if (part->GetMother(0) >= 0 && fabs(((TParticle*)plist->At(part->GetMother(0)))->GetPdgCode()) == 421) // reject D0 daughters to prevent double counting
	continue;

      if (abs(part->GetPdgCode()) == 211) {
	hEtaPtPion->Fill(part->Eta(),part->Pt());
      }
      new (aliplist[i_part]) AliMCParticle(part);

      //else if (abs(part->GetPdgCode()) == 421) {

      //}
      // Only use charged tracks?

      JetFinderEvent.AddCalTrkTrackKine((AliMCParticle*)aliplist[i_part],1,1);
    }

    aod->ClearStd();
    FastJet->Reset();
    FastJet->SetCalTrkEvent(JetFinderEvent);
    FastJet->ProcessEvent();
    for (Int_t iJet = 0; iJet < aod->GetNJets(); iJet++) {
      jet = aod->GetJet(iJet);
      //cout << "\t jet " << iJet << " pt " << jet->Pt() << " eta " << jet->Eta() << " phi " << jet->Phi() << endl; 
      if (TMath::Abs(jet->Eta()) < maxEtaJet) {
	hJetPtEtaPhi->Fill(jet->Pt(),jet->Eta(),jet->Phi());
	// loop over constituents
	Int_t ntrk = jet->GetRefTracks()->GetEntries();
	for (Int_t itrk = 0; itrk < ntrk; itrk++) {
	  AliMCParticle *trk = dynamic_cast<AliMCParticle*>(jet->GetTrack(itrk));
	  hEtaPtJetPtCh->Fill(jet->Eta(), jet->Pt(), trk->Pt());

	  if (abs(trk->PdgCode()) == 22) {
	    hEtaPtJetPtGamma->Fill(jet->Eta(), jet->Pt(), trk->Pt());
	    hEtaPtJetZGamma->Fill(jet->Eta(), jet->Pt(), trk->Pt()/jet->Pt());
	  }
	  else if (abs(trk->PdgCode()) == 211) {
	    hEtaPtJetPtPion->Fill(jet->Eta(), jet->Pt(), trk->Pt());
	    hEtaPtJetZPion->Fill(jet->Eta(), jet->Pt(), trk->Pt()/jet->Pt());
	  }
	  else if (abs(trk->PdgCode()) == 2212) {
	    hEtaPtJetPtProton->Fill(jet->Eta(), jet->Pt(), trk->Pt());
	  }
	  else if (abs(trk->PdgCode()) == 321) {
	    hEtaPtJetPtKaon->Fill(jet->Eta(), jet->Pt(), trk->Pt());
	  }
	  else if (abs(trk->PdgCode()) == 421) {
	    hEtaPtJetPtD->Fill(jet->Eta(), jet->Pt(), trk->Pt());
	    hEtaPtJetZD->Fill(jet->Eta(), jet->Pt(), trk->Pt()/jet->Pt());
	  }
	}
      }
    }
  }

  hXSec->Fill(0.5,pythia->GetPARI(1));

  fout->Write();
}
