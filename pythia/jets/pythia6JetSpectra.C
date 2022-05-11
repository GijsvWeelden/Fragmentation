Int_t is_charged(TParticle *part) {
  Int_t abs_kf = abs(part->GetPdgCode());
  
  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16 && abs_kf != 3124)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

Int_t is_charged(AliMCParticle *part) {
  Int_t abs_kf = abs(part->PdgCode());
  
  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16 && abs_kf != 3124)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

void pythia6JetSpectra(Int_t nEvent = 50, Char_t const *foutname = "pythiaJetSpectra.root", Float_t pthard_min=10, Float_t pthard_max=100, Float_t e_cms = 7000) {

  // example macro for on-the-fly generation of PYTHIA6 events
  // and analysis with new reader interface and FastJet
  // M. van Leeuwen

  gSystem->AddDynamicPath("/usr/lib64");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libsiscone");
  gSystem->Load("libsiscone_spherical");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");

  gSystem->Load("libpythia6.so");
  gSystem->Load("libEGPythia6.so");
  gSystem->Load("libAliPythia6.so");


  const Float_t eta_jet_max = 1.0;
  const Float_t maxEtaPart = 1.0; // particles
  const Int_t nEtaBins = 20;

  const Int_t chargedOnly = 0;  // charged only or full jets

  const Int_t nPtJetBins = 250;
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 250;
  const Int_t nEtaBins = 20;
  const Float_t minEta = -1.0;
  const Float_t maxEta = 1.0;
  TFile *fout = new TFile(foutname,"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);
  TH1F *hPtHard = new TH1F("hPtHard","pThard;p_{T,hard} (GeV)",nPtJetBins,minPtJet,maxPtJet);

  TH1F *hPtPion = new TH1F("hPtPion","charged pion pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  TH2F *hEtaPtCharged = new TH2F("hEtaPtCharged","charged particle pt spectrum;p_{T} (GeV)",nEtaBins,minEta,maxEta,nPtJetBins,minPtJet,maxPtJet);

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  const Float_t Rvals[]={0.2,0.3,0.4};
  //static const Int_t nR = sizeof(*Rvals)/sizeof(Float_t);
  const Int_t nR = 3;
  TH3F *hJetPtEtaPhi[nR] = {0};
  TH3F *hLeadJetPtEtaPhi[nR] = {0};
  TH2F *hLeadJetPtJetPt[nR] = {0};
  TH2F *hLeadJetPtChargedPt[nR] = {0};
  TH2F *hLeadJetPtPi0Pt[nR] = {0};
  for (Int_t iR = 0; iR < nR; iR++) {
    hJetPtEtaPhi[iR] = new TH3F(Form("hJetPtEtaPhi_R%.0f",10*Rvals[iR]),Form("hJetPtEtaPhi R=%.1f;p_{T,jet} (GeV);#eta;#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,32,0,2*TMath::Pi());
    hLeadJetPtEtaPhi[iR] = new TH3F(Form("hLeadJetPtEtaPhi_R%.0f",10*Rvals[iR]),Form("hJetPtEtaPhi R=%.1f;p_{T,lead jet} (GeV);#eta;#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,32,0,2*TMath::Pi());
    hLeadJetPtJetPt[iR] = new TH2F(Form("hLeadJetPtJetPt_R%.0f",10*Rvals[iR]),Form("hLeadJetPtJetPt R=%.1f;p_{T,lead jet} (GeV);p_{T, jet}",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
    hLeadJetPtChargedPt[iR] = new TH2F(Form("hLeadJetPtChargedPt_R%.0f",10*Rvals[iR]),Form("hLeadJetPtChargedPt R=%.1f;p_{T,lead jet} (GeV);p_{T}",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
    hLeadJetPtPi0Pt[iR] = new TH2F(Form("hLeadJetPtPi0Pt_R%.0f",10*Rvals[iR]),Form("hLeadJetPtPi0Pt R=%.1f;p_{T,lead jet} (GeV);p_{T, #pi^{0}}",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,nPtJetBins,minPtJet,maxPtJet);
  }

  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

  AliPythia6 *pythia=AliPythia6::Instance();

  if (pthard_min > 0) {
    pythia->SetCKIN(3,pthard_min);  // minimum hard pt
    pythia->SetCKIN(4,pthard_max);  // maximum hard pt
  }

  pythia->SetMDCY(pythia->Pycomp(111),1,0);  // switch off pi0 decay

  pythia->Initialize("CMS","p","p",e_cms);

  AliFastJetFinder *FastJet = new AliFastJetFinder;
  AliAODEvent *aod = new AliAODEvent();
  aod->CreateStdContent();
  FastJet->ConnectAOD(aod);

  AliJetCalTrkEvent JetFinderEvent(0,1);
  TClonesArray *plist = new TClonesArray("TParticle");
  TClonesArray aliplist("AliMCParticle",1000);

  AliFastJetHeaderV1 *header = new AliFastJetHeaderV1;
  header->SetBGMode(0);
  header->SetAlgorithm(2); // antikt_algorithm = 2, kt = 0 (see fastjet/fastjet/JetDefinition.hh
  header->SetJetEtaMin(-maxEtaPart);
  header->SetJetEtaMax(maxEtaPart);
  header->SetPtMin(0);

  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet

    pythia->GenerateEvent();
    
    hNEvent->Fill(0.5);
    hPtHard->Fill(pythia->GetPtHard());

    pythia->GetParticles(plist);

    aliplist.Clear();
    JetFinderEvent.Clear();

    Int_t n_part = plist->GetEntries();
    Int_t i_part_sel = 0;
    for (Int_t i_part = 0; i_part < n_part; i_part++) {
      part=(TParticle*) plist->At(i_part);

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
      if (chargedOnly || is_charged(part))
	hEtaPtCharged->Fill(part->Eta(),part->Pt());
      if (TMath::Abs(part->Eta()) > maxEtaPart)
        continue;
      new (aliplist[i_part_sel]) AliMCParticle(part);

      //else if (abs(part->GetPdgCode()) == 421) {

      //}
      // Only use charged tracks?

      JetFinderEvent.AddCalTrkTrackKine((AliMCParticle*)aliplist[i_part_sel],1,1);
      i_part_sel++;
    }

    aod->ClearStd();
    FastJet->Reset();
    FastJet->SetCalTrkEvent(JetFinderEvent);
    for (Int_t iR = 0; iR < nR; iR++) {
      header->SetRparam(Rvals[iR]); 
      //header->SetGhostEtaMax(2);
      //header->SetGhostArea(0.05);
      
      FastJet->SetJetHeader(header);
      
      FastJet->ProcessEvent();
      Float_t ptLead = 0, phiLead = 0, etaLead = 0;
      for (Int_t iJet = 0; iJet < aod->GetNJets(); iJet++) {
	jet = aod->GetJet(iJet);
	//cout << "\t jet " << iJet << " pt " << jet->Pt() << " eta " << jet->Eta() << " phi " << jet->Phi() << endl; 
	if (TMath::Abs(jet->Eta()) < eta_jet_max) {
	  hJetPtEtaPhi[iR]->Fill(jet->Pt(),jet->Eta(),jet->Phi());
          if (jet->Pt() > ptLead) {
            ptLead = jet->Pt();
            phiLead = jet->Phi();
            etaLead = jet->Eta();
          }
	}
        else {
          cout << "Found jet outside acceptance, eta " << jet->Eta() << endl;
        }
      }
      hLeadJetPtEtaPhi[iR]->Fill(ptLead, phiLead, etaLead);
      if (ptLead == 0)
        cout << "No leading jet; total jets in event: " << aod->GetNJets() << endl;
      for (Int_t iJet = 0; iJet < aod->GetNJets(); iJet++) {
	jet = aod->GetJet(iJet);
	if (TMath::Abs(jet->Eta()) < eta_jet_max) {
          hLeadJetPtJetPt[iR]->Fill(ptLead,jet->Pt());
        }
      } 
      Int_t nPart = aliplist.GetEntriesFast();
      for (Int_t iPart = 0; iPart < nPart; iPart++) {
         AliMCParticle *apart = (AliMCParticle*) aliplist[iPart];
         if (chargedOnly || is_charged(apart)) {
           if (apart->Pt()>10 && apart->Pt() > 1.1*ptLead) {
             cout << "found track with pt > ptLead: " << apart->Pt()  << " eta " << apart->Eta() << " phi " << apart->Phi() << " ptLead " << ptLead << " etaLead " << etaLead << " phiLead " << phiLead << endl;
             cout << "Jet list" << endl;
             for (Int_t iJet = 0; iJet < aod->GetNJets(); iJet++) {
	       jet = aod->GetJet(iJet);
	       cout << "\t jet " << iJet << " pt " << jet->Pt() << " eta " << jet->Eta() << " phi " << jet->Phi() << endl; 
             }
           }
           hLeadJetPtChargedPt[iR]->Fill(ptLead, apart->Pt());
         }
         if (apart->PdgCode() == 112) {
           cout << "Fill Pi0" << endl;
           hLeadJetPtPi0Pt[iR]->Fill(ptLead,apart->Pt());
         }
      }
    }
  }
  hXSec->Fill(0.5,pythia->GetPARI(1));

  fout->Write();
  pythia->Pystat(1);
}
