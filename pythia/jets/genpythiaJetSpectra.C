//
// similar to pythia6jetspectra, but now with aligenpythia
//

Int_t is_charged(TParticle *part) {
  Int_t abs_kf = abs(part->GetPdgCode());
  
  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf != 310 && abs_kf != 3122 && abs_kf != 3112  && abs_kf != 3222  && abs_kf != 3312 && abs_kf != 3322 && abs_kf != 3334 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16 && abs_kf != 3124)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

void genpythiaJetSpectra(Int_t nEvent = 50, Char_t const *foutname = "pythiaJetSpectra.root", Float_t pthard_min=10, Float_t pthard_max=100, Float_t e_cms = 7000) {

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

  gSystem->Load("libpythia6_4_25");
  //gSystem->Load("libpythia6_4_28");
  gSystem->Load("libEGPythia6");
  gSystem->Load("liblhapdf");
  gSystem->Load("libAliPythia6");


  const Float_t eta_jet_max = 0.5;
  const Float_t maxEta = 1.0; // particles
  const Int_t nEtaBins = 20;

  const Int_t chargedOnly = 0;  // charged only or full jets

  const Int_t nPtJetBins = 250;
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 250;
  TFile *fout = new TFile(foutname,"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; ;N",1,0,1);
  TH1F *hNTrials = new TH1F("hNTrials","number of trials; ;N",1,0,1);
  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);
  TH1F *hPtHard = new TH1F("hPtHard","pThard;p_{T,hard} (GeV)",nPtJetBins,minPtJet,maxPtJet);

  TH1F *hPtPion = new TH1F("hPtPion","charged pion pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  TH2F *hEtaPtCharged = new TH2F("hEtaPtCharged","charged particle pt spectrum;p_{T} (GeV)",nEtaBins,-maxEta,maxEta,nPtJetBins,minPtJet,maxPtJet);

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  const Float_t Rvals[]={0.2,0.3,0.4};
  //static const Int_t nR = sizeof(*Rvals)/sizeof(Float_t);
  const Int_t nR = 3;
  TH3F *hJetPtEtaPhi[nR] = {0};
  for (Int_t iR = 0; iR < nR; iR++) {
    hJetPtEtaPhi[iR] = new TH3F(Form("hJetPtEtaPhi_R%.0f",10*Rvals[iR]),Form("hJetPtEtaPhi R=%.1f;p_{T,jet} (GeV);#eta;#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,20,-eta_jet_max,eta_jet_max,32,0,2*TMath::Pi());
  }

  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

 AliGenPythia *pythia=new AliGenPythia(1);
  //pythia->SetProcess(kPyZgamma+1); //kPyJetsHardOnly); not sure why enum does not work..
  pythia->SetProcess(kPyJets);
  if (pthard_min > 0) {
    pythia->SetPtHard(pthard_min, pthard_max);
  }
  pythia->SetEnergyCMS(e_cms);
    //    Tune
    //        //    350     Perugia 11
    pythia->SetStrucFunc(kCTEQ5L);
    pythia->SetTune(350);
    pythia->UseNewMultipleInteractionsScenario();
    //                 
  pythia->Init();

  AliStack *stack = new AliStack();
  // There is a big mess with the connection between runloader, AliRun and the gAlice pointer. 
  // This order of things seems to work...
  AliRunLoader *dummyrl = new AliRunLoader("dummyEvtFolder");
  dummyrl->MakeHeader();
  dummyrl->SetEventNumber(0);
  gAlice->SetRunLoader(dummyrl);

  pythia->SetStack(stack);

  AliFastJetFinder *FastJet = new AliFastJetFinder;
  AliAODEvent *aod = new AliAODEvent();
  aod->CreateStdContent();
  FastJet->ConnectAOD(aod);

  AliJetCalTrkEvent JetFinderEvent(0,1);
  TClonesArray *plist = new TClonesArray("TParticle");
  TClonesArray aliplist("AliMCParticle",1000);

  AliFastJetHeaderV1 *fjheader = new AliFastJetHeaderV1;
  fjheader->SetBGMode(0);
  fjheader->SetAlgorithm(2); // antikt_algorithm = 2, kt = 0 (see fastjet/fastjet/JetDefinition.hh


  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    stack->Reset();
    TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet

    pythia->Generate();
    
    hNEvent->Fill(0.5);
    AliGenPythiaEventHeader *header = dynamic_cast<AliGenPythiaEventHeader*>(dummyrl->GetHeader()->GenEventHeader());
    hNTrials->Fill(0.5,header->Trials());
    hPtHard->Fill(header->GetPtHard());

    aliplist.Clear();
    JetFinderEvent.Clear();

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
      hEtaPtCharged->Fill(part->Eta(),part->Pt());
      if (TMath::Abs(part->Eta()) < maxEta) {
	new (aliplist[i_part]) AliMCParticle(part);

      //else if (abs(part->GetPdgCode()) == 421) {

      //}
      // Only use charged tracks?

	JetFinderEvent.AddCalTrkTrackKine((AliMCParticle*)aliplist[i_part],1,1);
      }
    }

    aod->ClearStd();
    FastJet->Reset();
    FastJet->SetCalTrkEvent(JetFinderEvent);
    for (Int_t iR = 0; iR < nR; iR++) {
      fjheader->SetRparam(Rvals[iR]); 
      //header->SetGhostEtaMax(2);
      //header->SetGhostArea(0.05);
      
      FastJet->SetJetHeader(fjheader);
      
      FastJet->ProcessEvent();
      for (Int_t iJet = 0; iJet < aod->GetNJets(); iJet++) {
	jet = aod->GetJet(iJet);
	//cout << "\t jet " << iJet << " pt " << jet->Pt() << " eta " << jet->Eta() << " phi " << jet->Phi() << endl; 
	if (TMath::Abs(jet->Eta()) < eta_jet_max) {
	  hJetPtEtaPhi[iR]->Fill(jet->Pt(),jet->Eta(),jet->Phi());
	}
      }
    }
  }
  
  hXSec->Fill(0.5,pythia->GetXsection());//GetPARI(1));

  fout->Write();
  TPythia6::Instance()->Pystat(0);
}
