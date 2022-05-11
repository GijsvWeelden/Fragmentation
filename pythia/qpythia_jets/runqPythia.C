void runqPythia(Int_t nEvent = 50, Float_t pthard_min=10, Float_t pthard_max=100, Float_t qHat = 50., Float_t e_cms = 2760) {
  // NB: qHat is mean qhat in units 0.1 GeV^2/fm

  // example macro for on-the-fly generation of PYTHIA6 events
  // and analysis with new reader interface and FastJet
  // M. van Leeuwen

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  //gSystem->Load("libTENDER");
  gSystem->Load("libTender");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGTools");

  gSystem->Load("libqpythia");
  gSystem->Load("libEGPythia6");
  gSystem->Load("liblhapdf");
  gSystem->Load("libAliPythia6");


  const Float_t eta_jet_max = 1.0;

  const Int_t chargedOnly = 1;  // charged only or full jets

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

  AliGenPythia *pythia=new AliGenPythia(1);
  pythia->SetProcess(kPyJets);
  pythia->SetPtHard(pthard_min, pthard_max);
  pythia->SetEnergyCMS(e_cms);

  pythia->SetTune(103); //tune DW, standard choice for Q2 showers

  if (qHat != 0) { // Switch on quenching
    // Set qhat : <qhat> ~ k and <qhat>=1.7 GeV^2/fm for k=0.6*10^6 (default AliPythia)
    Float_t myk=qHat*0.6e6/17.;
    pythia->SetQuench(4); // q-PYTHIA
    pythia->SetQhat(myk);
  }

  pythia->Init();

  TFile *fout = new TFile("pythia_xsec.root","RECREATE");
  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);

  // This order of things seems to work...
  AliRunLoader *rl = rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"RECREATE");
  rl->MakeTree("E");
  rl->LoadKinematics("RECREATE");
  //AliStack *stack = rl->Stack();
  rl->SetNumberOfEventsPerFile(nEvent+1);
  gAlice->SetRunLoader(rl);


  TObjArray plistSel(1000);

  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    rl->SetEventNumber(iEvent);
    //gAlice->SetEventNrInRun(iEvent); // suppress Pylist in AliGenPythia ?
    rl->MakeHeader();
    if (rl->Stack())
      rl->Stack()->Reset();
    else 
      rl->MakeStack();

    AliStack *stack = rl->Stack();
    pythia->SetStack(stack);
    rl->MakeTree("K");
    //stack->ConnectTree(rl->TreeK());

    if ((iEvent%1000)==0) 
      cout << "Event " << iEvent << endl;

    pythia->Generate();
    
    //stack->PurifyKine();

    stack->FinishEvent();

    AliHeader *header = rl->GetHeader();
    header->SetNprimary(stack->GetNprimary());
    header->SetNtrack(stack->GetNtrack());

    TTree* treeE = rl->TreeE();
    rl->GetHeader()->SetStack(stack);
    treeE->Fill();

    hNEvent->Fill(0.5);
    rl->WriteKinematics("OVERWRITE");
  }
  AliRunLoader::Instance()->WriteHeader("OVERWRITE");

  pythia->FinishRun();
  // Write AliRun info and all detectors parameters
  AliRunLoader::Instance()->CdGAFile();
  gAlice->Write(0,TObject::kOverwrite);//write AliRun
  AliRunLoader::Instance()->Write(0,TObject::kOverwrite);//write RunLoader itself
  AliRunLoader::Instance()->Synchronize();
  //
  //rl->WriteRunLoader();

  hXSec->Fill(0.5,pythia->GetXsection());
  fout->Write();
}
