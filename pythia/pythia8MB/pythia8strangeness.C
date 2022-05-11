//
// NB: need to run with aliroot, not root (to have gAlice)
//
R__LOAD_LIBRARY(libpythia6)
R__LOAD_LIBRARY(libpythia8)
// R__LOAD_LIBRARY(libAliPythia8)
AliGenerator*  CreateGenerator();

void pythia8strangeness(Int_t nev = 100, const char* foutname = "strangenessHists_MB7TeV_temp.root", const char* galicename="galice.root")
{
//  Runloader
    gSystem->Setenv("PYTHIA8DATA",gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8210/xmldoc"));
    gRandom->SetSeed(clock()+gSystem->GetPid());
    AliRunLoader* rl = AliRunLoader::Open(galicename,"FASTRUN","recreate");
    
    rl->SetCompressionLevel(2);
    rl->SetNumberOfEventsPerFile(10000);
    //rl->LoadKinematics("RECREATE");
    rl->MakeTree("E");
    gAlice->SetRunLoader(rl);

//  Create stack
    rl->MakeStack();
    AliStack* stack      = rl->Stack();
 
//  Header
    AliHeader* header = rl->GetHeader();
//
//  Create and Initialize Generator
    AliGenerator *gener = CreateGenerator();
    gener->Init();
    //(AliPythia8::Instance())->PrintDecayTable();
    gener->SetStack(stack);
    
//
//                        Event Loop
//
    Int_t pdgCodes[] = {211, 321, 2212, 3122, 3312, 3334};
    const Char_t *pdgAbbrev[] = {"Pion","Kaon","proton","Lambda","Xi","Omega"};
    Int_t nPdg = sizeof(pdgCodes)/sizeof(pdgCodes[0]);

    TFile *fout = new TFile(foutname,"RECREATE");
    TH2F *NchFwdNchMid = new TH2F("NchFwdNchMid","multiplicity mid rap vs fwd;N_{ch}^{fwd};N_{ch}^{mid}",250,0,250,250,0,250);
    TH2F **NchMidPt = new TH2F*[nPdg];
    TH2F **NchFwdPt = new TH2F*[nPdg];
    for(Int_t iPdg = 0; iPdg < nPdg; iPdg++) {
      NchMidPt[iPdg] = new TH2F(Form("NchMidPt%s",pdgAbbrev[iPdg]),Form("Pt spectra vs mid rap. %s;N_{ch}^{mid};p_{T} (GeV/c)",pdgAbbrev[iPdg]),250,0,250,50,0,10);
      NchFwdPt[iPdg] = new TH2F(Form("NchFwdPt%s",pdgAbbrev[iPdg]),Form("Pt spectra vs fwd rap. %s;N_{ch}^{fwd};p_{T} (GeV/c)",pdgAbbrev[iPdg]),250,0,250,50,0,10);
    }
    Int_t iev;
     
    for (iev = 0; iev < nev; iev++) {

	//printf("\n \n Event number %d \n \n", iev);
	
//  Initialize event
	header->Reset(0,iev);
	rl->SetEventNumber(iev);
	stack->Reset();
	rl->MakeTree("K");
//	stack->ConnectTree();
    
//  Generate event
	gener->Generate();
//  Analysis
	Int_t npart = stack->GetNprimary();
	//printf("Analyse %d Particles\n", npart);
        Int_t NchMid = 0;
        Int_t NchFwd = 0;
        TDatabasePDG *dbpdg = TDatabasePDG::Instance();
	for (Int_t part=0; part<npart; part++) {
	    TParticle *MPart = stack->Particle(part);
	    Int_t mpart  = MPart->GetPdgCode();
            if (MPart->GetStatusCode() > 10)
               continue;

	    TParticlePDG *ppdg = dbpdg->GetParticle(mpart);
            if (!ppdg || ppdg->Charge() == 0)
	   	continue;

            Double_t peta = MPart->Eta();
            if (TMath::Abs(peta) < 0.9)
               NchMid++;
            else if (TMath::Abs(peta) > 2.5 && TMath::Abs(peta) < 4.0)
               NchFwd++;
//	    printf("Particle %d\n", mpart);
	}
        NchFwdNchMid->Fill(NchFwd, NchMid);
	//cout << "NchMid " << NchMid << " fwd " << NchFwd << endl;
	for (Int_t part=0; part<npart; part++) {
	  TParticle *MPart = stack->Particle(part);
	  Int_t mpart  = MPart->GetPdgCode();
	  if (TMath::Abs(MPart->Eta()) < 1.0) {
	    for(Int_t iPdg = 0; iPdg < nPdg; iPdg++) {   
	      if (TMath::Abs(mpart) == pdgCodes[iPdg]) {
                //cout << "Found " << pdgCodes[iPdg] << " status " << MPart->GetStatusCode() << endl;
		NchMidPt[iPdg]->Fill(NchMid,MPart->Pt());
		NchFwdPt[iPdg]->Fill(NchFwd,MPart->Pt());
	      }
            } 
	  }
	}
	
//  Finish event
	header->SetNprimary(stack->GetNprimary());
	header->SetNtrack(stack->GetNtrack());  

    } // event loop
//
//                         Termination
//  Generator
    fout ->Write();
    gener->FinishRun();
    
}


AliGenerator*  CreateGenerator()
{
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());

//
//
    //gener->SetProcess(kPyCharmppMNRwmi);
    //gener->SetForceDecay(kHadronicDWithout4Bodies);
    gener->SetProcess(kPyMbDefault);
//   Centre of mass energy 
    gener->SetEnergyCMS(13000.);
    //gener->SetEnergyCMS(7000.);
//   Initialize generator    
    gener->SetEventListRange(-1, 10);
    return gener;
}
