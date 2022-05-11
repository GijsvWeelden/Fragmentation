Int_t is_charged(TParticle *part) {
  Int_t abs_kf = abs(part->GetPdgCode());
  
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
  gSystem->Load("libfastjettools");
  //gSystem->Load("libfastjetcontribfragile");
  //gSystem->Load("libSISConePlugin");
  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");

  gSystem->Load("libpythia6.so");
  gSystem->Load("libEGPythia6.so");
  gSystem->Load("libAliPythia6.so");


  const Int_t n_eta_jet_bins = 28;
  const Float_t eta_jet_max = 0.7;
  const Float_t phi_emcal_min = 1.4;
  const Float_t phi_emcal_max = TMath::Pi();
  // EMCAL: 96x120 towers = 48x60 L0 patches (2x2 tower modules)
  const Int_t n_emcal_eta = 48;
  const Int_t n_emcal_phi = 60;

  // for geom, see: https://indico.cern.ch/event/343403/session/5/contribution/18/material/slides/0.pdf slide 9
  // 3 modules PHOS + DCAL
  const Float_t phi_dcal_min = 260./180*TMath::Pi();
  const Float_t phi_dcal_max = 320./180*TMath::Pi();
  const Float_t eta_phos_max = 0.12;
  const Float_t eta_dcal_min = 0.22;
  const Float_t eta_dcal_max = 0.70;
  // there are two 1/3 modules from DCAL: 320-327 degrees, with full length; those are ignored for now.
  // There is a half PHOS module starting 240-260 degrees; ignoring that for jets.
  //PHOS: 56x64 (etaxphi) crystals per SM
  //DCAL: TRU 8x12 modules (2x2 towers each) x 2= 16x12 on each side 
  
  const Int_t n_dcal_eta = 32;
  const Int_t n_dcal_phi = 36;

  const Int_t n_phos_eta = 6; // take roughly EMCAL tower size

  const Int_t nEtaPartBins = 40;
  const Float_t maxEtaPart = 2;

  const Int_t chargedOnly = 0;  // charged only or full jets

  const Int_t nPtJetBins = 100;
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 100;

  const Float_t hadRespFrac = 0.25;

  TH2F *hTrigTowEMCAL = new TH2F("hTrigTowEMCAL","trigger towers;#eta;#phi",n_emcal_eta,-0.7,0.7,n_emcal_phi,phi_emcal_min,phi_emcal_max);
  // eta phi axis are approximate; gap between PHOS and DCAL is taken into account when filling
  TH2F *hTrigTowDCAL = new TH2F("hTrigTowDCAL","trigger towers (eta,phi approximate);#eta;#phi",n_dcal_eta+n_phos_eta,-0.7,0.7,n_dcal_phi,phi_dcal_min,phi_dcal_max);

  // open file with E/p histos
  TFile *EopFile = new TFile("1112050815_lhc11d_V0_pionsEop_noEcut.root");
  TH2F *hEoverP2D = (TH2F*) EopFile->Get("hEopM"); // there is also hEop without M, probably inclusive?
  TObjArray *EopList = new TObjArray(100);

  TFile *fout = new TFile(foutname,"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);
 
  TH2F *hEtaPtPion = new TH2F("hEtaPtPion","charged pion pt spectrum;#eta;p_{T} (GeV)",nEtaPartBins,-maxEtaPart,maxEtaPart,nPtJetBins,minPtJet,maxPtJet);
  TH2F *hEtaPtCharged = new TH2F("hEtaPtCharged","charged particle pt spectrum;p_{T} (GeV)",nEtaPartBins,-maxEtaPart,maxEtaPart,nPtJetBins,minPtJet,maxPtJet);
  TH1F *hEtTrigEMCAL = new TH1F("hEtTrigEMCAL","trigger Et EMCAL;E_{T,trig} (GeV)",nPtJetBins*4,minPtJet,maxPtJet);
  TH1F *hEtTrigDCAL = new TH1F("hEtTrigDCAL","trigger Et DCAL;E_{T,trig} (GeV)",nPtJetBins*4,minPtJet,maxPtJet);

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  const Float_t Rvals[]={0.2,0.3,0.4};
  //static const Int_t nR = sizeof(*Rvals)/sizeof(Float_t);
  const Int_t nR = 3;
  TH3F *hJetPtEtaPhi[nR] = {0};
  TH3F *hJetEtaPtTrigPtEMCAL[nR] = {0};
  TH3F *hJetEtaPtTrigPtDCAL[nR] = {0};
  TH3F *hLeadJetEtaPtTrigPt[nR] = {0};
  //TH2F *hNearestJetPtTrigPt[nR] = {0};
  for (Int_t iR = 0; iR < nR; iR++) {
    hJetPtEtaPhi[iR] = new TH3F(Form("hJetPtEtaPhi_R%.0f",10*Rvals[iR]),Form("hJetPtEtaPhi R=%.1f;p_{T,jet} (GeV);#eta;#phi",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,n_eta_jet_bins,-eta_jet_max,eta_jet_max,32,0,2*TMath::Pi());
    hJetEtaPtTrigPtEMCAL[iR] = new TH3F(Form("hJetEtaPtTrigPtEMCAL_R%.0f",10*Rvals[iR]),Form("hJetEtaPtTrigPt EMCAL R=%.1f;#eta;p_{T,jet} (GeV);E_{t,trig} (GeV)",Rvals[iR]),n_eta_jet_bins,-eta_jet_max,eta_jet_max,nPtJetBins,minPtJet,maxPtJet,nPtJetBins*4,minPtJet,maxPtJet);
    hJetEtaPtTrigPtDCAL[iR] = new TH3F(Form("hJetEtaPtTrigPtDCAL_R%.0f",10*Rvals[iR]),Form("hJetEtaPtTrigPt DCAL R=%.1f;#eta;p_{T,jet} (GeV);E_{t,trig} (GeV)",Rvals[iR]),n_eta_jet_bins,-eta_jet_max,eta_jet_max,nPtJetBins,minPtJet,maxPtJet,nPtJetBins*4,minPtJet,maxPtJet);
    hLeadJetEtaPtTrigPt[iR] = new TH3F(Form("hLeadJetEtaPtTrigPt_R%.0f",10*Rvals[iR]),Form("hLeadJetEtaPtTrigPt R=%.1f;#eta;p_{T,jet} (GeV);E_{t,trig} (GeV)",Rvals[iR]),n_eta_jet_bins,-eta_jet_max,eta_jet_max,nPtJetBins,minPtJet,maxPtJet,nPtJetBins*4,minPtJet,maxPtJet);
    //hNearestJetPtTrigPt[iR] = new TH2F(Form("hNearestJetPtTrigPt_R%.0f",10*Rvals[iR]),Form("hNearestJetTrigPt R=%.1f;p_{T,jet} (GeV);E_{t,trig} (GeV)",Rvals[iR]),nPtJetBins,minPtJet,maxPtJet,nPtJetBins*4,minPtJet,maxPtJet);
  }

  // Create random number generator and set seed
  AliPythiaRndm::SetPythiaRandom(new TRandom3());
  AliPythiaRndm::GetPythiaRandom()->SetSeed(clock()+gSystem->GetPid());

  AliPythia6 *pythia=AliPythia6::Instance();

  pythia->SetCKIN(3,pthard_min);   // minimum hard pt
  pythia->SetCKIN(4,pthard_max);  // maximum hard pt

  //pythia->SetMDCY(pythia->Pycomp(111),1,0);  // switch off pi0 decay

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


  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    TProcessID::SetObjectCount(0); // Needed for TRefs in AliCalTrkTrack and AliAODJet

    pythia->GenerateEvent();
    
    hNEvent->Fill(0.5);

    pythia->GetParticles(plist);

    aliplist.Clear();
    JetFinderEvent.Clear();

    hTrigTowEMCAL->Reset();
    hTrigTowDCAL->Reset();

    Int_t n_part = plist->GetEntries();
    for (Int_t i_part = 0; i_part < n_part; i_part++) {
      part=(TParticle*) plist->At(i_part);

      if (part->GetStatusCode() >= 10)  // Not a final state particle 
	continue;

      if (part->GetPdgCode()==kPi0)
        cout << "Unexpected pi0; don't decay????" << endl;

      Int_t isEMCAL = -1;  // 0: PHOS/DCAL 1: EMCAL
      Int_t phi_bin = -1; // only used for DCAL+PHOS
      Int_t eta_bin = -1; 

      if (part->Phi() > phi_emcal_min && part->Phi() < phi_emcal_max) {
	isEMCAL = 1;
      }
      else if (part->Phi() > phi_dcal_min && part->Phi() < phi_dcal_max) {
	phi_bin = hTrigTowDCAL->GetYaxis()->FindBin(part->Phi());
	eta_bin = -1;
	if (TMath::Abs(part->Eta()) < eta_phos_max) {
	  eta_bin = 0.5*(n_phos_eta+n_dcal_eta) + part->Eta()/eta_phos_max*n_phos_eta;
	  if (part->Eta() > 0)
	    eta_bin++;
	}
	else if (TMath::Abs(part->Eta()) > eta_dcal_min && TMath::Abs(part->Eta()) < eta_dcal_max) {
	  eta_bin = (TMath::Abs(part->Eta())-eta_dcal_min)/(eta_dcal_max-eta_dcal_min)*(n_dcal_eta/2); 
	  if (part->Eta() > 0)
	    eta_bin += 0.5*(n_phos_eta+n_dcal_eta) + 1 + 0.5*n_phos_eta;
	  else
	    eta_bin = 0.5*(n_phos_eta+n_dcal_eta) - 0.5*n_phos_eta - eta_bin;
	  if (eta_bin < 1 || eta_bin > n_phos_eta+n_dcal_eta) 
	    cout << "Error in eta bin calculation; eta: " << part->Eta() <<
	      " bin: " << eta_bin << endl;
	}
	if (eta_bin > 0)
	  isEMCAL = 0;
      }
      
      if (TMath::Abs(part->GetPdgCode())==11 || part->GetPdgCode()==22) {
	// photon or electron
	//cout << "pdg " << part->GetPdgCode() << " parent " << part->GetMother(0) << endl;
	// Could do some gauss response
	if (isEMCAL==1)
	  hTrigTowEMCAL->Fill(part->Eta(), part->Phi(), part->Pt());       
	else if (isEMCAL==0) {
	  hTrigTowDCAL->AddBinContent(hTrigTowDCAL->GetBin(eta_bin,phi_bin), part->Pt());
	}
      }
      else {
	// sample E/p histo (using EMCAL response also for PHOS)
	Int_t ptbin = hEoverP2D->GetXaxis()->FindBin(part->Pt());
        if (ptbin == 0)
	  ptbin = 1;
	if (ptbin > hEoverP2D->GetXaxis()->GetNbins())
	  ptbin = hEoverP2D->GetXaxis()->GetNbins();
	if (EopList->At(ptbin-1)==0) {
          TH1 *hproj = hEoverP2D->ProjectionY(Form("EoverP_bin%d",ptbin),ptbin,ptbin);
          hproj->SetDirectory(0);
	  EopList->AddAt(hproj,ptbin-1);
	}
	// Need to check for overflows...
	TH1F *hEoverP1D = (TH1F*) EopList->At(ptbin-1);
	Float_t eoverp = hEoverP1D->GetRandom();
	if (isEMCAL==1)
	  hTrigTowEMCAL->Fill(part->Eta(), part->Phi(), eoverp*part->Pt());
	else if (isEMCAL==0)
	  hTrigTowDCAL->AddBinContent(hTrigTowDCAL->GetBin(eta_bin,phi_bin), eoverp*part->Pt());
      }
      if (chargedOnly) {
       if (!is_charged(part)) 
        continue;
      }
      else { // full jets
        if (part->GetPdgCode()==12 || part->GetPdgCode()==14 || part->GetPdgCode()==16 ) // reject neutrinos
          continue;
      }

      if (abs(part->GetPdgCode()) == 211) {
	hEtaPtPion->Fill(part->Eta(), part->Pt());
      }
      if (chargedOnly)
	hEtaPtCharged->Fill(part->Eta(), part->Pt());
      else if (is_charged(part))  // need to check charged
	hEtaPtCharged->Fill(part->Eta(), part->Pt());
      new (aliplist[i_part]) AliMCParticle(part);

      //else if (abs(part->GetPdgCode()) == 421) {

      //}
      // Only use charged tracks?

      JetFinderEvent.AddCalTrkTrackKine((AliMCParticle*)aliplist[i_part],1,1);
    }

    // Find trigger patch energy
    // trigger patch = 8x8 L0 patches // using 16x16 instead??
    Int_t maxx, maxy, maxz;
    hTrigTowEMCAL->GetMaximumBin(maxx,maxy,maxz);
    Int_t minbinx = TMath::Max(1,maxx - 8);
    Int_t minbiny = TMath::Max(1,maxy - 8);
    Int_t maxbinx = TMath::Min(hTrigTowEMCAL->GetXaxis()->GetNbins(),maxx + 7);
    Int_t maxbiny = TMath::Min(hTrigTowEMCAL->GetYaxis()->GetNbins(),maxy + 7);
    Float_t trig_eta = hTrigTowEMCAL->GetXaxis()->GetBinCenter(maxx);
    Float_t trig_phi = hTrigTowEMCAL->GetYaxis()->GetBinCenter(maxy);

    Float_t EtrigEMCAL = 0;
    for (Int_t ibinx = minbinx; ibinx <= maxbinx; ibinx++) {
      for (Int_t ibiny = minbiny; ibiny <= maxbiny; ibiny++) {
	EtrigEMCAL += hTrigTowEMCAL->GetBinContent(ibinx, ibiny);
      }
    }
    hEtTrigEMCAL->Fill(EtrigEMCAL);

    Int_t maxx, maxy, maxz;
    hTrigTowDCAL->GetMaximumBin(maxx,maxy,maxz);
    Int_t minbinx = TMath::Max(1,maxx - 8);
    Int_t minbiny = TMath::Max(1,maxy - 8);
    Int_t maxbinx = TMath::Min(hTrigTowDCAL->GetXaxis()->GetNbins(),maxx + 7);
    Int_t maxbiny = TMath::Min(hTrigTowDCAL->GetYaxis()->GetNbins(),maxy + 7);
    Float_t trig_eta = hTrigTowDCAL->GetXaxis()->GetBinCenter(maxx);
    Float_t trig_phi = hTrigTowDCAL->GetYaxis()->GetBinCenter(maxy);

    Float_t EtrigDCAL = 0;
    for (Int_t ibinx = minbinx; ibinx <= maxbinx; ibinx++) {
      for (Int_t ibiny = minbiny; ibiny <= maxbiny; ibiny++) {
	EtrigDCAL += hTrigTowDCAL->GetBinContent(ibinx, ibiny);
      }
    }
    hEtTrigDCAL->Fill(EtrigDCAL);

    aod->ClearStd();
    FastJet->Reset();
    FastJet->SetCalTrkEvent(JetFinderEvent);
    for (Int_t iR = 0; iR < nR; iR++) {
      header->SetRparam(Rvals[iR]); 
      //header->SetGhostEtaMax(2);
      //header->SetGhostArea(0.05);
      
      FastJet->SetJetHeader(header);
      
      FastJet->ProcessEvent();
      Float_t lead_pt = -1;
      Float_t lead_eta = -99;
      Float_t trig_dist = 100;
      Float_t nearest_pt = -1;
      for (Int_t iJet = 0; iJet < aod->GetNJets(); iJet++) {
	jet = aod->GetJet(iJet);
	//cout << "\t jet " << iJet << " pt " << jet->Pt() << " eta " << jet->Eta() << " phi " << jet->Phi() << endl;
	/*
        Float_t deta = jet->Eta()-trig_eta;
	if (TMath::Abs(deta) < trig_dist) {
          Float_t dphi = jet->Phi() - trig_phi;
          if (dphi > TMath::Pi())
            dphi -= 2*TMath::Pi();
          if (dphi < -TMath::Pi())
            dphi += 2*TMath::Pi();
          Float_t dist = TMath::Sqrt(deta*deta + dphi*dphi);
          if (dist < trig_dist) {
            trig_dist = dist;
            nearest_pt = jet->Pt();
          }
        } 
	*/
	// EMCAL trigger
	if (TMath::Abs(jet->Eta()) < eta_jet_max &&
            jet->Phi() > phi_emcal_min+Rvals[iR] && jet->Phi() < phi_emcal_max-Rvals[iR]) {
	  hJetPtEtaPhi[iR]->Fill(jet->Pt(),jet->Eta(),jet->Phi());
	  hJetEtaPtTrigPtEMCAL[iR]->Fill(jet->Eta(),jet->Pt(), EtrigEMCAL);

          if (jet->Pt() > lead_pt) {
            lead_pt = jet->Pt();
            lead_eta = jet->Eta();
          }
	}

	// DCAL trigger
	if (TMath::Abs(jet->Eta()) < eta_jet_max &&
            jet->Phi() > phi_dcal_min+Rvals[iR] && jet->Phi() < phi_dcal_max-Rvals[iR]) {
	  hJetPtEtaPhi[iR]->Fill(jet->Pt(),jet->Eta(),jet->Phi());
	  hJetEtaPtTrigPtDCAL[iR]->Fill(jet->Eta(), jet->Pt(), EtrigDCAL);

          if (jet->Pt() > lead_pt) {
            lead_pt = jet->Pt();
            lead_eta = jet->Eta();
          }
	}

      }
      hLeadJetEtaPtTrigPt[iR]->Fill(lead_eta, lead_pt, TMath::Max(EtrigDCAL,EtrigEMCAL));
      //hNearestJetPtTrigPt[iR]->Fill(nearest_pt, Etrig);
    }
  }
  
  hXSec->Fill(0.5,pythia->GetPARI(1));

  fout->Write();
}
