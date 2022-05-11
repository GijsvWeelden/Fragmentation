#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"

#include "TPDGCode.h"
#include "TParticle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "THn.h"

#include <iostream>
#include <vector>

#include "AliPDG.h"
#include "AliRunLoader.h"
#include "AliStack.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

Int_t is_charged(TParticle *part) {
  Int_t abs_kf = abs(part->GetPdgCode());
  
  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf != 310 && abs_kf != 3122 && abs_kf != 3112  && abs_kf != 3222  && abs_kf != 3312 && abs_kf != 3322 && abs_kf != 3334 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16 && abs_kf != 3124)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

// Jet shapes from AliJetShape.h
float angularity(const fastjet::PseudoJet &jet) {
  if (!jet.has_constituents())
    return 0;
  Double_t den=0.;
  Double_t num = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    Double_t dphi = constits[ic].phi()-jet.phi();
    if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
    if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t dr2 = (constits[ic].eta()-jet.eta())*(constits[ic].eta()-jet.eta()) + dphi*dphi;
    Double_t dr = TMath::Sqrt(dr2);
    num=num+constits[ic].perp()*dr;
    den=den+constits[ic].perp();
  }
  return num/den;
}

float pTD(const fastjet::PseudoJet &jet) {
  if (!jet.has_constituents())
    return 0;
  Double_t den=0;
  Double_t num = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    num=num+constits[ic].perp()*constits[ic].perp();
    den=den+constits[ic].perp();
  }
  return TMath::Sqrt(num)/den;
}

float sigma2(const fastjet::PseudoJet &jet) {
  if (!jet.has_constituents())
    return 0;
  Double_t mxx    = 0.;
  Double_t myy    = 0.;
  Double_t mxy    = 0.;
  int  nc     = 0;
  Double_t sump2  = 0.;

  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    Double_t ppt=constits[ic].perp();
    Double_t dphi = constits[ic].phi()-jet.phi();
    if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
    if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t deta = constits[ic].eta()-jet.eta();
    mxx += ppt*ppt*deta*deta;
    myy += ppt*ppt*dphi*dphi;
    mxy -= ppt*ppt*deta*dphi;
    nc++;
    sump2 += ppt*ppt;
  }
  if(nc<2) return 0;
  if(sump2==0) return 0;
  // Sphericity Matrix
  Double_t ele[4] = {mxx , mxy, mxy, myy };
  TMatrixDSym m0(2,ele);

  // Find eigenvectors
  TMatrixDSymEigen m(m0);
  TVectorD eval(2);
  TMatrixD evecm = m.GetEigenVectors();
  eval  = m.GetEigenValues();
  // Largest eigenvector
  int jev = 0;
  if (eval[0] < eval[1]) jev = 1;
  TVectorD evec0(2);
  // Principle axis
  evec0 = TMatrixDColumn(evecm, jev);
  Double_t compx=evec0[0];
  Double_t compy=evec0[1];
  TVector2 evec(compx, compy);
  Double_t sigma2=0;
  if(jev==1) sigma2=TMath::Sqrt(TMath::Abs(eval[0])/sump2);
  if(jev==0) sigma2=TMath::Sqrt(TMath::Abs(eval[1])/sump2);
  return sigma2;
}

float LeSub(const fastjet::PseudoJet &jet) {
  if (!jet.has_constituents())
    return 0;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  std::vector<fastjet::PseudoJet> sortedconstits=sorted_by_pt(constits);
  if(sortedconstits.size()<2) return 0;
  Double_t num=TMath::Abs(sortedconstits[0].perp()-sortedconstits[1].perp());
  return num;
}

void analysePythia(Char_t const *indir = "output/events_40_80_50", Char_t const *foutname = "qPythiaJetSpectra.root") {

  const Float_t eta_jet_max = 1.0;

  const Int_t chargedOnly = 1;  // charged only or full jets
  // Constituents are charged only

  const Float_t pt_cut_rho = 1;

  const Int_t nPtJetBins = 250;
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 250;

  const Float_t maxEta = 3.0;

  TFile *fout = new TFile(foutname,"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);
 
  TH1F *hPtPion = new TH1F("hPtPion","charged pion pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  TH1F *hPtCharged = new TH1F("hPtCharged","charged particle pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  TH1F *hPtAll = new TH1F("hPtAll","all particle pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  const Float_t Rvals[]={0.2,0.3,0.4,0.5};
  //const Float_t Rvals[]={0.4,0.6,0.8,1.0};
  //static const Int_t nR = sizeof(*Rvals)/sizeof(Float_t);
  const Int_t nR = 4;
  TH3F *hJetPtEtaPhi[nR] = {0};
  TH3F *hJetPtEtaMass[nR] = {0};
  TH3F *hJetPtEtaAng[nR]={0};
  TH3F *hJetPtEtaPtD[nR]={0};

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
  
    hJetPtEtaAng[iR] = new TH3F(Form("hJetPtEtaAng_R%.0f",10*Rvals[iR]),Form("Angularity R=%.1f;p_{T,Lead_jet} (GeV/c);#eta;Angularity",Rvals[iR]),150,0,150,50,-2.5,2.5,50,0,0.5);

    hJetPtEtaPtD[iR] = new TH3F(Form("hJetPtEtaPtD_R%.0f",10*Rvals[iR]),Form("p_{T,D} R=%.1f;p_{T,Lead_jet} (GeV/c);#eta;p_{T,D} (GeV)",Rvals[iR]),150,0,150,50,-2.5,2.5,60,0,1.5);

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

  void *inDirP = gSystem->OpenDirectory(indir);
  const Char_t *subdir = 0;
  while (subdir = gSystem->GetDirEntry(inDirP)) {
    if (strcmp(subdir,".")==0 || strcmp(subdir,"..") == 0)
      continue;

    TString curDir(indir);
    curDir += "/";
    curDir += subdir;

    /*
    Int_t num = atoi(indir);
    if ( num <= 50 )
       continue;
    */
    cout << "input directory: " << curDir << endl;
    
    AliRunLoader *rl = AliRunLoader::Open(curDir+"/galice.root");  
    if (rl == 0)
      continue;
    
    rl->LoadKinematics();

    cout << "fout.ls() after LoadKinematics" << endl;
    fout->ls();

    TObjArray plistSel(1000);
    
    fastjet::GhostedAreaSpec ghostSpec(maxEta,1,0.01);
    fastjet::Strategy               strategy = fastjet::Best;
    //fastjet::RecombinationScheme    recombScheme = fastjet::BIpt_scheme;
    fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme; // need E scheme for jet mass
    fastjet::AreaType areaType =   fastjet::active_area_explicit_ghosts;
    fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);

    Int_t iEvent = 0;
    while (rl->GetEvent(iEvent)==0) {
      AliStack *stack = rl->Stack();

      Int_t index = 0;
      std::vector <fastjet::PseudoJet> fjInputs;
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

	fastjet::PseudoJet jInp(part->Px(),part->Py(),part->Pz(),part->Energy());  // need masses for E-scheme
	jInp.set_user_index(index);
	fjInputs.push_back(jInp);
	index++;

	n_trk++;
      }
      //cout << n_trk << " particles filled" << endl;
      //cout << JetFinderEvent.GetNCalTrkTracks() << " in event" << endl;

      for (Int_t iR = 0; iR < nR; iR++) {

       	fastjet::RangeDefinition range(-maxEta+Rvals[iR], maxEta-Rvals[iR], 0, 2.*fastjet::pi);

	fastjet::JetDefinition jetDefCh(fastjet::antikt_algorithm, Rvals[iR], recombScheme, strategy);
	fastjet::ClusterSequenceArea clustSeqCh(fjInputs, jetDefCh, areaDef);

	vector <fastjet::PseudoJet> inclusiveJetsCh = clustSeqCh.inclusive_jets();


	Float_t lead_pt = -1;
 
	Float_t lead_phi = -999;
	Float_t lead_eta = -999;
	for (UInt_t iJet = 0;  iJet < inclusiveJetsCh.size(); iJet++) {
	  fastjet::PseudoJet &jet = inclusiveJetsCh[iJet];
	  //cout << "\t jet " << iJet << " pt " << jet->Pt() << " eta " << jet.eta() << " phi " << jet.phi() << endl; 
	  if (TMath::Abs(jet.eta()) < eta_jet_max) {
	    hJetPtEtaPhi[iR]->Fill(jet.perp(),jet.eta(),jet.phi());
	    hJetPtEtaMass[iR]->Fill(jet.perp(),jet.eta(),jet.m());
	    hJetPtEtaAng[iR]->Fill(jet.perp(),jet.eta(),angularity(jet));
	    hJetPtEtaPtD[iR]->Fill(jet.perp(),jet.eta(),pTD(jet));
	    if (jet.perp() > lead_pt) {
	      lead_pt = jet.perp();
	      lead_eta = jet.eta();
	      lead_phi = jet.phi();
	    }
	  }
	}
	//cout << "lead_pt " << lead_pt << endl;
	if (lead_pt > 0) {
	  hLeadJetPtEtaPhi[iR]->Fill(lead_pt,lead_eta,lead_phi);
	  //cout << "aliplist entries: " << aliplist.GetEntriesFast() << " fast " << aliplist.GetEntriesFast() << " size " << aliplist.GetSize() << endl;
	  //Int_t n_trk2 = 0;
	  //Float_t sum_pt = 0;
	  /* 
	    // loop over constituents would need recoding
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
          }
	  */
	
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
    TH1D *hNEventIn = (TH1D*) fin.Get("hNEvent");
    hNEvent->Fill(0.5,hNEventIn->GetBinContent(1));  
    TH1D *hXSecIn = (TProfile*) fin.Get("hXSec");
    hXSec->Fill(0.5,hXSecIn->GetBinContent(1));  

    cout << "fout.ls() after filling xsec" << endl;
    fout->ls();

  }
  fout->cd();
  cout << "fout.ls() after everything" << endl;
  fout->ls();


  fout->Write();
}
