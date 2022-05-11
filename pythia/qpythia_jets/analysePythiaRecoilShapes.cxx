#include "fastjet/PseudoJet.hh"
#ifndef __CINT__
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#endif

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
//using std::numeric_limits; // Needed for Nsubjettiness?

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

void Nsubjettiness(const fastjet::PseudoJet &ParentJet, Double_t &oneSubjettiness, Double_t &twoSubjettiness, Double_t &Rdist, Double_t &SubJet1LeadingTrackPt, Double_t &SubJet2LeadingTrackPt){
  // From Nima:
  //Once you have the jet you want to perform the jet shape analysis on (Lets call it ParentJet) you can run the following code:

  Float_t Beta = 1.0; // A parameter of Nsubjettienss
  Float_t fR=0.4; //Jet resolution parameter
#ifndef __CINT__
  oneSubjettiness=-1;
  if (ParentJet.constituents().size() > 0){ //if jet has more than one track you can do a tau1 calculation on it
    fastjet::contrib::Nsubjettiness nSub_1(1, fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(Beta,fR));
    oneSubjettiness= nSub_1.result(ParentJet); //Tau1 Result   Note:Normally I do tau1 and tau2 in two separate calls to the function....so each time a copy of ParentJet is sent in...I dont think the .result function should change parentjet if it is a pointer but its good to keep an eye on
  }
  else oneSubjettiness = -2; //if it fails we fill the tree with -2


  twoSubjettiness=-1;
  if (ParentJet.constituents().size() > 1){ //if jet has more than two tracks we can do Tau2 and DeltaR
    std::vector<fastjet::PseudoJet> SubJet_Axes; //Axis of reclustering
    fastjet::PseudoJet SubJet1_Axis;
    fastjet::PseudoJet SubJet2_Axis;
    std::vector<fastjet::PseudoJet> SubJets;    //Reclustered Subjets
    fastjet::PseudoJet SubJet1;
    fastjet::PseudoJet SubJet2;
    fastjet::contrib::Nsubjettiness nSub_2(2, fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(Beta,fR));
    twoSubjettiness= nSub_2.result(ParentJet); //Tau2 Result
    SubJet_Axes=nSub_2.currentAxes(); //Reclustered Axes reesult
    SubJets=nSub_2.currentSubjets(); //Resclustered Subjets Result
    
    SubJet1_Axis=SubJet_Axes[0];	
    Double_t SubJet1_Eta = SubJet1_Axis.pseudorapidity();
    Double_t SubJet2_Eta = -999;
    Double_t SubJet1_Phi=SubJet1_Axis.phi();
    if(SubJet1_Phi < -1*TMath::Pi()) SubJet1_Phi += (2*TMath::Pi());
    else if (SubJet1_Phi > TMath::Pi()) SubJet1_Phi -= (2*TMath::Pi());
    Double_t SubJet2_Phi = -999;
    Double_t DeltaPhi=-5;  
    if (SubJet_Axes.size()>1){ //this should always be true
      SubJet2_Axis=SubJet_Axes[1];
      SubJet2_Eta=SubJet2_Axis.pseudorapidity();
      SubJet2_Phi=SubJet2_Axis.phi();
      if(SubJet2_Phi < -1*TMath::Pi()) SubJet2_Phi += (2*TMath::Pi());
      else if (SubJet2_Phi > TMath::Pi()) SubJet2_Phi -= (2*TMath::Pi());
      DeltaPhi=SubJet1_Phi-SubJet2_Phi; //Delta Phi between the two axis of reclustering
      if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
      else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
    }

    SubJet1=SubJets[0];
    SubJet1LeadingTrackPt=-3.0;
    SubJet2LeadingTrackPt=-3.0;
    std::vector<fastjet::PseudoJet> SubJet1Tracks = SubJet1.constituents();
    for (Int_t i=0; i<SubJet1Tracks.size(); i++){
      if (SubJet1Tracks[i].perp() > SubJet1LeadingTrackPt) SubJet1LeadingTrackPt=SubJet1Tracks[i].perp(); //Leading Track Pt in Subjet 1
    }
    if (SubJets.size()>1){ //should always be true as well
      SubJet2=SubJets[1];
      std::vector<fastjet::PseudoJet> SubJet2Tracks = SubJet2.constituents();
      for (Int_t i=0; i<SubJet2Tracks.size(); i++){
	if (SubJet2Tracks[i].perp() > SubJet2LeadingTrackPt) SubJet2LeadingTrackPt=SubJet2Tracks[i].perp(); //Leading Track  Pt in Subjet2
      }
    }

    if (ParentJet.constituents().size()>1 && SubJet_Axes.size()>1) 
      Rdist = TMath::Sqrt(TMath::Power(SubJet1_Eta-SubJet2_Eta,2)+TMath::Power(DeltaPhi,2)); //Delta R
    else
      Rdist = -2;
  }
  else twoSubjettiness=-2; //if there is less than 2 tracks in the jet return -2 as the answer

  //These are the 5 variables we would like to return at the end....and of course also the parentjet Pt
  
  /*
  if (ParentJet.constituents().size>1 && SubJets.size()>0)return SubJet1LeadingTrackPt; //we only care about the leading tracks if we can measure tau2 and delta R as well
  else return -2;
  if (ParentJet.constituents().size>1 && SubJets.size()>1) return SubJet2LeadingTrackPt;
  else return -2;
  */
#endif
}


void analysePythia(Char_t const *indir = "output/events_40_80_50", Char_t const *foutname = "qPythiaRecoilJetShapes.root", Int_t maxEvent = 200000) {

  const Float_t eta_jet_max = 1.0;

  const Int_t chargedOnly = 0;  // charged only or full jets
  // Constituents are charged only

  const Float_t pt_cut_rho = 1;

  const Int_t nPtJetBins = 250;
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 250;

  const Float_t maxEta = 0.9;
  const Float_t maxEtaJet = 0.5;
  const Float_t minPtTrig = 5;
  const Float_t maxPtTrig = 50;
  const Int_t nPtTrigBins = 45;
  const Float_t maxDPhi = 0.6;

  TFile *fout = new TFile(foutname,"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  TProfile *hXSec = new TProfile("hXSec","cross section; #sigma",1,0,1);
 
  TH1F *hPtPion = new TH1F("hPtPion","charged pion pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  TH1F *hPtCharged = new TH1F("hPtCharged","charged particle pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);
  TH1F *hPtAll = new TH1F("hPtAll","all particle pt spectrum;p_{T} (GeV)",nPtJetBins,minPtJet,maxPtJet);

  TH1F *hPtTrig = new TH1F("hPtTrig","trigger particle pt spectrum;p_{T,trig} (GeV)",nPtTrigBins,minPtTrig,maxPtTrig);

  AliPDG::AddParticlesToPdgDataBase(); // to add some PDF codes to TDatabasePDG

  const Float_t Rval = 0.4;

  TH3F *hTrigPtJetPtDPhi = 0;
  TH3F *hTrigPtJetPtMass = 0;
  TH3F *hTrigPtJetPtAng= 0;
  TH3F *hTrigPtJetPtPtD= 0;
  TH3F *hTrigPtJetPtTau1 = 0;
  TH3F *hTrigPtJetPtTau2 = 0;
  TH3F *hTrigPtJetPtTau2Tau1 = 0;
  TH3F *hTrigPtJetPtRsubdist = 0;

  hTrigPtJetPtDPhi = new TH3F(Form("hTrigPtJetPtDPhi_R%.0f",10*Rval),Form("hTrigPtJetPtDPhi R=%.1f;p_{T,trig} (GeV);p_{T,jet} (GeV);#Delta#phi",Rval),nPtTrigBins,minPtTrig,maxPtTrig,nPtJetBins,minPtJet,maxPtJet,48,0,TMath::Pi());

  hTrigPtJetPtMass = new TH3F(Form("hTrigPtJetPtMass_R%.0f",10*Rval),Form("hTrigPtJetPtMass R=%.1f;p_{T,trig} (GeV);p_{T,jet} (GeV);jet mass (GeV/c^{2})",Rval),nPtTrigBins,minPtTrig,maxPtTrig,nPtJetBins,minPtJet,maxPtJet,200,0,50);
  
  hTrigPtJetPtAng = new TH3F(Form("hTrigPtJetPtAng_R%.0f",10*Rval),Form("Angularity R=%.1f;p_{T,Lead_jet} (GeV/c);Angularity",Rval),nPtTrigBins,minPtTrig,maxPtTrig,nPtJetBins,minPtJet,maxPtJet,50,0,0.5);

  hTrigPtJetPtPtD = new TH3F(Form("hTrigPtJetPtPtD_R%.0f",10*Rval),Form("p_{T,D} R=%.1f;p_{T,Lead_jet} (GeV/c);p_{T,D} (GeV)",Rval),nPtTrigBins,minPtTrig,maxPtTrig,nPtJetBins,minPtJet,maxPtJet,60,0,1.5);

  hTrigPtJetPtTau1 = new TH3F(Form("hTrigPtJetPtTau1_R%.0f",10*Rval),Form("p_{T,D} R=%.1f;p_{T,Lead_jet} (GeV/c);#tau_{1}",Rval),nPtTrigBins,minPtTrig,maxPtTrig,nPtJetBins,minPtJet,maxPtJet,60,0,1.5);

  hTrigPtJetPtTau2 = new TH3F(Form("hTrigPtJetPtTau2_R%.0f",10*Rval),Form("p_{T,D} R=%.1f;p_{T,Lead_jet} (GeV/c);#tau_{2}",Rval),nPtTrigBins,minPtTrig,maxPtTrig,nPtJetBins,minPtJet,maxPtJet,60,0,1.5);

  hTrigPtJetPtTau2Tau1 = new TH3F(Form("hTrigPtJetPtTau2Tau1_R%.0f",10*Rval),Form("#tau_{2}/#tau_{1} R=%.1f;p_{T,Lead_jet} (GeV/c);#tau_{2}/#tau_{1}",Rval),nPtTrigBins,minPtTrig,maxPtTrig,nPtJetBins,minPtJet,maxPtJet,60,0,1.5);

  hTrigPtJetPtRsubdist = new TH3F(Form("hTrigPtJetPtRsubdist_R%.0f",10*Rval),Form("#DeltaR_{subjets} R=%.1f;p_{T,Lead_jet} (GeV/c);#DeltaR_{subjets}",Rval),nPtTrigBins,minPtTrig,maxPtTrig,nPtJetBins,minPtJet,maxPtJet,60,0,0.6);


  cout << "fout.ls() after creating histos" << endl;
  fout->ls();

  void *inDirP = gSystem->OpenDirectory(indir);
  const Char_t *subdir = 0;
  Int_t iEventTot = 0;
  while ((subdir = gSystem->GetDirEntry(inDirP)) && iEventTot < maxEvent) {
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
    // Changed dd 20170110: we had 'explicit ghosts' and switched them off
    // This impacts the Nsubjettiness code and also the jet population (pure ghost jets). Need to check on old results as well...
    // Explicit ghosts only needed when subtracting bkg in shapes
    //fastjet::AreaType areaType =   fastjet::active_area_explicit_ghosts;
    fastjet::AreaType areaType =   fastjet::active_area;
    fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);

    Int_t iEvent = 0;
    while (rl->GetEvent(iEvent)==0 && iEventTot < maxEvent) {
      AliStack *stack = rl->Stack();

      Int_t index = 0;
      std::vector <fastjet::PseudoJet> fjInputs;
      Int_t n_trk = 0;
      Int_t n_part = stack->GetNtrack();
      Int_t n_trig = 0;
      TObjArray triggers;
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
	if (TMath::Abs(part->Eta()) > maxEta)
	  continue;
	if (abs(part->GetPdgCode()) == 211) {
	  hPtPion->Fill(part->Pt());
	}
	hPtAll->Fill(part->Pt());
	if (is_charged(part)) {
	  hPtCharged->Fill(part->Pt());
	  if (part->Pt() > minPtTrig) {
	    n_trig++;
	    triggers.AddLast(part);
	  }
	}

	fastjet::PseudoJet jInp(part->Px(),part->Py(),part->Pz(),part->Energy());  // need masses for E-scheme
	jInp.set_user_index(index);
	fjInputs.push_back(jInp);
	index++;

	n_trk++;
      }
      //cout << n_trk << " particles filled" << endl;
      //cout << JetFinderEvent.GetNCalTrkTracks() << " in event" << endl;

      if (n_trig > 0) {
	fastjet::RangeDefinition range(-maxEta+Rval, maxEta-Rval, 0, 2.*fastjet::pi);

	fastjet::JetDefinition jetDefCh(fastjet::antikt_algorithm, Rval, recombScheme, strategy);
	fastjet::ClusterSequenceArea clustSeqCh(fjInputs, jetDefCh, areaDef);

	vector <fastjet::PseudoJet> inclusiveJetsCh = clustSeqCh.inclusive_jets();

	for (Int_t iTrig = 0; iTrig < triggers.GetEntriesFast(); iTrig++) { 
	  TParticle *trigTrack = (TParticle*) triggers[iTrig];
          hPtTrig->Fill(trigTrack->Pt());
	  for (UInt_t iJet = 0;  iJet < inclusiveJetsCh.size(); iJet++) {
	    fastjet::PseudoJet &jet = inclusiveJetsCh[iJet];
	    //cout << "\t jet " << iJet << " pt " << jet->Pt() << " eta " << jet.eta() << " phi " << jet.phi() << endl; 
	    if (TMath::Abs(jet.eta()) < eta_jet_max) {
	      Float_t dphi = jet.phi() - trigTrack->Phi();
	      while (dphi < 0)
		dphi += 2*TMath::Pi();
	      while (dphi > 2*TMath::Pi())
		dphi -= 2*TMath::Pi();
              if (dphi > TMath::Pi())
                dphi = TMath::Pi() - dphi;
	      hTrigPtJetPtDPhi->Fill(trigTrack->Pt(),jet.perp(),dphi);
	      if (TMath::Abs(TMath::Pi() - dphi) < maxDPhi) {
		hTrigPtJetPtMass->Fill(trigTrack->Pt(),jet.perp(),jet.m());
		hTrigPtJetPtAng->Fill(trigTrack->Pt(),jet.perp(),angularity(jet));
		hTrigPtJetPtPtD->Fill(trigTrack->Pt(),jet.perp(),pTD(jet));
		if (fabs(Rval-0.4)<0.001) {
		  Double_t tau1, tau2, Rdistsubjets, ptsub1, ptsub2;
		  Nsubjettiness(jet, tau1, tau2, Rdistsubjets, ptsub1, ptsub2);
		  /*
		    if (tau1 < 1e-10 || tau2 < 1e-10) {
		    cout << "jet pt " << jet.perp() << " nconst " << jet.constituents().size() << " : tau1 " << tau1 << ", tau2 " << tau2 << " distance " << Rdistsubjets << endl;
		    }
		  */
		  hTrigPtJetPtTau1->Fill(trigTrack->Pt(),jet.perp(),tau1);
		  hTrigPtJetPtTau2->Fill(trigTrack->Pt(),jet.perp(),tau2);
		  if (tau1 > 1e-10 && tau2 > 1e-10) {
		    hTrigPtJetPtTau2Tau1->Fill(trigTrack->Pt(),jet.perp(),tau2/tau1);
		  }
		  hTrigPtJetPtRsubdist->Fill(trigTrack->Pt(),jet.perp(),Rdistsubjets);
		}
	      }
	    }
	  }
	}
      }
      iEvent++;
      iEventTot++;
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
