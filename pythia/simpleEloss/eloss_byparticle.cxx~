#include "fastjet/Selector.hh" //.......... Background Sutraction event by event
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"//.......... Background Sutraction event by event
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/Recluster.hh"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TList.h"
#include "TVector3.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TNtuple.h"
#include "TString.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"  
#include <ctime>
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <valarray>

#define nEvents 100000

using namespace Pythia8;

//_________________________________________________________________________ 
Double_t RelativePhi(Double_t mphi,Double_t vphi){
   //Get relative azimuthal angle of two particles -pi to pi
   if      (vphi < -TMath::Pi()) vphi += TMath::TwoPi();
   else if (vphi > TMath::Pi())  vphi -= TMath::TwoPi();

   if      (mphi < -TMath::Pi()) mphi += TMath::TwoPi();
   else if (mphi > TMath::Pi())  mphi -= TMath::TwoPi();

   Double_t dphi = mphi - vphi;
   if      (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
   else if (dphi > TMath::Pi())  dphi -= TMath::TwoPi();

   return dphi;//dphi in [-Pi, Pi]
}
//__________________________________________________________________________

double getDeltaE(double Ejet) {
  return gRandom->Exp(0.15*Ejet);
}


//__________________________________________________________________________


int main(int , char**)
{

  //double jetRloss = 0.8;
  double minPtJetLoss = 10;

  double jetR   = 0.4; //jet R
  double trackEtaCut     = 3;
  double trackLowPtCut   = 0.15; //GeV
//__________________________________________________________________________
//PYTHIA SETTINGS

  TString name;
  
  int mecorr=1;
  
  Float_t ptHatMin=20;
  Float_t ptHatMax=200;

  gRandom->SetSeed(clock() + gSystem->GetPid());
  
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:idA = 2212"); //beam 1 proton
	pythia.readString("Beams:idB = 2212"); //beam 2 proton
	pythia.readString("Beams:eCM = 13000.");
	pythia.readString("Tune:pp = 5");  //tune 1-13    5=defaulr TUNE4C,  6=Tune 4Cx, 7=ATLAS MB Tune A2-CTEQ6L1

	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = 0");

  
	pythia.readString("HardQCD:all = on");
	if(ptHatMin<0 || ptHatMax <0){     
	pythia.readString("PhaseSpace:pTHatMin = 0."); // <<<<<<<<<<<<<<<<<<<<<<<
	}else{
	name = Form("PhaseSpace:pTHatMin = %f", (Float_t) ptHatMin);
	pythia.readString(name.Data()); 
	name = Form("PhaseSpace:pTHatMax = %f", (Float_t) ptHatMax);
	pythia.readString(name.Data()); 
	}

	//if(nEvents==0){
	//pythia.readString("PartonLevel:MPI = off");
	//pythia.readString("PartonLevel:ISR = off");
	//}
  
   
	pythia.readString("310:mayDecay  = off"); //K0s
	pythia.readString("3122:mayDecay = off"); //labda0
	pythia.readString("3112:mayDecay = off"); //sigma-
	pythia.readString("3212:mayDecay = off"); //sigma0
	pythia.readString("3222:mayDecay = off"); //sigma+
	pythia.readString("3312:mayDecay = off"); //xi-
	pythia.readString("3322:mayDecay = off"); //xi+
	pythia.readString("3334:mayDecay = off"); //omega-

	//ME corrections
	//use of matrix corrections where available
	if(mecorr==0){ 
	pythia.readString("TimeShower:MECorrections=off");
	}
	pythia.init();
 
//_________________________________________________________________________________________________
//FASTJET  SETTINGS

	double etaminJet = - trackEtaCut + jetR; //signal jet eta range
	double etamaxJet = - etaminJet;


	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recombScheme = fastjet::E_scheme;

	fastjet::JetDefinition *jetDefAKT_Sig = NULL;

	jetDefAKT_Sig = new fastjet::JetDefinition(fastjet::antikt_algorithm, jetR, recombScheme, strategy);
 
        /*
        // This code is not used
	fastjet::GhostedAreaSpec ghostareaspec(trackEtaCut, 1, 0.05); //ghost 
	//max rap, repeat, ghostarea default 0.01
	fastjet::AreaType areaType = fastjet::active_area_explicit_ghosts;
	fastjet::AreaDefinition *areaDef = new fastjet::AreaDefinition(areaType, ghostareaspec);
        */
	// Fastjet input
	std::vector<fastjet::PseudoJet> fjInputs;

//___________________________________________________ 
//HISTOGRAMS


	TFile* outFile = new TFile("Output_jetbyjet_relDE015_expo_expo.root","RECREATE");

	// Initial distributions
	TH2D *hJetEtaPt = new TH2D("hJetEtaPt","eta pt of jets;#eta;p_{T} (GeV/c)",60,-3,3,250,0,500);

	hJetEtaPt->Sumw2();
	TH2D *hTrackEtaPt = new TH2D("hTrackEtaPt","et pt of tracks;#eta;p_{T} (GeV/c)",60,-3,3,250,0,500);

	hTrackEtaPt->Sumw2();

	/*
	TH1D *hJet_deltaR_A = new TH1D("hJet_deltaR_A", "hJet_deltaR_A;#DeltaR", 100, 0.0, 1.0);
	hJet_deltaR_A->Sumw2();

	TH1D *hJet_deltaRg_A = new TH1D("hJet_deltaRg_A", "hJet_deltaRg_A;#DeltaR_{g}", 100, 0.0, 1.0);
	hJet_deltaRg_A->Sumw2();
        */

	// Distributions after energy loss
	TH2D *hJetEtaPtLoss = new TH2D("hJetEtaPtLoss","eta pt of jets after loss;#eta;p_{T} (GeV/c)",60,-3,3,250,0,500);

	hJetEtaPtLoss->Sumw2();
	TH2D *hTrackEtaPtLoss = new TH2D("hTrackEtaPtLoss","et pt of tracks after loss;#eta;p_{T} (GeV/c)",60,-3,3,250,0,500);

	hTrackEtaPtLoss->Sumw2();


//Begin event loop

	for(int i = 0; i < nEvents; i++)
	{
		double fourvec[4];

		if(!pythia.next()) continue;
		
		fjInputs.resize(0);
                int nPart = 0;

		for(int j = 0; j < pythia.event.size(); j++)
		{
			if(pythia.event[j].isFinal())
			{
				//Apply cuts in the particles
				if(pythia.event[j].pT() < trackLowPtCut) continue;                 //pt cut
	     			if(TMath::Abs(pythia.event[j].eta()) > trackEtaCut) continue;      //eta cut


				fourvec[0]=pythia.event[j].px();
             			fourvec[1]=pythia.event[j].py();
             			fourvec[2]=pythia.event[j].pz();
             			fourvec[3]=pythia.event[j].e();

				fastjet::PseudoJet PythiaParticle(fourvec);
                                PythiaParticle.set_user_index(nPart);
				fjInputs.push_back(PythiaParticle);
                                nPart++;
			}
		}

		//Jet Reconstruction:
		std::vector<fastjet::PseudoJet> PythiaJets;//Declaration of vector for Reconstructed Jets

		fastjet::GhostedAreaSpec ghost_spec(1, 1, 0.05);//Ghosts to calculate the Jet Area

		fastjet::AreaDefinition fAreaDef(fastjet::passive_area,ghost_spec);//Area Definition

		fastjet::ClusterSequenceArea clustSeq1(fjInputs, *jetDefAKT_Sig, fAreaDef);//Cluster Sequence

		PythiaJets = sorted_by_pt(clustSeq1.inclusive_jets(1.));//Vector with the Reconstructed Jets in pT order
		//__

		if(PythiaJets.size()==0) continue;

		if(PythiaJets[0].pt() < minPtJetLoss) continue;

		// Fill track histo
		for (auto part : fjInputs) {
		  hTrackEtaPt->Fill(part.pseudorapidity(),part.pt());
		}
                // Loop over jets
                for (auto jet : PythiaJets) {
		  hJetEtaPt->Fill(jet.pseudorapidity(), jet.pt());
                  if (jet.pt() > minPtJetLoss) {
		    double eLossFact = 1. - getDeltaE(jet.pt())/jet.pt();
                    if (eLossFact <= 0.1) { // Limit energy loss to 90 %
                      eLossFact = 0.1;
                      // cout << "ERROR: energy loss > jet energy; not sure what todo" << endl;
                    } 
		    for (auto part : jet.constituents()) {
                      double elossFrac = 1. - eLossFact;
                      elossFrac = gRandom->Exp(elossFrac);
                      if (elossFrac > 0.95)
                        elossFrac = 0.95; 
		      fjInputs[part.user_index()].reset_momentum((1.-elossFrac)*part.px(),(1.-elossFrac)*part.py(),(1.-elossFrac)*part.pz(),(1.-elossFrac)*part.E()); //NB mass changes; could do more accurate...
                    }
		  }
                }

		//
		//   After energy loss
		//

		fastjet::ClusterSequenceArea clustSeq2(fjInputs, *jetDefAKT_Sig, fAreaDef);
		auto ElossJets = clustSeq2.inclusive_jets();

		// Fill track histo
		for (auto part : fjInputs) {
		  hTrackEtaPtLoss->Fill(part.pseudorapidity(),part.pt());
		}
                // Loop over jets
                for (auto jet : ElossJets) {
		  hJetEtaPtLoss->Fill(jet.pseudorapidity(), jet.pt());
                }

	}
//End event loop
   
    	outFile->Write();

	return 0;
}
