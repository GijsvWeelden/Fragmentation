#include <iostream>
#include <string>
#include <vector>

#include "Pythia8/Pythia.h"

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

#include "TPDGCode.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

#define nEvents 20000

using namespace Pythia8;

std::vector <fastjet::PseudoJet> do_jet_finding();

int main(int /*argc*/, char** /*argv*/)
{
//PYTHIA SETTINGS

	TString name;

	int mecorr=1;

	Float_t ptHatMin=-1; //80;
	Float_t ptHatMax=-1; //100;

	int nCharged = 3, nNeutral = 5;
	int PDG[nCharged + nNeutral] = {211, 321, 2212, 111, 130, 310, 311, 3122};
	string Hadrons[nCharged + nNeutral] = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};

	// Generator. Process selection. LHC initialization. Histogram.
	Pythia pythia;
	pythia.readString("Beams:idA = 2212"); //beam 1 proton
	pythia.readString("Beams:idB = 2212"); //beam 2 proton
	pythia.readString("Beams:eCM = 13000.");
	pythia.readString("Tune:pp = 5");  //tune 1-13    5=defaulr TUNE4C,  6=Tune 4Cx, 7=ATLAS MB Tune A2-CTEQ6L1

	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = 0");

	pythia.readString("HardQCD:all = on");
	if(ptHatMin>0 && ptHatMax >0){
  	name = Form("PhaseSpace:pTHatMin = %f", (Float_t) ptHatMin);
		pythia.readString(name.Data());
		name = Form("PhaseSpace:pTHatMax = %f", (Float_t) ptHatMax);
		pythia.readString(name.Data());
	}

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

	// Output histograms
	TFile* outFile = new TFile("PythiaResult.root","RECREATE");
	TH2F *hEtaPt = new TH2F("hEtaPt","Pt vs Eta for all particles;#eta;p_{T} (GeV/c)",40,-2,2,50,0,10);
	TH1F* hists[nCharged + nNeutral];
	for (int i = 0; i < nCharged + nNeutral; i++){
		TH1F *hPt = new TH2F(TString::Format("hPt_%s", Hadrons[i].c_str()).Data(),
												 TString::Format("%s Pt;p_{T} (GeV/c)", Hadrons[i]).Data(),
												 50,0,10);
		hists[i] = *hPt;
	}

	float pt_lead = -1;
	float phi_lead = -100;
	float max_eta_jet = 2.0;
	float max_eta_track = 2.4;
	float jetR = 0.4;

	//Begin event loop
	for (int iEvent = 0; iEvent < nEvents; iEvent++){
		double fourvec[4];
		if (!pythia.next()) continue;

		Double_t ptSumPythia = 0;
		Int_t nPartPythia = 0;
		int nPart = pythia.event.size();
		for (int iPart = 0; iPart < nPart; iPart++){
    	const Particle &part = pythia.event[iPart];
			if (!part.isFinal()) continue;
			hEtaPt->Fill(part.eta(),part.pT());
			for (int i = 0; i < nCharged + nNeutral; i++){
				if (part.pdg_id() == PDG[i]){
					hists[i]->Fill(part.pT());
				}
			}
			// if (part.eta() < max_eta_track && part.pT() > min_track_pt){
			// 	fastjet::PseudoJet jInp(p->momentum().x(),p->momentum().y(),p->momentum().z(),p->momentum().e());  // need masses for E-scheme
			// 	jInp.set_user_index(part.pdg_id());//index);
			// 	fjInputs.push_back(jInp);
			// }
			// */
			nPartPythia++;
		}
		if ((iEvent%1000)==0)
				cout << "Pythia event: " << nPartPythia << " particles" << endl;
	}
	//End event loop
	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();
}

std::vector <fastjet::PseudoJet> do_jet_finding(){ // Takes in pythia event
	std::vector <fastjet::PseudoJet> particlesSig;
	int nPart = pythia.event.size();
	for (int iPart = 0; iPart < nPart; iPart++){
		const Particle &part = pythia.event[iPart];
		if (!part.isFinal())
		if (part.eta() > max_eta_track)
			continue; // Do we want this? Lambdas
			hEtaPt->Fill(part.eta(),part.pT());
			for (int i = 0; i < nCharged + nNeutral; i++){
				if (part.pdg_id() == PDG[i]){
					hists[i]->Fill(part.pT());
				}
			}
			// if (part.eta() < max_eta_track && part.pT() > min_track_pt){
			// 	fastjet::PseudoJet jInp(p->momentum().x(),p->momentum().y(),p->momentum().z(),p->momentum().e());  // need masses for E-scheme
			// 	jInp.set_user_index(part.pdg_id());//index);
			// 	fjInputs.push_back(jInp);
			// }
			// */
			nPartPythia++;
	}


	const HepMC::GenParticle *p = *pit; // for each particle
		if ( fabs(p->momentum().eta()) < max_eta_track && p->momentum().perp() > ptcut){
			fastjet::PseudoJet jInp(p->momentum().x(),p->momentum().y(),p->momentum().z(),p->momentum().e());  // need masses for E-scheme
			jInp.set_user_index(p->pdg_id());//index);
			fjInputs.push_back(jInp);
			index++;
		}
	fastjet::ClusterSequenceArea csSig(particlesSig, jetDef, areaDef);
	jetCollection jetCollection_Sig( sorted_by_pt( jet_selector( csSig.inclusive_jets(130.) ) ) );
		fastjet::GhostedAreaSpec ghostSpec(max_eta_track,1,0.01);
    fastjet::Strategy strategy = fastjet::Best;
    //fastjet::RecombinationScheme    recombScheme = fastjet::BIpt_scheme;
    fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; // need E scheme for jet mass
    fastjet::AreaType areaType = fastjet::active_area;
    fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
    fastjet::AreaDefinition areaDefShape = fastjet::AreaDefinition(fastjet::active_area,ghostSpec);

    fastjet::RangeDefinition range(-max_eta_jet, max_eta_jet, 0, 2.*fastjet::pi);

    fastjet::JetDefinition jetDefCh(fastjet::antikt_algorithm, jetR, recombScheme, strategy);
    fastjet::ClusterSequenceArea clustSeqCh(fjInputs, jetDefCh, areaDef);

    vector <fastjet::PseudoJet> inclusiveJetsCh = clustSeqCh.inclusive_jets();
    if (inclusiveJetsCh.size() <= 0){
      // delete the created event from memory
      delete evt;
      // read the next event
      ascii_in >> evt;
      ievt++;
      continue; // Skip events without jets
    }
}
