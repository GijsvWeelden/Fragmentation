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

// std::vector <fastjet::PseudoJet> do_jet_finding();
int find_matriarch(const Event& event, const Particle& particle);

int main(int /*argc*/, char** /*argv*/)
{
//PYTHIA SETTINGS

	TString name;

	int mecorr=1;

	Float_t ptHatMin=-1; //80;
	Float_t ptHatMax=-1; //100;

	float max_eta_track = 2, min_pt_track = 0., max_pt_track = 10.;
	float max_eta_jet = 2.0, min_pt_jet = 10, max_pt_jet = 200, jetR = 0.4;
	int nBins_eta_track = 40, nBins_pt_track = 50, nBins_eta_jet = 40, nBins_pt_jet = 200;

	int nCharged = 3, nNeutral = 5;
	std::vector<int> PDG = {211, 321, 2212, 111, 130, 310, 311, 3122};
	std::vector<string> Hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};

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

	pythia.readString("111:mayDecay  = off"); //pi0
	pythia.readString("310:mayDecay  = off"); //K0s
	pythia.readString("311:mayDecay  = off"); //K0
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
	TH2F *hEtaPt = new TH2F("hEtaPt","Pt vs Eta for all particles;#eta;p_{T} (GeV/c)", 40, -2, 2, 50, 0, 10);
	TH1F* hists[nCharged + nNeutral];
	for (int i = 0; i < nCharged + nNeutral; i++){
		TH1F *hPt = new TH1F(TString::Format("hPt_%s", Hadrons[i].c_str()).Data(),
												 TString::Format("%s Pt;p_{T} (GeV/c)", Hadrons[i].c_str()).Data(),
												 50,0,10);
		hists[i] = hPt;
	}

	//Begin event loop
	for (int iEvent = 0; iEvent < nEvents; iEvent++){
		if (!pythia.next()) continue;
		double fourvec[4];
    std::vector<Particle> matriarchs; // Outgoing particles from initial hard scattering

		Double_t ptSumPythia = 0;
		Int_t nMatriarchs = 0;
		Int_t matriarchID = -1;
		Int_t nPartPythia = 0;
		int nPart = pythia.event.size();

		for (int iPart = 0; iPart < nPart; iPart++){
    	const Particle &part = pythia.event[iPart];
			if (part.status() == -23){ // TODO: Should also include 22 and 24?
				nMatriarchs++;
        if (iEvent%1000 == 0) cout << "Found Initial particle: " << part.index() << ": " << part.id() << endl;
				if (nMatriarchs > 2){
					cout << "Warning: More than 2 outgoing particles found from the initial hard scattering. We will ignore these." << endl;
					continue;
				}
        matriarchs.push_back(part);
				continue;
			}
			if (!part.isFinal()) continue; // No decays yet
			// if (part.eta() > max_eta_track || part.pT() < min_track_pt) continue;
			hEtaPt->Fill(part.eta(),part.pT());
			for (int i = 0; i < nCharged + nNeutral; i++){
				if (part.id() == PDG[i]){
					hists[i]->Fill(part.pT());
				}
			}
			matriarchID = find_matriarch(pythia.event, part);
			if (iEvent%1000 == 0 && nPartPythia%10 == 0) cout << "Matriarch found: " << matriarchID << endl;
			nPartPythia++;
		}
		if ((iEvent%1000)==0){
			cout << "Pythia event: " << nPartPythia << " particles" << endl;
			//cout << nMatriarchs << " outgoing partons from hard scattering:" << endl
			//	<< matriarchs[0].id() << " " << matriarchs[0].p() << endl
			//	<< matriarchs[0].id() << " " << matriarchs[1].p() << endl;
		}
	}
	//End event loop
	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();
}

int find_matriarch(const Event& event, const Particle& particle){
	Particle mother1, mother2;
	Particle part = particle;
	for (int i = 0; i < event.size(); i++){
	// while (true){ // Could this loop infinitely?
		mother1 = event[part.mother1()];
		mother2 = event[part.mother2()];
		if (abs(mother1.status()) == 23){
			return mother1.id();
		}
		else if (abs(mother2.status()) == 23){
			return mother2.id();
		}
		else if (mother1.id() == 1 || mother1.id() == 2 || mother2.id() == 1 || mother2.id() == 2){
			// Particle originates from beam
			return -1;
		}
		part = mother1; // Had to choose one to avoid branching. Is this smart?
	}
	return -2;
}

/*
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
*/
