#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

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
#include "TH1D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TH2F.h"

using namespace Pythia8;

/*
 * This macro is used to study the difference between ch. jets and ch.+V0 jets
 *
 * For ch jets, we consider ch jets that have a V0 in the jet cone
 * We compare the ch jets with ch jets that have had their momentum corrected for the V0
 * Lastly, we use ch+V0 jets
 */

double z(fastjet::PseudoJet jet, fastjet::PseudoJet pj)
{
	double z = pj.px() * jet.px() + pj.py() * jet.py() + pj.pz() * jet.pz();
	z /= (jet.px() * jet.px() + jet.py() * jet.py() + jet.pz() * jet.pz());
	return z;
}

int main(int argc, char** argv)
{
	int nEvents = 200;
  string outName = "chJetV0JetStudy";

	if (argc >= 2){
		nEvents = atoi(argv[1]);
		if (nEvents == 0){
			cerr << "ERROR: Zero events requested. Aborting.";
			return 1;
		}
	}
	if (argc >= 3) {
		outName = argv[2];
	}

	double ptHatMin = -20;
	double ptHatMax = -80;

	//PYTHIA SETTINGS
	TString name;
	int mecorr=1;

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
  	name = Form("PhaseSpace:pTHatMin = %f", (double) ptHatMin);
		pythia.readString(name.Data());
		name = Form("PhaseSpace:pTHatMax = %f", (double) ptHatMax);
		pythia.readString(name.Data());
	}

	pythia.readString("111:mayDecay  = off"); //pi0
	pythia.readString("310:mayDecay  = off"); //K0s
	// pythia.readString("311:mayDecay  = off"); //K0
	pythia.readString("3122:mayDecay = off"); //labda0

	//ME corrections
	//use of matrix corrections where available
	if(mecorr==0){
		pythia.readString("TimeShower:MECorrections=off");
	}
	pythia.init();

  // Eta cut to ignore events with partons too far away
  double max_eta_matriarch = 2.5;
	// Settings for tracks and jets
	double max_eta_track = 2;
	double jetR = 0.4;
	double trackptcut = 0.1;
	double jetptcut = 5.;

	// Hist binning
	int nBins_pt_jet = 200;
	double min_pt_jet = .0, max_pt_jet = 200;
	int nBins_eta_jet = 40;
	double max_eta_jet = max_eta_track - jetR;
	int nBins_phi_jet = 30;
	double min_phi_jet = -1.*fastjet::pi, max_phi_jet = 2.*fastjet::pi;

	int nBins_pt_track = 120;
	double min_pt_track = 0., max_pt_track = 60.;
	int nBins_eta_track = 40;
	int nBins_phi_track = 30;

	int nBins_z = 100;
	double min_z = 1e-3, max_z = 1.001;

	string sPtChJet = "#it{p}_{T, ch. jet} [GeV]";
	string sPtV0jet = "#it{p}_{T, ch.+V0 jet} [GeV]";
	string sEtaJet   = "#eta_{jet}";
	string sPhiJet   = "#phi_{jet}";

	string sv0pt  = "#it{p}_{T, V0}";
	string sv0eta = "#eta_{V0}";
	string sv0phi = "#phi_{V0}";
	string sz     = "#it{z}_{V0}";

	// Output histograms
	TFile* outFile = new TFile(TString::Format("%s.root", outName.c_str()).Data(), "RECREATE");

	TH1D* hNEvts = new TH1D("hNEvts", "hNEvts", 2, -0.5, 1.5);

	// Jets
	TH3D* hChJet = new TH3D("hChJet", TString::Format("hChJet; %s; %s; %s", sPtChJet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data(),
													nBins_pt_jet, min_pt_jet, max_pt_jet, nBins_eta_jet, -1.*max_eta_jet, max_eta_jet, nBins_phi_jet, min_phi_jet, max_phi_jet);

	TH3D* hChJetCorrected = (TH3D*) hChJet->Clone("hChJetCorrected");
	hChJetCorrected->SetTitle(TString::Format("hChJetCorrected; %s; %s; %s", sPtChJet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data());

	TH3D* hV0Jet = (TH3D*) hChJet->Clone("hV0Jet");
	hV0Jet->SetTitle(TString::Format("hV0Jet; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data());

	// V0s
	TH3D* hV0 = new TH3D("hV0", TString::Format("hV0; %s; %s; %s", sv0pt.c_str(), sv0eta.c_str(), sv0phi.c_str()).Data(),
											 nBins_pt_track, min_pt_track, max_pt_track, nBins_eta_track, -1.*max_eta_track, max_eta_track, nBins_phi_track, min_phi_jet, max_phi_jet);

	TH3D* hTrack = (TH3D*) hV0->Clone("hTrack");
	hTrack->SetTitle("hTrack; #it{p}_{T, track}; #eta; #phi");

	// V0s in jets
	THnSparseD* hzCh
		= new THnSparseD("hzCh", TString::Format("hzCh; %s; %s; %s", sPtChJet.c_str(), sEtaJet.c_str(), sPhiJet.c_str(), sz.c_str()).Data(),
										 4, new int[4]{nBins_pt_jet, nBins_eta_jet, nBins_phi_jet, nBins_z},
												new double[4]{min_pt_jet, -1*max_eta_jet, min_phi_jet, min_z},
												new double[4]{max_pt_jet, max_eta_jet, max_phi_jet, max_z});
	THnSparseD* hzChCorrected = (THnSparseD*) hzCh->Clone("hzChCorrected");
	hzChCorrected->SetTitle(TString::Format("hzChCorrected; %s; %s; %s", sPtChJet.c_str(), sEtaJet.c_str(), sPhiJet.c_str(), sz.c_str()).Data());
	THnSparseD* hzV0 = (THnSparseD*) hzCh->Clone("hzV0");
	hzV0->SetTitle(TString::Format("hzV0; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str(), sz.c_str()).Data());

  for (int iEvent = 0; iEvent < nEvents; iEvent++)
	{
		if (!pythia.next()) continue;
		hNEvts->Fill(0);
		double fourvec[4];

		int nPart = pythia.event.size();

		std::vector<fastjet::PseudoJet> chjetInputs;
		std::vector<fastjet::PseudoJet> v0jetInputs;
		std::vector<int> v0Indices;

		for (int iPart = 0; iPart < nPart; iPart++)
    {
    	const Particle &part = pythia.event[iPart];
			if (!part.isFinal()) continue;
			if (abs( part.eta() ) > max_eta_track || part.pT() < trackptcut) continue;

			fastjet::PseudoJet jInp(part.px(), part.py(), part.pz(), part.e());
			jInp.set_user_index(iPart);

			// Save K0S and (anti)Lambda0
			if ( (abs(part.id()) == 310) || (abs(part.id()) == 3122) ){
				v0Indices.push_back(iPart);
				v0jetInputs.push_back(jInp);
				hV0->Fill(jInp.pt(), jInp.eta(), jInp.phi());
			}
			else {
				if (!part.isCharged()) continue;
				chjetInputs.push_back(jInp);
				v0jetInputs.push_back(jInp);
			}
		} // Particle loop

		// Skip events without V0s
		int nV0s = v0Indices.size();
		if (nV0s == 0) { continue; }
		hNEvts->Fill(1);

		for (auto& part : chjetInputs) {
			hTrack->Fill(part.pt(), part.eta(), part.phi());
		}

		// Do jet finding and analysis here.
		fastjet::GhostedAreaSpec ghostSpec(max_eta_track, 1, 0.01);
		fastjet::Strategy strategy = fastjet::Best;
		fastjet::AreaType areaType = fastjet::active_area;
		fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::AreaDefinition areaDefShape = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::RangeDefinition range(-1. * max_eta_jet, max_eta_jet, 0, 2.*fastjet::pi);
		fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; // need E scheme for jet mass
		fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR, recombScheme, strategy);

		fastjet::ClusterSequenceArea clustSeqCh(chjetInputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveChJets = clustSeqCh.inclusive_jets();
		std::vector <fastjet::PseudoJet> chJets = sorted_by_pt(inclusiveChJets); // Sort jets from high to low pt

		fastjet::ClusterSequenceArea clustSeqV0(v0jetInputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveV0Jets = clustSeqV0.inclusive_jets();
		std::vector <fastjet::PseudoJet> v0Jets = sorted_by_pt(inclusiveV0Jets); // Sort jets from high to low pt

		bool isV0Used[nV0s];
    for (int i = 0; i < nV0s; i++) {
      isV0Used[i] = false;
    }

		for (auto& jet : chJets) {
      if (jet.pt() < jetptcut) { continue; }
			fastjet::PseudoJet newjet(jet);
			bool jetContainsV0s = false;

			for (int i = 0; i < nV0s; i++) {
				if (isV0Used[i]) { continue; }
				const Particle &v0 = pythia.event[v0Indices[i]];
				fastjet::PseudoJet pjV0(v0.px(), v0.py(), v0.pz(), v0.e());
				bool v0inJet = (pjV0.delta_R(jet) <= jetR);
				if (!v0inJet) { continue; }
				jetContainsV0s = true;
				newjet += pjV0;
			} // V0 loop: momentum correction

			if (!jetContainsV0s) { continue; }

			hChJet->Fill(jet.pt(), jet.eta(), jet.phi());
			hChJetCorrected->Fill(newjet.pt(), newjet.eta(), newjet.phi());

			for (int i = 0; i < nV0s; i++) {
				if (isV0Used[i]) { continue; }
				const Particle &v0 = pythia.event[v0Indices[i]];
				fastjet::PseudoJet pjV0(v0.px(), v0.py(), v0.pz(), v0.e());
				bool v0inJet = (pjV0.delta_R(jet) <= jetR);
				if (!v0inJet) { continue; }

				isV0Used[i] = true;
				hzCh->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, pjV0));
				hzChCorrected->Fill(newjet.pt(), newjet.eta(), newjet.phi(), z(newjet, pjV0));
			} // V0 loop: hist filling
		} // ch jet loop

		for (auto& jet : v0Jets) {
			if (jet.pt() < jetptcut) { continue; }
			bool jetContainsV0s = false;
			for (const fastjet::PseudoJet &constituent : jet.constituents()) {
				int idx = constituent.user_index();
				const Particle &part = pythia.event[idx];

				if ( (abs(part.id()) == 310) || (abs(part.id()) == 3122) ){
					jetContainsV0s = true;
					fastjet::PseudoJet pj(part.px(), part.py(), part.pz(), part.e());
					hzV0->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, pj));
				} // if V0
			} // constituent loop
			if (!jetContainsV0s) { continue; }
			hV0Jet->Fill(jet.pt(), jet.eta(), jet.phi());
		} // V0 jet loop
	}

	hNEvts->Write(hNEvts->GetName(), TObject::kOverwrite);
	hV0->Write(hV0->GetName(), TObject::kOverwrite);
	hTrack->Write(hTrack->GetName(), TObject::kOverwrite);

	hChJet->Write(hChJet->GetName(), TObject::kOverwrite);
	hChJetCorrected->Write(hChJetCorrected->GetName(), TObject::kOverwrite);
	hV0Jet->Write(hV0Jet->GetName(), TObject::kOverwrite);

	hzCh->Write(hzCh->GetName(), TObject::kOverwrite);
	hzChCorrected->Write(hzChCorrected->GetName(), TObject::kOverwrite);
	hzV0->Write(hzV0->GetName(), TObject::kOverwrite);

	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();

  return 0;
}
