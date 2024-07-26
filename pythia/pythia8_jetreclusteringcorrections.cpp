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
 * This macro is used to study the effect of treating Lambdas as K0S during jet clustering
 *
 */

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

double z(fastjet::PseudoJet jet, fastjet::PseudoJet pj)
{
	double z = pj.px() * jet.px() + pj.py() * jet.py() + pj.pz() * jet.pz();
	z /= (jet.px() * jet.px() + jet.py() * jet.py() + jet.pz() * jet.pz());
	return z;
}

int main(int argc, char** argv)
{
	int nEvents = 200;
  string outName = "jetreclusteringcorrections";
	bool doPtScheme = false;

	if (argc >= 2){
		nEvents = atoi(argv[1]);
		if (nEvents <= 0){
			cerr << "ERROR: Zero events requested. Aborting.";
			return 1;
		}
	}
	if (argc >= 3) {
		outName = argv[2];
	}
	if (argc >= 4) {
		doPtScheme = atoi(argv[3]);
	}

	cout << "Running " << nEvents << " events" << endl;
	if (doPtScheme) {
		cout << "Using pt scheme for jet clustering" << endl;
	} else {
		cout << "Using E scheme for jet clustering" << endl;
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

	string sPtV0jet  = "#it{p}_{T, Ch+V0 jet} [GeV/#it{c}]";
	string sPtCorjet = "#it{p}_{T, Ch+V0 jet - V0} [GeV/#it{c}]";
	string sPtChjet  = "#it{p}_{T, Ch jet} [GeV/#it{c}]";

	string sEtaJet  = "#eta_{jet}";
	string sPhiJet  = "#phi_{jet}";

	string sv0pt  = "#it{p}_{T, V0}";
	string sv0eta = "#eta_{V0}";
	string sv0phi = "#phi_{V0}";
	string szv0   = "#it{z}_{V0}";
	string szk0   = "#it{z}_{K0}";

	// Output histograms
	if (outName.find(".root") == string::npos) {
		outName += ".root";
	}
	TFile* outFile = new TFile(TString::Format("%s", outName.c_str()).Data(), "RECREATE");

	TH1D* hNEvts = new TH1D("hNEvts", "hNEvts", 2, -0.5, 1.5);
	hNEvts->Sumw2();

	// Jets
	TH3D* hV0Jet = new TH3D("hV0Jet", TString::Format("hV0Jet; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data(),
													nBins_pt_jet, min_pt_jet, max_pt_jet, nBins_eta_jet, -1.*max_eta_jet, max_eta_jet, nBins_phi_jet, min_phi_jet, max_phi_jet);
	hV0Jet->Sumw2();
	TH3D* hV0JetMatched = (TH3D*) hV0Jet->Clone("hV0JetMatched");

	TH3D* hChJet = (TH3D*) hV0Jet->Clone("hChJet");
	hChJet->SetTitle(TString::Format("hChJet; %s; %s; %s", sPtChjet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data());
	TH3D* hChJetMatched = (TH3D*) hChJet->Clone("hChJetMatched");

	TH3D* hEJet = (TH3D*) hV0Jet->Clone("hEJet"); // Corrected with E scheme
	hEJet->SetTitle(TString::Format("hEJet; %s; %s; %s", sPtCorjet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data());
	TH3D* hEJetMatched = (TH3D*) hEJet->Clone("hEJetMatched");

	TH3D* hPtJet = (TH3D*) hEJet->Clone("hPtJet"); // Corrected with pt scheme
	TH3D* hPtJetMatched = (TH3D*) hPtJet->Clone("hPtJetMatched");

	// Difference hists
	THnSparseD* hJetDiff = new THnSparseD("hJetDiff", TString::Format("hJetDiff; %s; %s - %s; %s; %s", sPtV0jet.c_str(), sPtV0jet.c_str(), sPtChjet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data(),
																				4, new int[4]{nBins_pt_jet, nBins_pt_jet+1, nBins_eta_jet, nBins_phi_jet},
																					 new double[4]{min_pt_jet, -100.5, -1.*max_eta_jet, min_phi_jet},
																					 new double[4]{max_pt_jet, 100.5, max_eta_jet, max_phi_jet});
	hJetDiff->Sumw2();

	THnSparseD* hEJetDiff = (THnSparseD*) hJetDiff->Clone("hEJetDiff");
	hEJetDiff->SetTitle(TString::Format("hEJetDiff; %s; %s - %s; %s; %s", sPtV0jet.c_str(), sPtV0jet.c_str(), sPtCorjet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data());
	THnSparseD* hEJetDiffMatched = (THnSparseD*) hJetDiff->Clone("hEJetDiffMatched");
	THnSparseD* hPtJetDiff = (THnSparseD*) hJetDiff->Clone("hPtJetDiff");
	THnSparseD* hPtJetDiffMatched = (THnSparseD*) hJetDiff->Clone("hPtJetDiffMatched");

	// V0s
	TH3D* hV0 = new TH3D("hV0", TString::Format("hV0; %s; %s; %s", sv0pt.c_str(), sv0eta.c_str(), sv0phi.c_str()).Data(),
											 nBins_pt_track, min_pt_track, max_pt_track, nBins_eta_track, -1.*max_eta_track, max_eta_track, nBins_phi_track, min_phi_jet, max_phi_jet);
	hV0->Sumw2();

	TH3D* hTrack = (TH3D*) hV0->Clone("hTrack");
	hTrack->SetTitle("hTrack; #it{p}_{T, track}; #eta; #phi");

  for (int iEvent = 0; iEvent < nEvents; iEvent++)
	{
		if (!pythia.next()) continue;
		hNEvts->Fill(0);
		double fourvec[4];
		int nPart = pythia.event.size();

		std::vector<fastjet::PseudoJet> chParticles;
		std::vector<fastjet::PseudoJet> v0Particles;
		std::vector<int> v0Indices;
		for (int iPart = 0; iPart < nPart; iPart++) {
			const Particle &part = pythia.event[iPart];
			if (!part.isFinal()) continue;
			if (abs( part.eta() ) > max_eta_track || part.pT() < trackptcut) continue;

			fastjet::PseudoJet jInp(part.px(), part.py(), part.pz(), part.e());
			jInp.set_user_index(iPart);

			// Save K0S and (anti)Lambda0
			if ( abs(part.id()) == 310 || abs(part.id()) == 3122 ){
				v0Indices.push_back(iPart);
				v0Particles.push_back(jInp);
				hV0->Fill(jInp.pt(), jInp.eta(), jInp.phi());
			}
			else if (part.isCharged()) {
				chParticles.push_back(jInp);
				hTrack->Fill(jInp.pt(), jInp.eta(), jInp.phi());
			}
		} // Particle loop

		// Skip events, without V0s
		int nV0s = v0Indices.size();
		if (nV0s == 0) { continue; }
		hNEvts->Fill(1);

		// Do jet finding and analysis here.
		fastjet::GhostedAreaSpec ghostSpec(max_eta_track, 1, 0.01);
		fastjet::Strategy strategy = fastjet::Best;
		fastjet::AreaType areaType = fastjet::active_area;
		fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::AreaDefinition areaDefShape = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::RangeDefinition range(-1. * max_eta_jet, max_eta_jet, 0, 2.*fastjet::pi);
		fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; // need E scheme for jet mass
		if (doPtScheme) { recombScheme = fastjet::pt_scheme; }
		fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR, recombScheme, strategy);

		std::vector <fastjet::PseudoJet> v0Inputs = chParticles;
		v0Inputs.insert(v0Inputs.end(), v0Particles.begin(), v0Particles.end());
		fastjet::ClusterSequenceArea v0CS(v0Inputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveV0Jets = v0CS.inclusive_jets();
		std::vector <fastjet::PseudoJet> v0Jets = sorted_by_pt(inclusiveV0Jets); // Sort jets from high to low pt

		for (auto& v0jet : v0Jets) {
			if (v0jet.pt() < jetptcut) { continue; }
			bool jetContainsV0s = false;
			for (const fastjet::PseudoJet &constituent : v0jet.constituents()) {
				int idx = constituent.user_index();
				const Particle &part = pythia.event[idx];
				if (abs(part.id()) == 310 || abs(part.id()) == 3122) {
					jetContainsV0s = true;
				}
			}
			if (!jetContainsV0s) { continue; }
			hV0Jet->Fill(v0jet.pt(), v0jet.eta(), v0jet.phi());

			fastjet::PseudoJet jetSubE(v0jet);
			double ptSub = v0jet.pt();

			for (const fastjet::PseudoJet &constituent : v0jet.constituents()) {
				int idx = constituent.user_index();
				const Particle &part = pythia.event[idx];
				if (abs(part.id()) == 310 || abs(part.id()) == 3122) {
					fastjet::PseudoJet pj(part.px(), part.py(), part.pz(), part.e());
					jetSubE -= pj;
					ptSub -= part.pT();
				}
			}
			hEJet->Fill(jetSubE.pt(), jetSubE.eta(), jetSubE.phi());
			hPtJet->Fill(ptSub, v0jet.eta(), v0jet.phi());
			hEJetDiff->Fill(v0jet.pt(), v0jet.pt() - jetSubE.pt(), v0jet.eta() - jetSubE.eta(), v0jet.delta_phi_to(jetSubE));
			hPtJetDiff->Fill(v0jet.pt(), v0jet.pt() - ptSub, v0jet.eta() - v0jet.eta(), v0jet.delta_phi_to(v0jet));

			// Recluster v0jet with anti-kt without the V0s and compare to ch+V0 jet
			std::vector<fastjet::PseudoJet> chjetInputs;
			for (fastjet::PseudoJet &constituent : v0jet.constituents()) {
				int idx = constituent.user_index();
				const Particle &part = pythia.event[idx];
				if (abs(part.id()) == 310 || abs(part.id()) == 3122) {
					continue;
				}
				chjetInputs.push_back(constituent);
			}
			fastjet::ClusterSequenceArea chCS(chjetInputs, jetDef, areaDef);
			std::vector <fastjet::PseudoJet> inclusiveChJets = chCS.inclusive_jets();
			std::vector <fastjet::PseudoJet> chJets = sorted_by_pt(inclusiveChJets); // Sort jets from high to low pt

			auto& chjet = chJets[0];
			hChJet->Fill(chjet.pt(), chjet.eta(), chjet.phi());
			hJetDiff->Fill(v0jet.pt(), v0jet.pt() - chjet.pt(), v0jet.eta() - chjet.eta(), v0jet.delta_phi_to(chjet));
			hPtJetDiff->Fill(v0jet.pt(), ptSub - chjet.pt(), v0jet.eta() - chjet.eta(), v0jet.delta_phi_to(chjet));
			hEJetDiffMatched->Fill(v0jet.pt(), jetSubE.pt() - v0jet.pt(), jetSubE.eta() - v0jet.eta(), jetSubE.delta_phi_to(v0jet));
		} // v0jet loop
	} // Event loop

	hNEvts->Write(hNEvts->GetName(), TObject::kOverwrite);

	hV0Jet->Write(hV0Jet->GetName(), TObject::kOverwrite);
	hChJet->Write(hChJet->GetName(), TObject::kOverwrite);
	hEJet->Write(hEJet->GetName(), TObject::kOverwrite);
	hPtJet->Write(hPtJet->GetName(), TObject::kOverwrite);

	hJetDiff->Write(hJetDiff->GetName(), TObject::kOverwrite);
	hEJetDiff->Write(hEJetDiff->GetName(), TObject::kOverwrite);
	hEJetDiffMatched->Write(hEJetDiffMatched->GetName(), TObject::kOverwrite);
	hPtJetDiff->Write(hPtJetDiff->GetName(), TObject::kOverwrite);
	hPtJetDiffMatched->Write(hPtJetDiffMatched->GetName(), TObject::kOverwrite);

	hV0->Write(hV0->GetName(), TObject::kOverwrite);
	hTrack->Write(hTrack->GetName(), TObject::kOverwrite);

	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();

  return 0;
}
