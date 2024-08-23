#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>

#include <random>
#include <chrono>

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
 * This macro is used to study a weighted jet finder
 * For every K0, we roll the dice whether to include it as input or not
 * In the limit of infinite events, we should find the same average jet properties
 */

const int pdgK0S = 310;
const int pdgLambda = 3122;

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

double z(fastjet::PseudoJet jet, fastjet::PseudoJet pj)
{
	double z = pj.px() * jet.px() + pj.py() * jet.py() + pj.pz() * jet.pz();
	z /= (jet.px() * jet.px() + jet.py() * jet.py() + jet.pz() * jet.pz());
	return z;
}

// Source: https://stackoverflow.com/questions/9878965/rand-between-0-and-1
double myRandom()
{
	std::mt19937_64 rng;
	// initialize the random number generator with time-dependent seed
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
	rng.seed(ss);
	// initialize a uniform distribution between 0 and 1
	std::uniform_real_distribution<double> unif(0, 1);
	// ready to generate random numbers
	return unif(rng);
}

bool randsig(double threshold)
{
	return (myRandom() < threshold);
}

bool isV0(int id, int setting = 1)
{
	switch (setting) {
		case 1:
			return (abs(id) == pdgK0S);
		case 2:
			return (abs(id) == pdgLambda);
		case 3:
			return (abs(id) == pdgK0S || abs(id) == pdgLambda);
		default:
			return false;
	}
}

int main(int argc, char** argv)
{
	int nEvents = 200;
  string outName = "weightedjetfinder";
	bool doPtScheme = false;
	double wSignal = 0.8;
	int v0setting = 1; // 1 = K0S, 2 = Lambda0, 3 = K0S + Lambda0

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
	if (argc >= 5) {
		wSignal = atof(argv[4]);
		if (wSignal < 0 || wSignal > 1) {
			cerr << "ERROR: Invalid signal fraction. Aborting.";
			return 1;
		}
	}
	if (argc >= 6) {
		v0setting = atoi(argv[5]);
		if (v0setting < 1 || v0setting > 3) {
			cerr << "ERROR: Invalid V0 setting. Aborting.";
			return 1;
		}
	}

	cout << "Running " << nEvents << " events" << endl;
	if (doPtScheme) {
		cout << "Using pt scheme for jet clustering" << endl;
	} else {
		cout << "Using E scheme for jet clustering" << endl;
	}

	cout << "Including a fraction " << wSignal << " of V0s" << endl;
	string sV0Setting = (1 == v0setting) ? "K0S" : ((2 == v0setting) ? "Lambda0" : "K0S and Lambda0");
	cout << "V0 setting = " << v0setting << " means we consider " << sV0Setting << endl;

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
	double jetptcut = 2.5;

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

	TH1D* hNEvts = new TH1D("hNEvts", "hNEvts", 3, -0.5, 2.5);
	hNEvts->Sumw2();

	TH1D* hNV0s = new TH1D("hNV0s", "hNV0s", 20, -0.5, 19.5);
	hNV0s->Sumw2();

	TH1D* hNW0s = new TH1D("hNW0s", "hNW0s", 20, -0.5, 19.5);
	hNW0s->Sumw2();

	// Jets
	TH3D* hV0Jet = new TH3D("hV0Jet", TString::Format("hV0Jet; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data(),
													nBins_pt_jet, min_pt_jet, max_pt_jet, nBins_eta_jet, -1.*max_eta_jet, max_eta_jet, nBins_phi_jet, min_phi_jet, max_phi_jet);
	hV0Jet->Sumw2();
	TH3D* hV0JetMatched = (TH3D*) hV0Jet->Clone("hV0JetMatched");

	TH3D* hW0Jet = (TH3D*) hV0Jet->Clone("hW0Jet");
	TH3D* hW0JetMatched = (TH3D*) hW0Jet->Clone("hW0JetMatched");

	TH3D* hChJet = (TH3D*) hV0Jet->Clone("hChJet");
	hChJet->SetTitle(TString::Format("hChJet; %s; %s; %s", sPtChjet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data());

	TH3D* hEJet = (TH3D*) hV0Jet->Clone("hEJet"); // Corrected with E scheme
	hEJet->SetTitle(TString::Format("hEJet; %s; %s; %s", sPtCorjet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data());
	TH3D* hEJetMatched = (TH3D*) hEJet->Clone("hEJetMatched");

	TH3D* hPtJet = (TH3D*) hEJet->Clone("hPtJet"); // Corrected with pt scheme
	TH3D* hPtJetMatched = (TH3D*) hPtJet->Clone("hPtJetMatched");

	// Difference hists
	THnSparseD* hJetDiff = new THnSparseD("hJetDiff", TString::Format("hJetDiff; %s - %s; %s; %s", sPtChjet.c_str(), sPtCorjet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data(),
																				4, new int[4]{nBins_pt_jet, nBins_pt_jet+1, nBins_eta_jet, nBins_phi_jet},
																					 new double[4]{min_pt_jet, -100.5, -1.*max_eta_jet, min_phi_jet},
																					 new double[4]{max_pt_jet, 100.5, max_eta_jet, max_phi_jet});
	hJetDiff->Sumw2();

	// THnSparseD* hEJetDiff = (THnSparseD*) hJetDiff->Clone("hEJetDiff");
	THnSparseD* hEJetDiffMatched = (THnSparseD*) hJetDiff->Clone("hEJetDiffMatched");
	// THnSparseD* hPtJetDiff = (THnSparseD*) hJetDiff->Clone("hPtJetDiff");
	THnSparseD* hPtJetDiffMatched = (THnSparseD*) hJetDiff->Clone("hPtJetDiffMatched");

	// V0 in jets
	TH3D* hV0inJet = new TH3D("hV0inJet", TString::Format("hV0inJet; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), sv0pt.c_str()).Data(),
														nBins_pt_jet, min_pt_jet, max_pt_jet, nBins_eta_jet, -1.*max_eta_jet, max_eta_jet, nBins_pt_track, min_pt_track, max_pt_track);
	hV0inJet->Sumw2();

	TH3D* hW0inJet = (TH3D*) hV0inJet->Clone("hW0inJet");

	// V0s
	TH3D* hV0 = new TH3D("hV0", TString::Format("hV0; %s; %s; %s", sv0pt.c_str(), sv0eta.c_str(), sv0phi.c_str()).Data(),
											 nBins_pt_track, min_pt_track, max_pt_track, nBins_eta_track, -1.*max_eta_track, max_eta_track, nBins_phi_track, min_phi_jet, max_phi_jet);
	hV0->Sumw2();

	TH3D* hW0 = (TH3D*) hV0->Clone("hW0");
	hW0->SetTitle("hV0 weighted; #it{p}_{T, W0}; #eta; #phi");

	TH3D* hTrack = (TH3D*) hV0->Clone("hTrack");
	hTrack->SetTitle("hTrack; #it{p}_{T, track}; #eta; #phi");

	// Z hists
	TH3D* hV0JetZ = new TH3D("hV0JetZ", TString::Format("hV0JetZ; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), szv0.c_str()).Data(),
													nBins_pt_jet, min_pt_jet, max_pt_jet, nBins_eta_jet, -1.*max_eta_jet, max_eta_jet, nBins_z, min_z, max_z);
	hV0JetZ->Sumw2();
	TH3D* hW0JetZ = (TH3D*) hV0JetZ->Clone("hW0JetZ");
	TH3D* hEJetZ = (TH3D*) hV0JetZ->Clone("hEJetZ");
	TH3D* hPtJetZ = (TH3D*) hV0JetZ->Clone("hPtJetZ");

  for (int iEvent = 0; iEvent < nEvents; iEvent++)
	{
		if (!pythia.next()) continue;
		hNEvts->Fill(0);
		double fourvec[4];
		int nPart = pythia.event.size();

		std::vector<fastjet::PseudoJet> chParticles;
		std::vector<fastjet::PseudoJet> v0Particles;
		std::vector<fastjet::PseudoJet> w0Particles;
		std::vector<int> v0Indices;
		for (int iPart = 0; iPart < nPart; iPart++) {
			const Particle &part = pythia.event[iPart];
			if (!part.isFinal()) continue;
			if (abs( part.eta() ) > max_eta_track || part.pT() < trackptcut) continue;

			fastjet::PseudoJet jInp(part.px(), part.py(), part.pz(), part.e());
			jInp.set_user_index(iPart);

			// Save K0S and (anti)Lambda0
			if (isV0(part.id(), v0setting)) {
				v0Indices.push_back(iPart);
				v0Particles.push_back(jInp);
				hV0->Fill(jInp.pt(), jInp.eta(), jInp.phi());

				if (randsig(wSignal)) {
					w0Particles.push_back(jInp);
					hW0->Fill(jInp.pt(), jInp.eta(), jInp.phi());
				}
			}
			else if (part.isCharged()) {
				chParticles.push_back(jInp);
				hTrack->Fill(jInp.pt(), jInp.eta(), jInp.phi());
			}
		} // Particle loop

		// Skip events, without V0s
		int nV0s = v0Indices.size();
		int nW0s = w0Particles.size();
		if (nV0s == 0) { continue; }
		hNEvts->Fill(1);
		hNV0s->Fill(nV0s);
		if (nW0s > 0) hNEvts->Fill(2);
		hNW0s->Fill(nW0s);

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

		fastjet::ClusterSequenceArea chCS(chParticles, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveChJets = chCS.inclusive_jets();
		std::vector <fastjet::PseudoJet> chJets = sorted_by_pt(inclusiveChJets); // Sort jets from high to low pt

		std::vector <fastjet::PseudoJet> v0Inputs = chParticles;
		v0Inputs.insert(v0Inputs.end(), v0Particles.begin(), v0Particles.end());
		fastjet::ClusterSequenceArea v0CS(v0Inputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveV0Jets = v0CS.inclusive_jets();
		std::vector <fastjet::PseudoJet> v0Jets = sorted_by_pt(inclusiveV0Jets); // Sort jets from high to low pt

		std::vector <fastjet::PseudoJet> w0Inputs = chParticles;
		w0Inputs.insert(w0Inputs.end(), w0Particles.begin(), w0Particles.end());
		fastjet::ClusterSequenceArea w0CS(w0Inputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveW0Jets = w0CS.inclusive_jets();
		std::vector <fastjet::PseudoJet> w0Jets = sorted_by_pt(inclusiveW0Jets); // Sort jets from high to low pt

		for (auto& chjet : chJets) {
			if (chjet.pt() < jetptcut) { continue; }
			hChJet->Fill(chjet.pt(), chjet.eta(), chjet.phi());
		} // chjet loop

		for (auto& v0jet : v0Jets) {
			if (v0jet.pt() < jetptcut) { continue; }
			hV0Jet->Fill(v0jet.pt(), v0jet.eta(), v0jet.phi());

			// Check correction as well
			fastjet::PseudoJet jetSubE(v0jet);
			double ptSub = v0jet.pt();

			bool jetContainsV0s = false;
			std::vector<double> zValues;
			for (const fastjet::PseudoJet &constituent : v0jet.constituents()) {
				int idx = constituent.user_index();
				const Particle &part = pythia.event[idx];
				if (isV0(part.id(), v0setting)) {
					jetContainsV0s = true;
					hV0inJet->Fill(v0jet.pt(), v0jet.eta(), constituent.pt());
					hV0JetZ->Fill(v0jet.pt(), v0jet.eta(), z(v0jet, constituent));

					// Subtract randomly
					if (randsig(wSignal)) {
						zValues.push_back(z(v0jet, constituent));
					} else {
						jetSubE -= constituent;
						ptSub -= constituent.pt();
					}
				}
			}
			hEJet->Fill(jetSubE.pt(), jetSubE.eta(), jetSubE.phi());
			hPtJet->Fill(ptSub, v0jet.eta(), v0jet.phi());
			if (!jetContainsV0s) { continue; }
			hV0JetMatched->Fill(v0jet.pt(), v0jet.eta(), v0jet.phi());
			hEJetMatched->Fill(jetSubE.pt(), jetSubE.eta(), jetSubE.phi());
			hPtJetMatched->Fill(ptSub, v0jet.eta(), v0jet.phi());
			hEJetDiffMatched->Fill(v0jet.pt(), v0jet.pt() - jetSubE.pt(), v0jet.eta() - jetSubE.eta(), v0jet.phi() - jetSubE.phi());
			hPtJetDiffMatched->Fill(v0jet.pt(), v0jet.pt() - ptSub, v0jet.eta() - v0jet.eta(), v0jet.phi() - v0jet.phi());

			// We don't know which V0s were randomly selected as signal, but we've saved their z values
			for (auto& z : zValues) {
				hEJetZ->Fill(jetSubE.pt(), jetSubE.eta(), z * (v0jet.pt() / jetSubE.pt()));
				hPtJetZ->Fill(ptSub, v0jet.eta(), z * (v0jet.pt() / ptSub));
			}
		}

		for (auto& w0jet : w0Jets) {
			if (w0jet.pt() < jetptcut) { continue; }
			hW0Jet->Fill(w0jet.pt(), w0jet.eta(), w0jet.phi());
			bool jetContainsV0s = false;
			for (const fastjet::PseudoJet &constituent : w0jet.constituents()) {
				int idx = constituent.user_index();
				const Particle &part = pythia.event[idx];
				if (isV0(part.id(), v0setting)) {
					jetContainsV0s = true;
					hW0inJet->Fill(w0jet.pt(), w0jet.eta(), constituent.pt());
					hW0JetZ->Fill(w0jet.pt(), w0jet.eta(), z(w0jet, constituent));
				}
			}
			if (!jetContainsV0s) { continue; }
			hW0JetMatched->Fill(w0jet.pt(), w0jet.eta(), w0jet.phi());
		}
	} // Event loop

	hNEvts->Write(hNEvts->GetName(), TObject::kOverwrite);

	hChJet->Write(hChJet->GetName(), TObject::kOverwrite);
	hV0Jet->Write(hV0Jet->GetName(), TObject::kOverwrite);
	hW0Jet->Write(hW0Jet->GetName(), TObject::kOverwrite);
	hEJet->Write(hEJet->GetName(), TObject::kOverwrite);
	hPtJet->Write(hPtJet->GetName(), TObject::kOverwrite);

	hV0JetMatched->Write(hV0JetMatched->GetName(), TObject::kOverwrite);
	hW0JetMatched->Write(hW0JetMatched->GetName(), TObject::kOverwrite);

	hEJetDiffMatched->Write(hEJetDiffMatched->GetName(), TObject::kOverwrite);
	hPtJetDiffMatched->Write(hPtJetDiffMatched->GetName(), TObject::kOverwrite);

	hV0inJet->Write(hV0inJet->GetName(), TObject::kOverwrite);
	hW0inJet->Write(hW0inJet->GetName(), TObject::kOverwrite);
	hV0JetZ->Write(hV0JetZ->GetName(), TObject::kOverwrite);
	hW0JetZ->Write(hW0JetZ->GetName(), TObject::kOverwrite);
	hEJetZ->Write(hEJetZ->GetName(), TObject::kOverwrite);
	hPtJetZ->Write(hPtJetZ->GetName(), TObject::kOverwrite);

	hTrack->Write(hTrack->GetName(), TObject::kOverwrite);
	hV0->Write(hV0->GetName(), TObject::kOverwrite);
	hW0->Write(hW0->GetName(), TObject::kOverwrite);

	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();

  return 0;
}
