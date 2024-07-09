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
  string outName = "studyLambdaAsK0S";
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
		cout << "Using pt recombination scheme" << endl;
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

	string sPtK0jet = "#it{p}_{T, Ch+K0 jet} [GeV/#it{c}]";
	string sPtV0jet = "#it{p}_{T, Ch+V0 jet} [GeV/#it{c}]";
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

	// Jets
	TH3D* hV0Jet = new TH3D("hV0Jet", TString::Format("hV0Jet; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data(),
													nBins_pt_jet, min_pt_jet, max_pt_jet, nBins_eta_jet, -1.*max_eta_jet, max_eta_jet, nBins_phi_jet, min_phi_jet, max_phi_jet);

	TH3D* hK0Jet = (TH3D*) hV0Jet->Clone("hK0Jet"); // Treating L as K
	hK0Jet->SetTitle(TString::Format("hK0Jet; %s; %s; %s", sPtK0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str()).Data());

	// V0s
	TH3D* hV0 = new TH3D("hV0", TString::Format("hV0; %s; %s; %s", sv0pt.c_str(), sv0eta.c_str(), sv0phi.c_str()).Data(),
											 nBins_pt_track, min_pt_track, max_pt_track, nBins_eta_track, -1.*max_eta_track, max_eta_track, nBins_phi_track, min_phi_jet, max_phi_jet);

	TH3D* hTrack = (TH3D*) hV0->Clone("hTrack");
	hTrack->SetTitle("hTrack; #it{p}_{T, track}; #eta; #phi");

	// V0s in jets
	THnSparseD* hzV0_K0S
		= new THnSparseD("hzV0_K0S", TString::Format("hzV0_K0S; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str(), szv0.c_str()).Data(),
										 4, new int[4]{nBins_pt_jet, nBins_eta_jet, nBins_phi_jet, nBins_z},
												new double[4]{min_pt_jet, -1*max_eta_jet, min_phi_jet, min_z},
												new double[4]{max_pt_jet, max_eta_jet, max_phi_jet, max_z});
	THnSparseD* hzV0_Lambda0 = (THnSparseD*) hzV0_K0S->Clone("hzV0_Lambda0");
	hzV0_Lambda0->SetTitle(TString::Format("hzV0_Lambda0; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str(), szv0.c_str()).Data());
	THnSparseD* hzV0_V0 = (THnSparseD*) hzV0_K0S->Clone("hzV0_V0");
	hzV0_V0->SetTitle(TString::Format("hzV0_V0; %s; %s; %s", sPtV0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str(), szv0.c_str()).Data());

	THnSparseD* hzK0_K0S = (THnSparseD*) hzV0_K0S->Clone("hzK0_K0S");
	hzK0_K0S->SetTitle(TString::Format("hzK0_K0S; %s; %s; %s", sPtK0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str(), szk0.c_str()).Data());
	THnSparseD* hzK0_Lambda0 = (THnSparseD*) hzV0_K0S->Clone("hzK0_Lambda0");
	hzK0_Lambda0->SetTitle(TString::Format("hzK0_Lambda0; %s; %s; %s", sPtK0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str(), szk0.c_str()).Data());
	THnSparseD* hzK0_V0 = (THnSparseD*) hzV0_K0S->Clone("hzK0_V0");
	hzK0_V0->SetTitle(TString::Format("hzK0_V0; %s; %s; %s", sPtK0jet.c_str(), sEtaJet.c_str(), sPhiJet.c_str(), szk0.c_str()).Data());

  for (int iEvent = 0; iEvent < nEvents; iEvent++)
	{
		if (!pythia.next()) continue;
		hNEvts->Fill(0);
		double fourvec[4];

		int nPart = pythia.event.size();

		std::vector<fastjet::PseudoJet> chParticles;
		std::vector<fastjet::PseudoJet> K0SParticles;
		std::vector<fastjet::PseudoJet> L0Particles;
		std::vector<fastjet::PseudoJet> L0AsK0SParticles;
		std::vector<int> v0Indices;

		for (int iPart = 0; iPart < nPart; iPart++)
    {
    	const Particle &part = pythia.event[iPart];
			if (!part.isFinal()) continue;
			if (abs( part.eta() ) > max_eta_track || part.pT() < trackptcut) continue;

			fastjet::PseudoJet jInp(part.px(), part.py(), part.pz(), part.e());
			jInp.set_user_index(iPart);

			// Save K0S and (anti)Lambda0
			if (abs(part.id()) == 310) {
				v0Indices.push_back(iPart);
				K0SParticles.push_back(jInp);
				hV0->Fill(jInp.pt(), jInp.eta(), jInp.phi());
			}
			else if (abs(part.id()) == 3122) {
				v0Indices.push_back(iPart);
				L0Particles.push_back(jInp);
				hV0->Fill(jInp.pt(), jInp.eta(), jInp.phi());

				double wrongEnergy = sqrt( part.pAbs2() + MassK0S*MassK0S );
				fastjet::PseudoJet jInpK0S(part.px(), part.py(), part.pz(), wrongEnergy);
				jInpK0S.set_user_index(iPart);
				L0AsK0SParticles.push_back(jInpK0S);
			}
			else if (part.isCharged()) {
				chParticles.push_back(jInp);
				hTrack->Fill(jInp.pt(), jInp.eta(), jInp.phi());
			}
		} // Particle loop

		// Skip events without V0s
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

		std::vector <fastjet::PseudoJet> V0Inputs = chParticles;
		V0Inputs.insert(V0Inputs.end(), K0SParticles.begin(), K0SParticles.end());
		V0Inputs.insert(V0Inputs.end(), L0Particles.begin(), L0Particles.end());
		fastjet::ClusterSequenceArea clustSeqV0(V0Inputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveV0Jets = clustSeqV0.inclusive_jets();
		std::vector <fastjet::PseudoJet> V0Jets = sorted_by_pt(inclusiveV0Jets); // Sort jets from high to low pt

		std::vector <fastjet::PseudoJet> K0inputs = chParticles;
		K0inputs.insert(K0inputs.end(), K0SParticles.begin(), K0SParticles.end());
		K0inputs.insert(K0inputs.end(), L0AsK0SParticles.begin(), L0AsK0SParticles.end());
		fastjet::ClusterSequenceArea clustSeqK0(K0inputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveK0Jets = clustSeqK0.inclusive_jets();
		std::vector <fastjet::PseudoJet> K0Jets = sorted_by_pt(inclusiveK0Jets); // Sort jets from high to low pt

		for (auto& jet : V0Jets) {
			if (jet.pt() < jetptcut) { continue; }
			bool jetContainsV0s = false;
			for (const fastjet::PseudoJet &constituent : jet.constituents()) {
				int idx = constituent.user_index();
				const Particle &part = pythia.event[idx];

				if (abs(part.id()) == 310) {
					jetContainsV0s = true;
					hzV0_K0S->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, constituent));
					hzV0_V0->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, constituent));
				}
				if (abs(part.id()) == 3122) {
					jetContainsV0s = true;
					hzV0_Lambda0->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, constituent));
					hzV0_V0->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, constituent));
				}
			} // constituent loop
			if (!jetContainsV0s) { continue; }
			hV0Jet->Fill(jet.pt(), jet.eta(), jet.phi());
		} // V0 jet loop

		for (auto& jet : K0Jets) {
			if (jet.pt() < jetptcut) { continue; }
			bool jetContainsV0s = false;
			for (const fastjet::PseudoJet &constituent : jet.constituents()) {
				int idx = constituent.user_index();
				const Particle &part = pythia.event[idx];

				if (abs(part.id()) == 310) {
					jetContainsV0s = true;
					hzK0_K0S->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, constituent));
					hzK0_V0->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, constituent));
				}
				if (abs(part.id()) == 3122) {
					jetContainsV0s = true;
					hzK0_Lambda0->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, constituent));
					hzK0_V0->Fill(jet.pt(), jet.eta(), jet.phi(), z(jet, constituent));
				}
			} // constituent loop
			if (!jetContainsV0s) { continue; }
			hK0Jet->Fill(jet.pt(), jet.eta(), jet.phi());
		}
	}

	hNEvts->Write(hNEvts->GetName(), TObject::kOverwrite);
	hV0->Write(hV0->GetName(), TObject::kOverwrite);
	hTrack->Write(hTrack->GetName(), TObject::kOverwrite);

	hV0Jet->Write(hV0Jet->GetName(), TObject::kOverwrite);
	hK0Jet->Write(hK0Jet->GetName(), TObject::kOverwrite);

	hzV0_K0S->Write(hzV0_K0S->GetName(), TObject::kOverwrite);
	hzV0_Lambda0->Write(hzV0_Lambda0->GetName(), TObject::kOverwrite);
	hzV0_V0->Write(hzV0_V0->GetName(), TObject::kOverwrite);
	hzK0_K0S->Write(hzK0_K0S->GetName(), TObject::kOverwrite);
	hzK0_Lambda0->Write(hzK0_Lambda0->GetName(), TObject::kOverwrite);
	hzK0_V0->Write(hzK0_V0->GetName(), TObject::kOverwrite);

	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();

  return 0;
}
