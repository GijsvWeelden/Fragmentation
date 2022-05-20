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
#include "TH2F.h"

#define nEvents 20000

using namespace Pythia8;

std::vector <fastjet::PseudoJet> do_jet_finding(std::vector <fastjet::PseudoJet> fjInputs,
																								double max_eta_track, double max_eta_jet, double jetR);
int do_matching(double eta, double eta_base, double phi, double phi_base, double matchDist);
void fill_fragmentation(double px, double py, double pz, int id,
												double px_base, double py_base, double pz_base, double p2_base,
												TH1F* hFrag, std::vector<int> PDG);

int main(int /*argc*/, char** /*argv*/)
{
//PYTHIA SETTINGS
	TString name;
	int mecorr=1;
	Float_t ptHatMin = 80;
	Float_t ptHatMax = 200;

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

	// Settings for tracks and jets
	float max_eta_track = 2, min_pt_track = 0., max_pt_track = 10.;
	float max_eta_jet = 2.0, min_pt_jet = 10, max_pt_jet = 200, jetR = 0.4;
	int nBins_eta_track = 40, nBins_pt_track = 50, nBins_eta_jet = 40, nBins_pt_jet = 200;
	int nCharged = 3, nNeutral = 5;
	std::vector<int> PDG = {211, 321, 2212, 111, 130, 310, 311, 3122};
	std::vector<string> Hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};

	// Output histograms
	TFile* outFile = new TFile("PythiaResult.root","RECREATE");
	TH2F *hEtaPt = new TH2F("hEtaPt","Pt vs Eta for all particles;#eta;p_{T} (GeV/c)",
													nBins_eta_track, -1 * max_eta_track, max_eta_track,
													nBins_pt_track, min_pt_track, max_pt_track);
													// 40, -2, 2, 50, 0, 10);
	TH2F *hJetEtaPt = new TH2F("hJetEtaPt","Jet Pt vs Eta;#eta;p^{jet}_{T} (GeV/c)",
														 nBins_eta_jet, -1 * max_eta_jet, max_eta_jet,
														 nBins_pt_jet, min_pt_jet, max_pt_jet);
														//  40, -2, 2, 200, 0, 200);
	TH1F* hists[nCharged + nNeutral];
	TH1F* frags[nCharged + nNeutral];
	TH1F* jetFrags[nCharged + nNeutral];
	for (int i = 0; i < nCharged + nNeutral; i++){
		TH1F* hPt = new TH1F(TString::Format("hPt_%s", Hadrons[i].c_str()).Data(),
												 TString::Format("%s Pt;p_{T} (GeV/c)", Hadrons[i].c_str()).Data(),
												 nBins_pt_track, min_pt_track, max_pt_track);
												// 50, 0, 10);
		hists[i] = hPt;
		TH1F* hFrag = new TH1F(TString::Format("hFrag_%s", Hadrons[i].c_str()).Data(),
													 TString::Format("D(z) for all %s;z", Hadrons[i].c_str()).Data(),
													 100, 0., 1.);
		frags[i] = hFrag;
		TH1F* hJetFrag = new TH1F(TString::Format("hFrag_%s", Hadrons[i].c_str()).Data(),
															TString::Format("D(z) for all %s;z", Hadrons[i].c_str()).Data(),
															100, 0., 1.);
		jetFrags[i] = hJetFrag;
	}

	//Begin event loop
	for (int iEvent = 0; iEvent < nEvents; iEvent++){
		if (!pythia.next()) continue;
		double fourvec[4];

		Double_t ptSumPythia = 0;
		Int_t nPartPythia = 0;
		int nPart = pythia.event.size();

    int a = 0; int b = 0; int c = 0; // For counting descendants
		Int_t nMatriarchs = 0;
		Int_t matriarch1Index = -1;
		Int_t matriarch2Index = -1;
    double pxM1 = 0, pyM1 = 0, pzM1 = 0, p2M1 = 0, etaM1 = 0, phiM1 = 0;
		double pxM2 = 0, pyM2 = 0, pzM2 = 0, p2M2 = 0, etaM2 = 0, phiM2 = 0;
		double matchDist = 1.;

		std::vector<fastjet::PseudoJet> fastjetInputs;

		for (int iPart = 0; iPart < nPart; iPart++){
    	const Particle &part = pythia.event[iPart];
			if (part.status() == -23){
				nMatriarchs++;
        if (nMatriarchs == 1){
          matriarch1Index = part.index();
					etaM1 = part.eta();
					phiM1 = part.phi();
          pxM1 = part.px();
          pyM1 = part.py();
          pzM1 = part.pz();
          p2M1 = part.pAbs2();
        }
        else if (nMatriarchs == 2){
          matriarch2Index = part.index();
					etaM2 = part.eta();
					phiM2 = part.phi();
          pxM2 = part.px();
          pyM2 = part.py();
          pzM2 = part.pz();
          p2M2 = part.pAbs2();
        }
        // else if (nMatriarchs > 2){
					// cout << "Warning: More than 2 outgoing particles found from the initial hard scattering. We will ignore these." << endl;
				// }
				continue;
			}
			if (!part.isFinal()) continue; // No decays yet
			if (abs( part.eta() ) > max_eta_track || part.pT() < min_pt_track) continue;
			hEtaPt->Fill(part.eta(),part.pT());
			for (int i = 0; i < nCharged + nNeutral; i++){
				if (part.id() == PDG[i]){
					hists[i]->Fill(part.pT());
				}
			}
			double px = part.px();
			double py = part.py();
			double pz = part.pz();
			double p2 = part.pAbs2();

			fastjet::PseudoJet jInp(px, py, pz, part.e());
			jInp.set_user_index(part.id());
			fastjetInputs.push_back(jInp);

			int match = do_matching(part.eta(), etaM1, etaM2, part.phi(), phiM1, phiM2);
			if (match == 1){
				fill_fragmentation(px, py, pz, part.id(), pxM1, pyM1, pzM1, p2M1);
			}
			else if (match == 2){
				fill_fragmentation(px, py, pz, part.id(), pxM2, pyM2, pzM2, p2M2);
			}

			// double deltaR1 = (part.eta() - etaM1) * (part.eta() - etaM1) + (part.phi() - phiM1) + (part.phi() - phiM1);
			// double deltaR2 = (part.eta() - etaM2) * (part.eta() - etaM2) + (part.phi() - phiM2) + (part.phi() - phiM2);
			// if (deltaR1 < matchDist && deltaR1 < deltaR2){
			// 	a++;
			// 	double z = (px * pxM1 + py * pyM1 + pz * pzM1)/p2M1;
			// 	for (int i = 0; i < nCharged + nNeutral; i++){
			// 		if (part.id() == PDG[i]){
			// 			frags[i]->Fill(z);
			// 		}
			// 	}
			// }
      // else if (deltaR2 < matchDist){
			// 	b++;
			// 	double z = (px * pxM2 + py * pyM2 + pz * pzM2)/p2M2;
			// 	for (int i = 0; i < nCharged + nNeutral; i++){
			// 		if (part.id() == PDG[i]){
			// 			frags[i]->Fill(z);
			// 		}
			// 	}
			// }
      // else if (deltaR1 > matchDist && deltaR2 > matchDist){
			// 	c++;
			// }
			nPartPythia++;
		}
		if ((iEvent%1000)==0){
			cout << "Pythia event: " << nPartPythia << " particles" << endl;
		}
		// Do jet finding and analysis here.
		std::vector <fastjet::PseudoJet> ptSortedJets = do_jet_finding(fjInputs, max_eta_track, max_eta_jet, jetR);
	}
	//End event loop
	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();
}

std::vector <fastjet::PseudoJet> do_jet_finding(std::vector <fastjet::PseudoJet> fjInputs, double max_eta_track, double max_eta_jet, double jetR){
	fastjet::GhostedAreaSpec ghostSpec(max_eta_track, 1, 0.01);
	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; // need E scheme for jet mass
	fastjet::AreaType areaType = fastjet::active_area;
	fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType, ghostSpec);
	fastjet::AreaDefinition areaDefShape = fastjet::AreaDefinition(areaType, ghostSpec);
		//fastjet::active_area, ghostSpec);
	fastjet::RangeDefinition range(-1. * max_eta_jet, max_eta_jet, 0, 2.*fastjet::pi);
	fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR, recombScheme, strategy);
	fastjet::ClusterSequenceArea clustSeq(fjInputs, jetDef, areaDef);

	std::vector <fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets();
	vector <fastjet::PseudoJet> ptSortedJets = sorted_by_pt(inclusiveJets);
	return ptSortedJets;
}

int do_matching(double eta, double etaM1, double etaM2, double phi, double phiM1, double phiM2, double matchDist){
	double deltaR1 = (eta() - etaM1) * (eta() - etaM1) + (phi() - phiM1) + (phi() - phiM1);
	double deltaR2 = (eta() - etaM2) * (eta() - etaM2) + (phi() - phiM2) + (phi() - phiM2);
	if (deltaR1 < matchDist && deltaR1 < deltaR2){
		return 1;
	}
	else if (deltaR2 < matchDist){
		return 2;
	}
	else return 0;
}

void fragmentation(double px, double py, double pz, int id, double px_base, double py_base, double pz_base, double p2_base, TH1F* hFrag, std::vector<int> PDG){
	double z = px * px_base + py * py_base + pz * pz_base;
	z /= p2_base;
	for (int i = 0; i < PDG.size(); i++){
		if (id == PDG[i]){
			hFrag->Fill(z);
			return;
		}
	}
}