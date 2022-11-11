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

// #define nEvents 200

using namespace Pythia8;

std::vector <fastjet::PseudoJet> do_jet_finding(std::vector <fastjet::PseudoJet> &fastjetInputs,
																								double max_eta_track, double max_eta_jet, double jetR);
int do_matching(double eta, double etaM1, double etaM2, double phi, double phiM1, double phiM2, double matchDist);
void fill_fragmentation(double px, double py, double pz, int id,
												double px_base, double py_base, double pz_base, double p2_base,
												std::vector<TH1F*> &frags, std::vector<int> &PDG);
void fill_fragmentation(const fastjet::PseudoJet &jet, std::vector<TH2F*> &jetFrags, std::vector<int> &PDG);

int main(int argc, char** argv)
{
	int nEvents = 200;
	string outName = "PythiaResult";
	Float_t ptHatMin = 80;
	Float_t ptHatMax = 200;

	if (argc >= 2){
		nEvents = atoi(argv[1]);
		if (nEvents == 0){
			cerr << "ERROR: Zero events requested. Aborting.";
			return 1;
		}
	}
	if (argc >= 3){
		outName = argv[2];
	}
	if (argc >= 4){
		ptHatMin = atof(argv[3]);
	}
	if (argc >= 5){
		ptHatMax = atof(argv[4]);
	}
	if (argc >= 6){
		cout << "Superfluous arguments: ";
		for (int i = 5; i < argc; i++) cout << argv[i] << " ";
		cout << endl;
	}

	//PYTHIA SETTINGS
	TString name;
	int mecorr=1;
	// Float_t ptHatMin = 80;
	// Float_t ptHatMax = 200;

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

	//ME corrections
	//use of matrix corrections where available
	if(mecorr==0){
		pythia.readString("TimeShower:MECorrections=off");
	}
	pythia.init();

	// Settings for tracks and jets
	float max_eta_track = 2, min_pt_track = 0., max_pt_track = 10.;
	float jetR = 0.4, min_pt_jet = 10, max_pt_jet = 200;
	float max_eta_jet = max_eta_track - jetR;
	float min_z = -1e-3, max_z = 1.001;
	int nBins_eta_track = 40, nBins_pt_track = 50, nBins_eta_jet = 40, nBins_pt_jet = 190, nBins_z = 100;
	int match_0 = 0, match_1 = 0, match_2 = 0;
	int nGluons = 0, nQuarks = 0;
	std::vector<int> PDG = {211, 321, 2212, 111, 130, 310, 311, 3122};
	std::vector<string> Hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};
	std::vector<string> Partons = {"g", "q"};

	// Output histograms
	TFile* outFile = new TFile(TString::Format("%s_pthat%.0f_%.0f.root",
																						 outName.c_str(), ptHatMin, ptHatMax).Data(),
																						 "RECREATE");

	TH2F *hEtaPt = new TH2F("hEtaPt","Pt vs Eta for all particles;#eta;p_{T} (GeV/c)",
													nBins_eta_track, -1 * max_eta_track, max_eta_track,
													nBins_pt_track, min_pt_track, max_pt_track);
	TH2F *hJetEtaPt = new TH2F("hJetEtaPt","Jet Pt vs Eta;#eta;p^{jet}_{T} (GeV/c)",
														 nBins_eta_jet, -1 * max_eta_jet, max_eta_jet,
														 nBins_pt_jet, min_pt_jet, max_pt_jet);
	TH1F *hDeltaPartonJet = new TH1F("hDeltaPartonJet","Distance between parton and jet",
																	 nBins_eta_jet, 0., 1.);
	TH1F *hNPartons = new TH1F("hNPartons","Matriarchs per eta;#eta",
														 3 * nBins_eta_jet, -3 * max_eta_jet, 3 * max_eta_jet);
	TH1F *hNJets = new TH1F("hNJets","Jets per eta;#eta",
													 3 * nBins_eta_jet, -3 * max_eta_jet, 3 * max_eta_jet);
	TH2F *hNJetTypes = new TH2F("hNJetTypes","Number of gluon/quark jets;;p^{jet}_{T}",
													 Partons.size(), -0.5, Partons.size()-0.5,
													 nBins_pt_jet, min_pt_jet, max_pt_jet);
	for (int i = 0; i < Partons.size(); i++){
		hNJetTypes->GetXaxis()->SetBinLabel(i+1, Partons[i].c_str());
	}
	std::vector<TH1F*> hists;
	std::vector<TH1F*> frags;
	std::vector<TH2F*> jetFrags;
	std::vector<vector<TH2F*>> partonFrags;

	// Hadron specific spectra, fragmentation, jet fragmentation
	for (int i = 0; i < Hadrons.size(); i++){
		TH1F* hPt = new TH1F(TString::Format("hPt_%s", Hadrons[i].c_str()).Data(),
												 TString::Format("%s Pt;p_{T} (GeV/c)", Hadrons[i].c_str()).Data(),
												 nBins_pt_track, min_pt_track, max_pt_track);
		hists.push_back(hPt);
		TH1F* hFrag = new TH1F(TString::Format("hFrag_%s", Hadrons[i].c_str()).Data(),
													 TString::Format("D^{%s}(z);z", Hadrons[i].c_str()).Data(),
													nBins_z, min_z, max_z);
		frags.push_back(hFrag);
		TH2F* hJetFrag = new TH2F(TString::Format("hJetFrag_%s", Hadrons[i].c_str()).Data(),
															TString::Format("D^{%s}(z);z;p_{T}", Hadrons[i].c_str()).Data(),
															nBins_z, min_z, max_z,
															nBins_pt_jet, min_pt_jet, max_pt_jet);
		jetFrags.push_back(hJetFrag);
	}
	// Parton specific hadron
	for (int j = 0; j < Partons.size(); j++){
		std::vector<TH2F*> tmp;
		for (int i = 0; i < Hadrons.size(); i++){
			TH2F* hPartonFrags = new TH2F(TString::Format("h%sFrags_%s", Partons[j].c_str(), Hadrons[i].c_str()).Data(),
																		TString::Format("D^{%s/%s}(z);z;p_{T}^{jet}",
																										Hadrons[i].c_str(), Partons[j].c_str()).Data(),
																		nBins_z, min_z, max_z,
																		nBins_pt_jet, min_pt_jet, max_pt_jet);
			tmp.push_back(hPartonFrags);
		}
		partonFrags.push_back(tmp);
	}

	//Begin event loop
	for (int iEvent = 0; iEvent < nEvents; iEvent++){
		if (!pythia.next()) continue;
		double fourvec[4];

		Double_t ptSumPythia = 0;
		Int_t nPartPythia = 0;
		int nPart = pythia.event.size();

    int a = 0; int b = 0; int c = 0; // For counting descendants
		Int_t nMatriarchs = 0; Int_t flavourM1 = -999; Int_t flavourM2 = -999;
		int nMatchedJets = 0; int jetMatch1 = 0; int jetMatch2 = 0;
    double pxM1 = 0, pyM1 = 0, pzM1 = 0, p2M1 = 0, etaM1 = 0, phiM1 = 0;
		double pxM2 = 0, pyM2 = 0, pzM2 = 0, p2M2 = 0, etaM2 = 0, phiM2 = 0;
		double matchDist = 1.;

		std::vector<fastjet::PseudoJet> fastjetInputs;

		for (int iPart = 0; iPart < nPart; iPart++){
    	const Particle &part = pythia.event[iPart];
			// Initial scattering partons
			if (part.status() == -23){
				nMatriarchs++;
        if (nMatriarchs == 1){
          flavourM1 = part.id();
					etaM1 = part.eta();
					phiM1 = part.phi();
          pxM1 = part.px();
          pyM1 = part.py();
          pzM1 = part.pz();
          p2M1 = part.pAbs2();
					hNPartons->Fill(etaM1);
        }
        else if (nMatriarchs == 2){
          flavourM2 = part.id();
					etaM2 = part.eta();
					phiM2 = part.phi();
          pxM2 = part.px();
          pyM2 = part.py();
          pzM2 = part.pz();
          p2M2 = part.pAbs2();
					hNPartons->Fill(etaM2);
        }
				continue;
			}
			if (!part.isFinal()) continue;
			if (abs( part.eta() ) > max_eta_track || part.pT() < min_pt_track) continue;
			hEtaPt->Fill(part.eta(),part.pT());
			for (int i = 0; i < Hadrons.size(); i++){
				if (abs(part.id()) == PDG[i]){
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

			//  Check if particle came from parton
			int partMatch = do_matching(part.eta(), etaM1, etaM2, part.phi(), phiM1, phiM2, matchDist);
			if (partMatch == 1){
				fill_fragmentation(px, py, pz, part.id(), pxM1, pyM1, pzM1, p2M1, frags, PDG);
			}
			else if (partMatch == 2){
				fill_fragmentation(px, py, pz, part.id(), pxM2, pyM2, pzM2, p2M2, frags, PDG);
			}
			nPartPythia++;
		} // Particle loop

		// Do jet finding and analysis here.
		fastjet::GhostedAreaSpec ghostSpec(max_eta_track, 1, 0.01);
		fastjet::Strategy strategy = fastjet::Best;
		fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; // need E scheme for jet mass
		fastjet::AreaType areaType = fastjet::active_area;
		fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::AreaDefinition areaDefShape = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::RangeDefinition range(-1. * max_eta_jet, max_eta_jet, 0, 2.*fastjet::pi);
		fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR, recombScheme, strategy);
		fastjet::ClusterSequenceArea clustSeq(fastjetInputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets();
		vector <fastjet::PseudoJet> ptSortedJets = sorted_by_pt(inclusiveJets); // Sort jets from high to low pt

		for (auto jet : ptSortedJets){
      if (jet.pt() < min_pt_jet) continue;
			if (nMatchedJets == 2) continue;
			hJetEtaPt->Fill(jet.eta(), jet.pt());
			// Check if jet came from parton
			int jetMatch = do_matching(jet.eta(), etaM1, etaM2, jet.phi(), phiM1, phiM2, matchDist, hDeltaPartonJet);
			if (jetMatch == 1){
				if (jetMatch1) continue;
				hNJets->Fill(jet.eta());
				fill_fragmentation(jet, jetFrags, PDG);
				if (flavourM1 == 21){ // Gluon
					fill_fragmentation(jet, partonFrags[0], PDG);
					hNJetTypes->Fill(1, jet.pt());
					nGluons++;
				}
				else{ // Quark
					fill_fragmentation(jet, partonFrags[1], PDG);
					hNJetTypes->Fill(2, jet.pt());
					nQuarks++;
				}
				jetMatch1 = 1;
				nMatchedJets++;
			}
			else if (jetMatch == 2){
				if (jetMatch2) continue;
				hNJets->Fill(jet.eta());
				fill_fragmentation(jet, jetFrags, PDG);
				if (flavourM2 == 21){ // Gluon
					fill_fragmentation(jet, partonFrags[0], PDG);
					hNJetTypes->Fill(1, jet.pt());
					nGluons++;
				}
				else{ // Quark
					fill_fragmentation(jet, partonFrags[1], PDG);
					hNJetTypes->Fill(2, jet.pt());
					nQuarks++;
				}
				jetMatch2 = 1;
				nMatchedJets++;
			}
		}
		if (nMatchedJets == 2){
			match_2++;
		}
		else if (nMatchedJets == 1){
			match_1++;
		}
		else if (nMatchedJets == 0){
			match_0++;
		}
	}
	// hNJetTypes->SetBinContent(1, 1.*nGluons);
	// hNJetTypes->SetBinContent(2, 1.*nQuarks);

	cout << "Number of events: " << nEvents << endl
		<< "Events with (2, 1, 0) matches:" << endl
		<< "2: " << match_2 << " ("	<< 1.*match_2/nEvents << ")" << endl
		<< "1: " << match_1 << " ("	<< 1.*match_1/nEvents << ")" << endl
		<< "0: " << match_0 << " ("	<< 1.*match_0/nEvents << ")" << endl
		<< "Number of gluon jets: " << nGluons << " (" << nGluons/(2. * match_2 + 1. * match_1) << ")" << endl
		<< "Number of quark jets: " << nQuarks << " (" << nQuarks/(2. * match_2 + 1. * match_1) << ")" << endl;
	//End event loop
	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();
}

std::vector <fastjet::PseudoJet> do_jet_finding(std::vector <fastjet::PseudoJet> &fastjetInputs, double max_eta_track, double max_eta_jet, double jetR){
	fastjet::GhostedAreaSpec ghostSpec(max_eta_track, 1, 0.01);
	fastjet::Strategy strategy = fastjet::Best;
	fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; // need E scheme for jet mass
	fastjet::AreaType areaType = fastjet::active_area;
	fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType, ghostSpec);
	fastjet::AreaDefinition areaDefShape = fastjet::AreaDefinition(areaType, ghostSpec);
		//fastjet::active_area, ghostSpec);
	fastjet::RangeDefinition range(-1. * max_eta_jet, max_eta_jet, 0, 2.*fastjet::pi);
	fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR, recombScheme, strategy);
	fastjet::ClusterSequenceArea clustSeq(fastjetInputs, jetDef, areaDef);

	std::vector <fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets();
	vector <fastjet::PseudoJet> ptSortedJets = sorted_by_pt(inclusiveJets);
	return ptSortedJets;
}

int do_matching(double eta, double etaM1, double etaM2, double phi, double phiM1, double phiM2, double matchDist, TH1F* hDeltaPartonJet = nullptr){
	double dphi1 = phi - phiM1;
	if (dphi1 < - fastjet::pi) dphi1 += fastjet::twopi;
	if (dphi1 > fastjet::pi) dphi1 -= fastjet::twopi;
	double dphi2 = phi - phiM2;
	if (dphi2 < - fastjet::pi) dphi2 += fastjet::twopi;
	if (dphi2 > fastjet::pi) dphi2 -= fastjet::twopi;
	double deltaR1 = (eta - etaM1) * (eta - etaM1) + dphi1 * dphi1;
	double deltaR2 = (eta - etaM2) * (eta - etaM2) + dphi2 * dphi2;
	deltaR1 = sqrt(deltaR1);
	deltaR2 = sqrt(deltaR2);
	if (deltaR1 < matchDist && deltaR1 < deltaR2){
		if (hDeltaPartonJet) hDeltaPartonJet->Fill(deltaR1);
		return 1;
	}
	else if (deltaR2 < matchDist){
		if (hDeltaPartonJet) hDeltaPartonJet->Fill(deltaR2);
		return 2;
	}
	else return 0;
}

void fill_fragmentation(double px, double py, double pz, int id, double px_base, double py_base, double pz_base, double p2_base, std::vector<TH1F*> &frags, std::vector<int> &PDG){
	double z = px * px_base + py * py_base + pz * pz_base;
	z /= p2_base;
	for (int i = 0; i < PDG.size(); i++){
		if (abs(id) == PDG[i]){
			frags[i]->Fill(z);
			return;
		}
	}
}

void fill_fragmentation(const fastjet::PseudoJet &jet, std::vector<TH2F*> &jetFrags, std::vector<Int_t> &PDG)
{
  if (!jet.has_constituents())
    return;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  Double_t jpx = jet.px(), jpy = jet.py(), jpz = jet.pz(), jp2 = jet.modp2(); Double_t z;
	double px, py, pz; int id;
  for(UInt_t ic = 0; ic < constits.size(); ++ic){
		px = constits[ic].px();
		py = constits[ic].py();
		pz = constits[ic].pz();
		id = constits[ic].user_index();
		z = px * jpx + py * jpy + pz * jpz;
		z /= jp2;
		for (int i = 0; i < PDG.size(); i++){
			if (abs(id) == PDG[i]){
				jetFrags[i]->Fill(z, jet.perp());
			}
		}
  }
  return;
}
