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
void withDecays(int nEvents);
void withoutDecays(int nEvents);

int main(int argc, char** argv)
{
	int decays = 1;
	int nEvents = 200;
	if (argc >= 2){
		decays = atoi(argv[1]);
	}
	if (argc >= 3){
		nEvents = atoi(argv[2]);
		if (nEvents == 0){
			cerr << "ERROR: Zero events requested. Aborting.";
			return 1;
		}
	}

  if (decays) withDecays(nEvents);
	else withoutDecays(nEvents);

  return 0;
}

// /*
void withDecays(int nEvents)
{
	string outName = "withdecays";
	double ptHatMin = 20;
	double ptHatMax = 80;

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
	// pythia.readString("310:mayDecay  = off"); //K0s
	// pythia.readString("311:mayDecay  = off"); //K0
	// pythia.readString("3122:mayDecay = off"); //labda0

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
	double jetptcut = 10;

	int nBins_eta_track = 40;
	int nBins_pt_track = 100;
	int nBins_delta_pt = 101;
	int nBins_ratio_pt = 100;
	int nBins_eta_jet = 40;
	int nBins_pt_jet = 200;
	int nBins_dR = 100;
	int nBins_z = 100;
  int nBins_dau = 3;
  int nBins_v0 = 21;
	double min_delta_pt = -50.5, max_delta_pt = 50.5;
	double min_pt_track = 0., max_pt_track = 20.;
	double min_ratio_pt = 0., max_ratio_pt = 5;
	double max_eta_jet = max_eta_track - jetR;
	double min_pt_jet = 0, max_pt_jet = 200;
	double min_nv0 = -0.5, max_nv0 = 20.5;
  double min_dau = -0.5, max_dau = 2.5;
	double min_z = 1e-3, max_z = 1.001;
	double min_dR = 0., max_dR = 2.;

  int nDim = 4;
  int hzBins[nDim]   = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_z};
  double hzMin[nDim] = {min_pt_jet, min_ratio_pt, min_pt_track, min_z};
  double hzMax[nDim] = {max_pt_jet, max_ratio_pt, max_pt_track, max_z};

  int hNdauBins[nDim]   = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_dau};
  double hNdauMin[nDim] = {min_pt_jet, min_ratio_pt, min_pt_track, min_dau};
  double hNdauMax[nDim] = {max_pt_jet, max_ratio_pt, max_pt_track, max_dau};

	int hNptv0Bins[nDim]   = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_v0};
  double hNptv0Min[nDim] = {min_pt_jet, min_ratio_pt, min_pt_track, min_nv0};
  double hNptv0Max[nDim] = {max_pt_jet, max_ratio_pt, max_pt_track, max_nv0};

  int hDauPtBins[nDim]   = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_pt_track};
  double hDauPtMin[nDim] = {min_pt_jet, min_ratio_pt, min_pt_track, min_pt_track};
  double hDauPtMax[nDim] = {max_pt_jet, max_ratio_pt, max_pt_track, max_pt_track};

  int hDeltaPtBins[nDim]   = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_delta_pt};
  double hDeltaPtMin[nDim] = {min_pt_jet, min_ratio_pt, min_pt_track, min_delta_pt};
  double hDeltaPtMax[nDim] = {max_pt_jet, max_ratio_pt, max_pt_track, max_delta_pt};

  int hdRBins[nDim]   = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_dR};
  double hdRMin[nDim] = {min_pt_jet, min_ratio_pt, min_pt_track, min_dR};
  double hdRMax[nDim] = {max_pt_jet, max_ratio_pt, max_pt_track, max_dR};

	int hMissDauPtBins[nDim]   = {nBins_pt_jet, nBins_pt_track, nBins_pt_track, nBins_pt_track};
  double hMissDauPtMin[nDim] = {min_pt_jet, min_pt_track, min_pt_track, min_pt_track};
  double hMissDauPtMax[nDim] = {max_pt_jet, max_pt_track, max_pt_track, max_pt_track};

	int hMissDauEtaBins[nDim]   = {nBins_pt_jet, nBins_pt_track, nBins_pt_track, nBins_pt_track};
  double hMissDauEtaMin[nDim] = {min_pt_jet, min_pt_track, -1 * max_eta_track, -1 * max_eta_track};
  double hMissDauEtaMax[nDim] = {max_pt_jet, max_pt_track, max_eta_track, max_eta_track};

	int hMissDaudRBins[nDim+1]   = {nBins_pt_jet, nBins_pt_track, nBins_dR, nBins_dR, nBins_dR};
  double hMissDaudRMin[nDim+1] = {min_pt_jet, min_pt_track, min_dR, min_dR, min_dR};
  double hMissDaudRMax[nDim+1] = {max_pt_jet, max_pt_track, max_dR, max_dR, max_dR};

	double hEventWeightBins[15] = {1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0};

	int match_0 = 0, match_1 = 0, match_2 = 0;
	int nGluons = 0, nQuarks = 0;
	std::vector<int> PDG = {211, 321, 2212, 111, 130, 310, 311, 3122};
	std::vector<string> Hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};
	std::vector<string> Partons = {"g", "q"};

	string sorig   = "#it{p}_{T, jet}^{orig.}";
	string snew    = "#it{p}_{T, jet}^{new}";
	string sratio  = "#it{p}_{T, jet}^{new} / #it{p}_{T, jet}^{orig.}";
	string sv0pt   = "#it{p}_{T, V0}";
	string sz      = "#it{z}_{V0}";
	string sndau   = "#it{N}(V0 daughters #in jet)";
	string sdaupt  = "#it{p}_{T}^{daughter #notin jet}";
	string sdaueta = "#eta^{daughter #notin jet}";
	string sdelta  = "#it{p}_{T, jet}^{orig.} - #it{p}_{T, jet}^{new}";
	string sdr     = "d #it{R}";
	string snv0    = "#it{N}(V0s #in jet)";

	// Output histograms
	TFile* outFile = new TFile(TString::Format("%s_pthat%.0f_%.0f.root",
																						 outName.c_str(), ptHatMin, ptHatMax).Data(),
																						 "RECREATE");

	TH1D* hNEvts = new TH1D("hNEvts", "hNEvts", 1, 0.5, 1.5);
	TH1D* hEvtWeights = new TH1D("hEvtWeights", "hEvtWeights", 14, hEventWeightBins);
	TH1D* hjetpt = new TH1D("hjetpt", "hjetpt; #it{p}_{T, jet}^{orig.}", nBins_pt_jet, min_pt_jet, max_pt_jet);
	TH1D* hnewjetpt = new TH1D("hnewjetpt", "hnewjetpt; #it{p}_{T, jet}^{new}", nBins_pt_jet, min_pt_jet, max_pt_jet);
	TH1D* hv0pt = new TH1D("hv0pt", "hv0pt; #it{p}_{T, v0}", nBins_pt_track, min_pt_track, max_pt_track);

	TH3D* hNv0 = new TH3D("hNv0",
												TString::Format("hNv0; %s; %s; %s",
														sorig.c_str(), sratio.c_str(), snv0.c_str()).Data(),
											  nBins_pt_jet, min_pt_jet, max_pt_jet,
												nBins_ratio_pt, min_ratio_pt, max_ratio_pt,
												nBins_v0, min_nv0, max_nv0
											 );
	THnSparseD* hNptv0
		= new THnSparseD("hNptv0",
										 TString::Format("hNptv0;%s;%s;%s;%s",
												sorig.c_str(), sratio.c_str(), sv0pt.c_str(), snv0.c_str()).Data(),
										 nDim, hNptv0Bins, hNptv0Min, hNptv0Max
		);
	THnSparseD* hzU // z calculated with uncorrected jet pt
		= new THnSparseD(TString::Format("hzU").Data(),
										 TString::Format("hzU;%s;%s;%s;%s",
												sorig.c_str(), sratio.c_str(), sv0pt.c_str(), sz.c_str()).Data(),
                     nDim, hzBins, hzMin, hzMax
    );
	THnSparseD* hzC // z calculated with corrected jet pt
		= new THnSparseD(TString::Format("hzC").Data(),
										 TString::Format("hzC;%s;%s;%s;%s",
												snew.c_str(), sratio.c_str(), sv0pt.c_str(), sz.c_str()).Data(),
                     nDim, hzBins, hzMin, hzMax
    );
	// V0 in jet
	THnSparseD* hNdauV0in
		= new THnSparseD("hNdauV0in",
										 TString::Format("hNdauV0in;%s;%s;%s;%s",
												sorig.c_str(), sratio.c_str(), sv0pt.c_str(), sndau.c_str()).Data(),
										 nDim, hNdauBins, hNdauMin, hNdauMax
										);
	THnSparseD* hDauPtV0in
		= new THnSparseD("hDauPtV0in",
										 TString::Format("hDauPtV0in;%s;%s;%s;%s",
									 		  sorig.c_str(), sratio.c_str(), sv0pt.c_str(), sdaupt.c_str()).Data(),
										 nDim, hDauPtBins, hDauPtMin, hDauPtMax
										);
	THnSparseD* hDeltaPtV0in
		= new THnSparseD("hDeltaPtV0in",
                     TString::Format("hDeltaPtV0in;%s;%s;%s;%s",
									 		  sorig.c_str(), sratio.c_str(), sv0pt.c_str(), sdelta.c_str()).Data(),
										 nDim, hDeltaPtBins, hDeltaPtMin, hDeltaPtMax
										);
	THnSparseD* hdRV0in
		= new THnSparseD("hdRV0in",
                     TString::Format("hdRV0in;%s;%s;%s;%s",
									 		  sorig.c_str(), sratio.c_str(), sv0pt.c_str(), sdr.c_str()).Data(),
										 nDim, hdRBins, hdRMin, hdRMax
										);
	THnSparseD* hMissDauPt
		= new THnSparseD("hMissDauPt",
                     TString::Format("hMissDauPt;%s;%s;%s;%s",
									 		  sorig.c_str(), sv0pt.c_str(), sdaupt.c_str(), sdaupt.c_str()).Data(),
										 nDim, hMissDauPtBins, hMissDauPtMin, hMissDauPtMax
										);
	THnSparseD* hMissDauEta
		= new THnSparseD("hMissDauEta",
                     TString::Format("hMissDauEta;%s;%s;%s;%s",
									 		  sorig.c_str(), sv0pt.c_str(), sdaueta.c_str(), sdaueta.c_str()).Data(),
										 nDim, hMissDauEtaBins, hMissDauEtaMin, hMissDauEtaMax
										);
	THnSparseD* hMissDaudR
		= new THnSparseD("hMissDaudR",
                     TString::Format("hMissDaudR;%s;%s;%s(V0);%s;%s",
									 		  sorig.c_str(), sv0pt.c_str(), sdr.c_str(), sdr.c_str(), sdr.c_str()).Data(),
										 nDim+1, hMissDaudRBins, hMissDaudRMin, hMissDaudRMax
										);

	// V0 out of jet, but daughter(s) in jet
	THnSparseD* hNdauV0out
		= new THnSparseD("hNdauV0out",
										 TString::Format("hNdauV0out;%s;%s;%s;%s",
												sorig.c_str(), sratio.c_str(), sv0pt.c_str(), sndau.c_str()).Data(),
										 nDim, hNdauBins, hNdauMin, hNdauMax
										);
	THnSparseD* hDauPtV0out
		= new THnSparseD("hDauPtV0out",
										 TString::Format("hDauPtV0out;%s;%s;%s;%s",
									 		  sorig.c_str(), sratio.c_str(), sv0pt.c_str(), sdaupt.c_str()).Data(),
										 nDim, hDauPtBins, hDauPtMin, hDauPtMax
										);
	THnSparseD* hDeltaPtV0out
		= new THnSparseD("hDeltaPtV0out",
                     TString::Format("hDeltaPtV0out;%s;%s;%s;%s",
									 		  sorig.c_str(), sratio.c_str(), sv0pt.c_str(), sdelta.c_str()).Data(),
										 nDim, hDeltaPtBins, hDeltaPtMin, hDeltaPtMax
										);
	THnSparseD* hdRV0out
		= new THnSparseD("hdRV0out",
                     TString::Format("hdRV0out;%s;%s;%s;%s",
									 		  sorig.c_str(), sratio.c_str(), sv0pt.c_str(), sdr.c_str()).Data(),
										 nDim, hdRBins, hdRMin, hdRMax
										);

  for (int iEvent = 0; iEvent < nEvents; iEvent++)
  {
		if (!pythia.next()) continue;
		hNEvts->Fill(1);
		double weight = pythia.info.weight();
		hEvtWeights->Fill(weight);
		double fourvec[4];

		double ptSumPythia = 0;
		int nPartPythia = 0;
		int nPart = pythia.event.size();

		std::vector<fastjet::PseudoJet> fastjetInputs;
		std::vector<int> v0Indices;

		for (int iPart = 0; iPart < nPart; iPart++)
    {
    	const Particle &part = pythia.event[iPart];
			if (abs( part.eta() ) > max_eta_track || part.pT() < trackptcut) continue;

			// Save K0S and (anti)Lambda0
			if ( (abs(part.id()) == 310) || (abs(part.id()) == 3122) ){
				vector<int> dList = part.daughterList();
				const Particle &d0 = pythia.event[dList[0]];
				const Particle &d1 = pythia.event[dList[1]];
				if (!d0.isFinal() || !d1.isFinal()) {
					pythia.event.list();
					cout << "V0 daughter is not final state particle!" << endl;
					cout << "V0 = " << iPart << ", d1 = " << dList[0] << ", d2 = " << dList[1] << endl;
					return;
				}
				if (abs( d0.eta() ) > max_eta_track || d0.pT() < trackptcut || !d0.isFinal()) continue;
				if (abs( d1.eta() ) > max_eta_track || d1.pT() < trackptcut || !d1.isFinal()) continue;

				v0Indices.push_back(iPart);
				fastjet::PseudoJet pjV0(part.px(), part.py(), part.pz(), part.e());
				hv0pt->Fill(pjV0.pt());
			}
			if (!part.isFinal()) continue;

			fastjet::PseudoJet jInp(part.px(), part.py(), part.pz(), part.e());
			jInp.set_user_index(iPart);
			fastjetInputs.push_back(jInp);
			nPartPythia++;
		} // Particle loop

		// Do jet finding and analysis here.
		fastjet::GhostedAreaSpec ghostSpec(max_eta_track, 1, 0.01);
		fastjet::Strategy strategy = fastjet::Best;
		fastjet::AreaType areaType = fastjet::active_area;
		fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::AreaDefinition areaDefShape = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::RangeDefinition range(-1. * max_eta_jet, max_eta_jet, 0, 2.*fastjet::pi);
		fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; // need E scheme for jet mass
		fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR, recombScheme, strategy);
		fastjet::ClusterSequenceArea clustSeq(fastjetInputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets();
		std::vector <fastjet::PseudoJet> ptSortedJets = sorted_by_pt(inclusiveJets); // Sort jets from high to low pt

		bool foundV0Match = false;
		for (auto& jet : ptSortedJets)
    {
      if (jet.pt() < jetptcut) {
        continue;
      }
			hjetpt->Fill(jet.pt());
			fastjet::PseudoJet newjet(jet);
			vector<int> particlesToAdd;
			vector<int> particlesToSubtract;
			vector<int> dauInJet;
			int nv0InJet = 0;

      for (int v0Index : v0Indices)
      {
        const Particle &v0 = pythia.event[v0Index];
				fastjet::PseudoJet pjV0(v0.px(), v0.py(), v0.pz(), v0.e());
				bool v0inJet = (pjV0.delta_R(jet) <= jetR);
				if (v0inJet) {
					nv0InJet++;
				}
        vector<int> dList = v0.daughterList();
        vector<bool> dInJet = { false, false };
				int ndinjet = 0;

        for (const fastjet::PseudoJet &constituent : jet.constituents()) {
          int idx = constituent.user_index();
          if (idx == dList[0]) {
            dInJet[0] = true;
						ndinjet++;
            // cout << "Daughter found! Particle " << idx << " is daughter of v0 " << v0Index << "(" << dList[0] << ", " << dList[1] << ")" << endl;
          }
          if (idx == dList[1]) {
            dInJet[1] = true;
						ndinjet++;
            // cout << "Daughter found! Particle " << idx << " is daughter of v0 " << v0Index << "(" << dList[0] << ", " << dList[1] << ")" << endl;
          }
        } // constituent loop
				dauInJet.push_back(ndinjet);

				// Correct jet for missing/superfluous daughters
        for (int i = 0; i < 2; i++) {
          int dIndex = dList[i];
          const Particle& dp = pythia.event[dIndex];
          fastjet::PseudoJet daughter(dp.px(), dp.py(), dp.pz(), dp.e());
          if (v0inJet && !dInJet[i]) {
            newjet += daughter;
          }
          if (dInJet[i] && !v0inJet) {
            newjet -= daughter;
          }
        }
			} // v0 loop

			// Jet momentum is corrected for missing/superfluous v0 daughters
			double ptratio = newjet.pt() / jet.pt();
			hnewjetpt->Fill(newjet.pt());
			hNv0->Fill(jet.pt(), ptratio, nv0InJet);

			for (int iv0 = 0; iv0 < v0Indices.size(); iv0++)
			{
        const Particle &v0 = pythia.event[v0Indices[iv0]];
				fastjet::PseudoJet pjV0(v0.px(), v0.py(), v0.pz(), v0.e());
				double dR = pjV0.delta_R(jet);
				bool v0inJet = dR <= jetR;
				double zU = pjV0.px() * jet.px() + pjV0.py() * jet.py() + pjV0.pz() * jet.pz();
				zU /= (jet.px() * jet.px() + jet.py() * jet.py() + jet.pz() * jet.pz());
				double zC = pjV0.px() * newjet.px() + pjV0.py() * newjet.py() + pjV0.pz() * newjet.pz();
				zC /= (newjet.px() * newjet.px() + newjet.py() * newjet.py() + newjet.pz() * newjet.pz());

				int ndinjet = dauInJet[iv0];

				if (v0inJet) {
					hzU->Fill(jet.pt(), ptratio, pjV0.pt(), zU);
					hzC->Fill(jet.pt(), ptratio, pjV0.pt(), zC);
					hNptv0->Fill(jet.pt(), ptratio, pjV0.pt(), nv0InJet);
					hNdauV0in->Fill(jet.pt(), ptratio, pjV0.pt(), ndinjet);
					hDeltaPtV0in->Fill(jet.pt(), ptratio, pjV0.pt(), jet.pt() - newjet.pt());
					hdRV0in->Fill(jet.pt(), ptratio, pjV0.pt(), dR);

					if (ndinjet == 0) {
						vector<int> dList = v0.daughterList();
						const Particle &d0 = pythia.event[dList[0]];
						const Particle &d1 = pythia.event[dList[1]];
						fastjet::PseudoJet pjd0(d0.px(), d0.py(), d0.pz(), d0.e());
						fastjet::PseudoJet pjd1(d1.px(), d1.py(), d1.pz(), d1.e());

						hMissDauPt->Fill(jet.pt(), pjV0.pt(), pjd0.pt(), pjd1.pt());
						hMissDauEta->Fill(jet.pt(), pjV0.pt(), pjd0.eta(), pjd1.eta());
						hMissDaudR->Fill(jet.pt(), pjV0.pt(), dR, pjd0.delta_R(jet), pjd1.delta_R(jet));

						// cout << "v0 " << v0Indices[iv0] << " has daughters " << dList[0] << ", " << dList[1] << endl;
						// cout << "jet: " << jet.px() << ", " << jet.py() << ", " << jet.pz() << endl;

						// cout << "(pt, eta, phi):" << endl
						// 		 << "jet (" << jet.pt() << ", " << jet.eta() << ", " << jet.phi() << ")" << endl
						// 		 << "v0 (" << pjV0.pt() << ", " << pjV0.eta() << ", " << pjV0.phi() << "), dR(jet) = " << pjV0.delta_R(jet) << endl
						// 		 << "d0 (" << pjd0.pt() << ", " << pjd0.eta() << ", " << pjd0.phi() << "), dR(jet) = " << pjd0.delta_R(jet) << endl
						// 		 << "d1 (" << pjd1.pt() << ", " << pjd1.eta() << ", " << pjd1.phi() << "), dR(jet) = " << pjd1.delta_R(jet) << endl;
						// pythia.event.list();
					}
				}
				if (!v0inJet && ndinjet > 0) {
					hNdauV0out->Fill(jet.pt(), ptratio, pjV0.pt(), ndinjet);
					hDeltaPtV0out->Fill(jet.pt(), ptratio, pjV0.pt(), jet.pt() - newjet.pt());
					hdRV0out->Fill(jet.pt(), ptratio, pjV0.pt(), dR);
				}
			} // v0 loop
		} // jet loop
	} //End event loop

	// hjetpt->Write(hjetpt->GetName(), TObject::kOverwrite);
	// hv0pt->Write(hv0pt->GetName(), TObject::kOverwrite);
	hNdauV0in->Write(hNdauV0in->GetName(), TObject::kOverwrite);
	hDeltaPtV0in->Write(hDeltaPtV0in->GetName(), TObject::kOverwrite);
	hdRV0in->Write(hdRV0in->GetName(), TObject::kOverwrite);
	hNdauV0out->Write(hNdauV0out->GetName(), TObject::kOverwrite);
	hDeltaPtV0out->Write(hDeltaPtV0out->GetName(), TObject::kOverwrite);
	hdRV0out->Write(hdRV0out->GetName(), TObject::kOverwrite);
	hzU->Write(hzU->GetName(), TObject::kOverwrite);
	hzC->Write(hzC->GetName(), TObject::kOverwrite);
	hMissDauPt->Write(hMissDauPt->GetName(), TObject::kOverwrite);
	hMissDauEta->Write(hMissDauEta->GetName(), TObject::kOverwrite);
	hMissDaudR->Write(hMissDaudR->GetName(), TObject::kOverwrite);
	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();
}
// */

// /*
void withoutDecays(int nEvents)
{
	string outName = "withoutdecays";
	double ptHatMin = 20;
	double ptHatMax = 80;

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
	double jetptcut = 10.;

	int nBins_eta_track = 40;
	int nBins_pt_track = 100;
	int nBins_delta_pt = 101;
	int nBins_ratio_pt = 100;
	int nBins_eta_jet = 40;
	int nBins_pt_jet = 190;
	int nBins_dR = 100;
	int nBins_z = 100;
  int nBins_dau = 3;
	int nBins_v0 = 20;
	double min_delta_pt = -5.05, max_delta_pt = 5.05;
	double min_pt_track = 0., max_pt_track = 20.;
	double min_ratio_pt = 0., max_ratio_pt = 5;
	double max_eta_jet = max_eta_track - jetR;
	double min_pt_jet = .0, max_pt_jet = 200;
	double min_nv0 = -0.5, max_nv0 = 19.5;
  double min_dau = -0.5, max_dau = 2.5;
	double min_z = 1e-3, max_z = 1.001;
	double min_dR = 0., max_dR = 1.;

  int nDim = 4;
  int hzBins[nDim]         = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_z};
  double hzMin[nDim]       = {min_pt_jet, min_ratio_pt, min_pt_track, min_z};
  double hzMax[nDim]       = {max_pt_jet, max_ratio_pt, max_pt_track, max_z};

  int hNdauBins[nDim]      = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_dau};
  double hNdauMin[nDim]    = {min_pt_jet, min_ratio_pt, min_pt_track, min_dau};
  double hNdauMax[nDim]    = {max_pt_jet, max_ratio_pt, max_pt_track, max_dau};

  int hDauPtBins[nDim]     = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_pt_track};
  double hDauPtMin[nDim]   = {min_pt_jet, min_ratio_pt, min_pt_track, min_pt_track};
  double hDauPtMax[nDim]   = {max_pt_jet, max_ratio_pt, max_pt_track, max_pt_track};

  int hDeltaPtBins[nDim]   = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_delta_pt};
  double hDeltaPtMin[nDim] = {min_pt_jet, min_ratio_pt, min_pt_track, min_delta_pt};
  double hDeltaPtMax[nDim] = {max_pt_jet, max_ratio_pt, max_pt_track, max_delta_pt};

  int hdRBins[nDim]        = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_dR};
  double hdRMin[nDim]      = {min_pt_jet, min_ratio_pt, min_pt_track, min_dR};
  double hdRMax[nDim]      = {max_pt_jet, max_ratio_pt, max_pt_track, max_dR};

	int hNv0Bins[nDim]       = {nBins_pt_jet, nBins_ratio_pt, nBins_pt_track, nBins_v0};
  double hNv0Min[nDim]     = {min_pt_jet, min_ratio_pt, min_pt_track, min_nv0};
  double hNv0Max[nDim]     = {max_pt_jet, max_ratio_pt, max_pt_track, max_nv0};

	int match_0 = 0, match_1 = 0, match_2 = 0;
	int nGluons = 0, nQuarks = 0;
	std::vector<int> PDG = {211, 321, 2212, 111, 130, 310, 311, 3122};
	std::vector<string> Hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};
	std::vector<string> Partons = {"g", "q"};

	string sorig  = "#it{p}_{T, jet}^{orig.}";
	string sratio = "#it{p}_{T, jet}^{new} / #it{p}_{T, jet}^{orig.}";
	string sv0pt  = "#it{p}_{T, V0}";
	string sz     = "#it{z}_{V0}";
	string sndau  = "#it{N}(V0 daughters #in jet)";
	string sdaupt = "#it{p}_{T}^{daughter #notin jet}";
	string sdelta = "#it{p}_{T, jet}^{orig.} - #it{p}_{T, jet}^{new}";
	string sdr    = "d#it{R}";

	// Output histograms
	TFile* outFile = new TFile(TString::Format("%s_pthat%.0f_%.0f.root",
																						 outName.c_str(), ptHatMin, ptHatMax).Data(),
																						 "RECREATE");

	TH1D* hNEvts = new TH1D("hNEvts", "hNEvts", 1, 0.5, 1.5);
	TH1F* hjetpt = new TH1F("hjetpt", "hjetpt; #it{p}_{T, jet}", nBins_pt_jet, min_pt_jet, max_pt_jet);
	TH1F* hv0pt = new TH1F("hv0pt", "hv0pt; #it{p}_{T, v0}", nBins_pt_track, min_pt_track, max_pt_track);

	TH3D* hz
		= new TH3D("hz", TString::Format("hz;%s;%s;%s", sorig.c_str(), sv0pt.c_str(), sz.c_str()).Data(),
                     nBins_pt_jet, min_pt_jet, max_pt_jet,
										 nBins_pt_track, min_pt_track, max_pt_track,
										 nBins_z, min_z, max_z
    );
	TH3D* hdRV0in
		= new TH3D("hdRV0in", TString::Format("hdRV0in;%s;%s;%s", sorig.c_str(), sv0pt.c_str(), sz.c_str()).Data(),
                     nBins_pt_jet, min_pt_jet, max_pt_jet,
										 nBins_pt_track, min_pt_track, max_pt_track,
										 nBins_z, min_z, max_z
    );

  for (int iEvent = 0; iEvent < nEvents; iEvent++)
  {
		if (!pythia.next()) continue;
		hNEvts->Fill(1);
		double fourvec[4];

		double ptSumPythia = 0;
		int nPartPythia = 0;
		int nPart = pythia.event.size();

		std::vector<fastjet::PseudoJet> fastjetInputs;
		std::vector<int> v0Indices;

		for (int iPart = 0; iPart < nPart; iPart++)
    {
    	const Particle &part = pythia.event[iPart];
			if (!part.isFinal()) continue;
			if (abs( part.eta() ) > max_eta_track || part.pT() < trackptcut) continue;
			fastjet::PseudoJet jInp(part.px(), part.py(), part.pz(), part.e());
			jInp.set_user_index(iPart);
			fastjetInputs.push_back(jInp);
			nPartPythia++;

			// Save K0S and (anti)Lambda0
			if ( (abs(part.id()) == 310) || (abs(part.id()) == 3122) ){
				v0Indices.push_back(iPart);
				hv0pt->Fill(jInp.pt());
			}
		} // Particle loop

		// Do jet finding and analysis here.
		fastjet::GhostedAreaSpec ghostSpec(max_eta_track, 1, 0.01);
		fastjet::Strategy strategy = fastjet::Best;
		fastjet::AreaType areaType = fastjet::active_area;
		fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::AreaDefinition areaDefShape = fastjet::AreaDefinition(areaType, ghostSpec);
		fastjet::RangeDefinition range(-1. * max_eta_jet, max_eta_jet, 0, 2.*fastjet::pi);
		fastjet::RecombinationScheme recombScheme = fastjet::E_scheme; // need E scheme for jet mass
		fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR, recombScheme, strategy);
		fastjet::ClusterSequenceArea clustSeq(fastjetInputs, jetDef, areaDef);
		std::vector <fastjet::PseudoJet> inclusiveJets = clustSeq.inclusive_jets();
		std::vector <fastjet::PseudoJet> ptSortedJets = sorted_by_pt(inclusiveJets); // Sort jets from high to low pt

		bool foundV0Match = false;
		for (auto& jet : ptSortedJets)
    {
      if (jet.pt() < jetptcut) {
        continue;
      }
			hjetpt->Fill(jet.pt());

			for (const fastjet::PseudoJet &constituent : jet.constituents())
			{
				int idx = constituent.user_index();
				const Particle &part = pythia.event[idx];
				if ( (abs(part.id()) == 310) || (abs(part.id()) == 3122) ){
					fastjet::PseudoJet pj(part.px(), part.py(), part.pz(), part.e());
					double dR = pj.delta_R(jet);
					double z = pj.px() * jet.px() + pj.py() * jet.py() + pj.pz() * jet.pz();
					z /= (jet.px() * jet.px() + jet.py() * jet.py() + jet.pz() * jet.pz());

					hz->Fill(jet.pt(), pj.pt(), z);
					hdRV0in->Fill(jet.pt(), pj.pt(), dR);
				}
			}
		} // jet loop
	} //End event loop

	// hjetpt->Write(hjetpt->GetName(), TObject::kOverwrite);
	// hv0pt->Write(hv0pt->GetName(), TObject::kOverwrite);
	// hNdauV0in->Write(hNdauV0in->GetName(), TObject::kOverwrite);
	// hDeltaPtV0in->Write(hDeltaPtV0in->GetName(), TObject::kOverwrite);
	hdRV0in->Write(hdRV0in->GetName(), TObject::kOverwrite);
	hz->Write(hz->GetName(), TObject::kOverwrite);
	outFile->Write();
	cout << "Histos written to file " << outFile->GetName() << endl;
	outFile->Close();
}
// */
