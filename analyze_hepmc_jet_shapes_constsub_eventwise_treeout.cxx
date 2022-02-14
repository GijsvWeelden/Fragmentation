/*
 *  Program to analyse hepmc files and calculate jet shape variables.
 *     Output is stored a ROOT Tree (jetprops below)
 *
 *  Usage: analyze_hepmc_jet_shapes_constsub_eventwise_treeout [--chargedjets|--fulljets] [--nobkg] <infile> <outfilebase>
 *
 *  Argumemts/switches:
 *    --nobkg: do not subtract background (for pp and JEWEL without recoil)
 *    --fulljets|--chargedjets: use all particles, or only charged particles for jet finding.
 *  Author: Marco van Leeuwen, Nikhef
 *  Edited: Gijs van Weelden, Nikhef
 *
 */

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

#include "fastjet/Selector.hh" //.......... Background Sutraction event by event
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"//.......... Background Sutraction event by event
//#include "include/tools/Subtractor.hh"

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
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "THn.h"

#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
#include "getopt.h"

// defaults, can be set with arguments --chargedjets --fulljets --nobkg --bkg
int do_bkg = 0; // 0: no subtraction; 1: only jet energy; 2: energy and shape
int charged_jets = 0;

static const int debug = 0;
//static const int charged_constituents = 1; // This does not work yet for jet shapes!
static const float ptcut = 0.0; // GeV
static const float min_jet_pt = 10; // GeV; for shape histos

int is_stable(const HepMC::GenParticle *part) {
  // copied from AliStack::IsStable()
  int pdg = abs(part->pdg_id());
  if(pdg>1000000000)return kTRUE;

  const Int_t kNstable = 18;
  Int_t i;


  Int_t pdgStable[kNstable] = {
			       kGamma,             // Photon
			       kElectron,          // Electron
			       kMuonPlus,          // Muon
			       kPiPlus,            // Pion
			       kKPlus,             // Kaon
			       kK0Short,           // K0s
			       kK0Long,            // K0l
			       kProton,            // Proton
			       kNeutron,           // Neutron
			       kLambda0,           // Lambda_0
			       kSigmaMinus,        // Sigma Minus
			       kSigmaPlus,         // Sigma Plus
			       3312,               // Xsi Minus
			       3322,               // Xsi
			       3334,               // Omega
			       kNuE,               // Electron Neutrino
			       kNuMu,              // Muon Neutrino
			       kNuTau              // Tau Neutrino
  };

  Bool_t isStable = kFALSE;
  for (i = 0; i < kNstable; i++) {
    if (pdg == abs(pdgStable[i])) {
      isStable = kTRUE;
      break;
    }
  }

  return isStable;
}

int is_stable_charged(const HepMC::GenParticle *part) {
  // copied from AliStack::IsStable()
  int pdg = abs(part->pdg_id());
  if(pdg>1000000000)return kTRUE;

  const Int_t kNstableCharged = 9;
  Int_t i;


  Int_t pdgStableCharged[kNstableCharged] = {
					     kElectron,          // Electron
					     kMuonPlus,          // Muon
					     kPiPlus,            // Pion
					     kKPlus,             // Kaon
					     kProton,            // Proton
					     kSigmaMinus,        // Sigma Minus
					     kSigmaPlus,         // Sigma Plus
					     3312,               // Xsi Minus
					     3334                // Omega
  };

  Bool_t isStable = kFALSE;
  for (i = 0; i < kNstableCharged; i++) {
    if (pdg == abs(pdgStableCharged[i])) {
      isStable = kTRUE;
      break;
    }
  }

  return isStable;
}

int is_charged(const HepMC::GenParticle *part) {
  int abs_kf = abs(part->pdg_id());

  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13) // pi+, K+, p, e-, mu-
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16) // gamma, pi0, K0L, n, K0, nu-e, nu-mu, nu-tau
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

float dphi(float phi1, float phi2) {
  float dphi=phi1-phi2;
  float pi = 3.14159;
  if (dphi < -pi)
    dphi+=2*pi;
  if (dphi > pi)
    dphi-=2*pi;
  return dphi;
}

// Jet shapes from AliJetShape.h
float angularity(const fastjet::PseudoJet &jet) {
  if (!jet.has_constituents())
    return 0;
  Double_t den=0.;
  Double_t num = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    Double_t dphi = constits[ic].phi()-jet.phi();
    if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
    if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t dr2 = (constits[ic].eta()-jet.eta())*(constits[ic].eta()-jet.eta()) + dphi*dphi;
    Double_t dr = TMath::Sqrt(dr2);
    num=num+constits[ic].perp()*dr;
    den=den+constits[ic].perp();
  }
  return num/den;
}

void getmassangularities(const fastjet::PseudoJet &jet, float &rm, float &r2m, float &zs, float &z2m, float &rz, float &r2z) {
  if (!jet.has_constituents())
    return;
  rm = 0;
  r2m = 0;
  zs = 0;
  z2m = 0;
  rz = 0;
  r2z = 0;
  if (!jet.has_constituents())
    return;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    Double_t dphi = constits[ic].phi()-jet.phi();
    if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
    if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t dr2 = (constits[ic].eta()-jet.eta())*(constits[ic].eta()-jet.eta()) + dphi*dphi;
    Double_t dr = TMath::Sqrt(dr2);
    Double_t zfrag = constits[ic].perp()/jet.perp();
    rm += dr;
    r2m += dr*dr;
    zs += zfrag;
    z2m += zfrag*zfrag;
    rz += dr*zfrag;
    r2z += dr*dr*zfrag;
  }
  rm /= constits.size();
  r2m /= constits.size();
  z2m /= constits.size();
}

float pTD(const fastjet::PseudoJet &jet) {
  if (!jet.has_constituents())
    return 0;
  Double_t den=0;
  Double_t num = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    num=num+constits[ic].perp()*constits[ic].perp();
    den=den+constits[ic].perp();
  }
  return TMath::Sqrt(num)/den;
}

float circularity(const fastjet::PseudoJet &jet) {
  if (!jet.has_constituents())
    return 0;
  Double_t mxx    = 0.;
  Double_t myy    = 0.;
  Double_t mxy    = 0.;
  int  nc     = 0;
  Double_t sump2  = 0.;
  Double_t pxjet=jet.px();
  Double_t pyjet=jet.py();
  Double_t pzjet=jet.pz();

  //2 general normalized vectors perpendicular to the jet
  TVector3  ppJ1(pxjet, pyjet, pzjet);
  TVector3  ppJ3(- pxjet* pzjet, - pyjet * pzjet, pxjet * pxjet + pyjet * pyjet);
  ppJ3.SetMag(1.);
  TVector3  ppJ2(-pyjet, pxjet, 0);
  ppJ2.SetMag(1.);

  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    TVector3 pp(constits[ic].px(), constits[ic].py(), constits[ic].pz());
    //local frame
    TVector3 pLong = pp.Dot(ppJ1) / ppJ1.Mag2() * ppJ1;
    TVector3 pPerp = pp - pLong;
    //projection onto the two perpendicular vectors defined above
    Float_t ppjX = pPerp.Dot(ppJ2);
    Float_t ppjY = pPerp.Dot(ppJ3);
    Float_t ppjT = TMath::Sqrt(ppjX * ppjX + ppjY * ppjY);
    if(ppjT<=0) return 0;
    mxx += (ppjX * ppjX / ppjT);
    myy += (ppjY * ppjY / ppjT);
    mxy += (ppjX * ppjY / ppjT);
    nc++;
    sump2 += ppjT;
  }
  if(nc<2) return 0;
  if(sump2==0) return 0;
  // Sphericity Matrix
  Double_t ele[4] = {mxx / sump2, mxy / sump2, mxy / sump2, myy / sump2};
  TMatrixDSym m0(2,ele);

  // Find eigenvectors
  TMatrixDSymEigen m(m0);
  TVectorD eval(2);
  TMatrixD evecm = m.GetEigenVectors();
  eval  = m.GetEigenValues();
  // Largest eigenvector
  int jev = 0;
  if (eval[0] < eval[1]) jev = 1;
  TVectorD evec0(2);
  // Principle axis
  evec0 = TMatrixDColumn(evecm, jev);
  Double_t compx=evec0[0];
  Double_t compy=evec0[1];
  TVector2 evec(compx, compy);
  Double_t circ=0;
  if(jev==1) circ=2*eval[0];
  if(jev==0) circ=2*eval[1];

  return circ;
}

float sigma2(const fastjet::PseudoJet &jet) {
  if (!jet.has_constituents())
    return 0;
  Double_t mxx    = 0.;
  Double_t myy    = 0.;
  Double_t mxy    = 0.;
  int  nc     = 0;
  Double_t sump2  = 0.;

  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    Double_t ppt=constits[ic].perp();
    Double_t dphi = constits[ic].phi()-jet.phi();
    if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
    if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t deta = constits[ic].eta()-jet.eta();
    mxx += ppt*ppt*deta*deta;
    myy += ppt*ppt*dphi*dphi;
    mxy -= ppt*ppt*deta*dphi;
    nc++;
    sump2 += ppt*ppt;
  }
  if(nc<2) return 0;
  if(sump2==0) return 0;
  // Sphericity Matrix
  Double_t ele[4] = {mxx , mxy, mxy, myy };
  TMatrixDSym m0(2,ele);

  // Find eigenvectors
  TMatrixDSymEigen m(m0);
  TVectorD eval(2);
  TMatrixD evecm = m.GetEigenVectors();
  eval  = m.GetEigenValues();
  // Largest eigenvector
  int jev = 0;
  if (eval[0] < eval[1]) jev = 1;
  TVectorD evec0(2);
  // Principle axis
  evec0 = TMatrixDColumn(evecm, jev);
  Double_t compx=evec0[0];
  Double_t compy=evec0[1];
  TVector2 evec(compx, compy);
  Double_t sigma2=0;
  if(jev==1) sigma2=TMath::Sqrt(TMath::Abs(eval[0])/sump2);
  if(jev==0) sigma2=TMath::Sqrt(TMath::Abs(eval[1])/sump2);
  return sigma2;
}

float LeSub(const fastjet::PseudoJet &jet) {
  if (!jet.has_constituents())
    return 0;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  std::vector<fastjet::PseudoJet> sortedconstits=sorted_by_pt(constits);
  if(sortedconstits.size()<2) return 0;
  Double_t num=TMath::Abs(sortedconstits[0].perp()-sortedconstits[1].perp());
  return num;
}

int main(int argc, char **argv) {
  //
  // Takes two arguments: infile (HEPMC format) outfile (base name, ROOT format)
  // additional options: --nobkg --chargedjet|--fulljet
  //

  if (argc < 2) {
    cerr << "Need two arguments: infile outfile" << endl << "infile is HEPMC ascii format; outfile will be root format" << endl;
    cerr << "further option arguments: [--chargedjets|--fulljets] [--nobkg]" << endl;
    return 1;
  }

  const double jetR = 0.4;
  const double zcut = 0.1;
  const double beta = 0;

  const double max_eta_jet = 2.0;
  const double max_eta_track = 2.4; // max_eta_track = max_eta_jet + jetR

  int c;

  int nopt_parsed = 0;
  while (1) {

    static struct option long_options[] =
      {
       /* These options set a flag. */
       {"chargedjets", no_argument, &charged_jets, 1},
       {"fulljets",    no_argument, &charged_jets, 0},
       {"nobkg",       no_argument, &do_bkg, 0},
       {"bkg",         no_argument, &do_bkg, 2},
       /* it is also possible to have options that do not directly set a flag
	      * Not used for now */
       {0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 1;
    c = getopt_long (argc, argv, "",
		     long_options, &option_index);
    //cout << "c " << c << " option_index " << option_index << endl;
    /* Detect the end of the options. */
    if (c == -1)
      break;
    nopt_parsed++;
  }

  /* Print any remaining command line arguments (not options). */
  nopt_parsed++;
  cout << "option_index " << nopt_parsed << endl;
  if (nopt_parsed + 2 > argc){
    cerr << "Need two more arguments: infile outfile" << endl << "infile is HEPMC ascii format; outfile will be root format" << endl;
    return 1;
  }

  char *inname = argv[nopt_parsed];
  // specify an input file
  HepMC::IO_GenEvent ascii_in(inname,std::ios::in);

  // Make histos

  string outname(argv[nopt_parsed+1]);
  if (charged_jets)
    outname.append("_charged");
  else
    outname.append("_full");
  if (do_bkg == 0)
    outname.append("_nobkg");
  outname.append(".root");

  cout << "Input: " << inname << ", output " << outname << endl;

  const Int_t nPtJetBins = 40;  // Was 150 bins from 0 to 150
  const Float_t minPtJet = 0;
  const Float_t maxPtJet = 200;
  const Int_t useHarryZgBinning = 1;

  TFile fout(outname.c_str(),"RECREATE");

  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  hNEvent->Sumw2();
  TH1F *hPtLead = new TH1F("hPtLead","leading hadron pt;p_{T}",100,0,100);
  hPtLead->Sumw2();
  TH1D *hPtAllTrack = new TH1D("hPtAllTrack","Distribution of projected Pt;p_{T}",1000,-150,150);
  hPtAllTrack->Sumw2();
  TH2D *hPtPartEta = new TH2D("hPtPartEta","particle pt,eta;p_{T};#eta",100,0,100,200,-10,10);
  hPtPartEta->Sumw2();
  TH2D *hPtChTrackEta = new TH2D("hPtChTrackEta","charged particle pt,eta;p_{T};#eta",100,0,100,200,-10,10);
  hPtChTrackEta->Sumw2();

  Int_t ievt = 0, ijet = 0;
  Float_t evwt = 0, jet_eta=0, jet_phi=0, jet_pt=0, jet_dphi=0;
  Float_t zg=0, Rg=0, mass=0, mz2 = 0, mr = 0, mr2 = 0, rz = 0, r2z = 0, ptD = 0, t2t1 = 0, t3t2 = 0;
  Float_t t3dist[3];
  Int_t nconst=0, nSD=0;

  TTree *jetprops = new TTree("jetprops","Jet properties");
  jetprops->Branch("ievt",&ievt,"ievt/I");
  jetprops->Branch("ijet",&ijet,"ijet/I");
  jetprops->Branch("evwt",&evwt,"evwt/F");
  jetprops->Branch("pt",&jet_pt,"pt/F");
  jetprops->Branch("eta",&jet_eta,"eta/F");
  jetprops->Branch("phi",&jet_phi,"phi/F");
  jetprops->Branch("dphi",&jet_dphi,"dphi/F");
  jetprops->Branch("nconst",&nconst,"nconst/I");
  jetprops->Branch("zg",&zg,"zg/F");
  jetprops->Branch("Rg",&Rg,"Rg/F");
  jetprops->Branch("nSD",&nSD,"nSD/I");
  jetprops->Branch("mass",&mass,"mass/F");
  jetprops->Branch("mz2",&mz2,"mz2/F");
  jetprops->Branch("mr",&mr,"mr/F");
  jetprops->Branch("mr2",&mr2,"mr2/F");
  jetprops->Branch("rz",&rz,"rz/F");
  jetprops->Branch("r2z",&r2z,"r2z/F");
  jetprops->Branch("ptD",&ptD,"ptD/F");
  jetprops->Branch("t2t1",&t2t1,"t2t1/F");
  jetprops->Branch("t3t2",&t3t2,"t3t2/F");
  jetprops->Branch("t2dist",&t2dist,"t2dist/F");
  jetprops->Branch("t3dist",&t3dist,"t3dist[3]/F");

  // get the first event
  HepMC::GenEvent* evt = ascii_in.read_next_event();
  if (!evt)
    cerr << "Input file not found " << inname << endl;

  // loop until we run out of events
  while ( evt ) {
    // analyze the event
    if (debug)
      cout << "Event " << endl;

    evwt = evt->weights()[0]; // set event weight to fill in tree
    hNEvent->Fill(0.5,evwt); // count events
    // from example_UsingIterators.cc

    float pt_lead = -1;
    float phi_lead = -100;
    float max_eta_jet = 2.0;
    float max_eta_track = 2.4;
    float jetR = 0.4;

    int index = 0;
    std::vector <fastjet::PseudoJet> fjInputs;
    for (HepMC::GenEvent::particle_iterator pit = evt->particles_begin(); pit != evt->particles_end(); ++pit){
      const HepMC::GenParticle *p = *pit;
      if (!p->end_vertex() && p->status()==1 && (!charged_jets || is_charged(p))){ // final state charged particle
        hPtPartEta->Fill(p->momentum().perp(), p->momentum().eta(),evt->weights()[0]);
        if (is_charged(p))
          hPtChTrackEta->Fill(p->momentum().perp(), p->momentum().eta(),evt->weights()[0]);

        if ( fabs(p->momentum().eta()) < max_eta_track && p->momentum().perp() > ptcut){
          if (p->momentum().perp() > pt_lead){
            pt_lead = p->momentum().perp();
            phi_lead = p->momentum().phi();
          }
          double mom = sqrt(p->momentum().x()*p->momentum().x() +
                            p->momentum().y()*p->momentum().y() +
                            p->momentum().z()*p->momentum().z());
          //fastjet::PseudoJet jInp(p->momentum().x(),p->momentum().y(),p->momentum().z(),mom);
          fastjet::PseudoJet jInp(p->momentum().x(),p->momentum().y(),p->momentum().z(),p->momentum().e());  // need masses for E-scheme
          jInp.set_user_index(index);
          fjInputs.push_back(jInp);
          index++;
        }
      }
    }

    // Do jet finding
    // Need R =0.2 and R=0.4 later on...
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

    fastjet::JetMedianBackgroundEstimator bge;  //.......... Background Sutraction event by event
    fastjet::ClusterSequenceArea *clustSeqBG = 0;
    fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, jetR, recombScheme, strategy);

    vector<fastjet::PseudoJet> corrected_jets;
    fastjet::ClusterSequenceArea *clust_seq_corr = 0;

    if (do_bkg) {
      fastjet::Selector  BGSelector = fastjet::SelectorStrip(2*jetR);  //.......... Background Sutraction event by event
      clustSeqBG = new fastjet::ClusterSequenceArea(fjInputs, jetDefBG,areaDef); //............
      vector <fastjet::PseudoJet> BGJets = clustSeqBG->inclusive_jets();

      bge.set_selector(BGSelector);
      bge.set_jets(BGJets);

      // Check for zero area jets
      Int_t nZeroArea = 0;
      for (vector<fastjet::PseudoJet>::const_iterator it=BGJets.begin(); it!=BGJets.end(); it++) {
        if (it->area() == 0) {
          cout << "Found zero area jet pt " << it->perp() << " eta " << it->eta() << endl;
          nZeroArea++;
        }
      }
      if (nZeroArea != 0)
	      cout << "Found " << nZeroArea << " jets with zero area" << endl;

      fastjet::contrib::ConstituentSubtractor subtractor(&bge);
      // this sets the same background estimator to be used for deltaMass density, rho_m, as for pt density, rho:
      subtractor.set_common_bge_for_rho_and_rhom(true); // for massless input particles it does not make any difference (rho_m is always zero)
      //     cout << subtractor.description() << endl;
      subtractor.set_max_standardDeltaR(jetR);
      vector<fastjet::PseudoJet> corrected_event = subtractor.subtract_event(fjInputs,max_eta_track);
      if (debug > 0) {
        cout << "Event had " << fjInputs.size() << " tracks; " << corrected_event.size() << " after bkg sub" << endl;
      }
      clust_seq_corr = new fastjet::ClusterSequenceArea(corrected_event, jetDefCh, areaDefShape);
      corrected_jets = clust_seq_corr->inclusive_jets();
    }
    else
      corrected_jets = inclusiveJetsCh;

    // Sort corrected_jets by pt (vector<fastjet::PseudoJet>)
    vector <fastjet::PseudoJet> pt_sorted_jets = sorted_by_pt(corrected_jets);

    // modification 12/05/2014
    jet_eta= 0;
    jet_phi= 0;
    jet_pt= 0;
    float_t jet_phi_leading = 0;
    int n_jets_selected = 0; // Label for saved jet

    if (debug > 0)
      cout << corrected_jets.size() << " jets found" << endl;

    for (unsigned int iJet = 0; iJet < pt_sorted_jets.size(); iJet++){
      // Checks if jet in eta, phi range
      if (!range.is_in_range(pt_sorted_jets[iJet]))
        continue;

      jet_pt = pt_sorted_jets[iJet].perp();
      jet_eta = pt_sorted_jets[iJet].eta();
      float dphi_jh = dphi(pt_sorted_jets[iJet].phi(),phi_lead);
      jet_phi = pt_sorted_jets[iJet].phi();

      if (n_jets_selected == 0) // Jets are pt-sorted
        jet_phi_leading = jet_phi;

      jet_dphi = dphi(jet_phi,jet_phi_leading);

      if (jet_pt > min_jet_pt) {
        fastjet::PseudoJet &jet = pt_sorted_jets[iJet];
        float eta_jet = jet.eta();

        //float rm, rs, r2m, r2s, zs, rz, r2z;
        float zs;
        getmassangularities(jet, mr, mr2, zs, mz2, rz, r2z);
        mass = jet.m();
        ptD = pTD(jet);

        fastjet::contrib::SoftDrop sd(beta,zcut,jetR);
        fastjet::contrib::SoftDrop sd_kt(beta,zcut,jetR);
        fastjet::contrib::Recluster reclust_kt(fastjet::kt_algorithm, fastjet::JetDefinition::max_allowable_R);
        sd_kt.set_reclustering(true, &reclust_kt);
        fastjet::contrib::SoftDrop sd_akt(beta,zcut,jetR);
        fastjet::contrib::Recluster reclust_akt(fastjet::antikt_algorithm, fastjet::JetDefinition::max_allowable_R);
        sd_akt.set_reclustering(true, &reclust_akt);
        fastjet::contrib::ModifiedMassDropTagger mMDT(zcut);

        if (debug > 1)
          cout << "Softdrop " << jet.constituents().size() << " constituents " << endl;

        zg = -0.1;
        Rg = -0.1;
        nSD = 0;

        if (jet.has_associated_cluster_sequence()) {
          if (jet.has_pieces()) {
            fastjet::PseudoJet groomed_jet = sd(jet);
            if (groomed_jet.has_structure_of<fastjet::contrib::SoftDrop>()) {
              zg = groomed_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry();  // or mu() or delta_R()
              Rg = groomed_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R();  // or mu() or delta_R()
              // iterate to get nSD
              fastjet::PseudoJet j1 = groomed_jet, j2;
              while (j1.has_parents(j1,j2)) {
                if (j1.perp() < j2.perp()) std::swap(j1,j2);
                Double_t zg_cur = j2.perp()/(j1.perp()+j2.perp());
                if (zg_cur >= 0.1)
                  nSD++;
              }
            }
            else
              cout << "No groomed jet structure for jet with  pt " << jet.perp() <<  " E " << jet.E() << " eta " << jet.eta() << " :  " << jet.constituents().size() << " constituents; jet.has_structure(): " << jet.has_structure() << endl;
            if (debug > 2)
              cout << "z_g " << zg << " nSD " << nSD << endl;

          }
        }
        else cout << "No substructure stored with jet" << endl;

        const Float_t Beta = 1;
        fastjet::contrib::Nsubjettiness shape_tau1(1,fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(Beta,jetR)); // tau1
        fastjet::contrib::Nsubjettiness shape_tau2(2,fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(Beta,jetR)); // tau2
        fastjet::contrib::Nsubjettiness shape_tau3(3,fastjet::contrib::KT_Axes(), fastjet::contrib::NormalizedMeasure(Beta,jetR)); // tau3
        Double_t tau1 = shape_tau1(jet);
        Double_t tau2 = shape_tau2(jet);
        Double_t tau3 = shape_tau3(jet);

        // Distance between subjets
        std::vector<fastjet::PseudoJet>  SubJet_Axes=shape_tau2.currentAxes(); //Reclustered Axes result
        //SubJets=shape_tau2.currentSubjets(); //Resclustered Subjets Result
        Double_t R2subdist = -5;
        if (SubJet_Axes.size()>1){ //this should always be true
          fastjet::PseudoJet SubJet1_Axis=SubJet_Axes[0];
          Double_t SubJet1_Eta=SubJet1_Axis.pseudorapidity();
          Double_t SubJet1_Phi=SubJet1_Axis.phi();
          fastjet::PseudoJet SubJet2_Axis=SubJet_Axes[1];
          Double_t SubJet2_Eta=SubJet2_Axis.pseudorapidity();
          Double_t SubJet2_Phi=SubJet2_Axis.phi();
          Double_t DeltaPhi=SubJet1_Phi-SubJet2_Phi; //Delta Phi between the two axis of reclustering
          if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
          else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
          R2subdist = TMath::Sqrt(DeltaPhi*DeltaPhi + (SubJet2_Eta-SubJet1_Eta)*(SubJet2_Eta-SubJet1_Eta));
        }

        std::vector<fastjet::PseudoJet> SubJets_Axes = shape_tau3.currentAxes();
        Double_t R3subdist = {-5, -5, -5};
        if (SubJets_Axes.size() > 1){
          fastjet::PseudoJet SubJet1_Axis = SubJet_Axes[0];
          Double_t SubJet1_Eta=SubJet1_Axis.pseudorapidity();
          Double_t SubJet1_Phi=SubJet1_Axis.phi();
          fastjet::PseudoJet SubJet2_Axis = SubJet_Axes[1];
          Double_t SubJet2_Eta=SubJet2_Axis.pseudorapidity();
          Double_t SubJet2_Phi=SubJet2_Axis.phi();
          fastjet::PseudoJet SubJet3_Axis = SubJet_Axes[2];
          Double_t SubJet3_Eta=SubJet3_Axis.pseudorapidity();
          Double_t SubJet3_Phi=SubJet3_Axis.phi();

          Double_t DeltaPhi12 = SubJet1_Phi - SubJet2_Phi;
          if(DeltaPhi12 < -1*TMath::Pi()) DeltaPhi12 += (2*TMath::Pi());
          else if (DeltaPhi12 > TMath::Pi()) DeltaPhi12 -= (2*TMath::Pi());
          R3subdist[0] = TMath::Sqrt(DeltaPhi12*DeltaPhi12 + (SubJet2_Eta-SubJet1_Eta)*(SubJet2_Eta-SubJet1_Eta));

          Double_t DeltaPhi13 = SubJet1_Phi - SubJet3_Phi;
          if(DeltaPhi13 < -1*TMath::Pi()) DeltaPhi13 += (2*TMath::Pi());
          else if (DeltaPhi13 > TMath::Pi()) DeltaPhi13 -= (2*TMath::Pi());
          R3subdist[1] = TMath::Sqrt(DeltaPhi13*DeltaPhi13 + (SubJet2_Eta-SubJet1_Eta)*(SubJet2_Eta-SubJet1_Eta));

          Double_t DeltaPhi23 = SubJet2_Phi - SubJet3_Phi;
          if(DeltaPhi23 < -1*TMath::Pi()) DeltaPhi23 += (2*TMath::Pi());
          else if (DeltaPhi23 > TMath::Pi()) DeltaPhi23 -= (2*TMath::Pi());
          R3subdist[2] = TMath::Sqrt(DeltaPhi23*DeltaPhi23 + (SubJet2_Eta-SubJet1_Eta)*(SubJet2_Eta-SubJet1_Eta));
        }

        if (tau1 > 0) t2t1 = tau2 / tau1;
        else t2t1 = -1.;
        if (tau2 > 0){
          t3t2 = tau3 / tau2;
          t2dist = Rsubdist;
        }
        else t3t2 = -1.;
        if (tau3 > 0){
          for (int i = 0; i < t3dist.size(); ++i){
            t3dist[i] = R3subdist[i];
          }
        }


        ijet = n_jets_selected;
        n_jets_selected++;
        nconst = jet.constituents().size();
        jetprops->Fill();
      }
    }
    if (clustSeqBG) {
      delete clustSeqBG;
      clustSeqBG = 0;
    }
    if (clust_seq_corr) {
      delete clust_seq_corr;
      clust_seq_corr = 0;
    }

    // delete the created event from memory
    delete evt;
    // read the next event
    ascii_in >> evt;
    ievt++;
  }

  fout.Write();
  cout << "Simulation finished";
  fout.Close();
  return 0;
}


