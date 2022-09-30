#include <iostream>
#include <cmath>
#include <string>
#include <unistd.h>
#include "Pythia8/Pythia.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TH3F.h"

#define nEvents 5000000

using namespace Pythia8;

int main(int argc, char**argv)
{

  string foutname("PythiaDiCharmResult_pth20.root");

  if (argc > 1) foutname = argv[1];
  cout << "Output file name: " << foutname << endl;
  //PYTHIA SETTINGS

  TString name;

  int mecorr=1;

  Float_t ptHatMin=-1; //20; //-1; //80;
  Float_t ptHatMax=-1; //100;

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readString("Beams:idA = 2212"); //beam 1 proton
  pythia.readString("Beams:idB = 2212"); //beam 2 proton
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("Tune:pp = 5");  //tune 1-13    5=defaulr TUNE4C,  6=Tune 4Cx, 7=ATLAS MB Tune A2-CTEQ6L1

  pythia.readString("Random:setSeed = on");
  int seed = getpid() + clock();
  Char_t tmpstr[256];
  sprintf(tmpstr,"Random:seed = %d",seed);
  pythia.readString(tmpstr);
  //pythia.readString("Random:seed = 0");

  pythia.readString("SoftQCD:all = on");
  //pythia.readString("HardQCD:all = on");
  //pythia.readString("HardQCD:hardccbar = on");
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
  pythia.readString("421:onIfAll = 321 211"); // D0 -> K pi

  //ME corrections
  //use of matrix corrections where available
  if(mecorr==0){ 
    pythia.readString("TimeShower:MECorrections=off");
  }
  pythia.init();

  const Float_t etaMax = 4; // for delta eta and delta phi plots
  const Float_t dphiNearSide = 1;
  // Output histograms
  Int_t nPtBins = 100;
  Float_t ptMax = 50;
  TFile* outFile = new TFile(foutname.c_str(),"RECREATE");
  TH1F *hNEvent = new TH1F("hNEvent","n Event",1 ,0 ,2);
  TProfile *hXSec = new TProfile("hXSec","cross section",1 ,0 ,2);
  TH2F *hEtaPt = new TH2F("hEtaPt","Pt vs Eta for all particles;#eta;p_{T} (GeV/c)",40,-5,5,50,0,25);
  TH2F *hEtaPtCharm = new TH2F("hEtaPtCharm","Pt vs Eta for D^{0}, #bar{D}^{0};#eta;p_{T} (GeV/c)",40,-5,5,nPtBins, 0, ptMax);
  TH3F *hPt1Pt2EtaDiCharm = new TH3F("hPt1Pt2EtaCharm","Pt1 vs Pt2 vs max eta for di-charm;p_{T,1} (GeV/c);p_{T,2} (GeV/c);#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2EtaDiCharmNS = new TH3F("hPt1Pt2EtaCharmNS","Pt1 vs Pt2 vs max eta for di-charm (near side);p_{T,1} (GeV/c);p_{T,2} (GeV/c);#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2EtaDiCharmAS = new TH3F("hPt1Pt2EtaCharmAS","Pt1 vs Pt2 vs max eta for di-charm (away side);p_{T,1} (GeV/c);p_{T,2} (GeV/c);#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2EtaDiCharmDecay = new TH3F("hPt1Pt2EtaCharmDecay","Pt1 vs Pt2 vs max eta for di-charm decay daughters;p_{T,1} (GeV/c);p_{T,2} (GeV/c);#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2EtaDiCharmDecayNS = new TH3F("hPt1Pt2EtaCharmDecayNS","Pt1 vs Pt2 vs max eta for di-charm decay daughters (near side);p_{T,1} (GeV/c);p_{T,2} (GeV/c);#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2EtaDiCharmDecayAS = new TH3F("hPt1Pt2EtaCharmDecayAS","Pt1 vs Pt2 vs max eta for di-charm decay daughters (away side);p_{T,1} (GeV/c);p_{T,2} (GeV/c);#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);

  TH3F *hPt1Pt2DEtaDiCharm = new TH3F("hPt1Pt2DEtaCharm","Pt1 vs Pt2 vs delta eta for di-charm;p_{T,1} (GeV/c);p_{T,2} (GeV/c);#Delta#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2DEtaDiCharmNS = new TH3F("hPt1Pt2DEtaCharmNS","Pt1 vs Pt2 vs delta eta for di-charm (near side);p_{T,1} (GeV/c);p_{T,2} (GeV/c);#Delta#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2DEtaDiCharmAS = new TH3F("hPt1Pt2DEtaCharmAS","Pt1 vs Pt2 vs delta eta for di-charm (away side);p_{T,1} (GeV/c);p_{T,2} (GeV/c);#Delta#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2DEtaDiCharmDecay = new TH3F("hPt1Pt2DEtaCharmDecay","Pt1 vs Pt2 vs delta eta for di-charm decay products;p_{T,1} (GeV/c);p_{T,2} (GeV/c);#Delta#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2DEtaDiCharmDecayNS = new TH3F("hPt1Pt2DEtaCharmDecayNS","Pt1 vs Pt2 vs delta eta for di-charm decay products (near side);p_{T,1} (GeV/c);p_{T,2} (GeV/c);#Delta#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPt1Pt2DEtaDiCharmDecayAS = new TH3F("hPt1Pt2DEtaCharmDecayAS","Pt1 vs Pt2 vs delta eta for di-charm decay products (away side);p_{T,1} (GeV/c);p_{T,2} (GeV/c);#Delta#eta",nPtBins, 0, ptMax,nPtBins, 0, ptMax,40,-5,5);

  TH3F *hPt1Pt2DPhiDiCharm = new TH3F("hPt1Pt2DPhiCharm","Pt1 vs Pt2 vs delta phi for di-charm;p_{T,1} (GeV/c);p_{T,2} (GeV/c);#Delta#phi",nPtBins, 0, ptMax,nPtBins, 0, ptMax,48,-0.5*TMath::Pi(),1.5*TMath::Pi());

  TH3F *hPtLeadEtaDEtaCharm = new TH3F("hPtLeadEtaDEtaCharm","delta eta, eta vs Pt leading for di-charm;p_{T,lead} (GeV/c);#eta;#Delta#eta",nPtBins, 0, ptMax,40,-5,5,40,-5,5);
  TH3F *hPtLeadDEtaDPhiCharm = new TH3F("hPtLeadDEtaDPhiCharm","delta phi, delta eta vs Pt leading for di-charm;p_{T,lead} (GeV/c);#Delta#eta;#Delta#phi",nPtBins, 0, ptMax,40,-5,5,48,-0.5*TMath::Pi(),1.5*TMath::Pi());
  TH3F *hPtLeadDPtDPhiCharm = new TH3F("hPtLeadDPtDPhiCharm","delta phi, delta eta vs Pt leading for di-charm;p_{T,lead} (GeV/c);#Deltap_{T} (GeV/c);#Delta#phi",nPtBins, 0, ptMax, nPtBins, 0, ptMax,48,-0.5*TMath::Pi(),1.5*TMath::Pi());

  TH3F *hPt1Pt2EtaQuarks = new TH3F("hPt1Pt2EtaQuarks","Pt1 vs Pt2 vs max eta for ccbar quarks;p_{T,1} (GeV/c);p_{T,2} (GeV/c);#eta",nPtBins, 0, ptMax, nPtBins, 0, ptMax,40,-5,5);
  TH3F *hPtLeadDPtDPhiQuarks = new TH3F("hPtLeadDPtDPhiQuarks","delta phi, delta eta vs Pt leading for c-cbar quarks;p_{T,lead} (GeV/c);#Deltap_{T} (GeV/c);#Delta#phi",nPtBins, 0, ptMax, nPtBins, 0, ptMax,48,-0.5*TMath::Pi(),1.5*TMath::Pi());

  Int_t nCharm = 0;
  Int_t nDiCharm = 0;
  //Begin event loop
  for(int iEvent = 0; iEvent < nEvents; iEvent++)
    {
      double fourvec[4];

      if(!pythia.next()) continue;
      
      hNEvent->Fill(1);

      Double_t ptSumPythia = 0;
      Int_t nPartPythia = 0;
      int nPart = pythia.event.size();
      Particle const *charmQ = 0;
      Particle const *antiCharmQ = 0;
      Particle const *charmQI = 0;
      Particle const *antiCharmQI = 0;
      Particle const *charmPart = 0;
      Particle const *antiCharmPart = 0;
      Float_t maxEtaDght = 0;
      for(int iPart = 0; iPart < nPart; iPart++) {		  
	const Particle &part = pythia.event[iPart];
	if(part.isFinal()) {
	  hEtaPt->Fill(part.eta(),part.pT());
	  nPartPythia++;
	}
	else if (abs(part.id()) == 421) { // Charm is decayed in Pythia...
	    nCharm++;
	    hEtaPtCharm->Fill(part.eta(), part.pT());
	    if (part.id() > 0) {
	      charmPart = &part;
	    }
	    else {
	      antiCharmPart = &part;
	    }
            Float_t eta_dght = pythia.event[part.daughter1()].eta();
            if (fabs(eta_dght) > fabs(maxEtaDght)) maxEtaDght = eta_dght;
            eta_dght = pythia.event[part.daughter2()].eta();
            if (fabs(eta_dght) > fabs(maxEtaDght)) maxEtaDght = eta_dght;
	  }
        else if (part.id()==4) {
	  if (abs(part.status()) < 30)
	    charmQI = &part;
	  else
	    charmQ = &part;
        }
        else if (part.id()==-4) {
	  if (abs(part.status()) < 30)
	    antiCharmQI = &part;
	  else
	    antiCharmQ = &part;
        }
      }
		
      if (charmPart && antiCharmPart) {
	nDiCharm++;
	Float_t eta2 = charmPart->eta();
	if (fabs(antiCharmPart->eta()) > fabs(eta2)) eta2 = antiCharmPart->eta(); 
	hPt1Pt2EtaDiCharm->Fill(charmPart->pT(), antiCharmPart->pT(), eta2); 
	hPt1Pt2EtaDiCharmDecay->Fill(charmPart->pT(), antiCharmPart->pT(), maxEtaDght); 
	Float_t dphi = antiCharmPart->phi() - charmPart->phi();
	if (dphi < 0.5*TMath::Pi())
	   dphi += 2*TMath::Pi();
	if (dphi > 1.5*TMath::Pi())
	   dphi -= 2*TMath::Pi();
        if (fabs(dphi) < dphiNearSide) {
          	hPt1Pt2EtaDiCharmNS->Fill(charmPart->pT(), antiCharmPart->pT(), eta2); 
        	hPt1Pt2EtaDiCharmDecayNS->Fill(charmPart->pT(), antiCharmPart->pT(), maxEtaDght); 
        }
        else {
        	hPt1Pt2EtaDiCharmAS->Fill(charmPart->pT(), antiCharmPart->pT(), eta2); 
        	hPt1Pt2EtaDiCharmDecayAS->Fill(charmPart->pT(), antiCharmPart->pT(), maxEtaDght); 
        }
	if (fabs(charmPart->eta()) < etaMax && 
	    fabs(antiCharmPart->eta()) < etaMax) {
	  hPt1Pt2DEtaDiCharm->Fill(charmPart->pT(), antiCharmPart->pT(), antiCharmPart->eta() - charmPart->eta());
	  hPt1Pt2DPhiDiCharm->Fill(charmPart->pT(), antiCharmPart->pT(), dphi);
          if (fabs(dphi) < dphiNearSide) {
	    hPt1Pt2DEtaDiCharmNS->Fill(charmPart->pT(), antiCharmPart->pT(), antiCharmPart->eta() - charmPart->eta());
          }
          else {
	    hPt1Pt2DEtaDiCharmAS->Fill(charmPart->pT(), antiCharmPart->pT(), antiCharmPart->eta() - charmPart->eta());
          }
          if (fabs(charmPart->id()) == 421 && fabs(antiCharmPart->id()) == 421) {
            Float_t max_deta = 0;
            Float_t deta_dght = pythia.event[antiCharmPart->daughter1()].eta()-pythia.event[charmPart->daughter1()].eta();
            if (fabs(deta_dght) > fabs(max_deta)) max_deta = deta_dght;
            deta_dght = pythia.event[antiCharmPart->daughter1()].eta()-pythia.event[charmPart->daughter2()].eta();
            if (fabs(deta_dght) > fabs(max_deta)) max_deta = deta_dght;
            deta_dght = pythia.event[antiCharmPart->daughter2()].eta()-pythia.event[charmPart->daughter2()].eta();
            if (fabs(deta_dght) > fabs(max_deta)) max_deta = deta_dght;
            deta_dght = pythia.event[antiCharmPart->daughter2()].eta()-pythia.event[charmPart->daughter1()].eta();
            if (fabs(deta_dght) > fabs(max_deta)) max_deta = deta_dght;
	    hPt1Pt2DEtaDiCharmDecay->Fill(charmPart->pT(), antiCharmPart->pT(), max_deta);
            if (fabs(dphi) < dphiNearSide) {
	       hPt1Pt2DEtaDiCharmDecayNS->Fill(charmPart->pT(), antiCharmPart->pT(), max_deta);
            }
            else {
	       hPt1Pt2DEtaDiCharmDecayAS->Fill(charmPart->pT(), antiCharmPart->pT(), max_deta);
            }
          }
          Float_t ptLead = charmPart->pT();
          Float_t etaLead = charmPart->eta();
          if (antiCharmPart->pT() > charmPart->pT()) {
	    ptLead = antiCharmPart->pT();
	    etaLead = antiCharmPart->eta();
          }
          hPtLeadEtaDEtaCharm->Fill(ptLead, etaLead, charmPart->eta()-antiCharmPart->eta());
          hPtLeadDEtaDPhiCharm->Fill(ptLead, etaLead, charmPart->eta()-antiCharmPart->eta(), dphi);
          hPtLeadDPtDPhiCharm->Fill(ptLead, etaLead, fabs(charmPart->pT()-antiCharmPart->pT()), dphi);
 
	}
      }
      if (charmQ && antiCharmQ) {
	/*
	pythia.readString("Stat:showPartonLevel = on");
	cout << "found ccbar: " << charmQ->id() << " " << antiCharmQ->id() << " status " << charmQ->status() <<" " << antiCharmQ->status() << " phi " << charmQ->phi() << " " << antiCharmQ->phi() << " pt " << charmQ->pT() << " " << antiCharmQ->pT() << endl;
	if (charmQI && antiCharmQI) {
	  cout << "initial ccbar: " << charmQI->id() << " " << antiCharmQI->id() << " status " << charmQI->status() <<" " << antiCharmQI->status() << " phi " << charmQI->phi() << " " << antiCharmQI->phi() << " pt " << charmQI->pT() << " " << antiCharmQI->pT() << endl;
	}
	else
	  cout << "Initial charm not found" << endl;
	cout << "parton 7 " << pythia.event[7].id() << " phi " << pythia.event[7].phi() << endl;
	cout << "parton 8 " << pythia.event[8].id() << " phi " << pythia.event[8].phi() << endl;
	*/
        Float_t eta2 = charmQ->eta();
	if (fabs(antiCharmQ->eta()) > fabs(eta2)) eta2 = antiCharmQ->eta(); 
	hPt1Pt2EtaQuarks->Fill(charmQ->pT(), antiCharmQ->pT(), eta2);
        if (fabs(eta2) < etaMax) {
          Float_t ptLead = charmQ->pT();
          if (antiCharmQ->pT() > ptLead) ptLead = antiCharmQ->pT();
          Float_t dphi = charmQ->phi()-antiCharmQ->phi();
	  if (dphi < 0.5*TMath::Pi())
	    dphi += 2*TMath::Pi();
	  if (dphi > 1.5*TMath::Pi())
	    dphi -= 2*TMath::Pi();

	  hPtLeadDPtDPhiQuarks->Fill(ptLead, fabs(charmQ->pT()-antiCharmQ->pT()), dphi);
        }
      }
      if ((iEvent%1000)==0) {
	cout << "Pythia event: " << nPartPythia << " particles" << endl;
	cout << "   total inclusive charm: " << nCharm << " dicharm "<< nDiCharm << endl;
      }
    }
  //End event loop

  hXSec->Fill(1,pythia.info.sigmaGen()); // This only gives the xsec for the main process?
  hXSec->SetBinError(1,pythia.info.sigmaErr()); // This only gives the xsec for the main process?

  outFile->Write();
  cout << "Histos written to file " << outFile->GetName() << endl;
  outFile->Close();
}

