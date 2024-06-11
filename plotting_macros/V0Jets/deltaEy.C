
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TString.h"
#include "TLegend.h"

#include "../histUtils.C"

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

double energy(double m, double p)
{ return sqrt(m * m + p * p); }

// double rapidity(double m, double p, double pz)
// { return 0.5 * TMath::Log((energy(m, p) + pz) / (energy(m, p) - pz)); }
double rapidity(double m, double p, double eta)
{ return 0.5 * TMath::Log((energy(m, p) + p * TMath::TanH(eta)) / (energy(m, p) - p * TMath::TanH(eta))); }

void plotEnergyAndEta()
{
  double pMin = 0.; double pMax = 50.; int pBins = 100;
  double etaMin = -1.; double etaMax = 1.; int etaBins = 100;

  TH1D* eK0S = new TH1D("eK0S", "K^{0}_{S} energy;E (GeV);Counts", pBins, pMin, pMax);
  TH1D* eL0 = new TH1D("eL0", "#Lambda energy;E (GeV);Counts", pBins, pMin, pMax);
  TH1D* dE;
  TH2D* etaK0S = new TH2D("etaK0S", "K^{0}_{S} rapidity;E (GeV);#eta;Counts", pBins, pMin, pMax, etaBins, etaMin, etaMax);
  TH2D* etaL0 = new TH2D("etaL0", "#Lambda rapidity;E (GeV);#eta;Counts", pBins, pMin, pMax, etaBins, etaMin, etaMax);
  TH2D* dEta;

  for (int ip = 0; ip < pBins; ip++) {
    double p = eL0->GetBinCenter(ip);
    eK0S->Fill(p, energy(MassK0S, p));
    eL0->Fill(p, energy(MassLambda0, p));

    for (int ieta = 0; ieta < etaBins; ieta++) {
      double eta = etaL0->GetYaxis()->GetBinCenter(ieta);
      etaK0S->Fill(p, eta, rapidity(MassK0S, p, eta));
      etaL0->Fill(p, eta, rapidity(MassLambda0, p, eta));
    }
  }

  dE = (TH1D*) eL0->Clone("dE");
  dE->SetYTitle("E(#Lambda) - E(K^{0}_{S}) (GeV)");
  dE->Add(eK0S, -1);
  dE->SetLineWidth(3);

  TCanvas* cE = new TCanvas("cE", "cE", 800, 600);
  dE->Draw("hist L");
  cE->SaveAs("dE.pdf");

  dEta = (TH2D*) etaL0->Clone("dEta");
  dEta->SetZTitle("#Delta#eta");
  dEta->Add(etaK0S, -1);

  TCanvas* cEta = new TCanvas("cEta", "cEta", 800, 600);
  dEta->Draw("colz");
  cEta->SaveAs("dEta.pdf");

  // Loop over p and eta
}