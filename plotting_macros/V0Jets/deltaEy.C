
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

const double MassPi = 0.139570;
const double MassProton = 0.938272;
const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

double energy(double m, double p)
{ return sqrt(m * m + p * p); }

// double rapidity(double m, double p, double pz)
// { return 0.5 * TMath::Log((energy(m, p) + pz) / (energy(m, p) - pz)); }
double rapidity(double m, double p, double eta)
{ return 0.5 * TMath::Log((energy(m, p) + p * TMath::TanH(eta)) / (energy(m, p) - p * TMath::TanH(eta))); }

void plotEnergyAndRapidity()
{
  double pMin = 0.; double pMax = 50.; int pBins = 100;
  double etaMin = -1.; double etaMax = 1.; int etaBins = 100;

  TH1D* eK0S = new TH1D("eK0S", "K^{0}_{S} energy;|#bf{#it{p}}| (GeV/#it{c});E (GeV)", pBins, pMin, pMax);
  TH1D* eL0 = new TH1D("eL0", "#Lambda energy;|#bf{#it{p}}| (GeV/#it{c});E (GeV)", pBins, pMin, pMax);
  TH1D* dE;
  TH2D* yK0S = new TH2D("yK0S", "K^{0}_{S} rapidity;|#bf{#it{p}}| (GeV/#it{c});#eta;y", pBins, pMin, pMax, etaBins, etaMin, etaMax);
  TH2D* yL0 = new TH2D("yL0", "#Lambda rapidity;|#bf{#it{p}}| (GeV/#it{c});#eta;y", pBins, pMin, pMax, etaBins, etaMin, etaMax);
  TH2D* dy;

  for (int ip = 0; ip < pBins; ip++) {
    double p = eL0->GetBinCenter(ip);
    eK0S->Fill(p, energy(MassK0S, p));
    eL0->Fill(p, energy(MassLambda0, p));

    for (int ieta = 0; ieta < etaBins; ieta++) {
      double eta = yL0->GetYaxis()->GetBinCenter(ieta);
      yK0S->Fill(p, eta, rapidity(MassK0S, p, eta));
      yL0->Fill(p, eta, rapidity(MassLambda0, p, eta));
    }
  }

  TCanvas* cE = new TCanvas("cE", "cE", 800, 600);
  TH1F* Eframe = DrawFrame(0., 10., 0, 0.6, "|#bf{#it{p}}| (GeV/#it{c})", "E(#Lambda) - E(K^{0}_{S}) (GeV)");
  dE = (TH1D*) eL0->Clone("dE");
  dE->Add(eK0S, -1);
  dE->SetLineWidth(3);

  Eframe->Draw();
  dE->Draw("hist L same");
  cE->SaveAs("dE.pdf");

  TCanvas* cy = new TCanvas("cy", "cy", 800, 600);
  TH1F* yframe = DrawFrame(0., 10., -1., 1., "|#bf{#it{p}}| (GeV/#it{c})", "#eta");
  yframe->SetTitle("#it{y}(#Lambda) - #it{y}(K^{0}_{S})");
  dy = (TH2D*) yL0->Clone("dy");
  dy->Add(yK0S, -1);
  dy->SetZTitle("#it{y}(#Lambda) - #it{y}(K^{0}_{S})");

  yframe->Draw();
  dy->Draw("colz, same");
  cy->SaveAs("dy.pdf");
}

void plotRapidity(bool doStrange = true)
{
  double pMin = 0.; double pMax = 50.; int pBins = 100;
  double etaMin = -1.; double etaMax = 1.; int etaBins = 100;
  TH2D* yK0S = new TH2D("yK0S", "K^{0}_{S} rapidity;|#bf{#it{p}}| (GeV/#it{c});#eta;#it{y}", pBins, pMin, pMax, etaBins, etaMin, etaMax);
  TH2D* yL0 = new TH2D("yL0", "#Lambda rapidity;|#bf{#it{p}}| (GeV/#it{c});#eta;#it{y}", pBins, pMin, pMax, etaBins, etaMin, etaMax);
  TH2D* dy;

  for (int ip = 0; ip < pBins; ip++) {
    double p = yL0->GetXaxis()->GetBinCenter(ip);
    for (int ieta = 0; ieta < etaBins; ieta++) {
      double eta = yL0->GetYaxis()->GetBinCenter(ieta);
      yK0S->Fill(p, eta, rapidity(doStrange ? MassK0S : MassPi, p, eta));
      yL0->Fill(p, eta, rapidity(doStrange ? MassLambda0 : MassProton, p, eta));
      // yK0S->Fill(p, eta, rapidity(MassK0S, p, eta));
      // yL0->Fill(p, eta, rapidity(MassLambda0, p, eta));
    }
  }
  dy = (TH2D*) yL0->Clone("dy");
  dy->Add(yK0S, -1);

  vector<TH1D*> histVector;
  TLegend* leg = CreateLegend(0.5, 0.8, 0.2, 0.55);
  for (int i = 0; i < 8; i++) {
    TH2D* dyCopy = (TH2D*) dy->Clone(TString::Format("dyCopy_%d", i).Data());
    array<int, 2> etaProj = getProjectionBins(dyCopy->GetYaxis(), 0.1 * i, 0.1 * (i + 1));
    TH1D* dyProj = dyCopy->ProjectionX(TString::Format("dyProj_%d", i).Data(), etaProj[0], etaProj[1]);
    setStyle(dyProj, i);
    leg->AddEntry(dyProj, TString::Format("%.1f < #eta < %.1f", 0.1 * i, 0.1 * (i + 1)), "l");
    histVector.push_back(dyProj);
  }

  TCanvas* cdy = new TCanvas("cdy", "cdy", 800, 600);
  TH1F* dyframe = DrawFrame(0., 10., -1., 0.1, "|#bf{#it{p}}| (GeV)", "#it{y}(#Lambda) - #it{y}(K^{0}_{S})");
  if (!doStrange) { dyframe->SetYTitle("#it{y}(p) - #it{y}(#pi^{#pm})"); }
  dyframe->Draw();
  leg->Draw("same");
  for (auto dyProj : histVector) {
    dyProj->Draw("hist L same");
  }
  string saveName = "dy1D";
  saveName = TString::Format("%s_%s", saveName.c_str(), doStrange ? "LK" : "ppi").Data();
  saveName = TString::Format("%s.pdf", saveName.c_str());
  cdy->SaveAs(saveName.c_str());
}