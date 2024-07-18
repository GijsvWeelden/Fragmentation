
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

double MK(double p)
{ return 0.00281882007 + 0.00114057004 * p + 0.00172138005 * TMath::Exp(-1. * p / 0.500262022); }
double ML(double p)
{ return 0.00117517996 + 0.000124098995 * p + 0.00547936978 * TMath::Exp(-1. * p / 0.308008999); }

void plotMassWindow()
{
  double pMin = 0.; double pMax = 55.; int pBins = 110;

  TH1D* mK0S = new TH1D("mK0S", "K^{0}_{S} mass acceptance;|#bf{#it{p}}| (GeV/#it{c}); #Delta M (MeV/#it{c}^2)", pBins, pMin, pMax);
  setStyle(mK0S, kBlue);
  TH1D* mL0 = new TH1D("mL0", "#Lambda mass acceptance;|#bf{#it{p}}| (GeV/#it{c}); #Delta M (MeV/#it{c}^2)", pBins, pMin, pMax);
  setStyle(mL0, kRed);

  TLegend* legend = CreateLegend(0.3, 0.6, 0.6, 0.8);
  legend->AddEntry(mK0S, "#Delta M(K^{0}_{S})", "l");
  legend->AddEntry(mL0, "#Delta M(#Lambda)", "l");

  for (int ip = 0; ip < pBins; ip++) {
    double p = mL0->GetBinCenter(ip);
    mK0S->Fill(p, 1e3*MK(p));
    mL0->Fill(p, 1e3*ML(p));
  }

  TCanvas* cM = new TCanvas("cM", "cM", 800, 600);
  TH1F* mf = DrawFrame(pMin, 50., 0, 60., "|#bf{#it{p}}| (GeV/#it{c})", "#Delta M(#Lambda, K^{0}_{S}) (MeV/#it{c}^{2})");

  mf->Draw();
  mK0S->Draw("hist L same");
  mL0->Draw("hist L same");
  legend->Draw("same");
  cM->SaveAs("MassWindow.pdf");
}
