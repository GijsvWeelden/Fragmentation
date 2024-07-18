
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

// -------------------------------------------------------------------------------------------------
// Makes histogram to be used as a shaded region
template <typename T>
T* fillShadedRegion(T* data, int minBin, int maxBin)
{
  T* region = (T*)data->Clone("region");
  region->Reset();
  for (int i = minBin; i <= maxBin; i++) {
    region->SetBinContent(i, data->GetBinContent(i));
  }
  return region;
}

// -------------------------------------------------------------------------------------------------
// Initialise parameters and set bounds
void parsK0S(TF1* f, bool doubleGauss, int bkg)
{
  double gausPar[3] = {1., MassK0S, 0.01};
  double gausVar[6] = {0.4, 1., MassK0S - 5e-2, MassK0S + 5e-2, 1e-3, 0.02};
  double doubleGaussPar[3] = {0.1, MassK0S, 0.05};
  double doubleGaussVar[6] = {0.1, 0.6, MassK0S - 0.01, MassK0S + 0.01, 1e-3, 0.1};
  double bkgPar[3] = {0., 0., 0.};
  double bkgVar[6] = {-0.1, 0.1, -0.1, 0.1, -0.1, 0.1};

  f->SetParameter(0, gausPar[0]); f->SetParLimits(0, gausVar[0], gausVar[1]);
  f->SetParameter(1, gausPar[1]); f->SetParLimits(1, gausVar[2], gausVar[3]);
  f->SetParameter(2, gausPar[2]); f->SetParLimits(2, gausVar[4], gausVar[5]);
  if (doubleGauss) {
    f->SetParameter(3, doubleGaussPar[0]); f->SetParLimits(3, doubleGaussVar[0], doubleGaussVar[1]);
    f->SetParameter(4, doubleGaussPar[1]); f->SetParLimits(4, doubleGaussVar[2], doubleGaussVar[3]);
    f->SetParameter(5, doubleGaussPar[2]); f->SetParLimits(5, doubleGaussVar[4], doubleGaussVar[5]);
    f->SetParameter(6, bkgPar[0]);         f->SetParLimits(6, bkgVar[0], bkgVar[1]);
    f->SetParameter(7, bkgPar[1]);         f->SetParLimits(7, bkgVar[2], bkgVar[3]);
    if (2 >= bkg) {
      f->SetParameter(8, bkgPar[2]);       f->SetParLimits(8, bkgVar[4], bkgVar[5]);
    }
  }
  else {
    f->SetParameter(3, bkgPar[0]);   f->SetParLimits(3, bkgVar[0], bkgVar[1]);
    f->SetParameter(4, bkgPar[1]);   f->SetParLimits(4, bkgVar[2], bkgVar[3]);
    if (2 >= bkg) {
      f->SetParameter(5, bkgPar[2]); f->SetParLimits(5, bkgVar[4], bkgVar[5]);
    }
  }
}
void parsL0(TF1* f, bool doubleGauss, int bkg)
{
  double gausPar[3] = {1., MassLambda0, 0.01};
  double gausVar[6] = {0.4, 1., MassLambda0 - 0.01, MassLambda0 + 0.01, 1e-6, 0.02};
  double doubleGaussPar[3] = {0.1, MassLambda0, 0.05};
  double doubleGaussVar[6] = {1e-2, 0.6, MassLambda0 - 0.03, MassLambda0 + 0.03, 1e-3, 0.1};
  double bkgPar[3] = {0., 0., 0.};
  double bkgVar[6] = {-3., 3., -3., 3., -3., 3.};

  f->SetParameter(0, gausPar[0]); f->SetParLimits(0, gausVar[0], gausVar[1]);
  f->SetParameter(1, gausPar[1]); f->SetParLimits(1, gausVar[2], gausVar[3]);
  f->SetParameter(2, gausPar[2]); f->SetParLimits(2, gausVar[4], gausVar[5]);
  if (doubleGauss) {
    f->SetParameter(3, doubleGaussPar[0]); f->SetParLimits(3, doubleGaussVar[0], doubleGaussVar[1]);
    f->SetParameter(4, doubleGaussPar[1]); f->SetParLimits(4, doubleGaussVar[2], doubleGaussVar[3]);
    f->SetParameter(5, doubleGaussPar[2]); f->SetParLimits(5, doubleGaussVar[4], doubleGaussVar[5]);
    f->SetParameter(6, bkgPar[0]);         f->SetParLimits(6, bkgVar[0], bkgVar[1]);
    f->SetParameter(7, bkgPar[1]);         f->SetParLimits(7, bkgVar[2], bkgVar[3]);
    if (2 >= bkg) {
      f->SetParameter(8, bkgPar[2]);       f->SetParLimits(8, bkgVar[4], bkgVar[5]);
    }
  }
  else {
    f->SetParameter(3, bkgPar[0]);   f->SetParLimits(3, bkgVar[0], bkgVar[1]);
    f->SetParameter(4, bkgPar[1]);   f->SetParLimits(4, bkgVar[2], bkgVar[3]);
    if (2 >= bkg) {
      f->SetParameter(5, bkgPar[2]); f->SetParLimits(5, bkgVar[4], bkgVar[5]);
    }
  }
}
// Get fit function for signal and background
TF1* getSigBkgFit(string hadron, int bkg, bool doubleGauss)
{
  double fitRegion[2] = {1.08, 1.215};
  if ("K0S" == hadron) { fitRegion[0] = 0.44, fitRegion[1] = 0.555; }

  int arg = 3;
  string sFit = "gaus(0)";
  if (doubleGauss) {
    sFit += " + gaus(" + to_string(arg) + ")";
    arg += 3;
  }
  sFit += " + pol" + to_string(bkg) + "(" + to_string(arg) + ")";
  TF1* fit = new TF1("fit", sFit.c_str(), fitRegion[0], fitRegion[1]);

  if ("K0S" == hadron) {
    parsK0S(fit, doubleGauss, bkg);
  }
  else {
    parsL0(fit, doubleGauss, bkg);
  }
  return fit;
}
// Get shaded regions for signal and background
array<TH1*, 2> getSigBkgRegions(TH1* mass, double mean, double sigma, double nSigmaSignal, double nSigmaBkgMin, double nSigmaBkgMax)
{
  array<int, 2> sigBins      = getProjectionBins(mass->GetXaxis(), mean - nSigmaSignal * sigma, mean + nSigmaSignal * sigma);
  array<int, 2> leftBkgBins  = getProjectionBins(mass->GetXaxis(), mean - nSigmaBkgMax * sigma, mean - nSigmaBkgMin * sigma);
  array<int, 2> rightBkgBins = getProjectionBins(mass->GetXaxis(), mean + nSigmaBkgMin * sigma, mean + nSigmaBkgMax * sigma);
  // Prevent overlap of signal and background regions
  if (leftBkgBins[1] >= sigBins[0]) {
    int shift = leftBkgBins[1] - sigBins[0] + 1;
    leftBkgBins[0] -= shift;
    leftBkgBins[1] -= shift;
    cout << "Warning: Overlap of signal and background regions. Shifting left background region by " << shift << " bins." << endl;
  }
  if (rightBkgBins[0] <= sigBins[1]) {
    int shift = sigBins[1] - rightBkgBins[0] + 1;
    rightBkgBins[0] += shift;
    rightBkgBins[1] += shift;
    cout << "Warning: Overlap of signal and background regions. Shifting right background region by " << shift << " bins." << endl;
  }
  TH1D* sigRegion = (TH1D*)fillShadedRegion(mass, sigBins[0], sigBins[1]);
  TH1D* bkgRegion = (TH1D*)fillShadedRegion(mass, leftBkgBins[0], leftBkgBins[1]);
  bkgRegion->Add(fillShadedRegion(mass, rightBkgBins[0], rightBkgBins[1]));
  return {sigRegion, bkgRegion};
}

// -------------------------------------------------------------------------------------------------

// Extract purity for given hadron from histogram V0PtMass
// bkg = n: background is polynomial of order n
// doubleGauss = true: fit with 2 Gaussian
// flipGaussians = true: flip which Gaussian is treated as the signal peak
void v0Purity(string inName, string dataSet, string hadron, double ptmin = 0., double ptmax = 100., int bkg = 1, bool doubleGauss = true, bool flipGaussians = false, int rebinNumber = 1)
{
  int nSigmaSignal = 3, nSigmaBkgMin = 5, nSigmaBkgMax = 8;

  const int nDim = 4;
  const int ptAxis = 0;
  const int K0SmassAxis = 1;
  const int Lambda0massAxis = 2;
  const int AntiLambda0massAxis = 3;
  gStyle->SetNdivisions(505, "xy");

  int projectionAxis = ("K0S" == hadron)*K0SmassAxis + ("Lambda0" == hadron)*Lambda0massAxis + ("AntiLambda0" == hadron)*AntiLambda0massAxis;
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 1.07, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 1.1;
  if ("K0S" == hadron) { xMinFrame = 0.4, xMaxFrame = 0.6; }
  double xMinLegend = 0.6, xMaxLegend = 0.85, yMinLegend = 0.45, yMaxLegend = 0.7;
  double xLatex = 0.23, yLatex = 0.93;
  int xCanvas = 1800, yCanvas = 900;
  xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  yTitle = "arb. units";

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Line for hadron mass
  TLine* lineMass = new TLine(("K0S" == hadron) ? MassK0S : MassLambda0, yMinFrame, ("K0S" == hadron) ? MassK0S : MassLambda0, yMaxFrame);
  lineMass->SetLineColor(GetColor(0));
  lineMass->SetLineWidth(2);
  lineMass->SetLineStyle(9);

  histName = "V0PtMass";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
  TH1D* mass = (TH1D*)thn->Projection(projectionAxis);
  mass->Rebin(rebinNumber);

  double normalisationFactor = mass->GetBinContent(mass->GetMaximumBin());
  mass->Scale(1./normalisationFactor);
  setStyle(mass, 0);
  legend->AddEntry(mass, "data");

  canvas->cd();
  frame->Draw();
  mass->Draw("same");

  TF1* fit = getSigBkgFit(hadron, bkg, doubleGauss);
  TFitResultPtr fitResult = mass->Fit(fit, "RBS");
  double mean  = fit->GetParameter(1);
  double sigma = fit->GetParameter(2);
  double chiSq = fit->GetChisquare() / fit->GetNDF();
  if (doubleGauss && flipGaussians) {
    mean  = fit->GetParameter(4);
    sigma = fit->GetParameter(5);
  }

  array<TH1*, 2> regions = getSigBkgRegions(mass, mean, sigma, nSigmaSignal, nSigmaBkgMin, nSigmaBkgMax);
  TH1D* sigRegion = (TH1D*)regions[0];
  TH1D* bkgRegion = (TH1D*)regions[1];
  sigRegion->SetFillColorAlpha(kGreen, 0.3);
  bkgRegion->SetFillColorAlpha(kRed, 0.3);
  legend->AddEntry(sigRegion, "Signal+Background", "f");
  legend->AddEntry(bkgRegion, "Background", "f");
  sigRegion->Draw("bars same");
  bkgRegion->Draw("bars same");

  double sigPlusBkg = sigRegion->Integral();
  double background = bkgRegion->Integral();
  double purity = 1. - background/sigPlusBkg;
  cout << "Mean: " << mean << endl
       << "Sigma: " << sigma << endl
       << "Sig+Bkg: " << sigPlusBkg << endl
       << "Background: " << background << endl
       << "Purity: " << purity << endl;

  string sMean = TString::Format("#mu: %.3f GeV/#it{c}^{2}", mean).Data();
  if ("K0S" == hadron) {
    sMean = TString::Format("#mu: %.3f MeV/#it{c}^{2}", mean*1e3).Data();
  }
  string sSigma = TString::Format("#sigma: %.3f MeV/#it{c}^{2}", sigma*1e3).Data();
  TLatex* lData   = CreateLatex(0.55, 0.85, dataSet.c_str(), textSize);
  TLatex* lPt     = CreateLatex(0.55, 0.8, TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", lowpt, highpt), textSize);
  TLatex* lMean   = CreateLatex(0.25, 0.8, sMean.c_str(), textSize);
  TLatex* lSigma  = CreateLatex(0.25, 0.75, sSigma.c_str(), textSize);
  TLatex* lChiSq  = CreateLatex(0.25, 0.7, TString::Format("#chi^{2}/NDF: %.3f", chiSq), textSize);
  TLatex* lPurity = CreateLatex(0.25, 0.65, TString::Format("Purity: %.1f%%", purity*1e2), textSize);

  if (legend) legend->Draw("same");
  lineMass->Draw("same");
  lData->Draw("same");
  lPt->Draw("same");
  lPurity->Draw("same");
  lMean->Draw("same");
  lSigma->Draw("same");
  lChiSq->Draw("same");

  saveName = TString::Format("%sPurity", hadron.c_str());
  saveName = TString::Format("%s_pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax);
  if (doubleGauss) { saveName = TString::Format("%s_doubleGauss", saveName.c_str()); }
  saveName = TString::Format("%s_pol%d", saveName.c_str(), bkg);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  canvas->SaveAs(saveName.c_str());
}

void K0SPurity(string inName, string dataSet, double ptmin = 0., double ptmax = 100., int bkg = 1, bool doubleGauss = true, bool flipGaussians = false)
{ v0Purity(inName, dataSet, "K0S", ptmin, ptmax, bkg, doubleGauss, flipGaussians); }
void Lambda0Purity(string inName, string dataSet, double ptmin = 0., double ptmax = 100., int bkg = 1, bool doubleGauss = true, bool flipGaussians = false)
{ v0Purity(inName, dataSet, "Lambda0", ptmin, ptmax, bkg, doubleGauss, flipGaussians); }
void AntiLambda0Purity(string inName, string dataSet, double ptmin = 0., double ptmax = 100., int bkg = 1, bool doubleGauss = true, bool flipGaussians = false)
{ v0Purity(inName, dataSet, "AntiLambda0", ptmin, ptmax, bkg, doubleGauss, flipGaussians); }

// -------------------------------------------------------------------------------------------------

void purity225406(string hadron, double ptmin, double ptmax, int rebinNumber = 1, bool doubleGauss = true, bool flipGaussians = false)
{
  string inName = "~/cernbox/TrainOutput/225406/AnalysisResults.root";
  string dataSet = "LHC22o_pass6_minBias_small";
  v0Purity(inName, dataSet, hadron, ptmin, ptmax, 2, doubleGauss, flipGaussians, rebinNumber);
}

void purity233079(string hadron, double ptmin, double ptmax, int rebinNumber = 1, bool doubleGauss = true, bool flipGaussians = false)
{
  string inName = "~/cernbox/TrainOutput/233079/AnalysisResults.root";
  string dataSet = "LHC22o_pass6_minBias_small";
  v0Purity(inName, dataSet, hadron, ptmin, ptmax, 2, doubleGauss, flipGaussians, rebinNumber);
}