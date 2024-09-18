
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

void myPlotter(TCanvas* canvas, TH1F* frame, TH1* h, TLegend* legend, vector<TLatex*> latex)
{
  frame->Draw();
  h->Draw("same");
  if (legend) legend->Draw("same");
  for (auto l : latex) {
    l->Draw("same");
  }
  canvas->SaveAs(canvas->GetName());
}
void myPlotter(TCanvas* canvas, TH1* h, TLegend* legend, vector<TLatex*> latex)
{
  TH1F* frame = DrawFrame(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                          0, 1.1 * h->GetBinContent(h->GetMaximumBin()),
                          h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle());
  frame->SetTitle(h->GetTitle());
  myPlotter(canvas, frame, h, legend, latex);
}

// -------------------------------------------------------------------------------------------------

void plotConePt(vector<string> inputStrings)
{
  string fileName = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron = inputStrings[2];
  string canvasName = inputStrings[3]; // Also used as saveName
  vector<TLatex*> latex;

  TFile* file = TFile::Open(fileName.c_str());
  TH3D* h3D = (TH3D*)file->Get("jet-fragmentation/data/PC/ConePtEtaPhi");

  double etamin = -0.35;
  double etamax =  0.35;

  array<int, 2> etaBins = getProjectionBins(h3D->GetYaxis(), etamin, etamax);
  array<int, 2> phiBins = {0, h3D->GetNbinsZ() + 1};
  TH1D* hConePt = h3D->ProjectionX("hConePt", etaBins[0], etaBins[1], phiBins[0], phiBins[1]);
  setStyle(hConePt, 0);

  canvasName += ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
  myPlotter(canvas, hConePt, nullptr, latex);
}

// -------------------------------------------------------------------------------------------------

// void purity225406(string hadron, double ptmin, double ptmax, int rebinNumber = 1, bool doubleGauss = true, bool flipGaussians = false)
// {
//   string inName = "~/cernbox/TrainOutput/225406/AnalysisResults.root";
//   string dataSet = "LHC22o_pass6_minBias_small";
//   v0Purity(inName, dataSet, hadron, ptmin, ptmax, 2, doubleGauss, flipGaussians, rebinNumber);
// }

// void purity233079(string hadron, double ptmin, double ptmax, int rebinNumber = 1, bool doubleGauss = true, bool flipGaussians = false)
// {
//   string inName = "~/cernbox/TrainOutput/233079/AnalysisResults.root";
//   string dataSet = "LHC22o_pass6_minBias_small";
//   v0Purity(inName, dataSet, hadron, ptmin, ptmax, 2, doubleGauss, flipGaussians, rebinNumber);
// }