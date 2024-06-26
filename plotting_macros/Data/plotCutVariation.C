
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

const double MassK0S     = 0.497611;
const double MassLambda0 = 1.115683;

const int nDim                = 10;
const int ptAxis              =  0;
const int K0SmassAxis         =  1;
const int Lambda0massAxis     =  2;
const int AntiLambda0massAxis =  3;
const int RAxis               =  4;
const int ctauAxis            =  5;
const int cosPAAxis           =  6;
const int DCApAxis            =  7;
const int DCAnAxis            =  8;
const int DCAdAxis            =  9;

const vector<string> axisNames = {"pt", "K0Smass", "Lambda0mass", "AntiLambda0mass", "R", "ctau", "cosPA", "DCAp", "DCAn", "DCAd"};

// -------------------------------------------------------------------------------------------------

// Check if the cut is reversed (var < binUpEdge)
bool reverseCuts(int cutAxis)
{
  return ( (ctauAxis == cutAxis) || (DCAdAxis == cutAxis) );
}
// Format text for variable selections
string formatLatexText(int cutAxis, double binEdge)
{
  string latexText;
  switch (cutAxis) {
    case RAxis:
      latexText = TString::Format("%s > %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case ctauAxis:
      latexText = TString::Format("%s < %.0f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case cosPAAxis:
      latexText = TString::Format("%s > %.3f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case DCApAxis:
      latexText = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case DCAnAxis:
      latexText = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case DCAdAxis:
      latexText = TString::Format("%s < %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
  } // switch (cutAxis)
  return latexText;
}
// Format bin label for a given cut using formatLatexText
string formatBinLabel(THnSparseD* thn, int iBin, int cutAxis)
{
  if (0 == iBin) {
    return "No cut";
  }
  double binEdge = thn->GetAxis(cutAxis)->GetBinLowEdge(iBin);
  if (reverseCuts(cutAxis)) {
    binEdge = thn->GetAxis(cutAxis)->GetBinUpEdge(iBin);
  }
  return formatLatexText(cutAxis, binEdge);
}

// Returns a subset of a histogram from minBin to maxBin
template <typename T>
T* makeHistSubset(T* data, int minBin, int maxBin)
{
  T* region = (T*)data->Clone("region");
  region->Reset();
  for (int i = minBin; i <= maxBin; i++) {
    region->SetBinContent(i, data->GetBinContent(i));
  }
  return region;
}

// Initialise parameters and set bounds
void parsK0S(TF1* f, bool doubleGauss, int bkg, double scale = 1.)
{
  double gausPar[3] = {1. * scale, MassK0S, 0.01};
  double gausVar[6] = {0.4 * scale, 1. * scale, MassK0S - 5e-2, MassK0S + 5e-2, 1e-3, 0.02};
  double doubleGaussPar[3] = {0.1 * scale, MassK0S, 0.05};
  double doubleGaussVar[6] = {0.1 * scale, 0.6 * scale, MassK0S - 0.01, MassK0S + 0.01, 1e-3, 0.1};
  double bkgPar[3] = {0., 0., 0.};
  double bkgVar[6] = {-0.1 * scale, 0.1 * scale, -0.1 * scale, 0.1 * scale, -0.1 * scale, 0.1 * scale};

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
void parsL0(TF1* f, bool doubleGauss, int bkg, double scale = 1.)
{
  double gausPar[3] = {1. * scale, MassLambda0, 0.01};
  double gausVar[6] = {0.4 * scale, 1. * scale, MassLambda0 - 0.01, MassLambda0 + 0.01, 1e-4, 0.02};
  double doubleGaussPar[3] = {0.1 * scale, MassLambda0, 0.05};
  double doubleGaussVar[6] = {0.1 * scale, 0.6 * scale, MassLambda0 - 0.03, MassLambda0 + 0.03, 1e-3, 0.1};
  double bkgPar[3] = {0., 0., 0.};
  double bkgVar[6] = {-3. * scale, 3. * scale, -3. * scale, 3. * scale, -3. * scale, 3. * scale};

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
// Get fit function for signal and background, still need to initialise parameters
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
  return fit;
}
// Get shaded regions for signal and background, uses 2 background regions symmetric around the mean
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
  TH1D* sigRegion = (TH1D*)makeHistSubset(mass, sigBins[0], sigBins[1]);
  TH1D* bkgRegion = (TH1D*)makeHistSubset(mass, leftBkgBins[0], leftBkgBins[1]);
  bkgRegion->Add(makeHistSubset(mass, rightBkgBins[0], rightBkgBins[1]));
  return {sigRegion, bkgRegion};
}

// -------------------------------------------------------------------------------------------------
// Fit the signal+background for a mass spectrum. Returns mean and sigma of the signal peak
array<double, 2> doFitForRefBin(TH1* mass, string hadron, array<int, 3> fitSettings, TLegend* legend)
{
  double mean = -1., sigma = -1.;
  double scale = mass->GetBinContent(mass->GetMaximumBin());
  int bkg = fitSettings[0];
  bool doubleGauss = fitSettings[1], flipGaussians = fitSettings[2];
  TF1* fit = getSigBkgFit(hadron, bkg, doubleGauss);
  if ("K0S" == hadron) {
    parsK0S(fit, doubleGauss, bkg, scale);
  }
  else {
    parsL0(fit, doubleGauss, bkg, scale);
  }
  TFitResultPtr fitResult = mass->Fit(fit, "RBSQ");
  legend->AddEntry(fit, "Fit", "l");
  mean  = fit->GetParameter(1);
  sigma = fit->GetParameter(2);
  if (doubleGauss && flipGaussians) {
    mean  = fit->GetParameter(4);
    sigma = fit->GetParameter(5);
  }
  return {mean, sigma};
}
// Returns the mass spectrum for a given bin (i.e. cut) and projection axis
TH1D* getMassHist(THnSparseD* thn, array<int, 2> ptBins, int iBin, int cutAxis, int projectionAxis)
{
  THnSparseD* thn_copy = (THnSparseD*)thn->Clone("thn_copy");
  thn_copy->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);

  int overflowBin = 1 + thn->GetAxis(cutAxis)->GetNbins();
  int firstBin = 0, lastBin = overflowBin;
  if (reverseCuts(cutAxis)) {
    lastBin = overflowBin - iBin;
  }
  else {
    firstBin = iBin;
  }
  thn_copy->GetAxis(cutAxis)->SetRange(firstBin, lastBin);
  TH1D* mass = (TH1D*)thn_copy->Projection(projectionAxis);
  mass->SetName(TString::Format("mass_axis%d_bin%d", cutAxis, iBin).Data());
  return mass;
}
// -------------------------------------------------------------------------------------------------
// For a given bin (i.e. cut), calculate the purity and relative efficiency
array<double, 3> getPurityAndRelEfficiency(THnSparseD* thn, array<int, 2> ptBins,
                                           int iBin, int refBin, int cutAxis, int rebinNumber,
                                           string hadron, string dataSet,
                                           array<int, 3> fitSettings, array<int, 3> nSigma,
                                           array<double, 3> refValues,
                                           TH1D* hPurity, TH1D* hEfficiency)
{
  // Setup
  double textSize  = 0.04;
  double xMinFrame = 1.07, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 1.1;
  if ("K0S" == hadron) { xMinFrame = 0.4, xMaxFrame = 0.6; }
  double xMinLegend = 0.1, xMaxLegend = 0.7, yMinLegend = 0.1, yMaxLegend = 0.45;
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
  int xCanvas = 1800, yCanvas = 900;
  int projectionAxis = ("K0S" == hadron)*K0SmassAxis + ("Lambda0" == hadron)*Lambda0massAxis + ("AntiLambda0" == hadron)*AntiLambda0massAxis;
  string sPt = TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", lowpt, highpt).Data();
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "arb. units";
  string canvasName = TString::Format("canvas_%s_%s_%d", hadron.c_str(), axisNames[cutAxis].c_str(), iBin).Data();

  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), xCanvas, yCanvas);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, "", textSize);
  canvas->Divide(2,1);
  canvas->cd(1);
  gPad->SetRightMargin(0.05);

  TH1D* mass = getMassHist(thn, ptBins, iBin, cutAxis, projectionAxis);
  mass->Rebin(rebinNumber);
  double scale = mass->GetBinContent(mass->GetMaximumBin());
  setStyle(mass, 0);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, 1.1 * scale, xTitle, yTitle);
  frame->Draw();
  mass->Draw("same");
  legend->AddEntry(mass, "Data");

  // Fit mass peak only for the reference bin
  if (refBin == iBin) {
    array<double, 2> meanAndSigma = doFitForRefBin(mass, hadron, fitSettings, legend);
    refValues[0] = meanAndSigma[0];
    refValues[1] = meanAndSigma[1];
  }

  double mean         = refValues[0];
  double sigma        = refValues[1];
  double nSigmaSignal = nSigma[0];
  double nSigmaBkgMin = nSigma[1];
  double nSigmaBkgMax = nSigma[2];
  // Get Signal and Background regions
  array<TH1*, 2> regions = getSigBkgRegions(mass, mean, sigma, nSigmaSignal, nSigmaBkgMin, nSigmaBkgMax);
  TH1D* sigRegion = (TH1D*)regions[0];
  TH1D* bkgRegion = (TH1D*)regions[1];
  sigRegion->SetFillColorAlpha(kGreen, 0.3);
  bkgRegion->SetFillColorAlpha(kRed, 0.3);
  legend->AddEntry(sigRegion, "Signal+Background", "f");
  legend->AddEntry(bkgRegion, "Background", "f");
  sigRegion->Draw("bars same");
  bkgRegion->Draw("bars same");

  // Line for hadron mass
  TLine* lineMass = new TLine(("K0S" == hadron) ? MassK0S : MassLambda0, yMinFrame, ("K0S" == hadron) ? MassK0S : MassLambda0, yMaxFrame);
  lineMass->SetLineColor(GetColor(0));
  lineMass->SetLineWidth(2);
  lineMass->SetLineStyle(9);
  lineMass->Draw("same");

  // Calculate purity and relative efficiency
  double sigPlusBkg = sigRegion->Integral();
  double background = bkgRegion->Integral();
  double purity = 1. - background/sigPlusBkg;
  if (refBin == iBin) {
    refValues[2] = sigPlusBkg * purity;
  }
  double efficiency = sigPlusBkg * purity / refValues[2];

  // Format latex text, same as bin label
  string binLabel = formatBinLabel(thn, iBin, cutAxis);
  hPurity->SetBinContent(iBin + 1, purity);
  hPurity->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());
  hEfficiency->SetBinContent(iBin + 1, efficiency);
  hEfficiency->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());

  // On pad 2, we put important text and the legend
  canvas->cd(2);
  gPad->SetLeftMargin(0.);
  TLatex* lData   = CreateLatex(0.1, 0.85, dataSet.c_str(), textSize);
  TLatex* lPt     = CreateLatex(0.1, 0.8, sPt.c_str(), textSize);
  TLatex* lCut    = CreateLatex(0.1, 0.7, binLabel.c_str(), textSize);
  TLatex* lPurity = CreateLatex(0.1, 0.65, TString::Format("Purity: %.1f%%", purity*1e2), textSize);
  TLatex* lEff    = CreateLatex(0.1, 0.6, TString::Format("Rel. Efficiency: %.2f %%", 1e2*efficiency).Data(), textSize);
  if (legend) legend->Draw("same");
  lData->Draw("same");
  lPt->Draw("same");
  lCut->Draw("same");
  lPurity->Draw("same");
  lEff->Draw("same");

  string saveName = TString::Format("purity%s", hadron.c_str()).Data();
  saveName = TString::Format("%s_%s", saveName.c_str(), axisNames[cutAxis].c_str()).Data();
  saveName = TString::Format("%s_bin%d", saveName.c_str(), iBin).Data();
  saveName = TString::Format("%s_v0pt%.1f-%.1f", saveName.c_str(), lowpt, highpt).Data();
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(saveName.c_str());
  return refValues;
}
// -------------------------------------------------------------------------------------------------
// Makes plots and extracts purity and relative efficiency for a given cut axis
void cutVarPurity(string inName, string dataSet, string hadron, int cutAxis, double ptmin, double ptmax, int bkg, bool doubleGauss, bool flipGaussians)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if (cutAxis < 4 || cutAxis > 9) {
    cout << "Error: cutAxis must be specified (4-9)" << endl;
    return;
  }

  int nSigmaSignal = 3, nSigmaBkgMin = 5, nSigmaBkgMax = 8;
  int refBin = 0;
  int rebinNumber = 1;
  double mean = ("K0S" == hadron) ? MassK0S : MassLambda0, sigma = 0.03, refSignal = -1.;

  array<int, 3>    fitSettings = {bkg, doubleGauss, flipGaussians};
  array<int, 3>    nSigma      = {nSigmaSignal, nSigmaBkgMin, nSigmaBkgMax};
  array<double, 3> refValues   = {mean, sigma, refSignal};

  gStyle->SetNdivisions(505, "xy");
  string histName, histTitle, legendTitle, latexText;
  histName = "V0CutVariation";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);

  int underflowBin = 0, overflowBin = 1 + thn->GetAxis(cutAxis)->GetNbins();
  TH1D* hEfficiency = new TH1D("hEfficiency", "", overflowBin+1, thn->GetAxis(cutAxis)->GetBinLowEdge(0), thn->GetAxis(cutAxis)->GetBinUpEdge(overflowBin));
  TH1D* hPurity = (TH1D*)hEfficiency->Clone("hPurity");

  // Do reference bin first
  refValues = getPurityAndRelEfficiency(thn, ptBins, refBin, refBin, cutAxis, rebinNumber, hadron, dataSet, fitSettings, nSigma, refValues, hPurity, hEfficiency);
  for (int i = underflowBin; i <= overflowBin; i++) {
    if (i == refBin) { continue; }
    refValues = getPurityAndRelEfficiency(thn, ptBins, i, refBin, cutAxis, rebinNumber, hadron, dataSet, fitSettings, nSigma, refValues, hPurity, hEfficiency);
  }

  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);

  TCanvas* cEP   = new TCanvas("cEP", "cEP", 1800, 900);
  cEP->cd();
  TLatex*  lData = CreateLatex(0.25, 0.35, dataSet.c_str(), 0.04);
  TLatex*  lPt   = CreateLatex(0.25, 0.3, TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", lowpt, highpt).Data(), 0.04);
  TLegend* lEP   = CreateLegend(0.25, 0.5, 0.15, 0.25, "", 0.04);

  setStyle(hEfficiency, 0);
  setStyle(hPurity, 1);
  lEP->AddEntry(hEfficiency, "Relative Efficiency");
  lEP->AddEntry(hPurity, "Purity");

  hEfficiency->SetStats(0);
  hEfficiency->GetYaxis()->SetRangeUser(0., 1.3);
  hEfficiency->SetMarkerSize(2.5);
  hEfficiency->Draw();
  hEfficiency->Draw("same text00");

  hPurity->SetMarkerSize(2.5);
  hPurity->Draw("same");
  hPurity->Draw("same text00");

  lData->Draw("same");
  lPt->Draw("same");
  lEP->Draw("same");

  string saveName = "eff-pur";
  saveName += hadron;
  saveName += "_" + axisNames[cutAxis];
  saveName += TString::Format("_v0pt%.1f-%.1f", lowpt, highpt).Data();
  saveName += ".pdf";
  cEP->SaveAs(saveName.c_str());
}
