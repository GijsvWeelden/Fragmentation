
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

// Check if input is ok
bool inputIssue(string inName, string hadron, int cutAxis)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return true;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return true;
  }
  if (cutAxis < 4 || cutAxis > 9) {
    cout << "Error: cutAxis must be specified (4-9)" << endl;
    return true;
  }
  return false;
}

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
// THn version: Format bin label for a given cut using formatLatexText
string formatBinLabel(THnSparseD* thn, int iBin, int cutAxis)
{
  if (0 == iBin) {
    return "Default cut";
  }
  double binEdge = thn->GetAxis(cutAxis)->GetBinLowEdge(iBin);
  if (reverseCuts(cutAxis)) {
    binEdge = thn->GetAxis(cutAxis)->GetBinUpEdge(iBin);
  }
  return formatLatexText(cutAxis, binEdge);
}
// TH3 version: Format bin label for a given cut using formatLatexText
string formatBinLabel(TH3D* th3, int iBin, int cutAxis)
{
  if (0 == iBin) {
    return "Default cut";
  }
  double binEdge = th3->GetYaxis()->GetBinLowEdge(iBin);
  if (reverseCuts(cutAxis)) {
    binEdge = th3->GetYaxis()->GetBinUpEdge(iBin);
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
    if (bkg >= 2) {
      f->SetParameter(8, bkgPar[2]);       f->SetParLimits(8, bkgVar[4], bkgVar[5]);
    }
  }
  else {
    f->SetParameter(3, bkgPar[0]);   f->SetParLimits(3, bkgVar[0], bkgVar[1]);
    f->SetParameter(4, bkgPar[1]);   f->SetParLimits(4, bkgVar[2], bkgVar[3]);
    if (2 <= bkg) {
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
    if (bkg >= 2) {
      f->SetParameter(8, bkgPar[2]);       f->SetParLimits(8, bkgVar[4], bkgVar[5]);
    }
  }
  else {
    f->SetParameter(3, bkgPar[0]);   f->SetParLimits(3, bkgVar[0], bkgVar[1]);
    f->SetParameter(4, bkgPar[1]);   f->SetParLimits(4, bkgVar[2], bkgVar[3]);
    if (bkg >= 2) {
      f->SetParameter(5, bkgPar[2]); f->SetParLimits(5, bkgVar[4], bkgVar[5]);
    }
  }
}
// Get fit function for signal and background, still need to initialise parameters
TF1* getSigBkgFit(string hadron, int bkg, bool doubleGauss)
{
  double fitRegion[2] = {1.08, 1.2};
  if ("K0S" == hadron) { fitRegion[0] = 0.42, fitRegion[1] = 0.58; }

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
array<double, 3> doFitForRefBin(TH1* mass, string hadron, array<int, 3> fitSettings, TLegend* legend)
{
  double mean = -1., sigma = -1., chiSq = -1.;
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
  chiSq = fitResult->Chi2() / fitResult->Ndf();
  return {mean, sigma, chiSq};
}
// THn version: Returns the mass spectrum for a given bin (i.e. cut) and projection axis
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
// TH3 version: Returns the mass spectrum for a given bin (i.e. cut) and projection axis
TH1D* getMassHist(TH3D* th3, array<int, 2> ptBins, int iBin, int cutAxis)
{
  TH3D* th3_copy = (TH3D*)th3->Clone("th3_copy");
  int overflowBin = 1 + th3->GetYaxis()->GetNbins();
  int firstBin = 0, lastBin = overflowBin;
  if (reverseCuts(cutAxis)) {
    lastBin = overflowBin - iBin;
  }
  else {
    firstBin = iBin;
  }
  TH1D* mass = (TH1D*)th3_copy->ProjectionZ("mass", ptBins[0], ptBins[1], firstBin, lastBin);
  mass->SetName(TString::Format("mass_axis%d_bin%d", cutAxis, iBin).Data());
  return mass;
}
// -------------------------------------------------------------------------------------------------
// THn version: For a given bin (i.e. cut), calculate the purity and relative efficiency
array<double, 3> getPurityAndRelEfficiency(THnSparseD* thn, array<int, 2> ptBins,
                                           int iBin, int refBin, int cutAxis, int rebinNumber,
                                           string hadron, string dataSet,
                                           array<int, 3> fitSettings, array<int, 3> nSigma,
                                           array<double, 3> refValues,
                                           array<TH1D*, 4> hists)
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
  string yTitle = "counts";
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

  double chiSq = -1.;
  // Fit mass peak only for the reference bin
  if (refBin == iBin) {
    array<double, 3> fitResults = doFitForRefBin(mass, hadron, fitSettings, legend);
    refValues[0] = fitResults[0];
    refValues[1] = fitResults[1];
    chiSq = fitResults[2];
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
  legend->AddEntry(sigRegion, TString::Format("Signal+Background (#mu #pm %.0f #sigma)", nSigmaSignal).Data(), "f");
  legend->AddEntry(bkgRegion, TString::Format("Background (#mu #pm %.0f - %.0f #sigma)", nSigmaBkgMin, nSigmaBkgMax).Data(), "f");
  sigRegion->Draw("bars same");
  bkgRegion->Draw("bars same");

  // Line for hadron mass
  TLine* lineMass = new TLine(("K0S" == hadron) ? MassK0S : MassLambda0, yMinFrame, ("K0S" == hadron) ? MassK0S : MassLambda0, yMaxFrame * scale);
  lineMass->SetLineColor(GetColor(0));
  lineMass->SetLineWidth(2);
  lineMass->SetLineStyle(9);
  lineMass->Draw("same");

  // Calculate purity and relative efficiency
  double sigPlusBkg   = sigRegion->Integral();
  double background   = bkgRegion->Integral();
  double signal       = sigPlusBkg - background;
  double purity       = signal / sigPlusBkg;
  double significance = signal / TMath::Sqrt(sigPlusBkg);
  double sigOverBkg   = signal / ( (background > 0) ? background : -1.);

  if (refBin == iBin) {
    refValues[2] = signal;
  }
  double efficiency = signal / refValues[2];
  if (efficiency > 1) {
    cout << "Warning: relative efficiency > 1 (" << efficiency << "). Signal: " << signal << ", reference signal: " << refValues[2] << endl;
  }

  // Format latex text, same as bin label
  string binLabel = formatBinLabel(thn, iBin, cutAxis);
  hists[0]->SetBinContent(iBin + 1, purity);
  hists[0]->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());
  hists[1]->SetBinContent(iBin + 1, efficiency);
  hists[1]->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());
  hists[2]->SetBinContent(iBin + 1, significance);
  hists[2]->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());
  hists[3]->SetBinContent(iBin + 1, sigOverBkg);
  hists[3]->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());

  // On pad 2, we put important text and the legend
  canvas->cd(2);
  gPad->SetLeftMargin(0.);
  TLatex* lData   = CreateLatex(0.1, 0.9, dataSet.c_str(), textSize);
  TLatex* lPt     = CreateLatex(0.1, 0.85, sPt.c_str(), textSize);
  TLatex* lCut    = CreateLatex(0.1, 0.75, binLabel.c_str(), textSize);
  TLatex* lPurity = CreateLatex(0.1, 0.7, TString::Format("Purity: %.1f%%", purity*1e2), textSize);
  TLatex* lEff    = CreateLatex(0.1, 0.65, TString::Format("Rel. Efficiency: %.2f%%", 1e2*efficiency).Data(), textSize);
  if (legend) legend->Draw("same");
  lData->Draw("same");
  lPt->Draw("same");
  lCut->Draw("same");
  lPurity->Draw("same");
  lEff->Draw("same");

  if (refBin == iBin) {
    string sMean  = TString::Format("Mean: %.3f GeV/#it{c}^{2}", mean).Data();
    string sSigma = TString::Format("Sigma: %.3f GeV/#it{c}^{2}", sigma).Data();
    if ("K0S" == hadron) {
      sMean  = TString::Format("Mean: %.3f MeV/#it{c}^{2}", mean*1e3).Data();
      sSigma = TString::Format("Sigma: %.3f MeV/#it{c}^{2}", sigma*1e3).Data();
    }
    TLatex* lMean  = CreateLatex(0.1, 0.6, sMean.c_str(), textSize);
    TLatex* lSigma = CreateLatex(0.1, 0.55, sSigma.c_str(), textSize);
    TLatex* lChiSq = CreateLatex(0.1, 0.5, TString::Format("#chi^{2}/NDF: %.2f", chiSq).Data(), textSize);
    lMean->Draw("same");
    lSigma->Draw("same");
    lChiSq->Draw("same");
  }

  string saveName = hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowpt, highpt).Data();
  saveName += "_" + axisNames[cutAxis];
  saveName = TString::Format("%s_bin%d", saveName.c_str(), iBin).Data();
  saveName += ".pdf";
  canvas->SaveAs(saveName.c_str());
  return refValues;
}

// TH3 version: For a given bin (i.e. cut), calculate the purity and relative efficiency
// Allows you to put in a distribution that's been cut on another axis
array<double, 3> getPurityAndRelEfficiency(TH3D* th3, array<int, 2> ptBins,
                                           int iBin, int refBin, int cutAxis, int rebinNumber,
                                           string hadron, string dataSet, string cutString,
                                           array<int, 3> fitSettings, array<int, 3> nSigma,
                                           array<double, 3> refValues,
                                           array<TH1D*, 4> hists)
{
  // Input hist must be {pt, cutAxis, M}
  // Setup
  double textSize  = 0.04;
  double xMinFrame = 1.07, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 1.1;
  if ("K0S" == hadron) { xMinFrame = 0.4, xMaxFrame = 0.6; }
  double xMinLegend = 0.1, xMaxLegend = 0.7, yMinLegend = 0.1, yMaxLegend = 0.45;
  double lowpt = th3->GetXaxis()->GetBinLowEdge(ptBins[0]);
  double highpt = th3->GetXaxis()->GetBinUpEdge(ptBins[1]);
  int xCanvas = 1800, yCanvas = 900;
  string sPt = TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", lowpt, highpt).Data();
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";
  string canvasName = TString::Format("canvas_%s_%s_%d", hadron.c_str(), axisNames[cutAxis].c_str(), iBin).Data();

  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), xCanvas, yCanvas);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, "", textSize);
  canvas->Divide(2,1);
  canvas->cd(1);
  gPad->SetRightMargin(0.05);

  TH1D* mass = getMassHist(th3, ptBins, iBin, cutAxis);
  mass->Rebin(rebinNumber);
  double scale = mass->GetBinContent(mass->GetMaximumBin());
  setStyle(mass, 0);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, 1.1 * scale, xTitle, yTitle);
  frame->Draw();
  mass->Draw("same");
  legend->AddEntry(mass, "Data");

  double chiSq = -1.;
  // Fit mass peak only for the reference bin
  if (refBin == iBin) {
    array<double, 3> fitResults = doFitForRefBin(mass, hadron, fitSettings, legend);
    refValues[0] = fitResults[0];
    refValues[1] = fitResults[1];
    chiSq = fitResults[2];
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
  legend->AddEntry(sigRegion, TString::Format("Signal+Background (#mu #pm %.0f #sigma)", nSigmaSignal).Data(), "f");
  legend->AddEntry(bkgRegion, TString::Format("Background (#mu #pm %.0f - %.0f #sigma)", nSigmaBkgMin, nSigmaBkgMax).Data(), "f");
  sigRegion->Draw("bars same");
  bkgRegion->Draw("bars same");

  // Line for hadron mass
  TLine* lineMass = new TLine(("K0S" == hadron) ? MassK0S : MassLambda0, yMinFrame, ("K0S" == hadron) ? MassK0S : MassLambda0, yMaxFrame * scale);
  lineMass->SetLineColor(GetColor(0));
  lineMass->SetLineWidth(2);
  lineMass->SetLineStyle(9);
  lineMass->Draw("same");

  // Calculate purity and relative efficiency
  double sigPlusBkg = sigRegion->Integral();
  double background = bkgRegion->Integral();
  double signal     = sigPlusBkg - background;
  double purity     = signal / sigPlusBkg;
  double significance = signal / TMath::Sqrt(sigPlusBkg);
  double sigOverBkg   = signal / ( (background > 0) ? background : -1.);

  if (refBin == iBin) {
    refValues[2] = signal;
  }
  double efficiency = signal / refValues[2];
  if (efficiency > 1) {
    cout << "Warning: relative efficiency > 1 (" << efficiency << "). Signal: " << signal << ", reference signal: " << refValues[2] << endl;
  }

  // Format latex text, same as bin label
  string binLabel = formatBinLabel(th3, iBin, cutAxis);
  hists[0]->SetBinContent(iBin + 1, purity);
  hists[0]->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());
  hists[1]->SetBinContent(iBin + 1, efficiency);
  hists[1]->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());
  hists[2]->SetBinContent(iBin + 1, significance);
  hists[2]->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());
  hists[3]->SetBinContent(iBin + 1, sigOverBkg);
  hists[3]->GetXaxis()->SetBinLabel(iBin + 1, binLabel.c_str());

  // On pad 2, we put important text and the legend
  canvas->cd(2);
  gPad->SetLeftMargin(0.);
  TLatex* lData   = CreateLatex(0.1, 0.9, dataSet.c_str(), textSize);
  TLatex* lPt     = CreateLatex(0.1, 0.85, sPt.c_str(), textSize);
  TLatex* lCut    = CreateLatex(0.1, 0.75, binLabel.c_str(), textSize);
  TLatex* lPurity = CreateLatex(0.1, 0.7, TString::Format("Purity: %.1f%%", purity*1e2), textSize);
  TLatex* lEff    = CreateLatex(0.1, 0.65, TString::Format("Rel. Efficiency: %.2f%%", 1e2*efficiency).Data(), textSize);
  if (legend) legend->Draw("same");
  lData->Draw("same");
  lPt->Draw("same");
  lCut->Draw("same");
  lPurity->Draw("same");
  lEff->Draw("same");

  TLatex* lCutString = CreateLatex(0.6, 0.75, cutString.c_str(), textSize);
  lCutString->Draw("same");

  if (refBin == iBin) {
    string sMean  = TString::Format("Mean: %.3f GeV/#it{c}^{2}", mean).Data();
    string sSigma = TString::Format("Sigma: %.3f GeV/#it{c}^{2}", sigma).Data();
    if ("K0S" == hadron) {
      sMean  = TString::Format("Mean: %.3f MeV/#it{c}^{2}", mean*1e3).Data();
      sSigma = TString::Format("Sigma: %.3f MeV/#it{c}^{2}", sigma*1e3).Data();
    }
    TLatex* lMean  = CreateLatex(0.1, 0.6, sMean.c_str(), textSize);
    TLatex* lSigma = CreateLatex(0.1, 0.55, sSigma.c_str(), textSize);
    TLatex* lChiSq = CreateLatex(0.1, 0.5, TString::Format("#chi^{2}/NDF: %.2f", chiSq).Data(), textSize);
    lMean->Draw("same");
    lSigma->Draw("same");
    lChiSq->Draw("same");
  }

  string saveName = hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowpt, highpt).Data();
  saveName += "_" + axisNames[cutAxis];
  saveName = TString::Format("%s_bin%d", saveName.c_str(), iBin).Data();
  saveName += ".pdf";
  canvas->SaveAs(saveName.c_str());
  return refValues;
}
// -------------------------------------------------------------------------------------------------
// Makes plots and extracts purity and relative efficiency for a given cut axis
void cutVarPurity(string inName, string dataSet, string hadron, int cutAxis, double ptmin, double ptmax, int bkg, bool doubleGauss, bool flipGaussians, int rebinNumber)
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

  int refBin = 0;
  int nSigmaSignal = 3, nSigmaBkgMin = 5, nSigmaBkgMax = nSigmaBkgMin + nSigmaSignal;
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
  TH1D* hSignificance = (TH1D*)hEfficiency->Clone("hSignificance");
  TH1D* hSigOverBkg = (TH1D*)hEfficiency->Clone("hSigOverBkg");
  array<TH1D*, 4> hists = {hPurity, hEfficiency, hSignificance, hSigOverBkg};

  // Do reference bin first
  refValues = getPurityAndRelEfficiency(thn, ptBins, refBin, refBin, cutAxis, rebinNumber, hadron, dataSet, fitSettings, nSigma, refValues, hists);
  for (int i = underflowBin; i <= overflowBin; i++) {
    if (i == refBin) { continue; }
    refValues = getPurityAndRelEfficiency(thn, ptBins, i, refBin, cutAxis, rebinNumber, hadron, dataSet, fitSettings, nSigma, refValues, hists);
  }

  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);

  string canvasName = TString::Format("cEP_%s_%s_%.1f-%.1f", hadron.c_str(), axisNames[cutAxis].c_str(), lowpt, highpt).Data();
  TCanvas* cEP   = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1800, 900);
  cEP->cd();
  TLatex*  lData = CreateLatex(0.25, 0.35, dataSet.c_str(), 0.04);
  TLatex*  lPt   = CreateLatex(0.25, 0.3, TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", lowpt, highpt).Data(), 0.04);
  TLegend* lEP   = CreateLegend(0.25, 0.5, 0.15, 0.25, "", 0.04);

  hPurity = hists[0];
  hEfficiency = hists[1];
  hSignificance = hists[2];
  hSigOverBkg = hists[3];
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

  string saveName = hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowpt, highpt).Data();
  saveName += "_" + axisNames[cutAxis];
  saveName += "_EffPur";
  saveName += ".pdf";
  cEP->SaveAs(saveName.c_str());
}

void cutVarPurityWithExtraCuts(string inName, string dataSet, string hadron, int cutAxis, double ptmin, double ptmax, int bkg, bool doubleGauss, bool flipGaussians, int rebinNumber)
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

  int nSigmaSignal = 3, nSigmaBkgMin = 5, nSigmaBkgMax = nSigmaBkgMin + nSigmaSignal;
  if ("K0S" == hadron) {
    nSigmaSignal = 3, nSigmaBkgMin = 5, nSigmaBkgMax = nSigmaBkgMin + nSigmaSignal;
  }
  int refBin = 0;
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
  TH1D* hSignificance = (TH1D*)hEfficiency->Clone("hSignificance");
  TH1D* hSigOverBkg = (TH1D*)hEfficiency->Clone("hSigOverBkg");
  array<TH1D*, 4> hists = {hPurity, hEfficiency, hSignificance, hSigOverBkg};

  // Apply other cuts first: Change these to whatever you want
  // Set up cutString to write what you did
  string cutString = "", axisString = "";
  int axisBin = 0;

  axisBin = 0;
  axisString = formatBinLabel(thn, axisBin, RAxis);
  thn->GetAxis(RAxis)->SetRange(axisBin, 1 + thn->GetAxis(RAxis)->GetNbins());
  if (!("Default cut" == axisString)) {
    cutString = axisString;
  }
  axisBin = 0;
  axisString = formatBinLabel(thn, axisBin, ctauAxis);
  thn->GetAxis(ctauAxis)->SetRange(0, -1 * axisBin + 1 + thn->GetAxis(ctauAxis)->GetNbins());
  if (!("Default cut" == axisString)) {
    cutString = TString::Format("#splitline{%s}{%s}", cutString.c_str(), axisString.c_str()).Data();
  }
  axisBin = 1;
  axisString = formatBinLabel(thn, axisBin, cosPAAxis);
  thn->GetAxis(cosPAAxis)->SetRange(axisBin, 1 + thn->GetAxis(cosPAAxis)->GetNbins());
  if (!("Default cut" == axisString)) {
    cutString = TString::Format("#splitline{%s}{%s}", cutString.c_str(), axisString.c_str()).Data();
  }
  axisBin = 0;
  axisString = formatBinLabel(thn, axisBin, DCApAxis);
  thn->GetAxis(DCApAxis)->SetRange(0, -1 * axisBin + 1 + thn->GetAxis(DCApAxis)->GetNbins());
  if (!("Default cut" == axisString)) {
    cutString = TString::Format("#splitline{%s}{%s}", cutString.c_str(), axisString.c_str()).Data();
  }
  axisBin = 0;
  axisString = formatBinLabel(thn, axisBin, DCAnAxis);
  thn->GetAxis(DCAnAxis)->SetRange(axisBin, 1 + thn->GetAxis(DCAnAxis)->GetNbins());
  if (!("Default cut" == axisString)) {
    cutString = TString::Format("#splitline{%s}{%s}", cutString.c_str(), axisString.c_str()).Data();
  }
  axisBin = 0;
  axisString = formatBinLabel(thn, axisBin, DCAnAxis);
  thn->GetAxis(DCAdAxis)->SetRange(axisBin, 1 + thn->GetAxis(DCAdAxis)->GetNbins());
  if (!("Default cut" == axisString)) {
    cutString = TString::Format("#splitline{%s}{%s}", cutString.c_str(), axisString.c_str()).Data();
  }

  int projectionAxis = ("K0S" == hadron)*K0SmassAxis + ("Lambda0" == hadron)*Lambda0massAxis + ("AntiLambda0" == hadron)*AntiLambda0massAxis;
  TH3D* th3 = (TH3D*) thn->Projection(ptAxis, cutAxis, projectionAxis);
  // Do reference bin first
  refValues = getPurityAndRelEfficiency(th3, ptBins, refBin, refBin, cutAxis, rebinNumber, hadron, dataSet, cutString, fitSettings, nSigma, refValues, hists);
  for (int i = underflowBin; i <= overflowBin; i++) {
    if (i == refBin) { continue; }
    refValues = getPurityAndRelEfficiency(th3, ptBins, i, refBin, cutAxis, rebinNumber, hadron, dataSet, cutString, fitSettings, nSigma, refValues, hists);
  }

  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);

  string canvasName = TString::Format("cEP_%s_%s_%.1f-%.1f", hadron.c_str(), axisNames[cutAxis].c_str(), lowpt, highpt).Data();
  TCanvas* cEP   = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1800, 900);
  cEP->cd();
  TLatex*  lData      = CreateLatex(0.25, 0.35, dataSet.c_str(), 0.04);
  TLatex*  lPt        = CreateLatex(0.25, 0.3, TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", lowpt, highpt).Data(), 0.04);
  TLatex*  lCutString = CreateLatex(0.7, 0.35, cutString.c_str(), 0.04);
  TLegend* lEP        = CreateLegend(0.25, 0.5, 0.15, 0.25, "", 0.04);

  hPurity = hists[0];
  hEfficiency = hists[1];
  hSignificance = hists[2];
  hSigOverBkg = hists[3];
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
  lCutString->Draw("same");
  lEP->Draw("same");

  string saveName = hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowpt, highpt).Data();
  saveName += "_" + axisNames[cutAxis];
  saveName += "_EffPur";
  saveName += ".pdf";
  cEP->SaveAs(saveName.c_str());
}

void cutVarPurityWithCompetingMass(string inName, string dataSet, string hadron, int cutAxis, double ptmin, double ptmax, int bkg, bool doubleGauss, bool flipGaussians, int rebinNumber)
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

  int refBin = 0;
  int nSigmaSignal = 3, nSigmaBkgMin = 5, nSigmaBkgMax = nSigmaBkgMin + nSigmaSignal;
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
  THnSparseD* thn2 = (THnSparseD*)thn->Clone("thn2");
  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);

  int underflowBin = 0, overflowBin = 1 + thn->GetAxis(cutAxis)->GetNbins();
  TH1D* hEfficiency = new TH1D("hEfficiency", "", overflowBin+1, thn->GetAxis(cutAxis)->GetBinLowEdge(0), thn->GetAxis(cutAxis)->GetBinUpEdge(overflowBin));
  TH1D* hPurity = (TH1D*)hEfficiency->Clone("hPurity");
  TH1D* hSignificance = (TH1D*)hEfficiency->Clone("hSignificance");
  TH1D* hSigOverBkg = (TH1D*)hEfficiency->Clone("hSigOverBkg");
  array<TH1D*, 4> hists = {hPurity, hEfficiency, hSignificance, hSigOverBkg};

  // Competing mass cut
  double competingMassCut = 1e-2; // 10 MeV
  int b = thn2->GetAxis(K0SmassAxis)->FindBin(MassK0S);
  // array<int, 2> kBins = getProjectionBins(thn->GetAxis(K0SmassAxis), MassK0S - competingMassCut, MassK0S + competingMassCut);
  // thn->GetAxis(K0SmassAxis)->SetRange(kBins[0], kBins[1]);
  thn2->GetAxis(K0SmassAxis)->SetRange(b, b);
  cout << "K0S mass bin: " << b << " = " << thn2->GetAxis(K0SmassAxis)->GetBinLowEdge(b) << " GeV" << endl;
  // cout << "K0S mass range: " << kBins[0] << " = " << thn->GetAxis(K0SmassAxis)->GetBinLowEdge(kBins[0]) << " GeV" << endl;
  // cout << "K0S mass range: " << kBins[1] << " = " << thn->GetAxis(K0SmassAxis)->GetBinLowEdge(kBins[1]) << " GeV" << endl;

  // thn2->Scale(-1.);
  // thn->Add(thn2, -1.);

  TH1D* hMass = (TH1D*)thn->Projection(Lambda0massAxis);
  TH1D* hMass2 = (TH1D*)thn2->Projection(Lambda0massAxis);
  hMass->Add(hMass2);
  setStyle(hMass, 0);
  hMass->Draw();
  return;
  // Do reference bin first
  refValues = getPurityAndRelEfficiency(thn, ptBins, refBin, refBin, cutAxis, rebinNumber, hadron, dataSet, fitSettings, nSigma, refValues, hists);
  for (int i = underflowBin; i <= overflowBin; i++) {
    if (i == refBin) { continue; }
    refValues = getPurityAndRelEfficiency(thn, ptBins, i, refBin, cutAxis, rebinNumber, hadron, dataSet, fitSettings, nSigma, refValues, hists);
  }

  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);

  string canvasName = TString::Format("cEP_%s_%s_%.1f-%.1f", hadron.c_str(), axisNames[cutAxis].c_str(), lowpt, highpt).Data();
  TCanvas* cEP   = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1800, 900);
  cEP->cd();
  TLatex*  lData = CreateLatex(0.25, 0.35, dataSet.c_str(), 0.04);
  TLatex*  lPt   = CreateLatex(0.25, 0.3, TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", lowpt, highpt).Data(), 0.04);
  TLegend* lEP   = CreateLegend(0.25, 0.5, 0.15, 0.25, "", 0.04);

  hPurity = hists[0];
  hEfficiency = hists[1];
  hSignificance = hists[2];
  hSigOverBkg = hists[3];
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

  string saveName = hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowpt, highpt).Data();
  saveName += "_" + axisNames[cutAxis];
  saveName += "_EffPur";
  saveName += ".pdf";
  cEP->SaveAs(saveName.c_str());
}
