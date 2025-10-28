//
// This macro contains a set of useful functions for plotting histograms
// and performing standard operations
//

#include <vector>
#include <iostream>
#include <typeinfo>
#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TString.h"

#include "TLegend.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TROOT.h"

#ifndef HISTUTILS_H
#define HISTUTILS_H

namespace histutils {

// GLOBAL SETTINGS !!!
const double MassK0S = 0.497611;
const double MassLambda = 1.115683;
const double MassLambda0 = MassLambda;
gStyle->SetNdivisions(505, "xy");

// -------------------------------------------------------------------------------------------------
//
// Get the name of the dataset corresponding to a train number
//
// -------------------------------------------------------------------------------------------------

string getDataSet(int train) {
  switch (train) {
    case 210373: return "LHC24b1b";
    case 252064: return "LHC22o_pass6";
    case 271952: return "LHC24b1b";
    case 280432: return "LHC24g4";
    case 282430: return "LHC22o_pass7_small";
    case 350079: return "LHC22o_pass7_small";
    case 349872: return "LHC22o_pass7_small";
    case 349871: return "LHC24g4";
    case 419138: return "LHC22o_pass7_small";
    case 417810: return "LHC22o_pass7_small";
    case 419996: return "LHC22o_pass7_medium";
    case 420215: return "LHC22o_pass7_small";
    case 426828: return "LHC22o_pass7";
    case 428560: return "LHC22o_pass7_medium";
    case 436064: return "LHC25a2";
    case 436232: return "LHC22o_pass7_small";
    case 436233: return "LHC22o_pass7_small";
    case 439670: return "LHC25a2b";
    case 439671: return "LHC25a2b";
    default:     return "Could not find dataset";
  }
}

// -------------------------------------------------------------------------------------------------
//
// Get bins corresponding to a range on an axis
// Default behaviour is to include overflow and underflow bins
//
// -------------------------------------------------------------------------------------------------

std::array<int, 2> getProjectionBins(const TAxis* axis, const double min, const double max, double epsilon = 1e-5) {
  int nbins = axis->GetNbins();
  int firstBin = 0, lastBin = nbins + 1;

  if (min > -900)
    firstBin = axis->FindBin(min + epsilon);
  if (max > -900)
    lastBin = axis->FindBin(max - epsilon);

  string warningPrefix = "histutils::getProjectionBins() Warning: ";
  if (firstBin <= 0)
    cout << warningPrefix << "lower bound " << min << " is out of range of axis (" << axis->GetBinLowEdge(1) << ")!\n";

  if (firstBin > nbins)
    cout << warningPrefix << "lower bound " << min << " is out of range of axis (" << axis->GetBinUpEdge(nbins) << ")!\n";

  if (lastBin <= 0)
    cout << warningPrefix << "upper bound " << max << " is out of range of axis (" << axis->GetBinUpEdge(nbins) << ")!" << endl;

  if (lastBin > axis->GetNbins())
    cout << warningPrefix << "upper bound " << max << " is out of range of axis (" << axis->GetBinUpEdge(nbins) << ")!" << endl;

  return std::array{firstBin, lastBin};
}

std::array<double, 2> getProjectionEdges(const TAxis* axis, const std::array<int, 2>& bins) {
  return std::array{axis->GetBinLowEdge(bins[0]), axis->GetBinUpEdge(bins[1])};
}

// -------------------------------------------------------------------------------------------------
//
// Check if a histogram is empty in a given range
//
// -------------------------------------------------------------------------------------------------

// Check if a histogram is empty in a given range
bool isHistEmptyInRange(TH1* h, int low, int high, double threshold = 1e-10) {
  double integral = h->Integral(low, high);
  if (std::isnan(integral))
    return true;
  else
    return (integral < threshold);
}
bool isHistEmptyInRange(TH2* h, int xlow, int xhigh, int ylow, int yhigh, double threshold = 1e-10) {
  return isHistEmptyInRange(h->ProjectionX("px", ylow, yhigh), xlow, xhigh, threshold);
}

// -------------------------------------------------------------------------------------------------
//
// Get histogram bounds (1D only!)
//
// -------------------------------------------------------------------------------------------------

// Return the bin with the most extreme value in a range of the histogram
// Optionally accounts for error
int getExtremeBinInRange(TH1* h, int minBin, int maxBin, bool doError, bool doMin) {
  // Start scale at most extreme value
  int extremeBin = 0;
  double scale = std::numeric_limits<double>::lowest();
  if (doMin)
    scale = std::numeric_limits<double>::max();

  for (int i = minBin; i <= maxBin; i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    double s = bc;
    if (doError)
      doMin ? s -= be : s += be;

    bool isExtremeBin = (doMin && s < scale) || (!doMin && s > scale);
    if (isExtremeBin) {
      scale = s;
      extremeBin = i;
    }
  }
  return extremeBin;
}

// Return the upper or lower of a histogram in a bin range
double getScaleInRange(TH1* h, int minBin, int maxBin, bool doError, bool doMin) {
  int extremeBin = getExtremeBinInRange(h, minBin, maxBin, doError, doMin);
  double binContent = h->GetBinContent(extremeBin);
  double binError = h->GetBinError(extremeBin);

  if (doError)
    return doMin ? binContent - binError : binContent + binError;
  else
    return binContent;
}

// Return the scale of a number of histograms in a bin range
double getScaleInRange(vector<TH1*> v, int minBin, int maxBin, bool doError, bool doMin, bool doSum) {
  double scale = 0;
  for (auto h : v) {
    double s = getScaleInRange(h, minBin, maxBin, doError, doMin);

    if (doSum)
      scale = doMin ? scale - s : scale + s;
    else
      scale = doMin ? std::min(scale, s) : std::max(scale, s);
  }
  return scale;
}

// Lower bounds
int getLowerBoundBinInRange(TH1* h, int minBin, int maxBin, bool doError) {
  return getExtremeBinInRange(h, minBin, maxBin, doError, true);
}
int getLowerBoundBin(TH1* h, bool doError) {
  return getLowerBoundBinInRange(h, 1, h->GetNbinsX(), doError);
}

double getLowerBoundInRange(TH1* h, int minBin, int maxBin, bool doError) {
  return getScaleInRange(h, minBin, maxBin, doError, true);
}
double getLowerBound(TH1* h, bool doError) {
  return getLowerBoundInRange(h, 1, h->GetNbinsX(), doError);
}

double getLowerBoundInRange(vector<TH1*> v, int minBin, int maxBin, bool doError, bool doSum) {
  return getScaleInRange(v, minBin, maxBin, doError, true, doSum);
}
double getLowerBound(vector<TH1*> v, bool doError, bool doSum) {
  return getLowerBoundInRange(v, 1, v[0]->GetNbinsX(), doError, doSum);
}

// Upper bounds
int getUpperBoundBinInRange(TH1* h, int minBin, int maxBin, bool doError) {
  return getExtremeBinInRange(h, minBin, maxBin, doError, false);
}
int getUpperBoundBin(TH1* h, bool doError) {
  return getUpperBoundBinInRange(h, 1, h->GetNbinsX(), doError);
}

double getUpperBoundInRange(TH1* h, int minBin, int maxBin, bool doError) {
  return getScaleInRange(h, minBin, maxBin, doError, false);
}
double getUpperBound(TH1* h, bool doError) {
  return getUpperBoundInRange(h, 1, h->GetNbinsX(), doError);
}

double getUpperBoundInRange(vector<TH1*> v, int minBin, int maxBin, bool doError, bool doSum) {
  return getScaleInRange(v, minBin, maxBin, doError, false, doSum);
}
double getUpperBound(vector<TH1*> v, bool doError, bool doSum) {
  return getUpperBoundInRange(v, 1, v[0]->GetNbinsX(), doError, doSum);
}

// -------------------------------------------------------------------------------------------------
//
// Formatting text for V0s
//
// -------------------------------------------------------------------------------------------------

// Formats the hadron name to look nice (Greek letters, sub- and superscripts)
string formatHadronName(string hadron) {
  string had = hadron;
  if (hadron == "pi"){
    had = "#pi^{#pm}";
  }
  else if (hadron == "pi0"){
    had = "#pi^{0}";
  }
  else if (hadron == "K0S"){
    had = "K^{0}_{S}";
  }
  else if (hadron == "Lambda0" || hadron == "Lambda"){
    // had = "#Lambda^{0}";
    had = "#Lambda";
  }
  else if (hadron == "AntiLambda0" || hadron == "AntiLambda"){
    // had = "#bar{#Lambda}^{0}";
    had = "#bar{#Lambda}";
  }
  return had;
}

// Returns decay products given a hadron
string formatHadronDaughters(string hadron) {
  string daughters = hadron;
  if ("K0S" == hadron) {
    daughters = "#pi^{+}#pi^{-}";
  }
  else if ("Lambda0" == hadron || "Lambda" == hadron) {
    daughters = "p#pi^{-}";
  }
  else if ("AntiLambda0" == hadron || "AntiLambda" == hadron) {
    daughters = "#bar{p}#pi^{+}";
  }
  return daughters;
}

// -------------------------------------------------------------------------------------------------
//
// Normalise 2D histograms
//
// -------------------------------------------------------------------------------------------------

void normaliseHistRowByRow(TH2* hist) {
  int firstColBin = 1, lastColBin = hist->GetNbinsX();
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY();
  for (int iRow = 1; iRow <= lastRowBin; iRow++) {
    double integral = hist->Integral(firstColBin, lastColBin, iRow, iRow);
    if (integral < 1) { continue; }
    for (int iCol = 1; iCol <= lastColBin; iCol++) {
      double binContent = hist->GetBinContent(iCol, iRow);
      binContent /= integral;
      hist->SetBinContent(iCol, iRow, binContent);
    }
  }
}

void normaliseHistColByCol(TH2* hist) {
  int firstColBin = 1, lastColBin = hist->GetNbinsX();
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY();
  for (int iCol = 1; iCol <= lastColBin; iCol++) {
    double integral = hist->Integral(iCol, iCol, firstRowBin, lastRowBin);
    if (integral < 1) { continue; }
    for (int iRow = 1; iRow <= lastRowBin; iRow++) {
      double binContent = hist->GetBinContent(iCol, iRow);
      binContent /= integral;
      hist->SetBinContent(iCol, iRow, binContent);
    }
  }
}

// -------------------------------------------------------------------------------------------------
//
// Rounding for log plots
//
// -------------------------------------------------------------------------------------------------

double roundToNextPowerOfTen(double x) {
  return std::pow(10., std::ceil(std::log10(x)));
}

double roundToPrevPowerOfTen(double x) {
  return std::pow(10., std::floor(std::log10(x)));
}

// -------------------------------------------------------------------------------------------------
//
// Make a subset of a histogram
//
// -------------------------------------------------------------------------------------------------

TH1* makeHistSubset(TH1* data, int minBin, int maxBin, string name = "region") {
  TH1* region = (TH1*)data->Clone(name.c_str());
  region->Reset();
  for (int i = minBin; i <= maxBin; i++) {
    region->SetBinContent(i, data->GetBinContent(i));
    region->SetBinError(i, data->GetBinError(i));
  }
  return region;
}

// -------------------------------------------------------------------------------------------------
//
// Divide two histograms with protection against nan/null bin content
// Assumes uncorrelated errors
// Take caution: bins must match!
//
// -------------------------------------------------------------------------------------------------
TH1* divideWithProtection(TH1* base, TH1* divideBy, double threshold = 1e-25) {
  // TODO: Check if same binning
  TH1* result = (TH1*)base->Clone("result");
  result->Reset();

  for (int i = 0; i <= 1 + base->GetNbinsX(); i++) {
    double numerator = base->GetBinContent(i);
    double numError = base->GetBinError(i);
    double denominator = divideBy->GetBinContent(i);
    double denError = divideBy->GetBinError(i);

    if (std::isnan(numerator) || std::isnan(denominator))
      continue;
    else if (std::abs(numerator) < threshold || std::abs(denominator) < threshold)
      continue;

    double newBinContent = numerator / denominator;
    result->SetBinContent(i, newBinContent);

    if (std::isnan(numError) || std::isnan(denError))
      continue;

    double numRelError = numError / numerator;
    double denRelError = denError / denominator;
    double newBinError = newBinContent * std::sqrt((numRelError * numRelError) + (denRelError * denRelError));
    result->SetBinError(i, newBinError);
  }
  return result;
}

// -------------------------------------------------------------------------------------------------
//
// Make new rebinned histogram and propagates errors via sum of squares
// Take caution: bin edges must line up!
//
// -------------------------------------------------------------------------------------------------

TH1* rebinHist(const TH1* input, const TH1* output) {
  TH1* h = (TH1*)output->Clone(TString::Format("%s_rebinned", input->GetName()).Data());
  h->Reset();
  double newContents[h->GetNbinsX() + 2];
  double newErrorsSquared[h->GetNbinsX() + 2];
  for (int i = 0; i <= input->GetNbinsX()+1; i++) {
    double content = input->GetBinContent(i);
    double centre = input->GetBinCenter(i);
    double error = input->GetBinError(i);
    int newBin = h->FindBin(centre);
    newContents[newBin] += content;
    newErrorsSquared[newBin] += error * error;
  }
  for (int i = 0; i <= h->GetNbinsX()+1; i++) {
    double content = newContents[i];
    double error = std::sqrt(newErrorsSquared[i]);
    h->SetBinContent(i, content);
    h->SetBinError(i, error);
  }
  return h;
}

TH1* rebinnedV0PtHist(string hadron, string name) {
  if (hadron == "K0S") {
    const int nBinsK0S = 11;
    const double edgesK0S[nBinsK0S + 1] = {0., 1., 2., 3., 4., 5., 10., 15., 20., 25., 30., 40.};
    return new TH1D(name.c_str(), name.c_str(), nBinsK0S, edgesK0S);
  } else if (hadron == "Lambda" || hadron == "AntiLambda") {
    const int nBinsLAL = 10;
    const double edgesLAL[nBinsLAL + 1] = {0., 1., 2., 3., 4., 5., 10., 15., 20., 30., 40.};
    return new TH1D(name.c_str(), name.c_str(), nBinsLAL, edgesLAL);
  } else {
    cout << "Hadron not recognised for rebinned V0 pt hist: " << hadron << endl;
    return nullptr;
  }
}
TH1* rebinnedV0ZHist(string name) {
  const int nBins = 10;
  const double edges[nBins + 1] = {0.001, .101, .201, .301, .401, .501, .601, .701, .801, .901, 1.001};
  return new TH1D(name.c_str(), name.c_str(), nBins, edges);
}

// -------------------------------------------------------------------------------------------------
//
// Print the parameter names, values, and limits of a function
//
// -------------------------------------------------------------------------------------------------

void printParLimits(TF1* f) {
  for (int i = 0; i < f->GetNpar(); i++) {
    double min, max;
    f->GetParLimits(i, min, max);
    cout << f->GetName() << " (" << i << " = " << f->GetParName(i) << ") " << f->GetParameter(i) << " (" << min << ", " << max << ")" << endl;
  }
  cout << endl;
}

} // namespace histutils

#endif // HISTUTILS_H
