
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
gStyle->SetNdivisions(505, "xy");

string getDataSet(int train)
{
  if (252064 == train) return "LHC22o_pass6";
  if (282430 == train) return "LHC22o_pass7_small";
  return "Could not find dataset";
}

double getHistLowerBound(TH1* h, bool doError)
{
  if (!doError) { return h->GetBinContent(h->GetMinimumBin()); }

  double scale = 900.;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    scale = min(scale, bc - be);
  }
  return scale;
}
double getHistLowerBound(vector<TH1*> h, bool doError, bool doSum)
{
  double scale = 900.;
  for (auto i : h) {
    if (doSum)
      scale += getHistLowerBound(i, doError);
    else
      scale = min(scale, getHistLowerBound(i, doError));
  }
  return scale;
}
double getHistUpperBound(TH1* h, bool doError)
{
  if (!doError) { return h->GetBinContent(h->GetMaximumBin()); }

  double scale = 0.;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    scale = max(scale, bc + be);
  }
  return scale;
}
double getHistUpperBound(vector<TH1*> h, bool doError, bool doSum)
{
  double scale = -900.;
  for (auto i : h) {
    if (doSum)
      scale += getHistUpperBound(i, doError);
    else
      scale = max(scale, getHistUpperBound(i, doError));
  }
  return scale;
}

array<double, 2> getRange(TH1* h, string hadron, double ptmin, double ptmax, string range)
{
  // TAxis* a = h->GetXaxis();
  vector<double> ptBinEdges = {0., 1., 2., 3., 4., 5., 10., 15., 20.};
  int ptbin = 0;
  for (int i = 0; i < ptBinEdges.size(); i++) {
    if (ptmin < ptBinEdges[i] && ptmax >= ptBinEdges[i]) {
      ptbin = i;
      break;
    }
  }

  if ("K0S" == hadron) {
    if (ptmax < 10.) {
      if (range == "fit") {
        return {0.45, 0.55};
      }
      else if (range == "peak") {
        return {0.48, 0.53};
      }
      else if (range == "rtail") {
        return {0.53, 0.55};
      }
      else if (range == "lhs") {
        return {0.45, 0.5};
      }

    }
    return {0.45, 0.55};
  }
  else if ("Lambda0" == hadron || "AntiLambda0" == hadron) {
    return {0.45, 0.55};
  }
  else {
    cout << "Hadron " << hadron << " not recognized" << endl;
    return {-2., -1.};
  }
}

// Prints the parameter limits of a function
void printParLimits(TF1* f)
{
  for (int i = 0; i < f->GetNpar(); i++) {
    double min, max;
    f->GetParLimits(i, min, max);
    cout << f->GetName() << " (" << i << " = " << f->GetParName(i) << ") " << f->GetParameter(i) << " (" << min << ", " << max << ")" << endl;
  }
  cout << endl;
}
// Sorts a vector of functions alphabetically by name
// Important, because ROOT reshuffles the parameters when combining functions
vector<TF1*> sortAlphabetically(vector<TF1*> v)
{
  vector<TF1*> sorted;
  for (auto f : v) {
    string a = f->GetName();
    bool isLast = true;
    for (int i = 0; i < sorted.size(); i++) {
      TF1* s = sorted[i];
      string b = s->GetName();
      if (a.compare(b) < 0) {
        sorted.insert(sorted.begin() + i, f);
        isLast = false;
        break;
      }
    }
    if (isLast) sorted.push_back(f);
  }
  return sorted;
}
// Creates a new function that is the sum of all input functions
TF1* combineTFs(vector<TF1*> w, bool verbose = false)
{
  string fname = "";
  double xmin = w[0]->GetXmin();
  double xmax = w[0]->GetXmax();
  vector<double> pars;
  vector<vector<double> > vars;
  vector<TF1*> v = sortAlphabetically(w);

  for (int i = 0; i < v.size(); i++) {
    TF1* f = v[i];
    fname += f->GetName();
    if (i < v.size() - 1) {
      fname += "+";
    }

    if (f->GetXmin() < xmin) {
      xmin = f->GetXmin();
    }
    if (f->GetXmax() > xmax) {
      xmax = f->GetXmax();
    }

    for (int j = 0; j < f->GetNpar(); j++) {
      double pmin, pmax;
      f->GetParLimits(j, pmin, pmax);

      pars.push_back(f->GetParameter(j));
      vars.push_back({pmin, pmax});
    }
  }

  TF1* F = new TF1(fname.c_str(), fname.c_str(), xmin, xmax);
  for (int i = 0; i < pars.size(); i++) {
    F->SetParameter(i, pars[i]);
    F->SetParLimits(i, vars[i][0], vars[i][1]);
  }

  if (verbose) {
    cout << "unsorted: " << endl;
    for (auto i : w) { printParLimits(i); }
    cout << endl
        << "sorted: " << endl;
    for (auto i : v) { printParLimits(i); }
    cout << endl
        << "combined: " << endl;
    printParLimits(F);
    cout << endl;
  }
  return F;
}
// Update parameters of a function constructed with combineTFs
// !!! Assumes that the functions are sorted alphabetically !!!
// !!! If this is not true, the function will not work, as the parameters will be reshuffled !!!
void splitTFs(TF1* total, vector<TF1*> functions)
{
  vector<double> pars;
  // vector<int> nPars;
  for (int i = 0; i < total->GetNpar(); i++) {
    pars.push_back(total->GetParameter(i));
  }
  // for (auto f : functions) {
  //   nPars.push_back(f->GetNpar());
  // }
  for (int i = 0; i < pars.size(); i++) {
    int iLocal = i; // "local" index
    for (int j = 0; j < functions.size(); j++) {
      int np = functions[j]->GetNpar();
      if (iLocal < np) {
        functions[j]->SetParameter(iLocal, pars[i]);
        break;
      }
      else {
        iLocal -= np;
      }
    }
  }
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
// Make residual histogram of data and fit
TH1* makeResidual(TH1* data, TF1* fit)
{
  TH1* h = (TH1*)data->Clone("residual");
  TF1* f = (TF1*)fit->Clone("function");
  f->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()); // Make sure function is defined everywhere
  h->Add(f, -1.);
  return h;
}

// Calculates chisq over range, without fitting. Assumes exp error is sqrt(n)
double mychi2(TH1* h, TF1* f, int min, int max)
{
  double chi2 = 0.;
  for (int i = min; i <= max; i++) {
    double v = h->GetBinContent(i);
    double F = f->Eval(h->GetBinCenter(i));
    if (v) chi2 += (F - v) * (F - v) / v;
    // if (v) chi2 += TMath::Sq(f->Eval(h->GetBinCenter(i)) - v) / TMath::Abs(v);
  }
  return chi2;
}
double mychi2(TH1* h, TF1* f, double min, double max)
{
  array<int, 2> bins = getProjectionBins(h->GetXaxis(), min, max);
  return mychi2(h, f, bins[0], bins[1]);
}
// Printer for chi2 debugging
void printChi2(TH1* data, TF1* f)
{
  array<int, 2> bins = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  double NDF = 1 + bins[1] - bins[0] - f->GetNumberFreeParameters();
  cout << "0.45 - 0.55" << endl;
  cout << "chi2/NDF: " << mychi2(data, f, bins[0], bins[1]) / NDF << " = " << mychi2(data, f, bins[0], bins[1]) << " / " << NDF << endl;

  bins = getProjectionBins(data->GetXaxis(), 0.48, 0.525);
  NDF = 1 + bins[1] - bins[0] - f->GetNumberFreeParameters();
  cout << "0.48 - 0.525" << endl;
  cout << "chi2/NDF: " << mychi2(data, f, bins[0], bins[1]) / NDF << " = " << mychi2(data, f, bins[0], bins[1]) << " / " << NDF << endl;
  cout << endl;
}

// -------------------------------------------------------------------------------------------------
//
// Fit functions
//
// -------------------------------------------------------------------------------------------------

// Signal Gaussian, give a range for the peak to initialise the fit
TF1* sigGaus(TH1* h, string uniqueName, array<double, 3> vals = {MassLambda0, 3e-3, 0.01}, array<double, 2> fitRange = {1.08, 1.2})
{
  double massPole   = vals[0];
  double massRegion = vals[1];
  double massRange  = vals[2];
  array<int, 2> peakRegion = getProjectionBins(h->GetXaxis(), massPole - massRegion, massPole + massRegion);
  TH1* k = makeHistSubset(h, peakRegion[0], peakRegion[1]);
  double peakBin = k->GetMaximumBin();
  double peakVal = k->GetBinContent(peakBin);

  TF1* f = new TF1(uniqueName.c_str(), "gaus(0)", fitRange[0], fitRange[1]);
  f->SetParameter(0, peakVal);
  f->SetParameter(1, massPole);
  f->SetParameter(2, 0.5 * massRange);

  f->SetParLimits(0, 0., 1.5 * peakVal);
  f->SetParLimits(1, massPole - massRegion, massPole + massRegion);
  f->SetParLimits(2, 1e-6, massRange);

  for (int i = 0; i < f->GetNpar(); i++) {
    string parName = uniqueName + "_" + to_string(i);
    f->SetParName(i, parName.c_str());
  }
  return f;
}

// Gaussian background
TF1* gaussian(TH1* h, string uniqueName, array<double, 3> vals = {1e3, MassLambda0, 0.01}, array<double, 2> fitRange = {1.08, 1.2})
{
  double amp   = vals[0];
  double mean  = vals[1];
  double sigma = vals[2];

  TF1* f = new TF1(uniqueName.c_str(), "gaus(0)", fitRange[0], fitRange[1]);
  f->SetParameter(0, amp);
  f->SetParameter(1, mean);
  f->SetParameter(2, sigma);
  f->SetParLimits(0, 0.5 * amp, 2. * amp);
  f->SetParLimits(1, mean - 0.1, mean + 0.1);
  f->SetParLimits(2, 0.1 * sigma, 3. * sigma);

  for (int i = 0; i < f->GetNpar(); i++) {
    string parName = uniqueName + "_" + to_string(i);
    f->SetParName(i, parName.c_str());
  }
  return f;
}

// Linear background, give it two x values to which it will initialize the fit
TF1* linear(TH1* h, string uniqueName, array<double, 2> vals, array<double, 2> fitRange = {1.08, 1.2})
{
  TH1* k = (TH1*)h->Clone("k");
  array<int, 2> region = getProjectionBins(k->GetXaxis(), vals[0], vals[1]);
  double dx = k->GetBinCenter(region[1]) - k->GetBinCenter(region[0]);
  double dy = k->GetBinContent(region[1]) - k->GetBinContent(region[0]);

  double b = dy / dx;
  double a = k->GetBinContent(region[0]) - b * k->GetBinCenter(region[0]);

  // cout << "fit range: " << fitRange[0] << " - " << fitRange[1] << endl;
  // cout << "values: " << vals[0] << " - " << vals[1] << endl;
  // cout << "dx: " << dx << "(" << k->GetBinCenter() << ")" << endl;

  TF1* f = new TF1(uniqueName.c_str(), "pol1(0)", fitRange[0], fitRange[1]);
  f->SetParameter(0, a);
  f->SetParameter(1, b);
  f->SetParLimits(0, -2. * (abs(a) + 1), 2. * (abs(a) + 1));
  f->SetParLimits(1, -2. * (abs(b) + 1), 2. * (abs(b) + 1));

  for (int i = 0; i < f->GetNpar(); i++) {
    string parName = uniqueName + "_" + to_string(i);
    f->SetParName(i, parName.c_str());
  }
  return f;
}

// Sigmoid background
TF1* sigmoid(TH1* h, string uniqueName, array<double, 3> vals = {1., 200., 1.1}, array<double, 2> fitRange = {1.08, 1.2})
{
  double amp     = vals[0]; // Amplitude
  double slope   = vals[1]; // 4 * max slope; 4 * slope at x = halfway
  double halfway = vals[2]; // x value at which the function is halfway between 0 and amp

  TF1* f = new TF1(uniqueName.c_str(), "[0] / (1 + exp(-[1] * (x - [2])))", fitRange[0], fitRange[1]);
  f->SetParameter(0, amp);
  f->SetParameter(1, slope);
  f->SetParameter(2, halfway);
  f->SetParLimits(0, 0.1 * amp, 2. * amp);
  f->SetParLimits(1, 1e-2 * slope, 20. * slope);
  f->SetParLimits(2, fitRange[0], fitRange[1]);

  for (int i = 0; i < f->GetNpar(); i++) {
    string parName = uniqueName + "_" + to_string(i);
    f->SetParName(i, parName.c_str());
  }
  return f;
}

// Dilog background
TF1* dilog(TH1* h, string uniqueName, array<double, 3> vals = {2e3, 35., 1.09}, array<double, 2> fitRange = {1.08, 1.2})
{
  double amp   = vals[0]; // Amplitude
  double range = vals[1]; // Upper range of integral in dilog
  double zero  = vals[2]; // x value at which the function is zero

  TF1* f = new TF1(uniqueName.c_str(), "[0]*TMath::DiLog([1] * (x - [2]))", fitRange[0], fitRange[1]);
  f->SetParameter(0, amp);
  f->SetParameter(1, range);
  f->SetParameter(2, zero);
  f->SetParLimits(0, 0.1 * amp, 2. * amp);
  f->SetParLimits(1, 0., 2. * range);
  f->SetParLimits(2, zero - 0.05, zero + 0.05);

  for (int i = 0; i < f->GetNpar(); i++) {
    string parName = uniqueName + "_" + to_string(i);
    f->SetParName(i, parName.c_str());
  }
  return f;
}

// -------------------------------------------------------------------------------------------------
//
// Functions to use with setup
//
// -------------------------------------------------------------------------------------------------

// My "favourite" histogram
// TH1* getHist(double ptmin, double ptmax, string hadron)
TH1* getHist(double ptmin, double ptmax, string hadron, int train)
{
  // string inName = "~/cernbox/TrainOutput/252064/AnalysisResults.root";
  string inName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string histName = "jet-fragmentation/data/V0/V0CutVariation";

  int ptAxis = 0;
  int mAxis;
  if ("K0S" == hadron) {
    mAxis = 1;
  }
  else if ("Lambda0" == hadron) {
    mAxis = 2;
  }
  else if ("AntiLambda0" == hadron) {
    mAxis = 3;
  }
  else {
    cout << "Hadron " << hadron << " not recognized" << endl;
    return nullptr;
  }

  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*) inFile->Get(histName.c_str());

  array<int, 2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  return (thn->Projection(mAxis));
}
TH1* getHist(double ptmin, double ptmax, string hadron, string inName)
{
  string histName = "jet-fragmentation/data/V0/V0CutVariation";

  int ptAxis = 0;
  int mAxis;
  if ("K0S" == hadron) {
    mAxis = 1;
  }
  else if ("Lambda0" == hadron) {
    mAxis = 2;
  }
  else if ("AntiLambda0" == hadron) {
    mAxis = 3;
  }
  else {
    cout << "Hadron " << hadron << " not recognized" << endl;
    return nullptr;
  }

  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*) inFile->Get(histName.c_str());

  array<int, 2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  TH1* hist = thn->Projection(mAxis);
  hist->Sumw2();

  string name = TString::Format("pt%.1f-%.1f", ptmin, ptmax).Data();
  hist->SetName(name.c_str());
  return (hist);
}

void plotFitParts(TCanvas* canvas, TH1F* frame, TH1* h, TF1* fit, vector<TF1*> fitParts, TLegend* legend, vector<TLatex*> latex)
{
  frame->Draw();
  h->Draw("same");
  for (auto f : fitParts) {
    f->Draw("same");
  }
  legend->Draw("same");
  for (auto l : latex) {
    l->Draw("same");
  }
  fit->Draw("same");
  canvas->SaveAs(canvas->GetName());
}
void plotFitParts(TCanvas* canvas, TH1* h, TF1* fit, vector<TF1*> fitParts, TLegend* legend, vector<TLatex*> latex)
{
  TH1F* frame = DrawFrame(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                          0, 1.1 * h->GetBinContent(h->GetMaximumBin()),
                          h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle());
  frame->SetTitle(h->GetTitle());
  plotFitParts(canvas, frame, h, fit, fitParts, legend, latex);
}

TF1* getFitShape(string hadron, double min, double max, int pol = 1)
{
  if ("K0S" == hadron) {
    if (pol == 1)
      return new TF1("fSLP", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4]) + max(0., [5]*(x-[6]) * exp(-1. * (x-[6]) / [7]))", min, max);
    if (pol == 2)
      return new TF1("fSQP", "[0]+[1]*x+[2]*x*x + [3]*TMath::Gaus(x,[4],[5]) + max(0., [6]*(x-[7]) * exp(-1. * (x-[7]) / [8]))", min, max);
    else
      return nullptr;
  }
  else {
    cout << "Hadron " << hadron << " not recognized" << endl;
    return nullptr;
  }
}
// Can we set fit parameters in function?
void setFitParameters(TF1* f, string hadron, double ptmin, double ptmax, int pol = 1)
{
  if (pol < 1 || pol > 2) {
    cout << "Pol degree not recognised! Aborting" << endl;
    return;
  }
  string fName = f->GetName();
  for (int i = 0; i < f->GetNpar(); i++) {
    f->SetParName(i, (fName + "_" + to_string(i)).c_str());
  }

  if ("K0S" == hadron) {
    double mass = MassK0S;
    double signalWidth = 1e-2;
    double peakVal;

    double pOffset, pSlope, pQuad;
    double eAmplitude, eZero, eZeroToMax;

    if (ptmin > 3.9 && ptmax < 5.1) {
      pOffset = 186e3; pSlope = 183e3;
      peakVal = 4.9e6;
      eAmplitude = 1.9e8; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 2. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 1e3 * pSlope);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, 0., 10. * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 1e-3, 0.5);
    } else if (ptmin > 4.9 && ptmax < 10.1) {
      pOffset = 20e3; pSlope = -1e3;
      peakVal = 2.80851e+06;
      eAmplitude = 100e6; eZero = 0.5; eZeroToMax = 6e-3;
      if (1 == pol) {
        f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 15. * pOffset);
        f->SetParameter(1, pSlope);       f->SetParLimits(1, 1e3 * pSlope, 0.);
        f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
        f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
        f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
        f->SetParameter(5, eAmplitude);   f->SetParLimits(5, eAmplitude / 10., 2. * eAmplitude);
        f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
        f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 5e-3, 0.5);
      }
      if (2 == pol) {
        pOffset = -6.8e6; pSlope = 2.8e7; pQuad = -2.8e7;
        peakVal = 2.5e6; signalWidth = 7e-3;
        eAmplitude = 1.3e8;
        f->SetParameter(0, pOffset);      f->SetParLimits(0, 15. * pOffset, 0.);
        f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 1e3 * pSlope);
        f->SetParameter(2, pQuad);        f->SetParLimits(2, 1e3 * pQuad, 0.);
        f->SetParameter(3, peakVal);      f->SetParLimits(3, 0.8 * peakVal, 2 * peakVal);
        f->SetParameter(4, mass);         f->SetParLimits(4, 0.48, 0.51);
        f->SetParameter(5, signalWidth);  f->SetParLimits(5, 1e-3, 2e-2);
        f->SetParameter(6, eAmplitude);   f->SetParLimits(6, eAmplitude / 10., 2. * eAmplitude);
        f->SetParameter(7, eZero);        f->SetParLimits(7, 0.5, 0.55);
        f->SetParameter(8, eZeroToMax);   f->SetParLimits(8, 1e-3, 0.5);

      }
    } else if (ptmin > 9.9 && ptmax < 15.1) {
      pOffset = -7e3; pSlope = 22e3;
      peakVal = 70e3;
      eAmplitude = 2e6; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, -10e3, 0.);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 30e3);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, eAmplitude / 2., 2.5 * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 1e-3, 0.1);
    } else {
      cout << "Cannot determine fit parameters for " << hadron << " in pt range " << ptmin << " - " << ptmax << endl;
    }
  }
}
void setFitParametersPol1GausXex(TF1* f, string hadron, double ptmin, double ptmax)
{
  string fName = f->GetName();
  for (int i = 0; i < f->GetNpar(); i++) {
    f->SetParName(i, (fName + "_" + to_string(i)).c_str());
  }

  if ("K0S" == hadron) {
    double mass = MassK0S;
    double signalWidth = 1e-2;
    double peakVal;

    double pOffset, pSlope, pQuad;
    double eAmplitude, eZero, eZeroToMax;

    if (ptmin > 3.9 && ptmax < 5.1) {
      pOffset = 186e3; pSlope = 183e3;
      peakVal = 4.9e6;
      eAmplitude = 1.9e8; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 2. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 1e3 * pSlope);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, 0., 10. * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 1e-3, 0.5);
    } else if (ptmin > 4.9 && ptmax < 10.1) {
      pOffset = 20e3; pSlope = -1e3;
      peakVal = 2.80851e+06;
      eAmplitude = 100e6; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 15. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 1e3 * pSlope, 0.);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, eAmplitude / 10., 2. * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 5e-3, 0.5);
    } else if (ptmin > 9.9 && ptmax < 15.1) {
      pOffset = -7e3; pSlope = 22e3;
      peakVal = 70e3;
      eAmplitude = 2e6; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, -10e3, 0.);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 30e3);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, eAmplitude / 2., 2.5 * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 1e-3, 0.1);
    } else {
      cout << "Cannot determine fit parameters for " << hadron << " in pt range " << ptmin << " - " << ptmax << endl;
    }
  }
}
void setFitParametersPol2GausXex(TF1* f, string hadron, double ptmin, double ptmax)
{
  string fName = f->GetName();
  for (int i = 0; i < f->GetNpar(); i++) {
    f->SetParName(i, (fName + "_" + to_string(i)).c_str());
  }

  if ("K0S" == hadron) {
    double mass = MassK0S;
    double signalWidth = 1e-2;
    double peakVal;

    double pOffset, pSlope, pQuad;
    double eAmplitude, eZero, eZeroToMax;

    if (ptmin > 3.9 && ptmax < 5.1) {
      // FIXME:
    } else if (ptmin > 4.9 && ptmax < 10.1) {
      pOffset = -6.8e6; pSlope = 2.8e7; pQuad = -2.8e7;
      peakVal = 2.5e6; signalWidth = 7e-3;
      eAmplitude = 1.3e8;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 15. * pOffset, 0.);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 1e3 * pSlope);
      f->SetParameter(2, pQuad);        f->SetParLimits(2, 1e3 * pQuad, 0.);
      f->SetParameter(3, peakVal);      f->SetParLimits(3, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(4, mass);         f->SetParLimits(4, 0.48, 0.51);
      f->SetParameter(5, signalWidth);  f->SetParLimits(5, 1e-3, 2e-2);
      f->SetParameter(6, eAmplitude);   f->SetParLimits(6, eAmplitude / 10., 2. * eAmplitude);
      f->SetParameter(7, eZero);        f->SetParLimits(7, 0.5, 0.55);
      f->SetParameter(8, eZeroToMax);   f->SetParLimits(8, 1e-3, 0.5);
    // } else if (ptmin > 9.9 && ptmax < 15.1) {
      // FIXME:
    } else {
      cout << "Cannot determine fit parameters for " << hadron << " in pt range " << ptmin << " - " << ptmax << endl;
    }
  }
}
void setFitParametersPol1GausGaus(TF1* f, string hadron, double ptmin, double ptmax)
{
  string fName = f->GetName();
  for (int i = 0; i < f->GetNpar(); i++) {
    f->SetParName(i, (fName + "_" + to_string(i)).c_str());
  }

    if ("K0S" == hadron) {
    double mass = MassK0S;
    double signalWidth = 1e-2;
    double peakVal;

    double pOffset, pSlope;
    double gAmp, gWidth;

    if (ptmin > 3.9 && ptmax < 5.1) {
      pOffset = 186e3; pSlope = 183e3;
      peakVal = 4.9e6;
      gAmp = 1.9e8; gWidth = 2e-2;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 2. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 1e3 * pSlope);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, gAmp);         f->SetParLimits(5, 0., 10. * gAmp);
      f->SetParameter(6, gWidth);        f->SetParLimits(6, 1e-2, 0.1);
    } else if (ptmin > 4.9 && ptmax < 10.1) {
      pOffset = 20e3; pSlope = -1e3;
      peakVal = 2.80851e+06;
      gAmp = 100e6; gWidth = 2e-2;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 15. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 1e3 * pSlope, 0.);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, gAmp);         f->SetParLimits(5, gAmp / 10., 2. * gAmp);
      f->SetParameter(6, gWidth);        f->SetParLimits(6, 1e-2, 0.1);
    } else if (ptmin > 9.9 && ptmax < 15.1) {
      pOffset = -7e3; pSlope = 22e3;
      peakVal = 70e3;
      gAmp = 2e6; gWidth = 2e-2;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, -10e3, 0.);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 30e3);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, gAmp);         f->SetParLimits(5, gAmp / 2., 2.5 * gAmp);
      f->SetParameter(6, gWidth);        f->SetParLimits(6, 1e-2, 0.1);
    } else {
      cout << "Cannot determine fit parameters for " << hadron << " in pt range " << ptmin << " - " << ptmax << endl;
    }
  }
}
void setFitParametersPol1GausExp(TF1* f, string hadron, double ptmin, double ptmax)
{
  string fName = f->GetName();
  for (int i = 0; i < f->GetNpar(); i++) {
    f->SetParName(i, (fName + "_" + to_string(i)).c_str());
  }

  if ("K0S" == hadron) {
    double mass = MassK0S;
    double signalWidth = 1e-2;
    double peakVal;

    double pOffset, pSlope, pQuad;
    double eAmplitude, eZero, eZeroToMax;

    if (ptmin > 3.9 && ptmax < 5.1) {
      pOffset = 186e3; pSlope = 183e3;
      peakVal = 4.9e6;
      eAmplitude = 1.9e8; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 2. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 1e3 * pSlope);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, 0., 10. * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 1e-3, 0.5);
    } else if (ptmin > 4.9 && ptmax < 10.1) {
      pOffset = 20e3; pSlope = -1e3;
      peakVal = 2.80851e+06;
      eAmplitude = 100e6; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 15. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 1e3 * pSlope, 0.);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, eAmplitude / 10., 2. * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 5e-3, 0.5);
    } else if (ptmin > 9.9 && ptmax < 15.1) {
      pOffset = -7e3; pSlope = 22e3;
      peakVal = 70e3;
      eAmplitude = 2e6; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, -10e3, 0.);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 30e3);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, eAmplitude / 2., 2.5 * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 1e-3, 0.1);
    } else {
      cout << "Cannot determine fit parameters for " << hadron << " in pt range " << ptmin << " - " << ptmax << endl;
    }
  }
}
void setFitParametersPol2GausExp(TF1* f, string hadron, double ptmin, double ptmax)
{
  string fName = f->GetName();
  for (int i = 0; i < f->GetNpar(); i++) {
    f->SetParName(i, (fName + "_" + to_string(i)).c_str());
  }

  if ("K0S" == hadron) {
    double mass = MassK0S;
    double signalWidth = 1e-2;
    double peakVal;

    double pOffset, pSlope, pQuad;
    double eAmplitude, eZero, eZeroToMax;

    if (ptmin > 3.9 && ptmax < 5.1) {
      pOffset = 186e3; pSlope = 183e3;
      peakVal = 4.9e6;
      eAmplitude = 1.9e8; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 2. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 1e3 * pSlope);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, 0., 10. * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 1e-3, 0.5);
    } else if (ptmin > 4.9 && ptmax < 10.1) {
      pOffset = 20e3; pSlope = -1e3;
      peakVal = 2.80851e+06;
      eAmplitude = 100e6; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 15. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 1e3 * pSlope, 0.);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, eAmplitude / 10., 2. * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 5e-3, 0.5);
    } else if (ptmin > 9.9 && ptmax < 15.1) {
      pOffset = -7e3; pSlope = 22e3;
      peakVal = 70e3;
      eAmplitude = 2e6; eZero = 0.5; eZeroToMax = 6e-3;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, -10e3, 0.);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 0., 30e3);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, eAmplitude);   f->SetParLimits(5, eAmplitude / 2., 2.5 * eAmplitude);
      f->SetParameter(6, eZero);        f->SetParLimits(6, 0.5, 0.55);
      f->SetParameter(7, eZeroToMax);   f->SetParLimits(7, 1e-3, 0.1);
    } else {
      cout << "Cannot determine fit parameters for " << hadron << " in pt range " << ptmin << " - " << ptmax << endl;
    }
  }
}
// Input to add: rebin
// Plots multiple fits for signal+bkg
void plotBkgs(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  data->Sumw2();
  setStyle(data, 0);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  TCanvas* fitcanvas = new TCanvas((saveName + ".pdf").c_str(), (saveName + ".pdf").c_str(), 1800, 900);
  TCanvas* rescanvas = new TCanvas((saveName + "_residuals.pdf").c_str(), (saveName + "_residuals.pdf").c_str(), 1800, 900);
  vector<TF1*> fits;
  vector<TH1*> residuals;

  double fitmin = 0.45, fitmax = 0.55;
  TF1* fSLP = getFitShape(hadron, fitmin, fitmax);
  setFitParameters(fSLP, hadron, ptmin, ptmax);
  fits.push_back(fSLP);
  printParLimits(fSLP);
  // cout << fSLP->GetName() << " = " << fSLP->GetExpFormula() << endl;

  TLegend* legend = CreateLegend(0.25, 0.8, 0.7, 0.9);
  for (int i = 0; i < fits.size(); i++) {
    fits[i]->SetLineWidth(3);
    fits[i]->SetLineColor(GetColor(i+1));
    data->Fit(fits[i], "RSBQ0");
    string newName = TString::Format("%s (#chi^{2}/NDF = %.1f)", fits[i]->GetName(), fits[i]->GetChisquare() / fits[i]->GetNDF()).Data();
    fits[i]->SetName(newName.c_str());
    legend->AddEntry(fits[i], newName.c_str(), "l");
    printParLimits(fits[i]);
    // printChi2(data, fits[i]);
  }

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", ptmin, ptmax).Data();
  string xTitle = "#it{M}(" + formatHadronDaughters(hadron) + ") (GeV/#it{c}^{2})";
  string yTitle = "counts";
  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();

  fitcanvas->cd();
  double yMinFrame = 0., yMaxFrame = 1.5 * getHistScale(data, false);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());
  frame->Draw();
  legend->Draw("same");
  data->Draw("same");
  for (auto f : fits) {
    f->SetRange(xMinFrame, xMaxFrame);
    f->Draw("same");
    // printChi2(data, f);
  }
  fitcanvas->SaveAs(fitcanvas->GetName());

  rescanvas->cd();
  for (auto f : fits) {
    TH1* residual = makeResidual(data, f);
    residual->SetLineWidth(3);
    residual->SetLineColor(f->GetLineColor());
    residual->SetMarkerColor(f->GetLineColor());
    residuals.push_back(residual);
    // printChi2(data, f);
  }

  yMinFrame = 1.1 * getHistLowerBound(residuals, false, false);
  yMaxFrame = 1.1 * getHistScale(residuals, false, false);
  TH1F* resframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "residuals");
  resframe->SetTitle((dataSet + ", " + ptText).c_str());
  resframe->Draw();
  for (auto r : residuals) {
    r->Draw("same");
  }
  rescanvas->SaveAs(rescanvas->GetName());
}

// Test whether to fit in one go or in parts
void testFittingOrder()
{
  string inName  = "~/cernbox/TrainOutput/252064/AnalysisResults.root";
  string dataSet = "LHC22o_pass6";
  string hadron  = "K0S";
  double ptmin = 5., ptmax = 10.;

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  data->Sumw2();
  setStyle(data, 0);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_residuals";
  saveName += ".pdf";
  TCanvas* residualCanvas = new TCanvas(saveName.c_str(), saveName.c_str(), 1800, 2000);
  residualCanvas->Divide(1,2, 0.01, 0.01);

  TLegend* fitLegend = CreateLegend(0.2, 0.5, 0.65, 0.88);
  TLegend* residualLegend = CreateLegend(0.7, 0.9, 0.65, 0.88, "Residuals:");

  // Frame settings
  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", ptmin, ptmax).Data();
  string xTitle = "#it{M}(" + formatHadronDaughters(hadron) + ") (GeV/#it{c}^{2})";
  string yTitle = "counts";
  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistScale(data, false);

  // Fit settings
  double mass = MassK0S;
  double dM = 5e-2;
  double signalWidth = 1e-2;
  double peakVal = 2.80851e+06;
  // Functions to combine
  // array<int, 2> peakRegion = getProjectionBins(hist->GetXaxis(), mass - dM, mass + dM);
  // TH1* k = makeHistSubset(hist, peakRegion[0], peakRegion[1]);
  // double peakVal = k->GetBinContent(k->GetMaximumBin());

  // array<double, 2> sigRange = getRange(hist, hadron, ptmin, ptmax, "peak");
  array<double,2> sigRange = {0.48, 0.5};

  TF1* fSL = new TF1("fSL", "[0] + [1]*x + [2]*TMath::Gaus(x, [3], [4])", 0.48, 0.5);
  fSL->SetParameter(0, 200e3);       fSL->SetParLimits(0, 0., 300e3);                  fSL->SetParName(0, "fSL_0");
  fSL->SetParameter(1, -1e3);        fSL->SetParLimits(1, -1e6, 0.);                   fSL->SetParName(1, "fSL_1");
  fSL->SetParameter(2, peakVal);     fSL->SetParLimits(2, 0.8 * peakVal, 2 * peakVal); fSL->SetParName(2, "fSL_2");
  fSL->SetParameter(3, mass);        fSL->SetParLimits(3, 0.48, 0.51);                 fSL->SetParName(3, "fSL_3");
  fSL->SetParameter(4, signalWidth); fSL->SetParLimits(4, 1e-3, 2e-2);                 fSL->SetParName(4, "fSL_4");
  // printParLimits(fSL);
  fSL->SetLineWidth(3); fSL->SetLineColor(GetColor(1)); data->Fit(fSL, "RSBQ0");
  // printParLimits(fSL);
  fitLegend->AddEntry(fSL, "pol1+G(x)", "l");

  // Fit pol1 first. Then fit Gauss with fixed pol1
  double x0 = 0.45, x1 = 0.55;
  double y0 = data->GetBinContent(data->FindBin(x0+1e-3)); double y1 = data->GetBinContent(data->FindBin(x1-1e-3));
  TF1* l = new TF1("l", "pol1", x0, x1);
  double b = (y1 - y0) / (x1 - x0), a = (y1 + y0 - b * (x1 + x0)) / 2;
  l->SetParameter(1, b); l->SetParLimits(1, 100.*b, 0.); l->SetParName(1, "l_1");
  l->SetParameter(0, a); l->SetParLimits(0, 0., 100.*a); l->SetParName(0, "l_0");
  // data->Fit(l, "RSBQ0");
  l->SetLineColor(GetColor(2));
  fitLegend->AddEntry(l, "f(x)=pol1(x)", "l");
  // printParLimits(l);

  // Fix pol1 parameters and fit Gauss
  TF1* g = new TF1("g", "[0] + [1]*x + [2]*TMath::Gaus(x, [3], [4])", 0.48, 0.5);
  g->FixParameter(0, l->GetParameter(0));
  g->FixParameter(1, l->GetParameter(1));
  g->SetParameter(2, peakVal);            g->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
  g->SetParameter(3, mass);               g->SetParLimits(3, 0.48, 0.51);
  g->SetParameter(4, signalWidth);        g->SetParLimits(4, 1e-3, 2e-2);
  g->SetLineWidth(3); g->SetLineColor(GetColor(3)); data->Fit(g, "RSBQ0");
  // printParLimits(g);
  fitLegend->AddEntry(g, "f + G(x)", "l");

  // Residuals
  TH1* residual_fSL = makeResidual(data, fSL);
  residual_fSL->SetName("residual_p1+G");
  setStyle(residual_fSL, 1);

  TF1* fP = new TF1("fP", "max([0] * (x-[1]) * exp(-1. * (x-[1]) / [2]) - [3], 0.)", 0.505, 0.525);
  fP->SetParameter(0, 200e3 * 500); fP->SetParLimits(0, 100e5, 200e6); fP->SetParName(0, "fP_0");
  fP->SetParameter(1, 0.5);   fP->SetParLimits(1, 0.45, 0.55); fP->SetParName(1, "fP_1");
  fP->SetParameter(2, 6e-3);  fP->SetParLimits(2, 5e-3, 0.5);  fP->SetParName(2, "fP_2");
  fP->SetParameter(3, 200e3); fP->SetParLimits(3, 0., 300e3);  fP->SetParName(3, "fP_3");
  fP->SetLineWidth(3); fP->SetLineColor(GetColor(1)); residual_fSL->Fit(fP, "RSBQ0");
  residualLegend->AddEntry(fP, "Res(fSL): x e^{-x/b} -d", "l");
  // residualLegend->AddEntry(fP, "max(0, -d + a(x-b) e^{-#frac{(x-b)}{c}})", "l");

  TH1* residual_g = makeResidual(data, g);
  residual_g->SetName("residual_g");
  setStyle(residual_g, 3);

  TF1* fP2 = new TF1("fP2", "max([0] * (x-[1]) * exp(-1. * (x-[1]) / [2]), 0.)", 0.505, 0.525);
  fP2->SetParameter(0, 200e3 * 500); fP2->SetParLimits(0, 100e5, 200e6); fP2->SetParName(0, "fP2_0");
  fP2->SetParameter(1, 0.5);   fP2->SetParLimits(1, 0.45, 0.55); fP2->SetParName(1, "fP2_1");
  fP2->SetParameter(2, 6e-3);  fP2->SetParLimits(2, 5e-3, 0.5);  fP2->SetParName(2, "fP2_2");
  fP2->SetLineWidth(3); fP2->SetLineColor(GetColor(3)); residual_g->Fit(fP2, "RSBQ0");
  residualLegend->AddEntry(fP2, "Res(g): x e^{-x/b}", "l");

  // Plotting
  string frameTitle = (dataSet + ", " + ptText).c_str();
  residualCanvas->cd(1);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, "", "counts");
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->Draw();
  fitLegend->Draw("same");
  // residualLegend->Draw("same");
  data->Draw("same");
  fSL->Draw("same"); fSL->SetRange(xMinFrame, xMaxFrame);
  l->Draw("same");   l->SetRange(xMinFrame, xMaxFrame);
  g->Draw("same");   g->SetRange(xMinFrame, xMaxFrame);

  residualCanvas->cd(2);
  xMinFrame = 0.48; xMaxFrame = 0.55;
  yMinFrame = -200e3, yMaxFrame = 250e3;
  yMaxFrame = 1.1 * getHistScale({residual_fSL, residual_g}, false, false);
  TH1F* resframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "residual");

  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.05);
  resframe->Draw();
  resframe->GetYaxis()->SetTitleOffset(1.0);

  residual_fSL->Draw("same");
  residual_g->Draw("same");
  fP->Draw("same");  fP->SetRange(xMinFrame, xMaxFrame);
  fP2->Draw("same"); fP2->SetRange(xMinFrame, xMaxFrame);

  // printParLimits(fP);
  // printParLimits(fP2);

  residualCanvas->cd(0);
  residualCanvas->SaveAs(residualCanvas->GetName());

  // ------------------------------
  // Combined function fits
  // ------------------------------
  saveName  = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += ".pdf";
  TCanvas* fitcanvas = new TCanvas(saveName.c_str(), saveName.c_str(), 1800, 2000);
  fitcanvas->Divide(1, 2, 0.01, 0.01);
  TLegend* fulllegend = CreateLegend(0.175, 0.5, 0.55, 0.85, "pol(1) + G(x) + x e^{-x/b}");
  TLegend* fullreslegend = CreateLegend(0.175, 0.5, 0.55, 0.85, "#chi^{2}/NDF");

  double fitmin = 0.48, fitmax = 0.525;
  // pol1 + G(x) + xe-x/b
  TF1* fSLP = new TF1("fSLP", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4]) + max(0., [5]*(x-[6]) * exp(-1. * (x-[6]) / [7]))", fitmin, fitmax);
  fSLP->SetParameter(0, 20e3);        fSLP->SetParLimits(0, 0., 300e3);                  fSLP->SetParName(0, "fSLP_0");
  fSLP->SetParameter(1, -1e3);        fSLP->SetParLimits(1, -1e6, 0.);                   fSLP->SetParName(1, "fSLP_1");
  fSLP->SetParameter(2, peakVal);     fSLP->SetParLimits(2, 0.8 * peakVal, 2 * peakVal); fSLP->SetParName(2, "fSLP_2");
  fSLP->SetParameter(3, mass);        fSLP->SetParLimits(3, 0.48, 0.51);                 fSLP->SetParName(3, "fSLP_3");
  fSLP->SetParameter(4, signalWidth); fSLP->SetParLimits(4, 1e-3, 2e-2);                 fSLP->SetParName(4, "fSLP_4");
  fSLP->SetParameter(5, 200e3 * 500); fSLP->SetParLimits(5, 100e5, 200e6);               fSLP->SetParName(5, "fSLP_5");
  fSLP->SetParameter(6, 0.5);         fSLP->SetParLimits(6, 0.45, 0.55);                 fSLP->SetParName(6, "fSLP_6");
  fSLP->SetParameter(7, 6e-3);        fSLP->SetParLimits(7, 5e-3, 0.5);                  fSLP->SetParName(7, "fSLP_7");
  fSLP->SetLineWidth(3); fSLP->SetLineColor(GetColor(1)); data->Fit(fSLP, "RSBQ0");
  fulllegend->AddEntry(fSLP, "SLP: Fit 0.48 - 0.525", "l");
  fullreslegend->AddEntry(fSLP, TString::Format("%.2g = %.2g / %d", fSLP->GetChisquare() / fSLP->GetNDF(), fSLP->GetChisquare(), fSLP->GetNDF()).Data(), "l");

  // Wider fit range
  fitmin = 0.45, fitmax = 0.55;
  TF1* fSLP2 = new TF1("fSLP2", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4]) + max(0., [5]*(x-[6]) * exp(-1. * (x-[6]) / [7]))", fitmin, fitmax);
  fSLP2->SetParameter(0, 20e3);        fSLP2->SetParLimits(0, 0., 300e3);                  fSLP2->SetParName(0, "fSLP2_0");
  fSLP2->SetParameter(1, -1e3);        fSLP2->SetParLimits(1, -1e6, 0.);                   fSLP2->SetParName(1, "fSLP2_1");
  fSLP2->SetParameter(2, peakVal);     fSLP2->SetParLimits(2, 0.8 * peakVal, 2 * peakVal); fSLP2->SetParName(2, "fSLP2_2");
  fSLP2->SetParameter(3, mass);        fSLP2->SetParLimits(3, 0.48, 0.51);                 fSLP2->SetParName(3, "fSLP2_3");
  fSLP2->SetParameter(4, signalWidth); fSLP2->SetParLimits(4, 1e-3, 2e-2);                 fSLP2->SetParName(4, "fSLP2_4");
  fSLP2->SetParameter(5, 200e3 * 500); fSLP2->SetParLimits(5, 100e5, 200e6);               fSLP2->SetParName(5, "fSLP2_5");
  fSLP2->SetParameter(6, 0.5);         fSLP2->SetParLimits(6, 0.45, 0.55);                 fSLP2->SetParName(6, "fSLP2_6");
  fSLP2->SetParameter(7, 6e-3);        fSLP2->SetParLimits(7, 5e-3, 0.5);                  fSLP2->SetParName(7, "fSLP2_7");
  fSLP2->SetLineWidth(3); fSLP2->SetLineColor(GetColor(2)); data->Fit(fSLP2, "RSBQ0");
  fulllegend->AddEntry(fSLP2, "SLP2: Fit 0.45 - 0.55", "l");
  fullreslegend->AddEntry(fSLP2, TString::Format("%.2g = %.2g / %d", fSLP2->GetChisquare() / fSLP2->GetNDF(), fSLP2->GetChisquare(), fSLP2->GetNDF()).Data(), "l");
  printChi2(data, fSLP2);

  // Change fit range
  fitmin = 0.48, fitmax = 0.525;
  TF1* gex = new TF1("gex", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4]) + max(0., [5]*(x-[6]) * exp(-1. * (x-[6]) / [7]))", fitmin, fitmax);
  gex->FixParameter(0, g->GetParameter(0)); gex->SetParName(0, "gex_0");
  gex->FixParameter(1, g->GetParameter(1)); gex->SetParName(1, "gex_1");
  gex->FixParameter(2, g->GetParameter(2)); gex->SetParName(2, "gex_2");
  gex->FixParameter(3, g->GetParameter(3)); gex->SetParName(3, "gex_3");
  gex->FixParameter(4, g->GetParameter(4)); gex->SetParName(4, "gex_4");
  gex->SetParameter(5, 200e3 * 500); gex->SetParLimits(5, 100e5, 200e6);               gex->SetParName(5, "gex_5");
  gex->SetParameter(6, 0.5);         gex->SetParLimits(6, 0.45, 0.55);                 gex->SetParName(6, "gex_6");
  gex->SetParameter(7, 6e-3);        gex->SetParLimits(7, 5e-3, 0.5);                  gex->SetParName(7, "gex_7");
  gex->SetLineWidth(3); gex->SetLineColor(GetColor(3)); data->Fit(gex, "RSBQ0");
  fulllegend->AddEntry(gex, "g: Fixed pol1+G(x)", "l");
  fullreslegend->AddEntry(gex, TString::Format("%.2g = %.2g / %d", gex->GetChisquare() / gex->GetNDF(), gex->GetChisquare(), gex->GetNDF()).Data(), "l");

  double parmin, parmax;
  TF1* h1 = new TF1("h1", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4]) + max(0., [5]*(x-[6]) * exp(-1. * (x-[6]) / [7]))", fitmin, fitmax);
  h1->FixParameter(0, g->GetParameter(0)); h1->SetParName(0, "h1_0");
  h1->FixParameter(1, g->GetParameter(1)); h1->SetParName(1, "h1_1");
  g->GetParLimits(2, parmin, parmax);
  h1->SetParameter(2, g->GetParameter(2)); h1->SetParLimits(2, parmin, parmax);  h1->SetParName(2, "h1_2");
  g->GetParLimits(3, parmin, parmax);
  h1->SetParameter(3, g->GetParameter(3)); h1->SetParLimits(3, parmin, parmax);  h1->SetParName(3, "h1_3");
  g->GetParLimits(4, parmin, parmax);
  h1->SetParameter(4, g->GetParameter(4)); h1->SetParLimits(4, parmin, parmax);  h1->SetParName(4, "h1_4");
  h1->SetParameter(5, 200e3 * 500); h1->SetParLimits(5, 100e5, 200e6);               h1->SetParName(5, "h1_5");
  h1->SetParameter(6, 0.5);         h1->SetParLimits(6, 0.45, 0.55);                 h1->SetParName(6, "h1_6");
  h1->SetParameter(7, 6e-3);        h1->SetParLimits(7, 5e-3, 0.5);                  h1->SetParName(7, "h1_7");
  h1->SetLineWidth(3); h1->SetLineColor(GetColor(4)); data->Fit(h1, "RSBQ0");
  fulllegend->AddEntry(h1, "h1: Fixed pol1", "l");
  fullreslegend->AddEntry(h1, TString::Format("%.2g = %.2g / %d", h1->GetChisquare() / h1->GetNDF(), h1->GetChisquare(), h1->GetNDF()).Data(), "l");

  TF1* h2 = new TF1("h2", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4]) + max(0., [5]*(x-[6]) * exp(-1. * (x-[6]) / [7]))", fitmin, fitmax);
  h2->SetParameter(0, 22012.5); h2->SetParLimits(0, 0., 2. * 22012.5);  h2->SetParName(0, "h2_0");
  h2->SetParameter(1, 226650); h2->SetParLimits(1, 0., 2. * 226650);  h2->SetParName(1, "h2_1");
  g->GetParLimits(2, parmin, parmax);
  h2->SetParameter(2, g->GetParameter(2)); h2->SetParLimits(2, parmin, parmax);  h2->SetParName(2, "h2_2");
  g->GetParLimits(3, parmin, parmax);
  h2->SetParameter(3, g->GetParameter(3)); h2->SetParLimits(3, parmin, parmax);  h2->SetParName(3, "h2_3");
  g->GetParLimits(4, parmin, parmax);
  h2->SetParameter(4, g->GetParameter(4)); h2->SetParLimits(4, parmin, parmax);  h2->SetParName(4, "h2_4");
  h2->SetParameter(5, 200e3 * 500); h2->SetParLimits(5, 100e5, 200e6);               h2->SetParName(5, "h2_5");
  h2->SetParameter(6, 0.5);         h2->SetParLimits(6, 0.45, 0.55);                 h2->SetParName(6, "h2_6");
  h2->SetParameter(7, 6e-3);        h2->SetParLimits(7, 5e-3, 0.5);                  h2->SetParName(7, "h2_7");
  h2->SetLineWidth(3); h2->SetLineColor(GetColor(5)); data->Fit(h2, "RSBQ0");
  fulllegend->AddEntry(h2, "h2: Init with pol1+G", "l");
  fullreslegend->AddEntry(h2, TString::Format("%.2g = %.2g / %d", h2->GetChisquare() / h2->GetNDF(), h2->GetChisquare(), h2->GetNDF()).Data(), "l");

  TF1* h3 = new TF1("h3", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4]) + max(0., [5]*(x-[6]) * exp(-1. * (x-[6]) / [7]))", fitmin, fitmax);
  h3->SetParameter(0, 22012.5); h3->SetParLimits(0, 0., 10. * 22012.5);  h3->SetParName(0, "h3_0");
  h3->SetParameter(1, 226650); h3->SetParLimits(1, 0., 2. * 226650);  h3->SetParName(1, "h3_1");
  h3->FixParameter(2, g->GetParameter(2)); h3->SetParName(2, "h3_2");
  h3->FixParameter(3, g->GetParameter(3)); h3->SetParName(3, "h3_3");
  h3->FixParameter(4, g->GetParameter(4)); h3->SetParName(4, "h3_4");
  h3->SetParameter(5, 200e3 * 500); h3->SetParLimits(5, 100e5, 200e6);               h3->SetParName(5, "h3_5");
  h3->SetParameter(6, 0.5);         h3->SetParLimits(6, 0.45, 0.55);                 h3->SetParName(6, "h3_6");
  h3->SetParameter(7, 6e-3);        h3->SetParLimits(7, 5e-3, 0.5);                  h3->SetParName(7, "h3_7");
  h3->SetLineWidth(3); h3->SetLineColor(GetColor(6)); data->Fit(h3, "RSBQ0");
  fulllegend->AddEntry(h3, "h3: Fixed G", "l");
  fullreslegend->AddEntry(h3, TString::Format("%.2g = %.2g / %d", h3->GetChisquare() / h3->GetNDF(), h3->GetChisquare(), h3->GetNDF()).Data(), "l");

  // Residuals to compare
  TH1* residualSLP = makeResidual(data, fSLP);
  residualSLP->SetName("residualSLP");
  setStyle(residualSLP, 1);

  TH1* residualSLP2 = makeResidual(data, fSLP2);
  residualSLP2->SetName("residualSLP2");
  setStyle(residualSLP2, 2);
  printChi2(data, fSLP2);

  TH1* residualGex = makeResidual(data, gex);
  residualGex->SetName("residualGex");
  setStyle(residualGex, 3);

  TH1* residualH1 = makeResidual(data, h1);
  residualH1->SetName("residualH1");
  setStyle(residualH1, 4);

  TH1* residualH2 = makeResidual(data, h2);
  residualH2->SetName("residualH2");
  setStyle(residualH2, 5);

  TH1* residualH3 = makeResidual(data, h3);
  residualH3->SetName("residualH3");
  setStyle(residualH3, 6);

  // Plotting
  fitcanvas->cd(1);
  xMinFrame = data->GetXaxis()->GetXmin(); xMaxFrame = data->GetXaxis()->GetXmax();
  yMinFrame = 0.; yMaxFrame = 1.1 * getHistScale(data, false);
  TH1F* fullframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, "", "counts");
  fullframe->SetTitle((dataSet + ", " + ptText).c_str());

  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  fullframe->GetYaxis()->SetTitleOffset(1.0);
  fullframe->Draw();
  fulllegend->Draw("same");
  data->Draw("same");
  fSLP->Draw("same"); fSLP->SetRange(xMinFrame, xMaxFrame);
  fSLP2->Draw("same"); fSLP2->SetRange(xMinFrame, xMaxFrame);
  gex->Draw("same"); gex->SetRange(xMinFrame, xMaxFrame);
  h1->Draw("same");  h1->SetRange(xMinFrame, xMaxFrame);
  h2->Draw("same");  h2->SetRange(xMinFrame, xMaxFrame);
  h3->Draw("same");  h3->SetRange(xMinFrame, xMaxFrame);

  fitcanvas->cd(2);
  yMinFrame = -200e3, yMaxFrame = 250e3;
  TH1F* fullresframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "residual");
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.05);
  fullresframe->GetYaxis()->SetTitleOffset(1.0);
  fullresframe->Draw();
  fullreslegend->Draw("same");
  residualSLP->Draw("same");
  residualGex->Draw("same");
  residualH1->Draw("same");
  residualH2->Draw("same");
  residualH3->Draw("same");

  fitcanvas->cd(0);
  fitcanvas->SaveAs(fitcanvas->GetName());

  array<int, 2> myregion = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  array<int, 2> fitbins = getProjectionBins(data->GetXaxis(), fitmin, fitmax);
  double maxbin = data->GetNbinsX();

  // printParLimits(fSLP); printParLimits(gex); printParLimits(h1); printParLimits(h2); printParLimits(h3);
  // cout << "chi2:" << endl
  //      << "fSLP: " << mychi2(data, fSLP, fitbins[0], fitbins[1]) << ", full = " << mychi2(data, fSLP, myregion[0], myregion[1]) << endl
  //      << "gex: "  << mychi2(data, gex, fitbins[0], fitbins[1])  << ", full = " << mychi2(data, gex, myregion[0], myregion[1])  << endl
  //      << "h1: "   << mychi2(data, h1, fitbins[0], fitbins[1]) << ", full = " << mychi2(data, h1, myregion[0], myregion[1]) << endl
  //      << "h2: "   << mychi2(data, h2, fitbins[0], fitbins[1]) << ", full = " << mychi2(data, h2, myregion[0], myregion[1]) << endl
  //      << "h3: "   << mychi2(data, h3, fitbins[0], fitbins[1]) << ", full = " << mychi2(data, h3, myregion[0], myregion[1]) << endl;

  // cout << "nBins = " << fitbins[1] << " - " << fitbins[0] << " = " << fitbins[1] - fitbins[0] << endl;

  // cout << "fSLP: " << fSLP->GetNDF() << " =?= " << fitbins[1] - fitbins[0] << " - " << fSLP->GetNumberFreeParameters() << endl;
  // cout << "fSLP2: " << fSLP2->GetNDF() << " =?= " << fitbins[1] - fitbins[0] << " - " << fSLP2->GetNumberFreeParameters() << endl;
  // cout << "gex: " << gex->GetNDF() << " =?= " << fitbins[1] - fitbins[0] << " - " << gex->GetNumberFreeParameters() << endl;
  // cout << "h1: " << h1->GetNDF() << " =?= " << fitbins[1] - fitbins[0] << " - " << h1->GetNumberFreeParameters() << endl;
  // cout << "h2: " << h2->GetNDF() << " =?= " << fitbins[1] - fitbins[0] << " - " << h2->GetNumberFreeParameters() << endl;
  // cout << "h3: " << h3->GetNDF() << " =?= " << fitbins[1] - fitbins[0] << " - " << h3->GetNumberFreeParameters() << endl;

  //
  // Compare fit ranges: fSLP, fSLP2
  //
  string cName = (saveName + "_fitRange.pdf").c_str();
  TCanvas* canvas0 = new TCanvas(cName.c_str(), cName.c_str(), 1800, 2000);
  canvas0->Divide(1, 2, 0.01, 0.01);
  TLegend* legend0 = CreateLegend(0.175, 0.5, 0.55, 0.85, ""); // pol(1) + G(x) + x e^{-x/b}
  legend0->AddEntry(fSLP, "0.48 - 0.525", "l");
  legend0->AddEntry(fSLP2, "0.45 - 0.55", "l");
  TLegend* reslegend0 = CreateLegend(0.175, 0.5, 0.55, 0.85, "#chi^{2}/NDF: 0.48 - 0.525, 0.45 - 0.55");

  array<int, 2> bins = getProjectionBins(data->GetXaxis(), 0.48, 0.525);
  double chi2_narrow = mychi2(data, fSLP, 0.48, 0.525) / (1 + bins[1] - bins[0] - fSLP->GetNumberFreeParameters());
  bins = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  double chi2_wide   = mychi2(data, fSLP, 0.45, 0.55) / (1 + bins[1] - bins[0] - fSLP->GetNumberFreeParameters());
  reslegend0->AddEntry(fSLP, TString::Format("%.2g, %.2g", chi2_narrow, chi2_wide).Data(), "l");

  bins = getProjectionBins(data->GetXaxis(), 0.48, 0.525);
  chi2_narrow = mychi2(data, fSLP2, 0.48, 0.525) / (1 + bins[1] - bins[0] - fSLP2->GetNumberFreeParameters());
  bins = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  chi2_wide   = mychi2(data, fSLP2, 0.45, 0.55) / (1 + bins[1] - bins[0] - fSLP2->GetNumberFreeParameters());
  reslegend0->AddEntry(fSLP2, TString::Format("%.2g, %.2g", chi2_narrow, chi2_wide).Data(), "l");

  canvas0->cd(1);
  xMinFrame = data->GetXaxis()->GetXmin(); xMaxFrame = data->GetXaxis()->GetXmax();
  yMinFrame = 0.; yMaxFrame = 1.1 * getHistScale(data, false);
  TH1F* frame0 = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, "", "counts");
  frame0->SetTitle((dataSet + ", " + ptText).c_str());
  printChi2(data, fSLP2);

  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  frame0->GetYaxis()->SetTitleOffset(1.0);
  frame0->Draw();
  legend0->Draw("same");
  data->Draw("same");
  fSLP->Draw("same"); fSLP->SetRange(xMinFrame, xMaxFrame);
  fSLP2->Draw("same"); fSLP2->SetRange(xMinFrame, xMaxFrame);

  canvas0->cd(2);
  yMinFrame = -200e3, yMaxFrame = 250e3;
  TH1F* resframe0 = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "residual");
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.05);
  resframe0->GetYaxis()->SetTitleOffset(1.0);
  resframe0->Draw();
  reslegend0->Draw("same");
  residualSLP->Draw("same");
  residualSLP2->Draw("same");
  canvas0->SaveAs(canvas0->GetName());

  //
  // Compare pol1+G init/fix: fSLP, gex, h2
  //
  cName = (saveName + "_fitVSfixed.pdf").c_str();
  TCanvas* canvas1 = new TCanvas(cName.c_str(), cName.c_str(), 1800, 2000);
  canvas1->Divide(1, 2, 0.01, 0.01);
  TLegend* legend1 = CreateLegend(0.175, 0.5, 0.55, 0.85, "pol(1) + G(x) + x e^{-x/b}");
  legend1->AddEntry(fSLP2, "Fit in one go", "l");
  legend1->AddEntry(gex, "Fixed pol1+G(x)", "l");
  legend1->AddEntry(h2, "Init pol1+G(x)", "l");
  TLegend* reslegend1 = CreateLegend(0.175, 0.5, 0.55, 0.85, "#chi^{2}/NDF: 0.48 - 0.525, 0.45 - 0.55");

  bins = getProjectionBins(data->GetXaxis(), 0.48, 0.525);
  chi2_narrow = mychi2(data, fSLP2, 0.48, 0.525) / (1 + bins[1] - bins[0] - fSLP2->GetNumberFreeParameters());
  bins = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  chi2_wide   = mychi2(data, fSLP2, 0.45, 0.55) / (1 + bins[1] - bins[0] - fSLP2->GetNumberFreeParameters());
  reslegend1->AddEntry(fSLP2, TString::Format("%.2g, %.2g", chi2_narrow, chi2_wide).Data(), "l");

  bins = getProjectionBins(data->GetXaxis(), 0.48, 0.525);
  chi2_narrow = mychi2(data, gex, 0.48, 0.525) / (1 + bins[1] - bins[0] - gex->GetNumberFreeParameters());
  bins = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  chi2_wide   = mychi2(data, gex, 0.45, 0.55) / (1 + bins[1] - bins[0] - gex->GetNumberFreeParameters());
  reslegend1->AddEntry(gex, TString::Format("%.2g, %.2g", chi2_narrow, chi2_wide).Data(), "l");

  bins = getProjectionBins(data->GetXaxis(), 0.48, 0.525);
  chi2_narrow = mychi2(data, h2, 0.48, 0.525) / (1 + bins[1] - bins[0] - h2->GetNumberFreeParameters());
  bins = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  chi2_wide   = mychi2(data, h2, 0.45, 0.55) / (1 + bins[1] - bins[0] - h2->GetNumberFreeParameters());
  reslegend1->AddEntry(h2, TString::Format("%.2g, %.2g", chi2_narrow, chi2_wide).Data(), "l");

  canvas1->cd(1);
  xMinFrame = data->GetXaxis()->GetXmin(); xMaxFrame = data->GetXaxis()->GetXmax();
  yMinFrame = 0.; yMaxFrame = 1.1 * getHistScale(data, false);
  TH1F* frame1 = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, "", "counts");
  frame1->SetTitle((dataSet + ", " + ptText).c_str());

  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  frame1->GetYaxis()->SetTitleOffset(1.0);
  frame1->Draw();
  legend1->Draw("same");
  data->Draw("same");
  fSLP2->Draw("same"); fSLP2->SetRange(xMinFrame, xMaxFrame);
  gex->Draw("same"); gex->SetRange(xMinFrame, xMaxFrame);
  h2->Draw("same");  h2->SetRange(xMinFrame, xMaxFrame);

  canvas1->cd(2);
  yMinFrame = -200e3, yMaxFrame = 250e3;
  TH1F* resframe1 = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "residual");
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.05);
  resframe1->GetYaxis()->SetTitleOffset(1.0);
  resframe1->Draw();
  reslegend1->Draw("same");
  residualSLP2->Draw("same");
  residualGex->Draw("same");
  residualH2->Draw("same");
  canvas1->SaveAs(canvas1->GetName());

  //
  // Compare fit to partial fixed: fSLP, h1, h3
  //
  cName = (saveName + "_fitVSpartialfix.pdf").c_str();
  TCanvas* canvas2 = new TCanvas(cName.c_str(), cName.c_str(), 1800, 2000);
  canvas2->Divide(1, 2, 0.01, 0.01);
  TLegend* legend2 = CreateLegend(0.175, 0.5, 0.55, 0.85, "pol(1) + G(x) + x e^{-x/b}");
  legend2->AddEntry(fSLP2, "Fit in one go", "l");
  legend2->AddEntry(h1, "Fixed pol1", "l");
  legend2->AddEntry(h3, "Fixed G", "l");
  TLegend* reslegend2 = CreateLegend(0.175, 0.5, 0.55, 0.85, "#chi^{2}/NDF: 0.48 - 0.525, 0.45 - 0.55");

  bins = getProjectionBins(data->GetXaxis(), 0.48, 0.525);
  chi2_narrow = mychi2(data, fSLP2, 0.48, 0.525) / (1 + bins[1] - bins[0] - fSLP2->GetNumberFreeParameters());
  bins = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  chi2_wide   = mychi2(data, fSLP2, 0.45, 0.55) / (1 + bins[1] - bins[0] - fSLP2->GetNumberFreeParameters());
  reslegend2->AddEntry(fSLP2, TString::Format("%.2g, %.2g", chi2_narrow, chi2_wide).Data(), "l");

  bins = getProjectionBins(data->GetXaxis(), 0.48, 0.525);
  chi2_narrow = mychi2(data, h1, 0.48, 0.525) / (1 + bins[1] - bins[0] - h1->GetNumberFreeParameters());
  bins = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  chi2_wide   = mychi2(data, h1, 0.45, 0.55) / (1 + bins[1] - bins[0] - h1->GetNumberFreeParameters());
  reslegend2->AddEntry(h1, TString::Format("%.2g, %.2g", chi2_narrow, chi2_wide).Data(), "l");

  bins = getProjectionBins(data->GetXaxis(), 0.48, 0.525);
  chi2_narrow = mychi2(data, h3, 0.48, 0.525) / (1 + bins[1] - bins[0] - h3->GetNumberFreeParameters());
  bins = getProjectionBins(data->GetXaxis(), 0.45, 0.55);
  chi2_wide   = mychi2(data, h3, 0.45, 0.55) / (1 + bins[1] - bins[0] - h3->GetNumberFreeParameters());
  reslegend2->AddEntry(h3, TString::Format("%.2g, %.2g", chi2_narrow, chi2_wide).Data(), "l");
  canvas2->cd(1);
  xMinFrame = data->GetXaxis()->GetXmin(); xMaxFrame = data->GetXaxis()->GetXmax();
  yMinFrame = 0.; yMaxFrame = 1.1 * getHistScale(data, false);
  TH1F* frame2 = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, "", "counts");
  frame2->SetTitle((dataSet + ", " + ptText).c_str());

  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.1);
  frame2->GetYaxis()->SetTitleOffset(1.0);
  frame2->Draw();
  legend2->Draw("same");
  data->Draw("same");
  fSLP2->Draw("same"); fSLP2->SetRange(xMinFrame, xMaxFrame);
  h1->Draw("same");  h1->SetRange(xMinFrame, xMaxFrame);
  h3->Draw("same");  h3->SetRange(xMinFrame, xMaxFrame);

  canvas2->cd(2);
  yMinFrame = -200e3, yMaxFrame = 250e3;
  TH1F* resframe2 = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "residual");
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetTopMargin(0.05);
  resframe2->GetYaxis()->SetTitleOffset(1.0);
  resframe2->Draw();
  reslegend2->Draw("same");
  residualSLP2->Draw("same");
  residualH1->Draw("same");
  residualH3->Draw("same");
  canvas2->SaveAs(canvas2->GetName());
}

// Fit mass with pol1 + G(x) + xe-x
void plotPol1GausXex(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  double fitmin = 0.45, fitmax = 0.55;
  int pol = 1;
  // TF1* f = getFitShape(hadron, fitmin, fitmax, pol);
  TF1* f = new TF1("f", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4]) + max(0., [5]*(x-[6]) * exp(-1. * (x-[6]) / [7]))", fitmin, fitmax);
  // setFitParameters(f, hadron, ptmin, ptmax, pol);
  setFitParametersPol1GausXex(f, hadron, ptmin, ptmax);
  printParLimits(f);
  data->Fit(f, "RSBQ0");
  printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol1+G+tail";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";
  string ptText = TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", ptmin, ptmax).Data();
  TH1F* fitframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  fitframe->SetTitle((dataSet + ", " + ptText).c_str());

  fitframe->Draw();
  data->Draw("same");
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  fitcanvas->SaveAs(fitcanvas->GetName());

  // Plot parts
  TF1* bkg = new TF1("bkg", "[0] + [1]*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "[0] * TMath::Gaus(x, [1], [2])", fitmin, fitmax);
  TF1* tail = new TF1("tail", "max(0., [0] * (x - [1]) * exp(-1. * (x - [1]) / [2]))", fitmin, fitmax);
  vector<TF1*> functions = {bkg, sig, tail};
  setStyle(bkg, 2);
  setStyle(sig, 3);
  setStyle(tail, 4);

  bkg->SetParameter(0,  f->GetParameter(0));
  bkg->SetParameter(1,  f->GetParameter(1));
  sig->SetParameter(0,  f->GetParameter(2));
  sig->SetParameter(1,  f->GetParameter(3));
  sig->SetParameter(2,  f->GetParameter(4));
  tail->SetParameter(0, f->GetParameter(5));
  tail->SetParameter(1, f->GetParameter(6));
  tail->SetParameter(2, f->GetParameter(7));

  string partName = saveName + "_parts.pdf";
  TCanvas* partcanvas = new TCanvas(partName.c_str(), partName.c_str(), 1800, 900);
  partcanvas->cd();

  TH1F* partframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  partframe->SetTitle((dataSet + ", " + ptText).c_str());

  partframe->Draw();
  data->Draw("same");
  for (auto g : functions) {
    g->Draw("same");
  }
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  partcanvas->SaveAs(partcanvas->GetName());
}
// Fit mass with pol2 + G(x) + xe-x
void plotPol2GausXex(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  double fitmin = 0.45, fitmax = 0.55;
  int pol = 2;
  // TF1* f = getFitShape(hadron, fitmin, fitmax, pol);
  TF1* f = new TF1("f", "[0]+[1]*x+[2]*x*x + [3]*TMath::Gaus(x,[4],[5]) + max(0., [6]*(x-[7]) * exp(-1. * (x-[7]) / [8]))", fitmin, fitmax);
  // setFitParameters(f, hadron, ptmin, ptmax, pol);
  setFitParametersPol2GausXex(f, hadron, ptmin, ptmax);
  printParLimits(f);
  data->Fit(f, "RSBQ0");
  printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol2+G+tail";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";
  string ptText = TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", ptmin, ptmax).Data();
  TH1F* fitframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  fitframe->SetTitle((dataSet + ", " + ptText).c_str());

  fitframe->Draw();
  data->Draw("same");
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  fitcanvas->SaveAs(fitcanvas->GetName());

  // Plot parts
  TF1* bkg  = new TF1("bkg", "[0] + [1]*x + [2]*x*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "[0] * TMath::Gaus(x, [1], [2])", fitmin, fitmax);
  TF1* tail = new TF1("tail", "max(0., [0] * (x - [1]) * exp(-1. * (x - [1]) / [2]))", fitmin, fitmax);
  vector<TF1*> functions = {bkg, sig, tail};
  setStyle(bkg, 2);
  setStyle(sig, 3);
  setStyle(tail, 4);

  bkg->SetParameter(0,  f->GetParameter(0));
  bkg->SetParameter(1,  f->GetParameter(1));
  bkg->SetParameter(2,  f->GetParameter(2));
  sig->SetParameter(0,  f->GetParameter(3));
  sig->SetParameter(1,  f->GetParameter(4));
  sig->SetParameter(2,  f->GetParameter(5));
  tail->SetParameter(0, f->GetParameter(6));
  tail->SetParameter(1, f->GetParameter(7));
  tail->SetParameter(2, f->GetParameter(8));

  TLegend* legend = CreateLegend(0.5, 0.8, 0.7, 0.85, "");

  string partName = saveName + "_parts.pdf";
  TCanvas* partcanvas = new TCanvas(partName.c_str(), partName.c_str(), 1800, 900);
  partcanvas->cd();

  TH1F* partframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  partframe->SetTitle((dataSet + ", " + ptText).c_str());

  partframe->Draw();
  data->Draw("same");
  for (auto g : functions) {
    g->Draw("same");
  }
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  partcanvas->SaveAs(partcanvas->GetName());
}
// Fit mass with pol1 + G(x) + G(x)
void plotPol1GausGaus(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0] + [1]*x + [2] * TMath::Gaus(x, [3], [4]) + [5] * TMath::Gaus(x, [3], [6])", fitmin, fitmax);
  setFitParametersPol1GausGaus(f, hadron, ptmin, ptmax);
  printParLimits(f);
  data->Fit(f, "RSBQ0");
  printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol1+G+G";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";
  string ptText = TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", ptmin, ptmax).Data();
  TH1F* fitframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  fitframe->SetTitle((dataSet + ", " + ptText).c_str());

  fitframe->Draw();
  data->Draw("same");
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  fitcanvas->SaveAs(fitcanvas->GetName());

  // Plot parts
  TF1* bkg  = new TF1("bkg", "[0] + [1]*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "[0] * TMath::Gaus(x, [1], [2])", fitmin, fitmax);
  TF1* gaus = new TF1("gaus", "[0] * TMath::Gaus(x, [1], [2])", fitmin, fitmax);
  vector<TF1*> functions = {bkg, sig, gaus};
  setStyle(bkg, 2);
  setStyle(sig, 3);
  setStyle(gaus, 4);

  bkg->SetParameter(0,  f->GetParameter(0));
  bkg->SetParameter(1,  f->GetParameter(1));
  sig->SetParameter(0,  f->GetParameter(2));
  sig->SetParameter(1,  f->GetParameter(3));
  sig->SetParameter(2,  f->GetParameter(4));
  gaus->SetParameter(0, f->GetParameter(5));
  gaus->SetParameter(1, f->GetParameter(3));
  gaus->SetParameter(2, f->GetParameter(6));

  string partName = saveName + "_parts.pdf";
  TCanvas* partcanvas = new TCanvas(partName.c_str(), partName.c_str(), 1800, 900);
  partcanvas->cd();
  TH1F* partframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  partframe->SetTitle((dataSet + ", " + ptText).c_str());

  partframe->Draw();
  data->Draw("same");
  for (auto g : functions) {
    g->Draw("same");
  }
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  partcanvas->SaveAs(partcanvas->GetName());
}
void plotPol1GausExp(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  bool logplot = false;

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4]) * (x < ([3]+[4]*[5])) + [2]*TMath::Exp(-1.*(x - ([3]+[4]*[5]/2))/([4]/[5])) * (x > [3]+[4]*[5])", fitmin, fitmax);
  f->SetParameters(20e3, -1e3, 3e6, MassK0S, 1e-2, 2.); // a, b, A, mu, sigma, lambda
  // printParLimits(f);
  data->Fit(f, "RSBQ0");
  // printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol1+G+exp";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();
  if (logplot) fitcanvas->SetLogy(); // Easier to see exponentials

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = getHistLowerBound(data, false), yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame /= 2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";
  string ptText = TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", ptmin, ptmax).Data();
  TH1F* fitframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  fitframe->SetTitle((dataSet + ", " + ptText).c_str());

  fitframe->Draw();
  data->Draw("same");
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  fitcanvas->SaveAs(fitcanvas->GetName());

  // Plot parts
  double a = f->GetParameter(0), b = f->GetParameter(1), A = f->GetParameter(2), mu = f->GetParameter(3), sigma = f->GetParameter(4), lambda = f->GetParameter(5);
  TF1* bkg  = new TF1("bkg", "[0]+[1]*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "[0]*TMath::Gaus(x, [1], [2])", fitmin, mu + sigma * lambda);
  TF1* exp  = new TF1("exp", "[0]*TMath::Exp(-1.*(x-([1]+[2]*[3]/2))/([2]/[3]))", mu + sigma * lambda, fitmax);

  vector<TF1*> functions = {bkg, sig, exp};
  setStyle(bkg, 2);
  setStyle(sig, 3);
  setStyle(exp, 4);

  bkg->SetParameter(0, a);
  bkg->SetParameter(1, b);
  sig->SetParameter(0, A);
  sig->SetParameter(1, mu);
  sig->SetParameter(2, sigma);
  exp->SetParameter(0, A);
  exp->SetParameter(1, mu);
  exp->SetParameter(2, sigma);
  exp->SetParameter(3, lambda);
  // printParLimits(bkg); printParLimits(sig); printParLimits(exp);

  string partName = saveName + "_parts.pdf";
  TCanvas* partcanvas = new TCanvas(partName.c_str(), partName.c_str(), 1800, 900);
  partcanvas->cd();
  if (logplot) partcanvas->SetLogy();
  TH1F* partframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  partframe->SetTitle((dataSet + ", " + ptText).c_str());

  partframe->Draw();
  data->Draw("same");
  for (auto g : functions) {
    g->Draw("same"); //g->SetRange(xMinFrame, xMaxFrame);
  }
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  partcanvas->SaveAs(partcanvas->GetName());
}
void plotPol2GausExp(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  bool logplot = false;

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0]+[1]*x+[2]*x*x + [3]*TMath::Gaus(x,[4],[5]) * (x < ([4]+[5]*[6])) + [3]*TMath::Exp(-1.*(x - ([4]+[5]*[6]/2))/([5]/[6])) * (x > [4]+[5]*[6])", fitmin, fitmax);
  f->SetParameters(20e3, -1e3, -1e3, 3e6, MassK0S, 1e-2, 2.); // a, b, c, A, mu, sigma, lambda
  // printParLimits(f);
  data->Fit(f, "RSBQ0");
  // printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol1+G+exp";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();
  if (logplot) fitcanvas->SetLogy(); // Easier to see exponentials

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = getHistLowerBound(data, false), yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame /= 2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";
  string ptText = TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", ptmin, ptmax).Data();
  TH1F* fitframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  fitframe->SetTitle((dataSet + ", " + ptText).c_str());

  fitframe->Draw();
  data->Draw("same");
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  fitcanvas->SaveAs(fitcanvas->GetName());

  // Plot parts
  double a = f->GetParameter(0), b = f->GetParameter(1), c = f->GetParameter(2),
         A = f->GetParameter(3), mu = f->GetParameter(4), sigma = f->GetParameter(5), lambda = f->GetParameter(6);
  TF1* bkg  = new TF1("bkg", "[0]+[1]*x+[2]*x*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "[0]*TMath::Gaus(x, [1], [2])", fitmin, mu + sigma * lambda);
  TF1* exp  = new TF1("exp", "[0]*TMath::Exp(-1.*(x-([1]+[2]*[3]/2))/([2]/[3]))", mu + sigma * lambda, fitmax);

  vector<TF1*> functions = {bkg, sig, exp};
  setStyle(bkg, 2);
  setStyle(sig, 3);
  setStyle(exp, 4);

  bkg->SetParameter(0, a);
  bkg->SetParameter(1, b);
  bkg->SetParameter(2, c);
  sig->SetParameter(0, A);
  sig->SetParameter(1, mu);
  sig->SetParameter(2, sigma);
  exp->SetParameter(0, A);
  exp->SetParameter(1, mu);
  exp->SetParameter(2, sigma);
  exp->SetParameter(3, lambda);
  // printParLimits(bkg); printParLimits(sig); printParLimits(exp);

  string partName = saveName + "_parts.pdf";
  TCanvas* partcanvas = new TCanvas(partName.c_str(), partName.c_str(), 1800, 900);
  partcanvas->cd();
  if (logplot) partcanvas->SetLogy();
  TH1F* partframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  partframe->SetTitle((dataSet + ", " + ptText).c_str());

  partframe->Draw();
  data->Draw("same");
  for (auto g : functions) {
    g->Draw("same"); //g->SetRange(xMinFrame, xMaxFrame);
  }
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  partcanvas->SaveAs(partcanvas->GetName());
}
void plotBkgPartials(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  bool plotPol1 = false;
  bool plotPol2 = !plotPol1;

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  if (plotPol1)      saveName += "_fit=pol1+G_residuals";
  else if (plotPol2) saveName += "_fit=pol2+G_residuals";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 1800, 2000);
  canvas->Divide(1, 2, 0.01, 0.01);

  string ptText = TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", ptmin, ptmax).Data();
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistScale(data, false);

  if ("K0S" == hadron) {
    double signalWidth = 1e-2, peakVal = 3e6;
    double x0 = 0.45, x1 = 0.55, x2 = 0.57;
    double fitmin = 0.485, fitmax = 0.505;
    double y0 = data->GetBinContent(data->FindBin(x0 + 1e-3));
    double y1 = data->GetBinContent(data->FindBin(x1 + 1e-3));
    double y2 = data->GetBinContent(data->FindBin(x2 - 1e-3));

    TF1* f;
    if (plotPol1) {
      double b = (y1 - y0) / (x1 - x0);
      double a = y0 - b * x0;
      // TF1* f = new TF1("f", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4])", fitmin, fitmax);
      f = new TF1("f", "[0]+[1]*x + [2]*TMath::Gaus(x,[3],[4])", fitmin, fitmax);
      f->FixParameter(0, a);
      f->FixParameter(1, b);
      f->SetParameter(2, peakVal);     f->SetParLimits(2, 2e6, 6e6);
      f->SetParameter(3, MassK0S);     f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth); f->SetParLimits(4, 1e-3, 2e-2);
    }
    else if (plotPol2) {
      double c = ( (y2 - y0)/(x2 - x0) - (y1 - y0)/(x1 - x0) ) / (x2 - x1);
      double b = (y1 - y0) / (x1 - x0) - c * (x1 + x0);
      double a = y1 - b * x1 -c * x1 * x1;

      f = new TF1("f", "[0]+[1]*x+[2]*x*x + [3]*TMath::Gaus(x,[4],[5])", fitmin, fitmax);
      f->FixParameter(0, a);
      f->FixParameter(1, b);
      f->FixParameter(2, c);
      f->SetParameter(3, peakVal);     f->SetParLimits(3, 2e6, 6e6);
      f->SetParameter(4, MassK0S);     f->SetParLimits(4, 0.48, 0.51);
      f->SetParameter(5, signalWidth); f->SetParLimits(5, 1e-3, 2e-2);
    }

    printParLimits(f);
    data->Fit(f, "RSBQ0");
    printParLimits(f);
    setStyle(f, 1);

    canvas->cd(1);
    // gPad->SetLogy(); yMinFrame = 1e4; yMaxFrame *= 10.;
    TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "counts");
    frame->SetTitle((dataSet + ", " + ptText).c_str());
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.1);
    frame->GetYaxis()->SetTitleOffset(1.0);
    frame->Draw();
    data->Draw("same");
    f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);

    TF1* p;
    if (plotPol1) {
      p = new TF1("p", "pol1", f->GetXmin(), f->GetXmax());
      p->FixParameter(0, f->GetParameter(0));
      p->FixParameter(1, f->GetParameter(1));
      setStyle(p, 1, 7);
      p->Draw("same");
    } else if (plotPol2) {
      p = new TF1("p", "pol2", f->GetXmin(), f->GetXmax());
      p->FixParameter(0, f->GetParameter(0));
      p->FixParameter(1, f->GetParameter(1));
      p->FixParameter(2, f->GetParameter(2));
      setStyle(p, 1, 7);
      p->Draw("same");
    }

    canvas->cd(2);
    TH1* residual = makeResidual(data, f);
    residual->SetName("residual");
    setStyle(residual, 1);

    yMinFrame = 1.1 * getHistLowerBound(residual, false);
    yMaxFrame = 1.1 * getHistScale(residual, false);
    TH1F* resframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "residual");
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.05);
    resframe->GetYaxis()->SetTitleOffset(1.0);
    resframe->Draw();
    residual->Draw("same");

    // 0 at x = [1], max at x = [1]+[2], f([1]+[2]]) = [0]*[1]/e
    double gx0 = 0.505, gx1 = 0.55;
    if ("K0S" != hadron) {
      gx0 = 1.13; gx1 = 1.2;
    }
    TF1* g = new TF1("g", "max(0., [0]*(x-[1]) * exp(-1. * (x-[1]) / [2]))", gx0, gx1);
    double zeropoint = 0.505, maxpoint = 0.51, maxheight = 600e3;
    double amplitude = maxheight * exp(1) / (maxpoint - zeropoint);
    g->SetParameter(0, amplitude);  g->SetParLimits(0, 0., 10.*amplitude);
    g->SetParameter(1, zeropoint); g->SetParLimits(1, 0.5, 0.52);
    g->SetParameter(2, maxpoint - zeropoint); g->SetParLimits(2, 1e-3, 2e-2);
    printParLimits(g);
    residual->Fit(g, "RSBQ0");
    printParLimits(g);
    setStyle(g, 2);
    g->Draw("same"); g->SetRange(g->GetParameter(1), 0.55);

    TLegend* gLegend = CreateLegend(0.2, 0.4, 0.7, 0.85, "#chi^{2}/NDF");
    double gchi2 = g->GetChisquare();
    int gnfp     = g->GetNumberFreeParameters();
    int gnbins   = 1 + data->FindBin(gx1 - 1e-3) - data->FindBin(gx0 + 1e-3);
    gLegend->AddEntry(g, TString::Format("%.2g = %.2g / %d", gchi2 / (gnbins - gnfp), gchi2, gnbins - gnfp).Data());
    gLegend->Draw("same");
  } else { // (Anti)Lambda0
    xMinFrame = 1.08; xMaxFrame = 1.15;
    double x0 = 1.1, x1 = 1.13, x2 = 1.14;
    double fitmin = 1.1, fitmax = 1.13;

    double y0 = data->GetBinContent(data->FindBin(x0 + 1e-3));
    double y1 = data->GetBinContent(data->FindBin(x1 - 1e-3));
    double y2 = data->GetBinContent(data->FindBin(x2 - 1e-3));
    double c = ( (y2 - y0)/(x2 - x0) - (y1 - y0)/(x1 - x0) ) / (x2 - x1);
    double b = (y1 - y0) / (x1 - x0) - c * (x1 + x0);
    double a = y1 - b * x1 - c * x1 * x1;

    TF1* f = new TF1("f", "[0]+[1]*x+[2]*x*x + [3]*TMath::Gaus(x,[4],[5])", fitmin, fitmax);
    // f->SetParameter(0, a);     f->SetParLimits(0, 10.*a, -10.*a);
    // f->SetParameter(1, b);     f->SetParLimits(1, -10.*b, 10.*b);
    // f->SetParameter(2, c);     f->SetParLimits(2, -10.*c, 10.*c);
    f->FixParameter(0, a);
    f->FixParameter(1, b);
    f->FixParameter(2, c);
    // f->SetParameter(3, 15000e3);  f->SetParLimits(3, 900e3, 1500e3);
    // f->SetParameter(4, 1.116);   f->SetParLimits(4, 1.11, 1.12);
    // f->SetParameter(5, 5e-3);  f->SetParLimits(5, 1e-4, 1e-2);
    f->FixParameter(3, 10000e3);
    f->FixParameter(4, MassLambda0);
    f->FixParameter(5, 5e-4);

    printParLimits(f);
    data->Fit(f, "RSBQ0");
    printParLimits(f);
    setStyle(f, 1);

    canvas->cd(1);
    TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "counts");
    frame->SetTitle((dataSet + ", " + ptText).c_str());
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.05);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.1);
    frame->GetYaxis()->SetTitleOffset(1.0);
    frame->Draw();
    data->Draw("same");
    f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);

    canvas->cd(2);
    TH1* residual = makeResidual(data, f);
    residual->SetName("residual");
    setStyle(residual, 1);

    yMinFrame = 1.1 * getHistLowerBound(residual, false);
    yMaxFrame = 1.1 * getHistScale(residual, false);
    TH1F* resframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, "residual");
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.05);
    resframe->GetYaxis()->SetTitleOffset(1.0);
    resframe->Draw();
    residual->Draw("same");
  }

  canvas->SaveAs(canvas->GetName());
}

// -------------------------------------------------------------------------------------------------
//
// Interface
//
// -------------------------------------------------------------------------------------------------

void plotTrain(int train, string hadron, double v0min, double v0max, int setting)
{
  string inName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string dataSet = getDataSet(train);
  vector<string> inputStrings = {inName, dataSet, hadron};

  switch(setting) {
    case 0:
      plotPol1GausXex(inputStrings, v0min, v0max);
      break;
    case 1:
      plotPol2GausXex(inputStrings, v0min, v0max);
      break;
    case 2:
      plotPol1GausGaus(inputStrings, v0min, v0max);
      break;
    case 3:
      plotPol1GausExp(inputStrings, v0min, v0max);
      break;
    case 4:
      plotPol2GausExp(inputStrings, v0min, v0max);
      break;
  }
}

void plot252064(string hadron, double v0min, double v0max, int setting) { plotTrain(252064, hadron, v0min, v0max, setting); }
void plot282430(string hadron, double v0min, double v0max, int setting) { plotTrain(282430, hadron, v0min, v0max, setting); }
