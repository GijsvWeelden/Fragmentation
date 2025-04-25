
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

vector<double> polparams(vector<double> x, vector<double> y)
{
  if (x.size() != y.size()) {
    cout << "Error: polparams vectors must have the same size!" << endl;
    return {};
  }
  if (x.size() == 2) { // pol1
    double a = (y[1] - y[0]) / (x[1] - x[0]);
    double b = y[0] - a * x[0];
    return {b, a};
  }
  if (x.size() == 3) { // pol2
    double c = ( (y[2]-y[1])/(x[2]-x[1]) - (y[1]-y[0])/(x[1]-x[0]) ) / (x[2]-x[0]);
    double b = (y[1]-y[0])/(x[1]-x[0]) - c*(x[1]+x[0]);
    double a = y[0] - b*x[0] - c*x[0]*x[0];
    return {a, b, c};
  }
  return {};
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
// Functions to use with setup
//
// -------------------------------------------------------------------------------------------------

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
      f->SetParameter(6, gWidth);       f->SetParLimits(6, 1e-2, 0.1);
    } else if (ptmin > 4.9 && ptmax < 10.1) {
      pOffset = 20e3; pSlope = -1e3;
      peakVal = 2.80851e+06;
      gAmp = 1e6; gWidth = 2e-2;
      f->SetParameter(0, pOffset);      f->SetParLimits(0, 0., 15. * pOffset);
      f->SetParameter(1, pSlope);       f->SetParLimits(1, 1e3 * pSlope, 0.);
      f->SetParameter(2, peakVal);      f->SetParLimits(2, 0.8 * peakVal, 2 * peakVal);
      f->SetParameter(3, mass);         f->SetParLimits(3, 0.48, 0.51);
      f->SetParameter(4, signalWidth);  f->SetParLimits(4, 1e-3, 2e-2);
      f->SetParameter(5, gAmp);         f->SetParLimits(5, gAmp / 10., 2. * gAmp);
      f->SetParameter(6, gWidth);       f->SetParLimits(6, 1e-2, 0.1);
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
      f->SetParameter(6, gWidth);       f->SetParLimits(6, 1e-2, 0.1);
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
  double a = 20e3, b = -1e3, A = 3e6, mu = MassK0S, sigma = 1e-2, B = 3e5, rho = 2e-2;
  double x0 = 0.45, y0 = data->GetBinContent(data->FindBin(x0 + 1e-3));
  double x1 = 0.55, y1 = data->GetBinContent(data->FindBin(x1 - 1e-3));
  b = (y1 - y0) / (x1 - x0);
  a = y0 - b * x0;

  // f->SetParameters(a, b, A, mu, sigma, B, rho);
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
// Fit mass with pol2 + G(x) + G(x)
void plotPol2GausGaus(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0] + [1]*x + [2]*x*x +[3] * TMath::Gaus(x, [4], [5]) +[6] * TMath::Gaus(x, [4], [7])",
                   fitmin, fitmax);
  double a = 20e3, b = -1e3, c = -1e3, A = 3e6, mu = MassK0S, sigma = 1e-2, B = 3e5, rho = 2e-2;
  f->SetParameters(a, b, c, A, mu, sigma, B, rho);
  data->Fit(f, "RSBQ0");
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol2+G+G";
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
  TF1* gaus = new TF1("gaus", "[0] * TMath::Gaus(x, [1], [2])", fitmin, fitmax);
  vector<TF1*> functions = {bkg, sig, gaus};
  setStyle(bkg, 2);
  setStyle(sig, 3);
  setStyle(gaus, 4);

  bkg->SetParameter(0,  f->GetParameter(0));
  bkg->SetParameter(1,  f->GetParameter(1));
  bkg->SetParameter(2,  f->GetParameter(2));
  sig->SetParameter(0,  f->GetParameter(3));
  sig->SetParameter(1,  f->GetParameter(4));
  sig->SetParameter(2,  f->GetParameter(5));
  gaus->SetParameter(0, f->GetParameter(6));
  gaus->SetParameter(1, f->GetParameter(4));
  gaus->SetParameter(2, f->GetParameter(7));

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
// Fit mass with pol1 + [G(x) (x<a), e-x (x>a)]
void plotPol1GausExp(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  bool logplot = true;

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
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame = getHistLowerBound(data, false)/2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
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
// Fit mass with pol2 + [G(x) (x<a), e-x (x>a)]
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
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame = getHistLowerBound(data, false)/2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
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
// Fit mass with pol1 + G(x) + [G(x) (x<a), e-x (x>a)]
void plotPol1GausGausExp(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  bool logplot = true;

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0]+[1]*x +[2]*TMath::Gaus(x,[3],[4]) * (x < ([3]+[4]*[5])) +[2]*TMath::Exp(-1.*(x - ([3]+[4]*[5]/2))/([4]/[5])) * (x > [3]+[4]*[5]) +[6]*TMath::Gaus(x,[3],[7])",
                   fitmin, fitmax);

  double a = 20e3, b = -1e3, A = 3e6, mu = MassK0S, sigma = 1e-2, lambda = 2., B = 3e5, rho = 2e-2;
  f->SetParameters(a, b, A, mu, sigma, lambda, B, rho);
  // printParLimits(f);
  data->Fit(f, "RSBQ0");
  // printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol1+G+G+exp";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();
  if (logplot) fitcanvas->SetLogy(); // Easier to see exponentials

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame = getHistLowerBound(data, false)/2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
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
  a = f->GetParameter(0), b = f->GetParameter(1);
  A = f->GetParameter(2), mu = f->GetParameter(3), sigma = f->GetParameter(4), lambda = f->GetParameter(5);
  B = f->GetParameter(6), rho = f->GetParameter(7);
  TF1* bkg  = new TF1("bkg", "[0]+[1]*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "[0]*TMath::Gaus(x, [1], [2])", fitmin, mu + sigma * lambda);
  TF1* exp  = new TF1("exp", "[0]*TMath::Exp(-1.*(x-([1]+[2]*[3]/2))/([2]/[3]))", mu + sigma * lambda, fitmax);
  TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", fitmin, fitmax);

  vector<TF1*> functions = {bkg, sig, exp, gaus};
  setStyle(bkg, 2);
  setStyle(sig, 3);
  setStyle(exp, 3);
  setStyle(gaus, 4);

  bkg->SetParameter(0, a);
  bkg->SetParameter(1, b);
  sig->SetParameter(0, A);
  sig->SetParameter(1, mu);
  sig->SetParameter(2, sigma);
  exp->SetParameter(0, A);
  exp->SetParameter(1, mu);
  exp->SetParameter(2, sigma);
  exp->SetParameter(3, lambda);
  gaus->SetParameter(0, B);
  gaus->SetParameter(1, mu);
  // printParLimits(bkg); printParLimits(sig); printParLimits(exp); printParLimits(gaus);

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
// Fit mass with pol2 + G(x) + [G(x) (x<a), e-x (x>a)]
void plotPol2GausGausExp(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  bool logplot = true;

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0]+[1]*x+[2]*x*x +[3]*TMath::Gaus(x,[4],[5]) * (x < ([4]+[5]*[6])) +[3]*TMath::Exp(-1.*(x - ([4]+[5]*[6]/2))/([5]/[6])) * (x > [4]+[5]*[6]) +[7]*TMath::Gaus(x,[4],[8])",
                   fitmin, fitmax);

  double a = 20e3, b = -1e3, c = -1e3, A = 3e6, mu = MassK0S, sigma = 1e-2, lambda = 2., B = 3e5, rho = 2e-2;
  f->SetParameters(a, b, c, A, mu, sigma, lambda, B, rho);
  // printParLimits(f);
  data->Fit(f, "RSBQ0");
  // printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol2+G+G+exp";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();
  if (logplot) fitcanvas->SetLogy(); // Easier to see exponentials

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame = getHistLowerBound(data, false)/2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
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
  a = f->GetParameter(0), b = f->GetParameter(1), c = f->GetParameter(2);
  A = f->GetParameter(3), mu = f->GetParameter(4), sigma = f->GetParameter(5), lambda = f->GetParameter(6);
  B = f->GetParameter(7), rho = f->GetParameter(8);
  TF1* bkg  = new TF1("bkg", "[0]+[1]*x+[2]*x*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "[0]*TMath::Gaus(x, [1], [2])", fitmin, mu + sigma * lambda);
  TF1* exp  = new TF1("exp", "[0]*TMath::Exp(-1.*(x-([1]+[2]*[3]/2))/([2]/[3]))", mu + sigma * lambda, fitmax);
  TF1* gaus = new TF1("gaus", "[0]*TMath::Gaus(x, [1], [2])", fitmin, fitmax);

  vector<TF1*> functions = {bkg, sig, exp, gaus};
  setStyle(bkg, 2);
  setStyle(sig, 3);
  setStyle(exp, 3);
  setStyle(gaus, 4);

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
  gaus->SetParameter(0, B);
  gaus->SetParameter(1, mu);
  // printParLimits(bkg); printParLimits(sig); printParLimits(exp); printParLimits(gaus);

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
// Fit mass with pol1 + Breit-Wigner peak
void plotPol1BreitWigner(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  bool logplot = false;

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0]+[1]*x + breitwigner(x,[2],[3],[4])",
                   fitmin, fitmax);

  double a = 20e3, b = -1e3, A = 3e6, mu = MassK0S, Gamma = 1e-2;
  f->SetParameters(a, b, A, mu, Gamma);
  printParLimits(f);
  data->Fit(f, "RSBQ0");
  printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol1+BW";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();
  if (logplot) fitcanvas->SetLogy(); // Easier to see exponentials

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame = getHistLowerBound(data, false)/2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
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
  a = f->GetParameter(0), b = f->GetParameter(1);
  A = f->GetParameter(2), mu = f->GetParameter(3), Gamma = f->GetParameter(4);
  TF1* bkg  = new TF1("bkg", "[0]+[1]*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "breitwigner(x,[0],[1],[2])", fitmin, fitmax);

  vector<TF1*> functions = {bkg, sig};
  setStyle(bkg, 2);
  setStyle(sig, 3);

  bkg->SetParameter(0, a);
  bkg->SetParameter(1, b);
  sig->SetParameter(0, A);
  sig->SetParameter(1, mu);
  sig->SetParameter(2, Gamma);
  // printParLimits(bkg); printParLimits(sig); printParLimits(exp); printParLimits(gaus);

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
// Fit mass with pol2 + Breit-Wigner peak
void plotPol2BreitWigner(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  bool logplot = false;

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0]+[1]*x+[2]*x*x + breitwigner(x,[3],[4],[5])",
                   fitmin, fitmax);

  double a = 20e3, b = -1e3, c = -1e3, A = 3e6, mu = MassK0S, Gamma = 1e-2;
  f->SetParameters(a, b, c, A, mu, Gamma);
  printParLimits(f);
  data->Fit(f, "RSBQ0");
  printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol2+BW";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();
  if (logplot) fitcanvas->SetLogy(); // Easier to see exponentials

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame = getHistLowerBound(data, false)/2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
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
  a = f->GetParameter(0), b = f->GetParameter(1), c = f->GetParameter(2);
  A = f->GetParameter(3), mu = f->GetParameter(4), Gamma = f->GetParameter(5);
  TF1* bkg  = new TF1("bkg", "[0]+[1]*x+[2]*x*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "breitwigner(x,[0],[1],[2])", fitmin, fitmax);

  vector<TF1*> functions = {bkg, sig};
  setStyle(bkg, 2);
  setStyle(sig, 3);

  bkg->SetParameter(0, a);
  bkg->SetParameter(1, b);
  bkg->SetParameter(2, c);
  sig->SetParameter(0, A);
  sig->SetParameter(1, mu);
  sig->SetParameter(2, Gamma);
  // printParLimits(bkg); printParLimits(sig); printParLimits(exp); printParLimits(gaus);

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
// Fit mass with pol1 + Voigt peak
void plotPol1Voigt(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  bool logplot = false;

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0]+[1]*x + [2]*TMath::Voigt(x-[3],[4],[5])",
                   fitmin, fitmax);

  double a = 20e3, b = -1e3, A = 3e6, mu = MassK0S, sigma = 1e-2, gamma = 1e-2;
  f->SetParameters(a, b, A, mu, sigma, gamma);
  printParLimits(f);
  data->Fit(f, "RSBQ0");
  printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol1+V";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();
  if (logplot) fitcanvas->SetLogy(); // Easier to see exponentials

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame = getHistLowerBound(data, false)/2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
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
  a = f->GetParameter(0), b = f->GetParameter(1);
  A = f->GetParameter(2), mu = f->GetParameter(3), sigma = f->GetParameter(4), gamma = f->GetParameter(5);
  TF1* bkg  = new TF1("bkg", "[0]+[1]*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "[0]*TMath::Voigt(x-[1],[2],[3])", fitmin, fitmax);

  vector<TF1*> functions = {bkg, sig};
  setStyle(bkg, 2);
  setStyle(sig, 3);

  bkg->SetParameter(0, a);
  bkg->SetParameter(1, b);
  sig->SetParameter(0, A);
  sig->SetParameter(1, mu);
  sig->SetParameter(2, sigma);
  sig->SetParameter(3, gamma);
  // printParLimits(bkg); printParLimits(sig); printParLimits(exp); printParLimits(gaus);

  string partName = saveName + "_parts.pdf";
  TCanvas* partcanvas = new TCanvas(partName.c_str(), partName.c_str(), 1800, 900);
  partcanvas->cd();
  if (logplot) partcanvas->SetLogy();
  TH1F* partframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  partframe->SetTitle((dataSet + ", " + ptText).c_str());

  partframe->Draw();
  data->Draw("same");
  for (auto g : functions) {
    printParLimits(g);
    g->Draw("same"); //g->SetRange(xMinFrame, xMaxFrame);
  }
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  partcanvas->SaveAs(partcanvas->GetName());
}
// Fit mass with pol1 + Voigt peak
void plotPol2Voigt(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  bool logplot = false;

  TH1D* data = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  setStyle(data, 0);

  // Enforce same mean for both Gaussians
  double fitmin = 0.45, fitmax = 0.55;
  TF1* f = new TF1("f", "[0]+[1]*x+[2]*x*x + [3]*TMath::Voigt(x-[4],[5],[6])",
                   fitmin, fitmax);

  double a = 20e3, b = 1e3, c = -1e3, A = 3e6, mu = MassK0S, sigma = 1e-2, gamma = 1e-2;

  double x0 = 0.4, y0 = data->GetBinContent(data->FindBin(x0));
  double x1 = 0.45, y1 = data->GetBinContent(data->FindBin(x1));
  double x2 = 0.55,  y2 = data->GetBinContent(data->FindBin(x2));
  c = ( (y2-y1)/(x2-x1) - (y1-y0)/(x1-x0) ) / (x2-x0);
  b = (y1-y0)/(x1-x0) - c*(x1+x0);
  a = y0 - b*x0 - c*x0*x0;

  f->SetParameters(a, b, c, A, mu, sigma, gamma);
  printParLimits(f);
  data->Fit(f, "RSBQ0");
  printParLimits(f);
  setStyle(f, 1);

  string saveName = hadron;
  saveName += "_";
  saveName += data->GetName(); // hist name contains pt range
  saveName += "_fit=pol2+V";
  string fitName = saveName + ".pdf";
  TCanvas* fitcanvas = new TCanvas(fitName.c_str(), fitName.c_str(), 1800, 900);
  fitcanvas->cd();
  if (logplot) fitcanvas->SetLogy(); // Easier to see exponentials

  double xMinFrame = data->GetXaxis()->GetXmin(), xMaxFrame = data->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistUpperBound(data, false);
  if (logplot) yMinFrame = getHistLowerBound(data, false)/2., yMaxFrame = pow(10., ceil(log10(yMaxFrame))); // round up yMaxFrame to next power of 10
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
  a = f->GetParameter(0), b = f->GetParameter(1), c = f->GetParameter(2);
  A = f->GetParameter(3), mu = f->GetParameter(4), sigma = f->GetParameter(5), gamma = f->GetParameter(6);
  TF1* bkg  = new TF1("bkg", "[0]+[1]*x+[2]*x*x", fitmin, fitmax);
  TF1* sig  = new TF1("sig", "[0]*TMath::Voigt(x-[1],[2],[3])", fitmin, fitmax);

  vector<TF1*> functions = {bkg, sig};
  setStyle(bkg, 2);
  setStyle(sig, 3);

  bkg->SetParameter(0, a);
  bkg->SetParameter(1, b);
  bkg->SetParameter(2, c);
  sig->SetParameter(0, A);
  sig->SetParameter(1, mu);
  sig->SetParameter(2, sigma);
  sig->SetParameter(3, gamma);
  // printParLimits(bkg); printParLimits(sig); printParLimits(exp); printParLimits(gaus);

  string partName = saveName + "_parts.pdf";
  TCanvas* partcanvas = new TCanvas(partName.c_str(), partName.c_str(), 1800, 900);
  partcanvas->cd();
  if (logplot) partcanvas->SetLogy();
  TH1F* partframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  partframe->SetTitle((dataSet + ", " + ptText).c_str());

  partframe->Draw();
  data->Draw("same");
  for (auto g : functions) {
    printParLimits(g);
    g->Draw("same"); //g->SetRange(xMinFrame, xMaxFrame);
  }
  f->Draw("same"); f->SetRange(xMinFrame, xMaxFrame);
  partcanvas->SaveAs(partcanvas->GetName());
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
      plotPol2GausGaus(inputStrings, v0min, v0max);
      break;
    case 4:
      plotPol1GausExp(inputStrings, v0min, v0max);
      break;
    case 5:
      plotPol2GausExp(inputStrings, v0min, v0max);
      break;
    case 6:
      plotPol1GausGausExp(inputStrings, v0min, v0max);
      break;
    case 7:
      plotPol2GausGausExp(inputStrings, v0min, v0max);
      break;
    case 8:
      plotPol1BreitWigner(inputStrings, v0min, v0max);
      break;
    case 9:
      plotPol2BreitWigner(inputStrings, v0min, v0max);
      break;
    case 10:
      plotPol1Voigt(inputStrings, v0min, v0max);
      break;
    case 11:
      plotPol2Voigt(inputStrings, v0min, v0max);
      break;
    default:
      cout << "Invalid setting" << endl;
      return;
  }
}

void plot252064(string hadron, double v0min, double v0max, int setting) { plotTrain(252064, hadron, v0min, v0max, setting); }
void plot282430(string hadron, double v0min, double v0max, int setting) { plotTrain(282430, hadron, v0min, v0max, setting); }

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// EVERYTHING ABOVE THIS IS DEPRECATED
// USE ONLY FOR ADDING NEW FITTING FUNCTIONS TO THE STRUCTS BELOW
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

#ifndef __PLOTMASSFITS_H__
#define __PLOTMASSFITS_H__

// Return hist scale in a range of the histogram
double getHistScaleInRange(TH1* h, int minBin, int maxBin, bool doError, bool doMin) {
  double scale = doMin ? 1e12 : -1e12;
  for (int i = minBin; i <= maxBin; i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    double s = bc;
    if (doError)
      doMin ? s -= be : s += be;

    if (doMin)
      scale = min(scale, s);
    else
      scale = max(scale, s);
  }
  return scale;
}

int getExtremeBinInRange(TH1* h, int minBin, int maxBin, bool doError, bool doMin) {
  int BIN = -1;
  double scale = doMin ? 1e12 : -1e12;
  for (int i = minBin; i <= maxBin; i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    double s = bc;
    if (doError)
      doMin ? s -= be : s += be;

    if (doMin && s < scale) {
      scale = s;
      BIN = i;
    }
    else if (s > scale) {
      scale = s;
      BIN = i;
    }
  }
  return BIN;
}

double chisqInRange(TH1* hist, TF1* fit, int minBin, int maxBin) {
  double chisq = 0;
  for (int i = minBin; i <= maxBin; i++) {
    double bc = hist->GetBinContent(i);
    double be = hist->GetBinError(i);
    double fx = fit->Eval(hist->GetBinCenter(i));
    chisq += (fx - bc) * (fx - bc) / (be * be);
  }
  return chisq;
}

double chisqInRange(TH1* hist, TF1* fit, double xmin, double xmax) {
  array<int, 2> bins = getProjectionBins(hist->GetXaxis(), xmin, xmax);
  return chisqInRange(hist, fit, bins[0], bins[1]);
}

double ndfInRange(TF1* fit, int minBin, int maxBin) {
  return (1 + maxBin - minBin - fit->GetNpar());
}

double ndfInRange(TH1* data, TF1* fit, double xmin, double xmax) {
  array<int, 2> bins = getProjectionBins(data->GetXaxis(), xmin, xmax);
  return ndfInRange(fit, bins[0], bins[1]);
}

double roundToNextPowerOfTen(double x) {
  return std::pow(10., std::ceil(std::log10(x)));
}
double roundToPrevPowerOfTen(double x) {
  return std::pow(10., std::floor(std::log10(x)));
}

// Prints the parameter limits of a function
void printParLimits(TF1* f) {
  for (int i = 0; i < f->GetNpar(); i++) {
    double min, max;
    f->GetParLimits(i, min, max);
    cout << f->GetName() << " (" << i << " = " << f->GetParName(i) << ") " << f->GetParameter(i) << " (" << min << ", " << max << ")" << endl;
  }
  cout << endl;
}

// -------------------------------------------------------------------------------------------------
//
// Plotting struct
//
// -------------------------------------------------------------------------------------------------

struct InputSettings{
  private:
  public:
    int train = 0;
    string hadron = "";
    string fitName = "";
    string histName = "";
    string inputFileName = "";
    string outputFileName = "";
    double ptmin = -1e3, ptmax = -1e3, lowpt = -1e3, highpt = -1e3;
    double fitmin = -1e3, fitmax = -1e3; // Fit range
    bool logplot = false;
    bool normaliseData = false;
    vector<vector<double>> ptBinEdges = {};

    double massWindowMin = -1., massWindowMax = -1.;
    double polInitx0 = -1., polInitx1 = -1., polInitx2 = -1.;
    double sigalRegionMin = -1., sigalRegionMax = -1.;

    double getMass();
    string getFitExpression();
    string getSaveNameFromPt(string prefix, string suffix);
    int inputIssue(string home, string obj);
    void setFitX(double a, double b);
    string setInputFileNameFromTrain();
    string setInputFileNameFromFit();
    void setMassWindow(double a, double b);
    void setMassWindowDiff(double a, double b);
    void setPt(double a, double b);
    void setPtBinEdgesSorted(vector<vector<double>> x);
    void setPolInitX(double x0, double x1, double x2);
    void setSignalRegion(double a, double b);
    template <typename T> int writeOutputToFile(T* obj);
};

double InputSettings::getMass() {
  if (this->hadron == "K0S")
    return MassK0S;
  if (this->hadron == "Lambda" || this->hadron == "AntiLambda")
    return MassLambda;

  return -1.;
}

string InputSettings::getFitExpression() {
  string s = "";
  if (this->fitName == "pol2GausExp") {
    s = "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - ([1]+[2]/[3])/2.)/([2]/[3])) * (x > [1]+[2]*[3]) + [4]+[5]*x+[6]*x*x";
  } else if (this->fitName == "pol1BreitWigner") {
    s = "breitwigner(x,[0],[1],[2]) + [3]+[4]*x";
  } else
    cout << "InputSettings::getFitExpression() Error: requested unknown function" << endl;

  return s;
}

string InputSettings::getSaveNameFromPt(string prefix, string suffix = "") {
  string s = TString::Format("%s_pt%.1f-%.1f%s", prefix.c_str(), this->lowpt, this->highpt, suffix.c_str()).Data();
  return s;
}

int InputSettings::inputIssue(string home, string obj) {
  if (obj == "hadron") {
    if (this->hadron == "K0S" || this->hadron == "Lambda" || this->hadron == "AntiLambda") return 0;
  } else {
    cout << "InputSettings::inputIssue() Error: requested unknown object " << obj << endl;
    return 1;
  }

  cout << "InputSettings::" << home << "() Error: " << obj << " not initialised" << endl;
  return 2;
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(this->train) + "/AnalysisResults.root";
  this->inputFileName = s;
  return s;
}

string InputSettings::setInputFileNameFromFit() {
  string s = "./" + this->hadron + "_" + this->fitName + ".root";
  this->inputFileName = s;
  return s;
}

void InputSettings::setFitX(double a, double b) {
  if (a > b) {
    cout << "InputSettings::setFitX() Error: fitmin > fitmax" << endl;
    return;
  }
  this->fitmin = a;
  this->fitmax = b;
}

void InputSettings::setMassWindow(double a, double b) {
  if (a > b) {
    cout << "InputSettings::setMassWindow() Error: massWindowMin > massWindowMax" << endl;
    return;
  }
  this->massWindowMin = a;
  this->massWindowMax = b;
}

void InputSettings::setMassWindowDiff(double a, double b) {
  if (a > b) {
    cout << "InputSettings::setMassWindowDiff() Error: massWindowMin > massWindowMax" << endl;
    return;
  }
  double mass = this->getMass();
  this->setMassWindow(mass - a, mass + b);
}

void InputSettings::setPolInitX(double x0, double x1 = 0, double x2 = 0) {
  if (x0 > x1 || x0 > x2 || x1 > x2) {
    cout << "InputSettings::setPolInitX() Error: should have polInitx0 < polInitx1 < polInitx2" << endl;
    return;
  }
  this->polInitx0 = x0;
  this->polInitx1 = x1;
  this->polInitx2 = x2;
}

void InputSettings::setPt(double a, double b) {
  if (a > b) {
    cout << "InputSettings::setPt() Error: ptmin > ptmax" << endl;
    return;
  }
  this->ptmin = a;
  this->ptmax = b;
  this->lowpt = a;
  this->highpt = b;
}

void InputSettings::setPtBinEdgesSorted(vector<vector<double>> x) {
  sort(x.begin(), x.end());
  this->ptBinEdges = x;
}

void InputSettings::setSignalRegion(double a, double b) {
  if (a > b) {
    cout << "InputSettings::setSignalRegion() Error: sigalRegionMin > sigalRegionMax" << endl;
    return;
  }
  this->sigalRegionMin = a;
  this->sigalRegionMax = b;
}

template <typename T>
int InputSettings::writeOutputToFile(T* obj) {
  if (!obj)
    return 1;

  TFile* file = TFile::Open(this->outputFileName.c_str(), "UPDATE");
  obj->Write(obj->GetName(), TObject::kOverwrite);
  file->Close();
  return 0;
}

struct MassFitter {
  private:
  public:
    InputSettings inputs;
    TH1* data = nullptr;
    TF1* fit = nullptr;
    vector<TF1*> fitParts = {};
    TH1* fitResults = nullptr;
    TH1* fitParams = nullptr;
    TH1* pull = nullptr;
    TH1* residual = nullptr;

    double signal = -1.;
    double background = -1.;
    double signalPlusBackground = -1.;

    MassFitter() { InputSettings(); }
    MassFitter(InputSettings& x) { this->inputs = x; }

    void calcSigBkg();
    TF1* loadFitFunction();
    TH1* loadMassHist();
    void setFitInitialValues();
    vector<TF1*> loadFitParts();
    TH1* loadFitParams();
    TH1* loadFitResults();
    TH1* loadResidualHist();
    TH1* loadPullHist();
    string setHistNameFromTrain();
    void setSignalRegionFromSigma(double nSigma);
    void writeOutputsToFile();
};

void MassFitter::calcSigBkg() {
  if (!this->data || !this->fit || this->fitParts.size() == 0)
    return;

  int iBkg = 0;
  if (this->inputs.fitName == "pol2GausExp") {
    iBkg = 2;
  } else if (this->inputs.fitName == "pol1BreitWigner") {
    iBkg = 1;
  } else {
    cout << "MassFitter::calcSigBkg() Error: do not know how to set parameters for this function" << endl;
    return;
  }

  if (this->inputs.sigalRegionMin < 0 || this->inputs.sigalRegionMax < 0) {
    cout << "MassFitter::calcSigBkg() Automatically setting signal region to mu +/- 3 sigma" << endl;
    this->setSignalRegionFromSigma(3.);
  }

  array<int, 2> sigBins = getProjectionBins(this->data->GetXaxis(), this->inputs.sigalRegionMin, this->inputs.sigalRegionMax);
  double xmin = this->data->GetXaxis()->GetBinLowEdge(sigBins[0]);
  double xmax = this->data->GetXaxis()->GetBinUpEdge(sigBins[1]);

  TF1* bkgFit = this->fitParts[iBkg];
  bkgFit->SetRange(xmin, xmax);
  this->background = bkgFit->Integral(xmin, xmax);
  this->signalPlusBackground = this->data->Integral(sigBins[0], sigBins[1], "width");
  this->signal = this->signalPlusBackground - this->background;
}

TF1* MassFitter::loadFitFunction() {
  string saveName = this->inputs.getSaveNameFromPt("fit");
  string expression = this->inputs.getFitExpression();
  TF1* f = new TF1(saveName.c_str(), expression.c_str(), this->inputs.fitmin, this->inputs.fitmax);
  this->fit = (TF1*)f->Clone();
  return f;
}

TH1* MassFitter::loadMassHist() {
  TFile* f = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!f) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return nullptr;
  }

  if (this->inputs.histName == "") {
    this->setHistNameFromTrain();
    cout << "MassFitter::loadMassHist(): histName not set, setting to: " << this->inputs.histName << endl;
  }
  THnSparse* hist = (THnSparse*)f->Get(this->inputs.histName.c_str());
  if (!hist) {
    cout << "Could not find histogram " << this->inputs.histName << " in file " << this->inputs.inputFileName << endl;
    return nullptr;
  }

  int projectionAxis = 1;
  if (this->inputs.hadron == "Lambda")
    projectionAxis = 2;
  if (this->inputs.hadron == "AntiLambda")
    projectionAxis = 3;

  array<int, 2> bins = getProjectionBins(hist->GetAxis(0), this->inputs.ptmin, this->inputs.ptmax);
  int minBin = bins[0], maxBin = bins[1];
  this->inputs.lowpt = hist->GetAxis(0)->GetBinLowEdge(minBin);
  this->inputs.highpt = hist->GetAxis(0)->GetBinUpEdge(maxBin);

  hist->GetAxis(0)->SetRange(minBin, maxBin);
  TH1* h = hist->Projection(projectionAxis);

  if (this->inputs.normaliseData) {
    bins = getProjectionBins(h->GetXaxis(), this->inputs.massWindowMin, this->inputs.massWindowMax);
    minBin = bins[0], maxBin = bins[1];
    h->Scale(1./getHistScaleInRange(h, minBin, maxBin, false, false));
  }

  string hName = this->inputs.getSaveNameFromPt("data");
  this->data = (TH1*)h->Clone(hName.c_str());
  return h;
}

void MassFitter::setFitInitialValues() {
  array<int, 2> massBins = getProjectionBins(this->data->GetXaxis(), this->inputs.massWindowMin, this->inputs.massWindowMax);
  int minBin = massBins[0], maxBin = massBins[1];
  int extremeBin = getExtremeBinInRange(this->data, minBin, maxBin, false, false);

  int b0 = this->data->GetXaxis()->FindBin(this->inputs.polInitx0);
  int b1 = this->data->GetXaxis()->FindBin(this->inputs.polInitx1);
  int b2 = this->data->GetXaxis()->FindBin(this->inputs.polInitx2);
  double x0 = this->data->GetBinLowEdge(b0);
  double x1 = this->data->GetBinLowEdge(b1);
  double x2 = this->data->GetBinLowEdge(b2);
  double y0 = this->data->GetBinContent(b0);
  double y1 = this->data->GetBinContent(b1);
  double y2 = this->data->GetBinContent(b2);

  if (this->inputs.fitName == "pol2GausExp") {
    double A, mu, sigma, lambda, a, b, c;
    A = this->data->GetBinContent(extremeBin);
    mu = this->data->GetXaxis()->GetBinCenter(extremeBin);

    c = ( (y2-y1)/(x2-x1) - (y1-y0)/(x1-x0) ) / (x2-x0);
    b = (y1-y0)/(x1-x0) - c*(x1+x0);
    a = y0 - b*x0 - c*x0*x0;

    if (this->inputs.hadron == "K0S") {
      sigma = 1e-2;
      lambda = 2.;
    } else {
      sigma = 0.001;
      lambda = 2.;
    }
    this->fit->SetParameters(A, mu, sigma, lambda, a, b, c);
    this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
  } else if (this->inputs.fitName == "pol1BreitWigner") {
    double A, mu, Gamma, a, b;
    A = this->data->GetBinContent(extremeBin);
    mu = this->data->GetXaxis()->GetBinCenter(extremeBin);

    b = (y1-y0)/(x1-x0);
    a = y0 - b*x0;

    if (this->inputs.hadron == "K0S") {
      Gamma = 1e-2;
    } else {
      Gamma = 0.001;
    }

    this->fit->SetParameters(A, mu, Gamma, a, b);
    this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
    this->fit->SetParLimits(3, 0, 1);
    this->fit->SetParLimits(4, -1, 1);
  } else {
    cout << "InputSettings::setFitInitialValues() Error: do not have initial values for this function" << endl;
  }
}

vector<TF1*> MassFitter::loadFitParts() {
  if (!this->fit)
    return {};

  string fName = this->fit->GetName();
  vector<TF1*> parts;
  if (this->inputs.fitName == "pol2GausExp") {
    double A = this->fit->GetParameter(0);
    double mu = this->fit->GetParameter(1);
    double sigma = this->fit->GetParameter(2);
    double lambda = this->fit->GetParameter(3);
    double a = this->fit->GetParameter(4);
    double b = this->fit->GetParameter(5);
    double c = this->fit->GetParameter(6);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3]))", this->inputs.fitmin, mu+sigma*lambda);
    f1->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Exp(-(x - ([1]+[2]/[3])/2.)/([2]/[3])) * (x > [1]+[2]*[3])", mu+sigma*lambda, this->inputs.fitmax);
    f2->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs.fitmin, this->inputs.fitmax);
    f3->SetParameters(a, b, c);
    parts.push_back(f3);
  } else if (this->inputs.fitName == "pol1BreitWigner") {
    double A = this->fit->GetParameter(0);
    double mu = this->fit->GetParameter(1);
    double Gamma = this->fit->GetParameter(2);
    double a = this->fit->GetParameter(3);
    double b = this->fit->GetParameter(4);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "breitwigner(x,[0],[1],[2])", this->inputs.fitmin, this->inputs.fitmax);
    f1->SetParameters(A, mu, Gamma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs.fitmin, this->inputs.fitmax);
    f2->SetParameters(a, b);
    parts.push_back(f2);
  } else {
    cout << "InputSettings::setFitParts() Error: do not know the parts of this function" << endl;
  }
  this->fitParts = parts;
  return parts;
}

TH1* MassFitter::loadFitParams() {
  if (!this->fit)
    return nullptr;

  string saveName = this->inputs.getSaveNameFromPt("fitParams");

  int nPars = this->fit->GetNpar();
  TH1* h = new TH1D(saveName.c_str(), "fitParams", nPars, -0.5, 1.*nPars - 0.5);
  h->SetTitle("Fit parameters");
  h->GetXaxis()->SetTitle("");
  h->GetYaxis()->SetTitle("Value");
  for (int i = 0; i < nPars; i++) {
    h->GetXaxis()->SetBinLabel(i + 1, this->fit->GetParName(i));
    h->SetBinContent(i + 1, this->fit->GetParameter(i));
  }
  this->fitParams = (TH1*)h->Clone();
  return h;
}

TH1* MassFitter::loadFitResults() {
  if (!this->fit)
    return nullptr;

  string saveName = this->inputs.getSaveNameFromPt("fitResults");

  int nBins = 12;
  TH1* h = new TH1D(saveName.c_str(), "fitResults", nBins, -0.5, 1.*nBins - 0.5);
  h->SetTitle("Fit results");
  h->GetXaxis()->SetTitle("");
  h->GetXaxis()->SetBinLabel(1, "#chi^{2} (fit)");
  h->GetXaxis()->SetBinLabel(2, "NDF (fit)");
  h->GetXaxis()->SetBinLabel(3, "#chi^{2}/NDF (fit)");
  h->GetXaxis()->SetBinLabel(4, "#chi^{2} (all)");
  h->GetXaxis()->SetBinLabel(5, "NDF (all)");
  h->GetXaxis()->SetBinLabel(6, "#chi^{2}/NDF (all)");
  h->GetXaxis()->SetBinLabel(7, "Sig+Bkg");
  h->GetXaxis()->SetBinLabel(8, "Sig");
  h->GetXaxis()->SetBinLabel(9, "Bkg");
  h->GetXaxis()->SetBinLabel(10, "Sig/(Sig+Bkg)");
  h->GetXaxis()->SetBinLabel(11, "Sig/Bkg");
  h->GetXaxis()->SetBinLabel(12, "Sig/sqrt(Sig+Bkg)");

  double chisqFit = chisqInRange(this->data, this->fit, this->inputs.fitmin, this->inputs.fitmax);
  double ndfFit = ndfInRange(this->data, this->fit, this->inputs.fitmin, this->inputs.fitmax);
  double chisqAll = chisqInRange(this->data, this->fit, this->data->GetXaxis()->GetXmin(), this->data->GetXaxis()->GetXmax());
  double ndfAll = ndfInRange(this->data, this->fit, this->data->GetXaxis()->GetXmin(), this->data->GetXaxis()->GetXmax());
  this->calcSigBkg();

  h->SetBinContent(1, chisqFit);
  h->SetBinContent(2, ndfFit);
  h->SetBinContent(3, chisqFit / ndfFit);
  h->SetBinContent(4, chisqAll);
  h->SetBinContent(5, ndfAll);
  h->SetBinContent(6, chisqAll / ndfAll);
  h->SetBinContent(7, this->signalPlusBackground);
  h->SetBinContent(8, this->signal);
  h->SetBinContent(9, this->background);
  h->SetBinContent(10, this->signal / this->signalPlusBackground);
  h->SetBinContent(11, this->signal / this->background);
  h->SetBinContent(12, this->signal / sqrt(this->signalPlusBackground));

  this->fitResults = (TH1*)h->Clone(saveName.c_str());
  return h;
}

TH1* MassFitter::loadResidualHist() {
  if (!this->data || !this->fit)
    return nullptr;

  string saveName = this->inputs.getSaveNameFromPt("residual");

  TH1* h = (TH1*)this->data->Clone(saveName.c_str());
  TF1* f = (TF1*)this->fit->Clone("function");
  f->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()); // Make sure function is defined everywhere
  h->Add(f, -1.);
  this->residual = (TH1*)h->Clone();
  return h;
}

TH1* MassFitter::loadPullHist() {
  if (!this->residual)
    this->loadResidualHist();

  string saveName = this->inputs.getSaveNameFromPt("pull");

  TH1* h = (TH1*)this->residual->Clone(saveName.c_str());
  h->Reset();
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double bc = this->residual->GetBinContent(i);
    double be = this->data->GetBinError(i);
    double p = bc / be;
    h->SetBinContent(i, p);
  }
  this->pull = (TH1*)h->Clone();
  return h;
}

string MassFitter::setHistNameFromTrain() {
  if (this->inputs.inputIssue("setHistNameFromTrain", "hadron"))
    return "";

  string s = "jet-fragmentation";
  if (this->inputs.train == 287744) s += "_id12406";
  s += "/data/V0/V0PtMass";

  this->inputs.histName = s;
  return s;
}

void MassFitter::setSignalRegionFromSigma(double nSigma) {
  if (this->fitParts.size() == 0) {
    cout << "MassFitter::setSignalRegionFromSigma() Error: fit parts not set" << endl;
    return;
  }

  TF1* f = this->fitParts[0];
  double mu = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  this->inputs.setSignalRegion(mu - nSigma * sigma, mu + nSigma * sigma);
}

void MassFitter::writeOutputsToFile() {
  this->inputs.writeOutputToFile(this->data);
  this->inputs.writeOutputToFile(this->fit);
  for (auto f : this->fitParts) {
    this->inputs.writeOutputToFile(f);
  }
  this->inputs.writeOutputToFile(this->fitResults);
  this->inputs.writeOutputToFile(this->fitParams);
  this->inputs.writeOutputToFile(this->pull);
  this->inputs.writeOutputToFile(this->residual);
}

struct FitSummariser {
  private:
  public:
    InputSettings inputs;

    FitSummariser() { InputSettings(); }
    FitSummariser(InputSettings& x) { this->inputs = x; }

    void summariseFitInfo();
    void summariseFitInfo(string inputHistName);
};

void FitSummariser::summariseFitInfo() {
  this->summariseFitInfo(this->inputs.histName);
}

void FitSummariser::summariseFitInfo(string inputHistName) {
  int nPtBins = this->inputs.ptBinEdges.size();
  int nFitInfos = -1;
  if (inputHistName == "fitParams") {
    TF1* f = new TF1("f", this->inputs.getFitExpression().c_str(), 0., 1.);
    nFitInfos = f->GetNpar();
  } else if (inputHistName == "fitResults") {
    nFitInfos = 12;
  } else {
    cout << "writeFitInfo: unknown hist name " << inputHistName << endl;
    return;
  }

  TH2* fitInfo = new TH2D(inputHistName.c_str(), inputHistName.c_str(), nPtBins, -0.5, nPtBins - 0.5, nFitInfos, -0.5, nFitInfos - 0.5);
  for (int ix = 0; ix < nPtBins; ix++) {
    fitInfo->GetXaxis()->SetBinLabel(ix + 1, TString::Format("%.1f < #it{p}_{T, V0} < %.1f", this->inputs.ptBinEdges[ix][0], this->inputs.ptBinEdges[ix][1]).Data());
  }
  bool binLabelsSet = false;

  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!file) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return;
  }

  for (int iPt = 0; iPt < nPtBins; iPt++) {
    this->inputs.setPt(this->inputs.ptBinEdges[iPt][0], this->inputs.ptBinEdges[iPt][1]);
    string histName = this->inputs.getSaveNameFromPt(inputHistName);
    TH1* hist = (TH1*)file->Get(histName.c_str());
    if (!hist) {
      cout << "Could not find histogram " << histName << " in file " << this->inputs.inputFileName << endl;
      continue;
    }

    for (int iInfo = 0; iInfo < nFitInfos; iInfo++) {
      double parValue = hist->GetBinContent(iInfo + 1);
      fitInfo->SetBinContent(iPt + 1, iInfo + 1, parValue);

      if (binLabelsSet)
        continue;

      fitInfo->GetYaxis()->SetBinLabel(iInfo + 1, hist->GetXaxis()->GetBinLabel(iInfo + 1));
    }
    binLabelsSet = true;
  }

  file->Close();
  this->inputs.writeOutputToFile(fitInfo);
}

struct FitPlotter {
  private:
  public:
    InputSettings inputs;
    MassFitter massFitter;

    TCanvas* canvas = nullptr;
    TH1F* frame = nullptr;
    TLegend* legend = nullptr;

    FitPlotter() { InputSettings(); MassFitter(); }
    FitPlotter(InputSettings& x) { this->inputs = x; this->massFitter = MassFitter(x); }
    FitPlotter(MassFitter& x) { this->massFitter = x; this->inputs = x.inputs; }

    TH1* loadMassHist();
    vector<TF1*> loadFitParts();
    TF1* loadFit();
    void plotFitInfo();
    int plotFitInfo(int infoToPlot);
    void plotFitParts();
};

vector<TF1*> FitPlotter::loadFitParts() {
  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!file) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return {};
  }

  vector<TF1*> parts = {};
  int maxParts = 10; // To prevent infinite loop
  for (int i = 1; i < maxParts; i++) {
    string fName = this->inputs.getSaveNameFromPt("fit", "_f" + to_string(i));
    TF1* f = (TF1*)file->Get(fName.c_str());
    if (!f)
      break;

    parts.push_back(f);
  }

  this->massFitter.fitParts = parts;
  return parts;
}

TF1* FitPlotter::loadFit() {
  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!file) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return nullptr;
  }

  string fName = this->inputs.getSaveNameFromPt("fit");
  TF1* f = (TF1*)file->Get(fName.c_str());
  this->massFitter.fit = (TF1*)f->Clone();
  return f;
}

TH1* FitPlotter::loadMassHist() {
  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!file) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return nullptr;
  }

  string histName = this->inputs.getSaveNameFromPt("data");
  TH1* hist = (TH1*)file->Get(histName.c_str());
  if (!hist) {
    cout << "Could not find histogram " << histName << " in file " << this->inputs.inputFileName << endl;
    return nullptr;
  }
  this->massFitter.data = (TH1*)hist->Clone(histName.c_str());
  return hist;
}

void FitPlotter::plotFitInfo() {
  int maxInfos = 20; // To prevent infinite loop
  for (int i = 0; i < maxInfos; i++) {
    if (this->plotFitInfo(i) != 0)
      break;
  }
}

int FitPlotter::plotFitInfo(int infoToPlot)
{
  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");

  TH2* fitInfo = (TH2*)file->Get(this->inputs.histName.c_str());
  if (!fitInfo) {
    cout << "FitPlotter::plotFitInfo() Error: Could not find histogram " << this->inputs.histName << " in file " << this->inputs.inputFileName << endl;
    return 1;
  }

  if (infoToPlot < 0 || infoToPlot >= fitInfo->GetNbinsY()) {
    cout << "FitPlotter::plotFitInfo() Error: Invalid column " << infoToPlot << " selected for plotting " << this->inputs.histName << endl;
    cout << "Valid range is 0 to " << fitInfo->GetNbinsY() - 1 << endl;
    return 2;
  }

  TH1* hist = (TH1*)fitInfo->ProjectionX("hist", infoToPlot+1, infoToPlot+1);
  setStyle(hist, 0);
  if (!this->canvas) {
    cout << "FitPlotter::plotFitIfo() Warning: Canvas not set, creating automatically" << endl;
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetLogy(this->inputs.logplot);
    this->canvas = canvas;
  }

  this->canvas->cd();
  string histTitle;
  double yMinFrame, yMaxFrame;
  if (this->frame) {
    histTitle = this->frame->GetTitle();
    yMinFrame = this->frame->GetMinimum();
    yMaxFrame = this->frame->GetMaximum();
  } else {
    histTitle = TString::Format("%s, %s, ", getDataSet(this->inputs.train).c_str(), this->inputs.fitName.c_str()).Data();
    if (this->inputs.histName == "fitParams")
      histTitle += "par " + to_string(infoToPlot);
    if (this->inputs.histName == "fitResults")
      histTitle += fitInfo->GetYaxis()->GetBinLabel(infoToPlot + 1);

    yMinFrame = 0.;
    yMaxFrame = 1.1 * getHistUpperBound(hist, false);
  }

  hist->SetTitle(histTitle.c_str());
  hist->SetStats(0);
  hist->SetMinimum(yMinFrame);
  hist->SetMaximum(yMaxFrame);

  string outputFileName = this->inputs.outputFileName;
  if (outputFileName == "") {
    outputFileName = this->inputs.hadron + "_" + this->inputs.fitName + "_" + this->inputs.histName + to_string(infoToPlot) + ".pdf";
    cout << "FitPlotter::plotFitInfo() Warning: outputFileName not set, using: " << outputFileName << endl;
  }
  hist->Draw();
  canvas->SaveAs(outputFileName.c_str());
  return 0;
}

void FitPlotter::plotFitParts() {
  if (this->inputs.outputFileName == "") {
    this->inputs.outputFileName = this->inputs.getSaveNameFromPt(this->inputs.hadron + "_" + this->inputs.fitName, ".pdf");
    cout << "FitPlotter::plotFitParts() Warning: outputFileName not set, using: " << this->inputs.outputFileName << endl;
  }

  if (!this->canvas) {
    cout << "FitPlotter::plotFitParts() Warning: Canvas not set, creating automatically" << endl;
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetLogy(this->inputs.logplot);
    this->canvas = canvas;
  }

  this->canvas->cd();
  if (!this->frame) {
    cout << "FitPlotter::plotFitParts() Warning: Frame not set, creating automatically" << endl;
    double xMinFrame = this->massFitter.data->GetXaxis()->GetXmin();
    double xMaxFrame = this->massFitter.data->GetXaxis()->GetXmax();
    double yMinFrame = 0.;
    double yMaxFrame = 1.1 * getHistUpperBound(this->massFitter.data, false);
    string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(this->inputs.hadron).c_str()).Data();
    string yTitle = "";
    string histTitle = TString::Format("%s, (%.1f < #it{p}_{T, V0} < %.1f)", getDataSet(this->inputs.train).c_str(), this->inputs.lowpt, this->inputs.highpt).Data();

    TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
    frame->SetTitle(histTitle.c_str());
    this->frame = frame;
  }

  this->frame->Draw();
  if (this->legend)  this->legend->Draw("same");

  for (int i = 0; i < this->massFitter.fitParts.size(); i++) {
    TF1* f = this->massFitter.fitParts[i];
    setStyle(f, i + 2);
    f->Draw("same");
  }
  // Only do this if auto setup is requested
  setStyle(this->massFitter.data, 0);
  setStyle(this->massFitter.fit, 1);
  this->massFitter.fit->SetRange(this->massFitter.data->GetXaxis()->GetXmin(), this->massFitter.data->GetXaxis()->GetXmax());

  this->massFitter.data->Draw("same");
  this->massFitter.fit->Draw("same");
  this->canvas->SaveAs(this->inputs.outputFileName.c_str());
}

// -------------------------------------------------------------------------------------------------
//
// Interface
//
// -------------------------------------------------------------------------------------------------

// This runs the entire workflow in order
void quickRun() {
  double nSigma = 3.;
  InputSettings x;
  x.train = 252064;
  x.hadron = "K0S";
  x.fitName = "pol1BreitWigner";
  x.normaliseData = true;

  x.setFitX(0.45, 0.55);
  x.setPolInitX(0.45, 0.55, 0.57);
  x.setMassWindowDiff(2e-2, 2e-2);
  x.setInputFileNameFromTrain();
  x.outputFileName = x.hadron + "_" + x.fitName + ".root";

  vector<vector<double>> ptBinEdges;
  ptBinEdges.push_back({1., 2.});
  ptBinEdges.push_back({2., 3.});
  ptBinEdges.push_back({3., 4.});
  ptBinEdges.push_back({4., 5.});
  ptBinEdges.push_back({5., 10.});

  // For each pt bin, do fit, and plot fit parts
  for (int iPt = 0; iPt < ptBinEdges.size(); iPt++) {
    x.setPt(ptBinEdges[iPt][0], ptBinEdges[iPt][1]);
    MassFitter m(x);
    m.loadMassHist();
    m.loadFitFunction();
    m.setFitInitialValues();
    m.data->Fit(m.fit, "R");
    m.loadFitParts();
    m.loadFitParams();
    m.setSignalRegionFromSigma(nSigma);
    m.loadFitResults();
    m.loadResidualHist();
    m.loadPullHist();
    m.writeOutputsToFile();

    FitPlotter p(m);
    p.inputs.outputFileName = x.getSaveNameFromPt(x.hadron + "_" + x.fitName, ".pdf");
    p.plotFitParts();
  }

  // Compile the information into summary histograms
  x.setInputFileNameFromFit();
  x.outputFileName = x.inputFileName;
  x.setPtBinEdgesSorted(ptBinEdges);
  FitSummariser f(x);
  f.summariseFitInfo("fitParams");
  f.summariseFitInfo("fitResults");

  // Plot the evolution of fit parameters and results with pt
  x.outputFileName = "";
  FitPlotter p(x);
  p.inputs.histName = "fitParams";
  p.plotFitInfo();
  p.inputs.histName = "fitResults";
  p.plotFitInfo();
}

InputSettings quickSetup(double ptmin, double ptmax) {
  InputSettings x;
  x.train = 252064;
  x.hadron = "K0S";
  x.fitName = "pol1BreitWigner";
  x.normaliseData = true;
  x.setPt(ptmin, ptmax);
  x.setFitX(0.45, 0.55);
  x.setPolInitX(0.45, 0.55, 0.57);
  x.setMassWindowDiff(2e-2, 2e-2);

  x.setInputFileNameFromTrain();
  x.outputFileName = x.hadron + "_" + x.fitName + ".root";

  return x;
}

InputSettings quickSetup() {
  vector<vector<double>> ptBinEdges;
  ptBinEdges.push_back({1., 2.});
  ptBinEdges.push_back({5., 10.});
  ptBinEdges.push_back({3., 4.});
  ptBinEdges.push_back({2., 3.});
  ptBinEdges.push_back({4., 5.});

  InputSettings x;
  x.train = 252064;
  x.hadron = "K0S";
  x.fitName = "pol1BreitWigner";
  x.setInputFileNameFromFit();
  x.setPtBinEdgesSorted(ptBinEdges);

  cout << "Input file: " << x.inputFileName << endl;

  return x;
}

MassFitter doFitting(InputSettings& x, double nSigma = 3.) {
  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.setSignalRegionFromSigma(nSigma);
  m.loadFitResults();
  m.loadResidualHist();
  m.loadPullHist();
  m.writeOutputsToFile();
  return m;
}

void plotFitParts(InputSettings& x) {
  FitPlotter p(x);
  p.loadMassHist();
  p.loadFit();
  p.loadFitParts();
  p.plotFitParts();
}

void plotFitParts(MassFitter& m) {
  FitPlotter p(m);
  p.inputs.outputFileName = p.inputs.getSaveNameFromPt(p.inputs.hadron + "_" + p.inputs.fitName, ".pdf");
  p.plotFitParts();
}

void summariseFitInfo() {
  InputSettings x = quickSetup();
  x.outputFileName = x.inputFileName;
  cout << "Input file: " << x.inputFileName << "\nOutput file: " << x.outputFileName << endl;
  FitSummariser f(x);
  f.summariseFitInfo("fitParams");
  f.summariseFitInfo("fitResults");
}

void plotFitInfo() {
  InputSettings x = quickSetup();
  x.outputFileName = "";
  FitPlotter p(x);
  p.inputs.histName = "fitParams";
  p.plotFitInfo();
  p.inputs.histName = "fitResults";
  p.plotFitInfo();
}

#endif
