
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

// Crystal ball function, from: https://en.wikipedia.org/wiki/Crystal_Ball_function
double CrystalBall(double *x, double *p)
{
  double amp   = p[0];
  double alpha = p[1];
  double x0    = p[2];
  double sigma = p[3];
  double n     = p[4];

  double t = (x[0] - x0) / sigma;
  if (t > -1. * alpha) {
    return amp * exp(-0.5 * t * t);
  }
  else {
    double a = pow(n / abs(alpha), n) * exp(-0.5 * alpha * alpha);
    double b = n / abs(alpha) - abs(alpha);
    return amp * a / pow(b - t, n);
  }
}

// Prints the parameter limits of a function
void printParLimits(TF1* f)
{
  for (int i = 0; i < f->GetNpar(); i++) {
    double min, max;
    f->GetParLimits(i, min, max);
    cout << f->GetName() << "(" << i << " = " << f->GetParName(i) << ") " << f->GetParameter(i) << " (" << min << ", " << max << ")" << endl;
  }
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
TF1* combineTFs(vector<TF1*> w)
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
  return F;
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

// -------------------------------------------------------------------------------------------------
//
// Functions to use with setup
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
  f->SetParLimits(2, 0.5 * sigma, 2. * sigma);

  for (int i = 0; i < f->GetNpar(); i++) {
    string parName = uniqueName + "_" + to_string(i);
    f->SetParName(i, parName.c_str());
  }
  return f;
}

// Linear background, give it two x values to which it will initialize the fit
TF1* linear(TH1* h, string uniqueName, array<double, 2> vals, array<double, 2> fitRange = {1.08, 1.2})
{
  array<int, 2> region = getProjectionBins(h->GetXaxis(), vals[0], vals[1]);
  double dx = h->GetBinCenter(region[1]) - h->GetBinCenter(region[0]);
  double dy = h->GetBinContent(region[1]) - h->GetBinContent(region[0]);

  double b = dy / dx;
  double a = h->GetBinContent(region[0]) - b * h->GetBinCenter(region[0]);

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

// TF1* landau(TH1* h, string uniqueName, )

// -------------------------------------------------------------------------------------------------
//
// Functions to use with setup
//
// -------------------------------------------------------------------------------------------------

// My "favourite" histogram
TH1* getHist(double ptmin, double ptmax, string hadron)
{
  string inName = "~/cernbox/TrainOutput/252064/AnalysisResults.root";
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

void myDraw(TCanvas* canvas, TH1* hist, vector<TF1*> fits, TLegend* legend, vector<TLatex*> latex)
{
  canvas->Divide(2,1);
  canvas->cd(1);
  gPad->SetRightMargin(0.05);

  TH1F* frame = DrawFrame(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(),
                          0, 1.1 * hist->GetBinContent(hist->GetMaximumBin()),
                          hist->GetXaxis()->GetTitle(), hist->GetYaxis()->GetTitle());
  frame->Draw();

  hist->Draw("same");
  legend->AddEntry(hist, hist->GetName(), "l");

  for (int i = 0 ; i < fits.size(); i++) {
    legend->AddEntry(fits[i], fits[i]->GetName(), "l");
    fits[i]->Draw("same");
  }

  canvas->cd(2);
  gPad->SetLeftMargin(0.);
  legend->Draw("same");
  for (auto l : latex) {
    l->Draw("same");
  }

  string saveName = canvas->GetName();
  canvas->SaveAs(saveName.c_str());
}

// Input to add: hadron, rebin
// Plots multiple fits for signal+bkg
void plotBkgs(double ptmin, double ptmax, string hadron)
{
  gStyle->SetNdivisions(505, "xy");
  // string hadron = "Lambda0";
  string dataSet = "LHC22o_pass6";


  double textSize = 0.04;
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";

  TH1D* h = (TH1D*)getHist(ptmin, ptmax, hadron);
  h->SetStats(0);
  setStyle(h, 0);
  h->SetXTitle(xTitle.c_str());
  h->SetYTitle(yTitle.c_str());
  h->SetName("data");

  array<int, 2> ptBins = getProjectionBins(h->GetXaxis(), ptmin, ptmax);
  double lowpt  = h->GetXaxis()->GetBinLowEdge(ptBins[0]);
  double highpt = h->GetXaxis()->GetBinUpEdge(ptBins[1]);

  string canvasName = hadron;
  canvasName += TString::Format("_pt%.1f-%.1f", ptmin, ptmax).Data();
  canvasName += ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1800, 900);

  double xMinLegend = 0., xMaxLegend = 0.7, yMinLegend = 0.25, yMaxLegend = 0.45;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend);

  // array<double, 2> fitRange = {h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()};
  array<double, 2> fitRange = {1.08, 1.2};
  double mass = MassLambda0;
  double dM = 1e-2;
  double signalWidth = 1e-3;
  if ("K0S" == hadron) {
    fitRange = {0.4, 0.6};
    mass = MassK0S;
    dM = 5e-2;
    signalWidth = 1e-2;
  }

  // Functions to combine
  array<int, 2> peakRegion = getProjectionBins(h->GetXaxis(), mass - dM, mass + dM);
  TH1* k = makeHistSubset(h, peakRegion[0], peakRegion[1]);
  double peakVal = k->GetBinContent(k->GetMaximumBin());
  array<double, 3> sigGausVals = {peakVal / 1.5, mass, signalWidth};
  array<double, 2> linVals     = {sigGausVals[0] - 3*sigGausVals[2], sigGausVals[0] + 3*sigGausVals[2]};
  array<double, 3> sigmoidVals = {h->GetBinContent(h->GetNbinsX()), h->GetBinContent(h->GetNbinsX())/20., sigGausVals[0]};
  array<double, 3> dilogVals   = {h->GetBinContent(h->GetNbinsX()), 35., linVals[0]};
  array<double, 3> gausVals    = {sigGausVals[0] / 10., mass, signalWidth * (10. - ("K0S" == hadron) * 5.)};

  // TF1* fSigGaus = sigGaus(h, "fSigGaus", sigGausVals, fitRange);
  TF1* fSigGaus = gaussian(h, "fSigGaus", sigGausVals, fitRange);
  TF1* fLin = linear(h, "fLin", linVals, fitRange);
  TF1* fSigmoid = sigmoid(h, "fSigmoid", sigmoidVals, fitRange);
  TF1* fDilog = dilog(h, "fDilog", dilogVals, fitRange);
  TF1* fGaus = gaussian(h, "fGaus", gausVals, fitRange);

  vector<TF1*> fits;
  TF1* fGL = combineTFs({fSigGaus, fLin});
  fGL->SetName("+pol1");
  fGL->SetLineColor(GetColor(1));
  fGL->SetLineWidth(3);
  fits.push_back(fGL);

  TF1* fLan = new TF1("fLan", "[0]*TMath::Landau(x, [1], [2])", fitRange[0], fitRange[1]);
  fLan->SetParNames("Amp", "Mean", "Sigma");
  fLan->SetParameters(peakVal, mass, signalWidth);
  fLan->SetParLimits(0, 0.1 * peakVal, 2. * peakVal);
  fLan->SetParLimits(1, mass - dM, mass + dM);
  fLan->SetParLimits(2, 1e-6, 3. * signalWidth);
  fLan->SetName("Landau");
  fLan->SetLineWidth(3);
  fLan->SetLineColor(GetColor(2));
  fits.push_back(fLan);

  // TF1* fGS = combineTFs({fSigGaus, fSigmoid});
  // fGS->SetName("+#sigma(x)");
  // fGS->SetLineColor(GetColor(1));
  // fits.push_back(fGS);

  TF1* fGGS = combineTFs({fSigGaus, fGaus, fSigmoid});
  fGGS->SetName("+G+#sigma(x)");
  fGGS->SetLineColor(GetColor(2));
  // fits.push_back(fGGS);

  TF1* fGD = combineTFs({fSigGaus, fDilog});
  fGD->SetName("+Dilog");
  fGD->SetLineColor(GetColor(3));
  fits.push_back(fGD);

  TF1* fGGD = combineTFs({fSigGaus, fGaus, fDilog});
  fGGD->SetName("+G+Dilog");
  fGGD->SetLineColor(GetColor(4));
  // fits.push_back(fGGD);

  // TF1* fCb1 = new TF1("fCb1", CrystalBall, fitRange[0],fitRange[1], 5);
  // fCb1->SetParNames("Amp", "Alpha", "Mean", "Sigma", "N");
  // fCb1->SetParameters(peakVal, -1., MassLambda0, 0.01, 1.);
  // fCb1->SetLineColor(GetColor(1));
  // for (int i = 0; i < fCb1->GetNpar(); i++) {
  //   double p = abs(fCb1->GetParameter(i));
  //   fCb1->SetParLimits(i, -20 * p, 20 * p);
  // }
  // fits.push_back(fCb1);

  // TF1* fCb2 = new TF1("fCb2", CrystalBall, fitRange[0],fitRange[1], 5);
  // fCb2->SetParNames("Amp", "Alpha", "Mean", "Sigma", "N");
  // fCb2->SetParameters(peakVal, -2., MassLambda0, 0.01, 1.);
  // fCb2->SetLineColor(GetColor(2));
  // for (int i = 0; i < fCb2->GetNpar(); i++) {
  //   double p = abs(fCb2->GetParameter(i));
  //   fCb2->SetParLimits(i, -20 * p, 20 * p);
  // }
  // fits.push_back(fCb2);

  // TF1* fCb3 = new TF1("fCb3", CrystalBall, fitRange[0],fitRange[1], 5);
  // fCb3->SetParNames("Amp", "Alpha", "Mean", "Sigma", "N");
  // fCb3->SetParameters(peakVal, -0.5, MassLambda0, 0.01, 1.);
  // fCb3->SetLineColor(GetColor(3));
  // for (int i = 0; i < fCb3->GetNpar(); i++) {
  //   double p = abs(fCb3->GetParameter(i));
  //   fCb3->SetParLimits(i, -20 * p, 20 * p);
  // }
  // fits.push_back(fCb3);

  // TF1* fCb4 = new TF1("fCb4", CrystalBall, fitRange[0],fitRange[1], 5);
  // fCb4->SetParNames("Amp", "Alpha", "Mean", "Sigma", "N");
  // fCb4->SetParameters(peakVal, -1.5, MassLambda0, 0.01, 1.);
  // fCb4->SetLineColor(GetColor(4));
  // for (int i = 0; i < fCb4->GetNpar(); i++) {
  //   double p = abs(fCb4->GetParameter(i));
  //   fCb4->SetParLimits(i, -20 * p, 20 * p);
  // }
  // fits.push_back(fCb4);


  for (auto f : fits) {
    f->SetLineWidth(3);
    h->Fit(f, "RSBQ0");
    string newName = TString::Format("%s (#chi^{2}/NDF = %.1f)", f->GetName(), f->GetChisquare() / f->GetNDF()).Data();
    f->SetName(newName.c_str());

    printParLimits(f);
    cout << endl;
  }

  vector<TLatex*> latex;
  latex.push_back(
    CreateLatex(0.1, 0.9, dataSet.c_str(), textSize));
  latex.push_back(
    CreateLatex(0.1, 0.85, TString::Format("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptmin, ptmax).Data(), textSize));

  myDraw(canvas, h, fits, legend, latex);
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

void plotBkgParts(double ptmin, double ptmax, string hadron = "Lambda0")
{
  gStyle->SetNdivisions(505, "xy");
  string dataSet = "LHC22o_pass6";

  double textSize = 0.04;
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";

  TH1D* h = (TH1D*)getHist(ptmin, ptmax, hadron);
  h->SetStats(0);
  setStyle(h, 0);
  h->SetXTitle(xTitle.c_str());
  h->SetYTitle(yTitle.c_str());
  h->SetName("data");
  h->SetTitle("");

  array<int, 2> ptBins = getProjectionBins(h->GetXaxis(), ptmin, ptmax);
  double lowpt  = h->GetXaxis()->GetBinLowEdge(ptBins[0]);
  double highpt = h->GetXaxis()->GetBinUpEdge(ptBins[1]);

  string canvasName = hadron;
  canvasName += TString::Format("_pt%.1f-%.1f", ptmin, ptmax).Data();

  double xMinLegend = 0.5, xMaxLegend = 0.8, yMinLegend = 0.7, yMaxLegend = 0.85;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, "" , textSize);

  array<double, 2> fitRange = {1.09, 1.15};
  double mass = MassLambda0;
  double dM = 1e-2;
  double signalWidth = 1e-3;
  if ("K0S" == hadron) {
    fitRange = {0.4, 0.6};
    mass = MassK0S;
    dM = 5e-2;
    signalWidth = 1e-2;
  }

  // Functions to combine
  array<int, 2> peakRegion = getProjectionBins(h->GetXaxis(), mass - dM, mass + dM);
  TH1* k = makeHistSubset(h, peakRegion[0], peakRegion[1]);
  double peakVal = k->GetBinContent(k->GetMaximumBin());
  array<double, 3> sigGausVals = {peakVal / 1.5, mass, signalWidth};
  array<double, 2> linVals     = {sigGausVals[0] - 3*sigGausVals[2], sigGausVals[0] + 3*sigGausVals[2]};
  array<double, 3> sigmoidVals = {h->GetBinContent(h->GetNbinsX()), h->GetBinContent(h->GetNbinsX())/20., sigGausVals[0]};
  array<double, 3> dilogVals   = {peakVal / 10., 35., 1.1};
  array<double, 3> gausVals    = {sigGausVals[0] / 10., mass, signalWidth * (10. - ("K0S" == hadron) * 5.)}; // Change width here?

  // Background functions start with fx to ensure signal is first in alphabetical order
  TF1* fSigGaus = gaussian(h, "fSigGaus", sigGausVals, fitRange);
  TF1* fLin     = linear(h,   "fxLin", linVals, fitRange);
  TF1* fSigmoid = sigmoid(h,  "fxSigmoid", sigmoidVals, fitRange);
  TF1* fDilog   = dilog(h,    "fxDilog", dilogVals, fitRange);
  TF1* fGaus    = gaussian(h, "fxGaus", gausVals, fitRange);

  // vector<TF1*> functions = {fSigGaus, fSigmoid}; string fName = "+#sigma(x)";
  // vector<TF1*> functions = {fSigGaus, fGaus, fSigmoid}; string fName = "+G+#sigma(x)";
  vector<TF1*> functions = {fSigGaus, fDilog}; string fName = "+Dilog";
  // vector<TF1*> functions = {fSigGaus, fGaus, fDilog}; string fName = "+G+Dilog";

  functions = sortAlphabetically(functions);
  TF1* f = combineTFs(functions);
  canvasName += TString::Format("_fit=%s", f->GetName()).Data();

  f->SetLineWidth(3);
  h->Fit(f, "RSBQ0");
  string newName = TString::Format("%s (#chi^{2}/NDF = %.1f)", fName.c_str(), f->GetChisquare() / f->GetNDF()).Data();
  f->SetName(newName.c_str());
  legend->AddEntry(f, "Fit", "l");
  legend->SetHeader(TString::Format("#chi^{2}/NDF = %.1f", f->GetChisquare() / f->GetNDF()).Data());
  splitTFs(f, functions);

  f->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  for (int i = 0; i < functions.size(); i++) {
    functions[i]->SetLineColor(GetColor(i+2));
    functions[i]->SetLineWidth(3);
    legend->AddEntry(functions[i], functions[i]->GetName(), "l");
    functions[i]->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  }

  vector<TLatex*> latex;
  latex.push_back(
    CreateLatex(0., 0.95, dataSet.c_str(), textSize));
  latex.push_back(
    CreateLatex(0.5, 0.95, TString::Format("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", ptmin, ptmax).Data(), textSize));

  canvasName += ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);
  plotFitParts(canvas, h, f, functions, legend, latex);
}
