
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
    cout << f->GetName() << " (" << i << " = " << f->GetParName(i) << ") " << f->GetParameter(i) << " (" << min << ", " << max << ")" << endl;
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
  h->Add(fit, -1.);
  return h;
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

  cout << "fit range: " << fitRange[0] << " - " << fitRange[1] << endl;
  cout << "values: " << vals[0] << " - " << vals[1] << endl;
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

// Input to add: hadron, rebin
// Plots multiple fits for signal+bkg
void plotBkgs(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  TH1D* hist = (TH1D*)getHist(ptmin, ptmax, hadron, inName);
  hist->Sumw2();
  setStyle(hist, 0);
  // hist->Scale(1./hist->Integral(), "width");

  string saveName = hadron;
  saveName += "_";
  saveName += hist->GetName(); // hist name contains pt range
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 1800, 900);
  // hist->SetName("data");

  array<double, 2> fitRange = {1.08, 1.2};
  double mass = MassLambda0;
  double dM = 1e-2;
  double signalWidth = 1e-3;
  if ("K0S" == hadron) {
    fitRange = {0.45, 0.55};
    mass = MassK0S;
    dM = 5e-2;
    signalWidth = 1e-2;
  }

  // Functions to combine
  array<int, 2> peakRegion = getProjectionBins(hist->GetXaxis(), mass - dM, mass + dM);
  TH1* k = makeHistSubset(hist, peakRegion[0], peakRegion[1]);
  double peakVal = k->GetBinContent(k->GetMaximumBin());
  array<double, 3> sigGausVals = {peakVal / 1.5, mass, signalWidth};
  array<double, 2> linVals     = {sigGausVals[0] - 3*sigGausVals[2], sigGausVals[0] + 3*sigGausVals[2]};
  array<double, 3> sigmoidVals = {hist->GetBinContent(hist->GetNbinsX()), hist->GetBinContent(hist->GetNbinsX())/20., sigGausVals[0]};
  array<double, 3> dilogVals   = {hist->GetBinContent(hist->GetNbinsX()), 35., linVals[0]};
  array<double, 3> gausVals    = {sigGausVals[0] / 10., mass, signalWidth * (10. - ("K0S" == hadron) * 5.)};

  // TF1* fSigGaus = sigGaus(h, "fSigGaus", sigGausVals, fitRange);
  TF1* fSigGaus = gaussian(hist, "fSigGaus", sigGausVals, fitRange);
  // TF1* fLin = linear(hist, "fLin", linVals, fitRange);
  TF1* fSigmoid = sigmoid(hist, "fSigmoid", sigmoidVals, fitRange);
  TF1* fDilog = dilog(hist, "fDilog", dilogVals, fitRange);
  TF1* fGaus = gaussian(hist, "fGaus", gausVals, fitRange);

  vector<TF1*> fits;

  TF1* fS = fSigGaus;
  fS->SetName("SigGaus");
  // fits.push_back(fS);

  // TF1* fG = fSigGaus;
  TF1* fG = fGaus;
  fG->SetName("G");
  // fits.push_back(fG);

  TF1* fL = new TF1("fL", "[0] + [1]*x", fitRange[0], fitRange[1]);
  array<int, 2> vals = getProjectionBins(hist->GetXaxis(), fitRange[0], fitRange[1]);
  double dx = hist->GetBinCenter(vals[1]) - hist->GetBinCenter(vals[0]);
  double dy = hist->GetBinContent(vals[1]) - hist->GetBinContent(vals[0]);
  cout << "fit range: " << fitRange[0] << " " << fitRange[1] << endl;
  cout << "peak bins: " << peakRegion[0] << " " << peakRegion[1] << endl;
  cout << "bins: " << vals[0] << " " << vals[1] << endl;
  cout << "x: " << hist->GetBinCenter(vals[0]) << " " << hist->GetBinCenter(vals[1]) << endl;
  cout << "y: " << hist->GetBinContent(vals[0]) << " " << hist->GetBinContent(vals[1]) << endl;
  cout << "dx: " << dx << ", dy: " << dy << endl;
  double b = 1.;
  double a = 100.;
  // double b = dy / dx;
  // double a = hist->GetBinContent(vals[0]) - b * hist->GetBinCenter(vals[0]);
  fL->SetParameter(0, a);
  fL->SetParameter(1, b);
  fL->SetParLimits(0, -2. * (abs(a) + 1), 2. * (abs(a) + 1));
  fL->SetParLimits(1, -2. * (abs(b) + 1), 2. * (abs(b) + 1));

  for (int i = 0; i < fL->GetNpar(); i++) {
    string parName = "fL_" + to_string(i);
    fL->SetParName(i, parName.c_str());
  }
  fL->SetName("L");
  // fits.push_back(fL);

  // TF1* fGL = combineTFs({fG, fL});
  TF1* fGL = combineTFs({fL, fG});
  fGL->SetName("G+L");
  fits.push_back(fGL);
  // printParLimits(fGL);
  // return;

  // TF1* fSL = combineTFs({fS, fL});
  // fSL->SetName("+pol1");
  // fits.push_back(fSL);
  // printParLimits(fSL);

  // TF1* fLan = new TF1("fLan", "[0]*TMath::Landau(x, [1], [2])", fitRange[0], fitRange[1]);
  // fLan->SetParNames("Amp", "Mean", "Sigma");
  // fLan->SetParameters(peakVal, mass, signalWidth);
  // fLan->SetParLimits(0, 0.1 * peakVal, 2. * peakVal);
  // fLan->SetParLimits(1, mass - dM, mass + dM);
  // fLan->SetParLimits(2, 1e-6, 3. * signalWidth);
  // fLan->SetName("Landau");
  // fLan->SetLineWidth(3);
  // fLan->SetLineColor(GetColor(2));
  // fits.push_back(fLan);

  // TF1* fGS = combineTFs({fSigGaus, fSigmoid});
  // fGS->SetName("+#sigma(x)");
  // fGS->SetLineColor(GetColor(1));
  // fits.push_back(fGS);

  // TF1* fGGS = combineTFs({fSigGaus, fGaus, fSigmoid});
  // fGGS->SetName("+G+#sigma(x)");
  // fGGS->SetLineColor(GetColor(2));
  // fits.push_back(fGGS);

  // TF1* fGD = combineTFs({fSigGaus, fDilog});
  // fGD->SetName("+Dilog");
  // fGD->SetLineColor(GetColor(3));
  // fits.push_back(fGD);

  // TF1* fGGD = combineTFs({fSigGaus, fGaus, fDilog});
  // fGGD->SetName("+G+Dilog");
  // fGGD->SetLineColor(GetColor(4));
  // fits.push_back(fGGD);

  TLegend* legend = CreateLegend(0.25, 0.8, 0.7, 0.9);
  for (int i = 0; i < fits.size(); i++) {
    fits[i]->SetLineWidth(3);
    fits[i]->SetLineColor(GetColor(i+1));
    // printParLimits(fits[i]);
    hist->Fit(fits[i], "RSBQ0");
    fits[i]->SetRange(hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    string newName = TString::Format("%s (#chi^{2}/NDF = %.1f)", fits[i]->GetName(), fits[i]->GetChisquare() / fits[i]->GetNDF()).Data();
    fits[i]->SetName(newName.c_str());
    legend->AddEntry(fits[i], newName.c_str(), "l");

    // printParLimits(fits[i]);
    cout << endl;
  }

  // TF1* fGL = combineTFs({fits[0], fits[1]});
  // fGL->SetName("+pol1");
  // fits.push_back(fGL);
  // printParLimits(fits[0]);
  // printParLimits(fits[1]);

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", ptmin, ptmax).Data();
  string xTitle = "#it{M}(" + formatHadronDaughters(hadron) + ") (GeV/#it{c}^{2})";
  string yTitle = "counts";
  double xMinFrame = hist->GetXaxis()->GetXmin(), xMaxFrame = hist->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.5 * getHistScale(hist, false);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  canvas->cd();
  frame->Draw();
  legend->Draw("same");
  hist->Draw("same");
  for (auto f : fits) {
    f->Draw("same");
  }
  canvas->SaveAs(saveName.c_str());

  TCanvas* canvas2 = new TCanvas("residuals", "residuals", 1800, 900);
  canvas2->cd();
  TH1* residual = makeResidual(hist, fits[0]);
  residual->Fit(fS, "RSBQ0");
  residual->SetStats(0);
  residual->SetXTitle(xTitle.c_str());
  residual->SetTitle("Residual of fit");
  residual->Draw();
  fS->Draw("same");
  canvas2->SaveAs("residuals.pdf");
}

void plotBkgParts(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron = inputStrings[2];

  double textSize = 0.04;
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";

  TH1D* h = (TH1D*)getHist(ptmin, ptmax, hadron, 252064);
  h->SetStats(0);
  setStyle(h, 0);
  h->SetXTitle(xTitle.c_str());
  h->SetYTitle(yTitle.c_str());
  h->SetName("data");
  h->SetTitle("");

  array<int, 2> ptBins = getProjectionBins(h->GetXaxis(), ptmin, ptmax);
  double lowpt  = h->GetXaxis()->GetBinLowEdge(ptBins[0]);
  double highpt = h->GetXaxis()->GetBinUpEdge(ptBins[1]);

  string saveName = hadron;
  saveName += TString::Format("_pt%.1f-%.1f", ptmin, ptmax).Data();

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
  saveName += TString::Format("_fit=%s", f->GetName()).Data();

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

  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);
  plotFitParts(canvas, h, f, functions, legend, latex);
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
      plotBkgs(inputStrings, v0min, v0max);
      break;
    case 1:
      plotBkgParts(inputStrings, v0min, v0max);
      break;
  }
}

void plot252064(string hadron, double v0min, double v0max, int setting) { plotTrain(252064, hadron, v0min, v0max, setting); }
void plot282430(string hadron, double v0min, double v0max, int setting) { plotTrain(282430, hadron, v0min, v0max, setting); }

void test() {
  gROOT->SetBatch(kTRUE);
  plot252064("K0S", 5., 10., 0);
}