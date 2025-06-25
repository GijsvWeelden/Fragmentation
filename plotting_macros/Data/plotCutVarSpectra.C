
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

// Cut variation: See how the mass spectrum changes when varying the V0 variable cuts

#ifndef CUTVARIATION_C
#define CUTVARIATION_C

namespace MyEnums {
  enum Verbosity {kErrors, kWarnings, kInfo, kDebug};
  enum FitType   {kExpGausExp, kGausGaus, kGausGausExp};
  enum Variable  {kR, kCtau, kCosPA, kDCAp, kDCAn, kDCAd};
  enum CutType   {kLower, kUpper, kBoth};
}

using namespace MyEnums;

struct InputSettings {
  private:
  public:
    int train, rebinNumber, nHists;
    int fitType = -1, var = -1, verbosity = kWarnings, cutType = -1;
    string inputFileNameHist, inputFileNameFit, outputFileName;
    string hadron, fitName, histName, varName, varSymbol;

    vector<vector<double>> ptBinEdges = {};
    vector<double> varLimits = {};

    double ptmin, ptmax, lowpt, highpt;
    double varmin, varmax, lowvar, highvar;
    double fitmin, fitmax;
    double massWindowMin, massWindowMax, signalRegionMin, signalRegionMax;

    bool logplot = false;

    double getMass();
    string getSaveNameFromPt(string prefix, string suffix);
    string getSaveNameFromVar(string prefix, string suffix);
    string getVarSymbol(int xMaxLegend);
    string printLog(string message, int verbThreshold);
    string setFitName(int fit);
    int setFitType(string fit);
    void setFitX(double a, double b);
    string setInputFileNameFromTrain();
    string setInputFileNameFromFit();
    void setMassWindow(double a, double b);
    void setPt(double a, double b);
    vector<vector<double>> setPtBinEdges();
    void setSignalRegion(); // Auto-setup based on table. To be implemented
    void setSignalRegion(double a, double b);
    string setVar(int v);
    void setVarLimits(double a, double b);
    template <typename T> int writeOutputToFile(T* obj);
};


double InputSettings::getMass() {
  if (hadron == "K0S")
    return MassK0S;
  if (hadron == "Lambda" || hadron == "AntiLambda")
    return MassLambda;

  return -1.;
}

string InputSettings::getSaveNameFromPt(string prefix, string suffix = "") {
  string s = TString::Format("%s_pt%.1f-%.1f%s", prefix.c_str(), lowpt, highpt, suffix.c_str()).Data();
  return s;
}

string InputSettings::getSaveNameFromVar(string prefix, string suffix = "") {
  string s = TString::Format("%s_%s%.1f-%.1f%s", prefix.c_str(), varName.c_str(), lowvar, highvar, suffix.c_str()).Data();
  return s;
}

string InputSettings::getVarSymbol(int x = -1) {
  if (x < 0) x = var;
  string v = "";
  switch (x) {
    case kR:
      v = "R";
      break;
    case kCtau:
      v = "c#tau";
      break;
    case kCosPA:
      v = "cos(PA)";
      break;
    case kDCAp:
      v = "DCA pos";
      break;
    case kDCAn:
      v = "DCA neg";
      break;
    case kDCAd:
      v = "DCA dau";
      break;
    default:
      string s = "InputSettings::getVarSymbol() Error: requested unknown variable";
      printLog(s, kErrors);
  }
  return v;
}

string InputSettings::printLog(string message, int verbThreshold) {
  if (verbosity < verbThreshold)
    return "";

  cout << message << endl;
  return message;
}

string InputSettings::setFitName(int fit = -1) {
  if (fit >= 0) fitType = fit;

  switch (fitType) {
    case kExpGausExp:
      fitName = "ExpGausExp";
      break;
    case kGausGaus:
      fitName = "GausGaus";
      break;
    case kGausGausExp:
      fitName = "GausGausExp";
      break;
    default:
      string s = "InputSettings::setFitName() Error: requested unknown function";
      printLog(s, kErrors);
  }
  return fitName;
}

int InputSettings::setFitType(string fit = "") {
  if (fit != "") fitName = fit;

  if (fitName == "ExpGausExp") {
    fitType = kExpGausExp;
  } else if (fitName == "GausGaus") {
    fitType = kGausGaus;
  } else if (fitName == "GausGausExp") {
    fitType = kGausGausExp;
  } else {
    string s = "InputSettings::getFitType() Error: requested unknown function " + fitName;
    printLog(s, kErrors);
  }

  return fitType;
}

void InputSettings::setFitX(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setFitX() Error: fitmin > fitmax";
    printLog(s, kErrors);
    return;
  }
  fitmin = a;
  fitmax = b;
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  inputFileNameHist = s;
  return s;
}

string InputSettings::setInputFileNameFromFit() {
  string s = hadron + "_pol1" + fitName + ".root";
  inputFileNameFit = s;
  return s;
}

void InputSettings::setMassWindow(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setMassWindow() Error: massWindowMin > massWindowMax";
    printLog(s, kErrors);
    return;
  }
  massWindowMin = a;
  massWindowMax = b;
}

void InputSettings::setPt(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setPt() Error: ptmin > ptmax";
    printLog(s, kErrors);
    return;
  }
  ptmin = a;
  ptmax = b;
  lowpt = a;
  highpt = b;
}

vector<vector<double>> InputSettings::setPtBinEdges() {
  if (fitType == kGausGaus) {
    vector<vector<double>> x = {{0., 1.}, {1., 2.}, {2., 3.}, {3., 4.}, {4., 5.}, {5., 10.}};
    ptBinEdges = x;
  } else if (fitType == kGausGausExp || fitType == kExpGausExp) {
    vector<vector<double>> x = {{10., 15.}, {15., 20.}, {20., 25.}, {25., 30.}, {30., 40.}};
    ptBinEdges = x;
  }
  return ptBinEdges;
}

void InputSettings::setSignalRegion(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setSignalRegion() Error: signalRegionMin > signalRegionMax";
    printLog(s, kErrors);
    return;
  }
  signalRegionMin = a;
  signalRegionMax = b;
}

string InputSettings::setVar(int v) {
  string name;
  switch (v) {
    case kR:
      name = "Radius";
      break;
    case kCtau:
      name = "Ctau";
      break;
    case kCosPA:
      name = "CosPA";
      break;
    case kDCAp:
      name = "DCApos";
      break;
    case kDCAn:
      name = "DCAneg";
      break;
    case kDCAd:
      name = "DCAd";
      break;
    default:
      {
        string s = "InputSettings::setVar() Error: requested unknown var " + to_string(v);
        printLog(s, kErrors);
      }
  }
  var = v;
  varName = name;
  varSymbol = getVarSymbol();
  return varName;
}

void InputSettings::setVarLimits(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setVarLimits() Error: varmin > varmax";
    printLog(s, kErrors);
    return;
  }
  varmin = a;
  varmax = b;
  lowvar = a;
  highvar = b;
}

template <typename T>
int InputSettings::writeOutputToFile(T* obj) {
  if (!obj)
    return 1;

  TFile* file = TFile::Open(outputFileName.c_str(), "UPDATE");
  obj->Write(obj->GetName(), TObject::kOverwrite);
  file->Close();
  return 0;
}

// ---------------------------------------------------------------

struct FitLoader {
  private:
  public:
    InputSettings* inputs;
    TH3* data = nullptr;
    TH1* mass = nullptr;
    TF1* fit = nullptr;
    vector<TF1*> fitParts = {};

    vector<double> signal = {};
    vector<double> background = {};
    vector<double> signalPlusBackground = {};

    FitLoader() { inputs = new InputSettings(); }
    FitLoader(InputSettings& x) { inputs = &x; }

    TH3* loadData();
    TF1* loadFitFunction();
    vector<TF1*> loadFitParts();
    TH1* loadMass();

    string setHistNameFromTrain();
};

TF1* FitLoader::loadFitFunction() {
  TFile* file = TFile::Open(inputs->inputFileNameFit.c_str(), "READ");
  if (!file) {
    string s = "Could not open file " + inputs->inputFileNameFit;
    inputs->printLog(s, kErrors);
    return nullptr;
  }

  string fitName = inputs->getSaveNameFromPt("fit");
  TF1* f = (TF1*)file->Get(fitName.c_str());
  if (!f) {
    string s = "Could not find fit function " + fitName + " in file " + inputs->inputFileNameFit;
    inputs->printLog(s, kErrors);
    return nullptr;
  }

  fit = (TF1*)f->Clone();
  return f;
}

TH3* FitLoader::loadData() {
  TFile* f = TFile::Open(inputs->inputFileNameHist.c_str(), "READ");
  if (!f) {
    string s = "Could not open file " + inputs->inputFileNameHist;
    inputs->printLog(s, kErrors);
    return nullptr;
  }

  if (inputs->histName == "") {
    setHistNameFromTrain();
    string s = "FitLoader::loadData(): histName not set, setting to: " + inputs->histName;
    inputs->printLog(s, MyEnums::kInfo);
  }
  TH3* hist = (TH3*)f->Get(inputs->histName.c_str());
  if (!hist) {
    string s = "Could not find histogram " + inputs->histName + " in file " + inputs->inputFileNameHist;
    inputs->printLog(s, kErrors);
    return nullptr;
  }

  data = (TH3*)hist->Clone("data");
  return hist;
}

vector<TF1*> FitLoader::loadFitParts() {
  TFile* file = TFile::Open(inputs->inputFileNameFit.c_str(), "READ");
  if (!file) {
    string s = "Could not open file " + inputs->inputFileNameFit;
    inputs->printLog(s, kErrors);
    return {};
  }

  vector<TF1*> parts = {};
  int maxParts = 10; // To prevent infinite loop
  for (int i = 0; i < maxParts; i++) {
    string fName = inputs->getSaveNameFromPt("fit", "_f" + to_string(i));
    TF1* f = (TF1*)file->Get(fName.c_str());
    if (!f) break;

    parts.push_back(f);
  }

  fitParts = parts;
  return parts;
}

TH1* FitLoader::loadMass() {
  TH3* hist = (TH3*)data->Clone("hist");

  array<int, 2> ptBins = getProjectionBins(hist->GetXaxis(), inputs->ptmin, inputs->ptmax);
  array<int, 2> varBins = getProjectionBins(hist->GetYaxis(), inputs->varmin, inputs->varmax);

  inputs->lowpt = hist->GetXaxis()->GetBinLowEdge(ptBins[0]);
  inputs->highpt = hist->GetXaxis()->GetBinUpEdge(ptBins[1]);
  inputs->lowvar = hist->GetYaxis()->GetBinLowEdge(varBins[0]);
  inputs->highvar = hist->GetYaxis()->GetBinUpEdge(varBins[1]);

  string histName = inputs->getSaveNameFromPt("mass");
  histName = inputs->getSaveNameFromVar(histName);
  TH1* h = hist->ProjectionZ(histName.c_str(), ptBins[0], ptBins[1], varBins[0], varBins[1]);

  string s = TString::Format("FitLoader::loadMass(): \npt: %.1f-%.1f (%d, %d), \nvar: %.1f-%.1f (%d, %d), \nname: %s", inputs->lowpt, inputs->highpt, ptBins[0], ptBins[1], inputs->lowvar, inputs->highvar, varBins[0], varBins[1], histName.c_str()).Data();
  inputs->printLog(s, kDebug);
  return h;
}

string FitLoader::setHistNameFromTrain() {
  string s = "jet-fragmentation";
  if (inputs->train == 287744) s += "_id12406";
  if (inputs->train == 419996) s += "_id28293";
  if (inputs->train == 426828) s += "_id28293";

  s += "/data/V0/" + inputs->hadron + "Pt" + inputs->varName + "Mass";
  inputs->histName = s;
  return s;
}

// ---------------------------------------------------------------

struct SignalFinder {
  private:
  public:
    InputSettings* inputs;
    FitLoader* fl;

    double hSig, hBkg, hSigPlusBkg; // Calculated from hist
    double fSig, fBkg, fSigPlusBkg; // Calculated from fit

    vector<int> bkgFits = {};

    SignalFinder() { inputs = new InputSettings(); fl = new FitLoader(*inputs); }
    SignalFinder(InputSettings& x) { inputs = &x; fl = new FitLoader(x); }
    SignalFinder(FitLoader& x) { fl = &x; inputs = x.inputs; }

    // Select data within signal region
    void calcSigBkg();
};

void SignalFinder::calcSigBkg() {
  // Calculate signal and background with data hist and bkg fit function
  array<int, 2> sigBins = getProjectionBins(fl->data->GetXaxis(), inputs->signalRegionMin, inputs->signalRegionMax);
  double xmin = fl->data->GetXaxis()->GetBinLowEdge(sigBins[0]);
  double xmax = fl->data->GetXaxis()->GetBinUpEdge(sigBins[1]);

  hSigPlusBkg = fl->mass->Integral(sigBins[0], sigBins[1]);
  fSigPlusBkg = fl->fit->Integral(xmin, xmax);

  hBkg = 0;
  for (int i : bkgFits) {
    if (i <= 0 || i > bkgFits.size()) {
      string s = "SignalFinder::calcSigBkg Error: index " + to_string(i)
                 + " is out of range of fitParts (" + to_string(1)
                 + ", " + to_string(bkgFits.size()) + ")";
      inputs->printLog(s, kErrors);
      return;
    }
    TF1* f = (TF1*)fl->fitParts[i - 1]->Clone();
    f->SetRange(xmin, xmax);
    hBkg += f->Integral(xmin, xmax);
  }
  fBkg = hBkg; // Both calculated from fit

  hSig = hSigPlusBkg - hBkg;
  fSig = fSigPlusBkg - fBkg;
}

// ---------------------------------------------------------------

struct Plotter {
  private:
  public:
    InputSettings* inputs;
    FitLoader* fl;
    SignalFinder* sf;

    TCanvas* canvas = nullptr;
    TH1F* frame = nullptr;
    TLegend* legend = nullptr;
    vector<TLatex*> latex = {};

    Plotter() { inputs = new InputSettings(); fl = new FitLoader(*inputs); sf = new SignalFinder(*inputs); }
    Plotter(InputSettings& x) { inputs = &x; fl = new FitLoader(x); sf = new SignalFinder(x); }
    Plotter(FitLoader& x) { fl = &x; inputs = x.inputs; sf = new SignalFinder(*inputs); }
    Plotter(SignalFinder& x) { sf = &x; fl = x.fl; inputs = x.inputs; }

    // FIXME: What do we want to plot?
    // Fits are done to normalised histograms, but we need total integral for relative efficiency

    // Can read out peak heights from TH3/1
    // Multiply fit parameters by peak height
    // Peak height factors out of integrals
    // Plot spectrum with cuts as bar graph?
    vector<double> setVarLimits(int n);
    void plotCutVariationSpectrum();
};

vector<double> Plotter::setVarLimits(int n) {
  int nBins = fl->data->GetNbinsY();
  if (n <= 0 || n > nBins) {
    string s = "Plotter::setVarLimits() Error: n = " + to_string(n)
               + " is out of range of data histogram (" + to_string(1)
               + ", " + to_string(nBins) + ")";
    inputs->printLog(s, kErrors);
    return {};
  }

  if (nBins % n != 0) {
    string s = "Plotter::setVarLimits() Error: nBins = " + to_string(nBins)
               + " is not divisible by n = " + to_string(n);
    inputs->printLog(s, kErrors);
    return {};
  }

  // n hists to plot, each with a different cut
  int chunk = nBins / n;
  double xmin = fl->data->GetYaxis()->GetBinLowEdge(1);
  double xmax = fl->data->GetYaxis()->GetBinUpEdge(nBins);

  vector<double> limits;
  for (int i = 0; i < n; i++) {
    if (inputs->cutType == kUpper) {
      int b = nBins - i * chunk;
      limits.push_back(fl->data->GetYaxis()->GetBinUpEdge(b));
    }
    if (inputs->cutType == kLower) {
      int b = i * chunk + 1;
      limits.push_back(fl->data->GetYaxis()->GetBinLowEdge(b));
    }
  }
  inputs->varLimits = limits;
  return inputs->varLimits;
}

void Plotter::plotCutVariationSpectrum() {
  vector<TH1*> hists;

  double minVal = fl->data->GetYaxis()->GetBinLowEdge(1);
  double maxVal = fl->data->GetYaxis()->GetBinUpEdge(fl->data->GetNbinsY());

  if (!canvas) canvas = new TCanvas("canvas", "canvas", 800, 600);
  if (!legend) legend = CreateLegend(0.55, 0.85, 0.55, 0.85, "", 0.04);
  for (int i = 0; i < inputs->varLimits.size(); i++) {
    if (inputs->cutType == kUpper)
      inputs->setVarLimits(minVal, inputs->varLimits[i]);
    else
      inputs->setVarLimits(inputs->varLimits[i], maxVal);

    TH1* h = fl->loadMass();
    setStyle(h, i, 1);
    legend->AddEntry(h, TString::Format("%.1f < %s < %.1f", inputs->lowvar, inputs->varSymbol.c_str(), inputs->highvar), "f");
    hists.push_back(h);
  }

  if (!frame) {
    double xMinFrame = hists[0]->GetXaxis()->GetXmin();
    double xMaxFrame = hists[0]->GetXaxis()->GetXmax();
    double yMinFrame = 0., yMaxFrame = hists[0]->GetMaximum() * 1.2;
    string xTitle = TString::Format("#it{M} (%s) [GeV/#it{c}^{2}]", formatHadronDaughters(inputs->hadron).c_str()).Data();
    string yTitle = "Counts";
    frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  }

  frame->Draw();
  for (auto h : hists) h->Draw("same bar");

  legend->Draw("same");
  canvas->SaveAs(inputs->outputFileName.c_str());
}

// ---------------------------------------------------------------

// Plots the spectrum for defined variations of the reconstruction parameters
void plotCutVariationSpectrum(InputSettings x) {
  string outputFileName = x.varName;
  if (x.cutType == kLower) outputFileName += "_Lower";
  if (x.cutType == kUpper) outputFileName += "_Upper";
  x.outputFileName = x.getSaveNameFromPt(outputFileName, ".pdf");

  FitLoader m(x);
  m.setHistNameFromTrain();
  m.loadData();
  m.loadFitFunction();
  m.loadFitParts();

  Plotter p(m);
  if (x.varLimits.empty())
    p.setVarLimits(x.nHists); // Cut off an equal part of the hist each time

  p.legend = CreateLegend(0.55, 0.85, 0.55, 0.85, "", 0.04);
  p.plotCutVariationSpectrum();
}

// Plot spectra for 1 pt bin
void plotCutVariationSpectrumOneBin() {
  InputSettings x; x.verbosity = kWarnings;
  x.hadron = "K0S";
  x.train = 426828;
  x.setFitType("GausGaus");
  x.cutType = kLower; // Set cut type for variation
  x.setVar(kCtau);
  x.setPt(1., 2.);
  x.nHists = 5;

  x.setInputFileNameFromTrain();
  x.inputFileNameFit = x.hadron + "_pol1" + x.fitName + "/" + x.hadron + "_pol1" + x.fitName + ".root";

  plotCutVariationSpectrum(x);
}

// Plot spectra for all pt bins
void plotCutVariationSpectrumAllBins() {
  gROOT->SetBatch();
  InputSettings x; x.verbosity = kWarnings;
  x.hadron = "K0S";
  x.train = 426828;
  x.setFitType("GausGaus");
  x.cutType = kLower; // Set cut type for variation
  // Radius, Ctau, CosPA, DCAp, DCAn, DCAd
  x.setVar(kDCAd);
  x.nHists = 5;

  x.setInputFileNameFromTrain();
  x.inputFileNameFit = x.hadron + "_pol1" + x.fitName + "/" + x.hadron + "_pol1" + x.fitName + ".root";
  x.setPtBinEdges();

  for (int i = 0; i < x.ptBinEdges.size(); i++) {
    x.setPt(x.ptBinEdges[i][0], x.ptBinEdges[i][1]);
    plotCutVariationSpectrum(x);
  }
}

#endif // CUTVARIATION_C
