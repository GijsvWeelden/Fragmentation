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

#ifndef __PLOTPERPCONE_H__
#define __PLOTPERPCONE_H__

// -------------------------------------------------------------------------------------------------
//
// Struct for input settings
//
// -------------------------------------------------------------------------------------------------

// Divide two histograms with protection against nan/null bin content
// Assumes uncorrelated errors
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

    double numRelError = numError / numerator;
    double denRelError = denError / denominator;

    double newBinContent = numerator / denominator;
    double newBinError = newBinContent * std::sqrt((numRelError * numRelError) + (denRelError * denRelError));

    result->SetBinContent(i, newBinContent);
    result->SetBinError(i, newBinError);
  }
  return result;
}

namespace MyStrings {
  string getPtString(string subscript);
  string getZString(string subscript);
  string getRatioString(string num, string den);
  string getOneOverString(string s);
  string getdYdXString(string y, string x);
  string getdYdPtString(string y);
  string getdYdZString(string y);
  string getVarRangeString(string var, double high);
  string getVarRangeString(double low, string var, double high);
  string getPtJetRangeString(double ptmin, double ptmax, bool addUnits);
  string getPtV0RangeString(double ptmin, double ptmax, bool addUnits);

  const string sALICE      = "ALICE";
  const string sAntikt     = "Anti-#it{k}_{T}";
  const string sCharged    = "ch";
  const string sCounts     = "Counts";
  const string sEta        = "#eta";
  const string sGevC       = "GeV/#it{c}";
  const string sGevCC      = "GeV/#it{c}^{2}";
  const string sJet        = "jet";
  const string sJets       = "jets";
  const string sMass       = "#it{M}";
  const string sNumber     = "#it{N}";
  const string sRadius     = "#it{R} = 0.4";
  const string sRatio      = "Ratio";
  const string sSigma      = "#sigma";
  const string sSqrtS      = "#sqrt{s} = 13.6 TeV";
  const string sPpData     = "pp data";
  const string sPythia     = "PYTHIA";
  const string sThisThesis = "This Thesis";

  const string sV0         = "V0";
  const string sK0S        = formatHadronName("K0S");
  const string sLambda     = formatHadronName("Lambda");
  const string sAntiLambda = formatHadronName("AntiLambda");

  // Strings derived from the ones above
  const string sChJets      = sCharged + " " + sJets;
  const string sChV0Jets    = sCharged + "+" + sV0 + " " + sJets;
  const string sAliceData   = sALICE + " " + sPpData;
  const string sAlicePythia = sALICE + " " + sPythia;

  const string sNjets       = sNumber + "_{" + sJets + "}";
  const string sNevts       = sNumber + "_{evts}";
  const string sNV0         = sNumber + "_{" + sV0 + "}";
  const string sNK0S        = sNumber + "_{" + sK0S + "}";

  const string sEtaJet      = sEta + "_{" + sJet + "}";
  const string sEtaV0       = sEta + "_{" + sV0 + "}";
  const string sEtaK0S      = sEta + "_{" + sK0S + "}";

  const string sEtaJetRange = "|" + sEtaJet + "| < 0.5";
  const string sEtaV0Range  = "|" + sEtaV0 + "| < 0.9";
  const string sEtaK0SRange = "|" + sEtaK0S + "| < 0.9";

  const string sPtJet       = getPtString(sJet);
  const string sPtV0        = getPtString(sV0);
  const string sPtK0S       = getPtString(sK0S);
  const string sZV0         = getZString(sV0);
  const string sZK0S        = getZString(sK0S);

  const string sJetsPerEvent     = getOneOverString(sNevts) + " " + getdYdXString(sNjets, sPtJet);
  const string sV0PtPerEvt       = getOneOverString(sNevts) + " " + getdYdXString(sNV0, sPtV0);
  const string sK0SPtPerEvt      = getOneOverString(sNevts) + " " + getdYdXString(sNK0S, sPtK0S);
  const string sV0ZPerEvt       = getOneOverString(sNevts) + " " + getdYdXString(sNV0, sZV0);
  const string sK0SZPerEvt      = getOneOverString(sNevts) + " " + getdYdXString(sNK0S, sZK0S);

  const string sV0PtPerJet       = getOneOverString(sNjets) + " " + getdYdXString(sNV0, sPtV0);
  const string sK0SPtPerJet      = getOneOverString(sNjets) + " " + getdYdXString(sNK0S, sPtK0S);
  const string sV0ZPerJet        = getOneOverString(sNjets) + " " + getdYdZString(sV0);
  const string sK0SZPerJet       = getOneOverString(sNjets) + " " + getdYdXString(sNK0S, sZK0S);
  const string sLambdaPerJet     = getOneOverString(sNjets) + " " + getdYdZString(sLambda);
  const string sAntiLambdaPerJet = getOneOverString(sNjets) + " " + getdYdZString(sAntiLambda);
}

string MyStrings::getPtString(string subscript) {
  if (subscript.empty())
    return TString::Format("#it{p}_{T}").Data();
  else
    return TString::Format("#it{p}_{T, %s}", subscript.c_str()).Data();
}
string MyStrings::getZString(string subscript) {
  if (subscript.empty())
    return TString::Format("#it{z}").Data();
  else
    return TString::Format("#it{z}_{%s}", subscript.c_str()).Data();
}
string MyStrings::getRatioString(string num, string den) {
  return TString::Format("#frac{%s}{%s}", num.c_str(), den.c_str()).Data();
}
string MyStrings::getOneOverString(string s) {
  return (getRatioString("1", s));
}
string MyStrings::getdYdXString(string y, string x) {
  return getRatioString("d" + y, "d" + x);
}
string MyStrings::getdYdPtString(string y) {
  return getdYdXString(y, getPtString(y));
}
string MyStrings::getdYdZString(string y) {
  return getdYdXString(y, getZString(y));
}
string MyStrings::getVarRangeString(double low, string var, double high) {
  string sLow  = TString::Format("%.0f", low).Data();
  string sHigh = TString::Format("%.0f", high).Data();
  return TString::Format("%s < %s < %s", sLow.c_str(), var.c_str(), sHigh.c_str()).Data();
}
string MyStrings::getVarRangeString(string var, double high) {
  string sHigh = TString::Format("%.0f", high).Data();
  return TString::Format("%s < %s", var.c_str(), sHigh.c_str()).Data();
}
string MyStrings::getPtJetRangeString(double ptmin, double ptmax, bool addUnits = true) {
  string s = TString::Format("%.f < %s < %.f", ptmin, sPtJet.c_str(), ptmax).Data();
  if (addUnits)
    s += TString::Format(" %s", sGevC.c_str()).Data();

  return s;
}
string MyStrings::getPtV0RangeString(double ptmin, double ptmax, bool addUnits = true) {
  string s = TString::Format("%.1f < %s < %.1f", ptmin, sPtV0.c_str(), ptmax).Data();
  if (addUnits)
    s += TString::Format(" %s", sGevC.c_str()).Data();

  return s;
}

namespace VerbosityLevels {
  enum Verbosity {kErrors, kWarnings, kInfo, kDebug, kDebugMax};
}
namespace HistogramTypes {
  enum HistType {kGetJets, kGetJetV0s, kGetPCV0s};
}

using namespace MyStrings;
using namespace VerbosityLevels;
using namespace HistogramTypes;

struct InputSettings{
  private:
  public:
    const int conesPerJet = 2;
    int train = 0;
    int rebinNumber = -1;
    string hadron = "";
    string histName = "";
    string inputFileName = "";
    string outputFileName = "";
    double ptmin = -1e3, ptmax = -1e3, lowpt = -1e3, highpt = -1e3;
    double ptjetmin = -1e3, ptjetmax = -1e3, lowptjet = -1e3, highptjet = -1e3;
    bool logplot = false;
    bool ratioplot = false;
    vector<vector<double>> ptBinEdges = {};

    VerbosityLevels::Verbosity verbosity = VerbosityLevels::kWarnings;

    double massWindowMin = -1., massWindowMax = -1.;
    double polInitx0 = -1., polInitx1 = -1., polInitx2 = -1.;
    double signalRegionMin = -1., signalRegionMax = -1.;
    double nSigma = -1., nSigmaL = -1., nSigmaR = -1.;

    string getHistName(HistogramTypes::HistType x);
    double getMass();
    string getNameFromJetPt(string prefix, string suffix);
    string getNameFromPt(string prefix, string suffix);
    bool passVerbosityCheck(VerbosityLevels::Verbosity verbThreshold);
    string printLog(string message, VerbosityLevels::Verbosity verbThreshold);
    int setHadron(string h);
    string setInputFileNameFromTrain();
    void setLowHighFromAxis(const TAxis* axis, double& low, double& high);
    void setJetPt(double a, double b);
    void setPt(double a, double b);
    vector<vector<double>> setPtBinEdgesFromHadron();
    vector<vector<double>> setPtBinEdgesSorted(vector<vector<double>> x);
    template <typename T> int writeOutputToFile(T* obj);
};

string InputSettings::getHistName(HistogramTypes::HistType x) {
  string s = "jet-fragmentation";
  if (train == 436232) {
    s += (x == HistogramTypes::kGetPCV0s) ? "_id24581" : "_id24580";
  }
  s += "/data/";
  switch (x) {
    case HistogramTypes::kGetJets:
      s += "jets/jetPtEtaPhi";
      break;
    case HistogramTypes::kGetJetV0s:
      s += "jets/V0/jetPtK0SPtMass";
      break;
    case HistogramTypes::kGetPCV0s:
      s += "PC/JetPtK0SPtMass";
      break;
    default:
      printLog("InputSettings::getHistName() Error: invalid x", VerbosityLevels::kErrors);
      s = "";
  }
  return s;
}

double InputSettings::getMass() {
  if (this->hadron == "K0S")
    return MassK0S;
  if (this->hadron == "Lambda" || this->hadron == "AntiLambda")
    return MassLambda;

  return -1.;
}

string InputSettings::getNameFromJetPt(string prefix, string suffix = "") {
  string s = TString::Format("%s_jetpt%.f-%.f%s", prefix.c_str(), lowptjet, highptjet, suffix.c_str()).Data();
  return s;
}

string InputSettings::getNameFromPt(string prefix, string suffix = "") {
  string s = TString::Format("%s_pt%.1f-%.1f%s", prefix.c_str(), lowpt, highpt, suffix.c_str()).Data();
  return s;
}

bool InputSettings::passVerbosityCheck(VerbosityLevels::Verbosity verbThreshold) {
  return (verbosity >= verbThreshold);
}

string InputSettings::printLog(string message, VerbosityLevels::Verbosity verbThreshold) {
  if (!passVerbosityCheck(verbThreshold))
    return "";

  cout << message << endl;
  return message;
}

int InputSettings::setHadron(string h) {
  if (h == "K0S" || h == "Lambda" || h == "AntiLambda") {
    hadron = h;
    return 0;
  } else {
    printLog(TString::Format("InputSettings::setHadron() Error: requested invalid hadron %s", h.c_str()).Data(), kErrors);
    return 1;
  }
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(this->train) + "/AnalysisResults.root";
  this->inputFileName = s;
  return s;
}

void InputSettings::setPt(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setPt() Error: ptmin > ptmax", kErrors);
    return;
  }
  this->ptmin = a;
  this->ptmax = b;
  this->lowpt = a;
  this->highpt = b;
}

void InputSettings::setJetPt(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setJetPt() Error: ptjetmin > ptjetmax", kErrors);
    return;
  }
  this->ptjetmin = a;
  this->ptjetmax = b;
  this->lowptjet = a;
  this->highptjet = b;
}

vector<vector<double>> InputSettings::setPtBinEdgesFromHadron() {
  vector<vector<double>> x;

  x.push_back({0., 1.});
  x.push_back({1., 2.});
  x.push_back({2., 3.});
  x.push_back({3., 4.});
  x.push_back({4., 5.});
  x.push_back({5., 10.});
  x.push_back({10., 15.});
  x.push_back({15., 20.});
  if (this->hadron == "K0S") {
    x.push_back({20., 25.});
    x.push_back({25., 30.});
    x.push_back({30., 40.});
  } else {
    x.push_back({20., 30.});
    x.push_back({30., 40.});
  }
  x = this->setPtBinEdgesSorted(x);
  return x;
}

vector<vector<double>> InputSettings::setPtBinEdgesSorted(vector<vector<double>> x) {
  sort(x.begin(), x.end());
  this->ptBinEdges = x;
  return x;
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

// -------------------------------------------------------------------------------------------------
//
// Struct for plotting
//
// -------------------------------------------------------------------------------------------------

struct Plotter {
  private:
    string drawOption; // Private because it requires caution with spaces

  public:
    InputSettings* inputs;

    TCanvas* canvas = nullptr;
    TH1F* frame = nullptr;
    TLegend* legend = nullptr;
    vector<TH1*> hists = {};
    vector<TObject*> objects = {};

    double textSize = 0.04;
    const string sGevC = "GeV/#it{c}";
    const string sGevCC = "GeV/#it{c}^{2}";
    const string sNjets = "#it{N}_{jets}";
    const string sPtJet = "#it{p}_{T, jet} (GeV/#it{c})";
    const string sRatio = "Ratio";
    const string sZV0 = "#it{z}_{V0}";

    Plotter() { inputs = new InputSettings(); }
    Plotter(InputSettings& x) { inputs = &x; }

    void addLatex(double x, double y, string s);
    void addLine(double x0, double y0, double x1, double y1, int styleNumber, int lineStyle, int lineWidth);
    void makeCanvas(string s, double x, double y);
    void makeFrame(string sx, string sy);
    void makeFrame(double x0, double x1, double y0, double y1, string sx, string sy);
    void makeLegend(double x0, double x1, double y0, double y1, string s);
    void plot();
    void reset();
    void setHistStyles();

    string getMassString();
    string getdNdXString(string var);
    string setDrawOption(string s);
};

void Plotter::addLatex(double x, double y, string s) {
  TLatex* l = CreateLatex(x, y, s.c_str(), textSize);
  objects.push_back(l);
}
void Plotter::addLine(double x0, double y0, double x1, double y1, int styleNumber, int lineStyle = 9, int lineWidth = 3) {
  TLine* l = new TLine(x0, y0, x1, y1);
  setStyle(l, styleNumber, lineStyle, lineWidth);
  objects.push_back(l);
}
void Plotter::makeCanvas(string s = "canvas", double x = 800, double y = 600) {
  canvas = new TCanvas(s.c_str(), s.c_str(), x, y);
  canvas->SetLogy(inputs->logplot);
}
void Plotter::makeFrame(string sx, string sy) {
  if (hists.empty()) {
    inputs->printLog("Plotter::makeFrame() hists vector is empty! Aborting", kErrors);
    return;
  }
  if (!canvas) makeCanvas();

  double xMinFrame = hists[0]->GetXaxis()->GetXmin();
  double xMaxFrame = hists[0]->GetXaxis()->GetXmax();
  double yMinFrame = getLowerBound(hists, 0, 0) * 0.9;
  double yMaxFrame = getUpperBound(hists, 0, 0) * 1.2;

  inputs->printLog(TString::Format("Plotter::makeFrame() x: %.2f, %.2f \ny: %.2f, %.2f", xMinFrame, xMaxFrame, yMinFrame, yMaxFrame).Data(), VerbosityLevels::kDebug);

  if (inputs->logplot) {
    roundToNextPowerOfTen(yMaxFrame);
    if (yMinFrame < 1e-12) {
      int xBin = hists[0]->FindLastBinAbove(0.);
      yMinFrame = hists[0]->GetBinContent(xBin) * 0.5;
      yMaxFrame *= 1.75;
      roundToPrevPowerOfTen(yMinFrame);
    }
  }
  inputs->printLog(TString::Format("Plotter::makeFrame() x: %.2f, %.2f \ny: %.2f, %.2f", xMinFrame, xMaxFrame, yMinFrame, yMaxFrame).Data(), VerbosityLevels::kDebug);
  makeFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, sx, sy);
}
void Plotter::makeFrame(double x0, double x1, double y0, double y1, string sx, string sy) {
  if (!canvas) makeCanvas();
  frame = DrawFrame(x0, x1, y0, y1, sx, sy);
}
void Plotter::makeLegend(double x0, double x1, double y0, double y1, string s) {
  legend = CreateLegend(x0, x1, y0, y1, s.c_str(), textSize);
}

void Plotter::plot() {
  if (hists.empty())
    inputs->printLog("Plotter::plot(): Hist vector is empty!", kWarnings);

  if (inputs->ratioplot) {
    TH1* baseCopy = (TH1*)hists[0]->Clone("baseCopy");
    // for (auto& h : hists) h->Divide(baseCopy);
    for (auto& h : hists) h = divideWithProtection(h, baseCopy);
  }

  if (!canvas) makeCanvas();
  if (!frame) makeFrame("x", "y");

  frame->Draw();
  if (legend) legend->Draw("same");
  for (auto o : objects)
    o->Draw("same");

  for (auto h : hists) {
    h->Draw(("same" + drawOption).c_str());
    if (inputs->passVerbosityCheck(Verbosity::kDebugMax))
      h->Print("all");
    else if (inputs->passVerbosityCheck(Verbosity::kDebug))
      h->Print();
  }
  canvas->SaveAs(inputs->outputFileName.c_str());
}

void Plotter::reset() {
  canvas = nullptr;
  frame = nullptr;
  legend = nullptr;
  objects.clear();
  hists.clear();
}

void Plotter::setHistStyles() {
  for (unsigned int i = 0; i < hists.size(); i++)
    setStyle(hists[i], i);
}

string Plotter::getMassString() {
  return TString::Format("#it{M}(%s) (%s)", formatHadronDaughters(inputs->hadron).c_str(), sGevCC.c_str()).Data();
}
string Plotter::getdNdXString(string var) {
  return TString::Format("#frac{d#it{N}_{V0}}{d%s}", var.c_str()).Data();
}
string Plotter::setDrawOption(string s) {
  if (s.at(0) == ' ') { // Single quotes, because checking for a char
    drawOption = s;
  } else {
    drawOption = " " + s;
  }
  return drawOption;
}

// -------------------------------------------------------------------------------------------------
//
// Interface
//
// -------------------------------------------------------------------------------------------------

struct PerpCone {
  InputSettings* inputs;
  Plotter* plotter;

  PerpCone() { inputs = new InputSettings(); plotter = new Plotter(*inputs); }
  PerpCone(InputSettings& i) { inputs = &i; plotter = new Plotter(*inputs); }

  TH1* getMassHist(HistogramTypes::HistType type);
  double getNjets();
  array<TH1*, 2> getPerpConePtHists(bool rebin);
  TH1* getPtHist(HistogramTypes::HistType type);

  void plotPerpConeMass();
  void plotPerpConePt();
  TH1* rebinPtHist(TH1* hist);

  bool greaterThan(double a, double b, double epsilon = 1e-5) { return a - b > epsilon; }
  bool lessThan(double a, double b, double epsilon = 1e-5) { return greaterThan(b, a, epsilon); }
};

TH1* PerpCone::getMassHist(HistogramTypes::HistType type) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    inputs->printLog("Could not open file " + inputs->inputFileName, kErrors);
    return nullptr;
  }

  TH3* h3 = (TH3*)file->Get(inputs->getHistName(type).c_str());

  array<int, 2> jetBins = getProjectionBins(h3->GetXaxis(), inputs->ptjetmin, inputs->ptjetmax);
  array<int, 2> v0Bins = getProjectionBins(h3->GetYaxis(), inputs->ptmin, inputs->ptmax);
  inputs->lowptjet = h3->GetXaxis()->GetBinLowEdge(jetBins[0]);
  inputs->highptjet = h3->GetXaxis()->GetBinUpEdge(jetBins[1]);
  inputs->lowpt = h3->GetYaxis()->GetBinLowEdge(v0Bins[0]);
  inputs->highpt = h3->GetYaxis()->GetBinUpEdge(v0Bins[1]);

  string name = "m" + inputs->hadron;
  if (type == HistogramTypes::kGetPCV0s)
    name += "inPCs";
  if (type == HistogramTypes::kGetJetV0s)
    name += "inJets";
  name = inputs->getNameFromJetPt(name);
  name = inputs->getNameFromPt(name);

  TH1* hist = h3->ProjectionZ(name.c_str(), jetBins[0], jetBins[1], v0Bins[0], v0Bins[1]);
  return hist;
}

TH1* PerpCone::getPtHist(HistogramTypes::HistType type) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    inputs->printLog("Could not open file " + inputs->inputFileName, kErrors);
    return nullptr;
  }

  TH3* h3 = (TH3*)file->Get(inputs->getHistName(type).c_str());

  array<int, 2> jetBins = getProjectionBins(h3->GetXaxis(), inputs->ptjetmin, inputs->ptjetmax);
  inputs->lowptjet = h3->GetXaxis()->GetBinLowEdge(jetBins[0]);
  inputs->highptjet = h3->GetXaxis()->GetBinUpEdge(jetBins[1]);

  string name = "pt" + inputs->hadron;
  if (type == HistogramTypes::kGetPCV0s)
    name += "inPCs";
  if (type == HistogramTypes::kGetJetV0s)
    name += "inJets";
  name = inputs->getNameFromJetPt(name);

  TH1* hist = h3->ProjectionY(name.c_str(), jetBins[0], jetBins[1], 1, h3->GetNbinsZ());
  return hist;
}

double PerpCone::getNjets() {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    inputs->printLog("Could not open file " + inputs->inputFileName, kErrors);
    return -1.;
  }

  TH3* hJetPtEtaPhi = (TH3*)file->Get(inputs->getHistName(kGetJets).c_str());
  TH1* hJetPt = hJetPtEtaPhi->ProjectionX("jetPt");
  array<int, 2> jetBins = getProjectionBins(hJetPt->GetXaxis(), inputs->ptjetmin, inputs->ptjetmax);
  double nJets = hJetPt->Integral(jetBins[0], jetBins[1]);
  return nJets;
}

void PerpCone::plotPerpConeMass() {
  plotter->hists.clear();
  inputs->outputFileName = inputs->getNameFromPt(inputs->getNameFromJetPt("pcMass"), ".pdf");

  plotter->makeLegend(0.65, 0.75, 0.7, 0.8, "");
  array<int, 2> jetBins, v0Bins;

  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");

  double nJets = getNjets();
  if (nJets <= 0.) {
    inputs->printLog(TString::Format("PerpCone::plotPerpConeMass() Error: could not get number of jets for pt (%.1f, %.1f)", inputs->ptjetmin, inputs->ptjetmax).Data(), kErrors);
    return;
  }

  TH1* hJet = getMassHist(kGetJetV0s);
  setStyle(hJet, 0);
  hJet->Scale(1. / nJets, "width");
  plotter->hists.push_back(hJet);
  plotter->legend->AddEntry(hJet, "V0s in jet cone");

  TH1* hPC = getMassHist(kGetPCV0s);
  setStyle(hPC, 1);
  hPC->Scale(1. / (inputs->conesPerJet * nJets), "width");
  plotter->hists.push_back(hPC);
  plotter->legend->AddEntry(hPC, "V0s in UE");

  if (inputs->ratioplot) {
    plotter->makeFrame(0.4, 0.6, 0., 1., plotter->getMassString(), sRatio);
  } else {
    plotter->makeFrame(plotter->getMassString(), TString::Format("#frac{1}{%s} %s", sNjets.c_str(), plotter->getdNdXString("#it{M}").c_str()).Data());
  }

  plotter->addLatex(0.25, 0.80, sThisThesis + " " + sAliceData);
  plotter->addLatex(0.25, 0.75, sSqrtS);
  plotter->addLatex(0.25, 0.70, sAntikt + " " + sChV0Jets);
  plotter->addLatex(0.25, 0.65, getVarRangeString(inputs->lowptjet, sPtJet, inputs->highptjet) + " " + sGevC);
  plotter->addLatex(0.25, 0.60, getVarRangeString(inputs->lowpt, sPtK0S, inputs->highpt) + " " + sGevC);

  plotter->plot();
}

array<TH1*, 2> PerpCone::getPerpConePtHists(bool rebin = true) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  double nJets = getNjets();
  if (nJets <= 0.) {
    inputs->printLog("PerpCone::plotPerpConeMass() Error: could not get number of jets", kErrors);
    return {nullptr, nullptr};
  }

  TH1* hJet = getPtHist(kGetJetV0s);
  TH1* hPC = getPtHist(kGetPCV0s);
  if (rebin) {
    hJet = rebinPtHist(hJet);
    hPC = rebinPtHist(hPC);
  }

  hJet->Scale(1. / nJets, "width");
  hPC->Scale(1. / (inputs->conesPerJet * nJets), "width");

  return {hJet, hPC};
}

void PerpCone::plotPerpConePt() {
  plotter->hists.clear();
  inputs->outputFileName = inputs->getNameFromJetPt("pc");
  if (inputs->ratioplot)
    inputs->outputFileName += "ratio";
  inputs->outputFileName += ".pdf";

  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  double nJets = getNjets();
  if (nJets <= 0.) {
    inputs->printLog("PerpCone::plotPerpConeMass() Error: could not get number of jets", kErrors);
    return;
  }
  plotter->makeLegend(0.45, 0.6, 0.3, 0.4, "");

  const bool doRebin = true;
  array<TH1*, 2> ptHists = getPerpConePtHists(doRebin);
  // TH1* hJet = getPtHist(HistogramTypes::kGetJetV0s);
  // hJet = rebinPtHist(hJet);
  // hJet->Scale(1. / nJets, "width");
  TH1* hJet = ptHists[0];

  setStyle(hJet, 0);
  plotter->hists.push_back(hJet);
  plotter->legend->AddEntry(hJet, "V0s in jet cone");

  TH1* hPC = ptHists[1];
  setStyle(hPC, 1);
  plotter->hists.push_back(hPC);
  plotter->legend->AddEntry(hPC, "V0s in UE");

  if (inputs->ratioplot) {
    plotter->makeFrame(0., inputs->ptjetmax, 1e-5, 1., sPtK0S, sRatio);
  } else {
    string xTitle = sPtK0S;
    string yTitle = getOneOverString(sNumber + "_{jets, cones}") + getdYdXString(sNK0S, sPtK0S);
    plotter->makeFrame(0., inputs->ptjetmax, 1e-5, 0.2, xTitle, yTitle);
  }

  plotter->addLatex(0.47, 0.85, sThisThesis + ", " + sAliceData);
  plotter->addLatex(0.47, 0.80, sSqrtS);
  plotter->addLatex(0.47, 0.75, sAntikt + " " + sChV0Jets);
  plotter->addLatex(0.47, 0.70, sRadius + ", " + sEtaJetRange);
  plotter->addLatex(0.47, 0.65, getPtJetRangeString(inputs->lowptjet, inputs->highptjet, true));

  plotter->plot();
}

TH1* PerpCone::rebinPtHist(TH1* hist) {
  const int nBinsK0S = 11;
  const double edgesK0S[nBinsK0S + 1] = {0., 1., 2., 3., 4., 5., 10., 15., 20., 25., 30., 40.};

  const int nBinsLAL = 10;
  const double edgesLAL[nBinsLAL + 1] = {0., 1., 2., 3., 4., 5., 10., 15., 20., 30., 40.};

  string name = hist->GetName();
  name += "_rebinned";

  TH1D* h;
  if (inputs->hadron == "K0S") {
    h = new TH1D(name.c_str(), hist->GetTitle(), nBinsK0S, edgesK0S);
  } else {
    h = new TH1D(name.c_str(), hist->GetTitle(), nBinsLAL, edgesLAL);
  }

  // +2 to include underflow and overflow bins
  double newContents[h->GetNbinsX() + 2];
  double newErrorsSquared[h->GetNbinsX() + 2];
  for (int ib = 0; ib <= hist->GetNbinsX() + 1; ib++) {
    double binContent = hist->GetBinContent(ib);
    double binCentre = hist->GetBinCenter(ib);
    double binError = hist->GetBinError(ib);

    int newBin = h->FindBin(binCentre);
    newContents[newBin] += binContent;
    newErrorsSquared[newBin] += binError * binError;
  }

  // Loop over h to set the content and errors
  for (int ib = 0; ib <= h->GetNbinsX() + 1; ib++) {
    double content = newContents[ib];
    double errorSquared = newErrorsSquared[ib];

    h->SetBinContent(ib, content);
    h->SetBinError(ib, std::sqrt(errorSquared));
  }

  return h;
}

PerpCone setupPerpCone(double ptjetmin, double ptjetmax) {
  PerpCone p;
  p.inputs->verbosity = kDebug;
  p.inputs->hadron = "K0S";
  p.inputs->train = 436232;
  p.inputs->setInputFileNameFromTrain();
  p.inputs->setJetPt(ptjetmin, ptjetmax);
  p.inputs->setPtBinEdgesFromHadron();
  return p;
}

void perpConeMass(bool doRatio) {
  PerpCone p = setupPerpCone(20., 30.);
  p.inputs->ratioplot = doRatio;

  if (p.inputs->passVerbosityCheck(VerbosityLevels::kDebug)) {
    string s = TString::Format("Ptjet: (%.f, %.f)\n", p.inputs->ptjetmin, p.inputs->ptjetmax).Data();
    s += "Pt bin edges: ";
    for (int ipt = 0; ipt < p.inputs->ptBinEdges.size(); ipt++) {
      s += TString::Format("(%.1f, %.1f) ", p.inputs->ptBinEdges[ipt][0], p.inputs->ptBinEdges[ipt][1]).Data();
    }
    p.inputs->printLog(s, VerbosityLevels::kDebug);
  }

  for (int ipt = 0; ipt < p.inputs->ptBinEdges.size(); ipt++) {
    p.inputs->setPt(p.inputs->ptBinEdges[ipt][0], p.inputs->ptBinEdges[ipt][1]);
    p.inputs->printLog("Pt: " + to_string(p.inputs->ptmin) + " - " + to_string(p.inputs->ptmax), VerbosityLevels::kDebug);

    if (p.greaterThan(p.inputs->ptmin, p.inputs->ptjetmax)) {
      p.inputs->printLog("Stopping plotting because ptmin >= ptjetmax", VerbosityLevels::kDebug);
      break;
    }

    p.plotPerpConeMass();
  }
}

void perpconept(double ptjetmin, double ptjetmax, bool doRatio) {
  PerpCone p = setupPerpCone(ptjetmin, ptjetmax);
  p.inputs->verbosity = VerbosityLevels::kInfo;
  p.inputs->logplot = true;
  p.inputs->ratioplot = doRatio;
  p.inputs->outputFileName = p.inputs->getNameFromJetPt("pc");
  if (doRatio)
    p.inputs->outputFileName += "_ratio";
  p.inputs->outputFileName += ".pdf";

  array<TH1*, 2> ptHists = p.getPerpConePtHists(true);
  if (p.inputs->passVerbosityCheck(VerbosityLevels::kDebug)) {
    for (auto h : ptHists) {
      h->Print("all");
    }
  }
  p.plotter->makeLegend(0.45, 0.6, 0.3, 0.4, "");

  TH1* hJet = ptHists[0];
  setStyle(hJet, 0);
  p.plotter->legend->AddEntry(hJet, "V0s in jet cone");

  TH1* hPC = ptHists[1];
  setStyle(hPC, 1);
  p.plotter->legend->AddEntry(hPC, "V0s in UE");

  string xTitle = sPtK0S;
  string yTitle = getOneOverString(sNumber + "_{jets, cones}") + " " + getdYdXString(sNK0S, sPtK0S);
  double xMinFrame = 0., xMaxFrame = p.inputs->ptjetmax;
  double yMinFrame = 1e-5;
  double yMaxFrame = 0.2;
  if (p.lessThan(p.inputs->ptjetmin, 11.)) {
    xMaxFrame = 25.;
    yMinFrame = 1e-8;
  }
  p.plotter->makeFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  if (p.inputs->ratioplot)
    p.plotter->makeFrame(xMinFrame, xMaxFrame, yMinFrame, 1., sPtK0S, sRatio);

  p.plotter->addLatex(0.47, 0.85, sThisThesis + ", " + sAliceData);
  p.plotter->addLatex(0.47, 0.80, sSqrtS);
  p.plotter->addLatex(0.47, 0.75, sAntikt + " " + sChV0Jets);
  p.plotter->addLatex(0.47, 0.70, sRadius + ", " + sEtaJetRange);
  p.plotter->addLatex(0.47, 0.65, getPtJetRangeString(p.inputs->lowptjet, p.inputs->highptjet, true));

  if (doRatio) {
    TH1* baseCopy = (TH1*)hJet->Clone("baseCopy");
    hJet = divideWithProtection(hJet, baseCopy);
    hPC = divideWithProtection(hPC, baseCopy);
  }
  if (p.inputs->passVerbosityCheck(VerbosityLevels::kInfo)) {
    hJet->Print("all");
    hPC->Print("all");
  }

  p.plotter->hists.push_back(hJet);
  p.plotter->hists.push_back(hPC);
  p.plotter->setDrawOption("p");

  p.inputs->ratioplot = false;

  p.plotter->plot();
}

#endif
