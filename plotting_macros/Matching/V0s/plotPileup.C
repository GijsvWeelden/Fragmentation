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

#include "../../histUtils.C"
#include "../../plotUtils.C"

#ifndef __PLOTPERPCONE_H__
#define __PLOTPERPCONE_H__

// -------------------------------------------------------------------------------------------------
//
// Struct for input settings
//
// -------------------------------------------------------------------------------------------------

// Divide two histograms with protection against nan/null bin content
// Assumes uncorrelated errors
// TH1* divideWithProtection(TH1* base, TH1* divideBy, double threshold = 1e-25) {
//   // TODO: Check if same binning
//   TH1* result = (TH1*)base->Clone("result");
//   result->Reset();

//   for (int i = 0; i <= 1 + base->GetNbinsX(); i++) {
//     double numerator = base->GetBinContent(i);
//     double numError = base->GetBinError(i);
//     double denominator = divideBy->GetBinContent(i);
//     double denError = divideBy->GetBinError(i);

//     if (std::isnan(numerator) || std::isnan(denominator))
//       continue;
//     else if (std::abs(numerator) < threshold || std::abs(denominator) < threshold)
//       continue;

//     double newBinContent = numerator / denominator;

//     if (std::isnan(numError) || std::isnan(denError))
//       continue;

//     double numRelError = numError / numerator;
//     double denRelError = denError / denominator;
//     double newBinError = newBinContent * std::sqrt((numRelError * numRelError) + (denRelError * denRelError));

//     result->SetBinContent(i, newBinContent);
//     result->SetBinError(i, newBinError);
//   }
//   return result;
// }

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
  bool is_valid(int v) { return (v >= kErrors && v <= kDebugMax); }
  string to_string(Verbosity v) {
    switch (v) {
      case kErrors:   return "kErrors";
      case kWarnings: return "kWarnings";
      case kInfo:     return "kInfo";
      case kDebug:    return "kDebug";
      case kDebugMax: return "kDebugMax";
      default:        return "Unknown";
    }
  }
}
namespace HistogramTypes {
  enum HistType {kInclusive, kInclusiveWrongCollision, kJetPt, kJetPtWrongCollision, kJetZ, kJetZWrongCollision, kMatchedJetPt, kMatchedJetPtWrongCollision, kMatchedJetZ, kMatchedJetZWrongCollision};
  bool is_valid(int x) { return (x >= kInclusive && x <= kMatchedJetZWrongCollision); }
  string to_string(HistType x) {
    switch (x) {
      case kInclusive:                  return "kInclusive";
      case kInclusiveWrongCollision:    return "kInclusiveWrongCollision";
      case kJetPt:                      return "kJetPt";
      case kJetPtWrongCollision:        return "kJetPtWrongCollision";
      case kJetZ:                       return "kJetZ";
      case kJetZWrongCollision:         return "kJetZWrongCollision";
      case kMatchedJetPt:               return "kMatchedJetPt";
      case kMatchedJetPtWrongCollision: return "kMatchedJetPtWrongCollision";
      case kMatchedJetZ:                return "kMatchedJetZ";
      case kMatchedJetZWrongCollision:  return "kMatchedJetZWrongCollision";
      default: return "Unknown";
    }
  }
}
namespace ProjectionTypes {
  enum ProjType {kPtV0, kZV0, kMass};
  bool is_valid(int x) { return (x >= kPtV0 && x <= kMass); }
  string to_string(ProjType x) {
    switch (x) {
      case kPtV0: return "kPtV0";
      case kZV0:  return "kZV0";
      case kMass: return "kMass";
      default:    return "Unknown";
    }
  }
}

using namespace MyStrings;
// using namespace VerbosityLevels;
// using namespace HistogramTypes;

struct InputSettings{
  private:
  public:
    int train = 0;
    int rebinNumber = -1;
    string hadron = "";
    string histName = "";
    string inputFileName = "";
    string outputFileName = "";
    double ptmin = -1e3, ptmax = -1e3, lowpt = -1e3, highpt = -1e3;
    double ptjetmin = -1e3, ptjetmax = -1e3, lowptjet = -1e3, highptjet = -1e3;
    double etamin = -0.9, etamax = 0.9;
    bool logplot = false;
    bool ratioplot = false;
    vector<vector<double>> ptBinEdges = {};

    VerbosityLevels::Verbosity verbosity = VerbosityLevels::kWarnings;

    double massWindowMin = -1., massWindowMax = -1.;
    double polInitx0 = -1., polInitx1 = -1., polInitx2 = -1.;
    double signalRegionMin = -1., signalRegionMax = -1.;
    double nSigma = -1., nSigmaL = -1., nSigmaR = -1.;

    double getMass();
    string getNameFromJetPt(string prefix, string suffix);
    string getNameFromPt(string prefix, string suffix);
    bool passVerbosityCheck(VerbosityLevels::Verbosity verbThreshold);
    string printLog(string message, VerbosityLevels::Verbosity verbThreshold);
    int setHadron(string h);
    string setInputFileNameFromTrain();
    void setLowHighFromAxis(const TAxis* axis, double& low, double& high);
    void setEta(double a, double b);
    void setJetPt(double a, double b);
    void setPt(double a, double b);
    vector<vector<double>> setPtBinEdgesFromHadron();
    vector<vector<double>> setPtBinEdgesSorted(vector<vector<double>> x);
    template <typename T> int writeOutputToFile(T* obj);
};

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
    printLog(TString::Format("InputSettings::setHadron() Error: requested invalid hadron %s", h.c_str()).Data(), VerbosityLevels::kErrors);
    return 1;
  }
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(this->train) + "/AnalysisResults.root";
  this->inputFileName = s;
  return s;
}

void InputSettings::setEta(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setEta() Error: etamin > etamax", VerbosityLevels::kErrors);
    return;
  }
  this->etamin = a;
  this->etamax = b;
}

void InputSettings::setPt(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setPt() Error: ptmin > ptmax", VerbosityLevels::kErrors);
    return;
  }
  this->ptmin = a;
  this->ptmax = b;
  this->lowpt = a;
  this->highpt = b;
}

void InputSettings::setJetPt(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setJetPt() Error: ptjetmin > ptjetmax", VerbosityLevels::kErrors);
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
    string setDrawOption(string s);
};

void Plotter::addLatex(double x, double y, string s) {
  TLatex* l = plotutils::CreateLatex(x, y, s.c_str(), textSize);
  objects.push_back(l);
}
void Plotter::addLine(double x0, double y0, double x1, double y1, int styleNumber, int lineStyle = 9, int lineWidth = 3) {
  TLine* l = new TLine(x0, y0, x1, y1);
  plotutils::setStyle(l, styleNumber, lineStyle, lineWidth);
  objects.push_back(l);
}
void Plotter::makeCanvas(string s = "canvas", double x = 800, double y = 600) {
  canvas = new TCanvas(s.c_str(), s.c_str(), x, y);
  canvas->SetLogy(inputs->logplot);
}
void Plotter::makeFrame(string sx, string sy) {
  if (hists.empty()) {
    inputs->printLog("Plotter::makeFrame() hists vector is empty! Aborting", VerbosityLevels::kErrors);
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
  frame = plotutils::DrawFrame(x0, x1, y0, y1, sx, sy);
}
void Plotter::makeLegend(double x0, double x1, double y0, double y1, string s) {
  legend = plotutils::CreateLegend(x0, x1, y0, y1, s.c_str(), textSize);
}

void Plotter::plot() {
  if (hists.empty())
    inputs->printLog("Plotter::plot(): Hist vector is empty!", VerbosityLevels::kWarnings);

  if (inputs->ratioplot) {
    TH1* baseCopy = (TH1*)hists[0]->Clone("baseCopy");
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
    if (inputs->passVerbosityCheck(VerbosityLevels::kDebugMax))
      h->Print("all");
    else if (inputs->passVerbosityCheck(VerbosityLevels::kDebug))
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
    plotutils::setStyle(hists[i], i);
}

string Plotter::getMassString() {
  return TString::Format("#it{M}(%s) (%s)", formatHadronDaughters(inputs->hadron).c_str(), sGevCC.c_str()).Data();
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

struct PileUp {
  InputSettings* inputs;
  Plotter* plotter;

  PileUp() { inputs = new InputSettings(); plotter = new Plotter(*inputs); }
  PileUp(InputSettings& i) { inputs = &i; plotter = new Plotter(*inputs); }

  void autoFrame(HistogramTypes::HistType type, ProjectionTypes::ProjType proj);
  void autoLatex(HistogramTypes::HistType type, ProjectionTypes::ProjType proj);

  string getHistName(HistogramTypes::HistType x);
  template <typename T> T getHistFromFile(HistogramTypes::HistType type);
  TH1* getHistProjection(HistogramTypes::HistType type, ProjectionTypes::ProjType proj, string name);

  double getNevts();
  double getNjets(double ptjetmin, double ptjetmax);
  bool greaterThan(double a, double b, double epsilon = 1e-5) { return a - b > epsilon; }
  bool lessThan(double a, double b, double epsilon = 1e-5) { return greaterThan(b, a, epsilon); }

  void plotInclusivePt(bool doRatio);
  void plotInJet(bool doRatio, bool doZ);
};

void PileUp::autoFrame(HistogramTypes::HistType type, ProjectionTypes::ProjType proj) {
  plotter->frame = nullptr;
  switch (type) {
    case HistogramTypes::kInclusive:
    case HistogramTypes::kInclusiveWrongCollision:
      if (proj == ProjectionTypes::kPtV0) {
        if (inputs->ratioplot)
          plotter->makeFrame(0., 40., 1e-2, 1., sPtV0, sRatio);
        else
          plotter->makeFrame(0., 40., 1e-9, 1e-1, sPtV0, sV0PtPerEvt);
      }
      break;
    case HistogramTypes::kJetPt:
    case HistogramTypes::kJetPtWrongCollision:
    case HistogramTypes::kJetZ:
    case HistogramTypes::kJetZWrongCollision:
      if (proj == ProjectionTypes::kPtV0) {
        if (inputs->ratioplot)
          plotter->makeFrame(0., inputs->ptjetmax, 1e-3, 1., sPtV0, sRatio);
        else
          plotter->makeFrame(0., inputs->ptjetmax, 1e-6, 1e1, sPtV0, sV0PtPerJet);
      } else if(proj == ProjectionTypes::kZV0) {
        if (inputs->ratioplot) {
          plotter->makeFrame(0., 1., 1e-3, 1., sZV0, sRatio);
        } else {
          plotter->makeFrame(0., 1., 5e-4, 1e1, sZV0, sV0ZPerJet);
        }
      }
    default: break;
  }
  if (plotter->frame)
    return;

  inputs->printLog("PileUp::autoFrame() Error: invalid combination of histogram type " + to_string(type) + " and projection type " + to_string(proj), VerbosityLevels::kErrors);
}

void PileUp::autoLatex(HistogramTypes::HistType type, ProjectionTypes::ProjType proj) {
  bool success = false;
  switch (type) {
    case HistogramTypes::kInclusive:
    case HistogramTypes::kInclusiveWrongCollision:
      if (proj == ProjectionTypes::kPtV0) {
        success = true;
        if (inputs->ratioplot) {
          plotter->addLatex(0.25, 0.80, sThisThesis + ", " + sAliceData);
          plotter->addLatex(0.25, 0.75, sSqrtS);
        } else {
          plotter->addLatex(0.47, 0.80, sThisThesis + ", " + sAliceData);
          plotter->addLatex(0.47, 0.75, sSqrtS);
        }
      }
      break;
    case HistogramTypes::kJetPt:
    case HistogramTypes::kJetPtWrongCollision:
    case HistogramTypes::kJetZ:
    case HistogramTypes::kJetZWrongCollision:
      if (proj == ProjectionTypes::kPtV0 || proj == ProjectionTypes::kZV0) {
        success = true;
        if (inputs->ratioplot) {
          plotter->addLatex(0.25, 0.80, sThisThesis + ", " + sAliceData);
          plotter->addLatex(0.25, 0.75, sSqrtS);
          plotter->addLatex(0.25, 0.70, sAntikt + " " + sChV0Jets);
          plotter->addLatex(0.25, 0.65, sRadius + ", " + sEtaJetRange);
          plotter->addLatex(0.25, 0.60, getPtJetRangeString(inputs->ptjetmin, inputs->ptjetmax));
        } else {
          plotter->addLatex(0.47, 0.83, sThisThesis + ", " + sAliceData);
          plotter->addLatex(0.47, 0.78, sSqrtS);
          plotter->addLatex(0.47, 0.73, sAntikt + " " + sChV0Jets);
          plotter->addLatex(0.47, 0.68, sRadius + ", " + sEtaJetRange);
          plotter->addLatex(0.47, 0.63, getPtJetRangeString(inputs->ptjetmin, inputs->ptjetmax));
        }
      }
      break;
    default: break;
  }
  if (!success)
    inputs->printLog("PileUp::autoLatex() Error: invalid combination of histogram type " + to_string(type) + "(" + HistogramTypes::to_string(type) + ")" + " and projection type " + to_string(proj) + "(" + ProjectionTypes::to_string(proj) + ")", VerbosityLevels::kErrors);
}

string PileUp::getHistName(HistogramTypes::HistType x) {
  string s = "jet-v0qa/";
  switch (x) {
    case HistogramTypes::kInclusive:
      s += "inclusive/" + inputs->hadron + "PtEtaMass";
      break;
    case HistogramTypes::kInclusiveWrongCollision:
      s += "inclusive/" + inputs->hadron + "PtEtaMassWrongCollision";
      break;
    case HistogramTypes::kJetPt:
      s += "jets/JetPtEta" + inputs->hadron + "Pt";
      break;
    case HistogramTypes::kJetPtWrongCollision:
      s += "jets/JetPtEta" + inputs->hadron + "PtWrongCollision";
      break;
    case HistogramTypes::kJetZ:
      s += "jets/JetPtEta" + inputs->hadron + "Z";
      break;
    case HistogramTypes::kJetZWrongCollision:
      s += "jets/JetPtEta" + inputs->hadron + "ZWrongCollision";
      break;
    case HistogramTypes::kMatchedJetPt:
      s += "jets/JetsPtEta" + inputs->hadron + "Pt";
      break;
    case HistogramTypes::kMatchedJetPtWrongCollision:
      s += "jets/JetsPtEta" + inputs->hadron + "PtWrongCollision";
      break;
    case HistogramTypes::kMatchedJetZ:
      s += "jets/JetsPtEta" + inputs->hadron + "Z";
      break;
    case HistogramTypes::kMatchedJetZWrongCollision:
      s += "jets/JetsPtEta" + inputs->hadron + "ZWrongCollision";
      break;
    default:
      inputs->printLog("PileUp::getHistName() Error: invalid histogram type " + to_string(x), VerbosityLevels::kErrors);
      s = "";
  }
  return s;
}

template <typename T>
T PileUp::getHistFromFile(HistogramTypes::HistType type) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    inputs->printLog("PileUp::getHistFromFile() Error: could not open file " + inputs->inputFileName, VerbosityLevels::kErrors);
    return nullptr;
  }
  T hist = (T)file->Get(getHistName(type).c_str());
  if (!hist) {
    inputs->printLog("PileUp::getHistFromFile() Error: could not find histogram " + getHistName(type) + " in file " + inputs->inputFileName, VerbosityLevels::kErrors);
    return nullptr;
  }
  return hist;
}

TH1* PileUp::getHistProjection(HistogramTypes::HistType type, ProjectionTypes::ProjType proj, string name) {
  if (name == "") {
    inputs->printLog("PileUp::getHistProjection() Error: histogram name is empty", VerbosityLevels::kErrors);
    return nullptr;
  }
  const int thnAxisMcpJetPt = 0, thnAxisMcdJetPt = 1, thnAxisMcdJetEta = 2, thnAxisV0 = 3;
  array<int, 2> ptBins, etaBins;
  TH1* h = nullptr;

  switch (type) {
    case HistogramTypes::kInclusive:
    case HistogramTypes::kInclusiveWrongCollision: {
      TH3* h3 = (TH3*)getHistFromFile<TH3*>(type);
      ptBins = getProjectionBins(h3->GetXaxis(), inputs->ptmin, inputs->ptmax);
      etaBins = getProjectionBins(h3->GetYaxis(), inputs->etamin, inputs->etamax);

      if (proj == ProjectionTypes::kPtV0)
        h = (TH1*)h3->ProjectionX(name.c_str(), etaBins[0], etaBins[1], 0, 1 + h3->GetNbinsZ());
      else if (proj == ProjectionTypes::kMass)
        h = (TH1*)h3->ProjectionZ(name.c_str(), ptBins[0], ptBins[1], etaBins[0], etaBins[1]);
    } break;
    case HistogramTypes::kJetPt:
    case HistogramTypes::kJetPtWrongCollision:
    case HistogramTypes::kJetZ:
    case HistogramTypes::kJetZWrongCollision: {
      TH3* h3 = (TH3*)getHistFromFile<TH3*>(type);
      ptBins = getProjectionBins(h3->GetXaxis(), inputs->ptjetmin, inputs->ptjetmax);
      etaBins = getProjectionBins(h3->GetYaxis(), inputs->etamin, inputs->etamax);

      if (proj == ProjectionTypes::kPtV0 || proj == ProjectionTypes::kZV0)
        h = (TH1*)h3->ProjectionZ(name.c_str(), ptBins[0], ptBins[1], etaBins[0], etaBins[1]);
    } break;
    case HistogramTypes::kMatchedJetPt:
    case HistogramTypes::kMatchedJetPtWrongCollision:
    case HistogramTypes::kMatchedJetZ:
    case HistogramTypes::kMatchedJetZWrongCollision: {
      THnSparse* hn = (THnSparse*)getHistFromFile<THnSparse*>(type);
      ptBins = getProjectionBins(hn->GetAxis(thnAxisMcdJetPt), inputs->ptjetmin, inputs->ptjetmax);
      etaBins = getProjectionBins(hn->GetAxis(thnAxisMcdJetEta), inputs->etamin, inputs->etamax);

      hn->GetAxis(thnAxisMcdJetPt)->SetRange(ptBins[0], ptBins[1]);
      hn->GetAxis(thnAxisMcdJetEta)->SetRange(etaBins[0], etaBins[1]);
      if (proj == ProjectionTypes::kPtV0 || proj == ProjectionTypes::kZV0) {
        h = (TH1*)hn->Projection(thnAxisV0);
        h->SetName(name.c_str());
      }
    } break;
    default:
      inputs->printLog("PileUp::getHistProjection() Error: invalid histogram type " + to_string(type), VerbosityLevels::kErrors);
      return nullptr;
  }
  if (!h)
    inputs->printLog("PileUp::getHistProjection() Error: could not create projection " + to_string(proj) + " for histogram type " + to_string(type), VerbosityLevels::kErrors);
  return h;
}

double PileUp::getNevts() {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    inputs->printLog("PileUp::getNevts() Error: could not open file " + inputs->inputFileName, VerbosityLevels::kErrors);
    return -1;
  }
  TH1* hist = (TH1*)file->Get("jet-v0qa/inclusive/hEvents");
  if (!hist) {
    inputs->printLog("PileUp::getNevts() Error: could not find histogram inclusive/hEvents in file " + inputs->inputFileName, VerbosityLevels::kErrors);
    return -1;
  }
  return hist->GetBinContent(2); // Should be bin 3 for sum of weights, i.e. xsec
}

double PileUp::getNjets(double ptjetmin, double ptjetmax) {
  // TODO: Now only accounts for events with V0s. Need V0 spectra task for proper normalisation.
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    inputs->printLog("PileUp::getNjets() Error: could not open file " + inputs->inputFileName, VerbosityLevels::kErrors);
    return -1;
  }
  THnSparse* thn = (THnSparse*)file->Get("jet-finder-v0-mcd-charged/hJet");
  if (!thn) {
    inputs->printLog("PileUp::getNjets() Error: could not find histogram jet-finder-v0-mcd-charged/hJet in file " + inputs->inputFileName, VerbosityLevels::kErrors);
    return -1;
  }
  const int ptAxis = 1, etaAxis = 2;
  array<int, 2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptjetmin, ptjetmax);
  array<int, 2> etaBins = getProjectionBins(thn->GetAxis(etaAxis), inputs->etamin, inputs->etamax);
  thn->GetAxis(etaAxis)->SetRange(etaBins[0], etaBins[1]);
  TH1* hist = (TH1*)thn->Projection(ptAxis);
  return hist->Integral(ptBins[0], ptBins[1]);
}

void PileUp::plotInclusivePt(bool doRatio) {
  inputs->logplot = true;
  inputs->ratioplot = doRatio;
  inputs->setEta(-0.9, 0.9);

  inputs->outputFileName = "pileUp";
  if (doRatio)
    inputs->outputFileName += "_ratio";
  inputs->outputFileName += ".pdf";

  TH1* inclPt   = (TH1*)getHistProjection(HistogramTypes::kInclusive, ProjectionTypes::kPtV0, "hPt");
  TH1* inclPtPU = (TH1*)getHistProjection(HistogramTypes::kInclusiveWrongCollision, ProjectionTypes::kPtV0, "hPtPU");

  inclPt = rebinHist(inclPt, rebinnedV0PtHist(inputs->hadron, "hPtRebinned"));
  inclPtPU = rebinHist(inclPtPU, rebinnedV0PtHist(inputs->hadron, "hPtPURebinned"));

  double nEvts = getNevts();
  if (nEvts < 0) {
    inputs->printLog("PileUp::plotpileup() Error: could not get number of events", VerbosityLevels::kErrors);
    return;
  }
  inclPt->Scale(1.0 / nEvts, "width");
  inclPtPU->Scale(1.0 / nEvts, "width");

  plotter->hists.push_back(inclPt);
  plotter->hists.push_back(inclPtPU);
  plotter->setHistStyles();

  if (!plotter->frame)
    autoFrame(HistogramTypes::kInclusive, ProjectionTypes::kPtV0);
  if (plotter->objects.empty())
    autoLatex(HistogramTypes::kInclusive, ProjectionTypes::kPtV0);

  if (!inputs->ratioplot) {
    plotter->makeLegend(0.7, 0.8, 0.5, 0.6, "");
    plotter->legend->AddEntry(inclPt, "Total");
    plotter->legend->AddEntry(inclPtPU, "Pile-up");
  }
  plotter->plot();
}

void PileUp::plotInJet(bool doRatio, bool doZ) {
  inputs->logplot = true;
  inputs->ratioplot = doRatio;
  inputs->setEta(-0.5, 0.5);

  inputs->outputFileName = inputs->getNameFromJetPt("pileUp");
  if (doZ)
    inputs->outputFileName += "z";
  if (doRatio)
    inputs->outputFileName += "_ratio";
  inputs->outputFileName += ".pdf";

  ProjectionTypes::ProjType proj = doZ ? ProjectionTypes::kZV0 : ProjectionTypes::kPtV0;
  HistogramTypes::HistType histType = doZ ? HistogramTypes::kJetZ : HistogramTypes::kJetPt;
  HistogramTypes::HistType histTypePU = doZ ? HistogramTypes::kJetZWrongCollision : HistogramTypes::kJetPtWrongCollision;
  TH1* inclPt   = (TH1*)getHistProjection(histType, proj, "hPt");
  TH1* inclPtPU = (TH1*)getHistProjection(histTypePU, proj, "hPtPU");

  if (doZ) {
    inclPt = rebinHist(inclPt, rebinnedV0ZHist("hPtRebinned"));
    inclPtPU = rebinHist(inclPtPU, rebinnedV0ZHist("hPtPURebinned"));
  } else {
    inclPt = rebinHist(inclPt, rebinnedV0PtHist(inputs->hadron, "hPtRebinned"));
    inclPtPU = rebinHist(inclPtPU, rebinnedV0PtHist(inputs->hadron, "hPtPURebinned"));
  }

  double nJets = getNjets(inputs->ptjetmin, inputs->ptjetmax);
  if (nJets < 0) {
    inputs->printLog("PileUp::plotInJet() Error: could not get number of jets", VerbosityLevels::kErrors);
    return;
  }
  inclPt->Scale(1.0 / nJets, "width");
  inclPtPU->Scale(1.0 / nJets, "width");

  plotter->hists.push_back(inclPt);
  plotter->hists.push_back(inclPtPU);
  plotter->setHistStyles();

  if (!plotter->frame)
    autoFrame(histType, proj);
  if (plotter->objects.empty())
    autoLatex(histType, proj);

  // if (inputs->ratioplot) {
  //   plotter->makeFrame(0., inputs->ptjetmax, 1e-3, 1., sPtV0, sRatio);
  //   if (doZ)
  //     plotter->makeFrame(0., 1., 1e-3, 1., sZV0, sRatio);

  //   plotter->addLatex(0.25, 0.80, sThisThesis + ", " + sAliceData);
  //   plotter->addLatex(0.25, 0.75, sSqrtS);
  //   plotter->addLatex(0.25, 0.70, sAntikt + " " + sChV0Jets);
  //   plotter->addLatex(0.25, 0.65, sRadius + ", " + sEtaJetRange);
  //   plotter->addLatex(0.25, 0.60, getPtJetRangeString(inputs->ptjetmin, inputs->ptjetmax));
  // } else {
  //   plotter->makeFrame(0., inputs->ptjetmax, 1e-6, 1e1, sPtV0, getOneOverString(sNjets) + " " + getdYdXString(sNV0, sPtV0));
  //   if (doZ)
  //     plotter->makeFrame(0., 1., 5e-4, 10, sZV0, getOneOverString(sNjets) + " " + getdYdXString(sNV0, sZV0));

  //   plotter->makeLegend(0.7, 0.8, 0.5, 0.6, "");
  //   plotter->legend->AddEntry(inclPt, "Total");
  //   plotter->legend->AddEntry(inclPtPU, "Pile-up");

  //   plotter->addLatex(0.47, 0.83, sThisThesis + ", " + sAliceData);
  //   plotter->addLatex(0.47, 0.78, sSqrtS);
  //   plotter->addLatex(0.47, 0.73, sAntikt + " " + sChV0Jets);
  //   plotter->addLatex(0.47, 0.68, sRadius + ", " + sEtaJetRange);
  //   plotter->addLatex(0.47, 0.63, getPtJetRangeString(inputs->ptjetmin, inputs->ptjetmax));
  // }
  plotter->plot();
}

PileUp setuppu() {
  PileUp pu;
  pu.inputs->verbosity = VerbosityLevels::kDebug;
  pu.inputs->hadron = "K0S";
  pu.inputs->train = 496209;
  pu.inputs->setInputFileNameFromTrain();
  pu.inputs->logplot = true;
  return pu;
}

void printOptions() {
  cout << "Histogram types:\n";
  for (int i = 0; i < 20; i++) {
    if (HistogramTypes::is_valid(i))
      cout << (HistogramTypes::HistType)i << ": " << HistogramTypes::to_string((HistogramTypes::HistType)i) << endl;
  }
  cout << "Projection types:\n";
  for (int i = 0; i < 10; i++) {
    if (ProjectionTypes::is_valid(i))
      cout << (ProjectionTypes::ProjType)i << ": " << ProjectionTypes::to_string((ProjectionTypes::ProjType)i) << endl;
  }
}
printOptions();

void printHistContents(TH1* h, string s) {
  cout << s << "\n";
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    cout << TString::Format("%.1f - %.1f: %f (%f)\n", h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinUpEdge(i), 1 - h->GetBinContent(i), h->GetBinContent(i)).Data();
  }
}

void plotpileup() {
  const bool doRatio = true;
  const bool doZ = true;

  PileUp pu = setuppu();
  pu.inputs->verbosity = VerbosityLevels::kErrors;
  pu.plotter->makeFrame(0., 40., 1e-12, 1e-1, sPtV0, sV0PtPerEvt);
  pu.plotInclusivePt(!doRatio);
  pu.plotter->reset();
  pu.plotInclusivePt(doRatio);
  printHistContents(pu.plotter->hists[1], "Pile-up inclusive");
  pu.plotter->reset();

  pu.inputs->setJetPt(10., 20.);
  // pu.plotter->makeFrame(0., pu.inputs->ptjetmax, 1e-7, 10, sPtV0, sV0PtPerJet);
  // pu.plotInJet(!doRatio, !doZ);
  // pu.plotter->reset();
  pu.plotter->makeFrame(0., pu.inputs->ptjetmax, 1e-3, 1., sPtV0, sRatio);
  pu.plotInJet(doRatio, !doZ);
  printHistContents(pu.plotter->hists[1], "Pile-up in jets 10-20 GeV/c");
  pu.plotter->reset();
  // pu.plotter->makeFrame(0., 1., 1e-5, 10, sZV0, sV0ZPerJet);
  // pu.plotInJet(!doRatio, doZ);
  // pu.plotter->reset();
  // pu.plotter->makeFrame(0., 1., 1e-3, 1., sZV0, sRatio);
  // pu.plotInJet(doRatio, doZ);
  // pu.plotter->reset();

  pu.inputs->setJetPt(20., 30.);
  // pu.plotter->makeFrame(0., pu.inputs->ptjetmax, 1e-7, 10, sPtV0, sV0PtPerJet);
  // pu.plotInJet(!doRatio, !doZ);
  // pu.plotter->reset();
  pu.plotter->makeFrame(0., pu.inputs->ptjetmax, 1e-3, 1., sPtV0, sRatio);
  pu.plotInJet(doRatio, !doZ);
  printHistContents(pu.plotter->hists[1], "Pile-up in jets 20-30 GeV/c");
  pu.plotter->reset();
  // pu.plotter->makeFrame(0., 1., 1e-5, 10, sZV0, sV0ZPerJet);
  // pu.plotInJet(!doRatio, doZ);
  // pu.plotter->reset();
  // pu.plotInJet(doRatio, doZ);
  // pu.plotter->reset();

  pu.inputs->setJetPt(30., 40.);
  // pu.plotter->makeFrame(0., pu.inputs->ptjetmax, 1e-7, 10, sPtV0, sV0PtPerJet);
  // pu.plotInJet(!doRatio, !doZ);
  // pu.plotter->reset();
  pu.plotter->makeFrame(0., pu.inputs->ptjetmax, 1e-3, 1., sPtV0, sRatio);
  pu.plotInJet(doRatio, !doZ);
  printHistContents(pu.plotter->hists[1], "Pile-up in jets 30-40 GeV/c");
  pu.plotter->reset();
  // pu.plotter->makeFrame(0., 1., 1e-5, 10, sZV0, sV0ZPerJet);
  // pu.plotInJet(!doRatio, doZ);
  // pu.plotter->reset();
  // pu.plotter->makeFrame(0., 1., 1e-3, 1., sZV0, sRatio);
  // pu.plotInJet(doRatio, doZ);
  // pu.plotter->reset();
}

#endif
