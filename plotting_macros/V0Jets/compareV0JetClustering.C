
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TString.h"
#include "TLegend.h"

#include "../histUtils.C"
#include "../plotUtils.C"
#include "../myStrings.C"

// This script compares jets when we assign Lambdas the K0S mass or the Lambda mass:
// * Jet pt spectrum
// * K0S z spectrum
// * Lambda z spectrum

using namespace histutils;

namespace verbosityutilities {
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
  bool passVerbosityCheck(Verbosity level, Verbosity threshold) {
    return (is_valid(level) && level <= threshold);
  }
  void printLog(string message) {
    cout << message << endl;
  }
} // namespace verbosityutilities

namespace analysisutilities {
  enum AnaType {kJetLasK, kJetSchemes, kZK0SLasK, kZK0SSchemesV0, kZK0SSchemesK0, kZLambdaLasK, kZLambdaSchemesV0, kZLambdaSchemesK0};
}

struct Names {
  public:
    string inputFileName;
    string histogramName;
    string legendEntry;
    string jetHistName;

    Names() {}
    Names(string file, string hist, string legend) : inputFileName(file), histogramName(hist), legendEntry(legend) {}
    Names(string file, string hist, string legend, string jetHist) : inputFileName(file), histogramName(hist), legendEntry(legend), jetHistName(jetHist) {}
};

struct InputSettings{
  private:
    verbosityutilities::Verbosity verbosity = verbosityutilities::kInfo;

    string getNameFromVar(string prefix, string varstring, string suffix) {
      return prefix + varstring + suffix;
    }
    template <typename T> bool setVariable(T a, T b, T &x, T &y);
  public:
    // Analysis settings
    vector<Names> names = {};
    string hadron = "";
    string histName = "";
    string inputFileName = "";
    string outputFileName = "";

    double ptmin = -1e3, ptmax = -1e3, zmin = -1e3, zmax = -1e3;
    double ptjetmin = -1e3, ptjetmax = -1e3;
    double etamin = -0.9, etamax = 0.9;
    vector<vector<double>> ptBinEdges = {};

    // Plotting settings
    double textSize = 0.04;
    bool logplot = false;
    bool ratioplot = false;

    // Setters and getters
    verbosityutilities::Verbosity getVerbosity() { return verbosity; }

    void setVerbosity(verbosityutilities::Verbosity v) {
      if (verbosityutilities::is_valid(v)) {
        verbosity = v;
      } else {
        verbosityutilities::printLog("InputSettings::setVerbosity() Error: invalid verbosity");
      }
    }
    bool setEta(double a, double b) {   return setVariable(a, b, etamin, etamax); }
    bool setJetPt(double a, double b) { return setVariable(a, b, ptjetmin, ptjetmax); }
    bool setPt(double a, double b) {    return setVariable(a, b, ptmin, ptmax); }

    // Functions
    void addNames(string file, string hist, string legend) { names.push_back({file, hist, legend}); }
    void addNames(string file, string hist, string legend, string jethist) { names.push_back({file, hist, legend, jethist}); }
    string getNameFromJetPt(string prefix, string suffix) {
      return getNameFromVar(prefix, TString::Format("ptjet%.f-%.f", ptjetmin, ptjetmax).Data(), suffix);
    }
    bool passVerbosityCheck(verbosityutilities::Verbosity v) {
      return verbosityutilities::passVerbosityCheck(v, verbosity);
    }
    void printLog(string message, verbosityutilities::Verbosity messageVerbLevel) {
      if (passVerbosityCheck(messageVerbLevel))
        verbosityutilities::printLog(message);
    }
};

template <typename T>
bool InputSettings::setVariable(T a, T b, T &x, T &y) {
  if (a > b) {
    printLog("InputSettings::setVariable() Error: min > max", verbosityutilities::kErrors);
    return false;
  }
  x = a;
  y = b;
  return true;
}

double getNevts(Names names, bool withV0s = false) {
  TFile* f = TFile::Open(names.inputFileName.c_str());
  TH1D* hNEvents = (TH1D*)f->Get("hNEvts");
  return hNEvents->GetBinContent(1 + (int)withV0s);
}

double getNjets(InputSettings& inputs, Names names) {
  TFile* f = TFile::Open(names.inputFileName.c_str());
  if (!f) {
    inputs.printLog(TString::Format("getNjets() Error: Could not open file %s", names.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return -1;
  }
  TH3* jets = (TH3*)f->Get(names.jetHistName.c_str());
  if (!jets) {
    inputs.printLog(TString::Format("getNjets() Error: Could not find histogram %s in file %s", names.jetHistName.c_str(), names.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return -1;
  }
  array<int, 2> ptbins  = getProjectionBins(jets->GetXaxis(), inputs.ptjetmin, inputs.ptjetmax);
  array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), inputs.etamin, inputs.etamax);
  TH1* jetpt = jets->ProjectionX("jetpt", etabins[0], etabins[1]);
  return jetpt->Integral(ptbins[0], ptbins[1]);
}

std::array<double, 2> getNjetsAndError(InputSettings& inputs, Names names) {
  TFile* f = TFile::Open(names.inputFileName.c_str());
  if (!f) {
    inputs.printLog(TString::Format("getNjetsAndError() Error: Could not open file %s", names.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return {-1, -1};
  }
  TH3* jets = (TH3*)f->Get(names.jetHistName.c_str());
  if (!jets) {
    inputs.printLog(TString::Format("getNjetsAndError() Error: Could not find histogram %s in file %s", names.jetHistName.c_str(), names.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return {-1, -1};
  }
  array<int, 2> ptbins  = getProjectionBins(jets->GetXaxis(), inputs.ptjetmin, inputs.ptjetmax);
  array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), inputs.etamin, inputs.etamax);
  TH1* jetpt = jets->ProjectionX("jetpt", etabins[0], etabins[1]);
  double error;
  double njets = jetpt->IntegralAndError(ptbins[0], ptbins[1], error);
  return {njets, error};
}

TH1* getJetPtHist(InputSettings& inputs, Names name, const TH1* rebinTemplate = nullptr) {
  TFile *inFile = TFile::Open(name.inputFileName.c_str());
  if (!inFile) {
    inputs.printLog(TString::Format("getJetPtHist() Error: Could not open file %s", name.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return nullptr;
  }
  TH3* jets = (TH3*)inFile->Get(name.histogramName.c_str());
  if (!jets) {
    inputs.printLog(TString::Format("getJetPtHist() Error: Could not find histogram %s in file %s", name.histogramName.c_str(), name.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return nullptr;
  }

  array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), inputs.etamin, inputs.etamax);
  TH1* th1 = jets->ProjectionX("jetpt", etabins[0], etabins[1]);
  th1->SetName(TString::Format("jetpt_%s_%s", name.histogramName.c_str(), name.inputFileName.c_str()).Data());
  th1->Sumw2();
  if (rebinTemplate) {
    TH1* tmp = rebinHist(th1, rebinTemplate);
    th1 = tmp;
  }
  return th1;
}

TH1* getV0ZHist(InputSettings& inputs, Names name, const TH1* rebinTemplate = nullptr) {
  TFile *inFile = TFile::Open(name.inputFileName.c_str());
  if (!inFile) {
    inputs.printLog(TString::Format("getV0ZHist() Error: Could not open file %s", name.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return nullptr;
  }
  THnSparseD* thn = (THnSparseD*)inFile->Get(name.histogramName.c_str());
  if (!thn) {
    inputs.printLog(TString::Format("getV0ZHist() Error: Could not find histogram %s in file %s", name.histogramName.c_str(), name.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return nullptr;
  }

  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int jetetaAxis = 1;
  const int jetphiAxis = 2;
  const int zAxis      = 3;

  array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), inputs.ptjetmin, inputs.ptjetmax);
  array<int, 2> etabins   = getProjectionBins(thn->GetAxis(jetetaAxis), inputs.etamin, inputs.etamax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(jetetaAxis)->SetRange(etabins[0], etabins[1]);

  TH1* th1 = (TH1*)thn->Projection(zAxis);
  th1->SetName(TString::Format("z_%s_%s", name.histogramName.c_str(), name.inputFileName.c_str()).Data());
  th1->Sumw2();
  if (rebinTemplate) {
    TH1* tmp = rebinHist(th1, rebinTemplate);
    th1 = tmp;
  }
  return th1;
}

// --------------------------------------------------------
// Compare jet pt between E and pt recombination schemes 
// K0: All V0s have K0 mass
// V0: V0s have correct mass
// --------------------------------------------------------
void compareJetPtSchemesK0() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.setEta(-0.35, 0.35);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string jetHistName = "hK0Jet";
  x.addNames(eFileName, jetHistName, mystrings::sEscheme);
  x.addNames(ptFileName, jetHistName, mystrings::sPtscheme);

  x.outputFileName = "schemesK0-jetpt_Comparison.pdf";
  x.ratioplot = false;
  x.logplot = true;

  TH1* jetPtTemplate = new TH1D("jetpttemplate", "jetpttemplate", 12, 0., 60.);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);

  p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  for (int i = 0; i < x.names.size(); i++) {
    TH1* h = getJetPtHist(x, x.names[i], jetPtTemplate);
    h->Scale(1./getNevts(x.names[i]), "width");
    p.addHistogram(h);
  }

  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";
  p.addLatex(0.60, 0.80, mystrings::sThisThesis);
  p.addLatex(0.60, 0.75, mystrings::sPythiaSim);
  p.addLatex(0.60, 0.70, mystrings::sSqrtS + ", " + lambdaMassString);
  p.addLatex(0.60, 0.65, mystrings::sAntikt + " jets");
  p.addLatex(0.60, 0.60, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);

  string xTitle = mystrings::sPtJetWithUnits;
  string yTitle = mystrings::sJetsPerXsec;
  p.makeFrame(0., 60., 1e-8, 1e-3, xTitle, yTitle);

  p.setHistStyles();
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioJetPtSchemesK0() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.setEta(-0.35, 0.35);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string jetHistName = "hK0Jet";
  x.addNames(eFileName, jetHistName, mystrings::sEscheme);
  x.addNames(ptFileName, jetHistName, mystrings::sPtscheme);

  x.outputFileName = "schemesK0-jetpt_Ratio.pdf";
  x.ratioplot = true;
  x.logplot = false;

  TH1* jetPtTemplate = new TH1D("jetpttemplate", "jetpttemplate", 12, 0., 60.);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);

  p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  for (int i = 0; i < x.names.size(); i++) {
    TH1* h = getJetPtHist(x, x.names[i], jetPtTemplate);
    h->Scale(1./getNevts(x.names[i]), "width");
    p.addHistogram(h);
  }

  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";
  p.addLatex(0.275, 0.85, mystrings::sThisThesis);
  p.addLatex(0.275, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.75, mystrings::sSqrtS + ", " + lambdaMassString);
  p.addLatex(0.275, 0.70, mystrings::sAntikt + " jets");
  p.addLatex(0.275, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);

  string xTitle = mystrings::sPtJetWithUnits;
  string yTitle = mystrings::sRatio;
  p.makeFrame(0., 60., 0.8, 1.2, xTitle, yTitle);

  p.setHistStyles();
  p.makeRatios();
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void compareJetPtSchemesV0() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.setEta(-0.35, 0.35);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, jetHistName, mystrings::sEscheme);
  x.addNames(ptFileName, jetHistName, mystrings::sPtscheme);

  x.outputFileName = "schemesV0-jetpt_Comparison.pdf";
  x.ratioplot = false;
  x.logplot = true;

  TH1* jetPtTemplate = new TH1D("jetpttemplate", "jetpttemplate", 12, 0., 60.);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);

  for (int i = 0; i < x.names.size(); i++) {
    TH1* h = getJetPtHist(x, x.names[i], jetPtTemplate);
    h->Scale(1./getNevts(x.names[i]), "width");
    p.addHistogram(h);
  }

  p.addLatex(0.60, 0.80, mystrings::sThisThesis);
  p.addLatex(0.60, 0.75, mystrings::sPythiaSim);
  p.addLatex(0.60, 0.70, mystrings::sSqrtS);
  p.addLatex(0.60, 0.65, mystrings::sAntikt + " jets");
  p.addLatex(0.60, 0.60, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);

  string xTitle = mystrings::sPtJetWithUnits;
  string yTitle = mystrings::sJetsPerXsec;
  p.makeFrame(0., 60., 1e-8, 1e-3, xTitle, yTitle);

  p.setHistStyles();
  p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioJetPtSchemesV0() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.setEta(-0.35, 0.35);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, jetHistName, mystrings::sEscheme);
  x.addNames(ptFileName, jetHistName, mystrings::sPtscheme);

  x.outputFileName = "schemesV0-jetpt_Ratio.pdf";
  x.ratioplot = true;
  x.logplot = false;

  TH1* jetPtTemplate = new TH1D("jetpttemplate", "jetpttemplate", 12, 0., 60.);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);

  p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  for (int i = 0; i < x.names.size(); i++) {
    TH1* h = getJetPtHist(x, x.names[i], jetPtTemplate);
    h->Scale(1./getNevts(x.names[i]), "width");
    p.addHistogram(h);
  }

  p.addLatex(0.275, 0.85, mystrings::sThisThesis);
  p.addLatex(0.275, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.75, mystrings::sSqrtS);
  p.addLatex(0.275, 0.70, mystrings::sAntikt + " jets");
  p.addLatex(0.275, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);

  string xTitle = mystrings::sPtJetWithUnits;
  string yTitle = mystrings::sRatio;
  p.makeFrame(0., 60., 0.8, 1.2, xTitle, yTitle);

  p.setHistStyles();
  p.makeRatios();
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void jetPtSchemes() {
  gROOT->SetBatch(true);
  compareJetPtSchemesK0();
  ratioJetPtSchemesK0();
  compareJetPtSchemesV0();
  ratioJetPtSchemesV0();
}

// --------------------------------------------------------
// Compare V0 z between E and pt recombination schemes 
// K0: All V0s have K0 mass
// V0: V0s have correct mass
// --------------------------------------------------------

// TODO: Add manual error calculations!!!
// K0S
void compareZKSchemesV01020() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(10., 20.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Comparison.pdf");
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    std::array<double, 2> njetsAndError = getNjetsAndError(x, ns);
    double njets = njetsAndError[0];
    double njetsError = njetsAndError[1];
    double njetsRelError = njetsError / njets;

    TH1* hRelError = histutils::getHistRelErrors(h);
    for (int bin = 1; bin <= hRelError->GetNbinsX(); bin++) {
      double relError = hRelError->GetBinContent(bin);
      double newRelError = std::sqrt(relError * relError + njetsRelError * njetsRelError);
      hRelError->SetBinContent(bin, newRelError);
    }
    
    h->Scale(1./njets, "width");
    histutils::setHistErrors(h, hRelError);
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sK0SZPerJetXsec;
  p.makeFrame(1e-3, 1.+1e-3, 1e-2, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZKSchemesV01020() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(10., 20.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    // double njets = getNjets(x, ns);
    std::array<double, 2> njetsAndError = getNjetsAndError(x, ns);
    double njets = njetsAndError[0];
    double njetsError = njetsAndError[1];
    double njetsRelError = njetsError / njets;

    TH1* hRelError = histutils::getHistRelErrors(h);
    for (int bin = 1; bin <= hRelError->GetNbinsX(); bin++) {
      double relError = hRelError->GetBinContent(bin);
      double newRelError = std::sqrt(relError * relError + njetsRelError * njetsRelError);
      hRelError->SetBinContent(bin, newRelError);
    }
    
    h->Scale(1./njets, "width");
    histutils::setHistErrors(h, hRelError);
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sThisThesis);
  p.addLatex(0.275, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.75, mystrings::sSqrtS);
  p.addLatex(0.275, 0.70, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.60, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.65, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void compareZKSchemesV02030() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(20., 30.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Comparison.pdf");
  x.ratioplot = false;
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sK0SZPerJetXsec;
  p.makeFrame(1e-3, 1.+1e-3, 1e-2, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZKSchemesV02030() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(20., 30.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sThisThesis);
  p.addLatex(0.275, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.75, mystrings::sSqrtS);
  p.addLatex(0.275, 0.70, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.60, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.65, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void compareZKSchemesV03040() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(30., 40.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Comparison.pdf");
  x.ratioplot = false;
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sK0SZPerJetXsec;
  p.makeFrame(1e-3, 1.+1e-3, 1e-3, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZKSchemesV03040() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(30., 40.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0., 4., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sThisThesis);
  p.addLatex(0.275, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.75, mystrings::sSqrtS);
  p.addLatex(0.275, 0.70, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.60, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.65, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void zKSchemes() {
  gROOT->SetBatch(true);
  compareZKSchemesV01020();
  ratioZKSchemesV01020();
  compareZKSchemesV02030();
  ratioZKSchemesV02030();
  compareZKSchemesV03040();
  ratioZKSchemesV03040();
}

// Lambda
void compareZLSchemesV01020() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(10., 20.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Comparison.pdf");
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sK0SZPerJetXsec;
  p.makeFrame(1e-3, 1.+1e-3, 1e-2, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZLSchemesV01020() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(10., 20.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sThisThesis);
  p.addLatex(0.275, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.75, mystrings::sSqrtS);
  p.addLatex(0.275, 0.70, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.60, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.65, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void compareZLSchemesV02030() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(20., 30.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Comparison.pdf");
  x.ratioplot = false;
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sK0SZPerJetXsec;
  p.makeFrame(1e-3, 1.+1e-3, 1e-2, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZLSchemesV02030() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(20., 30.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sThisThesis);
  p.addLatex(0.275, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.75, mystrings::sSqrtS);
  p.addLatex(0.275, 0.70, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.60, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.65, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void compareZLSchemesV03040() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(30., 40.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Comparison.pdf");
  x.ratioplot = false;
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sK0SZPerJetXsec;
  p.makeFrame(1e-3, 1.+1e-3, 1e-3, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZLSchemesV03040() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(30., 40.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  string zHistName = "hzV0_" + x.hadron;
  string jetHistName = "hV0Jet";
  x.addNames(eFileName, zHistName, mystrings::sEscheme, jetHistName);
  x.addNames(ptFileName, zHistName, mystrings::sPtscheme, jetHistName);
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";

  x.outputFileName = x.getNameFromJetPt("schemesV0-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0., 4., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sThisThesis);
  p.addLatex(0.275, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.75, mystrings::sSqrtS);
  p.addLatex(0.275, 0.70, mystrings::sAntikt + " jets, " + lambdaMassString);
  p.addLatex(0.275, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.60, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.65, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void zLSchemes() {
  gROOT->SetBatch(true);
  compareZLSchemesV01020();
  ratioZLSchemesV01020();
  compareZLSchemesV02030();
  ratioZLSchemesV02030();
  compareZLSchemesV03040();
  ratioZLSchemesV03040();
}

// --------------------------------------------------------
// Compare jet pt when assigning K0 mass to Lambda
// This comparison only makes sense for E scheme
// K0: All V0s have K0 mass
// V0: V0s have correct mass
// --------------------------------------------------------

void compareJetPtLasK() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string fileName = filePath + "v0jetclustering-Escheme.root";
  string histNameV0 = "hV0Jet";
  string histNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";
  string kaonMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sK0S + ")";

  x.addNames(fileName, histNameV0, lambdaMassString);
  x.addNames(fileName, histNameK0, kaonMassString);

  x.setEta(-0.35, 0.35);

  x.outputFileName = "LasK-jetpt_Comparison.pdf";
  x.logplot = true;

  TH1* jetPtTemplate = new TH1D("jetpttemplate", "jetpttemplate", 12, 0., 60.);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);

  for (auto ns : x.names) {
    TH1* h = getJetPtHist(x, ns, jetPtTemplate);
    h->Scale(1./getNevts(ns), "width");
    p.addHistogram(h);
  }

  p.addLatex(0.50, 0.85, mystrings::sThisThesis);
  p.addLatex(0.50, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.50, 0.75, mystrings::sSqrtS);
  p.addLatex(0.50, 0.70, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.50, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);

  string xTitle = mystrings::sPtJetWithUnits;
  string yTitle = mystrings::sJetsPerXsec;
  p.makeFrame(0., 60., 1e-8, 1e-3, xTitle, yTitle);

  p.setHistStyles();
  p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  for (int i = 0; i < p.getHists().size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }

  p.plot();
}

void ratioJetPtLasK() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string fileName = filePath + "v0jetclustering-Escheme.root";
  string histNameV0 = "hV0Jet";
  string histNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sLambda + ")";
  string kaonMassString = mystrings::sLambda + " mass = #it{m}(" + mystrings::sK0S + ")";

  x.addNames(fileName, histNameV0, lambdaMassString);
  x.addNames(fileName, histNameK0, kaonMassString);

  x.setEta(-0.35, 0.35);

  x.outputFileName = "LasK-jetpt_Ratio.pdf";
  x.logplot = false;

  TH1* jetPtTemplate = new TH1D("jetpttemplate", "jetpttemplate", 12, 0., 60.);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);

  for (auto ns : x.names) {
    TH1* h = getJetPtHist(x, ns, jetPtTemplate);
    h->Scale(1./getNevts(ns), "width");
    p.addHistogram(h);
  }

  p.addLatex(0.275, 0.85, mystrings::sThisThesis);
  p.addLatex(0.275, 0.80, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.75, mystrings::sSqrtS);
  p.addLatex(0.275, 0.70, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.65, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);

  string xTitle = mystrings::sPtJetWithUnits;
  string yTitle = mystrings::sRatio;
  p.makeFrame(0., 60., 0.9, 1.1, xTitle, yTitle);

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  for (int i = 0; i < p.getHists().size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }

  p.plot();
}

void jetPtLasK() {
  gROOT->SetBatch(true);
  compareJetPtLasK();
  ratioJetPtLasK();
}

// --------------------------------------------------------
// Compare V0 z when assigning K0 mass to Lambda
// This comparison only makes sense for E scheme
// K0: All V0s have K0 mass
// V0: V0s have correct mass
// --------------------------------------------------------

// K0S
void compareZKLasK1020() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(10., 20.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Comparison.pdf");
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sK0SZPerJetXsec;
  p.makeFrame(1e-3, 1.+1e-3, 1e-2, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZKLasK1020() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(10., 20.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.80, mystrings::sSqrtS);
  p.addLatex(0.275, 0.75, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.70, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.65, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void compareZKLasK2030() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(20., 30.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Comparison.pdf");
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sK0SZPerJetXsec;
  p.makeFrame(1e-3, 1.+1e-3, 1e-2, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZKLasK2030() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(20., 30.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.80, mystrings::sSqrtS);
  p.addLatex(0.275, 0.75, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.70, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.65, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void compareZKLasK3040() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(30., 40.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Comparison.pdf");
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sK0SZPerJetXsec;
  p.makeFrame(1e-3, 1.+1e-3, 1e-3, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZKLasK3040() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "K0S";
  x.setJetPt(30., 40.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZK0S;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.80, mystrings::sSqrtS);
  p.addLatex(0.275, 0.75, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.70, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.65, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void zKLasK() {
  gROOT->SetBatch(true);
  compareZKLasK1020();
  ratioZKLasK1020();
  compareZKLasK2030();
  ratioZKLasK2030();
  compareZKLasK3040();
  ratioZKLasK3040();
}

// Lambda
void compareZLLasK1020() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(10., 20.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Comparison.pdf");
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZLambda;
  string yTitle = mystrings::sLambdaZPerJet;
  p.makeFrame(1e-3, 1.+1e-3, 1e-2, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZLLasK1020() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(10., 20.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZLambda;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.80, mystrings::sSqrtS);
  p.addLatex(0.275, 0.75, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.70, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.65, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void compareZLLasK2030() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(20., 30.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Comparison.pdf");
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZLambda;
  string yTitle = mystrings::sLambdaZPerJet;
  p.makeFrame(1e-3, 1.+1e-3, 1e-3, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZLLasK2030() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(20., 30.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZLambda;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.80, mystrings::sSqrtS);
  p.addLatex(0.275, 0.75, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.70, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.65, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void compareZLLasK3040() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(30., 40.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Comparison.pdf");
  x.logplot = true;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZLambda;
  string yTitle = mystrings::sLambdaZPerJet;
  p.makeFrame(1e-3, 1.+1e-3, 1e-3, 10, xTitle, yTitle);

  p.addLatex(0.275, 0.45, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.40, mystrings::sSqrtS);
  p.addLatex(0.275, 0.35, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.30, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.25, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void ratioZLLasK3040() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.hadron = "Lambda0";
  x.setJetPt(30., 40.);
  x.setEta(-0.35, 0.35);
  
  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string zHistNameV0 = "hzV0_" + x.hadron;
  string zHistNameK0 = "hzK0_" + x.hadron;
  string jetHistNameV0 = "hV0Jet";
  string jetHistNameK0 = "hK0Jet";
  string lambdaMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sLambda);
  string kaonMassString = mystrings::sLambda + " mass = " + mystrings::massOf(mystrings::sK0S);

  x.addNames(eFileName, zHistNameV0, lambdaMassString, jetHistNameV0);
  x.addNames(eFileName, zHistNameK0, kaonMassString, jetHistNameK0);

  x.outputFileName = x.getNameFromJetPt("LasK-z" + x.hadron + "_", "_Ratio.pdf");
  x.logplot = false;

  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  plotutils::Plotter p(x.outputFileName, x.logplot, 0.04);
  for (auto ns : x.names) {
    TH1* h = getV0ZHist(x, ns, zTemplate);
    double njets = getNjets(x, ns);
    h->Scale(1./njets, "width");
    p.addHistogram(h);
  }

  string xTitle = mystrings::sZLambda;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1.+1e-3, 0.5, 2., xTitle, yTitle);

  p.addLatex(0.275, 0.85, mystrings::sPythiaSim);
  p.addLatex(0.275, 0.80, mystrings::sSqrtS);
  p.addLatex(0.275, 0.75, mystrings::sAntiktJets + ", " + mystrings::sEscheme);
  p.addLatex(0.275, 0.70, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(0.275, 0.65, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax, true));

  p.setHistStyles();
  p.makeRatios();
  p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
  for (int i = 0; i < x.names.size(); i++) {
    p.addLegendEntry(p.getHists()[i], x.names[i].legendEntry.c_str());
  }
  p.plot();
}

void zLLasK() {
  gROOT->SetBatch(true);
  compareZLLasK1020();
  ratioZLLasK1020();
  compareZLLasK2030();
  ratioZLLasK2030();
  compareZLLasK3040();
  ratioZLLasK3040();
}
