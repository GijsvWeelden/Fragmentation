
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

#ifndef SUBTRACTION_COMPARISON_C
#define SUBTRACTION_COMPARISON_C

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
  string getPtJetRangeString(double ptmin, double ptmax);

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

  const string sPtJet       = getPtString(sJet);
  const string sPtV0        = getPtString(sV0);
  const string sPtK0S       = getPtString(sK0S);
  const string sZV0         = getZString(sV0);
  const string sZK0S        = getZString(sK0S);

  const string sJetsPerEvent     = getOneOverString(sNevts) + " " + getdYdXString(sNjets, sPtJet);
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
string MyStrings::getPtJetRangeString(double ptmin, double ptmax) {
  return getVarRangeString(ptmin, sPtJet, ptmax);
}

namespace MyEnums {
  enum Verbosity {kErrors, kWarnings, kInfo, kDebug};
  enum Subtraction {kNone, kSubtract, kAddTracks, kUseFlags};
  enum HistType {kEvents, kJet, kV0Pt, kV0InJetPt, kV0InJetZ};
}

using namespace MyStrings;
using namespace MyEnums;

struct InputSettings {
  private:
  public:
    int train, rebinNumber, verbosity = kWarnings, subtraction = -1;
    string inputFileName, outputFileName;
    string hadron;

    double jetptmin, jetptmax, jetptlow, jetpthigh;
    double ptmin, ptmax, ptlow, pthigh;
    double etaMin = -0.9, etaMax = 0.9;

    bool logplot = false, ratioplot = false;
    bool normalisePerJet = true;

    double getMass();
    string getSaveNameFromPt(string prefix, string suffix);
    string getSaveNameFromJetPt(string prefix, string suffix);
    string printLog(string message, int verbThreshold);
    string setInputFileNameFromTrain();
    void setJetPt(double a, double b);
    void setPt(double a, double b);
    void setEta(double a, double b);
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
  string s = TString::Format("%s_pt%.1f-%.1f%s", prefix.c_str(), ptlow, pthigh, suffix.c_str()).Data();
  return s;
}

string InputSettings::getSaveNameFromJetPt(string prefix, string suffix = "") {
  string s = TString::Format("%s_jetpt%.1f-%.1f%s", prefix.c_str(), jetptlow, jetpthigh, suffix.c_str()).Data();
  return s;
}

string InputSettings::printLog(string message, int verbThreshold) {
  if (verbosity < verbThreshold)
    return "";

  cout << message << endl;
  return message;
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  inputFileName = s;
  return s;
}

void InputSettings::setJetPt(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setJetPt() Error: jetptmin > jetptmax";
    printLog(s, kErrors);
    return;
  }
  jetptmin = a;
  jetptmax = b;
  jetptlow = a;
  jetpthigh = b;
}

void InputSettings::setPt(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setPt() Error: ptmin > ptmax";
    printLog(s, kErrors);
    return;
  }
  ptmin = a;
  ptmax = b;
  ptlow = a;
  pthigh = b;
}

void InputSettings::setEta(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setEta() Error: etaMin > etaMax";
    printLog(s, kErrors);
    return;
  }
  etaMin = a;
  etaMax = b;
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

    Plotter() { inputs = new InputSettings(); }
    Plotter(InputSettings& x) { inputs = &x; }
    Plotter(InputSettings* x) { inputs = x; }

    void addLatex(double x, double y, string s);
    void addLine(double x0, double y0, double x1, double y1, int styleNumber, int lineStyle, int lineWidth);
    void makeCanvas(string s, double x, double y);
    void makeFrame(string sx, string sy);
    void makeFrame(double x0, double x1, double y0, double y1, string sx, string sy);
    void makeLegend(double x0, double x1, double y0, double y1, string s);
    void plotHists();
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

  inputs->printLog(TString::Format("x: %.2f, %.2f \ny: %.2f, %.2f", xMinFrame, xMaxFrame, yMinFrame, yMaxFrame).Data(), MyEnums::kDebug);

  if (inputs->logplot) {
    roundToNextPowerOfTen(yMaxFrame);
    if (yMinFrame < 1e-12) {
      int xBin = hists[0]->FindLastBinAbove(0.);
      yMinFrame = hists[0]->GetBinContent(xBin) * 0.5;
      yMaxFrame *= 1.75;
      roundToPrevPowerOfTen(yMinFrame);
    }
  }
  inputs->printLog(TString::Format("x: %.2f, %.2f \ny: %.2f, %.2f", xMinFrame, xMaxFrame, yMinFrame, yMaxFrame).Data(), MyEnums::kDebug);
  makeFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, sx, sy);
}
void Plotter::makeFrame(double x0, double x1, double y0, double y1, string sx, string sy) {
  frame = DrawFrame(x0, x1, y0, y1, sx, sy);
}
void Plotter::makeLegend(double x0, double x1, double y0, double y1, string s) {
  legend = CreateLegend(x0, x1, y0, y1, s.c_str(), textSize);
}

void Plotter::plotHists() {
  if (hists.empty()) {
    string s = "Plotter::plotHists(): Hist vector is empty! Aborting";
    inputs->printLog(s, kErrors);
    return;
  }

  if (inputs->ratioplot) {
    TH1* baseCopy = (TH1*)hists[0]->Clone("baseCopy");
    for (auto& h : hists) h->Divide(baseCopy);
  } // Self-normalise otherwise?

  if (!canvas) makeCanvas();
  if (!frame) {
    string xTitle = "x";
    string yTitle = "y";
    makeFrame(xTitle, yTitle);
  }

  frame->Draw();
  if (legend) legend->Draw("same");
  for (auto o : objects)
    o->Draw("same");

  int iStyle = 0;
  for (auto h : hists) {
    setStyle(h, iStyle++);
    h->Draw(("same" + drawOption).c_str());
    if (inputs->verbosity >= kDebug) h->Print();
  }
  canvas->SaveAs(inputs->outputFileName.c_str());
}

string Plotter::setDrawOption(string s) {
  if (s.at(0) == ' ') { // Single quotes, because checking for a char
    drawOption = s;
  } else {
    drawOption = " " + s;
  }
  return drawOption;
}

// ---------------------------------------------------------------

// Comparison of V0 subtraction with V0 subtraction + track adding

// ---------------------------------------------------------------

struct SubtractionComparison {
  string getFragHistNameFromTrain(int train, int x);
  string getJetHistName(int train, int x);

  void plotJetPt(bool doBase, bool doRatio);
  void plotJetPt();
  void plotV0InJetZ(bool doBase, bool doRatio);
  void plotV0InJetZ();
};

string SubtractionComparison::getFragHistNameFromTrain(int train, int x) {
  string s = "jet-fragmentation";
  if (train == 436232) {
    if (x == MyEnums::kNone || x == kSubtract)
      s += "_id30952";
    else
      s += "_id24580";

    s += "/data/jets/";

    if (x == kSubtract || x == kAddTracks)
      s += "weighted/";

    s += "V0/jetPtV0TrackProjMass";
  }
  return s;
}
string SubtractionComparison::getJetHistName(int train, int x) {
  string s = "jet-fragmentation";
  if (train == 436232) {
    if (x == MyEnums::kNone || x == kSubtract)
      s += "_id30952";
    else
      s += "_id24580";

    s += "/data/jets/";

    if (x == kSubtract || x == kAddTracks)
      s += "weighted/";

    s += "jetPtEtaPhi";
  }
  return s;
}

// Plot jet spectrum
void SubtractionComparison::plotJetPt(bool doBase, bool doRatio) {
  InputSettings x; //x.verbosity = kDebug;
  x.train = 436232;

  x.setInputFileNameFromTrain();
  x.ratioplot = false;
  x.logplot = true;

  Plotter p(x);
  if (!doBase && doRatio) {
    p.setDrawOption("hist");
    p.makeFrame(0., 200., 0.8, 1.2, sPtJet, sRatio);
  }

  // Spectrum - Base vs V0 Sub vs Track-Adding
  x.outputFileName = "subcomp-";
  x.outputFileName += "jetspectrum";
  if (doRatio) {
    x.ratioplot = true;
    x.logplot = false;
    x.outputFileName = "jetratio";
  }
  if (doBase) x.outputFileName += "-all";
  x.outputFileName += ".pdf";
  p.legend = CreateLegend(0.25, 0.45, 0.17, 0.37, "", 0.04);

  if (doRatio) {
    p.addLatex(0.25, 0.8, "This Thesis");
  } else {
    p.addLatex(0.6, 0.8, "This Thesis");
  }

  TFile* file = TFile::Open(x.inputFileName.c_str(), "READ");
  if (!file) {
    string s = "Could not open file " + x.inputFileName;
    x.printLog(s, kErrors);
    return;
  }

  TH3* h3Base = (TH3*)file->Get(getJetHistName(x.train, MyEnums::kNone).c_str());
  TH1* base = (TH1*)h3Base->ProjectionX("base");
  if (doBase) {
    p.hists.push_back(base);
    p.legend->AddEntry(base, "No subtraction");
  }

  TH3* h3Sub = (TH3*)file->Get(getJetHistName(x.train, MyEnums::kSubtract).c_str());
  TH1* sub = (TH1*)h3Sub->ProjectionX("sub");
  p.hists.push_back(sub);
  p.legend->AddEntry(sub, "Subtract V0s");

  TH3* h3TrackAdd = (TH3*)file->Get(getJetHistName(x.train, MyEnums::kAddTracks).c_str());
  TH1* trackAdd = (TH1*)h3TrackAdd->ProjectionX("trackAdd");
  p.hists.push_back(trackAdd);
  p.legend->AddEntry(trackAdd, "Subtract V0s, add tracks");

  if (x.verbosity >= MyEnums::kInfo) {
    // Print integrals. Should have baseInt > subInt
    double baseInt = p.hists[0]->Integral();
    string baseName = p.hists[0]->GetName();
    string s = "Jet pt hist integrals:";

    for (int i = 0; i < p.hists.size(); i++) {
      TH1* h = p.hists[i];
      double integral = h->Integral();
      string name = h->GetName();

      s = TString::Format("%s\n%s: %.3g, diff% .3g, ratio: %.3g", s.c_str(), name.c_str(), integral, baseInt - integral, integral / baseInt).Data();
    }
    x.printLog(s, MyEnums::kInfo);
  }

  p.plotHists();
}
void SubtractionComparison::plotJetPt() {
  cout << "Plotting all jet pt figures" << endl;
  plotJetPt(0, 0);
  plotJetPt(0, 1);
  plotJetPt(1, 0);
  plotJetPt(1, 1);
}

// Plot V0 fragmentation spectrum (currently not per-jet normalised)
void SubtractionComparison::plotV0InJetZ(bool doBase, bool doRatio) {
  InputSettings x;
  x.hadron = "K0S";
  x.train = 436232;

  x.setInputFileNameFromTrain();
  x.ratioplot = false;
  x.logplot = true;
  x.setJetPt(20., 30.);
  const int projectionAxis = 1;

  x.outputFileName = "subcomp-";
  x.outputFileName += "v0z";
  if (doRatio) {
    x.ratioplot = true;
    x.logplot = false;
    x.outputFileName += "ratio";
  }
  if (doBase) x.outputFileName += "-all";
  x.outputFileName += ".pdf";

  if (x.verbosity >= MyEnums::kDebug) {
    cout << "in: " << x.inputFileName << "\nout: " << x.outputFileName
    << "\nBase: " << getFragHistNameFromTrain(x.train, MyEnums::kNone)
    << "\nSub: " << getFragHistNameFromTrain(x.train, MyEnums::kSubtract)
    << "\nAdd: " << getFragHistNameFromTrain(x.train, MyEnums::kAddTracks)
    << endl;
  }

  Plotter p(x);
  p.makeCanvas();
  p.makeLegend(0.25, 0.45, 0.17, 0.37, "");

  if (doRatio && doBase) {
    p.makeFrame(1e-3, 1+1e-3, 0.2, 1.1, sZK0S, sRatio);
    p.addLatex(0.4, 0.70, "This Thesis");
    p.addLatex(0.4, 0.65, "Anti-#it{k}_{T} ch+V0 jets, R = 0.4");
    p.addLatex(0.4, 0.60, TString::Format("%.0f < #it{p}_{T, jet} < %.0f GeV/#it{c}", x.jetptlow, x.jetpthigh).Data());
  } else if (doRatio) {
    p.makeFrame(1e-3, 1+1e-3, 0.5, 1.7, sZK0S, sRatio);
    p.addLatex(0.25, 0.80, "This Thesis");
    p.addLatex(0.25, 0.75, "Anti-#it{k}_{T} ch+V0 jets, R = 0.4");
    p.addLatex(0.25, 0.70, TString::Format("%.0f < #it{p}_{T, jet} < %.0f GeV/#it{c}", x.jetptlow, x.jetpthigh).Data());
  } else {
    p.makeFrame(1e-3, 1+1e-3, 10, 2e4, sZK0S, "\"Counts\"");
    p.addLatex(0.5, 0.80, "This Thesis");
    p.addLatex(0.5, 0.75, "Anti-#it{k}_{T} ch+V0 jets, R = 0.4");
    p.addLatex(0.5, 0.70, TString::Format("%.0f < #it{p}_{T, jet} < %.0f GeV/#it{c}", x.jetptlow, x.jetpthigh).Data());
  }

  TFile* file = TFile::Open(x.inputFileName.c_str(), "READ");
  if (!file) {
    string s = "Could not open file " + x.inputFileName;
    x.printLog(s, kErrors);
    return;
  }

  THnSparse* hnBase = (THnSparse*)file->Get(getFragHistNameFromTrain(x.train, MyEnums::kNone).c_str());
  array<int, 2> jetptbins = getProjectionBins(hnBase->GetAxis(0), x.jetptmin, x.jetptmax);
  x.jetptlow = hnBase->GetAxis(0)->GetBinLowEdge(jetptbins[0]);
  x.jetpthigh = hnBase->GetAxis(0)->GetBinUpEdge(jetptbins[1]);

  hnBase->GetAxis(0)->SetRange(jetptbins[0], jetptbins[1]);
  TH1* base = (TH1*)hnBase->Projection(projectionAxis);
  base->SetName(TString::Format("base_%s_jetpt%.0f-%.0f", x.hadron.c_str(), x.jetptlow, x.jetpthigh).Data());
  if (x.verbosity >= MyEnums::kDebug) base->Print();

  if (doBase) {
    p.hists.push_back(base);
    p.legend->AddEntry(base, "No subtraction");
  }

  THnSparse* hnSub = (THnSparse*)file->Get(getFragHistNameFromTrain(x.train, MyEnums::kSubtract).c_str());
  jetptbins = getProjectionBins(hnSub->GetAxis(0), x.jetptmin, x.jetptmax);
  hnSub->GetAxis(0)->SetRange(jetptbins[0], jetptbins[1]);
  TH1* sub = (TH1*)hnSub->Projection(projectionAxis);
  sub->SetName(TString::Format("sub_%s_jetpt%.0f-%.0f", x.hadron.c_str(), x.jetptlow, x.jetpthigh).Data());
  if (x.verbosity >= MyEnums::kDebug) sub->Print();
  p.hists.push_back(sub);
  p.legend->AddEntry(sub, "Subtract V0s");

  THnSparse* hnAdd = (THnSparse*)file->Get(getFragHistNameFromTrain(x.train, MyEnums::kAddTracks).c_str());
  jetptbins = getProjectionBins(hnAdd->GetAxis(0), x.jetptmin, x.jetptmax);
  hnAdd->GetAxis(0)->SetRange(jetptbins[0], jetptbins[1]);
  TH1* add = (TH1*)hnAdd->Projection(projectionAxis);
  add->SetName(TString::Format("add_%s_jetpt%.0f-%.0f", x.hadron.c_str(), x.jetptlow, x.jetpthigh).Data());
  if (x.verbosity >= MyEnums::kDebug) add->Print();
  p.hists.push_back(add);
  p.legend->AddEntry(add, "Subtract V0s, add tracks");

  if (doRatio && doBase)
    p.addLine(p.hists[0]->GetXaxis()->GetXmin(), 0.5, p.hists[0]->GetXaxis()->GetXmax(), 0.5, 6);

  p.plotHists();
}
void SubtractionComparison::plotV0InJetZ() {
  cout << "Plotting all V0 z figures" << endl;
  plotV0InJetZ(0, 0);
  plotV0InJetZ(0, 1);
  plotV0InJetZ(1, 0);
  plotV0InJetZ(1, 1);
}

// ---------------------------------------------------------------

// Comparison of V0 subtraction vs excluding them beforehand

// ---------------------------------------------------------------

struct JetFinderComparison {
  InputSettings* inputs;
  int trainSub = -1, trainFlags = -1;

  JetFinderComparison(const int tSub, const int tFlags) {
    trainSub = tSub;
    trainFlags = tFlags;
    inputs = new InputSettings();
  }
  JetFinderComparison(const int tSub, const int tFlags, const double jetmin, const double jetmax) {
    trainSub = tSub;
    trainFlags = tFlags;
    inputs = new InputSettings();
    inputs->setPt(jetmin, jetmax);
  }

  string getHistName(Subtraction iSub, HistType iHist, string hadron = "");
  double getNevts(Subtraction iSub);
  double getNjets(double ptmin, double ptmax, Subtraction iSub);
  void plotJetFinderJetPt(bool doBase, bool doRatio);
  void plotJetPt(bool doBase, bool doRatio);
  void plotV0Pt(bool doBase, bool doRatio);
  void plotV0InJetPt(bool doBase, bool doRatio);
  void plotV0InJetZ(bool doBase, bool doRatio);
};

string JetFinderComparison::getHistName(Subtraction iSub, HistType iHist, string hadron = "") {
  string base = "jet-v0qa/tests/";
  string dir, hist;
  bool useDir = true;

  switch (iSub) {
    case MyEnums::kUseFlags:
      dir = "weighted/";
      break;
    case MyEnums::kNone:
      dir = "nosub/";
      break;
    case MyEnums::kSubtract:
      dir = "sub/";
      break;
    default:
      dir = "unknown-sub/";
  }

  switch (iHist) {
    case MyEnums::kEvents:
      useDir = (iSub == MyEnums::kUseFlags);
      hist = "hEvents";
      break;
    case MyEnums::kJet:
      hist = "JetPtEtaPhi";
      break;
    case MyEnums::kV0Pt:
      hist = hadron + "PtEtaPhi";
      break;
    case MyEnums::kV0InJetPt:
      hist = "JetPtEta" + hadron + "Pt";
      break;
    case MyEnums::kV0InJetZ:
      hist = "JetPtEta" + hadron + "Z";
      break;
    default:
      hist = "unknown-hist";
  }

  if (useDir) return base + dir + hist;
  return base + hist;
}

double JetFinderComparison::getNevts(Subtraction iSub) {
  // Final bin contains selected events
  TFile* file = TFile::Open(inputs->inputFileName.c_str());
  TH1* hist = (TH1*)file->Get(getHistName(iSub, MyEnums::kEvents).c_str());
  return hist->GetBinContent(hist->GetNbinsX());
}

double JetFinderComparison::getNjets(double ptmin, double ptmax, Subtraction iSub) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str());
  // TH3* jetPtEtaPhi = (TH3*)file->Get(getHistName(iSub, MyEnums::kJet).c_str());
  THn* jetPtEtaPhi = (THn*)file->Get("jet-finder-v0-data-charged/hJet");
  const int ptAxis = 1, etaAxis = 2;
  array<int, 2> jetPtBins = getProjectionBins(jetPtEtaPhi->GetAxis(ptAxis), ptmin, ptmax);
  array<int, 2> etaBins = getProjectionBins(jetPtEtaPhi->GetAxis(etaAxis), inputs->etaMin, inputs->etaMax);
  TH1* jetPt = (TH1*)jetPtEtaPhi->Projection(ptAxis, "jetPt");

  // if (iSub == MyEnums::kSubtract) {
  //   // Need to correct for the fact that the jet pt spectrum has changed
  //   TH3* h3NoSub = (TH3*)file->Get(getHistName(MyEnums::kNone, MyEnums::kJet).c_str());
  //   etaBins = getProjectionBins(h3NoSub->GetYaxis(), inputs->etaMin, inputs->etaMax);
  //   TH1* hNoSub = (TH1*)h3NoSub->ProjectionX("noSub");

  //   TH3* h3Sub = (TH3*)file->Get(getHistName(MyEnums::kSubtract, MyEnums::kJet).c_str());
  //   etaBins = getProjectionBins(h3Sub->GetYaxis(), inputs->etaMin, inputs->etaMax);
  //   TH1* hSub = (TH1*)h3Sub->ProjectionX("sub");

  //   jetPt->Add(hNoSub, -1);
  //   jetPt->Add(hSub);
  // }

  return jetPt->Integral(jetPtBins[0], jetPtBins[1]);
}

void JetFinderComparison::plotJetFinderJetPt(bool doBase, bool doRatio) {
  InputSettings inputSubtract;
  inputSubtract.train = trainSub;
  InputSettings inputFlags;
  inputFlags.train = trainFlags;

  // Spectrum - Base vs V0 Sub vs Track-Adding
  inputs->outputFileName = "jetfindercomp-jetfinderjetpt";
  if (doRatio) {
    inputs->ratioplot = true;
    inputs->outputFileName += "ratio";
  } else {
    inputs->ratioplot = false;
    inputs->logplot = true;
  }
  if (doBase) inputs->outputFileName += "-all";
  inputs->outputFileName += ".pdf";

  Plotter p(inputs);
  p.makeCanvas();

  if (doRatio) {
    p.setDrawOption("hist");
    p.makeFrame(0., 50., 0.6, 1.15, getPtString(sJet), sRatio);
    p.legend = CreateLegend(0.30, 0.45, 0.20, 0.40, "", 0.04);
    p.addLatex(0.30, 0.83, sThisThesis + ", " + sAliceData + ", " + sSqrtS);
    p.addLatex(0.30, 0.78, sAntikt + " " + sChV0Jets + ", " + sRadius);
  } else {
    p.makeFrame(0., 50., 1e-8, 1e-2, getPtString(sJet), sJetsPerEvent);
    p.legend = CreateLegend(0.30, 0.50, 0.17, 0.37, "", 0.04);
    p.addLatex(0.40, 0.80, sThisThesis + ", " + sAliceData);
    p.addLatex(0.40, 0.75, sSqrtS);
    p.addLatex(0.40, 0.70, sAntikt + " " + sChV0Jets + ", " + sRadius);
  }

  TFile* fileX = TFile::Open(inputSubtract.setInputFileNameFromTrain().c_str(), "READ");
  TFile* fileY = TFile::Open(inputFlags.setInputFileNameFromTrain().c_str(), "READ");

  if (!fileX) {
    string s = "Could not open file " + inputSubtract.inputFileName;
    inputSubtract.printLog(s, kErrors);
    return;
  }
  if (!fileY) {
    string s = "Could not open file " + inputFlags.inputFileName;
    inputFlags.printLog(s, kErrors);
    return;
  }

  // inputs->printLog(TString::Format("From file %s, get %s \nFrom file %s, get %s \nFrom file %s, get %s \nSave to %s",
  //   inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kNone, MyEnums::kJet).c_str(),
  //   inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kSubtract, MyEnums::kJet).c_str(),
  //   inputFlags.inputFileName.c_str(), getHistName(MyEnums::kUseFlags, MyEnums::kJet).c_str(),
  //   inputs->outputFileName.c_str()).Data(),
  //   MyEnums::kDebug);

  inputs->inputFileName = inputSubtract.inputFileName;
  double nEvtsBase = getNevts(MyEnums::kNone);
  THn* h3Base = (THn*)fileX->Get("jet-finder-v0-data-charged/hJet");
  array<int, 2> etaBins = getProjectionBins(h3Base->GetAxis(2), inputs->etaMin, inputs->etaMax);
  h3Base->GetAxis(2)->SetRange(etaBins[0], etaBins[1]);
  TH1* base = (TH1*)h3Base->Projection(1, "base");
  base->Scale(1.0 / nEvtsBase, "width");
  if (doBase) {
    p.hists.push_back(base);
    p.legend->AddEntry(base, "No subtraction");
  }

  double nEvtsSub = getNevts(MyEnums::kSubtract);
  THn* h3Sub = (THn*)fileX->Get("jet-finder-v0-data-charged/hJet");
  h3Sub->GetAxis(2)->SetRange(etaBins[0], etaBins[1]);
  TH1* sub = (TH1*)h3Sub->Projection(1, "sub");
  sub->Scale(1.0 / nEvtsSub, "width");
  p.hists.push_back(sub);
  p.legend->AddEntry(sub, "Subtract V0s");

  inputs->inputFileName = inputFlags.inputFileName;
  double nEvtsFlags = getNevts(MyEnums::kUseFlags);
  THn* h3Flags = (THn*)fileY->Get("jet-finder-v0-data-charged/hJet");
  etaBins = getProjectionBins(h3Flags->GetAxis(2), inputs->etaMin, inputs->etaMax);
  h3Flags->GetAxis(2)->SetRange(etaBins[0], etaBins[1]);
  TH1* flags = (TH1*)h3Flags->Projection(1, "flags");
  flags->Scale(1.0 / nEvtsFlags, "width");
  p.hists.push_back(flags);
  p.legend->AddEntry(flags, "Exclude V0s before clustering");

  if (inputs->verbosity >= MyEnums::kInfo) {
    // Print integrals. Should have baseInt > subInt
    double baseInt = p.hists[0]->Integral();
    string baseName = p.hists[0]->GetName();
    string s = "Jet pt hist integrals:";

    for (int i = 0; i < p.hists.size(); i++) {
      TH1* h = p.hists[i];
      double integral = h->Integral();
      string name = h->GetName();

      s = TString::Format("%s\n%s: %.3g, diff% .3g, ratio: %.3g", s.c_str(), name.c_str(), integral, baseInt - integral, integral / baseInt).Data();
    }
    inputs->printLog(s, MyEnums::kInfo);
  }

  p.plotHists();
}

void JetFinderComparison::plotJetPt(bool doBase, bool doRatio) {
  InputSettings inputSubtract;
  inputSubtract.train = trainSub;
  InputSettings inputFlags;
  inputFlags.train = trainFlags;

  // Spectrum - Base vs V0 Sub vs Track-Adding
  inputs->outputFileName = "jetfindercomp-jetpt";
  if (doRatio) {
    inputs->ratioplot = true;
    inputs->outputFileName += "ratio";
  } else {
    inputs->ratioplot = false;
    inputs->logplot = true;
  }
  if (doBase) inputs->outputFileName += "-all";
  inputs->outputFileName += ".pdf";

  Plotter p(inputs);
  p.makeCanvas();

  if (doRatio) {
    p.setDrawOption("hist");
    p.makeFrame(0., 50., 0.6, 1.15, getPtString(sJet), sRatio);
    p.legend = CreateLegend(0.30, 0.45, 0.20, 0.40, "", 0.04);
    p.addLatex(0.30, 0.83, sThisThesis + ", " + sAliceData + ", " + sSqrtS);
    p.addLatex(0.30, 0.78, sAntikt + " " + sChV0Jets + ", " + sRadius);
  } else {
    p.makeFrame(0., 50., 1e-8, 1e-2, getPtString(sJet), sJetsPerEvent);
    p.legend = CreateLegend(0.30, 0.50, 0.17, 0.37, "", 0.04);
    p.addLatex(0.40, 0.80, sThisThesis + ", " + sAliceData);
    p.addLatex(0.40, 0.75, sSqrtS);
    p.addLatex(0.40, 0.70, sAntikt + " " + sChV0Jets + ", " + sRadius);
  }

  TFile* fileX = TFile::Open(inputSubtract.setInputFileNameFromTrain().c_str(), "READ");
  TFile* fileY = TFile::Open(inputFlags.setInputFileNameFromTrain().c_str(), "READ");

  if (!fileX) {
    string s = "Could not open file " + inputSubtract.inputFileName;
    inputSubtract.printLog(s, kErrors);
    return;
  }
  if (!fileY) {
    string s = "Could not open file " + inputFlags.inputFileName;
    inputFlags.printLog(s, kErrors);
    return;
  }

  inputs->printLog(TString::Format("From file %s, get %s \nFrom file %s, get %s \nFrom file %s, get %s \nSave to %s",
    inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kNone, MyEnums::kJet).c_str(),
    inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kSubtract, MyEnums::kJet).c_str(),
    inputFlags.inputFileName.c_str(), getHistName(MyEnums::kUseFlags, MyEnums::kJet).c_str(),
    inputs->outputFileName.c_str()).Data(),
    MyEnums::kDebug);

  inputs->inputFileName = inputSubtract.inputFileName;
  double nEvtsBase = getNevts(MyEnums::kNone);
  TH3* h3Base = (TH3*)fileX->Get(getHistName(MyEnums::kNone, MyEnums::kJet).c_str());
  TH1* base = (TH1*)h3Base->ProjectionX("base");
  base->Scale(1.0 / nEvtsBase, "width");
  if (doBase) {
    p.hists.push_back(base);
    p.legend->AddEntry(base, "No subtraction");
  }

  double nEvtsSub = getNevts(MyEnums::kSubtract);
  TH3* h3Sub = (TH3*)fileX->Get(getHistName(MyEnums::kSubtract, MyEnums::kJet).c_str());
  TH1* sub = (TH1*)h3Sub->ProjectionX("sub");
  sub->Scale(1.0 / nEvtsSub, "width");
  p.hists.push_back(sub);
  p.legend->AddEntry(sub, "Subtract V0s");

  inputs->inputFileName = inputFlags.inputFileName;
  double nEvtsFlags = getNevts(MyEnums::kUseFlags);
  TH3* h3Flags = (TH3*)fileY->Get(getHistName(MyEnums::kUseFlags, MyEnums::kJet).c_str());
  TH1* flags = (TH1*)h3Flags->ProjectionX("flags");
  flags->Scale(1.0 / nEvtsFlags, "width");
  p.hists.push_back(flags);
  p.legend->AddEntry(flags, "Exclude V0s before clustering");

  if (inputs->verbosity >= MyEnums::kInfo) {
    // Print integrals. Should have baseInt > subInt
    double baseInt = p.hists[0]->Integral();
    string baseName = p.hists[0]->GetName();
    string s = "Jet pt hist integrals:";

    for (int i = 0; i < p.hists.size(); i++) {
      TH1* h = p.hists[i];
      double integral = h->Integral();
      string name = h->GetName();

      s = TString::Format("%s\n%s: %.3g, diff% .3g, ratio: %.3g", s.c_str(), name.c_str(), integral, baseInt - integral, integral / baseInt).Data();
    }
    inputs->printLog(s, MyEnums::kInfo);
  }

  p.plotHists();
}

void JetFinderComparison::plotV0Pt(bool doBase, bool doRatio) {
  InputSettings inputSubtract;
  inputSubtract.train = trainSub;
  InputSettings inputFlags;
  inputFlags.train = trainFlags;

  // Spectrum - Base vs V0 Sub vs Track-Adding
  inputs->outputFileName = "jetfindercomp-v0pt";
  if (doRatio) {
    inputs->ratioplot = true;
    inputs->logplot = false;
    inputs->outputFileName += "ratio";
  } else {
    inputs->ratioplot = false;
    inputs->logplot = false;
  }
  if (doBase) inputs->outputFileName += "-all";
  inputs->outputFileName += ".pdf";

  Plotter p(inputs);
  p.makeCanvas();
  p.makeLegend(0.25, 0.45, 0.17, 0.37, "");

  if (doRatio) {
    p.setDrawOption("hist");
    p.makeFrame(0., 50., 0., 2., sPtK0S, sRatio);

    p.addLatex(0.25, 0.80, sThisThesis + ", " + sAliceData);
    p.addLatex(0.25, 0.75, sSqrtS);
    p.addLatex(0.25, 0.70, sAntikt + " " + sChV0Jets);
    p.addLatex(0.25, 0.65, sRadius + ", " + TString::Format("|%s| < %.1f", sEtaJet.c_str(), inputs->etaMax).Data());
    p.addLatex(0.25, 0.60, getPtJetRangeString(inputs->jetptmin, inputs->jetptmax).c_str());
  } else {
    if (inputs->logplot) {
      p.canvas->SetLogy();
      p.makeFrame(0., 50., 1e-1, 1e3, sPtK0S, sK0SPtPerJet);
    } else {
      p.makeFrame(0., 50., 0., 350., sPtK0S, sK0SPtPerJet);
    }

    p.addLatex(0.45, 0.80, sThisThesis + ", " + sAliceData);
    p.addLatex(0.45, 0.75, sSqrtS);
    p.addLatex(0.45, 0.70, sAntikt + " " + sChV0Jets);
    p.addLatex(0.45, 0.65, sRadius + ", " + TString::Format("|%s| < %.1f", sEtaJet.c_str(), inputs->etaMax).Data());
    p.addLatex(0.45, 0.60, getPtJetRangeString(inputs->jetptmin, inputs->jetptmax).c_str());
  }

  TFile* fileX = TFile::Open(inputSubtract.setInputFileNameFromTrain().c_str(), "READ");
  TFile* fileY = TFile::Open(inputFlags.setInputFileNameFromTrain().c_str(), "READ");

  if (!fileX) {
    string s = "Could not open file " + inputSubtract.inputFileName;
    inputSubtract.printLog(s, kErrors);
    return;
  }
  if (!fileY) {
    string s = "Could not open file " + inputFlags.inputFileName;
    inputFlags.printLog(s, kErrors);
    return;
  }

  inputs->printLog(TString::Format("From file %s, get %s \nFrom file %s, get %s \nFrom file %s, get %s \nSave to %s",
    inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kNone, MyEnums::kV0Pt, inputs->hadron).c_str(),
    inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kSubtract, MyEnums::kV0Pt, inputs->hadron).c_str(),
    inputFlags.inputFileName.c_str(), getHistName(MyEnums::kUseFlags, MyEnums::kV0Pt, inputs->hadron).c_str(),
    inputs->outputFileName.c_str()).Data(),
    MyEnums::kDebug);

  inputs->inputFileName = inputSubtract.inputFileName;
  double nEvtsBase = getNevts(MyEnums::kNone);
  TH3* h3Base = (TH3*)fileX->Get(getHistName(MyEnums::kNone, MyEnums::kV0Pt, inputs->hadron).c_str());
  TH1* base = (TH1*)h3Base->ProjectionX("base");
  base->Scale(1.0 / nEvtsBase, "width");
  if (doBase) {
    p.hists.push_back(base);
    p.legend->AddEntry(base, "No subtraction");
  }

  double nEvtsSub = getNevts(MyEnums::kSubtract);
  TH3* h3Sub = (TH3*)fileX->Get(getHistName(MyEnums::kSubtract, MyEnums::kV0Pt, inputs->hadron).c_str());
  TH1* sub = (TH1*)h3Sub->ProjectionX("sub");
  sub->Scale(1.0 / nEvtsSub, "width");
  p.hists.push_back(sub);
  p.legend->AddEntry(sub, "Subtract V0s");

  inputs->inputFileName = inputFlags.inputFileName;
  double nEvtsFlags = getNevts(MyEnums::kUseFlags);
  TH3* h3Flags = (TH3*)fileY->Get(getHistName(MyEnums::kUseFlags, MyEnums::kV0Pt, inputs->hadron).c_str());
  TH1* flags = (TH1*)h3Flags->ProjectionX("flags");
  flags->Scale(1.0 / nEvtsFlags, "width");
  p.hists.push_back(flags);
  p.legend->AddEntry(flags, "Exclude V0s before clustering");

  if (inputs->verbosity >= MyEnums::kInfo) {
    // Print integrals. Should have baseInt > subInt
    double baseInt = p.hists[0]->Integral();
    string baseName = p.hists[0]->GetName();
    string s = "V0 pt hist integrals:";

    for (int i = 0; i < p.hists.size(); i++) {
      TH1* h = p.hists[i];
      double integral = h->Integral();
      string name = h->GetName();

      s = TString::Format("%s\n%s: %.3g, diff: % .3g, ratio: %.3g", s.c_str(), name.c_str(), integral, baseInt - integral, integral / baseInt).Data();
    }
    s += TString::Format("\nN evts: base: %.0f, sub: %.0f, flags: %.0f", nEvtsBase, nEvtsSub, nEvtsFlags).Data();
    inputs->printLog(s, MyEnums::kInfo);
  }

  p.plotHists();
}

void JetFinderComparison::plotV0InJetPt(bool doBase, bool doRatio) {
  InputSettings inputSubtract;
  inputSubtract.train = trainSub;
  InputSettings inputFlags;
  inputFlags.train = trainFlags;

  // Spectrum - Base vs V0 Sub vs Track-Adding
  inputs->outputFileName = "jetfindercomp-v0pt";
  if (doRatio) {
    inputs->ratioplot = true;
    inputs->logplot = false;
    inputs->outputFileName += "ratio";
  } else {
    inputs->ratioplot = false;
    inputs->logplot = true;
  }
  if (doBase) inputs->outputFileName += "-all";
  inputs->outputFileName += inputs->getSaveNameFromJetPt("", "");
  inputs->outputFileName += ".pdf";

  Plotter p(inputs);
  p.makeCanvas();
  p.makeLegend(0.25, 0.45, 0.17, 0.37, "");

  if (doRatio) {
    p.setDrawOption("hist");
    p.makeFrame(0., inputs->jetptmax, 0., 1.2, sPtK0S, sRatio);

    p.addLatex(0.25, 0.80, sThisThesis + ", " + sAliceData);
    p.addLatex(0.25, 0.75, sSqrtS);
    p.addLatex(0.25, 0.70, sAntikt + " " + sChV0Jets);
    p.addLatex(0.25, 0.65, sRadius + ", " + TString::Format("|%s| < %.1f", sEtaJet.c_str(), inputs->etaMax).Data());
    p.addLatex(0.25, 0.60, getPtJetRangeString(inputs->jetptmin, inputs->jetptmax).c_str());
  } else {
    if (inputs->logplot) {
      p.canvas->SetLogy();
      p.makeFrame(0., inputs->jetptmax, 1e-6, 1e1, sPtK0S, sK0SPtPerJet);
    } else {
      p.makeFrame(0., inputs->jetptmax, 0., 50., sPtK0S, sK0SPtPerJet);
    }

    p.addLatex(0.45, 0.80, sThisThesis + ", " + sAliceData);
    p.addLatex(0.45, 0.75, sSqrtS);
    p.addLatex(0.45, 0.70, sAntikt + " " + sChV0Jets);
    p.addLatex(0.45, 0.65, sRadius + ", " + TString::Format("|%s| < %.1f", sEtaJet.c_str(), inputs->etaMax).Data());
    p.addLatex(0.45, 0.60, getPtJetRangeString(inputs->jetptmin, inputs->jetptmax).c_str());
  }

  TFile* fileX = TFile::Open(inputSubtract.setInputFileNameFromTrain().c_str(), "READ");
  TFile* fileY = TFile::Open(inputFlags.setInputFileNameFromTrain().c_str(), "READ");

  if (!fileX) {
    string s = "Could not open file " + inputSubtract.inputFileName;
    inputSubtract.printLog(s, kErrors);
    return;
  }
  if (!fileY) {
    string s = "Could not open file " + inputFlags.inputFileName;
    inputFlags.printLog(s, kErrors);
    return;
  }

  inputs->printLog(TString::Format("From file %s, get %s \nFrom file %s, get %s \nFrom file %s, get %s \nSave to %s",
    inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kNone, MyEnums::kV0InJetPt, inputs->hadron).c_str(),
    inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kSubtract, MyEnums::kV0InJetPt, inputs->hadron).c_str(),
    inputFlags.inputFileName.c_str(), getHistName(MyEnums::kUseFlags, MyEnums::kV0InJetPt, inputs->hadron).c_str(),
    inputs->outputFileName.c_str()).Data(),
    MyEnums::kDebug);

  inputs->inputFileName = inputSubtract.inputFileName;
  double nJetsBase = getNjets(inputs->jetptmin, inputs->jetptmax, MyEnums::kNone);
  double nEvtsBase = getNevts(MyEnums::kNone);
  TH3* h3Base = (TH3*)fileX->Get(getHistName(MyEnums::kNone, MyEnums::kV0InJetPt, inputs->hadron).c_str());
  array<int, 2> jetBins = getProjectionBins(h3Base->GetXaxis(), inputs->jetptmin, inputs->jetptmax);
  TH1* base = (TH1*)h3Base->ProjectionZ("base", jetBins[0], jetBins[1], 1, h3Base->GetNbinsY());
  if (inputs->normalisePerJet)
    base->Scale(1.0 / nJetsBase, "width");
  else
    base->Scale(1.0 / nEvtsBase, "width");

  if (doBase) {
    p.hists.push_back(base);
    p.legend->AddEntry(base, "No subtraction");
  }

  double nJetsSub = getNjets(inputs->jetptmin, inputs->jetptmax, MyEnums::kSubtract);
  double nEvtsSub = getNevts(MyEnums::kSubtract);
  TH3* h3Sub = (TH3*)fileX->Get(getHistName(MyEnums::kSubtract, MyEnums::kV0InJetPt, inputs->hadron).c_str());
  jetBins = getProjectionBins(h3Sub->GetXaxis(), inputs->jetptmin, inputs->jetptmax);
  TH1* sub = (TH1*)h3Sub->ProjectionZ("sub", jetBins[0], jetBins[1], 1, h3Sub->GetNbinsY());
  if (inputs->normalisePerJet)
    sub->Scale(1.0 / nJetsSub, "width");
  else
    sub->Scale(1.0 / nEvtsSub, "width");
  p.hists.push_back(sub);
  p.legend->AddEntry(sub, "Subtract V0s");

  inputs->inputFileName = inputFlags.inputFileName;
  double nJetsFlags = getNjets(inputs->jetptmin, inputs->jetptmax, MyEnums::kUseFlags);
  double nEvtsFlags = getNevts(MyEnums::kUseFlags);
  TH3* h3Flags = (TH3*)fileY->Get(getHistName(MyEnums::kUseFlags, MyEnums::kV0InJetPt, inputs->hadron).c_str());
  jetBins = getProjectionBins(h3Flags->GetXaxis(), inputs->jetptmin, inputs->jetptmax);
  TH1* flags = (TH1*)h3Flags->ProjectionZ("flags", jetBins[0], jetBins[1], 1, h3Flags->GetNbinsY());
  if (inputs->normalisePerJet)
    flags->Scale(1.0 / nJetsFlags, "width");
  else
    flags->Scale(1.0 / nEvtsFlags, "width");

  p.hists.push_back(flags);
  p.legend->AddEntry(flags, "Exclude V0s before clustering");

  if (inputs->verbosity >= MyEnums::kInfo) {
    // Print integrals. Should have baseInt > subInt
    double baseInt = p.hists[0]->Integral();
    string baseName = p.hists[0]->GetName();
    string s = "V0 z hist integrals:";

    for (int i = 0; i < p.hists.size(); i++) {
      TH1* h = p.hists[i];
      double integral = h->Integral();
      string name = h->GetName();

      s = TString::Format("%s\n%s: %.3g, diff: % .3g, ratio: %.3g", s.c_str(), name.c_str(), integral, baseInt - integral, integral / baseInt).Data();
    }
    s += TString::Format("\nN jets (%.0f < pT < %.0f): base: %.0f, sub: %.0f, flags: %.0f", inputs->jetptmin, inputs->jetpthigh, nJetsBase, nJetsSub, nJetsFlags).Data();
    inputs->printLog(s, MyEnums::kInfo);
  }

  p.plotHists();
 }

void JetFinderComparison::plotV0InJetZ(bool doBase, bool doRatio) {
  InputSettings inputSubtract;
  inputSubtract.train = trainSub;
  InputSettings inputFlags;
  inputFlags.train = trainFlags;

  // Spectrum - Base vs V0 Sub vs Track-Adding
  inputs->outputFileName = "jetfindercomp-v0z";
  if (doRatio) {
    inputs->ratioplot = true;
    inputs->logplot = false;
    inputs->outputFileName += "ratio";
  } else {
    inputs->ratioplot = false;
    inputs->logplot = false;
  }
  if (doBase) inputs->outputFileName += "-all";
  inputs->outputFileName += inputs->getSaveNameFromJetPt("", "");
  inputs->outputFileName += ".pdf";

  Plotter p(inputs);
  p.makeCanvas();
  p.makeLegend(0.25, 0.45, 0.17, 0.37, "");

  if (doRatio) {
    p.setDrawOption("hist");
    p.makeFrame(0., 1., 0., 2., sZK0S, sRatio);

    p.addLatex(0.25, 0.80, sThisThesis + ", " + sAliceData);
    p.addLatex(0.25, 0.75, sSqrtS);
    p.addLatex(0.25, 0.70, sAntikt + " " + sChV0Jets);
    p.addLatex(0.25, 0.65, sRadius + ", " + TString::Format("|%s| < %.1f", sEtaJet.c_str(), inputs->etaMax).Data());
    p.addLatex(0.25, 0.60, getPtJetRangeString(inputs->jetptmin, inputs->jetptmax).c_str());
  } else {
    if (inputs->logplot) {
      p.canvas->SetLogy();
      p.makeFrame(0., 1., 1e-1, 1e3, sZK0S, sK0SZPerJet);
    } else {
      p.makeFrame(0., 1., 0., 350., sZK0S, sK0SZPerJet);
    }

    p.addLatex(0.45, 0.80, sThisThesis + ", " + sAliceData);
    p.addLatex(0.45, 0.75, sSqrtS);
    p.addLatex(0.45, 0.70, sAntikt + " " + sChV0Jets);
    p.addLatex(0.45, 0.65, sRadius + ", " + TString::Format("|%s| < %.1f", sEtaJet.c_str(), inputs->etaMax).Data());
    p.addLatex(0.45, 0.60, getPtJetRangeString(inputs->jetptmin, inputs->jetptmax).c_str());
  }

  TFile* fileX = TFile::Open(inputSubtract.setInputFileNameFromTrain().c_str(), "READ");
  TFile* fileY = TFile::Open(inputFlags.setInputFileNameFromTrain().c_str(), "READ");

  if (!fileX) {
    string s = "Could not open file " + inputSubtract.inputFileName;
    inputSubtract.printLog(s, kErrors);
    return;
  }
  if (!fileY) {
    string s = "Could not open file " + inputFlags.inputFileName;
    inputFlags.printLog(s, kErrors);
    return;
  }

  inputs->printLog(TString::Format("From file %s, get %s \nFrom file %s, get %s \nFrom file %s, get %s \nSave to %s",
    inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kNone, MyEnums::kV0InJetZ, inputs->hadron).c_str(),
    inputSubtract.inputFileName.c_str(), getHistName(MyEnums::kSubtract, MyEnums::kV0InJetZ, inputs->hadron).c_str(),
    inputFlags.inputFileName.c_str(), getHistName(MyEnums::kUseFlags, MyEnums::kV0InJetZ, inputs->hadron).c_str(),
    inputs->outputFileName.c_str()).Data(),
    MyEnums::kDebug);

  inputs->inputFileName = inputSubtract.inputFileName;
  double nJetsBase = getNjets(inputs->jetptmin, inputs->jetptmax, MyEnums::kNone);
  double nEvtsBase = getNevts(MyEnums::kNone);
  TH3* h3Base = (TH3*)fileX->Get(getHistName(MyEnums::kNone, MyEnums::kV0InJetZ, inputs->hadron).c_str());
  TH1* base = (TH1*)h3Base->ProjectionZ("base");
  if (inputs->normalisePerJet)
    base->Scale(1.0 / nJetsBase, "width");
  else
    base->Scale(1.0 / nEvtsBase, "width");

  if (doBase) {
    p.hists.push_back(base);
    p.legend->AddEntry(base, "No subtraction");
  }

  double nJetsSub = getNjets(inputs->jetptmin, inputs->jetptmax, MyEnums::kSubtract);
  double nEvtsSub = getNevts(MyEnums::kSubtract);
  TH3* h3Sub = (TH3*)fileX->Get(getHistName(MyEnums::kSubtract, MyEnums::kV0InJetZ, inputs->hadron).c_str());
  TH1* sub = (TH1*)h3Sub->ProjectionZ("sub");
  if (inputs->normalisePerJet)
    sub->Scale(1.0 / nJetsSub, "width");
  else
    sub->Scale(1.0 / nEvtsSub, "width");

  p.hists.push_back(sub);
  p.legend->AddEntry(sub, "Subtract V0s");

  inputs->inputFileName = inputFlags.inputFileName;
  double nJetsFlags = getNjets(inputs->jetptmin, inputs->jetptmax, MyEnums::kUseFlags);
  double nEvtsFlags = getNevts(MyEnums::kUseFlags);
  TH3* h3Flags = (TH3*)fileY->Get(getHistName(MyEnums::kUseFlags, MyEnums::kV0InJetZ, inputs->hadron).c_str());
  TH1* flags = (TH1*)h3Flags->ProjectionZ("flags");
  if (inputs->normalisePerJet)
    flags->Scale(1.0 / nJetsFlags, "width");
  else
    flags->Scale(1.0 / nEvtsFlags, "width");

  p.hists.push_back(flags);
  p.legend->AddEntry(flags, "Exclude V0s before clustering");

  if (inputs->verbosity >= MyEnums::kInfo) {
    // Print integrals. Should have baseInt > subInt
    double baseInt = p.hists[0]->Integral();
    string baseName = p.hists[0]->GetName();
    string s = "V0 z hist integrals:";

    for (int i = 0; i < p.hists.size(); i++) {
      TH1* h = p.hists[i];
      double integral = h->Integral();
      string name = h->GetName();

      s = TString::Format("%s\n%s: %.3g, diff: % .3g, ratio: %.3g", s.c_str(), name.c_str(), integral, baseInt - integral, integral / baseInt).Data();
    }
    s += TString::Format("\nN jets (%.0f < pT < %.0f): base: %.0f, sub: %.0f, flags: %.0f", inputs->jetptmin, inputs->jetpthigh, nJetsBase, nJetsSub, nJetsFlags).Data();
    inputs->printLog(s, MyEnums::kInfo);
  }

  p.plotHists();
}

// ---------------------------------------------------------------

void plotjetpt() {
  int trainSub = 486343;
  int trainFlags = 486347;
  JetFinderComparison jfc(trainSub, trainFlags);
  jfc.inputs->verbosity = MyEnums::kDebug;
  jfc.inputs->logplot = true;
  jfc.inputs->setEta(-0.5, 0.5);
  jfc.plotJetPt(1, 0);
  jfc.inputs->logplot = false;
  jfc.plotJetPt(1, 1);
}
void plotv0pt() {
  int trainSub = 486343;
  int trainFlags = 486347;
  JetFinderComparison jfc(trainSub, trainFlags);
  jfc.inputs->hadron = "K0S";
  jfc.inputs->verbosity = MyEnums::kDebug;
  jfc.inputs->setEta(-0.9, 0.9);
  jfc.plotV0Pt(1, 0);
  jfc.plotV0Pt(1, 1);
}
void plotv0injet() {
  gROOT->SetBatch();
  int trainSub = 486343;
  int trainFlags = 486347;
  JetFinderComparison jfc(trainSub, trainFlags);
  jfc.inputs->hadron = "K0S";
  jfc.inputs->verbosity = MyEnums::kDebug;
  jfc.inputs->setEta(-0.5, 0.5);

  jfc.inputs->setJetPt(10., 20.);
  jfc.plotV0InJetZ(1, 0);
  jfc.plotV0InJetZ(1, 1);
  jfc.plotV0InJetPt(1, 0);
  jfc.plotV0InJetPt(1, 1);
  jfc.inputs->setJetPt(20., 30.);
  jfc.plotV0InJetZ(1, 0);
  jfc.plotV0InJetZ(1, 1);
  jfc.plotV0InJetPt(1, 0);
  jfc.plotV0InJetPt(1, 1);
}

JetFinderComparison setupJFC(double jetptmin = 10., double jetptmax = 20., bool normalisePerJet = true) {
  int trainSub = 486343;
  int trainFlags = 486347;
  JetFinderComparison jfc(trainSub, trainFlags, jetptmin, jetptmax);
  jfc.inputs->hadron = "K0S";
  jfc.inputs->verbosity = MyEnums::kDebug;
  jfc.inputs->normalisePerJet = normalisePerJet;
  return jfc;
}

// Compare all jets (from jet finder thn) with ch+v0 jets
void plotJetPt_Ch_ChV0() {
  InputSettings x; x.verbosity = kDebug;
  x.train = 476306;
  x.setInputFileNameFromTrain();
  // TODO: Finish
}

#endif // SUBTRACTION_COMPARISON_C
