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

#ifndef __PLOTTRIGGERCORRELATION_H__
#define __PLOTTRIGGERCORRELATION_H__

namespace MyStrings {
  const string sALICE      = "ALICE";
  const string sAntiktJets = "Anti-#it{k}_{T} jets";
  const string sCounts     = "Counts";
  const string sGevC       = "GeV/#it{c}";
  const string sGevCC      = "GeV/#it{c}^{2}";
  const string sMass       = "#it{M}";
  const string sNumber     = "#it{N}";
  const string sRadius     = "#it{R} = 0.4";
  const string sRatio      = "Ratio";
  const string sSigma      = "#sigma";
  const string sSqrtS      = "#sqrt{s} = 13.6 TeV";
  const string sThisThesis = "This Thesis";
  const string sV0         = "V0";
  const string sZ          = "#it{z}";

  const string sData       = sALICE + " pp data";
  const string sPythia     = sALICE + " PYTHIA";

  string getPtString(string subscript) {
    if (subscript.empty())
      return TString::Format("#it{p}_{T}").Data();
    else
      return TString::Format("#it{p}_{T, %s}", subscript.c_str()).Data();
  }

  // Triggers
  const string sEMCalReadout = "EMCalReadout";
  const string sJetFullHighPt = "JetFullHighPt";
  const string sJetFullLowPt = "JetFullLowPt";
  const string sGammaHighPtEMCal = "GammaHighPtEMCal";
  const string sGammaLowPtEMCal = "GammaLowPtEMCal";
}

namespace MyEnums {
  enum Verb {kErrors, kWarnings, kInfo, kDebug, kDebugMax};
  enum TriggerBins {kEMCalReadout = 9, kJetFullHighPt = 10, kJetFullLowPt = 11, kGammaHighPtEMCal = 16, kGammaLowPtEMCal = 18, kGammaVeryLowPtEMCal = 20};
}

string getTriggerName(int trigger) {
  switch (trigger) {
    case MyEnums::kEMCalReadout: return MyStrings::sEMCalReadout;
    case MyEnums::kJetFullHighPt: return MyStrings::sJetFullHighPt;
    case MyEnums::kJetFullLowPt: return MyStrings::sJetFullLowPt;
    case MyEnums::kGammaHighPtEMCal: return MyStrings::sGammaHighPtEMCal;
    case MyEnums::kGammaLowPtEMCal: return MyStrings::sGammaLowPtEMCal;
    default: return "Unknown Trigger";
  }
}

using namespace MyStrings;
using namespace MyEnums;

// -------------------------------------------------------------------------------------------------
//
// Struct for input settings
//
// -------------------------------------------------------------------------------------------------

struct InputSettings{
  private:
  public:
    const int conesPerJet = 2;
    int train = 0;
    int rebinNumber = -1;
    string hadron = "";
    string histName = "trigger-correlations/triggerCorrelations";
    string inputFileName = "";
    string outputFileName = "";
    double ptmin = -1e3, ptmax = -1e3, lowpt = -1e3, highpt = -1e3;
    double ptminjet = -1e3, ptmaxjet = -1e3, lowptjet = -1e3, highptjet = -1e3;
    bool logplot = false;
    bool ratioplot = false;
    vector<vector<double>> ptBinEdges = {};

    Verb verbosity = kWarnings;

    double massWindowMin = -1., massWindowMax = -1.;
    double polInitx0 = -1., polInitx1 = -1., polInitx2 = -1.;
    double signalRegionMin = -1., signalRegionMax = -1.;
    double nSigma = -1., nSigmaL = -1., nSigmaR = -1.;

    string getHistName(int x);
    double getMass();
    string getNameFromJetPt(string prefix, string suffix);
    string getNameFromPt(string prefix, string suffix);
    string printLog(string message, int verbThreshold);
    int setHadron(string h);
    string setInputFileNameFromTrain();
    void setLowHighFromAxis(const TAxis* axis, double& low, double& high);
    void setJetPt(double a, double b);
    void setPt(double a, double b);
    vector<vector<double>> setPtBinEdgesFromHadron();
    vector<vector<double>> setPtBinEdgesSorted(vector<vector<double>> x);
    template <typename T> int writeOutputToFile(T* obj);
};

string InputSettings::getNameFromJetPt(string prefix, string suffix = "") {
  string s = TString::Format("%s_jetpt%.f-%.f%s", prefix.c_str(), lowptjet, highptjet, suffix.c_str()).Data();
  return s;
}

string InputSettings::getNameFromPt(string prefix, string suffix = "") {
  string s = TString::Format("%s_pt%.1f-%.1f%s", prefix.c_str(), lowpt, highpt, suffix.c_str()).Data();
  return s;
}

string InputSettings::printLog(string message, int verbThreshold) {
  if (this->verbosity < verbThreshold)
    return "";

  cout << message << endl;
  return message;
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(this->train) + "/AnalysisResults.root";
  this->inputFileName = s;
  return s;
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
    const string sCounts = "Counts";
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
    void makeCanvas(double x, double y, string s);
    void makeFrame(string sx, string sy);
    void makeFrame(double x0, double x1, double y0, double y1, string sx, string sy);
    void makeLegend(double x0, double x1, double y0, double y1, string s);
    void plot();
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

  inputs->printLog(TString::Format("Plotter::makeFrame() x: %.2f, %.2f \ny: %.2f, %.2f", xMinFrame, xMaxFrame, yMinFrame, yMaxFrame).Data(), MyEnums::kDebug);

  if (inputs->logplot) {
    roundToNextPowerOfTen(yMaxFrame);
    if (yMinFrame < 1e-12) {
      int xBin = hists[0]->FindLastBinAbove(0.);
      yMinFrame = hists[0]->GetBinContent(xBin) * 0.5;
      yMaxFrame *= 1.75;
      roundToPrevPowerOfTen(yMinFrame);
    }
  }
  inputs->printLog(TString::Format("Plotter::makeFrame() x: %.2f, %.2f \ny: %.2f, %.2f", xMinFrame, xMaxFrame, yMinFrame, yMaxFrame).Data(), MyEnums::kDebug);
  makeFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, sx, sy);
}
void Plotter::makeCanvas(double x, double y, string s = "canvas") {
  makeCanvas(s, x, y);
}
void Plotter::makeFrame(double x0, double x1, double y0, double y1, string sx, string sy) {
  if (!canvas) makeCanvas();
  frame = DrawFrame(x0, x1, y0, y1, sx, sy);
}
void Plotter::makeLegend(double x0, double x1, double y0, double y1, string s) {
  legend = CreateLegend(x0, x1, y0, y1, s.c_str(), textSize);
}

void Plotter::plot() {
  if (inputs->ratioplot) {
    TH1* baseCopy = (TH1*)hists[0]->Clone("baseCopy");
    for (auto& h : hists) h->Divide(baseCopy);
  }

  if (!canvas) makeCanvas();
  if (!frame) makeFrame("x", "y");

  frame->Draw();
  if (legend) legend->Draw("same");
  for (auto o : objects)
    o->Draw("same");

  for (auto h : hists) {
    h->Draw(("same" + drawOption).c_str());
    if (inputs->verbosity >= kDebug) h->Print();
  }
  canvas->SaveAs(inputs->outputFileName.c_str());
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

cout << "Trigger options: "
     << "k" << sEMCalReadout << ", "
     << "k" << sJetFullHighPt << ", "
     << "k" << sJetFullLowPt << ", "
     << "k" << sGammaHighPtEMCal << ", "
     << "k" << sGammaLowPtEMCal << endl;

void plotTriggerCorrelationAll() {
  InputSettings x; x.verbosity = kDebug;
  x.hadron = "K0S";
  x.train = 461466;

  enum NewTriggerBins {kFJHigh = 1, kFJLow, kGHigh, kGLow, kFJHighAndGHigh, kFJLowAndGHigh, kFJHighAndGLow, kFJLowAndGLow};

  x.setInputFileNameFromTrain();
  x.outputFileName = "triggerCorrelation.pdf";

  TFile* file = TFile::Open(x.inputFileName.c_str(), "READ");
  TH2* hTriggerCorrelation = (TH2*)file->Get(x.histName.c_str());

  int nTriggers = kFJLowAndGLow;
  TH1* trigCor = new TH1D("trigCor", "Trigger Correlations", nTriggers, 0, nTriggers);

  trigCor->GetXaxis()->SetBinLabel(kFJHigh, "JetFullHighPt");
  trigCor->GetXaxis()->SetBinLabel(kFJLow, "JetFullLowPt");
  trigCor->GetXaxis()->SetBinLabel(kGHigh, "GammaHighPtEMCal");
  trigCor->GetXaxis()->SetBinLabel(kGLow, "GammaLowPtEMCal");
  trigCor->GetXaxis()->SetBinLabel(kFJHighAndGHigh, "JetFullHighPtAndGammaHighPtEMCal");
  trigCor->GetXaxis()->SetBinLabel(kFJLowAndGHigh, "JetFullLowPtAndGammaHighPtEMCal");
  trigCor->GetXaxis()->SetBinLabel(kFJHighAndGLow, "JetFullHighPtAndGammaLowPtEMCal");
  trigCor->GetXaxis()->SetBinLabel(kFJLowAndGLow, "JetFullLowPtAndGammaLowPtEMCal");

  trigCor->SetBinContent(kFJHigh, hTriggerCorrelation->GetBinContent(kJetFullHighPt, kJetFullHighPt));
  trigCor->SetBinContent(kFJLow, hTriggerCorrelation->GetBinContent(kJetFullLowPt, kJetFullLowPt));
  trigCor->SetBinContent(kGHigh, hTriggerCorrelation->GetBinContent(kGammaHighPtEMCal, kGammaHighPtEMCal));
  trigCor->SetBinContent(kGLow, hTriggerCorrelation->GetBinContent(kGammaLowPtEMCal, kGammaLowPtEMCal));
  trigCor->SetBinContent(kFJHighAndGHigh, hTriggerCorrelation->GetBinContent(kJetFullHighPt, kGammaHighPtEMCal));
  trigCor->SetBinContent(kFJLowAndGHigh, hTriggerCorrelation->GetBinContent(kJetFullLowPt, kGammaHighPtEMCal));
  trigCor->SetBinContent(kFJHighAndGLow, hTriggerCorrelation->GetBinContent(kJetFullHighPt, kGammaLowPtEMCal));
  trigCor->SetBinContent(kFJLowAndGLow, hTriggerCorrelation->GetBinContent(kJetFullLowPt, kGammaLowPtEMCal));

  trigCor->Draw();
}

TH1* getTriggerCorrelation(InputSettings& x, TriggerBins trigger, TriggerBins associate) {
  if (trigger == associate) {
    x.printLog("getTriggerCorrelation() Error: trigger and associate are the same", kErrors);
    return nullptr;
  }

  TFile* file = TFile::Open(x.inputFileName.c_str(), "READ");
  TH2* hTriggerCorrelation = (TH2*)file->Get(x.histName.c_str());

  enum Bins {kTrigger = 1, kAssociate, kBothBin};

  int nTriggersPerPlot = kBothBin;
  TH1* triggers = new TH1D(("triggerCorrelation_" + getTriggerName(trigger) + "-" + getTriggerName(associate)).c_str(), "", nTriggersPerPlot, 0, nTriggersPerPlot);

  triggers->GetXaxis()->SetBinLabel(kTrigger, getTriggerName(trigger).c_str());
  triggers->GetXaxis()->SetBinLabel(kAssociate, getTriggerName(associate).c_str());
  triggers->GetXaxis()->SetBinLabel(kBothBin, "Overlap");

  triggers->SetBinContent(kTrigger, hTriggerCorrelation->GetBinContent(trigger, trigger));
  triggers->SetBinContent(kAssociate, hTriggerCorrelation->GetBinContent(associate, associate));
  triggers->SetBinContent(kBothBin, hTriggerCorrelation->GetBinContent(trigger, associate));
  return triggers;
}

void plotTriggerCorrelation(TriggerBins trigger, TriggerBins associate) {
  InputSettings x; x.verbosity = kDebug;
  x.train = 462605;
  x.setInputFileNameFromTrain();
  x.outputFileName = "triggerCorrelation_" + getTriggerName(trigger) + "-" + getTriggerName(associate) + ".pdf";
  TH1* triggers = getTriggerCorrelation(x, trigger, associate);
  if (!triggers) {
    x.printLog(TString::Format("plotTriggerCorrelation() Error: could not get trigger correlation for %d and %d", trigger, associate).Data(), kErrors);
    return;
  }

  setStyle(triggers, 0);
  triggers->SetStats(0);
  triggers->Scale(1. / triggers->GetBinContent(1)); // Normalise to primary trigger
  triggers->SetMinimum(0.);
  triggers->SetMaximum(1.1);

  Plotter p(x);
  p.frame = (TH1F*)triggers->Clone("frame");
  p.hists.push_back(triggers);

  p.addLatex(0.45, 0.83, TString::Format("%s, %s", sThisThesis.c_str(), sData.c_str()).Data());
  p.addLatex(0.45, 0.78, sSqrtS.c_str());
  p.plot();
}

void plotTriggerCorrelationText(TriggerBins trigger, TriggerBins associate) {
  InputSettings x; x.verbosity = kDebug;
  x.train = 462605;
  x.setInputFileNameFromTrain();
  x.outputFileName = "triggerCorrelation_" + getTriggerName(trigger) + "-" + getTriggerName(associate) + ".pdf";
  TH1* triggers = getTriggerCorrelation(x, trigger, associate);
  if (!triggers) {
    x.printLog(TString::Format("plotTriggerCorrelation() Error: could not get trigger correlation for %d and %d", trigger, associate).Data(), kErrors);
    return;
  }

  setStyle(triggers, 0);
  triggers->SetStats(0);
  triggers->Scale(1. / triggers->GetBinContent(1)); // Normalise to primary trigger
  triggers->SetMinimum(0.);
  triggers->SetMaximum(1.1);

  Plotter p(x);
  p.makeCanvas();
  p.frame = (TH1F*)triggers->Clone("frame");

  p.addLatex(0.45, 0.83, TString::Format("%s, %s", sThisThesis.c_str(), sData.c_str()).Data());
  p.addLatex(0.45, 0.78, sSqrtS.c_str());

  // p.addLatex(0.45, triggers->GetBinContent(2), TString::Format("%.2f", triggers->GetBinContent(2)).Data());
  // p.addLatex(0.45, triggers->GetBinContent(3), TString::Format("%.2f", triggers->GetBinContent(3)).Data());

  p.frame->Draw("hist text0");
  for (auto l : p.objects) l->Draw();
  p.canvas->SaveAs(x.outputFileName.c_str());
}

void plotTriggerCorrelation() {
  gROOT->SetBatch(true);
  plotTriggerCorrelation(kGammaHighPtEMCal, kJetFullHighPt);
  plotTriggerCorrelation(kGammaHighPtEMCal, kJetFullLowPt);
  plotTriggerCorrelation(kGammaLowPtEMCal, kJetFullHighPt);
  plotTriggerCorrelation(kGammaLowPtEMCal, kJetFullLowPt);

  plotTriggerCorrelation(kJetFullHighPt, kJetFullLowPt);
  plotTriggerCorrelation(kGammaHighPtEMCal, kGammaLowPtEMCal);
}

void jetGammaCorrelation1D(bool drawText = false) {
  vector<TriggerBins> triggers = {kJetFullLowPt, kJetFullHighPt, kGammaHighPtEMCal, kGammaLowPtEMCal};
  enum Bins {kJetHighGammaHigh = 1, kJetLowGammaHigh, kJetHighGammaLow, kJetLowGammaLow, kNbinsPlusOne};
  InputSettings x; x.verbosity = kDebug;
  x.train = 462605;
  x.setInputFileNameFromTrain();
  double textSize = 0.05;

  string histName = "triggerCorrelation_";
  for (auto trig : triggers) {
    histName += getTriggerName(trig) + "-";
  }
  x.outputFileName = histName + "1D.pdf";

  TFile* file = TFile::Open(x.inputFileName.c_str(), "READ");
  TH2* hTriggerCorrelation = (TH2*)file->Get(x.histName.c_str());

  int nTriggersPerPlot = kNbinsPlusOne - 1;
  TH1* hTriggers = new TH1D(histName.c_str(), ";;Fraction of jet trigger rate", nTriggersPerPlot, 0, nTriggersPerPlot);
  hTriggers->GetXaxis()->SetBinLabel(kJetHighGammaHigh, TString::Format("#splitline{%s &}{ %s}", getTriggerName(kJetFullHighPt).c_str(), getTriggerName(kGammaHighPtEMCal).c_str()).Data());
  hTriggers->GetXaxis()->SetBinLabel(kJetLowGammaHigh, TString::Format("#splitline{%s &}{ %s}", getTriggerName(kJetFullLowPt).c_str(), getTriggerName(kGammaHighPtEMCal).c_str()).Data());
  hTriggers->GetXaxis()->SetBinLabel(kJetHighGammaLow, TString::Format("#splitline{%s &}{ %s}", getTriggerName(kJetFullHighPt).c_str(), getTriggerName(kGammaLowPtEMCal).c_str()).Data());
  hTriggers->GetXaxis()->SetBinLabel(kJetLowGammaLow, TString::Format("#splitline{%s &}{ %s}", getTriggerName(kJetFullLowPt).c_str(), getTriggerName(kGammaLowPtEMCal).c_str()).Data());

  double jetLowCorrection = 10, gammaLowCorrection = 20.;
  double nJetFullHigh          = hTriggerCorrelation->GetBinContent(kJetFullHighPt, kJetFullHighPt);
  double nJetFullLow           = hTriggerCorrelation->GetBinContent(kJetFullLowPt, kJetFullLowPt) * jetLowCorrection;
  double nGammaHigh            = hTriggerCorrelation->GetBinContent(kGammaHighPtEMCal, kGammaHighPtEMCal);
  double nGammaLow             = hTriggerCorrelation->GetBinContent(kGammaLowPtEMCal, kGammaLowPtEMCal) * gammaLowCorrection;
  double nJetFullHighGammaHigh = hTriggerCorrelation->GetBinContent(kJetFullHighPt, kGammaHighPtEMCal);
  double nJetFullHighGammaLow  = hTriggerCorrelation->GetBinContent(kJetFullHighPt, kGammaLowPtEMCal) * gammaLowCorrection;
  double nJetFullLowGammaHigh  = hTriggerCorrelation->GetBinContent(kJetFullLowPt, kGammaHighPtEMCal) * jetLowCorrection;
  double nJetFullLowGammaLow   = hTriggerCorrelation->GetBinContent(kJetFullLowPt, kGammaLowPtEMCal) * jetLowCorrection * gammaLowCorrection;

  hTriggers->SetBinContent(kJetHighGammaHigh, nJetFullHighGammaHigh / nJetFullHigh);
  hTriggers->SetBinContent(kJetHighGammaLow, nJetFullHighGammaLow / nJetFullHigh);
  hTriggers->SetBinContent(kJetLowGammaHigh, nJetFullLowGammaHigh / nJetFullLow);
  hTriggers->SetBinContent(kJetLowGammaLow, nJetFullLowGammaLow / nJetFullLow);

  setStyle(hTriggers, 0);
  hTriggers->SetStats(0);
  hTriggers->SetMinimum(0.);
  hTriggers->SetMaximum(1.);
  hTriggers->GetXaxis()->SetLabelSize(textSize);
  hTriggers->GetYaxis()->SetLabelSize(textSize);
  hTriggers->GetYaxis()->SetTitleSize(textSize);
  hTriggers->SetMarkerSize(2.);
  // hTriggers->SetBarOffset(0.6);

  Plotter p(x);
  p.textSize = textSize;
  if (drawText) p.setDrawOption("text0");
  // p.setDrawOption("bar");
  p.makeCanvas(1400, 600);
  p.frame = (TH1F*)hTriggers->Clone("frame");
  // p.frame->Reset();
  p.hists.push_back(hTriggers);

  p.addLatex(0.25, 0.80, sThisThesis.c_str());
  p.addLatex(0.25, 0.75, sData.c_str());
  p.addLatex(0.25, 0.70, sSqrtS.c_str());
  p.plot();
}

void jetGammaCorrelation2D(bool drawText = true) {
  vector<TriggerBins> triggers = {kJetFullLowPt, kJetFullHighPt, kGammaHighPtEMCal, kGammaLowPtEMCal};
  enum xBins {kJetLow = 1, kJetHigh, kNxPlusOne};
  enum yBins {kGammaLow = 1, kGammaHigh, kNyPlusOne};
  InputSettings x; x.verbosity = kDebug;
  x.train = 462605;
  x.setInputFileNameFromTrain();

  string histName = "triggerCorrelation_";
  for (auto trig : triggers) {
    histName += getTriggerName(trig) + "-";
  }
  x.outputFileName = histName + "2D.pdf";

  TFile* file = TFile::Open(x.inputFileName.c_str(), "READ");
  TH2* hTriggerCorrelation = (TH2*)file->Get(x.histName.c_str());

  TH2* hTriggers = new TH2D(histName.c_str(), "Fraction of jet trigger rate", kNxPlusOne - 1, 0, kNxPlusOne - 1, kNyPlusOne - 1, 0, kNyPlusOne - 1);
  hTriggers->GetXaxis()->SetBinLabel(kJetLow, getTriggerName(kJetFullLowPt).c_str());
  hTriggers->GetXaxis()->SetBinLabel(kJetHigh, getTriggerName(kJetFullHighPt).c_str());
  hTriggers->GetYaxis()->SetBinLabel(kGammaLow, getTriggerName(kGammaLowPtEMCal).c_str());
  hTriggers->GetYaxis()->SetBinLabel(kGammaHigh, getTriggerName(kGammaHighPtEMCal).c_str());

  double jetLowCorrection = 10, gammaLowCorrection = 20.;
  double nJetFullHigh          = hTriggerCorrelation->GetBinContent(kJetFullHighPt, kJetFullHighPt);
  double nJetFullLow           = hTriggerCorrelation->GetBinContent(kJetFullLowPt, kJetFullLowPt) * jetLowCorrection;
  double nGammaHigh            = hTriggerCorrelation->GetBinContent(kGammaHighPtEMCal, kGammaHighPtEMCal);
  double nGammaLow             = hTriggerCorrelation->GetBinContent(kGammaLowPtEMCal, kGammaLowPtEMCal) * gammaLowCorrection;
  double nJetFullHighGammaHigh = hTriggerCorrelation->GetBinContent(kJetFullHighPt, kGammaHighPtEMCal);
  double nJetFullHighGammaLow  = hTriggerCorrelation->GetBinContent(kJetFullHighPt, kGammaLowPtEMCal) * gammaLowCorrection;
  double nJetFullLowGammaHigh  = hTriggerCorrelation->GetBinContent(kJetFullLowPt, kGammaHighPtEMCal) * jetLowCorrection;
  double nJetFullLowGammaLow   = hTriggerCorrelation->GetBinContent(kJetFullLowPt, kGammaLowPtEMCal) * jetLowCorrection * gammaLowCorrection;

  hTriggers->SetBinContent(kJetLow, kGammaLow, nJetFullLowGammaLow / nJetFullLow);
  hTriggers->SetBinContent(kJetLow, kGammaHigh, nJetFullLowGammaHigh / nJetFullLow);
  hTriggers->SetBinContent(kJetHigh, kGammaLow, nJetFullHighGammaLow / nJetFullHigh);
  hTriggers->SetBinContent(kJetHigh, kGammaHigh, nJetFullHighGammaHigh / nJetFullHigh);

  setStyle(hTriggers, 0);
  hTriggers->SetStats(0);
  hTriggers->SetMinimum(0.);
  hTriggers->SetMaximum(1.);
  hTriggers->SetMarkerSize(2.);

  Plotter p(x);
  if (drawText) p.setDrawOption("text0");
  else p.setDrawOption("colz");

  p.makeCanvas(800, 600);
  gPad->SetLeftMargin(0.25);

  p.frame = (TH1F*)hTriggers->Clone("frame");
  p.frame->GetXaxis()->SetLabelSize(0.05);
  p.frame->GetYaxis()->SetLabelSize(0.05);
  p.frame->Reset();
  p.hists.push_back(hTriggers);

  p.addLatex(0.03, 0.90, sThisThesis.c_str());
  p.addLatex(0.03, 0.85, sData.c_str());
  p.addLatex(0.03, 0.80, sSqrtS.c_str());
  p.plot();
}

#endif
