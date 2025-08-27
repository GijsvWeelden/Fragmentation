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

#ifndef __PLOTEFFICIENCY_H__
#define __PLOTEFFICIENCY_H__

// -------------------------------------------------------------------------------------------------
//
// Struct for input settings
//
// -------------------------------------------------------------------------------------------------

namespace MyEnums {
  enum Verb {kErrors, kWarnings, kInfo, kDebug, kDebugMax};
  enum Level {kGenerator, kReconstructed, kGeneratorJet, kReconstructedJet, kGeneratorJetMatched, kReconstructedJetMatched, kReconstructedPileUpOnly, kReconstructedJetPileUpOnly, kReconstructedJetMatchedPileUpOnly};
  enum PlotType {kPlotPt, kPlotEff, kPlotPtJets, kPlotEffJets, kPlotPtJetsMatched, kPlotEffJetsMatched};
}
using namespace MyEnums;

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
    double ptminjet = -1e3, ptmaxjet = -1e3, lowptjet = -1e3, highptjet = -1e3;
    double etamin = -0.9, etamax = 0.9, loweta = etamin, higheta = etamax;
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
    void setEta(double a, double b);
    void setJetPt(double a, double b);
    void setJetPt(vector<double> bin);
    void setPt(double a, double b);
    vector<vector<double>> setPtBinEdgesFromHadron();
    vector<vector<double>> setPtBinEdgesSorted(vector<vector<double>> x);
    template <typename T> int writeOutputToFile(T* obj);
};

string InputSettings::getHistName(int x) {
  string s = "jet-v0qa";

  switch (train) {
    case 436064:
      s += "_id13556";
      break;
    case 439670:
      if (x == kReconstructedPileUpOnly)
        s += "_id31538";
      else
        s += "_id31540";
      break;
  }

  switch (x) {
    case kGenerator:
      s += "/inclusive/Generated" + hadron;
      break;
    case kReconstructed:
      s += "/inclusive/" + hadron + "PtEtaMass";
      break;
    case kGeneratorJet:
    case kGeneratorJetMatched:
      s += "/jets/GeneratedJet" + hadron;
      break;
    case kReconstructedJet:
      s += "/jets/JetPtEta" + hadron + "Pt";
      break;
    case kReconstructedJetMatched:
      s += "/jets/JetsPtEta" + hadron + "Pt";
      break;
    case kReconstructedPileUpOnly:
      s += "/inclusive/" + hadron + "PtEtaMassWrongCollision";
      break;
    case kReconstructedJetPileUpOnly:
      s += "/jets/JetPtEta" + hadron + "PtWrongCollision";
      break;
    case kReconstructedJetMatchedPileUpOnly:
      s += "/jets/JetsPtEta" + hadron + "PtWrongCollision";
      break;
    default:
      printLog("InputSettings::getHistName() Error: invalid setting " + to_string(x), kErrors);
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

string InputSettings::printLog(string message, int verbThreshold) {
  if (this->verbosity < verbThreshold)
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

void InputSettings::setEta(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setEta() Error: etamin > etamax", kErrors);
    return;
  }
  this->etamin = a;
  this->etamax = b;
  this->loweta = a;
  this->higheta = b;
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
    printLog("InputSettings::setJetPt() Error: ptminjet > ptmaxjet", kErrors);
    return;
  }
  this->ptminjet = a;
  this->ptmaxjet = b;
  this->lowptjet = a;
  this->highptjet = b;
}

void InputSettings::setJetPt(vector<double> bin) {
  setJetPt(bin[0], bin[1]);
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

    Plotter() { inputs = new InputSettings(); }
    Plotter(InputSettings& x) { inputs = &x; }

    void addLatex(double x, double y, string s);
    void addLine(double x0, double y0, double x1, double y1, int styleNumber, int lineStyle, int lineWidth);
    void makeCanvas(string s, double x, double y);
    void makeFrame(string sx, string sy);
    void makeFrame(double x0, double x1, double y0, double y1, string sx, string sy);
    void makeLegend(double x0, double x1, double y0, double y1, string s);
    void plot();
    void setHistStyles();

    // String formatting
    string addUnits(string s, string unit);
    string getdYdXString(string y, string x);
    string getPtString(string subscript = "");
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

string Plotter::addUnits(string s, string unit) {
  return TString::Format("%s (%s)", s.c_str(), unit.c_str()).Data();
}

string Plotter::getdYdXString(string y, string x) {
  return TString::Format("#frac{d%s}{d%s}", y.c_str(), x.c_str()).Data();
}

string Plotter::getPtString(string subscript) {
  if (subscript.empty())
    return TString::Format("#it{p}_{T}").Data();
  else
    return TString::Format("#it{p}_{T, %s}", subscript.c_str()).Data();
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
// Helper struct for efficiencies
//
// -------------------------------------------------------------------------------------------------

struct EfficiencyFinder {
  InputSettings* inputs;
  Plotter* plotter;

  EfficiencyFinder() { inputs = new InputSettings(); plotter = new Plotter(*inputs); }
  EfficiencyFinder(InputSettings& i) { inputs = &i; plotter = new Plotter(*inputs); }

  TH1* getPtHist(int type, string name);
  std::array<TH1*, 2> getPtHistRebinnedAndNormalised(int type);

  void plotIncl(bool plotSpectra, bool subPileUp, bool rebinHists);
  void plotPtIncl(bool subPileUp, bool rebinHists);
  void plotEffIncl(bool subPileUp, bool rebinHists);

  void plotJets(double ptmin, double ptmax, bool plotSpectra, bool subPileUp, bool plotMatched, bool rebinHists);
  void plotPtJets(double ptmin, double ptmax, bool subPileUp, bool rebinHists);
  void plotPtJetsMatched(double ptmin, double ptmax, bool subPileUp, bool rebinHists);
  void plotEffJets(double ptmin, double ptmax, bool subPileUp, bool rebinHists);
  void plotEffJetsMatched(double ptmin, double ptmax, bool subPileUp, bool rebinHists);
};

TH1* EfficiencyFinder::getPtHist(int type, string name) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    inputs->printLog("EfficiencyFinder::getPtHist() Error: could not open file " + inputs->inputFileName, kErrors);
    return nullptr;
  }
  inputs->printLog("Loading histogram " + inputs->getHistName(type), kDebug);

  switch (type) {
    case kGenerator:
    case kReconstructed:
    case kReconstructedPileUpOnly:
    {
      TH3* h3 = (TH3*)file->Get(inputs->getHistName(type).c_str());
      array<int, 2> etaBins = getProjectionBins(h3->GetYaxis(), inputs->etamin, inputs->etamax);
      TH1* h = h3->ProjectionX(name.c_str(), etaBins[0], etaBins[1], 1, h3->GetNbinsZ());
      return h;
    }
    case kGeneratorJet:
    case kReconstructedJet:
    case kGeneratorJetMatched:
    case kReconstructedJetPileUpOnly:
    {
      TH3* h3 = (TH3*)file->Get(inputs->getHistName(type).c_str());
      array<int, 2> jetBins = getProjectionBins(h3->GetXaxis(), inputs->ptminjet, inputs->ptmaxjet);
      array<int, 2> etaBins = getProjectionBins(h3->GetYaxis(), inputs->etamin, inputs->etamax);
      TH1* h = h3->ProjectionZ(name.c_str(), jetBins[0], jetBins[1], etaBins[0], etaBins[1]);
      return h;
    }
    case kReconstructedJetMatched:
    case kReconstructedJetMatchedPileUpOnly:
    {
      THnSparse* hn = (THnSparse*)file->Get(inputs->getHistName(type).c_str());
      const int ptPartAxis = 0, ptDetAxis = 1, etaDetAxis = 2, ptV0Axis = 3;
      array<int, 2> jetBins = getProjectionBins(hn->GetAxis(ptDetAxis), inputs->ptminjet, inputs->ptmaxjet);
      array<int, 2> etaBins = getProjectionBins(hn->GetAxis(etaDetAxis), inputs->etamin, inputs->etamax);
      hn->GetAxis(ptDetAxis)->SetRange(jetBins[0], jetBins[1]);
      hn->GetAxis(etaDetAxis)->SetRange(etaBins[0], etaBins[1]);
      TH1* h = hn->Projection(ptV0Axis);
      h->SetName(name.c_str());
      return h;
    }
    default:
      inputs->printLog("EfficiencyFinder::getPtHist() Error: invalid type " + to_string(type), kErrors);
      return nullptr;
  }
}

std::array<TH1*, 2> EfficiencyFinder::getPtHistRebinnedAndNormalised(int type) {
  TH1* hRec = nullptr;
  TH1* hGen = nullptr;

  switch (type) {
    case kReconstructed:
    case kGenerator:
      hRec = getPtHist(kReconstructed, "hRec");
      hGen = getPtHist(kGenerator, "hGen");
      break;
    case kReconstructedJet:
    case kGeneratorJet:
      hRec = getPtHist(kReconstructedJet, "hRec");
      hGen = getPtHist(kGeneratorJet, "hGen");
      break;
    case kReconstructedJetMatched:
    case kGeneratorJetMatched:
      hRec = getPtHist(kReconstructedJetMatched, "hRec");
      hGen = getPtHist(kGeneratorJetMatched, "hGen");
      break;
    default:
      inputs->printLog("EfficiencyFinder::getPtHistRebinnedAndNormalised() Error: invalid type " + to_string(type), kErrors);
      return std::array{hRec, hGen};
  }

  TH1* hRecRebinned = rebinHist(hRec, rebinnedV0PtHist(inputs->hadron, "hRecRebinned"));
  TH1* hGenRebinned = rebinHist(hGen, rebinnedV0PtHist(inputs->hadron, "hGenRebinned"));

  // FIXME: Am I doing this right?
  double crossSection = hGenRebinned->Integral(1, hGenRebinned->GetNbinsX(), "width");
  hGenRebinned->Scale(1. / crossSection, "width");
  hRecRebinned->Scale(1. / crossSection, "width");

  return std::array{hRecRebinned, hGenRebinned};
}

void EfficiencyFinder::plotIncl(bool plotSpectra, bool subPileUp, bool rebinHists) {
  plotter->hists.clear();
  TH1* hRec = getPtHist(kReconstructed, "hRec");
  TH1* hGen = getPtHist(kGenerator, "hGen");

  if (subPileUp) {
    TH1* hRecPileUp = getPtHist(kReconstructedPileUpOnly, "hRecPileUp");
    hRec->Add(hRecPileUp, -1);
  }

  if (rebinHists) {
    TH1* tmp = rebinHist(hRec, rebinnedV0PtHist(inputs->hadron, "hRecRebinned"));
    hRec = (TH1*)tmp->Clone();
    tmp = rebinHist(hGen, rebinnedV0PtHist(inputs->hadron, "hGenRebinned"));
    hGen = (TH1*)tmp->Clone();
  }

  double xSec = hGen->Integral(1, hGen->GetNbinsX(), "width");
  hGen->Scale(1. / xSec, "width");
  hRec->Scale(1. / xSec, "width");

  if (plotSpectra) {
    inputs->logplot = true;
    plotter->hists.push_back(hRec);
    plotter->hists.push_back(hGen);
    plotter->setHistStyles();

    plotter->makeLegend(0.25, 0.45, 0.25, 0.35, "");
    plotter->legend->AddEntry(hGen, "Generated");
    plotter->legend->AddEntry(hRec, "Reconstructed");

    string xTitle = plotter->getPtString(plotter->sV0);
    string yTitle = TString::Format("#frac{1}{%s} %s", plotter->sSigma.c_str(), plotter->getdYdXString(plotter->sSigma, xTitle).c_str()).Data();
    plotter->makeFrame(0., 40., 1e-9, 10., xTitle, yTitle);

    plotter->addLatex(0.45, 0.83, TString::Format("%s, %s", plotter->sThisThesis.c_str(), plotter->sPythia.c_str()).Data());
    plotter->addLatex(0.45, 0.78, plotter->sSqrtS.c_str());
    plotter->addLatex(0.45, 0.73, TString::Format("|#eta| < %.1f", inputs->etamax).Data());
  } else { // Efficiency
    TH1* hEff = (TH1*)hRec->Clone("hEff");
    hEff->Divide(hGen);
    plotter->hists.push_back(hEff);
    plotter->setHistStyles();

    string xTitle = plotter->getPtString(plotter->sV0);
    string yTitle = "Efficiency";

    plotter->makeFrame(0., 40., 0., 1.1, xTitle, yTitle);
    plotter->addLatex(0.63, 0.80, plotter->sThisThesis.c_str());
    plotter->addLatex(0.63, 0.75, plotter->sPythia.c_str());
    plotter->addLatex(0.63, 0.70, plotter->sSqrtS.c_str());
    plotter->addLatex(0.63, 0.65, TString::Format("|#eta| < %.1f", inputs->etamax).Data());
  }

  plotter->plot();
}

void EfficiencyFinder::plotPtIncl(bool subPileUp, bool rebinHists) {
  plotIncl(true, subPileUp, rebinHists);
}

void EfficiencyFinder::plotEffIncl(bool subPileUp, bool rebinHists) {
  plotIncl(false, subPileUp, rebinHists);
}

void EfficiencyFinder::plotJets(double ptmin, double ptmax, bool plotSpectra, bool subPileUp, bool plotMatched, bool rebinHists) {
  plotter->hists.clear();
  int jetType = plotMatched ? kReconstructedJetMatched : kReconstructedJet;
  TH1* hRec = getPtHist(jetType, "hRec");
  jetType = plotMatched ? kGeneratorJetMatched : kGeneratorJet;
  TH1* hGen = getPtHist(jetType, "hGen");

  if (subPileUp) {
    jetType = plotMatched ? kReconstructedJetMatchedPileUpOnly : kReconstructedJetPileUpOnly;
    TH1* hRecPileUp = getPtHist(jetType, "hRecPileUp");
    hRec->Add(hRecPileUp, -1);
  }

  if (rebinHists) {
    TH1* tmp = rebinHist(hRec, rebinnedV0PtHist(inputs->hadron, "hRecRebinned"));
    hRec = (TH1*)tmp->Clone();
    tmp = rebinHist(hGen, rebinnedV0PtHist(inputs->hadron, "hGenRebinned"));
    hGen = (TH1*)tmp->Clone();
  }

  // FIXME: Need to divide by njets instead!
  double xSec = hGen->Integral(1, hGen->GetNbinsX(), "width");
  hGen->Scale(1. / xSec, "width");
  hRec->Scale(1. / xSec, "width");

  if (plotSpectra) {
    inputs->logplot = true;
    plotter->hists.push_back(hRec);
    plotter->hists.push_back(hGen);
    plotter->setHistStyles();

    plotter->makeLegend(0.25, 0.45, 0.25, 0.35, "");
    plotter->legend->AddEntry(hGen, "Generated");
    plotter->legend->AddEntry(hRec, "Reconstructed");

    string xTitle = plotter->getPtString(plotter->sV0);
    string yTitle = TString::Format("#frac{1}{%s} %s", plotter->sSigma.c_str(), plotter->getdYdXString(plotter->sSigma, xTitle).c_str()).Data();
    plotter->makeFrame(0., ptmax, 1e-5, 1., xTitle, yTitle);

    plotter->addLatex(0.45, 0.83, TString::Format("%s, %s", plotter->sThisThesis.c_str(), plotter->sPythia.c_str()).Data());
    plotter->addLatex(0.45, 0.78, plotter->sSqrtS.c_str());
    plotter->addLatex(0.45, 0.73, TString::Format("%s, %s, |#eta| < %.1f", plotter->sAntiktJets.c_str(), plotter->sRadius.c_str(), inputs->etamax).Data());
    plotter->addLatex(0.45, 0.68, TString::Format("%.1f < %s < %.1f %s", inputs->ptminjet, plotter->getPtString("jet").c_str(), inputs->ptmaxjet, plotter->sGevC.c_str()).Data());
  } else { // Efficiency
    TH1* hEff = (TH1*)hRec->Clone("hEff");
    hEff->Divide(hGen);
    plotter->hists.push_back(hEff);
    plotter->setHistStyles();

    string xTitle = plotter->getPtString(plotter->sV0);
    string yTitle = "Efficiency";
    plotter->makeFrame(0., ptmax, 0., 1.1, xTitle, yTitle);

    plotter->addLatex(0.63, 0.80, plotter->sThisThesis.c_str());
    plotter->addLatex(0.63, 0.75, plotter->sPythia.c_str());
    plotter->addLatex(0.63, 0.70, plotter->sSqrtS.c_str());
    plotter->addLatex(0.63, 0.65, TString::Format("%s, %s, |#eta| < %.1f", plotter->sAntiktJets.c_str(), plotter->sRadius.c_str(), inputs->etamax).Data());
    plotter->addLatex(0.63, 0.60, TString::Format("%.1f < %s < %.1f %s", inputs->ptminjet, plotter->getPtString("jet").c_str(), inputs->ptmaxjet, plotter->sGevC.c_str()).Data());
  }
  plotter->plot();
}

void EfficiencyFinder::plotPtJets(double ptmin, double ptmax, bool subPileUp, bool rebinHists) {
  bool plotSpectra = true, plotMatched = false;
  plotJets(ptmin, ptmax, plotSpectra, subPileUp, plotMatched, rebinHists);
}
void EfficiencyFinder::plotPtJetsMatched(double ptmin, double ptmax, bool subPileUp, bool rebinHists) {
  bool plotSpectra = true, plotMatched = true;
  plotJets(ptmin, ptmax, plotSpectra, subPileUp, plotMatched, rebinHists);
}
void EfficiencyFinder::plotEffJets(double ptmin, double ptmax, bool subPileUp, bool rebinHists) {
  bool plotSpectra = false, plotMatched = false;
  plotJets(ptmin, ptmax, plotSpectra, subPileUp, plotMatched, rebinHists);
}
void EfficiencyFinder::plotEffJetsMatched(double ptmin, double ptmax, bool subPileUp, bool rebinHists) {
  bool plotSpectra = false, plotMatched = true;
  plotJets(ptmin, ptmax, plotSpectra, subPileUp, plotMatched, rebinHists);
}

// -------------------------------------------------------------------------------------------------
//
// Interface
//
// -------------------------------------------------------------------------------------------------

void doPlotting(PlotType mode, bool subPileUp, bool rebinHists, double ptmin = -1., double ptmax = -1.) {
  InputSettings x; x.verbosity = kDebug;
  x.hadron = "K0S";
  x.train = 476307;
  // 476306: MB
  // 475816: Jet-Jet

  x.setJetPt(ptmin, ptmax);
  x.setInputFileNameFromTrain();
  EfficiencyFinder ef(x);

  switch (mode) {
    case kPlotPt:
      x.outputFileName = x.hadron + "_pt" + (subPileUp ? "_subPU" : "") + ".pdf";
      x.setEta(-0.9, 0.9);
      ef.plotPtIncl(subPileUp, rebinHists);
      break;
    case kPlotEff:
      x.outputFileName = x.hadron + "_eff" + (subPileUp ? "_subPU" : "") + ".pdf";
      x.setEta(-0.9, 0.9);
      ef.plotEffIncl(subPileUp, rebinHists);
      break;
    case kPlotPtJets:
      x.outputFileName = x.getNameFromJetPt(x.hadron + "_pt", (subPileUp ? "_subPU" : "")) + ".pdf";
      x.setEta(-0.5, 0.5);
      ef.plotPtJets(ptmin, ptmax, subPileUp, rebinHists);
      break;
    case kPlotEffJets:
      x.outputFileName = x.getNameFromJetPt(x.hadron + "_eff", (subPileUp ? "_subPU" : "")) + ".pdf";
      x.setEta(-0.5, 0.5);
      ef.plotEffJets(ptmin, ptmax, subPileUp, rebinHists);
      break;
    case kPlotPtJetsMatched:
      x.outputFileName = x.getNameFromJetPt(x.hadron + "_pt", "_matched") + (subPileUp ? "_subPU" : "") + ".pdf";
      x.setEta(-0.5, 0.5);
      ef.plotPtJetsMatched(ptmin, ptmax, subPileUp, rebinHists);
      break;
    case kPlotEffJetsMatched:
      x.outputFileName = x.getNameFromJetPt(x.hadron + "_eff", "_matched") + (subPileUp ? "_subPU" : "") + ".pdf";
      x.setEta(-0.5, 0.5);
      ef.plotEffJetsMatched(ptmin, ptmax, subPileUp, rebinHists);
      break;
    default:
      x.printLog("doPlotting() Error: invalid mode " + to_string(mode), kErrors);
  }
}

void plotpt(bool subPileUp, bool rebinHists) { doPlotting(kPlotPt, subPileUp, rebinHists); }
void ploteff(bool subPileUp, bool rebinHists) { doPlotting(kPlotEff, subPileUp, rebinHists); }

void plotptjets(double ptmin, double ptmax, bool subPileUp, bool rebinHists) {
  doPlotting(kPlotPtJets, subPileUp, rebinHists, ptmin, ptmax);
}
void ploteffjets(double ptmin, double ptmax, bool subPileUp, bool rebinHists) {
  doPlotting(kPlotEffJets, subPileUp, rebinHists, ptmin, ptmax);
}
void plotptjetsmatched(double ptmin, double ptmax, bool subPileUp, bool rebinHists) {
  doPlotting(kPlotPtJetsMatched, subPileUp, rebinHists, ptmin, ptmax);
}
void ploteffjetsmatched(double ptmin, double ptmax, bool subPileUp, bool rebinHists) {
  doPlotting(kPlotEffJetsMatched, subPileUp, rebinHists, ptmin, ptmax);
}

void allmyplots(bool subPileUp = true, bool rebinHists = true) {
  gROOT->SetBatch();

  plotpt(subPileUp, rebinHists);
  ploteff(subPileUp, rebinHists);

  plotptjets(10., 20., subPileUp, rebinHists);
  plotptjets(20., 30., subPileUp, rebinHists);

  ploteffjets(10., 20., subPileUp, rebinHists);
  ploteffjets(20., 30., subPileUp, rebinHists);

  plotptjetsmatched(10., 20., subPileUp, rebinHists);
  plotptjetsmatched(20., 30., subPileUp, rebinHists);

  ploteffjetsmatched(10., 20., subPileUp, rebinHists);
  ploteffjetsmatched(20., 30., subPileUp, rebinHists);
}

// void plotEffComparison(Level plotLevel, bool subPileUp = true, bool rebinHists = true) {
//   vector<int> trains = {476306, 475816};
//   vector<vector<double>> jetptbins = { {10., 20.}, {20., 30.} };
//   // Plot eff for both trains or both pt bins
//   InputSettings x; x.verbosity = kDebug;
//   x.train = trains[0];
//   x.setInputFileNameFromTrain();
//   x.setJetPt(jetptbins[0]);

//   Level lRec =
//   TH1* hRec = getPtHist(kReconstructed, "hRec");
//   TH1* hGen = getPtHist(kGenerator, "hGen");
// }



//   if (subPileUp) {
//     TH1* hRecPileUp = getPtHist(kReconstructedPileUpOnly, "hRecPileUp");
//     hRec->Add(hRecPileUp, -1);
//   }

//   if (rebinHists) {
//     TH1* tmp = rebinHist(hRec, rebinnedV0PtHist(inputs->hadron, "hRecRebinned"));
//     hRec = (TH1*)tmp->Clone();
//     tmp = rebinHist(hGen, rebinnedV0PtHist(inputs->hadron, "hGenRebinned"));
//     hGen = (TH1*)tmp->Clone();
//   }

//   int jetType = plotMatched ? kReconstructedJetMatched : kReconstructedJet;
//   TH1* hRec = getPtHist(jetType, "hRec");
//   jetType = plotMatched ? kGeneratorJetMatched : kGeneratorJet;
//   TH1* hGen = getPtHist(jetType, "hGen");

//   if (subPileUp) {
//     jetType = plotMatched ? kReconstructedJetMatchedPileUpOnly : kReconstructedJetPileUpOnly;
//     TH1* hRecPileUp = getPtHist(jetType, "hRecPileUp");
//     hRec->Add(hRecPileUp, -1);
//   }

//   if (rebinHists) {
//     TH1* tmp = rebinHist(hRec, rebinnedV0PtHist(inputs->hadron, "hRecRebinned"));
//     hRec = (TH1*)tmp->Clone();
//     tmp = rebinHist(hGen, rebinnedV0PtHist(inputs->hadron, "hGenRebinned"));
//     hGen = (TH1*)tmp->Clone();
//   }


#endif
