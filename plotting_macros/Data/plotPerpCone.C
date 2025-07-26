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

namespace MyEnums {
  enum Verb {kErrors, kWarnings, kInfo, kDebug, kDebugMax};
  enum Foo {kGetJets, kGetJetV0s, kGetPCV0s};
}
using namespace MyEnums;

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

string InputSettings::getHistName(int x) {
  string s = "jet-fragmentation";
  if (train == 436232) {
    s += (x == kGetPCV0s) ? "_id24581" : "_id24580";
  }
  s += "/data/";
  switch (x) {
    case kGetJets:
      s += "jets/jetPtEtaPhi";
      break;
    case kGetJetV0s:
      s += "jets/V0/jetPtK0SPtMass";
      break;
    case kGetPCV0s:
      s += "PC/JetPtK0SPtMass";
      break;
    default:
      printLog("InputSettings::getHistName() Error: invalid x", kErrors);
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
void Plotter::makeFrame(double x0, double x1, double y0, double y1, string sx, string sy) {
  if (!canvas) makeCanvas();
  frame = DrawFrame(x0, x1, y0, y1, sx, sy);
}
void Plotter::makeLegend(double x0, double x1, double y0, double y1, string s) {
  legend = CreateLegend(x0, x1, y0, y1, s.c_str(), textSize);
}

void Plotter::plot() {
  if (hists.empty()) {
    string s = "Plotter::plot(): Hist vector is empty! Aborting";
    inputs->printLog(s, kErrors);
    return;
  }

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

struct PerpCone {
  InputSettings* inputs;
  Plotter* plotter;

  PerpCone() { inputs = new InputSettings(); plotter = new Plotter(*inputs); }
  PerpCone(InputSettings& i) { inputs = &i; plotter = new Plotter(*inputs); }

  TH1* getMassHist(int type);
  double getNjets();
  TH1* getPtHist(int type);

  void plotPerpConeMass();
  TH1* rebinPtHist(TH1* hist);
};

TH1* PerpCone::getMassHist(int type) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    inputs->printLog("Could not open file " + inputs->inputFileName, kErrors);
    return nullptr;
  }

  TH3* h3 = (TH3*)file->Get(inputs->getHistName(type).c_str());

  array<int, 2> jetBins = getProjectionBins(h3->GetXaxis(), inputs->ptminjet, inputs->ptmaxjet);
  array<int, 2> v0Bins = getProjectionBins(h3->GetYaxis(), inputs->ptmin, inputs->ptmax);
  inputs->lowptjet = h3->GetXaxis()->GetBinLowEdge(jetBins[0]);
  inputs->highptjet = h3->GetXaxis()->GetBinUpEdge(jetBins[1]);
  inputs->lowpt = h3->GetYaxis()->GetBinLowEdge(v0Bins[0]);
  inputs->highpt = h3->GetYaxis()->GetBinUpEdge(v0Bins[1]);

  string name = "m" + inputs->hadron;
  if (type == MyEnums::kGetPCV0s)
    name += "inPCs";
  if (type == MyEnums::kGetJetV0s)
    name += "inJets";
  name = inputs->getNameFromJetPt(name);
  name = inputs->getNameFromPt(name);

  TH1* hist = h3->ProjectionZ(name.c_str(), jetBins[0], jetBins[1], v0Bins[0], v0Bins[1]);
  return hist;
}

TH1* PerpCone::getPtHist(int type) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    inputs->printLog("Could not open file " + inputs->inputFileName, kErrors);
    return nullptr;
  }

  TH3* h3 = (TH3*)file->Get(inputs->getHistName(type).c_str());

  array<int, 2> jetBins = getProjectionBins(h3->GetXaxis(), inputs->ptminjet, inputs->ptmaxjet);
  inputs->lowptjet = h3->GetXaxis()->GetBinLowEdge(jetBins[0]);
  inputs->highptjet = h3->GetXaxis()->GetBinUpEdge(jetBins[1]);

  string name = "pt" + inputs->hadron;
  if (type == MyEnums::kGetPCV0s)
    name += "inPCs";
  if (type == MyEnums::kGetJetV0s)
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
  array<int, 2> jetBins = getProjectionBins(hJetPt->GetXaxis(), inputs->ptminjet, inputs->ptmaxjet);
  double nJets = hJetPt->Integral(jetBins[0], jetBins[1]);
  return nJets;
}

void PerpCone::plotPerpConeMass() {
  plotter->hists.clear();
  inputs->setInputFileNameFromTrain();
  inputs->outputFileName = inputs->getNameFromPt(inputs->getNameFromJetPt("pcMass"), ".pdf");

  plotter->makeLegend(0.65, 0.75, 0.7, 0.8, "");
  array<int, 2> jetBins, v0Bins;

  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");

  double nJets = getNjets();
  if (nJets < 0.) {
    inputs->printLog("PerpCone::plotPerpConeMass() Error: could not get number of jets", kErrors);
    return;
  }

  TH1* hJet = getMassHist(kGetJetV0s);
  setStyle(hJet, 0);
  hJet->Scale(1. / nJets, "width");
  plotter->hists.push_back(hJet);
  plotter->legend->AddEntry(hJet, "V0s in jets");

  TH1* hPC = getMassHist(kGetPCV0s);
  setStyle(hPC, 1);
  hPC->Scale(1. / (inputs->conesPerJet * nJets), "width");
  plotter->hists.push_back(hPC);
  plotter->legend->AddEntry(hPC, "V0s in UE");

  if (inputs->ratioplot) {
    plotter->makeFrame(0.4, 0.6, 0., 1., plotter->getMassString(), plotter->sRatio);
  } else {
    plotter->makeFrame(plotter->getMassString(), TString::Format("#frac{1}{%s} %s", plotter->sNjets.c_str(), plotter->getdNdXString("#it{M}").c_str()).Data());
  }

  plotter->addLatex(0.25, 0.8, "This Thesis, ALICE pp data");
  plotter->addLatex(0.25, 0.70, "#sqrt{s} = 13.6 TeV");
  plotter->addLatex(0.25, 0.65, TString::Format("Anti-#it{k}_{T} ch+V0 jets").Data());
  plotter->addLatex(0.25, 0.60, TString::Format("%.0f < #it{p}_{T, jet} < %.0f %s", inputs->lowptjet, inputs->highptjet, plotter->sGevC.c_str()).Data());
  plotter->addLatex(0.25, 0.55, TString::Format("%.1f < #it{p}_{T, V0} < %.1f %s", inputs->lowpt, inputs->highpt, plotter->sGevC.c_str()).Data());

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

void perpConeMass(bool doRatio) {
  InputSettings x; x.verbosity = kDebug;
  x.hadron = "K0S";
  x.train = 436232;
  x.ratioplot = doRatio;

  // Set jet pt
  x.setJetPt(20., 30.);
  x.setPtBinEdgesFromHadron();

  for (int ipt = 0; ipt < x.ptBinEdges.size(); ipt++) {
    x.setPt(x.ptBinEdges[ipt][0], x.ptBinEdges[ipt][1]);
    PerpCone p(x);
    p.plotPerpConeMass();
  }
}

void perpConePt() {
  InputSettings x; x.verbosity = kDebug;
  x.hadron = "K0S";
  x.train = 436232;
  x.setJetPt(20., 30.);
  x.logplot = true;

  // plotter->hists.clear();
  x.setInputFileNameFromTrain();
  x.outputFileName = x.getNameFromJetPt("pc", ".pdf");

  PerpCone p(x);
  p.plotter->makeLegend(0.45, 0.6, 0.3, 0.4, "");
  array<int, 2> jetBins;

  TFile* file = TFile::Open(x.inputFileName.c_str(), "READ");
  double nJets = p.getNjets();
  if (nJets < 0.) {
    x.printLog("PerpCone::plotPerpConeMass() Error: could not get number of jets", kErrors);
    return;
  }

  TH1* hJet = p.getPtHist(kGetJetV0s);
  hJet = p.rebinPtHist(hJet);
  setStyle(hJet, 0);
  hJet->Scale(1. / nJets, "width");
  p.plotter->hists.push_back(hJet);
  p.plotter->legend->AddEntry(hJet, "V0s in jets");

  TH1* hPC = p.getPtHist(kGetPCV0s);
  hPC = p.rebinPtHist(hPC);
  setStyle(hPC, 1);
  hPC->Scale(1. / (x.conesPerJet * nJets), "width");
  p.plotter->hists.push_back(hPC);
  p.plotter->legend->AddEntry(hPC, "V0s in UE");

  if (x.ratioplot) {
    p.plotter->makeFrame(0.4, 0.6, 0., 1., p.plotter->getMassString(), p.plotter->sRatio);
  } else {
    string xTitle = p.plotter->sPtJet;
    string yTitle = TString::Format("#frac{1}{#it{N}_{jets, cones}} %s", p.plotter->getdNdXString("#it{p}_{T}").c_str()).Data();
    p.plotter->makeFrame(0., 30., 1e-5, 0.2, xTitle, yTitle);
    // p.plotter->makeFrame(0., 30., 1e-1, 1e5, xTitle, yTitle);
  }

  p.plotter->addLatex(0.55, 0.85, "This Thesis");
  p.plotter->addLatex(0.55, 0.80, "ALICE pp data");
  p.plotter->addLatex(0.55, 0.75, "#sqrt{s} = 13.6 TeV");
  p.plotter->addLatex(0.55, 0.70, TString::Format("Anti-#it{k}_{T} ch+V0 jets").Data());
  p.plotter->addLatex(0.55, 0.65, TString::Format("#it{R} = 0.4, |#eta| < 0.5").Data());
  p.plotter->addLatex(0.55, 0.60, TString::Format("%.0f < #it{p}_{T, jet} < %.0f %s", x.lowptjet, x.highptjet, p.plotter->sGevC.c_str()).Data());

  p.plotter->plot();
}

#endif
