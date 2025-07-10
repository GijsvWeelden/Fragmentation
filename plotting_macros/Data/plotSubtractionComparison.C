
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

#ifndef SUBTRACTION_COMPARISON_C
#define SUBTRACTION_COMPARISON_C

int getHistScaleBinInRange(TH1* h, int minbin, int maxbin, bool doMin = false, bool doError = false) {
  // Initialise with max/min bin content to ensure scale is overwritten when checking min/max
  int bin = doMin ? h->GetMaximumBin() : h->GetMinimumBin();
  for (int i = 1; i < h->GetNbinsX(); i++) {
    double scale = h->GetBinContent(bin);
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);

    double s = bc;
    if (doError)  {
      if (doMin) s -= be;
      else s += be;
    }

    if (doMin && scale < s) bin = i;
    else if (scale > s) bin = i;
  }
  return bin;
}
double getHistScaleInRange(TH1* h, int minbin, int maxbin, bool doMin = false, bool doError = false) {
  // Initialise with max/min bin content to ensurescale is overwritten when checking min/max
  int bin = getHistScaleBinInRange(h, minbin, maxbin, doMin, doError);
  double bc = h->GetBinContent(bin);
  double scale = bc;
  if (doError) {
    double be = h->GetBinError(bin);
    scale = doMin ? bc - be : bc + be;
  }

  return scale;
}
int getHistLowBinInRange(TH1* hist, int minbin, int maxbin, bool doError = false) {
  return getHistScaleBinInRange(hist, minbin, maxbin, true, doError);
}
double gethistLowBoundInRange(TH1* hist, int minbin, int maxbin, bool doError = false) {
  return getHistScaleBinInRange(hist, minbin, maxbin, true, doError);
}

double roundToPowerOfTen(double x, int round = 0) {
  double y = std::log10(x);
  if (round == 1)
    y = std::ceil(y);
  else if (round == 2)
    y = std::floor(y);
  else
    y = std::round(y);

  return std::pow(10., y);
}
double roundToNextPowerOfTen(double x) {
  return roundToPowerOfTen(x, 1);
}
double roundToPrevPowerOfTen(double x) {
  return roundToPowerOfTen(x, 2);
}

namespace MyEnums {
  enum Verbosity {kErrors, kWarnings, kInfo, kDebug};
  enum Subtraction {kNone, kSubtract, kAddTracks};
}

using namespace MyEnums;

struct InputSettings {
  private:
  public:
    int train, rebinNumber, verbosity = kWarnings, subtraction = -1;
    string inputFileName, outputFileName;
    string hadron, histNameSub, histNameTrackAdd;

    double jetptmin, jetptmax, jetptlow, jetpthigh;
    double ptmin, ptmax, ptlow, pthigh;

    bool logplot = false, ratioplot = false;

    string getFragHistNameFromTrain(int x);
    string getJetHistNameFromTrain(int x);
    double getMass();
    string getSaveNameFromPt(string prefix, string suffix);
    string getSaveNameFromJetPt(string prefix, string suffix);
    string printLog(string message, int verbThreshold);
    string setInputFileNameFromTrain();
    void setJetPt(double a, double b);
    void setPt(double a, double b);
    template <typename T> int writeOutputToFile(T* obj);
};

string InputSettings::getFragHistNameFromTrain(int x) {
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
string InputSettings::getJetHistNameFromTrain(int x) {
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
    const string sPtJet = "#it{p}_{T, jet} (GeV/#it{c})";
    const string sCounts = "Counts";
    const string sZV0 = "#it{z}_{V0}";
    const string sRatio = "Ratio";

    Plotter() { inputs = new InputSettings(); }
    Plotter(InputSettings& x) { inputs = &x; }

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
  // latex.push_back(l);
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
  double yMinFrame = getHistLowerBound(hists, 0, 0) * 0.9;
  double yMaxFrame = getHistUpperBound(hists, 0, 0) * 1.2;

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
    string xTitle = sPtJet;
    string yTitle = sCounts;
    if (inputs->ratioplot) yTitle = sRatio;
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

// Plot jet spectrum
void plotJetPt(bool doBase, bool doRatio) {
  InputSettings x; //x.verbosity = kDebug;
  x.train = 436232;

  x.setInputFileNameFromTrain();
  x.ratioplot = false;
  x.logplot = true;

  Plotter p(x);
  if (!doBase && doRatio) {
    p.setDrawOption("hist");
    p.makeFrame(0., 200., 0.8, 1.2, p.sPtJet, p.sRatio);
  }

  // Spectrum - Base vs V0 Sub vs Track-Adding
  x.outputFileName = "jetspectrum";
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

  TH3* h3Base = (TH3*)file->Get(x.getJetHistNameFromTrain(MyEnums::kNone).c_str());
  TH1* base = (TH1*)h3Base->ProjectionX("base");
  if (doBase) {
    p.hists.push_back(base);
    p.legend->AddEntry(base, "No subtraction");
  }

  TH3* h3Sub = (TH3*)file->Get(x.getJetHistNameFromTrain(MyEnums::kSubtract).c_str());
  TH1* sub = (TH1*)h3Sub->ProjectionX("sub");
  p.hists.push_back(sub);
  p.legend->AddEntry(sub, "Subtract V0s");

  TH3* h3TrackAdd = (TH3*)file->Get(x.getJetHistNameFromTrain(MyEnums::kAddTracks).c_str());
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
void plotJetPt() {
  cout << "Plotting all jet pt figures" << endl;
  plotJetPt(0, 0);
  plotJetPt(0, 1);
  plotJetPt(1, 0);
  plotJetPt(1, 1);
}

void plotV0z(bool doBase, bool doRatio) {
  InputSettings x;
  x.hadron = "K0S";
  x.train = 436232;

  x.setInputFileNameFromTrain();
  x.ratioplot = false;
  x.logplot = true;
  x.setJetPt(20., 30.);
  const int projectionAxis = 1;

  x.outputFileName = "v0z";
  if (doRatio) {
    x.ratioplot = true;
    x.logplot = false;
    x.outputFileName += "ratio";
  }
  if (doBase) x.outputFileName += "-all";
  x.outputFileName += ".pdf";

  if (x.verbosity >= MyEnums::kDebug) {
    cout << "in: " << x.inputFileName << "\nout: " << x.outputFileName
    << "\nBase: " << x.getFragHistNameFromTrain(MyEnums::kNone)
    << "\nSub: " << x.getFragHistNameFromTrain(MyEnums::kSubtract)
    << "\nAdd: " << x.getFragHistNameFromTrain(MyEnums::kAddTracks)
    << endl;
  }

  Plotter p(x);
  p.makeCanvas();
  p.makeLegend(0.25, 0.45, 0.17, 0.37, "");

  if (doRatio && doBase) {
    p.makeFrame(1e-3, 1+1e-3, 0.2, 1.1, p.sZV0, p.sRatio);
    p.addLatex(0.4, 0.70, "This Thesis");
    p.addLatex(0.4, 0.65, "Anti-#it{k}_{T} ch+V0 jets, R = 0.4");
    p.addLatex(0.4, 0.60, TString::Format("%.0f < #it{p}_{T, jet} < %.0f GeV/#it{c}", x.jetptlow, x.jetpthigh).Data());
  } else if (doRatio) {
    p.makeFrame(1e-3, 1+1e-3, 0.5, 1.7, p.sZV0, p.sRatio);
    p.addLatex(0.25, 0.80, "This Thesis");
    p.addLatex(0.25, 0.75, "Anti-#it{k}_{T} ch+V0 jets, R = 0.4");
    p.addLatex(0.25, 0.70, TString::Format("%.0f < #it{p}_{T, jet} < %.0f GeV/#it{c}", x.jetptlow, x.jetpthigh).Data());
  } else {
    p.makeFrame(1e-3, 1+1e-3, 10, 2e4, p.sZV0, "\"Counts\"");
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

  THnSparse* hnBase = (THnSparse*)file->Get(x.getFragHistNameFromTrain(MyEnums::kNone).c_str());
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

  THnSparse* hnSub = (THnSparse*)file->Get(x.getFragHistNameFromTrain(MyEnums::kSubtract).c_str());
  jetptbins = getProjectionBins(hnSub->GetAxis(0), x.jetptmin, x.jetptmax);
  hnSub->GetAxis(0)->SetRange(jetptbins[0], jetptbins[1]);
  TH1* sub = (TH1*)hnSub->Projection(projectionAxis);
  sub->SetName(TString::Format("sub_%s_jetpt%.0f-%.0f", x.hadron.c_str(), x.jetptlow, x.jetpthigh).Data());
  if (x.verbosity >= MyEnums::kDebug) sub->Print();
  p.hists.push_back(sub);
  p.legend->AddEntry(sub, "Subtract V0s");

  THnSparse* hnAdd = (THnSparse*)file->Get(x.getFragHistNameFromTrain(MyEnums::kAddTracks).c_str());
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
void plotV0z() {
  cout << "Plotting all V0 z figures" << endl;
  plotV0z(0, 0);
  plotV0z(0, 1);
  plotV0z(1, 0);
  plotV0z(1, 1);
}

#endif // SUBTRACTION_COMPARISON_C
