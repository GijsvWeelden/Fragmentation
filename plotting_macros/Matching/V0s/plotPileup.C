
#include "../../histUtils.C"
#include "../../plotUtils.C"
#include "../../myStrings.C"

#ifndef __PLOTPERPCONE_H__
#define __PLOTPERPCONE_H__

namespace verbosityutils {
  enum Verbosity {kErrors, kWarnings, kInfo, kDebug, kDebugMax};
  // Name of the verbosity as written in messages
  string to_string(Verbosity v) {
    switch (v) {
      case kErrors:   return "Error";
      case kWarnings: return "Warning";
      case kInfo:     return "Info";
      case kDebug:    return "Debug";
      case kDebugMax: return "DebugMax";
      default:        return "Unknown";
    }
  }
  bool passVerbosityCheck(Verbosity level, Verbosity threshold) { return ( level <= threshold); }
  void printLog(string message, Verbosity level) { cout << to_string(level) << ": " << message << endl; }
} // namespace verbosityutils

namespace typeutils {
  enum histtypes { kNevts, kNjets, kInclusiveV0s, kInclusiveV0sWrongCollision, kInJetsV0Pt, kInJetsV0PtWrongCollision, kInJetsV0Z, kInJetsV0ZWrongCollision, kInMatchedJetsV0Pt, kInMatchedJetsV0PtWrongCollision, kInMatchedJetsV0Z, kInMatchedJetsV0ZWrongCollision };
} // namespace typeutils

struct InputSettings {
  private:
    int _train;
    double _etamin, _etamax;
    double _ptjetmin, _ptjetmax;
    string _inputFileName, _outputFileName;
    verbosityutils::Verbosity _verbosity = verbosityutils::kInfo;
    typeutils::histtypes _histtype;

    template <typename T>
    void setVar(T a, T b, T& _x, T& _y, string name);
  public:
    bool passVerbosityCheck(verbosityutils::Verbosity level) {
      return verbosityutils::passVerbosityCheck(level, _verbosity);
    }
    void printLog(string message, verbosityutils::Verbosity messageVerbLevel) {
      if (passVerbosityCheck(messageVerbLevel))
        verbosityutils::printLog(message, messageVerbLevel);
    }

    // Getters and setters
    double getEtaMin() { return _etamin; }
    double getEtaMax() { return _etamax; }
    array<double, 2> getEtaRange() { return std::array<double, 2>{_etamin, _etamax}; }
    typeutils::histtypes getHistType() { return _histtype; }
    string getInputFileName() { return _inputFileName; }
    string getOutputFileName() { return _outputFileName; }
    double getPtJetMin() { return _ptjetmin; }
    double getPtJetMax() { return _ptjetmax; }
    array<double, 2> getPtJetRange() { return std::array<double, 2>{_ptjetmin, _ptjetmax}; }
    int getTrain() { return _train; }
    verbosityutils::Verbosity getVerbosity() { return _verbosity; }

    void setEta(double a, double b) { setVar(a, b, _etamin, _etamax, "setEta()"); }
    void setEta(array<double, 2> x) { setEta(x[0], x[1]); }
    void setHistType(typeutils::histtypes h) { _histtype = h; }
    void setInputFileName(string s) { _inputFileName = s; }
    void setOutputFileName(string s) { _outputFileName = s; }
    void setPtJet(double a, double b) { setVar(a, b, _ptjetmin, _ptjetmax, "setPtJet()"); }
    void setPtJet(array<double, 2> x) { setPtJet(x[0], x[1]); }
    void setTrain(int t) { _train = t; }
    void setVerbosity(verbosityutils::Verbosity v) {  _verbosity = v; }

    // Utilities
    double GetEvtXsec();
    TFile* GetFile();
    template <typename T> T* GetHist(string name);
    template <typename T> T* GetHist(typeutils::histtypes htype);
    string GetHistName(typeutils::histtypes htype);
    double GetNevts();
    double GetNjets(double ptmin, double ptmax);
    void SetInputFileNameFromTrain();
    string GetNameFromPtJet(string prefix, string suffix);
};

TFile* InputSettings::GetFile() {
  TFile* f = TFile::Open(_inputFileName.c_str());
  if (!f) {
    printLog("InputSettings::GetFile() Could not open file " + _inputFileName, verbosityutils::kErrors);
    return nullptr;
  }
  return f;
}

template <typename T>
T* InputSettings::GetHist(string name) {
  TFile* file = GetFile();
  if (!file)
    return nullptr;

  T* h = (T*)file->Get(name.c_str());
  if (!h) {
    printLog(TString::Format("InputSetting::GetHist() Could not find histogram %s in file %s", name.c_str(), _inputFileName.c_str()).Data(), verbosityutils::kErrors);
    return nullptr;
  }
  return h;
}

template <typename T>
T* InputSettings::GetHist(typeutils::histtypes htype) {
  T* h = GetHist<T>(GetHistName(htype));
  return h;
}

string InputSettings::GetHistName(typeutils::histtypes htype) {
  string n = "jet-v0qa/";

  switch (_train) {
    case 745299:
      if (htype == typeutils::kNevts)
        n = "jet-fragmentation_id38322/";
      if (htype == typeutils::kNjets)
        n = "jet-fragmentation_id38322/";
  }

  switch (htype) {
    case typeutils::kNevts:
      n += "matching/hEvents";
      break;
    case typeutils::kNjets:
      n += "mcd/jets/inclDetJetPtEtaPhi"; //FIXME: is this the right jet collection?
      break;
    case typeutils::kInclusiveV0s:
      n += "collisions/K0SPtEtaMass";
      break;
    case typeutils::kInclusiveV0sWrongCollision:
      n += "collisions/K0SPtEtaMassWrongColl";
      break;
    case typeutils::kInJetsV0Pt:
      n += "collisions/JetPtEtaK0SPtMass";
      break;
    case typeutils::kInJetsV0PtWrongCollision:
      n += "collisions/JetPtEtaK0SPtMassWrongColl";
      break;
    case typeutils::kInJetsV0Z:
      n += "collisions/JetPtEtaK0SFragMass";
      break;
    case typeutils::kInJetsV0ZWrongCollision:
      n += "collisions/JetPtEtaK0SFragMassWrongColl";
      break;
    case typeutils::kInMatchedJetsV0Pt:
      n += "collisions/JetsPtEtaK0SPt";
      break;
    case typeutils::kInMatchedJetsV0PtWrongCollision:
      n += "collisions/JetsPtEtaK0SPtWrongColl";
      break;
    case typeutils::kInMatchedJetsV0Z:
      n += "collisions/JetsPtEtaK0SZ";
      break;
    case typeutils::kInMatchedJetsV0ZWrongCollision:
      n += "collisions/JetsPtEtaK0SZWrongColl";
      break;
    default:
      printLog("InputSettings::getHistName() invalid histtype", verbosityutils::kErrors);
  }
  return n;
}


void InputSettings::SetInputFileNameFromTrain() {
  _inputFileName = "~/cernbox/TrainOutput/" + to_string(_train) + "/AnalysisResults.root";
}

string InputSettings::GetNameFromPtJet(string prefix, string suffix) {
  return TString::Format("%s_ptjet%.f-%.f%s", prefix.c_str(), _ptjetmin, _ptjetmax, suffix.c_str()).Data();
}

double InputSettings::GetNevts() {
  TH1D* h = GetHist<TH1D>(typeutils::kNevts);
  if (!h)
    return -1;

  return h->GetBinContent(2); // All reconstructed events
}

double InputSettings::GetEvtXsec() {
  TH1D* h = GetHist<TH1D>(typeutils::kNevts);
  if (!h)
    return -1;

  return h->GetBinContent(3);
}
double InputSettings::GetNjets(double ptmin, double ptmax) {
  TH3D* h = GetHist<TH3D>(typeutils::kNjets);
  if (!h)
    return -1.;

  array<int, 2> bins = histutils::getProjectionBins(h->GetYaxis(), _etamin, _etamax);
  TH1D* hpt = h->ProjectionX("hpt", bins[0], bins[1], 0, 1 + h->GetNbinsZ());

  bins = histutils::getProjectionBins(hpt->GetXaxis(), ptmin, ptmax);
  return hpt->Integral(bins[0], bins[1]); // FIXME: Should this have option width???
}

template <typename T>
void InputSettings::setVar(T a, T b, T& _x, T& _y, string name) {
  if (a > b) {
    printLog(TString::Format("%s: min > max", name.c_str()).Data(), verbosityutils::kErrors);
    return;
  }
  _x = a;
  _y = b;
}

// ------------------------------------------------------------------------------------
//
// Plot spectra of inclusive V0s matched with the wrong collision
//
// ------------------------------------------------------------------------------------

array<TH1D*, 2> getHistsInclusive(InputSettings& inputs) {
  inputs.printLog(TString::Format("getHistsInclusive() Getting ptV0 histograms for eta in [%.2f, %.2f]", inputs.getEtaMin(), inputs.getEtaMax()).Data(), verbosityutils::kInfo);
  array<TH1D*, 2> errResult = { nullptr, nullptr };

  TFile* f = inputs.GetFile();
  if (!f)
    return errResult;

  TH3D* h3All = inputs.GetHist<TH3D>(typeutils::kInclusiveV0s);
  TH3D* h3WrongColl = inputs.GetHist<TH3D>(typeutils::kInclusiveV0sWrongCollision);
  if (!h3All || !h3WrongColl)
    return errResult;

  array<int, 2> etaBins = histutils::getProjectionBins(h3All->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hAll = h3All->ProjectionX("hAll", etaBins[0], etaBins[1], 0, h3All->GetNbinsZ() + 1);

  etaBins  = histutils::getProjectionBins(h3WrongColl->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hWrongColl = (TH1D*)h3WrongColl->ProjectionX("hWrongColl", etaBins[0], etaBins[1], 0, h3WrongColl->GetNbinsZ()+1);

  hAll = (TH1D*)histutils::rebinHist(hAll, histutils::rebinnedV0PtHist("K0S", "hAllRebinned"));
  hWrongColl = (TH1D*)histutils::rebinHist(hWrongColl, histutils::rebinnedV0PtHist("K0S", "hWrongCollRebinned"));
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    inputs.printLog("getHistsInclusive() Printing histograms", verbosityutils::kDebug);
    hAll->Print("all");
    hWrongColl->Print("all");
  }

  double nEvts = inputs.GetNevts();
  if (nEvts <= 0)
    return errResult;

  hAll->Scale(1. / nEvts, "width");
  hWrongColl->Scale(1. / nEvts, "width");

  return array<TH1D*, 2>{hAll, hWrongColl};
}

void plotInclusivePt(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.75, 0.75);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = getHistsInclusive(inputs);
  TH1D* hInclusive = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p("pileUp.pdf", true, 0.04);
  string xTitle = mystrings::sPtV0WithUnits;
  string yTitle = mystrings::sV0PtPerXsec;
  p.makeFrame(0., inputs.getPtJetMax(), 1e-10, 0.1, xTitle, yTitle);

  p.setHists({hInclusive, hPileUp});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInclusive, "Inclusive V0s");
  p.addLegendEntry(hPileUp, "V0 pile-up");

  double xLatex = 0.45, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sEtaV0Range075);
  p.plot();
}

void plotInclusivePtRatio(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.75, 0.75);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = getHistsInclusive(inputs);
  TH1D* hInclusive = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p("pileUp_ratio.pdf", true, 0.04);
  string xTitle = mystrings::sPtV0WithUnits;
  string yTitle = mystrings::sRatio;
  p.makeFrame(0., inputs.getPtJetMax(), 1e-3, 2., xTitle, yTitle);

  p.setHists({hInclusive, hPileUp});
  p.setHistStyles();
  p.makeRatios();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInclusive, "Inclusive V0s");
  p.addLegendEntry(hPileUp, "V0 pile-up");

  double xLatex = 0.25, yLatex = 0.70;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sEtaV0Range075);
  p.plot();
}

void plotincl() {
  int train = 745299;
  plotInclusivePt(train);
  plotInclusivePtRatio(train);
}

// ------------------------------------------------------------------------------------
//
// Plot pt spectra of in-jet V0s matched with the wrong collision
//
// ------------------------------------------------------------------------------------

array<TH1D*, 2> getPtHistsInJets(InputSettings& inputs) {
  inputs.printLog(TString::Format("getPtHistsInJets() Getting ptV0 histograms for jet pT in [%.f, %.f] GeV/c and eta in [%.2f, %.2f]", inputs.getPtJetMin(), inputs.getPtJetMax(), inputs.getEtaMin(), inputs.getEtaMax()).Data(), verbosityutils::kInfo);
  array<TH1D*, 2> errResult = { nullptr, nullptr };

  TFile* f = inputs.GetFile();
  if (!f)
    return errResult;

  const int axisJetPt = 0;
  const int axisJetEta = 1;
  const int axisV0Pt = 2;
  const int axisV0Mass = 3;
  THnSparse* hnAll = inputs.GetHist<THnSparse>(typeutils::kInJetsV0Pt);
  THnSparse* hnWrongColl = inputs.GetHist<THnSparse>(typeutils::kInJetsV0PtWrongCollision);
  if (!hnAll || !hnWrongColl)
    return errResult;

  array<int, 2> ptBins = histutils::getProjectionBins(hnAll->GetAxis(axisJetPt), inputs.getPtJetMin(), inputs.getPtJetMax());
  array<int, 2> etaBins = histutils::getProjectionBins(hnAll->GetAxis(axisJetEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnAll->GetAxis(axisJetPt)->SetRange(ptBins[0], ptBins[1]);
  hnAll->GetAxis(axisJetEta)->SetRange(etaBins[0], etaBins[1]);
  TH1D* hAll = (TH1D*)hnAll->Projection(axisV0Pt);
  hAll->SetName("hAll");

  ptBins = histutils::getProjectionBins(hnWrongColl->GetAxis(axisJetPt), inputs.getPtJetMin(), inputs.getPtJetMax());
  etaBins  = histutils::getProjectionBins(hnWrongColl->GetAxis(axisJetEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnWrongColl->GetAxis(axisJetPt)->SetRange(ptBins[0], ptBins[1]);
  hnWrongColl->GetAxis(axisJetEta)->SetRange(etaBins[0], etaBins[1]);
  TH1D* hWrongColl = (TH1D*)hnWrongColl->Projection(axisV0Pt);
  hWrongColl->SetName("hWrongColl");

  hAll = (TH1D*)histutils::rebinHist(hAll, histutils::rebinnedV0PtHist("K0S", "hAllRebinned"));
  hWrongColl = (TH1D*)histutils::rebinHist(hWrongColl, histutils::rebinnedV0PtHist("K0S", "hWrongCollRebinned"));
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    inputs.printLog("getPtHistsInJets() Printing histograms", verbosityutils::kDebug);
    hAll->Print("all");
    hWrongColl->Print("all");
  }

  double nJets = inputs.GetNjets(inputs.getPtJetMin(), inputs.getPtJetMax());
  if (nJets <= 0)
    return errResult;

  hAll->Scale(1. / nJets, "width");
  hWrongColl->Scale(1. / nJets, "width");

  return array<TH1D*, 2>{hAll, hWrongColl};
}

void plotInJetPt1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = getPtHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp", ".pdf"), true, 0.04);
  string xTitle = mystrings::sPtV0WithUnits;
  string yTitle = mystrings::sV0PtPerJetXsecWithUnits;
  p.makeFrame(0., inputs.getPtJetMax(), 1e-7, 0.1, xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.45, yLatex = 0.875;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetPtRatio1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = getPtHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp", "_ratio.pdf"), true, 0.04);
  string xTitle = mystrings::sPtV0WithUnits;
  string yTitle = mystrings::sRatio;
  p.makeFrame(0., inputs.getPtJetMax(), 1e-3, 2., xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();
  p.makeRatios();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "Inclusive V0s");
  p.addLegendEntry(hPileUp, "V0 pile-up");

  double xLatex = 0.25, yLatex = 0.70;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetPt2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = getPtHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    hInJets->Print("all");
    hPileUp->Print("all");
  }

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp", ".pdf"), true, 0.04);
  string xTitle = mystrings::sPtV0WithUnits;
  string yTitle = mystrings::sV0PtPerJetXsecWithUnits;
  p.makeFrame(0., inputs.getPtJetMax(), 1e-7, 0.1, xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.45, yLatex = 0.875;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetPtRatio2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = getPtHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    hInJets->Print("all");
    hPileUp->Print("all");
  }

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp", "_ratio.pdf"), true, 0.04);
  string xTitle = mystrings::sPtV0WithUnits;
  string yTitle = mystrings::sRatio;
  p.makeFrame(0., inputs.getPtJetMax(), 1e-3, 2., xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();
  p.makeRatios();
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    hInJets->Print("all");
    hPileUp->Print("all");
  }

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.25, yLatex = 0.70;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetPt3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = getPtHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    hInJets->Print("all");
    hPileUp->Print("all");
  }

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp", ".pdf"), true, 0.04);
  string xTitle = mystrings::sPtV0WithUnits;
  string yTitle = mystrings::sV0PtPerJetXsecWithUnits;
  p.makeFrame(0., inputs.getPtJetMax(), 1e-7, 0.1, xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.45, yLatex = 0.875;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetPtRatio3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = getPtHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    inputs.printLog("plotInJetPtRatio3040() Printing histograms", verbosityutils::kDebug);
    hInJets->Print("all");
    hPileUp->Print("all");
  }

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp", "_ratio.pdf"), true, 0.04);
  string xTitle = mystrings::sPtV0WithUnits;
  string yTitle = mystrings::sRatio;
  p.makeFrame(0., inputs.getPtJetMax(), 1e-3, 2., xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();
  p.makeRatios();
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    inputs.printLog("plotInJetPtRatio3040() Printing histograms after taking ratios", verbosityutils::kDebug);
    for (auto h : p.getHists()) {
      h->Print("all");
    }
  }

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.25, yLatex = 0.70;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotinjetpt() {
  int train = 745299;
  plotInJetPt1020(train);
  plotInJetPtRatio1020(train);
  plotInJetPt2030(train);
  plotInJetPtRatio2030(train);
  plotInJetPt3040(train);
  plotInJetPtRatio3040(train);
}

// ------------------------------------------------------------------------------------
//
// Plot z spectra of in-jet V0s matched with the wrong collision
//
// ------------------------------------------------------------------------------------

array<TH1D*, 2> getZHistsInJets(InputSettings& inputs) {
  inputs.printLog(TString::Format("getZHistsInJets() Getting zV0 histograms for jet pT in [%.f, %.f] GeV/c and eta in [%.2f, %.2f]", inputs.getPtJetMin(), inputs.getPtJetMax(), inputs.getEtaMin(), inputs.getEtaMax()).Data(), verbosityutils::kInfo);
  array<TH1D*, 2> errResult = { nullptr, nullptr };

  TFile* f = inputs.GetFile();
  if (!f)
    return errResult;

  const int axisJetPt = 0;
  const int axisJetEta = 1;
  const int axisV0Z = 2;
  const int axisV0Mass = 3;
  THnSparse* hnAll = inputs.GetHist<THnSparse>(typeutils::kInJetsV0Z);
  THnSparse* hnWrongColl = inputs.GetHist<THnSparse>(typeutils::kInJetsV0ZWrongCollision);
  if (!hnAll || !hnWrongColl)
    return errResult;

  array<int, 2> ptBins = histutils::getProjectionBins(hnAll->GetAxis(axisJetPt), inputs.getPtJetMin(), inputs.getPtJetMax());
  array<int, 2> etaBins = histutils::getProjectionBins(hnAll->GetAxis(axisJetEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnAll->GetAxis(axisJetPt)->SetRange(ptBins[0], ptBins[1]);
  hnAll->GetAxis(axisJetEta)->SetRange(etaBins[0], etaBins[1]);
  TH1D* hAll = (TH1D*)hnAll->Projection(axisV0Z);
  hAll->SetName("hAll");

  ptBins = histutils::getProjectionBins(hnWrongColl->GetAxis(axisJetPt), inputs.getPtJetMin(), inputs.getPtJetMax());
  etaBins  = histutils::getProjectionBins(hnWrongColl->GetAxis(axisJetEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnWrongColl->GetAxis(axisJetPt)->SetRange(ptBins[0], ptBins[1]);
  hnWrongColl->GetAxis(axisJetEta)->SetRange(etaBins[0], etaBins[1]);
  TH1D* hWrongColl = (TH1D*)hnWrongColl->Projection(axisV0Z);
  hWrongColl->SetName("hWrongColl");

  hAll = (TH1D*)histutils::rebinHist(hAll, histutils::rebinnedV0ZHist("hAllRebinned"));
  hWrongColl = (TH1D*)histutils::rebinHist(hWrongColl, histutils::rebinnedV0ZHist("hWrongCollRebinned"));
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    inputs.printLog("getZHistsInJets() Printing histograms", verbosityutils::kDebug);
    hAll->Print("all");
    hWrongColl->Print("all");
  }

  double nJets = inputs.GetNjets(inputs.getPtJetMin(), inputs.getPtJetMax());
  if (nJets <= 0)
    return errResult;

  hAll->Scale(1. / nJets, "width");
  hWrongColl->Scale(1. / nJets, "width");

  return array<TH1D*, 2>{hAll, hWrongColl};
}

void plotInJetZ1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = getZHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp_z", ".pdf"), true, 0.04);
  string xTitle = mystrings::sZV0;
  string yTitle = mystrings::sV0ZPerJetXsec;
  p.makeFrame(1e-3, 1 + 1e-3, 1e-5, 1., xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.45, yLatex = 0.875;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetZRatio1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = getZHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp_z", "_ratio.pdf"), true, 0.04);
  string xTitle = mystrings::sZV0;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1 + 1e-3, 1e-3, 2., xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();
  p.makeRatios();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.25, yLatex = 0.70;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetZ2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = getZHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp_z", ".pdf"), true, 0.04);
  string xTitle = mystrings::sZV0;
  string yTitle = mystrings::sV0ZPerJetXsec;
  p.makeFrame(1e-3, 1 + 1e-3, 1e-5, 1., xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.45, yLatex = 0.875;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetZRatio2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = getZHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp_z", "_ratio.pdf"), true, 0.04);
  string xTitle = mystrings::sZV0;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1 + 1e-3, 1e-3, 2., xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();
  p.makeRatios();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.25, yLatex = 0.70;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetZ3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = getZHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp_z", ".pdf"), true, 0.04);
  string xTitle = mystrings::sZV0;
  string yTitle = mystrings::sV0ZPerJetXsec;
  p.makeFrame(1e-3, 1 + 1e-3, 1e-5, 1., xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.45, yLatex = 0.875;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotInJetZRatio3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = getZHistsInJets(inputs);
  TH1D* hInJets = hists[0];
  TH1D* hPileUp = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pileUp_z", "_ratio.pdf"), true, 0.04);
  string xTitle = mystrings::sZV0;
  string yTitle = mystrings::sRatio;
  p.makeFrame(1e-3, 1 + 1e-3, 1e-3, 2., xTitle, yTitle);

  p.setHists({hInJets, hPileUp});
  p.setHistStyles();
  p.makeRatios();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jets");
  p.addLegendEntry(hPileUp, "V0 pile-up in jets");

  double xLatex = 0.25, yLatex = 0.70;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets + ", " + mystrings::sJetR04Eta035);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotinjetz() {
  int train = 745299;
  // plotInJetZ1020(train);
  // plotInJetZRatio1020(train);
  // plotInJetZ2030(train);
  // plotInJetZRatio2030(train);
  plotInJetZ3040(train);
  plotInJetZRatio3040(train);
}

#endif