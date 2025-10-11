
#include "../histUtils.C"
#include "../plotUtils.C"
#include "../myStrings.C"

#ifndef __PLOTDAUGHTERSHARING_H__
#define __PLOTDAUGHTERSHARING_H__

namespace verbosityutils {
  enum Verbosity {kErrors, kWarnings, kInfo, kDebug, kDebugMax};
  bool is_valid(int v) {
    bool b = (v >= kErrors && v <= kDebugMax);
    if (!b)
      cout << "verbosityutils Error: invalid verbosity level " << v << endl;
    return b;
  }
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
  bool passVerbosityCheck(Verbosity level, Verbosity threshold) { return (is_valid(level) && is_valid(threshold) && level <= threshold); }
  void printLog(string message, Verbosity level, string prefix = "") {
    string s;
    if (level == kErrors || level == kWarnings)
      s = to_string(level) + ": ";

    s += prefix + " " + message;
    cout << s << endl;
  }
} // namespace verbosityutils

namespace histutils {
  enum histtypes {kV0Incl, kV0InclShared, kV0InJet, kV0InJetShared, kJetIncl, kJetV0, kJetShared, kEvents };
  bool is_valid(histtypes h) {
    bool v = (h >= kV0Incl && h <= kEvents);
    if (!v)
      cout << "histutils Error: invalid histtype " << h << endl;
    return v;
  }
} // namespace histutils

struct InputSettings {
  private:
    int _train;
    double _etamin, _etamax;
    double _ptjetmin, _ptjetmax;
    string _inputFileName, _outputFileName;
    verbosityutils::Verbosity _verbosity = verbosityutils::kInfo;
    histutils::histtypes _histtype;

    template <typename T>
    void setVar(T a, T b, T& x, T& y, string name = "setVar()");
  public:
    bool passVerbosityCheck(verbosityutils::Verbosity level) {
      return verbosityutils::passVerbosityCheck(level, _verbosity);
    }
    void printLog(string prefix, string message, verbosityutils::Verbosity messageVerbLevel) {
      if (passVerbosityCheck(messageVerbLevel))
        verbosityutils::printLog(message, messageVerbLevel, prefix);
    }

    // Getters and setters
    double getEtaMin() { return _etamin; }
    double getEtaMax() { return _etamax; }
    array<double, 2> getEtaRange() { return std::array<double, 2>{_etamin, _etamax}; }
    histutils::histtypes getHistType() { return _histtype; }
    string getInputFileName() { return _inputFileName; }
    string getOutputFileName() { return _outputFileName; }
    double getPtJetMin() { return _ptjetmin; }
    double getPtJetMax() { return _ptjetmax; }
    array<double, 2> getPtJetRange() { return std::array<double, 2>{_ptjetmin, _ptjetmax}; }
    int getTrain() { return _train; }
    verbosityutils::Verbosity getVerbosity() { return _verbosity; }

    void setEta(double a, double b) { setVar(a, b, _etamin, _etamax, "setEta()"); }
    void setEta(array<double, 2> x) { setEta(x[0], x[1]); }
    void setEtaMin(double x) { _etamin = x; }
    void setEtaMax(double x) { _etamax = x; }
    void setHistType(histutils::histtypes h) { if (histutils::is_valid(h)) { _histtype = h; } }
    void setInputFileName(string s) { _inputFileName = s; }
    void setOutputFileName(string s) { _outputFileName = s; }
    void setPtJet(double a, double b) { setVar(a, b, _ptjetmin, _ptjetmax, "setPtJet()"); }
    void setPtJet(array<double, 2> x) { setPtJet(x[0], x[1]); }
    void setPtJetMin(double x) { _ptjetmin = x; }
    void setPtJetMax(double x) { _ptjetmax = x; }
    void setTrain(int t) { _train = t; }
    void setVerbosity(verbosityutils::Verbosity v) { if (verbosityutils::is_valid(v)) { _verbosity = v; } }

    // Utilities
    TFile* GetFile();
    template <typename T> T* GetHist(histutils::histtypes htype);
    template <typename T> T* GetHist(string name);
    string GetHistName(histutils::histtypes htype);
    double GetNevts();
    double GetNjets(double ptmin, double ptmax);
    void SetInputFileNameFromTrain();
};

TFile* InputSettings::GetFile() {
  TFile* f = TFile::Open(_inputFileName.c_str());
  if (!f) {
    printLog("GetFile()", "Could not open file " + _inputFileName, verbosityutils::kErrors);
    return nullptr;
  }
  return f;
}

string InputSettings::GetHistName(histutils::histtypes htype) {
  if (!is_valid(htype)) {
    return "";
  }
  string d = "jet-v0qa/sharing/", n = "";

  switch (htype) {
    case histutils::kV0Incl:
      n = "V0PtEtaPhi";
      break;
    case histutils::kV0InclShared:
      n = "V0PtEtaPt";
      break;
    case histutils::kV0InJet:
      n = "JetPtEtaV0Pt";
      break;
    case histutils::kV0InJetShared:
      n = "JetPtEtaV0PtPt";
      break;
    case histutils::kJetIncl:
      n = "JetPtEtaPhi";
      break;
    case histutils::kJetV0:
      n = ""; // Requires summation of histograms in GetJetPtHist()
      break;
    case histutils::kJetShared:
      n = "JetPtEtaPhiShared";
      break;
    case histutils::kEvents:
      n = "hEvents";
      break;
  }
  return d.append(n);
}

template <typename T>
T* InputSettings::GetHist(string name) {
  TFile* file = GetFile();
  if (!file)
    return nullptr;

  T* h = (T*)file->Get(name.c_str());
  if (!h) {
    printLog("GetHist()", "Could not find histogram " + name, verbosityutils::kErrors);
    return nullptr;
  }
  return h;
}

template <typename T>
T* InputSettings::GetHist(histutils::histtypes htype) {
  T* h = nullptr;
  if (htype == histutils::kJetV0) {
    TFile* file = GetFile();
    if (!file)
      return nullptr;

    string dir = GetHistName(htype);
    string nSingle = dir + "JetPtEtaPhiSingle";
    string nMultiple = dir + "JetPtEtaPhiMultiple";
    h = GetHist<T>(nSingle);
    if (!h)
      return nullptr;

    h->SetName("JetPtEtaPhiV0");

    T* g = GetHist<T>(nMultiple);
    if (!g)
      return nullptr;

    h->Add(g);
  } else {
    h = GetHist<T>(GetHistName(htype));
  }
  if (!h) {
    printLog("GetHist()", "Could not find histogram " + GetHistName(htype), verbosityutils::kErrors);
    return nullptr;
  }
  return h;
}

void InputSettings::SetInputFileNameFromTrain() {
  _inputFileName = "~/cernbox/TrainOutput/" + to_string(_train) + "/AnalysisResults.root";
}

double InputSettings::GetNevts() {
  TH1D* h = GetHist<TH1D>(histutils::kEvents);
  if (!h)
    return -1.;

  return h->GetBinContent(1);
}

double InputSettings::GetNjets(double ptmin, double ptmax) {
  TH3D* h = GetHist<TH3D>(histutils::kJetIncl);
  if (!h)
    return -1.;

  array<int, 2> bins = histutils::getProjectionBins(h->GetYaxis(), _etamin, _etamax);
  TH1D* hpt = h->ProjectionX("hpt", bins[0], bins[1], 0, 1 + h->GetNbinsZ());

  bins = histutils::getProjectionBins(hpt->GetXaxis(), ptmin, ptmax);
  return hpt->Integral(bins[0], bins[1]);
}

template <typename T>
void InputSettings::setVar(T a, T b, T& x, T& y, string name) {
  if (a > b) {
    printLog(name, "Error: min > max", verbosityutils::kErrors);
    return;
  }
  x = a;
  y = b;
}

// Get ptjet of jets with shared daughters vs all jets
// Get pt of V0s with shared daughters vs all V0s (in jets)

array<TH1D*, 2> gethistsincl(InputSettings& inputs) {
  THnSparse* hnShared = inputs.GetHist<THnSparse>(histutils::kV0InclShared);
  TH3D* h3All = inputs.GetHist<TH3D>(histutils::kV0Incl);

  array<int, 2> etaBins = histutils::getProjectionBins(h3All->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hAll = h3All->ProjectionX("hAll", etaBins[0], etaBins[1], 0, 1 + h3All->GetNbinsZ());

  etaBins = histutils::getProjectionBins(hnShared->GetAxis(1), inputs.getEtaMin(), inputs.getEtaMax()); // This only cuts on the eta of the trigger V0! Should be fixed in task
  hnShared->GetAxis(1)->SetRange(etaBins[0], etaBins[1]);
  TH1D* hSharedTrigger = hnShared->Projection(0);
  TH1D* hSharedAssoc   = hnShared->Projection(2);
  TH1D* hShared = (TH1D*)hSharedTrigger->Clone("hShared");
  hShared->Add(hSharedAssoc);

  hAll = (TH1D*)histutils::rebinHist(hAll, histutils::rebinnedV0PtHist("K0S", "hAllRebinned"));
  hShared = (TH1D*)histutils::rebinHist(hShared, histutils::rebinnedV0PtHist("K0S", "hSharedRebinned"));
  hAll->Scale(1. / inputs.GetNevts(), "width");
  hShared->Scale(1. / inputs.GetNevts(), "width");

  return array<TH1D*, 2>{hAll, hShared};
}

void plotincl() {
  InputSettings inputs;
  inputs.setTrain(515268); // Does not contain JetPtEtaPhiShared!
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.75, 0.75);

  array<TH1D*, 2> hists = gethistsincl(inputs);
  TH1D* hAll = hists[0];
  TH1D* hShared = hists[1];

  plotutils::Plotter plSpectra("sharing_inclusive.pdf", true, 0.04);
  plSpectra.makeLegend(0.45, 0.70, 0.60, 0.70, "");
  plSpectra.makeFrame(0., 40., 1e-10, 1., mystrings::sPtV0, mystrings::sV0PtPerEvt);

  plSpectra.setHists({hAll, hShared});
  plSpectra.setHistStyles();
  plSpectra.addLegendEntry(hAll, "All V0s");
  plSpectra.addLegendEntry(hShared, "V0s with shared daughters");

  plSpectra.addLatex(0.45, 0.85, "This Thesis");
  plSpectra.addLatex(0.45, 0.80, "ALICE pp data, #sqrt{s} = 13.6 TeV");
  plSpectra.addLatex(0.45, 0.75, "|#eta| < 0.75 Incorrectly applied!");

  plSpectra.plot();

  plotutils::Plotter plRatio("sharing_inclusive_ratio.pdf", true, 0.04);
  plRatio.makeFrame(0., 40., 1e-2, 2., mystrings::sPtV0, "Ratio");
  plRatio.makeLegend(0.50, 0.75, 0.20, 0.35, "");
  plRatio.setHists({hAll, hShared});
  plRatio.setHistStyles();
  plRatio.addLegendEntry(hAll, "All V0s");
  plRatio.addLegendEntry(hShared, "V0s with shared daughters");

  plRatio.addLatex(0.25, 0.75, "This Thesis");
  plRatio.addLatex(0.25, 0.70, "ALICE pp data, #sqrt{s} = 13.6 TeV");
  plRatio.addLatex(0.25, 0.65, "|#eta| < 0.75 Incorrectly applied!");
  plRatio.makeRatios(0);
  plRatio.plot();
}

array<TH1D*, 2> gethistsinjet(InputSettings& inputs) {
  inputs.printLog("gethistsinjet()", TString::Format("Getting histograms for jet pT in [%.f, %.f] GeV/c and eta in [%.2f, %.2f]", inputs.getPtJetMin(), inputs.getPtJetMax(), inputs.getEtaMin(), inputs.getEtaMax()).Data(), verbosityutils::kInfo);

  const int axisJetPt   = 0;
  const int axisJetEta  = 1;
  const int axisTrigger = 2;
  const int axisAssoc   = 3;

  THnSparse* hnShared = inputs.GetHist<THnSparse>(histutils::kV0InJetShared);
  TH3D* h3All = inputs.GetHist<TH3D>(histutils::kV0InJet);

  array<int, 2> ptBins = histutils::getProjectionBins(h3All->GetXaxis(), inputs.getPtJetMin(), inputs.getPtJetMax());
  array<int, 2> etaBins = histutils::getProjectionBins(h3All->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hAllV0sInJets = h3All->ProjectionZ("hAllV0sInJets", ptBins[0], ptBins[1], etaBins[0], etaBins[1]);

  ptBins  = histutils::getProjectionBins(hnShared->GetAxis(axisJetPt), inputs.getPtJetMin(), inputs.getPtJetMax());
  etaBins = histutils::getProjectionBins(hnShared->GetAxis(axisJetEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnShared->GetAxis(axisJetPt)->SetRange(ptBins[0], ptBins[1]);
  hnShared->GetAxis(axisJetEta)->SetRange(etaBins[0], etaBins[1]);
  TH1D* hSharedTrigger = hnShared->Projection(axisTrigger);
  TH1D* hSharedAssoc   = hnShared->Projection(axisAssoc);
  TH1D* hSharedV0sInJets = (TH1D*)hSharedTrigger->Clone("hSharedV0sInJets");

  hAllV0sInJets = (TH1D*)histutils::rebinHist(hAllV0sInJets, histutils::rebinnedV0PtHist("K0S", "hAllV0sInJetsRebinned"));
  hSharedV0sInJets = (TH1D*)histutils::rebinHist(hSharedV0sInJets, histutils::rebinnedV0PtHist("K0S", "hSharedV0sInJetsRebinned"));
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    hAllV0sInJets->Print("all");
    hSharedV0sInJets->Print("all");
  }

  double nJets = inputs.GetNjets(inputs.getPtJetMin(), inputs.getPtJetMax());
  hAllV0sInJets->Scale(1. / nJets, "width");
  hSharedV0sInJets->Scale(1. / nJets, "width");

  return array<TH1D*, 2>{hAllV0sInJets, hSharedV0sInJets};
}

void plotinjet1020() {
  InputSettings inputs;
  inputs.setTrain(515268); // Does not contain JetPtEtaPhiShared!
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = gethistsinjet(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];
  // FIXME: There is a bug in the task, the pTV0 axis of hAllV0sInJets goes up to 10 GeV/c only!
  // hAllV0sInJets->Print("all");
  // hSharedV0sInJets->Print("all");

  plotutils::Plotter plSpectra("sharing_injet.pdf", true, 0.04);
  plSpectra.makeFrame(0., 20., 1e-8, 0.1, mystrings::sPtV0, mystrings::sV0PtPerJet);

  plSpectra.setHists({hAllV0sInJets, hSharedV0sInJets});
  plSpectra.setHistStyles();

  plSpectra.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  plSpectra.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  plSpectra.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.55;
  plSpectra.addLatex(xLatex, yLatex, "This Thesis, ALICE pp data");
  plSpectra.addLatex(xLatex, yLatex - 0.05, "#sqrt{s} = 13.6 TeV");
  plSpectra.addLatex(xLatex, yLatex - 0.10, "Anti-#it{k}_{T} ch+V0 jets");
  plSpectra.addLatex(xLatex, yLatex - 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.35");
  plSpectra.addLatex(xLatex, yLatex - 0.20, TString::Format("%.f < #it{p}_{T,jet} < %.f GeV/#it{c}", inputs.getPtJetMin(), inputs.getPtJetMax()).Data());

  plSpectra.plot();

  plotutils::Plotter plRatio("sharing_injet_ratio.pdf", true, 0.04);
  plRatio.makeFrame(0., 20., 1e-2, 2., mystrings::sPtV0, "Ratio");
  plRatio.makeLegend(0.50, 0.75, 0.20, 0.35, "");
  plRatio.setHists({hSharedV0sInJets, hAllV0sInJets});
  plRatio.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");
  plRatio.addLegendEntry(hAllV0sInJets, "All V0s in jets");

  xLatex = 0.25, yLatex = 0.60;
  plSpectra.addLatex(xLatex, yLatex, "This Thesis, ALICE pp data");
  plSpectra.addLatex(xLatex, yLatex - 0.05, "#sqrt{s} = 13.6 TeV");
  plSpectra.addLatex(xLatex, yLatex - 0.10, "Anti-#it{k}_{T} ch+V0 jets");
  plSpectra.addLatex(xLatex, yLatex - 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.35");
  plSpectra.addLatex(xLatex, yLatex - 0.20, TString::Format("%.f < #it{p}_{T,jet} < %.f GeV/#it{c}", inputs.getPtJetMin(), inputs.getPtJetMax()).Data());

  plRatio.makeRatios(hAllV0sInJets);
  plRatio.plot();
}

array<TH1D*, 2> gethistsjet(InputSettings& inputs) {
  inputs.printLog("gethistjet()", TString::Format("Getting histograms for jet pT in [%.f, %.f] GeV/c and eta in [%.2f, %.2f]", inputs.getPtJetMin(), inputs.getPtJetMax(), inputs.getEtaMin(), inputs.getEtaMax()).Data(), verbosityutils::kInfo);

  TH3D* h3All     = inputs.GetHist<TH3D>(histutils::kJetIncl); // All jets
  TH3D* h3AllwV0s = inputs.GetHist<TH3D>(histutils::kJetV0); // All jets with V0 candidates
  TH3D* h3Shared  = inputs.GetHist<TH3D>(histutils::kJetShared); // All jets with V0s that share daughters

  array<int, 2> etaBins = histutils::getProjectionBins(h3All->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hAll = h3All->ProjectionX("hAll", etaBins[0], etaBins[1], 0, 1 + h3All->GetNbinsZ());

  etaBins = histutils::getProjectionBins(h3AllwV0s->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hAllwV0s = h3AllwV0s->ProjectionX("hAll", etaBins[0], etaBins[1], 0, 1 + h3AllwV0s->GetNbinsZ());

  etaBins = histutils::getProjectionBins(h3Shared->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hShared = h3Shared->ProjectionX("hShared", etaBins[0], etaBins[1], 0, 1 + h3Shared->GetNbinsZ());

  hAll->Scale(1. / inputs.GetNevts(), "width");
  hAllwV0s->Scale(1. / inputs.GetNevts(), "width");
  hShared->Scale(1. / inputs.GetNevts(), "width");

  return array<TH1D*, 2>{hAll, hAllwV0s, hShared};
}

void plotjet() {
  InputSettings inputs;
  inputs.setTrain(515268); // Does not contain JetPtEtaPhiShared!
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);

  array<TH1D*, 2> hists = gethistsjet(inputs);
  TH1D* hAll = hists[0];
  TH1D* hAllwV0s = hists[1];
  TH1D* hShared = hists[2];

  plotutils::Plotter plSpectra("sharing_jets.pdf", false, 0.04);
  // Which comparison do you want to make?
  plSpectra.setHists({hAll, hShared});
  // plSpectra.setHists({hAllwV0s, hShared});
  plSpectra.setHistStyles();

  plSpectra.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  plSpectra.addLegendEntry(hAll, "Inclusive jets");
  plSpectra.addLegendEntry(hAllwV0s, "Jets with V0s");
  plSpectra.addLegendEntry(hShared, "Jets with V0s that share daughters");

  double xLatex = 0.25, yLatex = 0.55;
  plSpectra.addLatex(xLatex, yLatex, "This Thesis, ALICE pp data");
  plSpectra.addLatex(xLatex, yLatex - 0.05, "#sqrt{s} = 13.6 TeV");
  plSpectra.addLatex(xLatex, yLatex - 0.10, "Anti-#it{k}_{T} ch+V0 jets");
  plSpectra.addLatex(xLatex, yLatex - 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.35");

  plSpectra.makeFrame(mystrings::sPtJet, mystrings::sJetsPerEvent);
  plSpectra.plot();

  plotutils::Plotter plRatio("sharing_jets_ratio.pdf", false, 0.04);
  plRatio.makeLegend(0.50, 0.75, 0.20, 0.35, "");

  plRatio.setHists(plSpectra.getHists()); // Should ensure the histograms are the same
  plRatio.setHistStyles();
  plRatio.makeRatios(0);

  plSpectra.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  plSpectra.addLegendEntry(hAll, "Inclusive jets");
  plSpectra.addLegendEntry(hAllwV0s, "Jets with V0s");
  plSpectra.addLegendEntry(hShared, "Jets with V0s that share daughters");

  double xLatex = 0.25, yLatex = 0.55;
  plSpectra.addLatex(xLatex, yLatex, "This Thesis, ALICE pp data");
  plSpectra.addLatex(xLatex, yLatex - 0.05, "#sqrt{s} = 13.6 TeV");
  plSpectra.addLatex(xLatex, yLatex - 0.10, "Anti-#it{k}_{T} ch+V0 jets");
  plSpectra.addLatex(xLatex, yLatex - 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.35");

  plRatio.makeFrame(mystrings::sPtJet, "Ratio");
  plRatio.plot();
}

#endif