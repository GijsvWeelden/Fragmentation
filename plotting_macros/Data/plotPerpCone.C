
#include "../histUtils.C"
#include "../plotUtils.C"
#include "../myStrings.C"

#ifndef __PLOTPERPCONE_H__
#define __PLOTPERPCONE_H__

namespace verbosityutils {
  enum Verbosity {kErrors, kWarnings, kInfo, kDebug, kDebugMax};
  // Name of the verbosity as written in code
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
  enum histtypes { kJets, kPerpConeV0s, kInJetV0s };
} // namespace typeutils

struct InputSettings {
  private:
    const int _conesPerJet = 2;
    int _train;
    double _etamin, _etamax;
    double _ptjetmin, _ptjetmax;
    string _inputFileName, _outputFileName;
    verbosityutils::Verbosity _verbosity = verbosityutils::kInfo;
    typeutils::histtypes _histtype;
    // hadron, histname, ptBinEdges

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
    const int getConesPerJet() { return _conesPerJet; }
    verbosityutils::Verbosity getVerbosity() { return _verbosity; }

    void setEta(double a, double b) { setVar(a, b, _etamin, _etamax, "setEta()"); }
    void setEta(array<double, 2> x) { setEta(x[0], x[1]); }
    void setEtaMin(double x) { _etamin = x; }
    void setEtaMax(double x) { _etamax = x; }
    void setHistType(typeutils::histtypes h) { _histtype = h; }
    void setInputFileName(string s) { _inputFileName = s; }
    void setOutputFileName(string s) { _outputFileName = s; }
    void setPtJet(double a, double b) { setVar(a, b, _ptjetmin, _ptjetmax, "setPtJet()"); }
    void setPtJet(array<double, 2> x) { setPtJet(x[0], x[1]); }
    void setPtJetMin(double x) { _ptjetmin = x; }
    void setPtJetMax(double x) { _ptjetmax = x; }
    void setTrain(int t) { _train = t; }
    void setVerbosity(verbosityutils::Verbosity v) {  _verbosity = v; }

    // Utilities
    TFile* GetFile();
    template <typename T> T* GetHist(string name);
    template <typename T> T* GetHist(typeutils::histtypes htype);
    string GetHistName(typeutils::histtypes htype);
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
  string n = "jet-fragmentation/data/";

  switch (htype) {
    case typeutils::kJets:
      n += "jets/inclJetPtEtaPhi";
      break;
    case typeutils::kPerpConeV0s:
      n += "PC/JetPtEtaK0SPt";
      break;
    case typeutils::kInJetV0s:
      n += "jets/V0/jetPtK0SPtMass";
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

double InputSettings::GetNjets(double ptmin, double ptmax) {
  TH3D* h = GetHist<TH3D>(typeutils::kJets);
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
// Plot spectra of V0s in the jet cone and the perp cone
//
// ------------------------------------------------------------------------------------

array<TH1D*, 2> getHistsPerpConePt(InputSettings& inputs) {
  inputs.printLog(TString::Format("getHistsPerpConePt() Getting ptV0 histograms for jet pT in [%.f, %.f] GeV/c and eta in [%.2f, %.2f]", inputs.getPtJetMin(), inputs.getPtJetMax(), inputs.getEtaMin(), inputs.getEtaMax()).Data(), verbosityutils::kInfo);
  array<TH1D*, 2> errResult = { nullptr, nullptr };

  TFile* f = inputs.GetFile();
  if (!f)
    return errResult;

  TH3D* h3PerpCone = inputs.GetHist<TH3D>(typeutils::kPerpConeV0s);
  TH3D* h3InJets = inputs.GetHist<TH3D>(typeutils::kInJetV0s); // Does not contain eta axis, but eta cut is applied in task
  if (!h3PerpCone || !h3InJets)
    return errResult;

  array<int, 2> ptBins = histutils::getProjectionBins(h3PerpCone->GetXaxis(), inputs.getPtJetMin(), inputs.getPtJetMax());
  array<int, 2> etaBins = histutils::getProjectionBins(h3PerpCone->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hPerpCone = h3PerpCone->ProjectionZ("hPerpCone", ptBins[0], ptBins[1], etaBins[0], etaBins[1]);

  ptBins  = histutils::getProjectionBins(h3InJets->GetXaxis(), inputs.getPtJetMin(), inputs.getPtJetMax());
  TH1D* hInJets = (TH1D*)h3InJets->ProjectionY("hInJets", ptBins[0], ptBins[1], 0, h3InJets->GetNbinsZ()+1);

  hPerpCone = (TH1D*)histutils::rebinHist(hPerpCone, histutils::rebinnedV0PtHist("K0S", "hPerpConeRebinned"));
  hInJets = (TH1D*)histutils::rebinHist(hInJets, histutils::rebinnedV0PtHist("K0S", "hInJetsRebinned"));
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    hPerpCone->Print("all");
    hInJets->Print("all");
  }

  double nJets = inputs.GetNjets(inputs.getPtJetMin(), inputs.getPtJetMax());
  if (nJets < 0)
    return errResult;

  hPerpCone->Scale(1. / nJets, "width");
  hInJets->Scale(1. / nJets, "width");

  return array<TH1D*, 2>{hPerpCone, hInJets};
}

void plotPerpConePtFromTrain1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = getHistsPerpConePt(inputs);
  TH1D* hPerpCone = hists[0];
  TH1D* hInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pc", ".pdf"), true, 0.04);
  string xTitle = mystrings::sPtK0SWithUnits;
  string sNJetsCones = mystrings::addSubscript(mystrings::sNumber, "jets, cones");
  string yTitle = mystrings::getOneOverString(sNJetsCones) + mystrings::getdYdXString(mystrings::sNK0S, mystrings::sPtK0S);
  yTitle = mystrings::addUnits(yTitle, mystrings::sCGeV, true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-7, 10, xTitle, yTitle);

  p.setHists({hInJets, hPerpCone});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jet cone");
  p.addLegendEntry(hPerpCone, "V0s in perp. cone");

  double xLatex = 0.25, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotPerpConePtRatioFromTrain1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = getHistsPerpConePt(inputs);
  TH1D* hPerpCone = hists[0];
  TH1D* hInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pc", "_ratio.pdf"), true, 0.04);
  string xTitle = mystrings::sPtK0SWithUnits;
  string sNJetsCones = mystrings::addSubscript(mystrings::sNumber, "jets, cones");
  string yTitle = mystrings::getOneOverString(sNJetsCones) + mystrings::getdYdXString(mystrings::sNK0S, mystrings::sPtK0S);
  yTitle = mystrings::addUnits(yTitle, mystrings::sCGeV, true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-3, 2., mystrings::sPtK0SWithUnits, mystrings::sRatio);

  p.setHists({hInJets, hPerpCone});
  p.setHistStyles();
  p.makeRatios();

  double xLegend = 0.25, yLegend = 0.20;
  p.makeLegend(xLegend, xLegend + 0.25, yLegend, yLegend + 0.10, "");
  p.addLegendEntry(hInJets, "V0s in jet cone");
  p.addLegendEntry(hPerpCone, "V0s in perp. cone");

  double xLatex = 0.40, yLatex = 0.75;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
} 

void plotPerpConePtFromTrain2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = getHistsPerpConePt(inputs);
  TH1D* hPerpCone = hists[0];
  TH1D* hInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pc", ".pdf"), true, 0.04);
  string xTitle = mystrings::sPtK0SWithUnits;
  string sNJetsCones = mystrings::addSubscript(mystrings::sNumber, "jets, cones");
  string yTitle = mystrings::getOneOverString(sNJetsCones) + mystrings::getdYdXString(mystrings::sNK0S, mystrings::sPtK0S);
  yTitle = mystrings::addUnits(yTitle, mystrings::sCGeV, true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-7, 10, xTitle, yTitle);

  p.setHists({hInJets, hPerpCone});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jet cone");
  p.addLegendEntry(hPerpCone, "V0s in perp. cone");

  double xLatex = 0.25, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotPerpConePtRatioFromTrain2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = getHistsPerpConePt(inputs);
  TH1D* hPerpCone = hists[0];
  TH1D* hInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pc", "_ratio.pdf"), true, 0.04);
  string xTitle = mystrings::sPtK0SWithUnits;
  string sNJetsCones = mystrings::addSubscript(mystrings::sNumber, "jets, cones");
  string yTitle = mystrings::getOneOverString(sNJetsCones) + mystrings::getdYdXString(mystrings::sNK0S, mystrings::sPtK0S);
  yTitle = mystrings::addUnits(yTitle, mystrings::sCGeV, true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-3, 2., mystrings::sPtK0SWithUnits, mystrings::sRatio);

  p.setHists({hInJets, hPerpCone});
  p.setHistStyles();
  p.makeRatios();

  double xLegend = 0.25, yLegend = 0.20;
  p.makeLegend(xLegend, xLegend + 0.25, yLegend, yLegend + 0.10, "");
  p.addLegendEntry(hInJets, "V0s in jet cone");
  p.addLegendEntry(hPerpCone, "V0s in perp. cone");

  double xLatex = 0.40, yLatex = 0.75;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
} 

void plotPerpConePtFromTrain3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = getHistsPerpConePt(inputs);
  TH1D* hPerpCone = hists[0];
  TH1D* hInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pc", ".pdf"), true, 0.04);
  string xTitle = mystrings::sPtK0SWithUnits;
  string sNJetsCones = mystrings::addSubscript(mystrings::sNumber, "jets, cones");
  string yTitle = mystrings::getOneOverString(sNJetsCones) + mystrings::getdYdXString(mystrings::sNK0S, mystrings::sPtK0S);
  yTitle = mystrings::addUnits(yTitle, mystrings::sCGeV, true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-7, 10, xTitle, yTitle);

  p.setHists({hInJets, hPerpCone});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hInJets, "V0s in jet cone");
  p.addLegendEntry(hPerpCone, "V0s in perp. cone");

  double xLatex = 0.25, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
}

void plotPerpConePtRatioFromTrain3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = getHistsPerpConePt(inputs);
  TH1D* hPerpCone = hists[0];
  TH1D* hInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("pc", "_ratio.pdf"), true, 0.04);
  string xTitle = mystrings::sPtK0SWithUnits;
  string sNJetsCones = mystrings::addSubscript(mystrings::sNumber, "jets, cones");
  string yTitle = mystrings::getOneOverString(sNJetsCones) + mystrings::getdYdXString(mystrings::sNK0S, mystrings::sPtK0S);
  yTitle = mystrings::addUnits(yTitle, mystrings::sCGeV, true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-3, 2., mystrings::sPtK0SWithUnits, mystrings::sRatio);

  p.setHists({hInJets, hPerpCone});
  p.setHistStyles();
  p.makeRatios();

  double xLegend = 0.25, yLegend = 0.20;
  p.makeLegend(xLegend, xLegend + 0.25, yLegend, yLegend + 0.10, "");
  p.addLegendEntry(hInJets, "V0s in jet cone");
  p.addLegendEntry(hPerpCone, "V0s in perp. cone");

  double xLatex = 0.40, yLatex = 0.75;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.plot();
} 

void plotinconept() {
  int train = 747926;
  // plotPerpConePtFromTrain1020(train);
  // plotPerpConePtRatioFromTrain1020(train);
  // plotPerpConePtFromTrain2030(train);
  // plotPerpConePtRatioFromTrain2030(train);
  // plotPerpConePtFromTrain3040(train);
  // plotPerpConePtRatioFromTrain3040(train);
}

#endif