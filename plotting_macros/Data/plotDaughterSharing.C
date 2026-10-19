
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
  void printLog(string message, Verbosity level) {
    string s;
    if (level == kErrors || level == kWarnings)
      s = to_string(level) + ": ";

    s += message;
    cout << s << endl;
  }
} // namespace verbosityutils

namespace typeutils {
  enum histtypes {kV0InclPt, kV0InclPtShared, kV0InclPtSharedDaughter, kV0PtInJet, kV0PtInJetShared, kV0ZInJet, kV0ZInJetShared, kJetPtIncl, kJetPtWithV0s, kJetPtShared, kEvents };
  bool is_valid(histtypes h) {
    bool v = (h >= kV0InclPt && h <= kEvents);
    if (!v)
      cout << "typeutils Error: invalid histtype " << h << endl;
    return v;
  }
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
    void setVar(T a, T b, T& x, T& y, string name = "setVar()");
  public:
    bool passVerbosityCheck(verbosityutils::Verbosity level) {
      return verbosityutils::passVerbosityCheck(level, _verbosity);
    }
    bool printLog(string message, verbosityutils::Verbosity messageVerbLevel) {
      if (passVerbosityCheck(messageVerbLevel))
        verbosityutils::printLog(message, messageVerbLevel);
      
      return passVerbosityCheck(messageVerbLevel);
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
    void setEtaMin(double x) { _etamin = x; }
    void setEtaMax(double x) { _etamax = x; }
    void setHistType(typeutils::histtypes h) { if (typeutils::is_valid(h)) { _histtype = h; } }
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
    template <typename T> T* GetHist(typeutils::histtypes htype);
    template <typename T> T* GetHist(string name);
    string GetHistName(typeutils::histtypes htype);
    double GetNevts();
    double GetNjets(double ptmin, double ptmax);
    void SetInputFileNameFromTrain();
    string GetNameFromPtJet(string prefix, string suffix);
};

TFile* InputSettings::GetFile() {
  TFile* f = TFile::Open(_inputFileName.c_str());
  if (!f) {
    printLog("GetFile() Could not open file " + _inputFileName, verbosityutils::kErrors);
    return nullptr;
  }
  return f;
}

string InputSettings::GetHistName(typeutils::histtypes htype) {
  if (!is_valid(htype)) {
    return "";
  }
  string d = "jet-v0qa", n = "";

  switch (_train) {
    case 538318:
      d += "_id38323"; // Different directory structure in this train
      break;
    default:
      break;
  }
  d += "/sharing/";

  switch (htype) {
    case typeutils::kV0InclPt:
      n = "V0PtEtaPhi";
      break;
    case typeutils::kV0InclPtShared:
      n = "V0PtEtaPt";
      break;
    case typeutils::kV0InclPtSharedDaughter:
      n = "V0PtEtaPtDaughterPt";
      break;
    case typeutils::kV0PtInJet:
      n = "JetPtEtaV0Pt";
      break;
    case typeutils::kV0PtInJetShared:
      n = "JetPtEtaV0PtPt";
      break;
    case typeutils::kV0ZInJet:
      n = "JetPtEtaV0Z";
      break;
    case typeutils::kV0ZInJetShared:
      n = "JetPtEtaV0ZZ";
      break;
    case typeutils::kJetPtIncl:
      n = "JetPtEtaPhi";
      break;
    case typeutils::kJetPtWithV0s:
      n = ""; // Requires summation of histograms in GetJetPtHist()
      break;
    case typeutils::kJetPtShared:
      n = "JetPtEtaPhiShared";
      break;
    case typeutils::kEvents:
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
    printLog("GetHist() Could not find histogram " + name, verbosityutils::kErrors);
    return nullptr;
  }
  return h;
}

template <typename T>
T* InputSettings::GetHist(typeutils::histtypes htype) {
  T* h = nullptr;
  if (htype == typeutils::kJetPtWithV0s) {
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
    // h->Print();

    T* g = GetHist<T>(nMultiple);
    if (!g)
      return nullptr;
    // g->Print();

    h->Add(g);
    // h->Print();
  } else {
    h = GetHist<T>(GetHistName(htype));
  }
  if (!h) {
    printLog("GetHist() Could not find histogram " + GetHistName(htype), verbosityutils::kErrors);
    return nullptr;
  }
  return h;
}

void InputSettings::SetInputFileNameFromTrain() {
  _inputFileName = "~/cernbox/TrainOutput/" + to_string(_train) + "/AnalysisResults.root";
}

string InputSettings::GetNameFromPtJet(string prefix, string suffix) {
  return TString::Format("%s_ptjet%.f-%.f%s", prefix.c_str(), _ptjetmin, _ptjetmax, suffix.c_str()).Data();
}

double InputSettings::GetNevts() {
  TH1D* h = GetHist<TH1D>(typeutils::kEvents);
  if (!h)
    return -1.;

  return h->GetBinContent(1);
}

double InputSettings::GetNjets(double ptmin, double ptmax) {
  TH3D* h = GetHist<TH3D>(typeutils::kJetPtIncl);
  if (!h)
    return -1.;

  array<int, 2> bins = histutils::getProjectionBins(h->GetYaxis(), _etamin, _etamax);
  TH1D* hpt = h->ProjectionX("hpt", bins[0], bins[1], 0, 1 + h->GetNbinsZ());

  bins = histutils::getProjectionBins(hpt->GetXaxis(), ptmin, ptmax);
  return hpt->Integral(bins[0], bins[1]); // FIXME: Should this have option width???
}

template <typename T>
void InputSettings::setVar(T a, T b, T& x, T& y, string name) {
  if (a > b) {
    printLog(name + " min > max", verbosityutils::kErrors);
    return;
  }
  x = a;
  y = b;
}


// ------------------------------------------------------------------------------------
//
// Plot spectra of inclusive V0s that share daughters
//
// ------------------------------------------------------------------------------------

// Get ptjet of jets with shared daughters vs all jets
// Get pt of V0s with shared daughters vs all V0s (in jets)

array<TH1D*, 2> gethistsincl(InputSettings& inputs) {
  THnSparse* hnShared = inputs.GetHist<THnSparse>(typeutils::kV0InclPtShared); // This is saved as a THn, even though D=3
  TH3D* h3All = inputs.GetHist<TH3D>(typeutils::kV0InclPt);

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

void plotInclFromTrain(int train) {
  InputSettings inputs;
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.75, 0.75);

  array<TH1D*, 2> hists = gethistsincl(inputs);
  TH1D* hAll = hists[0];
  TH1D* hShared = hists[1];

  plotutils::Plotter plSpectra("sharing_inclusive.pdf", true);
  double xLegend = 0.50, yLegend = 0.65;
  plSpectra.makeLegend(xLegend, xLegend + 0.15, yLegend, yLegend + 0.10, "");
  plSpectra.makeFrame(0., 40., 1e-12, 1., mystrings::sPtV0WithUnits, mystrings::sV0PtPerEvtWithUnits);

  plSpectra.setHists({hAll, hShared});
  plSpectra.setHistStyles();
  plSpectra.addLegendEntry(hAll, "All V0s");
  plSpectra.addLegendEntry(hShared, "V0s with shared daughters");

  double xLatex = 0.35, yLatex = 0.85;
  plSpectra.addLatex(xLatex, yLatex, mystrings::sThisThesis + ", " + mystrings::sAlicePpData + ", " + mystrings::sSqrtS);
  plSpectra.addLatex(xLatex, yLatex - 0.05, mystrings::sEtaV0Range075);

  plSpectra.plot();
}

void plotInclRatioFromTrain(int train) {
  InputSettings inputs;
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.75, 0.75);

  array<TH1D*, 2> hists = gethistsincl(inputs);
  TH1D* hAll = hists[0];
  TH1D* hShared = hists[1];

  plotutils::Plotter plRatio("sharing_inclusive_ratio.pdf", true);
  plRatio.makeFrame(0., 40., 1e-4, 2., mystrings::sPtV0, mystrings::sRatio);
  double xLegend = 0.30, yLegend = 0.50;
  plRatio.makeLegend(xLegend, xLegend + 0.20, yLegend, yLegend + 0.15, "");
  plRatio.setHists({hAll, hShared});
  plRatio.setHistStyles();
  plRatio.addLegendEntry(hAll, "All V0s");
  plRatio.addLegendEntry(hShared, "V0s with shared daughters");

  double xLatex = 0.30, yLatex = 0.75;
  plRatio.addLatex(xLatex, yLatex, mystrings::sThisThesis + ", " + mystrings::sAlicePpData + ", " + mystrings::sSqrtS);
  plRatio.addLatex(xLatex, yLatex - 0.05, mystrings::sEtaV0Range075);

  plRatio.makeRatios(0);
  plRatio.plot();
}

// 2d plot of the pt of the V0s that share daughters
void plotIncl2dFromTrain(int train) {
  InputSettings inputs;
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.75, 0.75);

  const int axisPtHard = 0;
  const int axisEta    = 1;
  const int axisPtSoft = 2;
  THnSparse* hnShared = inputs.GetHist<THnSparse>(typeutils::kV0InclPtShared);
  array<int, 2> etaBins = histutils::getProjectionBins(hnShared->GetAxis(axisEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnShared->GetAxis(axisEta)->SetRange(etaBins[0], etaBins[1]);
  TH2D* hTrigAssoc = (TH2D*)hnShared->Projection(axisPtSoft, axisPtHard);
  hTrigAssoc = (TH2D*)histutils::rebinHist2D(hTrigAssoc, histutils::rebinnedV0PtHist2D("K0S", "K0S", "hTrigAssocRebinned"));

  // hTrigAssoc->Scale(1. / inputs.GetNevts(), "width");
  histutils::normaliseHistColByCol(hTrigAssoc);

  string sPtHard = mystrings::addSuperscript(mystrings::sPtV0, "hard");
  string sPtSoft = mystrings::addSuperscript(mystrings::sPtV0, "soft");

  plotutils::Plotter p("sharing_inclusive_correlation.pdf", false, false, false);
  p.setDrawOption("colz");
  p.addHistogram(hTrigAssoc);
  p.setZAxisRange(1e-5, 0.75); // Eyeballed values
  p.setFrameTitleOffsetY(plotutils::defaultValues::defTitleOffset);

  string xTitle = sPtHard;
  string yTitle = sPtSoft;
  xTitle = mystrings::addUnits(xTitle, mystrings::sGevC, true);
  yTitle = mystrings::addUnits(yTitle, mystrings::sGevC, true);
  p.makeFrame(0., 40., 0., 40., xTitle, yTitle);
  // FIXME: Label size for Z axis!
  plotutils::SetPadMargins(0.15, 0.10, 0.15, 0.05);

  double xLatex = 0.20, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesis + ", " + mystrings::sAlicePpData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sEtaV0Range075);


  if (inputs.printLog("plotInJetZFromTrain1020() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

// Projection of correlation hist. For each column, project the dist
void plotIncl2dProjectionFromTrain(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.75, 0.75);

  const int axisPtHard = 0;
  const int axisEta    = 1;
  const int axisPtSoft = 2;
  THnSparse* hnShared = inputs.GetHist<THnSparse>(typeutils::kV0InclPtShared);
  array<int, 2> etaBins = histutils::getProjectionBins(hnShared->GetAxis(axisEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnShared->GetAxis(axisEta)->SetRange(etaBins[0], etaBins[1]);
  TH2D* hTrigAssoc = (TH2D*)hnShared->Projection(axisPtSoft, axisPtHard);
  hTrigAssoc = (TH2D*)histutils::rebinHist2D(hTrigAssoc, histutils::rebinnedV0PtHist2D("K0S", "K0S", "hTrigAssocRebinned"));

  // hTrigAssoc->Scale(1. / inputs.GetNevts(), "width");
  // histutils::normaliseHistColByCol(hTrigAssoc);

  string sPtHard = mystrings::addSuperscript(mystrings::sPtV0, "hard");
  string sPtSoft = mystrings::addSuperscript(mystrings::sPtV0, "soft");
  string sNSoft  = mystrings::addSuperscript(mystrings::sNV0, "soft");

  plotutils::Plotter p("sharing_inclusive_correlation_projection.pdf", true);

  double xLegend = 0.60, yLegend = 0.80;
  p.makeLegend(xLegend, xLegend + 0.15, yLegend - 0.40, yLegend, "");

  double xLatex = 0.20, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesis + ", " + mystrings::sAlicePpData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sEtaV0Range075);

  string xTitle = mystrings::addUnits(sPtSoft, mystrings::sGevC, true);
  string yTitle = sNSoft;
  // yTitle = mystrings::getdYdXString(sNSoft, sPtSoft);
  yTitle = mystrings::addUnits(yTitle, mystrings::sAUShort, true);
  p.makeFrame(0., 40., 1e-3, 10., xTitle, yTitle);
  // plotutils::SetPadMargins(0.15, 0.10, 0.15, 0.05);

  const int nCols = hTrigAssoc->GetNbinsX(), nRows = hTrigAssoc->GetNbinsY();
  for (int iCol = 1; iCol <= nCols; iCol++) {
    double ptHardMin = hTrigAssoc->GetXaxis()->GetBinLowEdge(iCol);
    double ptHardMax = hTrigAssoc->GetXaxis()->GetBinUpEdge(iCol);
    string name = TString::Format("hAssoc%.f_%.f", ptHardMin, ptHardMax).Data();
    TH1D* hAssoc = (TH1D*)hTrigAssoc->ProjectionY(name.c_str(), iCol, iCol);
    if (histutils::isHistEmptyInRange(hAssoc, 1, nRows, 1.0))
      continue;

    if (inputs.printLog("plotIncl2dProjectionFromTrain() Checking projected histograms before normalisation", verbosityutils::kDebug))
      hAssoc->Print("all");

    // FIXME: What normalisation to use?
    // double factor = hAssoc->GetBinContent(hAssoc->FindFirstBinAbove(0.1));
    double factor = histutils::getUpperBound(hAssoc, false);
    // hAssoc->Scale(1. / factor, "width");
    hAssoc->Scale(1. / factor);
    string legendEntry = mystrings::addUnits(mystrings::getVarRangeString(ptHardMin, sPtHard, ptHardMax), mystrings::sGevC, false);
    p.addHistogram(hAssoc);
    p.addLegendEntry(hAssoc, legendEntry);

    if (inputs.printLog(TString::Format("plotIncl2dProjectionFromTrain() Checking projected histograms after normalisation with upper bound %f", factor).Data(), verbosityutils::kDebug))
      hAssoc->Print("all");
  } 

  p.setHistStyles();

  if (inputs.printLog("plotInJetZFromTrain1020() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotinclSharedDaughterPtFromTrain(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kInfo);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.75, 0.75);

  vector<double> ptHardBins = { 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0 };
  const int nPtHardBins = ptHardBins.size() - 1;

  const int axisPtHard     = 0;
  const int axisEta        = 1;
  const int axisPtSoft     = 2;
  const int axisPtDaughter = 3;
  THnSparse* hnSharedDaughter = inputs.GetHist<THnSparse>(typeutils::kV0InclPtSharedDaughter);
  array<int, 2> etaBins = histutils::getProjectionBins(hnSharedDaughter->GetAxis(axisEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnSharedDaughter->GetAxis(axisEta)->SetRange(etaBins[0], etaBins[1]);

  string sPtDaughter = mystrings::getPtString("daughter");
  string sPtHard     = mystrings::addSuperscript(mystrings::sPtV0, "hard");
  string sPtSoft     = mystrings::addSuperscript(mystrings::sPtV0, "soft");

  plotutils::Plotter p("sharing_inclusive_daughter.pdf", true);
  double xLegend = 0.60, yLegend = 0.80;
  p.makeLegend(xLegend, xLegend + 0.15, yLegend - 0.05 * nPtHardBins, yLegend, "");

  array<int, 2> ptBins;
  for (int i = 0; i < nPtHardBins; i++) {
    double ptHardMin = ptHardBins[i], ptHardMax = ptHardBins[i + 1];
    ptBins = histutils::getProjectionBins(hnSharedDaughter->GetAxis(axisPtHard), ptHardMin, ptHardMax);
    hnSharedDaughter->GetAxis(axisPtHard)->SetRange(ptBins[0], ptBins[1]);

    TH1D* hDaughter = (TH1D*)hnSharedDaughter->Projection(axisPtSoft);
    // hDaughter = (TH1D*)histutils::rebinHist(hDaughter, histutils::rebinnedV0PtHist("K0S", TString::Format("hDaughterRebinned%.f_%.f", ptHardMin, ptHardMax).Data()));
    hDaughter->Scale(1. / inputs.GetNevts(), "width");
    if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
      hDaughter->Print("all");
    }

    string legendEntry = mystrings::addUnits(mystrings::getVarRangeString(ptHardMin, sPtHard, ptHardMax), mystrings::sGevC, false);
    p.addHistogram(hDaughter);
    p.addLegendEntry(hDaughter, legendEntry);
  }
  
  string xTitle = sPtDaughter;
  string yTitle = "y";
  xTitle = mystrings::addUnits(xTitle, mystrings::sGevC, true);
  // yTitle = mystrings::addUnits(yTitle, mystrings::sGevC, true);
  p.makeFrame(0., ptHardBins[nPtHardBins], 1e-11, 1e-4, xTitle, yTitle);

  double xLatex = 0.35, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesis + ", " + mystrings::sAlicePpData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sEtaV0Range075);

  p.setHistStyles();

  if (inputs.printLog("plotInJetZFromTrain1020() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotincl() {
  int train = 767540;
  // plotInclFromTrain(train);
  // plotInclRatioFromTrain(train);
  // plotIncl2dFromTrain(train);
  plotIncl2dProjectionFromTrain(train);
  // plotinclSharedDaughterPtFromTrain(train);
}

// ------------------------------------------------------------------------------------
//
// Plot spectra of in-jet V0s that share daughters
//
// ------------------------------------------------------------------------------------

array<TH1D*, 2> gethistsinjetPt(InputSettings& inputs) {
  inputs.printLog(TString::Format("gethistsinjetPt() Getting ptV0 histograms for jet pT in [%.f, %.f] GeV/c and eta in [%.2f, %.2f]", inputs.getPtJetMin(), inputs.getPtJetMax(), inputs.getEtaMin(), inputs.getEtaMax()).Data(), verbosityutils::kInfo);

  const int axisJetPt   = 0;
  const int axisJetEta  = 1;
  const int axisTrigger = 2;
  const int axisAssoc   = 3;

  THnSparse* hnShared = inputs.GetHist<THnSparse>(typeutils::kV0PtInJetShared);
  TH3D* h3All = inputs.GetHist<TH3D>(typeutils::kV0PtInJet);

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
  hSharedV0sInJets->Add(hSharedAssoc);

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

void plotInJetPtFromTrain1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = gethistsinjetPt(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_pt", ".pdf"), true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-9, 0.1, mystrings::sPtV0WithUnits, mystrings::sV0PtPerJetWithUnits);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.55;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.20, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));

  
  if (inputs.printLog("plotInJetPtFromTrain1020() Checking histogram entries", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  } 
  p.plot();
} 

void plotInJetPtRatioFromTrain1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = gethistsinjetPt(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_pt", "_ratio.pdf"), true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-5, 2., mystrings::sPtV0WithUnits, mystrings::sRatio);
  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  double xLegend = 0.25, yLegend = 0.2;
  p.makeLegend(xLegend, xLegend + 0.25, yLegend, yLegend + 0.15, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.makeRatios(hAllV0sInJets);

  
  if (inputs.printLog("plotInJetPtRatioFromTrain1020() Checking histogram entries", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  } 
  p.plot();
}

void plotInJetPtFromTrain2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = gethistsinjetPt(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_pt", ".pdf"), true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-8, 0.1, mystrings::sPtV0WithUnits, mystrings::sV0PtPerJetWithUnits);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();
  
  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");
  
  double xLatex = 0.25, yLatex = 0.55;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.20, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));

  if (inputs.printLog("plotInJetPtFromTrain2030() Checking histogram entries", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  } 
  p.plot();
}

void plotInJetPtRatioFromTrain2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = gethistsinjetPt(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_pt", "_ratio.pdf"), true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-4, 2., mystrings::sPtV0WithUnits, mystrings::sRatio);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  double xLegend = 0.25, yLegend = 0.2;
  p.makeLegend(xLegend, xLegend + 0.25, yLegend, yLegend + 0.15, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");
  
  double xLatex = 0.25, yLatex = 0.75;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.makeRatios();

  if (inputs.printLog("plotInJetPtRatioFromTrain2030() Checking histogram entries", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  } 
  p.plot();
}

void plotInJetPtFromTrain3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = gethistsinjetPt(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_pt", ".pdf"), true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-8, 0.1, mystrings::sPtV0WithUnits, mystrings::sV0PtPerJetWithUnits);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.55;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.20, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));

  if (inputs.printLog("plotInJetPtFromTrain3040() Checking histogram entries", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotInJetPtRatioFromTrain3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = gethistsinjetPt(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_pt", "_ratio.pdf"), true);
  p.makeFrame(0., inputs.getPtJetMax(), 1e-4, 2., mystrings::sPtV0WithUnits, mystrings::sRatio);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  double xLegend = 0.25, yLegend = 0.2;
  p.makeLegend(xLegend, xLegend + 0.25, yLegend, yLegend + 0.15, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.75;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  p.makeRatios();

  if (inputs.printLog("plotInJetPtRatioFromTrain3040() Checking histogram entries", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotinjetpt() {
  int train = 767540;
  // plotInJetPtFromTrain1020(train);
  // plotInJetPtRatioFromTrain1020(train);
  // plotInJetPtFromTrain2030(train);
  // plotInJetPtRatioFromTrain2030(train);
  plotInJetPtFromTrain3040(train);
  plotInJetPtRatioFromTrain3040(train);
}

array<TH1D*, 2> gethistsinjetZ(InputSettings& inputs) {
  inputs.printLog(TString::Format("gethistsinjetZ() Getting zV0 histograms for jet pT in [%.f, %.f] GeV/c and eta in [%.2f, %.2f]", inputs.getPtJetMin(), inputs.getPtJetMax(), inputs.getEtaMin(), inputs.getEtaMax()).Data(), verbosityutils::kInfo);

  const int axisJetPt   = 0;
  const int axisJetEta  = 1;
  const int axisTrigger = 2;
  const int axisAssoc   = 3;

  THnSparse* hnShared = inputs.GetHist<THnSparse>(typeutils::kV0ZInJetShared);
  TH3D* h3All = inputs.GetHist<TH3D>(typeutils::kV0ZInJet);

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
  hSharedV0sInJets->Add(hSharedAssoc);

  hAllV0sInJets = (TH1D*)histutils::rebinHist(hAllV0sInJets, histutils::rebinnedV0ZHist("hAllV0sInJetsRebinned"));
  hSharedV0sInJets = (TH1D*)histutils::rebinHist(hSharedV0sInJets, histutils::rebinnedV0ZHist("hSharedV0sInJetsRebinned"));
  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    hAllV0sInJets->Print("all");
    hSharedV0sInJets->Print("all");
  }

  double nJets = inputs.GetNjets(inputs.getPtJetMin(), inputs.getPtJetMax());
  hAllV0sInJets->Scale(1. / nJets, "width");
  hSharedV0sInJets->Scale(1. / nJets, "width");

  return array<TH1D*, 2>{hAllV0sInJets, hSharedV0sInJets};
}

void plotInJetZFromTrain1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = gethistsinjetZ(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_z", ".pdf"), true);
  p.makeFrame(1e-3, 1+1e-3, 1e-8, 1., mystrings::sZV0, mystrings::sV0ZPerJet);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.55;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.20, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));

  if (inputs.printLog("plotInJetZFromTrain1020() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotInJetZRatioFromTrain1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  array<TH1D*, 2> hists = gethistsinjetZ(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_z", "_ratio.pdf"), true);
  p.makeFrame(1e-3, 1+1e-3, 1e-5, 2., mystrings::sZV0, mystrings::sV0ZPerJet);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.55;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.20, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));

  p.makeRatios();

  if (inputs.printLog("plotInJetZRatioFromTrain1020() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotInJetZFromTrain2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = gethistsinjetZ(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_z", ".pdf"), true);
  p.makeFrame(1e-3, 1+1e-3, 1e-8, 1., mystrings::sZV0, mystrings::sV0ZPerJet);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.55;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.20, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));

  if (inputs.printLog("plotInJetZFromTrain2030() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotInJetZRatioFromTrain2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  array<TH1D*, 2> hists = gethistsinjetZ(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_z", "_ratio.pdf"), true);
  p.makeFrame(1e-3, 1+1e-3, 1e-4, 2., mystrings::sZV0, mystrings::sV0ZPerJet);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.55;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.20, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));

  p.makeRatios();

  if (inputs.printLog("plotInJetZRatioFromTrain2030() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotInJetZFromTrain3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = gethistsinjetZ(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_z", ".pdf"), true);
  p.makeFrame(1e-3, 1+1e-3, 1e-8, 1., mystrings::sZV0, mystrings::sV0ZPerJet);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.55;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.20, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));

  if (inputs.printLog("plotInJetZFromTrain3040() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotInJetZRatioFromTrain3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  array<TH1D*, 2> hists = gethistsinjetZ(inputs);
  TH1D* hAllV0sInJets = hists[0];
  TH1D* hSharedV0sInJets = hists[1];

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_z", "_ratio.pdf"), true);
  p.makeFrame(1e-3, 1+1e-3, 1e-4, 2., mystrings::sZV0, mystrings::sV0ZPerJet);

  p.setHists({hAllV0sInJets, hSharedV0sInJets});
  p.setHistStyles();

  p.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  p.addLegendEntry(hAllV0sInJets, "All V0s in jets");
  p.addLegendEntry(hSharedV0sInJets, "V0s with shared daughters in jets");

  double xLatex = 0.25, yLatex = 0.55;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntiktJets);
  p.addLatex(xLatex, yLatex - 0.15, mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.20, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));

  p.makeRatios();

  if (inputs.printLog("plotInJetZRatioFromTrain3040() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotInJetZCorrelationFromTrain1020(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(10., 20.);

  const int axisJetPt = 0;
  const int axisJetEta = 1;
  const int axisTrigger = 2;
  const int axisAssoc = 3;
  THnSparse* hnShared = inputs.GetHist<THnSparse>(typeutils::kV0ZInJetShared);
  array<int, 2> ptBins = histutils::getProjectionBins(hnShared->GetAxis(axisJetPt), inputs.getPtJetMin(), inputs.getPtJetMax());
  array<int, 2> etaBins = histutils::getProjectionBins(hnShared->GetAxis(axisJetEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnShared->GetAxis(axisJetPt)->SetRange(ptBins[0], ptBins[1]);
  hnShared->GetAxis(axisJetEta)->SetRange(etaBins[0], etaBins[1]);
  TH2D* hTrigAssoc = (TH2D*)hnShared->Projection(axisAssoc, axisTrigger);
  hTrigAssoc = (TH2D*)histutils::rebinHist2D(hTrigAssoc, histutils::rebinnedV0ZHist2D("hTrigAssocRebinned"));

  // hTrigAssoc->Scale(1. / inputs.GetNjets(inputs.getPtJetMin(), inputs.getPtJetMax()), "width");
  histutils::normaliseHistColByCol(hTrigAssoc);

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_z", "_correlation.pdf"), false, false, false);
  p.setDrawOption("colz");
  p.addHistogram(hTrigAssoc);
  p.setZAxisRange(1e-5, 1.);
  p.setFrameTitleOffsetY(plotutils::defaultValues::defTitleOffset);
  p.setPadMargins(0.15, 0.10, 0.15, 0.05);

  string xTitle = mystrings::addSuperscript(mystrings::sZV0, "hard");
  string yTitle = mystrings::addSuperscript(mystrings::sZV0, "soft");
  p.makeFrame(1e-3, 1+1e-3, 1e-8, 1., xTitle, yTitle);

  double xLatex = 0.25, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  
  if (inputs.printLog("plotInJetZCorrelationFromTrain1020() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotInJetZCorrelationFromTrain2030(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(20., 30.);

  const int axisJetPt = 0;
  const int axisJetEta = 1;
  const int axisTrigger = 2;
  const int axisAssoc = 3;
  THnSparse* hnShared = inputs.GetHist<THnSparse>(typeutils::kV0ZInJetShared);
  array<int, 2> ptBins = histutils::getProjectionBins(hnShared->GetAxis(axisJetPt), inputs.getPtJetMin(), inputs.getPtJetMax());
  array<int, 2> etaBins = histutils::getProjectionBins(hnShared->GetAxis(axisJetEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnShared->GetAxis(axisJetPt)->SetRange(ptBins[0], ptBins[1]);
  hnShared->GetAxis(axisJetEta)->SetRange(etaBins[0], etaBins[1]);
  TH2D* hTrigAssoc = (TH2D*)hnShared->Projection(axisAssoc, axisTrigger);
  hTrigAssoc = (TH2D*)histutils::rebinHist2D(hTrigAssoc, histutils::rebinnedV0ZHist2D("hTrigAssocRebinned"));

  // hTrigAssoc->Scale(1. / inputs.GetNjets(inputs.getPtJetMin(), inputs.getPtJetMax()), "width");
  histutils::normaliseHistColByCol(hTrigAssoc);

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_z", "_correlation.pdf"), false, false, false);
  p.setDrawOption("colz");
  p.addHistogram(hTrigAssoc);
  p.setZAxisRange(1e-5, 1.);
  p.setFrameTitleOffsetY(plotutils::defaultValues::defTitleOffset);
  p.setPadMargins(0.15, 0.10, 0.15, 0.05);

  string xTitle = mystrings::addSuperscript(mystrings::sZV0, "hard");
  string yTitle = mystrings::addSuperscript(mystrings::sZV0, "soft");
  p.makeFrame(1e-3, 1+1e-3, 1e-8, 1., xTitle, yTitle);

  double xLatex = 0.25, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  

  if (inputs.printLog("plotInJetZCorrelationFromTrain2030() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotInJetZCorrelationFromTrain3040(int train) {
  InputSettings inputs; inputs.setVerbosity(verbosityutils::kDebug);
  inputs.setTrain(train);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);
  inputs.setPtJet(30., 40.);

  const int axisJetPt = 0;
  const int axisJetEta = 1;
  const int axisTrigger = 2;
  const int axisAssoc = 3;
  THnSparse* hnShared = inputs.GetHist<THnSparse>(typeutils::kV0ZInJetShared);
  array<int, 2> ptBins = histutils::getProjectionBins(hnShared->GetAxis(axisJetPt), inputs.getPtJetMin(), inputs.getPtJetMax());
  array<int, 2> etaBins = histutils::getProjectionBins(hnShared->GetAxis(axisJetEta), inputs.getEtaMin(), inputs.getEtaMax());
  hnShared->GetAxis(axisJetPt)->SetRange(ptBins[0], ptBins[1]);
  hnShared->GetAxis(axisJetEta)->SetRange(etaBins[0], etaBins[1]);
  TH2D* hTrigAssoc = (TH2D*)hnShared->Projection(axisAssoc, axisTrigger);
  hTrigAssoc = (TH2D*)histutils::rebinHist2D(hTrigAssoc, histutils::rebinnedV0ZHist2D("hTrigAssocRebinned"));

  // hTrigAssoc->Scale(1. / inputs.GetNjets(inputs.getPtJetMin(), inputs.getPtJetMax()), "width");
  histutils::normaliseHistColByCol(hTrigAssoc);

  plotutils::Plotter p(inputs.GetNameFromPtJet("sharing_injet_z", "_correlation.pdf"), false, false, false);
  p.setDrawOption("colz");
  p.addHistogram(hTrigAssoc);
  p.setZAxisRange(1e-5, 1.);
  p.setFrameTitleOffsetY(plotutils::defaultValues::defTitleOffset);
  p.setPadMargins(0.15, 0.10, 0.15, 0.05);

  string xTitle = mystrings::addSuperscript(mystrings::sZV0, "hard");
  string yTitle = mystrings::addSuperscript(mystrings::sZV0, "soft");
  p.makeFrame(1e-3, 1+1e-3, 1e-8, 1., xTitle, yTitle);

  double xLatex = 0.25, yLatex = 0.85;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData + ", " + mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sAntiktJets + ", " + mystrings::sJetRadius04 + ", " + mystrings::sEtaJetRange035);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::getPtJetRangeString(inputs.getPtJetMin(), inputs.getPtJetMax(), true));
  

  if (inputs.printLog("plotInJetZCorrelationFromTrain3040() Checking histogram values", verbosityutils::kDebug)) {
    for (auto h : p.getHists())
      h->Print("all");
  }
  p.plot();
}

void plotinjetz() {
  int train = 767540;
  // gROOT->SetBatch(true);
  // plotInJetZFromTrain1020(train);
  // plotInJetZRatioFromTrain1020(train);
  // plotInJetZFromTrain2030(train);
  // plotInJetZRatioFromTrain2030(train);
  // plotInJetZFromTrain3040(train);
  // plotInJetZRatioFromTrain3040(train);
  // plotInJetZCorrelationFromTrain1020(train);
  // plotInJetZCorrelationFromTrain2030(train);
  plotInJetZCorrelationFromTrain3040(train);
}


// ------------------------------------------------------------------------------------
//
// Plot spectra of jets, jets with V0s, and jets with V0s that share daughters
//
// ------------------------------------------------------------------------------------

array<TH1D*, 3> gethistsjet(InputSettings& inputs) {
  inputs.printLog(TString::Format("gethistjet() Getting histograms for jet pT with eta in [%.2f, %.2f]", inputs.getEtaMin(), inputs.getEtaMax()).Data(), verbosityutils::kInfo);

  TH3D* h3All     = inputs.GetHist<TH3D>(typeutils::kJetPtIncl); // All jets
  TH3D* h3AllwV0s = inputs.GetHist<TH3D>(typeutils::kJetPtWithV0s); // All jets with V0 candidates
  TH3D* h3Shared  = inputs.GetHist<TH3D>(typeutils::kJetPtShared); // All jets with V0s that share daughters

  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    h3All->Print();
    h3AllwV0s->Print();
    h3Shared->Print();
  }

  array<int, 2> etaBins = histutils::getProjectionBins(h3All->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hAll = h3All->ProjectionX("hAll", etaBins[0], etaBins[1], 0, 1 + h3All->GetNbinsZ());

  etaBins = histutils::getProjectionBins(h3AllwV0s->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hAllwV0s = h3AllwV0s->ProjectionX("hAllwV0s", etaBins[0], etaBins[1], 0, 1 + h3AllwV0s->GetNbinsZ());

  etaBins = histutils::getProjectionBins(h3Shared->GetYaxis(), inputs.getEtaMin(), inputs.getEtaMax());
  TH1D* hShared = h3Shared->ProjectionX("hShared", etaBins[0], etaBins[1], 0, 1 + h3Shared->GetNbinsZ());

  if (inputs.passVerbosityCheck(verbosityutils::kDebug)) {
    h3All->Print();
    h3AllwV0s->Print();
    h3Shared->Print();
  }

  return array<TH1D*, 3>{hAll, hAllwV0s, hShared};
}

void plotjetall3() {
  TH1D* jetTemplate = new TH1D("hJetPtTemplate", "hJetPtTemplate", 10, 0., 50.);
  InputSettings inputs;
  inputs.setTrain(538318);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);

  array<TH1D*, 3> hists = gethistsjet(inputs);
  TH1D* hAll     = (TH1D*)histutils::rebinHist(hists[0], jetTemplate);
  TH1D* hAllwV0s = (TH1D*)histutils::rebinHist(hists[1], jetTemplate);
  TH1D* hShared  = (TH1D*)histutils::rebinHist(hists[2], jetTemplate);

  double nevts = inputs.GetNevts();
  hAll->Scale(1. / nevts, "width");
  hAllwV0s->Scale(1. / nevts, "width");
  hShared->Scale(1. / nevts, "width");

  plotutils::Plotter plSpectra("sharing_jets_all.pdf", false);
  plSpectra.setHists({hAll, hAllwV0s, hShared});
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

  plotutils::Plotter plRatio("sharing_jets_ratio_all.pdf", true);
  plRatio.makeFrame(0., 50., 1e-4, 2., mystrings::sPtJet, "Ratio");

  plRatio.setHists(plSpectra.getHists()); // Should ensure the histograms are the same
  plRatio.setHistStyles();
  plRatio.makeRatios(0);

  plRatio.makeLegend(0.25, 0.50, 0.40, 0.50, "");
  plRatio.addLegendEntry(hAll, "Inclusive jets");
  plRatio.addLegendEntry(hAllwV0s, "Jets with V0s");
  plRatio.addLegendEntry(hShared, "Jets with V0s that share daughters");

  xLatex = 0.25, yLatex = 0.75;
  plRatio.addLatex(xLatex, yLatex, "This Thesis, ALICE pp data");
  plRatio.addLatex(xLatex, yLatex - 0.05, "#sqrt{s} = 13.6 TeV");
  plRatio.addLatex(xLatex, yLatex - 0.10, "Anti-#it{k}_{T} ch+V0 jets");
  plRatio.addLatex(xLatex, yLatex - 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.35");

  plRatio.plot();
}

void plotjetallvsv0() {
  TH1D* jetTemplate = new TH1D("hJetPtTemplate", "hJetPtTemplate", 10, 0., 50.);
  InputSettings inputs;
  inputs.setTrain(538318);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);

  array<TH1D*, 3> hists = gethistsjet(inputs);
  TH1D* hAll     = (TH1D*)histutils::rebinHist(hists[0], jetTemplate);
  TH1D* hAllwV0s = (TH1D*)histutils::rebinHist(hists[1], jetTemplate);
  TH1D* hShared  = (TH1D*)histutils::rebinHist(hists[2], jetTemplate);

  double nevts = inputs.GetNevts();
  hAll->Scale(1. / nevts, "width");
  hAllwV0s->Scale(1. / nevts, "width");
  hShared->Scale(1. / nevts, "width");

  plotutils::Plotter plSpectra("sharing_jets_all-v0s.pdf", false);
  plSpectra.setHists({hAll, hAllwV0s});
  plSpectra.setHistStyles();

  plSpectra.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  plSpectra.addLegendEntry(hAll, "Inclusive jets");
  plSpectra.addLegendEntry(hAllwV0s, "Jets with V0s");

  double xLatex = 0.25, yLatex = 0.55;
  plSpectra.addLatex(xLatex, yLatex, "This Thesis, ALICE pp data");
  plSpectra.addLatex(xLatex, yLatex - 0.05, "#sqrt{s} = 13.6 TeV");
  plSpectra.addLatex(xLatex, yLatex - 0.10, "Anti-#it{k}_{T} ch+V0 jets");
  plSpectra.addLatex(xLatex, yLatex - 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.35");

  plSpectra.makeFrame(mystrings::sPtJet, mystrings::sJetsPerEvent);
  plSpectra.plot();

  plotutils::Plotter plRatio("sharing_jets_ratio_all-v0s.pdf", true);
  plRatio.makeFrame(0., 50., 1e-2, 2., mystrings::sPtJet, "Ratio");

  plRatio.setHists(plSpectra.getHists()); // Should ensure the histograms are the same
  plRatio.setHistStyles();
  plRatio.makeRatios(0);

  plRatio.makeLegend(0.25, 0.50, 0.40, 0.50, "");
  plRatio.addLegendEntry(hAll, "Inclusive jets");
  plRatio.addLegendEntry(hAllwV0s, "Jets with V0s");

  xLatex = 0.25, yLatex = 0.75;
  plRatio.addLatex(xLatex, yLatex, "This Thesis, ALICE pp data");
  plRatio.addLatex(xLatex, yLatex - 0.05, "#sqrt{s} = 13.6 TeV");
  plRatio.addLatex(xLatex, yLatex - 0.10, "Anti-#it{k}_{T} ch+V0 jets");
  plRatio.addLatex(xLatex, yLatex - 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.35");

  plRatio.plot();
}

void plotjetwv0vsshared() {
  TH1D* jetTemplate = new TH1D("hJetPtTemplate", "hJetPtTemplate", 10, 0., 50.);
  InputSettings inputs;
  inputs.setTrain(538318);
  inputs.SetInputFileNameFromTrain();
  inputs.setEta(-0.35, 0.35);

  array<TH1D*, 3> hists = gethistsjet(inputs);
  TH1D* hAll     = (TH1D*)histutils::rebinHist(hists[0], jetTemplate);
  TH1D* hAllwV0s = (TH1D*)histutils::rebinHist(hists[1], jetTemplate);
  TH1D* hShared  = (TH1D*)histutils::rebinHist(hists[2], jetTemplate);
  double nevts = inputs.GetNevts();
  hAll->Scale(1. / nevts, "width");
  hAllwV0s->Scale(1. / nevts, "width");
  hShared->Scale(1. / nevts, "width");

  plotutils::Plotter plSpectra("sharing_jets_shared.pdf", false);
  plSpectra.setHists({hAllwV0s, hShared});
  plSpectra.setHistStyles();

  plSpectra.makeLegend(0.25, 0.50, 0.20, 0.30, "");
  plSpectra.addLegendEntry(hAllwV0s, "Jets with V0s");
  plSpectra.addLegendEntry(hShared, "Jets with V0s that share daughters");

  double xLatex = 0.25, yLatex = 0.55;
  plSpectra.addLatex(xLatex, yLatex, "This Thesis, ALICE pp data");
  plSpectra.addLatex(xLatex, yLatex - 0.05, "#sqrt{s} = 13.6 TeV");
  plSpectra.addLatex(xLatex, yLatex - 0.10, "Anti-#it{k}_{T} ch+V0 jets");
  plSpectra.addLatex(xLatex, yLatex - 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.35");

  plSpectra.makeFrame(mystrings::sPtJet, mystrings::sJetsPerEvent);
  plSpectra.plot();

  plotutils::Plotter plRatio("sharing_jets_ratio_shared.pdf", true);
  plRatio.makeFrame(0., 50., 1e-4, 2., mystrings::sPtJet, "Ratio");

  plRatio.setHists(plSpectra.getHists()); // Should ensure the histograms are the same
  plRatio.setHistStyles();
  plRatio.makeRatios(0);

  plRatio.makeLegend(0.25, 0.50, 0.40, 0.50, "");
  plRatio.addLegendEntry(hAllwV0s, "Jets with V0s");
  plRatio.addLegendEntry(hShared, "Jets with V0s that share daughters");

  xLatex = 0.25, yLatex = 0.75;
  plRatio.addLatex(xLatex, yLatex, "This Thesis, ALICE pp data");
  plRatio.addLatex(xLatex, yLatex - 0.05, "#sqrt{s} = 13.6 TeV");
  plRatio.addLatex(xLatex, yLatex - 0.10, "Anti-#it{k}_{T} ch+V0 jets");
  plRatio.addLatex(xLatex, yLatex - 0.15, "#it{R} = 0.4, |#eta_{jet}| < 0.35");

  plRatio.plot();
}

#endif