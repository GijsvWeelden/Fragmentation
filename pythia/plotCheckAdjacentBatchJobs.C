
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

#include "../plotting_macros/histUtils.C"
#include "../plotting_macros/plotUtils.C"

// ------------------------------------------------
// Compare the spectra of adjacent jobs submitted with the condor batch system.
// The pythia seed is determined based on the time, but adjacent jobs may be submitted at the same time, thereby leading to duplicate seeds
// ------------------------------------------------

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
      case kErrors:   return "Error";
      case kWarnings: return "Warning";
      case kInfo:     return "Info";
      case kDebug:    return "Debug";
      case kDebugMax: return "DebugMax";
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

namespace histogramtypeutils {
  enum HistType { hTrack, hV0, hV0Jet, hzV0_K0S, hzV0_Lambda0, hzV0_V0 };
  bool is_valid(int t) {
    bool b = (t >= hTrack && t <= hzV0_V0);
    if (!b)
      cout << "histogramtypeutils Error: invalid HistType " << t << endl;
    return b;
  }
  string to_string(HistType t) {
    switch (t) {
      case hTrack:       return "hTrack";
      case hV0:          return "hV0";
      case hV0Jet:       return "hV0Jet";
      case hzV0_K0S:     return "hzV0_K0S";
      case hzV0_Lambda0: return "hzV0_Lambda0";
      case hzV0_V0:      return "hzV0_V0";
      default:           return "Unknown";
    }
  }
} // namespace histogramtypeutils

struct InputSettings{
  private:
    int _batch = 0;
    vector<int> _jobs = {}; // All jobs to compare
    histogramtypeutils::HistType _histtype = histogramtypeutils::hTrack;
    string _outputFileName = "";
    verbosityutils::Verbosity _verbosity = verbosityutils::kWarnings;

  public:
    // Getters and Setters
    int getBatch() { return _batch; }
    histogramtypeutils::HistType getHistType() { return _histtype; }
    int getJob(int i) { return _jobs[i]; }
    vector<int> getJobs() { return _jobs; }
    int getNJobs() { return _jobs.size(); }
    string getOutputFileName() { return _outputFileName; }
    verbosityutils::Verbosity getVerbosity() { return _verbosity; }
    
    void setBatch(int b) {_batch = b; }
    void setHistType(histogramtypeutils::HistType t) { _histtype = t; }
    void setJobs(vector<int> j) { _jobs = j; }
    void setOutputFileName(string s) { _outputFileName = s; }
    void setVerbosity(verbosityutils::Verbosity v) { _verbosity = v; }

    // Methods
    string formatHistName() {
      return histogramtypeutils::to_string(_histtype);
    }
    bool passVerbosityCheck(verbosityutils::Verbosity level) {
      return verbosityutils::passVerbosityCheck(level, _verbosity);
    }
    void printLog(string message, verbosityutils::Verbosity messageVerbLevel) {
      if (passVerbosityCheck(messageVerbLevel))
        verbosityutils::printLog(message, messageVerbLevel);
    }
};

string formatInputFileName(int batchNumber, int jobNumber) {
  string sb = to_string(batchNumber);
  string sj = to_string(jobNumber);
  string s = "V0JetClustering/Batch" + sb + "/" + sb + "_" + sj + ".root";
  return s;
}

string formatOutputFileName(string prefix, int batchNumber, vector<int> jobNumbers, bool ratio) {
  string s = prefix + "_b" + to_string(batchNumber) + "_j";
  for (int i = 0; i < jobNumbers.size(); i++) {
    s += to_string(jobNumbers[i]);
    if (i != jobNumbers.size() - 1)
      s += "-";
  }
  if (ratio)
    s += "_ratio";
  
  s += ".pdf";
  return s;
}

TFile* getFile(InputSettings &x, string fileName) {
  x.printLog("getFile() Opening file " + fileName, verbosityutils::kInfo);
  TFile* file = TFile::Open(fileName.c_str(), "READ");
  if (!file)
    x.printLog("Could not open file " + fileName, verbosityutils::kErrors);

  return file;
}

template <typename T> T* getHist(InputSettings &x, int index) {
  string inputFileName = formatInputFileName(x.getBatch(), x.getJob(index));
  TFile* inputFile = getFile(x, inputFileName);
  if (!inputFile)
    return nullptr;

  string histName = x.formatHistName();
  x.printLog("getHist() Retrieving histogram " + histName, verbosityutils::kInfo);
  
  T* hist = (T*)inputFile->Get(histName.c_str());
  if (!hist)
    x.printLog(TString::Format("Could not retrieve histogram %s from file %s", histName.c_str(), inputFileName.c_str()).Data(), verbosityutils::kErrors);

  return hist;
}

void plotV0s(bool makeRatios = false) {
  InputSettings x; x.setVerbosity(verbosityutils::kDebug);
  x.setBatch(3226091);
  // x.setJobs({0, 1, 2});
  x.setJobs({1, 2});
  x.setHistType(histogramtypeutils::hV0);
  x.setOutputFileName(formatOutputFileName("inclV0", x.getBatch(), x.getJobs(), makeRatios));

  bool logplot = true;
  if (makeRatios)
    logplot = false;
  plotutils::Plotter p(x.getOutputFileName(), logplot);

  p.makeLegend(0.5, 0.7, 0.7, 0.9, TString::Format("Batch %d", x.getBatch()).Data());

  double xLatex = 0.5, yLatex = 0.4;
  p.addLatex(xLatex, yLatex, "Pythia simulation pp");

  double xMinFrame = 0., xMaxFrame = 60.;
  double yMinFrame = 1e-1, yMaxFrame = 1e7;
  string xTitle = "#it{p}_{T, V0}", yTitle = "#it{N}_{V0}";
  if (makeRatios) {
    yMinFrame = 0., yMaxFrame = 2.;
    yTitle = "Ratio";
  }
  p.makeFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);

  for (int iJob = 0; iJob < x.getNJobs(); iJob++) {
    TH3D* h3V0 = getHist<TH3D>(x, iJob);
    TH1D* hV0 = (TH1D*)h3V0->ProjectionX(TString::Format("hV0Pt_%d", iJob).Data());
    p.addHistogram(hV0);
    p.addLegendEntry(hV0, TString::Format("Job %d", iJob).Data());
  }
  if (makeRatios)
    p.makeRatios();

  p.setHistStyles();
  p.plot();
}
