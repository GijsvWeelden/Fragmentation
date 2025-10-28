
#include "dounfolding.C"
#include "../plotting_macros/myStrings.C"

// !!! Global Variable !!!
string _fileNameSingle = "ClosureTest_527899_80_ptjetgen10-65-5_ptjetrec10-60-5.root";
string _fileNameNarrow = "ClosureTest_527899_80_ptjetgen10-65-5_ptjetrec10-60-5.root";
string _fileNameWide   = "ClosureTest_527899_80_ptjetgen5-80-5_ptjetrec10-60-5.root";

void SwapNames() {
  string tmp      = _fileNameNarrow;
  _fileNameNarrow = _fileNameWide;
  _fileNameWide   = tmp;
}

string RemovePrependingJunk(const string s, const string premarker, const string postmarker) {
  // Remove everything before the marker
  string t = s;
  size_t pos = t.find(premarker);
  if (pos != string::npos)
    t.erase(0, pos);

  pos = t.find(postmarker);
  if (pos != string::npos)
    t.erase(pos, string::npos);

  return t;
}

TFile* GetFile(InputSettings& inputs) {
  TFile* file = TFile::Open(inputs.inputFileName.c_str(), "READ");
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return nullptr;
  }
  return file;
}

template <typename T>
T* GetHist(InputSettings& inputs, string histName) {
  TFile* file = GetFile(inputs);
  if (!file)
    return nullptr;

  // Get the histogram from the file
  T* hist = dynamic_cast<T*>(file->Get(histName.c_str()));
  if (!hist) {
    inputs.printLog("Error: could not find histogram " + histName, verbosityutilities::kErrors);
    return nullptr;
  }

  return hist;
}

void GetJetIterations(InputSettings& inputs, plotutils::Plotter& plotter, vector<int> iterations, bool isUnfolded) {
  if (iterations.empty()) {
    cout << "GetJetIterations() Error: no iterations provided!" << endl;
    return;
  }
  if (!plotter.getLegend()) {
    cout << "GetJetIterations() Error: plotter does not have a legend!" << endl;
    return;
  }

  string histName = isUnfolded ? "unfoldedJets" : "refoldedJets";
  for (int i : iterations) {
    TH1* hist = GetHist<TH1>(inputs, histName + to_string(i));
    if (!hist)
      continue;

    plotutils::setStyle(hist, i);
    plotter.addLegendEntry(hist, to_string(i));
    plotter.addHistogram(hist);
  }
}

void GetJetIterations(InputSettings& inputs, plotutils::Plotter& plotter, bool isUnfolded) {
  if (inputs.maxIteration < inputs.minIteration) {
    inputs.printLog("GetJetIterations() Error: invalid iteration range " + to_string(inputs.minIteration) + " to " + to_string(inputs.maxIteration), verbosityutilities::kErrors);
    return;
  }

  vector<int> iterations;
  for (int i = inputs.minIteration; i <= inputs.maxIteration; i++)
    iterations.push_back(i);

  GetJetIterations(inputs, plotter, iterations, isUnfolded);
}

void GetK0SIterations(InputSettings& inputs, plotutils::Plotter& plotter, vector<int> iterations, bool isUnfolded, bool doZ) {
  if (iterations.empty()) {
    inputs.printLog("GetK0SIterations() Error: no iterations provided!", verbosityutilities::kErrors);
    return;
  }
  if (!plotter.getLegend()) {
    inputs.printLog("GetK0SIterations() Error: plotter does not have a legend!", verbosityutilities::kErrors);
    return;
  }

  if (doZ)
    inputs.printLog("GetK0SIterations() getting Z histograms", verbosityutilities::kInfo);
  else
    inputs.printLog("GetK0SIterations() getting Pt histograms", verbosityutilities::kInfo);

  string baseName = isUnfolded ? rmutilities::testing::nameGenK0SPt : rmutilities::testing::nameRecK0SPt;
  string histName = isUnfolded ? "unfoldedK0SPt" : "refoldedK0SPt";
  if (doZ) {
    baseName = isUnfolded ? rmutilities::testing::nameGenK0SZ : rmutilities::testing::nameRecK0SZ;
    histName = isUnfolded ? "unfoldedK0SZ" : "refoldedK0SZ";
  }

  for (auto i : iterations) {
    TH2* hist = GetHist<TH2>(inputs, histName + to_string(i));
    if (!hist)
      continue;

    array<int, 2> ptjetBins = histutils::getProjectionBins(hist->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
    TH1* hProj = (TH1*)hist->ProjectionX(inputs.getNameFromPtJetProjection(histName + to_string(i), "").c_str(), ptjetBins[0], ptjetBins[1]);
    plotutils::setStyle(hProj, i);
    plotter.addLegendEntry(hProj, to_string(i));
    plotter.addHistogram(hProj);
  }
}

void GetK0SIterations(InputSettings& inputs, plotutils::Plotter& plotter, bool isUnfolded, bool doZ) {
  if (inputs.maxIteration < inputs.minIteration) {
    inputs.printLog("GetK0SIterations() Error: invalid iteration range " + to_string(inputs.minIteration) + " to " + to_string(inputs.maxIteration), verbosityutilities::kErrors);
    return;
  }

  vector<int> iterations;
  for (int i = inputs.minIteration; i <= inputs.maxIteration; i++)
    iterations.push_back(i);

  GetK0SIterations(inputs, plotter, iterations, isUnfolded, doZ);
}

// -------------------------------------------------------------------------------------------------
//
// Jet Plotters
//
// -------------------------------------------------------------------------------------------------

void plotjetsRefolded(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  inputs.printLog("Plotting refolded spectra", verbosityutilities::kInfo);
  plotutils::Plotter p(inputs.getNameFromIterations("refoldedJets", ".pdf"), true, 0.04);
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  p.makeFrame(10., 60., 1e2, 1e6, mystrings::sPtJet, "Refolded");
  GetJetIterations(inputs, p, !isUnfolded);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotjetsRefoldedOverRec(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.inputFileName = _fileNameSingle;

  TH1* hRec = GetHist<TH1>(inputs, rmutilities::testing::nameRecJets);

  const bool isUnfolded = true;
  inputs.printLog("Plotting refolded ratios", verbosityutilities::kInfo);
  plotutils::Plotter p(inputs.getNameFromIterations("refoldedJetsRatio", ".pdf"), false, 0.04);
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  p.makeFrame(10., 60., 0.8, 1.2, mystrings::sPtJet, "#frac{Refolded}{Reconstructed}");
  GetJetIterations(inputs, p, !isUnfolded);
  p.makeRatios(hRec);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotjetsUnfolded(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  TH1* hGen = GetHist<TH1>(inputs, rmutilities::testing::nameGenJets);

  inputs.printLog("Plotting unfolded spectra", verbosityutilities::kInfo);
  plotutils::Plotter p(inputs.getNameFromIterations("unfoldedJets", ".pdf"), true, 0.04);
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  p.makeFrame(5., 80., 1e2, 1e6, mystrings::sPtJet, "Unfolded");
  GetJetIterations(inputs, p, isUnfolded);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotjetsUnfoldedOverGen(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs; inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.inputFileName = _fileNameSingle;

  TH1* hGen = GetHist<TH1>(inputs, rmutilities::testing::nameGenJets);

  const bool isUnfolded = true;
  inputs.printLog("Plotting unfolded ratios", verbosityutilities::kInfo);
  plotutils::Plotter p(inputs.getNameFromIterations("unfoldedJetsRatio", ".pdf"), false, 0.04);
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  p.makeFrame(5., 80., 0.8, 1.2, mystrings::sPtJet, "#frac{Unfolded}{Generated}");
  GetJetIterations(inputs, p, isUnfolded);
  p.makeRatios(hGen);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotjetsGenRecUnfolded(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs; inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.inputFileName = _fileNameSingle;

  TH1* hGen = GetHist<TH1>(inputs, rmutilities::testing::nameGenJets);
  plotutils::setStyle(hGen, 1);

  TH1* hRec = GetHist<TH1>(inputs, rmutilities::testing::nameRecJets);
  plotutils::setStyle(hRec, 2);

  const bool isUnfolded = true;
  inputs.printLog("Plotting unfolded ratios", verbosityutilities::kInfo);
  plotutils::Plotter p(inputs.getNameFromIterations("unfoldedJetsRatio", ".pdf"), true, 0.04);
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  p.makeFrame(5., 80., 1e-6, 1., mystrings::sPtJet, "");
  GetJetIterations(inputs, p, isUnfolded);
  p.addHistogram(hGen); p.addLegendEntry(hGen, "Generated");
  p.addHistogram(hRec); p.addLegendEntry(hRec, "Reconstructed");
  p.makeRatios(hGen);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotalljets(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kErrors) {
  gROOT->SetBatch(true);
  plotjetsUnfolded(minIteration, maxIteration, v);
  plotjetsUnfoldedOverGen(minIteration, maxIteration, v);
  plotjetsRefolded(minIteration, maxIteration, v);
  plotjetsRefoldedOverRec(minIteration, maxIteration, v);
}

// -------------------------------------------------------------------------------------------------
//
// K0S Pt Plotters
//
// -------------------------------------------------------------------------------------------------

void plotK0SPtUnfolded1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S pt spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SPt", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(0., 20., 1e-5, 1., mystrings::sPtK0S, "Unfolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, !doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtUnfoldedOverGen1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Gen = GetHist<TH2>(inputs, rmutilities::testing::nameGenK0SPt);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Gen->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hGen = (TH1*)h2Gen->ProjectionX("hGen", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S pt ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SPtRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(0., 20., 0.8, 1.2, mystrings::sPtK0S, "#frac{Unfolded}{Generated}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, !doZ);
  p.makeRatios(hGen);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtRefolded1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S pt spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SPt", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(0., 20., 1e-6, 1., mystrings::sPtK0S, "Refolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, !doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtRefoldedOverRec1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Rec = GetHist<TH2>(inputs, rmutilities::testing::nameRecK0SPt);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Rec->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hRec = (TH1*)h2Rec->ProjectionX("hRec", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S pt ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SPtRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(0., 20., 0.8, 1.2, mystrings::sPtK0S, "#frac{Refolded}{Reconstructed}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, !doZ);
  p.makeRatios(hRec);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtUnfolded2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S pt spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SPt", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(0., 30., 1e-5, 1., mystrings::sPtK0S, "Unfolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, !doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtUnfoldedOverGen2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Gen = GetHist<TH2>(inputs, rmutilities::testing::nameGenK0SPt);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Gen->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hGen = (TH1*)h2Gen->ProjectionX("hGen", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S pt ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SPtRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(0., 30., 0.8, 1.2, mystrings::sPtK0S, "#frac{Unfolded}{Generated}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, !doZ);
  p.makeRatios(hGen);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtRefolded2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S pt spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SPt", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(0., 30., 1e-5, 1., mystrings::sPtK0S, "Refolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, !doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtRefoldedOverRec2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Rec = GetHist<TH2>(inputs, rmutilities::testing::nameRecK0SPt);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Rec->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hRec = (TH1*)h2Rec->ProjectionX("hRec", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S pt ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SPtRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(0., 30., 0.8, 1.2, mystrings::sPtK0S, "#frac{Refolded}{Reconstructed}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, !doZ);
  p.makeRatios(hRec);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtUnfolded3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(30., 40.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S pt spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SPt", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(0., 40., 1e-5, 1., mystrings::sPtK0S, "Unfolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, !doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtUnfoldedOverGen3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(30., 40.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Gen = GetHist<TH2>(inputs, rmutilities::testing::nameGenK0SPt);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Gen->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hGen = (TH1*)h2Gen->ProjectionX("hGen", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S pt ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SPtRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(0., 40., 0.8, 1.2, mystrings::sPtK0S, "#frac{Unfolded}{Generated}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, !doZ);
  p.makeRatios(hGen);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtRefolded3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(30., 40.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S pt spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SPt", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(0., 40., 1e-5, 1., mystrings::sPtK0S, "Refolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, !doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SPtRefoldedOverRec3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(30., 40.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Rec = GetHist<TH2>(inputs, rmutilities::testing::nameRecK0SPt);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Rec->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hRec = (TH1*)h2Rec->ProjectionX("hRec", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S pt ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SPtRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(0., 40., 0.8, 1.2, mystrings::sPtK0S, "#frac{Refolded}{Reconstructed}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, !doZ);
  p.makeRatios(hRec);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotallK0SPt(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kErrors) {
  gROOT->SetBatch(true);
  plotK0SPtUnfolded1020(minIteration, maxIteration, v);
  plotK0SPtUnfoldedOverGen1020(minIteration, maxIteration, v);
  plotK0SPtRefolded1020(minIteration, maxIteration, v);
  plotK0SPtRefoldedOverRec1020(minIteration, maxIteration, v);

  plotK0SPtUnfolded2030(minIteration, maxIteration, v);
  plotK0SPtUnfoldedOverGen2030(minIteration, maxIteration, v);
  plotK0SPtRefolded2030(minIteration, maxIteration, v);
  plotK0SPtRefoldedOverRec2030(minIteration, maxIteration, v);

  plotK0SPtUnfolded3040(minIteration, maxIteration, v);
  plotK0SPtUnfoldedOverGen3040(minIteration, maxIteration, v);
  plotK0SPtRefolded3040(minIteration, maxIteration, v);
  plotK0SPtRefoldedOverRec3040(minIteration, maxIteration, v);
}

// -------------------------------------------------------------------------------------------------
//
// K0S Z Plotters
//
// -------------------------------------------------------------------------------------------------

void plotK0SZUnfolded1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZ", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 1e-2, 10., mystrings::sZK0S, "Unfolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZUnfoldedOverGen1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Gen = GetHist<TH2>(inputs, rmutilities::testing::nameGenK0SZ);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Gen->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hGen = (TH1*)h2Gen->ProjectionX("hGen", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0.8, 1.2, mystrings::sZK0S, "#frac{Unfolded}{Generated}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.makeRatios(hGen);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZUnfoldedRelError1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z relative errors for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZRelError", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0., 0.5, mystrings::sZK0S, "Relative Error");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.addLatex(0.3, 0.8, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data());

  plotutils::Plotter q(p);
  q.resetHists();

  vector<TH1*> distributions = p.getHists();
  vector<TH1*> relErrors;
  for (int i = 0; i < p.getHists().size(); i++) {
    TH1* h = (TH1*)p.getHists()[i];
    TH1* r = (TH1*)h->Clone(TString::Format("%s_RelError", distributions[i]->GetName()).Data());
    r->Reset();
    for (int bin = 1; bin <= h->GetNbinsX(); bin++) {
      double bc = h->GetBinContent(bin);
      double be = h->GetBinError(bin);
      if (std::isnan(bc) || std::isnan(be))
        continue;
      if (bc < 1e-25)
        continue;
      r->SetBinContent(bin, be / bc);
    }
    plotutils::setStyle(r, i + 1);
    q.addHistogram(r);
  }
  q.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (int i = 0; i < p.getHists().size(); i++) {
      p.getHists()[i]->Print("all");
      q.getHists()[i]->Print("all");
    }
  }
}

void plotK0SZRefoldedRelError1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z relative errors for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZRelError", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0., 0.5, mystrings::sZK0S, "Relative Error");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.addLatex(0.3, 0.8, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data());

  plotutils::Plotter q(p);
  q.resetHists();

  vector<TH1*> distributions = p.getHists();
  vector<TH1*> relErrors;
  for (int i = 0; i < p.getHists().size(); i++) {
    TH1* h = (TH1*)p.getHists()[i];
    TH1* r = (TH1*)h->Clone(TString::Format("%s_RelError", distributions[i]->GetName()).Data());
    r->Reset();
    for (int bin = 1; bin <= h->GetNbinsX(); bin++) {
      double bc = h->GetBinContent(bin);
      double be = h->GetBinError(bin);
      if (std::isnan(bc) || std::isnan(be))
        continue;
      if (bc < 1e-25)
        continue;
      r->SetBinContent(bin, be / bc);
    }
    plotutils::setStyle(r, i + 1);
    q.addHistogram(r);
  }
  q.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (int i = 0; i < p.getHists().size(); i++) {
      p.getHists()[i]->Print("all");
      q.getHists()[i]->Print("all");
    }
  }
}

void plotK0SZUnfoldedRelError2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z relative errors for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZRelError", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0., 0.5, mystrings::sZK0S, "Relative Error");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.addLatex(0.3, 0.8, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data());

  plotutils::Plotter q(p);
  q.resetHists();

  vector<TH1*> distributions = p.getHists();
  vector<TH1*> relErrors;
  for (int i = 0; i < p.getHists().size(); i++) {
    TH1* h = (TH1*)p.getHists()[i];
    TH1* r = (TH1*)h->Clone(TString::Format("%s_RelError", distributions[i]->GetName()).Data());
    r->Reset();
    for (int bin = 1; bin <= h->GetNbinsX(); bin++) {
      double bc = h->GetBinContent(bin);
      double be = h->GetBinError(bin);
      if (std::isnan(bc) || std::isnan(be))
        continue;
      if (bc < 1e-25)
        continue;
      r->SetBinContent(bin, be / bc);
    }
    plotutils::setStyle(r, i + 1);
    q.addHistogram(r);
  }
  q.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (int i = 0; i < p.getHists().size(); i++) {
      p.getHists()[i]->Print("all");
      q.getHists()[i]->Print("all");
    }
  }
}

void plotK0SZRefoldedRelError2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z relative errors for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZRelError", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0., 0.5, mystrings::sZK0S, "Relative Error");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.addLatex(0.3, 0.8, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data());

  plotutils::Plotter q(p);
  q.resetHists();

  vector<TH1*> distributions = p.getHists();
  vector<TH1*> relErrors;
  for (int i = 0; i < p.getHists().size(); i++) {
    TH1* h = (TH1*)p.getHists()[i];
    TH1* r = (TH1*)h->Clone(TString::Format("%s_RelError", distributions[i]->GetName()).Data());
    r->Reset();
    for (int bin = 1; bin <= h->GetNbinsX(); bin++) {
      double bc = h->GetBinContent(bin);
      double be = h->GetBinError(bin);
      if (std::isnan(bc) || std::isnan(be))
        continue;
      if (bc < 1e-25)
        continue;
      r->SetBinContent(bin, be / bc);
    }
    plotutils::setStyle(r, i + 1);
    q.addHistogram(r);
  }
  q.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (int i = 0; i < p.getHists().size(); i++) {
      p.getHists()[i]->Print("all");
      q.getHists()[i]->Print("all");
    }
  }
}

void plotK0SZRefolded1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZ", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 1e-2, 10., mystrings::sZK0S, "Refolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZRefoldedOverRec1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Rec = GetHist<TH2>(inputs, rmutilities::testing::nameRecK0SZ);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Rec->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hRec = (TH1*)h2Rec->ProjectionX("hRec", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0.8, 1.2, mystrings::sZK0S, "#frac{Refolded}{Reconstructed}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.makeRatios(hRec);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZUnfolded1525(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(15., 25.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZ", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 1e-2, 10., mystrings::sZK0S, "Unfolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZUnfoldedOverGen1525(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(15., 25.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Gen = GetHist<TH2>(inputs, rmutilities::testing::nameGenK0SZ);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Gen->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hGen = (TH1*)h2Gen->ProjectionX("hGen", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0.8, 1.2, mystrings::sZK0S, "#frac{Unfolded}{Generated}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.makeRatios(hGen);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZRefolded1525(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(15., 25.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZ", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 1e-3, 10., mystrings::sZK0S, "Refolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZRefoldedOverRec1525(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(15., 25.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Rec = GetHist<TH2>(inputs, rmutilities::testing::nameRecK0SZ);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Rec->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hRec = (TH1*)h2Rec->ProjectionX("hRec", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0.8, 1.2, mystrings::sZK0S, "#frac{Refolded}{Reconstructed}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.makeRatios(hRec);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZUnfolded2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZ", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 1e-2, 10., mystrings::sZK0S, "Unfolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZUnfoldedOverGen2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Gen = GetHist<TH2>(inputs, rmutilities::testing::nameGenK0SZ);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Gen->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hGen = (TH1*)h2Gen->ProjectionX("hGen", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0.8, 1.2, mystrings::sZK0S, "#frac{Unfolded}{Generated}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.makeRatios(hGen);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZRefolded2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZ", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 1e-3, 10., mystrings::sZK0S, "Refolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZRefoldedOverRec2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(20., 30.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Rec = GetHist<TH2>(inputs, rmutilities::testing::nameRecK0SZ);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Rec->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hRec = (TH1*)h2Rec->ProjectionX("hRec", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0.8, 1.2, mystrings::sZK0S, "#frac{Refolded}{Reconstructed}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.makeRatios(hRec);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZUnfolded3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(30., 40.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZ", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 1e-2, 10., mystrings::sZK0S, "Unfolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZUnfoldedOverGen3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(30., 40.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Gen = GetHist<TH2>(inputs, rmutilities::testing::nameGenK0SZ);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Gen->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hGen = (TH1*)h2Gen->ProjectionX("hGen", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting unfolded K0S z ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("unfoldedK0SZRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0.8, 1.2, mystrings::sZK0S, "#frac{Unfolded}{Generated}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, isUnfolded, doZ);
  p.makeRatios(hGen);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZRefolded3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(30., 40.);
  inputs.inputFileName = _fileNameSingle;

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z spectra for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZ", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), true, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 1e-3, 10., mystrings::sZK0S, "Refolded");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotK0SZRefoldedOverRec3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(30., 40.);
  inputs.inputFileName = _fileNameSingle;

  TH2* h2Rec = GetHist<TH2>(inputs, rmutilities::testing::nameRecK0SZ);
  array<int, 2> ptjetBins = histutils::getProjectionBins(h2Rec->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1* hRec = (TH1*)h2Rec->ProjectionX("hRec", ptjetBins[0], ptjetBins[1]);

  const bool isUnfolded = true;
  const bool doZ = true;
  inputs.printLog(TString::Format("Plotting refolded K0S z ratios for ptjet = %.f - %.f", inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);
  string outputFileName = inputs.getNameFromPtJetProjection("refoldedK0SZRatio", "");
  plotutils::Plotter p(inputs.getNameFromIterations(outputFileName, ".pdf"), false, 0.04);
  p.makeFrame(1e-3, 1. + 1e-3, 0.8, 1.5, mystrings::sZK0S, "#frac{Refolded}{Reconstructed}");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  GetK0SIterations(inputs, p, !isUnfolded, doZ);
  p.makeRatios(hRec);
  p.setDrawOption("hist");
  p.plot();

  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    for (auto& hist : p.getHists())
      hist->Print("all");
  }
}

void plotallK0SZ(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kErrors) {
  gROOT->SetBatch(true);
  plotK0SZUnfolded1020(minIteration, maxIteration, v);
  plotK0SZUnfoldedOverGen1020(minIteration, maxIteration, v);
  plotK0SZRefolded1020(minIteration, maxIteration, v);
  plotK0SZRefoldedOverRec1020(minIteration, maxIteration, v);

  plotK0SZUnfolded1525(minIteration, maxIteration, v);
  plotK0SZUnfoldedOverGen1525(minIteration, maxIteration, v);
  plotK0SZRefolded1525(minIteration, maxIteration, v);
  plotK0SZRefoldedOverRec1525(minIteration, maxIteration, v);

  plotK0SZUnfolded2030(minIteration, maxIteration, v);
  plotK0SZUnfoldedOverGen2030(minIteration, maxIteration, v);
  plotK0SZRefolded2030(minIteration, maxIteration, v);
  plotK0SZRefoldedOverRec2030(minIteration, maxIteration, v);

  plotK0SZUnfolded3040(minIteration, maxIteration, v);
  plotK0SZUnfoldedOverGen3040(minIteration, maxIteration, v);
  plotK0SZRefolded3040(minIteration, maxIteration, v);
  plotK0SZRefoldedOverRec3040(minIteration, maxIteration, v);
}

// -------------------------------------------------------------------------------------------------
//
// Compare results with different ranges
//
// -------------------------------------------------------------------------------------------------

void GetJetIterationsForComp(InputSettings& inputs, plotutils::Plotter& plotter, vector<int> iterations, bool isUnfolded) {
  if (iterations.empty()) {
    cout << "GetJetIterationsForComp() Error: no iterations provided!" << endl;
    return;
  }
  if (!plotter.getLegend()) {
    cout << "GetJetIterationsForComp() Error: plotter does not have a legend!" << endl;
    return;
  }

  string histName = isUnfolded ? "unfoldedJets" : "refoldedJets";
  for (int i : iterations) {
    TH1* hist = GetHist<TH1>(inputs, histName + to_string(i));
    if (!hist)
      continue;

    plotutils::setStyle(hist, i);
    // plotter.addLegendEntry(hist, to_string(i));
    plotter.addLegendEntry(hist, TString::Format("%d (#it{p}_{T,jet}^{gen} %.f-%.f GeV/c, #it{p}_{T,jet}^{rec} %.f-%.f GeV/c)", i, inputs.ptjetminGen, inputs.ptjetmaxGen, inputs.ptjetminRec, inputs.ptjetmaxRec).Data());
    plotter.addHistogram(hist);
  }
}

void GetJetIterationsForComp(InputSettings& inputs, plotutils::Plotter& plotter, bool isUnfolded) {
  if (inputs.maxIteration < inputs.minIteration) {
    inputs.printLog("GetJetIterationsForComp() Error: invalid iteration range " + to_string(inputs.minIteration) + " to " + to_string(inputs.maxIteration), verbosityutilities::kErrors);
    return;
  }

  vector<int> iterations;
  for (int i = inputs.minIteration; i <= inputs.maxIteration; i++)
    iterations.push_back(i);

  GetJetIterationsForComp(inputs, plotter, iterations, isUnfolded);
}

void GetK0SIterationsForComp(InputSettings& inputs, plotutils::Plotter& plotter, vector<int> iterations, bool isUnfolded, bool doZ) {
  if (iterations.empty()) {
    inputs.printLog("GetK0SIterationsForComp() Error: no iterations provided!", verbosityutilities::kErrors);
    return;
  }
  if (!plotter.getLegend()) {
    inputs.printLog("GetK0SIterationsForComp() Error: plotter does not have a legend!", verbosityutilities::kErrors);
    return;
  }

  if (doZ)
    inputs.printLog("GetK0SIterationsForComp() getting Z histograms", verbosityutilities::kInfo);
  else
    inputs.printLog("GetK0SIterationsForComp() getting Pt histograms", verbosityutilities::kInfo);

  string baseName = isUnfolded ? rmutilities::testing::nameGenK0SPt : rmutilities::testing::nameRecK0SPt;
  string histName = isUnfolded ? "unfoldedK0SPt" : "refoldedK0SPt";
  if (doZ) {
    baseName = isUnfolded ? rmutilities::testing::nameGenK0SZ : rmutilities::testing::nameRecK0SZ;
    histName = isUnfolded ? "unfoldedK0SZ" : "refoldedK0SZ";
  }

  for (auto i : iterations) {
    TH2* hist = GetHist<TH2>(inputs, histName + to_string(i));
    if (!hist)
      continue;

    if (inputs.passVerbosityCheck(verbosityutilities::kDebug))
      hist->Print();

    array<int, 2> ptjetBins = histutils::getProjectionBins(hist->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
    TH1* hProj = (TH1*)hist->ProjectionX(inputs.getNameFromPtJetProjection(histName + to_string(i), "").c_str(), ptjetBins[0], ptjetBins[1]);
    plotutils::setStyle(hProj, i);
    // plotter.addLegendEntry(hProj, TString::Format("%d (#it{p}_{T,jet}^{gen} %.f-%.f GeV/c, #it{p}_{T,jet}^{rec} %.f-%.f GeV/c)", i, inputs.ptjetminGen, inputs.ptjetmaxGen, inputs.ptjetminRec, inputs.ptjetmaxRec).Data());
    plotter.addLegendEntry(hProj, TString::Format("%d, %s", i, RemovePrependingJunk(inputs.inputFileName, "ptjetgen", ".root").c_str()).Data());
    plotter.addHistogram(hProj);
  }
}

void GetK0SIterationsForComp(InputSettings& inputs, plotutils::Plotter& plotter, bool isUnfolded, bool doZ) {
  if (inputs.maxIteration < inputs.minIteration) {
    inputs.printLog("GetK0SIterationsForComp() Error: invalid iteration range " + to_string(inputs.minIteration) + " to " + to_string(inputs.maxIteration), verbosityutilities::kErrors);
    return;
  }

  vector<int> iterations;
  for (int i = inputs.minIteration; i <= inputs.maxIteration; i++)
    iterations.push_back(i);

  GetK0SIterationsForComp(inputs, plotter, iterations, isUnfolded, doZ);
}

void PlotJetsComparison(InputSettings& inputsWide, InputSettings& inputsNarrow, plotutils::Plotter& plotter, bool isUnfolded, bool doRatio) {
  GetJetIterationsForComp(inputsWide, plotter, isUnfolded);
  plotter.addLegendEntry((TObject*)0, ""); // Plotter struct does not like this!
  GetJetIterationsForComp(inputsNarrow, plotter, isUnfolded);
  plotter.setHistStyles();
  if (doRatio) {
    TH1* hErrors = (TH1*)plotter.getHists()[1]->Clone("hErrors");
    plotter.makeRatios(0);
    for (int i = 1; i <= hErrors->GetNbinsX(); i++) {
      double be = hErrors->GetBinError(i);
      double bc = hErrors->GetBinContent(i);
      double relError = 0;
      if (!std::isnan(be) || bc < 1e-25)
        relError = be / bc;

      plotter.getHists()[1]->SetBinError(i, relError);
    }
  }
  plotter.plot();

  if (inputsWide.passVerbosityCheck(verbosityutilities::kDebug)) {
    inputsWide.printLog("K0S hists (wide, then narrow):", verbosityutilities::kDebug);
    for (auto& hist : plotter.getHists())
      hist->Print("all");
  }
}

void PlotK0SComparison(InputSettings& inputsWide, InputSettings& inputsNarrow, plotutils::Plotter& plotter, bool isUnfolded, bool doZ, bool doRatio) {
  plotutils::SetPadMargins(0.22, 0.10, 0.15, 0.05);
  GetK0SIterationsForComp(inputsWide, plotter, isUnfolded, doZ);
  plotter.addLegendEntry((TObject*)0, ""); // Plotter struct does not like this!
  GetK0SIterationsForComp(inputsNarrow, plotter, isUnfolded, doZ);
  plotter.setHistStyles();
  if (doRatio) {
    TH1* hErrors = (TH1*)plotter.getHists()[1]->Clone("hErrors");
    plotter.makeRatios(0);
    for (int i = 1; i <= hErrors->GetNbinsX(); i++) {
      double be = hErrors->GetBinError(i);
      double bc = hErrors->GetBinContent(i);
      double relError = 0;
      if (!std::isnan(be) || bc < 1e-25)
        relError = be / bc;

      plotter.getHists()[1]->SetBinError(i, relError);
    }
  }
  plotter.plot();

  if (inputsWide.passVerbosityCheck(verbosityutilities::kDebug)) {
    inputsWide.printLog("K0S hists (wide, then narrow):", verbosityutilities::kDebug);
    for (auto& hist : plotter.getHists())
      hist->Print("all");
  }
}

void comparejetsRefolded(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  inputsWide.printLog("Comparing refolded jet spectra for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = min(inputsWide.ptjetminRec, inputsNarrow.ptjetminRec);
  double xMaxFrame = max(inputsWide.ptjetmaxRec, inputsNarrow.ptjetmaxRec);
  plotutils::Plotter pSpectra(inputsWide.getNameFromIterations("refoldedJetsComp", ".pdf"), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1., 1e9, mystrings::sPtJet, "Refolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  PlotJetsComparison(inputsWide, inputsNarrow, pSpectra, !isUnfolded, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(inputsWide.getNameFromIterations("refoldedJetsRatio", ".pdf"), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sPtJet, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  PlotJetsComparison(inputsWide, inputsNarrow, pRatio, !isUnfolded, doRatio);
}

void comparejetsUnfolded(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  inputsWide.printLog("Comparing unfolded jet spectra for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = min(inputsWide.ptjetminGen, inputsNarrow.ptjetminGen);
  double xMaxFrame = max(inputsWide.ptjetmaxGen, inputsNarrow.ptjetmaxGen);
  plotutils::Plotter pSpectra(inputsWide.getNameFromIterations("unfoldedJetsComp", ".pdf"), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1., 1e9, mystrings::sPtJet, "Unfolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  PlotJetsComparison(inputsWide, inputsNarrow, pSpectra, isUnfolded, !doRatio);

  // Ratio doesn't work due to different binning!
  // plotutils::Plotter pRatio(inputsWide.getNameFromIterations("unfoldedJetsRatio", ".pdf"), false, 0.04);
  // pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sPtJet, "Ratio");
  // pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  // PlotJetsComparison(inputsWide, inputsNarrow, pRatio, isUnfolded, doRatio);
}

void compareK0SZRefolded1015(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(10., 15.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(10., 15.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing refolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("refoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Refolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, !isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("refoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, !isUnfolded, doZ, doRatio);
}

void compareK0SZUnfolded1015(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(10., 15.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(10., 15.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing unfolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("unfoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Unfolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("unfoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, isUnfolded, doZ, doRatio);
}

void compareK0SZRefolded1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(10., 20.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(10., 20.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing refolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("refoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Refolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, !isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("refoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, !isUnfolded, doZ, doRatio);
}

void compareK0SZUnfolded1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(10., 20.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(10., 20.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing unfolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("unfoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Unfolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("unfoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, isUnfolded, doZ, doRatio);
}

void compareK0SZRefolded1520(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(15., 20.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(15., 20.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing refolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("refoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Refolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, !isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("refoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, !isUnfolded, doZ, doRatio);
}

void compareK0SZUnfolded1520(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(15., 20.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(15., 20.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing unfolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("unfoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Unfolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("unfoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, isUnfolded, doZ, doRatio);
}

void compareK0SZRefolded1525(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(15., 25.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(15., 25.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing refolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("refoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Refolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, !isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("refoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, !isUnfolded, doZ, doRatio);
}

void compareK0SZUnfolded1525(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(15., 25.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(15., 25.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing unfolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("unfoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Unfolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("unfoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, isUnfolded, doZ, doRatio);
}

void compareK0SZRefolded2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(20., 30.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(20., 30.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing refolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("refoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Refolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, !isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("refoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, !isUnfolded, doZ, doRatio);
}

void compareK0SZUnfolded2030(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(20., 30.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(20., 30.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing unfolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("unfoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Unfolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("unfoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, isUnfolded, doZ, doRatio);
}

void compareK0SZRefolded3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(30., 40.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(30., 40.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing refolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("refoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Refolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, !isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("refoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, !isUnfolded, doZ, doRatio);
}

void compareK0SZUnfolded3040(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.setPtJetGen(5., 80.);
  inputsWide.setPtJetRec(10., 60.);
  inputsWide.setPtJetProjection(30., 40.);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow;
  inputsNarrow.setVerbosity(v);
  inputsNarrow.setIterations(minIteration, maxIteration);
  inputsNarrow.setPtJetGen(10., 65.);
  inputsNarrow.setPtJetRec(10., 60.);
  inputsNarrow.setPtJetProjection(30., 40.);
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  const bool doRatio = true;
  const bool doZ = true;
  inputsWide.printLog("Comparing unfolded K0S z for different unfolding jet pt ranges", verbosityutilities::kInfo);

  double xMinFrame = 1e-3;
  double xMaxFrame = 1. + 1e-3;
  string outputFileNameSuffix = inputsWide.getNameFromIterations(inputsWide.getNameFromPtJetProjection("", ""), ".pdf");
  plotutils::Plotter pSpectra(("unfoldedK0SZComp" + outputFileNameSuffix).c_str(), true, 0.04);
  pSpectra.makeFrame(xMinFrame, xMaxFrame, 1e2, 1e6, mystrings::sZK0S, "Unfolded");
  pSpectra.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pSpectra.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pSpectra, isUnfolded, doZ, !doRatio);
  pSpectra.reset();

  plotutils::Plotter pRatio(("unfoldedK0SZRatio" + outputFileNameSuffix).c_str(), false, 0.04);
  pRatio.makeFrame(xMinFrame, xMaxFrame, 0.5, 1.5, mystrings::sZK0S, "Ratio");
  pRatio.makeLegend(0.25, 0.55, 0.70, 0.9, "Iteration (ranges)");
  pRatio.addLatex(0.3, 0.3, TString::Format("%.f < #it{p}_{T, jet} < %.f", inputsNarrow.ptjetminProjection, inputsNarrow.ptjetmaxProjection).Data());
  PlotK0SComparison(inputsWide, inputsNarrow, pRatio, isUnfolded, doZ, doRatio);
}

void compareK0SZ(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  gROOT->SetBatch(true);
  compareK0SZRefolded1015(minIteration, maxIteration, v);
  compareK0SZUnfolded1015(minIteration, maxIteration, v);
  compareK0SZRefolded1020(minIteration, maxIteration, v);
  compareK0SZUnfolded1020(minIteration, maxIteration, v);
  compareK0SZRefolded1520(minIteration, maxIteration, v);
  compareK0SZUnfolded1520(minIteration, maxIteration, v);
  compareK0SZRefolded1525(minIteration, maxIteration, v);
  compareK0SZUnfolded1525(minIteration, maxIteration, v);
  compareK0SZRefolded2030(minIteration, maxIteration, v);
  compareK0SZUnfolded2030(minIteration, maxIteration, v);
  compareK0SZRefolded3040(minIteration, maxIteration, v);
  compareK0SZUnfolded3040(minIteration, maxIteration, v);
}

void compareK0SZ(int iteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  compareK0SZ(iteration, iteration, v);
}

void DrawRelativeErrorsK0SZ(InputSettings& inputsNarrow, InputSettings& inputsWide, bool isUnfolded) {
  string histName = isUnfolded ? "unfoldedK0SZ" : "refoldedK0SZ";

  string outputFileNameSuffix = inputsNarrow.getNameFromIterations("Errors", ".pdf");
  plotutils::Plotter p((histName + outputFileNameSuffix).c_str(), false, 0.04);
  p.makeLegend(0.3, 0.5, 0.65, 0.9, "Iteration (#it{p}_{T, jet} ranges)");

  plotutils::Plotter q(p);

  vector<vector<double>> ptjetRanges = { {10., 20.}, {15., 25.}, {20., 30.}, {30., 40.} };
  for (int i = inputsNarrow.minIteration; i <= inputsNarrow.maxIteration; i++) {
    TH2* histNarrow = GetHist<TH2>(inputsNarrow, (histName + to_string(i)).c_str());
    if (histNarrow)
      histNarrow->SetName(TString::Format("%s_narrow", histNarrow->GetName()).Data());

    TH2* histWide = GetHist<TH2>(inputsWide, (histName + to_string(i)).c_str());
    if (histWide)
      histWide->SetName(TString::Format("%s_wide", histWide->GetName()).Data());

    if (!histWide || !histNarrow)
      continue;

    histNarrow->Print();
    histWide->Print();
    for (int iRange = 0; iRange < ptjetRanges.size(); iRange++) {
      array<int, 2>    ptjetBins     = histutils::getProjectionBins(histNarrow->GetYaxis(), ptjetRanges[iRange][0], ptjetRanges[iRange][1]);
      array<double, 2> ptjetBinEdges = histutils::getProjectionEdges(histNarrow->GetYaxis(), ptjetBins);

      string hName = TString::Format("%s_ptjet%.f-%.f_iter%d", histName.c_str(), ptjetBinEdges[0], ptjetBinEdges[1], i).Data();
      TH1* hProj   = (TH1*)histNarrow->ProjectionX(hName.c_str(), ptjetBins[0], ptjetBins[1]);
      TH1* hErrors = (TH1*)hProj->Clone(TString::Format("%s_errors", hName.c_str()).Data());
      hErrors->Reset();

      for (int iBin = 1; iBin <= hProj->GetNbinsX(); iBin++) {
        double be = hProj->GetBinError(iBin);
        double bc = hProj->GetBinContent(iBin);
        double relError = 0;
        if (!std::isnan(be) && bc > 1e-25)
          relError = be / bc;

        hErrors->SetBinContent(iBin, relError);
      }

      plotutils::setStyle(hErrors, iRange);
      p.addHistogram(hErrors);
      p.addLegendEntry(hErrors, TString::Format("It. %d (%.f-%.f GeV/c)", i, ptjetBinEdges[0], ptjetBinEdges[1]).Data(), "l");

      hProj->Print("all");
      hErrors->Print("all");

      ptjetBins       = histutils::getProjectionBins(histWide->GetYaxis(), ptjetRanges[iRange][0], ptjetRanges[iRange][1]);
      ptjetBinEdges   = histutils::getProjectionEdges(histWide->GetYaxis(), ptjetBins);
      TH1* hProjWide  = (TH1*)histWide->ProjectionX(TString::Format("%s_wide", hName.c_str()).Data(), ptjetBins[0], ptjetBins[1]);
      TH1* hRatioWide = (TH1*)hProj->Clone(TString::Format("%s_ratioWide", hName.c_str()).Data());
      hRatioWide->Divide(hProjWide);

      TH1* hErrorsWide = (TH1*)hRatioWide->Clone(TString::Format("%s_ptjet%.f-%.f_iter%d_wide_errors", histName.c_str(), ptjetBinEdges[0], ptjetBinEdges[1], i).Data());
      hErrorsWide->Reset();
      for (int iBin = 1; iBin <= hRatioWide->GetNbinsX(); iBin++) {
        hErrorsWide->SetBinContent(iBin, hRatioWide->GetBinError(iBin));
      }
      double alpha = 0.3;
      double lineWidth = 0;
      plotutils::setStyle(hErrorsWide, iRange, alpha, lineWidth);
      q.addHistogram(hErrorsWide);

      hProjWide->Print("all");
      hErrorsWide->Print("all");
    }
  }
  for (int i = inputsNarrow.minIteration; i <= inputsNarrow.maxIteration; i++) {
    TH2* histNarrow = GetHist<TH2>(inputsNarrow, (histName + to_string(i)).c_str());
    if (histNarrow)
      histNarrow->SetName(TString::Format("%s_narrow", histNarrow->GetName()).Data());

    TH2* histWide = GetHist<TH2>(inputsWide, (histName + to_string(i)).c_str());
    if (histWide)
      histWide->SetName(TString::Format("%s_wide", histWide->GetName()).Data());

    if (!histWide || !histNarrow)
      continue;

    histNarrow->Print();
    histWide->Print();
    for (int iRange = 0; iRange < ptjetRanges.size(); iRange++) {
      array<int, 2>    ptjetBins     = histutils::getProjectionBins(histNarrow->GetYaxis(), ptjetRanges[iRange][0], ptjetRanges[iRange][1]);
      array<double, 2> ptjetBinEdges = histutils::getProjectionEdges(histNarrow->GetYaxis(), ptjetBins);

      string hName = TString::Format("%s_ptjet%.f-%.f_iter%d", histName.c_str(), ptjetBinEdges[0], ptjetBinEdges[1], i).Data();
      TH1* hProj   = (TH1*)histNarrow->ProjectionX(hName.c_str(), ptjetBins[0], ptjetBins[1]);
      TH1* hErrors = (TH1*)hProj->Clone(TString::Format("%s_errors", hName.c_str()).Data());
      hErrors->Reset();

      for (int iBin = 1; iBin <= hProj->GetNbinsX(); iBin++) {
        double be = hProj->GetBinError(iBin);
        double bc = hProj->GetBinContent(iBin);
        double relError = 0;
        if (!std::isnan(be) && bc > 1e-25)
          relError = be / bc;

        hErrors->SetBinContent(iBin, relError);
      }

      plotutils::setStyle(hErrors, iRange);
      p.addHistogram(hErrors);
      p.addLegendEntry(hErrors, TString::Format("It. %d (%.f-%.f GeV/c)", i, ptjetBinEdges[0], ptjetBinEdges[1]).Data(), "l");

      hProj->Print("all");
      hErrors->Print("all");

      ptjetBins       = histutils::getProjectionBins(histWide->GetYaxis(), ptjetRanges[iRange][0], ptjetRanges[iRange][1]);
      ptjetBinEdges   = histutils::getProjectionEdges(histWide->GetYaxis(), ptjetBins);
      TH1* hProjWide  = (TH1*)histWide->ProjectionX(TString::Format("%s_wide", hName.c_str()).Data(), ptjetBins[0], ptjetBins[1]);
      TH1* hRatioWide = (TH1*)hProj->Clone(TString::Format("%s_ratioWide", hName.c_str()).Data());
      hRatioWide->Divide(hProjWide);

      TH1* hErrorsWide = (TH1*)hRatioWide->Clone(TString::Format("%s_ptjet%.f-%.f_iter%d_wide_errors", histName.c_str(), ptjetBinEdges[0], ptjetBinEdges[1], i).Data());
      hErrorsWide->Reset();
      for (int iBin = 1; iBin <= hRatioWide->GetNbinsX(); iBin++) {
        hErrorsWide->SetBinContent(iBin, hRatioWide->GetBinError(iBin));
      }
      double alpha = 0.3;
      double lineWidth = 0;
      plotutils::setStyle(hErrorsWide, iRange, alpha, lineWidth);
      q.addHistogram(hErrorsWide);

      hProjWide->Print("all");
      hErrorsWide->Print("all");
    }
  }


  p.makeFrame(mystrings::sZK0S, "Relative Error");
  p.getFrame()->SetMinimum(0.);
  p.getFrame()->SetMaximum(0.3);
  p.getFrame()->Draw();
  for (auto& hist : p.getHists())
    hist->Draw("same hist");

  for (auto& hist : q.getHists())
    hist->Draw("same hist p");

  p.getLegend()->Draw("same");
  for (auto& obj : p.getObjects())
    obj->Draw("same");

  p.getCanvas()->SaveAs(p.getOutputFileName().c_str());

  if (inputsNarrow.passVerbosityCheck(verbosityutilities::kDebug)) {
    inputsNarrow.printLog("K0S relative error hists:", verbosityutilities::kDebug);
    for (auto& hist : p.getHists()) {
      if (inputsNarrow.passVerbosityCheck(verbosityutilities::kDebugMax))
        hist->Print("all");
      else
        hist->Print();
    }
  }
}

void RelativeErrorsK0SZUnfolded(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow = inputsWide;
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  // DrawRelativeErrorsK0SZ(inputsWide, inputsNarrow, isUnfolded);
  DrawRelativeErrorsK0SZ(inputsNarrow, inputsWide, isUnfolded);
}

void RelativeErrorsK0SZRefolded(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(minIteration, maxIteration);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow = inputsWide;
  inputsNarrow.inputFileName = _fileNameNarrow;

  const bool isUnfolded = true;
  DrawRelativeErrorsK0SZ(inputsNarrow, inputsWide, !isUnfolded);
  // DrawRelativeErrorsK0SZ(inputsNarrow, inputsWide, !isUnfolded);
}

void RelativeErrorsK0SZ(int iteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  RelativeErrorsK0SZUnfolded(iteration, iteration, v);
  RelativeErrorsK0SZRefolded(iteration, iteration, v);
}

void RelativeErrorsK0SZAlt(int iteration, bool isUnfolded, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputsWide;
  inputsWide.setVerbosity(v);
  inputsWide.setIterations(iteration, iteration);
  inputsWide.inputFileName = _fileNameWide;

  InputSettings inputsNarrow = inputsWide;
  inputsNarrow.inputFileName = _fileNameNarrow;

  string histName = isUnfolded ? "unfoldedK0SZ1" : "refoldedK0SZ1";
  TH2* hWide = GetHist<TH2>(inputsWide, histName.c_str());
  hWide->SetName("hWide");
  TH2* hNarrow = GetHist<TH2>(inputsNarrow, histName.c_str());
  hNarrow->SetName("hNarrow");

  hWide->Print();
  hNarrow->Print();

  double jetptmin = 10., jetptmax = 20.;
  array<int, 2> ptjetBins = histutils::getProjectionBins(hWide->GetYaxis(), jetptmin, jetptmax);
  TH1* hProjWide = (TH1*)hWide->ProjectionX("hProjWide", ptjetBins[0], ptjetBins[1]);
  TH1* hProjNarrow = (TH1*)hNarrow->ProjectionX("hProjNarrow", ptjetBins[0], ptjetBins[1]);

  TH1* hRatioWide = (TH1*)hProjWide->Clone("hRatioWide");
  hRatioWide->Divide(hProjWide);

  TH1* hErrorsWide = (TH1*)hProjWide->Clone("hErrorsWide");
  hErrorsWide->Reset();
  TH1* hErrorsNarrow = (TH1*)hProjNarrow->Clone("hErrorsNarrow");
  hErrorsNarrow->Reset();

  for (int iBin = 1; iBin <= hProjNarrow->GetNbinsX(); iBin++) {
    double be = hProjNarrow->GetBinError(iBin);
    double bc = hProjNarrow->GetBinContent(iBin);
    double relError = 0;
    if (!std::isnan(be) && bc > 1e-25)
      hErrorsNarrow->SetBinContent(iBin, be/bc);

    hErrorsWide->SetBinContent(iBin, hRatioWide->GetBinError(iBin));
  }

  cout << "\nNarrow:\n";
  hProjNarrow->Print("all");
  hErrorsNarrow->Print("all");

  cout << "\nWide:\n";
  hProjWide->Print("all");
  hErrorsWide->Print("all");
}
