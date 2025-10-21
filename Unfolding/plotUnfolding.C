
#include "dounfolding.C"
#include "../plotting_macros/myStrings.C"

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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

void plotK0SZRefolded1020(int minIteration, int maxIteration, verbosityutilities::Verbosity v = verbosityutilities::kDebug) {
  InputSettings inputs;
  inputs.setVerbosity(v);
  inputs.setIterations(minIteration, maxIteration);
  inputs.setPtJetProjection(10., 20.);
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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
  inputs.inputFileName = "ClosureTest_527899_80.root";

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

