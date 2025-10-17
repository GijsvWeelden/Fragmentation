
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

// Want:
// * Proper plotter for closure tests
// * Plot multiple iterations on one canvas
void FillJetPlotter(InputSettings& inputs, plotutils::Plotter& plotter, bool isUnfolded) {
  if (inputs.maxIteration < inputs.minIteration) {
    inputs.printLog("FillJetPlotter() Error: invalid iteration range " + to_string(inputs.minIteration) + " to " + to_string(inputs.maxIteration), verbosityutilities::kErrors);
    return;
  }
  if (!plotter.getLegend()) {
    inputs.printLog("FillJetPlotter() Error: plotter does not have a legend!", verbosityutilities::kErrors);
    return;
  }

  string baseName = isUnfolded ? rmutilities::testing::nameGenJets : rmutilities::testing::nameRecJets;
  string histName = isUnfolded ? "unfoldedJets" : "refoldedJets";

  TH1* hBase = GetHist<TH1>(inputs, baseName);
  hBase->Scale(1. / hBase->Integral(), "width");

  for (int i = inputs.minIteration; i <= inputs.maxIteration; i++) {
    TH1* hist = GetHist<TH1>(inputs, histName + to_string(i));
    if (!hist)
      continue;

    hist->Scale(1. / hist->Integral(), "width");
    plotutils::setStyle(hist, i);
    plotter.addLegendEntry(hist, to_string(i));
    plotter.addHistogram(hist);
    if (inputs.passVerbosityCheck(verbosityutilities::kDebug))
      hist->Print("all");
  }
  plotter.makeRatios(hBase);
}

void plotjets(int minIteration, int maxIteration, bool isUnfolded) {
  InputSettings inputs; inputs.setVerbosity(verbosityutilities::kDebug);
  inputs.setIterations(minIteration, maxIteration);
  inputs.inputFileName = "ClosureTest_520161_40.root";

  string outputPrefix = "refoldedJets";
  double xMinFrame = 10., xMaxFrame = 60.;
  string yTitle = "#frac{Refolded}{Reconstructed}";
  if (isUnfolded) {
    xMinFrame = 5.;
    xMaxFrame = 80.;
    yTitle = "#frac{Unfolded}{Generated}";
    outputPrefix = "unfoldedJets";
  }

  inputs.outputFileName = outputPrefix + "_iter" + to_string(minIteration) + "-" + to_string(maxIteration) + ".pdf";
  plotutils::Plotter p(inputs.outputFileName, false, 0.04);
  p.makeFrame(xMinFrame, xMaxFrame, 0.8, 1.2, mystrings::sPtJet, yTitle);
  p.setDrawOption("hist");
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "Iteration");
  FillJetPlotter(inputs, p, isUnfolded);
  p.plot();
}

