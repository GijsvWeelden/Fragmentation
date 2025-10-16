
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

void plotjetsrec() {}

void plotjetsgen(int minIteration, int maxIteration) {
  InputSettings inputs; inputs.setVerbosity(verbosityutilities::kDebug);
  inputs.setIterations(minIteration, maxIteration);

  inputs.inputFileName = "ClosureTest_520161_40.root";
  inputs.outputFileName = TString::Format("unfoldedJets_iter%d-%d.pdf", minIteration, maxIteration).Data();

  TH1* hGen = GetHist<TH1>(inputs, rmutilities::testing::nameGenJets);
  hGen->Scale(1. / hGen->Integral(), "width");

  vector<TH1*> unfHists;
  plotutils::Plotter p(inputs.outputFileName);
  p.makeLegend(0.6, 0.9, 0.7, 0.9, "");

  for (int i = minIteration; i <= maxIteration; i++) {
    TH1* h = GetHist<TH1>(inputs, "unfoldedJets" + to_string(i));
    if (!h)
      continue;

    h->Scale(1. / h->Integral(), "width");
    plotutils::setStyle(h, i);
    // h = (TH1*)histutils::divideWithProtection(h, hGen);
    unfHists.push_back(h);

    p.addLegendEntry(h, to_string(i));
    p.addHistogram(h);
    h->Print("all");
  }
  // unfHists[0]->Draw("hist");
  // p.setHists(unfHists);
  p.makeRatios(hGen);

  // p.makeFrame(5., 80., 0.8, 1.2, mystrings::sPtJet, "Unfolded / Generated");
  p.setDrawOption("hist");
  p.plot();
}

