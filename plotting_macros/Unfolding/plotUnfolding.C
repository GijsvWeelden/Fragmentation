
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

#include "/Users/gijsvanweelden/Documents/Fragmentation/plotting_macros/histUtils.C"

void checkAllBinsExactlyOne(TH1* hist);
void checkBinning(std::vector<string> histNameVector);
void checkBinning(void);
template <typename T>
T loadHist(string fileName, string histName);
template <typename T>
T loadHist(string fileName, string dirName, string histName);
template <typename T>
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<T> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);
void printHists(void);
template <typename T>
void setStyle(T hist, int styleNumber);

std::vector<string> nameVector = {
  "./closureTest_projection.root"
};
std::vector<string> legendVector = {
  ""
};

void plotUnfoldedTrainingTruth(bool doPt = true)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName = "zUnfoldedOverTrainingTruth";
  string numeratorName = "testZUnfolded";
  string denominatorName = "trainingZTruth";
  string histTitle = "";
  string xTitle = "#it{z}";
  string yTitle = "#frac{Unfolded}{Training truth}";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  if (doPt) {
    xTitle = "#it{p}_{T}";
    saveName = "ptUnfoldedOverTrainingTruth";
    numeratorName = "testPtUnfolded";
    denominatorName = "trainingPtTruth";
  }

  // Plotting stuff
  bool setLogY = false;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2;
  if (doPt) { xMaxFrame = 300; }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1D*> histVector;
  // std::vector<string> nameVector;
  // std::vector<string> legendVector;

  // nameVector.push_back("./closureTest_pt10-300_ptmin10_z0-1_nIter1.root");
  // legendVector.push_back("");


  for (unsigned int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string legendEntry = legendVector[i];
    TH1D* numerator = loadHist<TH1D*>(fileName, numeratorName);
    TH1D* denominator = loadHist<TH1D*>(fileName, denominatorName);
    TH1D* hist = (TH1D*)numerator->Clone("ratio");
    hist->Divide(denominator);
    setStyle(hist, i);
    if (nameVector.size() > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
    histVector.push_back(hist);
    checkAllBinsExactlyOne(hist);
  }

  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

// -------------------- Refolded, Detector --------------------
void plotRefoldedAndDetector(string inName, int setting = 1, int itmin = 3, int itmax = 3,
                             double ptmin = 10., double ptmax = 300., string drawOption = "")
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;

  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = -2, yMaxFrame = 2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;

  std::vector<string> ptSettings = { "ptRefoldedOverDetector", "ptRefoldedMinusDetector", "ptDetectorDiffRatio" };
  std::vector<string> zSettings  = { "zRefoldedOverDetector",  "zRefoldedMinusDetector",  "zDetectorDiffRatio" };
  std::vector<string> yTitles    = { "#frac{Refolded}{Detector}",  "Refolded - Detector",  "#frac{Refolded - Detector}{Detector}" };
  yTitle = yTitles[abs(setting) - 1];
  yMinFrame = -2 * (abs(setting) != 1); // 0 for the ratio only

  if (setting < 0) {
    xTitle = "#it{p}_{T}";
    saveName = ptSettings[abs(setting) - 1];
    histName = ptSettings[abs(setting) - 1];
    xMaxFrame = 300.;
  }
  else if (setting > 0) {
    xTitle = "#it{z}";
    saveName = zSettings[abs(setting) - 1];
    histName = zSettings[abs(setting) - 1];
    xMinFrame = 1e-3;
    xMaxFrame = 1. + 1e-3;
    latexText = TString::Format("#it{p}_{T,jet}: %.0f - %.0f GeV/#it{c}", ptmin, ptmax).Data();
  }
  else {
    cout << "Error: invalid setting " << setting << endl;
    return;
  }

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1D*> histVector;
  int nIter = itmax + 1 - itmin;
  // cout << nIter << endl;
  for (int iter = itmin; iter <= itmax; iter++) {
    string legendEntry = std::to_string(iter);
    TH1D* hist = loadHist<TH1D*>(inName, TString::Format("%s_iter%d", histName.c_str(), iter).Data());
    setStyle(hist, iter);
    if (nIter > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
    histVector.push_back(hist);
  }
  saveName = TString::Format("%s_iter%d", saveName.c_str(), itmin);
  if (nIter > 1) { saveName = TString::Format("%s-%d", saveName.c_str(), itmax); }
  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, drawOption, latexText);
  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}
void plotRefoldedMinusDetector(string inName, bool doPt = true, int itmin = 3, int itmax = 3, string drawOption = "")
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName = "zRefoldedMinusDetector";
  string histName = "zRefoldedMinusDetector";
  string histTitle = "";
  string xTitle = "#it{z}";
  string yTitle = "Refolded - Detector";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  if (doPt) {
    xTitle = "#it{p}_{T}";
    saveName = "ptRefoldedMinusDetector";
    histName = "ptRefoldedMinusDetector";
  }

  // Plotting stuff
  bool setLogY = false;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = -2, yMaxFrame = 2;
  if (doPt) { xMaxFrame = 300; }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1D*> histVector;

  int nIter = itmax + 1 - itmin;
  cout << nIter << endl;
  for (int iter = itmin; iter <= itmax; iter++) {
    string legendEntry = std::to_string(iter);
    TH1D* hist = loadHist<TH1D*>(inName, TString::Format("%s_iter%d", histName.c_str(), iter).Data());
    setStyle(hist, iter);
    if (nIter > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
    histVector.push_back(hist);
  }

  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, drawOption, latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}
void plotRefoldedOverDetector(string inName, bool doPt = true, int itmin = 3, int itmax = 3, string drawOption = "")
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName = "zRefoldedOverDetector";
  string histName = "zRefoldedOverDetector";
  string histTitle = "";
  string xTitle = "#it{z}";
  string yTitle = "#frac{Refolded}{Detector}";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  if (doPt) {
    xTitle = "#it{p}_{T}";
    saveName = "ptRefoldedOverDetector";
    histName = "ptRefoldedOverDetector";
  }

  // Plotting stuff
  bool setLogY = false;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2;
  if (doPt) { xMaxFrame = 300; }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1D*> histVector;

  int nIter = itmax + 1 - itmin;
  cout << nIter << endl;
  for (int iter = itmin; iter <= itmax; iter++) {
    string legendEntry = std::to_string(iter);
    TH1D* hist = loadHist<TH1D*>(inName, TString::Format("%s_iter%d", histName.c_str(), iter).Data());
    setStyle(hist, iter);
    if (nIter > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
    histVector.push_back(hist);
  }

  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, drawOption, latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

// -------------------- Unfolded, Truth --------------------
void plotUnfoldedAndTruth(string inName, int setting = 1, int itmin = 3, int itmax = 3,
                          double ptmin = 10., double ptmax = 300., string drawOption = "")
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;

  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = -2, yMaxFrame = 2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;

  std::vector<string> ptSettings = { "ptUnfoldedOverTruth", "ptUnfoldedMinusTruth", "ptTruthDiffRatio" };
  std::vector<string> zSettings  = { "zUnfoldedOverTruth",  "zUnfoldedMinusTruth",  "zTruthDiffRatio" };
  std::vector<string> yTitles    = { "#frac{Unfolded}{Truth}",  "Unfolded - Truth",  "#frac{Unfolded - Truth}{Truth}" };
  yTitle = yTitles[abs(setting) - 1];
  yMinFrame = -2 * (abs(setting) != 1); // 0 for the ratio only

  if (setting < 0) {
    saveName = ptSettings[abs(setting) - 1];
    histName = ptSettings[abs(setting) - 1];
    xTitle = "#it{p}_{T}";
    xMaxFrame = 300.;
  }
  else if (setting > 0) {
    saveName = zSettings[abs(setting) - 1];
    histName = zSettings[abs(setting) - 1];
    xTitle = "#it{z}";
    xMinFrame = 1e-3;
    xMaxFrame = 1. + 1e-3;
    latexText = TString::Format("#it{p}_{T,jet}: %.0f - %.0f GeV/#it{c}", ptmin, ptmax).Data();
  }
  else {
    cout << "Error: invalid setting " << setting << endl;
    return;
  }

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1D*> histVector;
  int nIter = itmax + 1 - itmin;
  // cout << nIter << endl;
  for (int iter = itmin; iter <= itmax; iter++) {
    string legendEntry = std::to_string(iter);
    TH1D* hist = loadHist<TH1D*>(inName, TString::Format("%s_iter%d", histName.c_str(), iter).Data());
    setStyle(hist, iter);
    if (nIter > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
    histVector.push_back(hist);
  }
  saveName = TString::Format("%s_iter%d", saveName.c_str(), itmin).Data();
  if (nIter > 1) { saveName = TString::Format("%s-%d", saveName.c_str(), itmax).Data(); }
  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, drawOption, latexText);
  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}
void plotUnfoldedTruth(bool doPt = true, int iteration = 3)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string histName = "testZUnfoldedOverTruth";
  string saveName = "testZUnfoldedOverTruth";
  string histTitle = "";
  string xTitle = "#it{z}";
  string yTitle = "#frac{Unfolded}{Truth}";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;

  if (doPt) {
    xTitle = "#it{p}_{T}";
    histName = "testPtUnfoldedOverTruth";
    saveName = "testPtUnfoldedOverTruth";
  }

  // Plotting stuff
  bool setLogY = false;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2;
  if (doPt) { xMaxFrame = 300; }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1D*> histVector;
  // std::vector<string> nameVector;
  // std::vector<string> legendVector;

  // nameVector.push_back("./closureTest_pt10-300_ptmin10_z0-1_nIter1.root");
  // legendVector.push_back("");

  for (unsigned int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string legendEntry = legendVector[i];
    TH1D* hist = loadHist<TH1D*>(fileName, histName);
    setStyle(hist, i);
    if (nameVector.size() > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
    histVector.push_back(hist);
    checkAllBinsExactlyOne(hist);
  }

  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void compareMeasured(bool doPt = true)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string histName = "zMeasuredComparison";
  string saveName = "";
  string histTitle = "";
  string xTitle = "#it{z}";
  string yTitle = "#frac{Measured (projected)}{Measured (filled)}";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;

  if (doPt) {
    xTitle = "#it{p}_{T}";
    saveName = "ptMeasuredComparison";
  }

  // Plotting stuff
  bool setLogY = false;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2;
  if (doPt) { xMaxFrame = 300; }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  // std::vector<string> nameVector;
  // std::vector<string> legendVector;

  // nameVector.push_back("./RooUnfoldResponse.root");
  // legendVector.push_back("");

  for (int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string legendEntry = legendVector[i];
    TH1F* hDetectorProjected; TH1F* hMeasFullProjected;

    // Projected from RM
    TH2F* hDetector = loadHist<TH2F*>(fileName, "hDetector");
    // Filled together with ruResponse
    TH2F* h2MeasFull = loadHist<TH2F*>(fileName, "h2MeasFull");

    if (doPt) {
      hDetectorProjected = (TH1F*)hDetector->ProjectionX();
      hMeasFullProjected = (TH1F*)h2MeasFull->ProjectionX();
    } else {
      hDetectorProjected = (TH1F*)hDetector->ProjectionY();
      hMeasFullProjected = (TH1F*)h2MeasFull->ProjectionY();
    }
    TH1F* ratio = (TH1F*)hDetectorProjected->Clone("ratio");
    ratio->Divide(hMeasFullProjected);
    setStyle(ratio, i);
    if (nameVector.size() > 1) { legend->AddEntry(ratio, legendEntry.c_str()); }
    histVector.push_back(ratio);
    checkAllBinsExactlyOne(ratio);
  }
  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}
void compareTruth(bool doPt = true)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string histName = "zTruthComparison";
  string saveName = "";
  string histTitle = "";
  string xTitle = "#it{z}";
  string yTitle = "#frac{Truth (projected)}{Truth (filled)}";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;

  if (doPt) {
    xTitle = "#it{p}_{T}";
    saveName = "ptTruthComparison";
  }

  // Plotting stuff
  bool setLogY = false;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2;
  if (doPt) { xMaxFrame = 300; }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  // std::vector<string> nameVector;
  // std::vector<string> legendVector;

  // nameVector.push_back("./RooUnfoldResponse.root");
  // legendVector.push_back("");

  for (int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string legendEntry = legendVector[i];
    TH1F* hTruthProjected; TH1F* hTrueFullProjected;

    // Projected from RM
    TH2F* hTruth = loadHist<TH2F*>(fileName, "hTruth");
    // Filled together with ruResponse
    TH2F* h2TrueFull = loadHist<TH2F*>(fileName, "h2TrueFull");

    if (doPt) {
      hTruthProjected = (TH1F*)hTruth->ProjectionX();
      hTrueFullProjected = (TH1F*)h2TrueFull->ProjectionX();
    } else {
      hTruthProjected = (TH1F*)hTruth->ProjectionY();
      hTrueFullProjected = (TH1F*)h2TrueFull->ProjectionY();
    }
    TH1F* ratio = (TH1F*)hTruthProjected->Clone("ratio");
    ratio->Divide(hTrueFullProjected);
    setStyle(ratio, i);
    if (nameVector.size() > 1) { legend->AddEntry(ratio, legendEntry.c_str()); }
    histVector.push_back(ratio);
    checkAllBinsExactlyOne(ratio);
  }
  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void printHists(void)
{
  // string fileName = "./closureTest_pt10-300_ptmin10_z0-1_nIter1.root";
  string fileName = nameVector[0];
  std::vector<string> histNameVector = {
    "testZRefoldedOverDetector",
    "testPtRefoldedOverDetector",
    "testZUnfoldedOverTruth",
    "testPtUnfoldedOverTruth"
  };
  for (string histName : histNameVector) {
    TH1D* hist = loadHist<TH1D*>(fileName, histName);
    hist->Print("all");
  }
}
void checkBinning(std::vector<string> histNameVector)
{
  // string fileName = "./closureTest_pt10-300_ptmin10_z0-1_nIter1.root";
  string fileName = nameVector[0];
  for (string histName : histNameVector) {
    TH1D* hist = loadHist<TH1D*>(fileName, histName);
    cout << histName << " : "
         << hist->GetXaxis()->GetBinLowEdge(1) << " - "
         << hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX())
         << " (" << hist->GetNbinsX() << " bins, width "
         << hist->GetXaxis()->GetBinWidth(1) << ")" << endl;
  }
  cout << endl;
}
void checkBinning(void)
{
  std::vector<string> ptHistNameVector = {
    "trainingPtDetector",
    "trainingPtTruth",
    "testPtDetector",
    "testPtTruth",
    "testPtUnfolded",
    "testPtRefolded",
    "testPtRefoldedOverDetector",
    "testPtUnfoldedOverTruth"
  };
  checkBinning(ptHistNameVector);

  std::vector<string> zHistNameVector = {
    "trainingZDetector",
    "trainingZTruth",
    "testZDetector",
    "testZTruth",
    "testZUnfolded",
    "testZRefolded",
    "testZRefoldedOverDetector",
    "testZUnfoldedOverTruth"
  };
  checkBinning(zHistNameVector);
}
void checkAllBinsExactlyOne(TH1* hist)
{
  cout << "Checking bins for histogram " << hist->GetName() << endl;
  cout << "(";
  for (int iBin = 1; iBin <= hist->GetNbinsX(); iBin++) {
    cout << hist->GetBinContent(iBin);
    if (iBin < hist->GetNbinsX()) { cout << ", "; }
  }
  cout << ")" << endl;
}

template <typename T>
T loadHist(string fileName, string histName)
{
  TFile *inFile = TFile::Open(TString::Format("./%s", fileName.c_str()).Data());
  T hist = (T)inFile->Get(histName.c_str());
  return hist;
}
template <typename T>
T loadHist(string fileName, string dirName, string histName)
{
  TFile *inFile = TFile::Open(TString::Format("./%s", fileName.c_str()).Data());
  TDirectory* dir = (TDirectory*)inFile->Get(TString::Format("%s", dirName.c_str()).Data());
  T hist = (T)dir->Get(histName.c_str());
  return hist;
}
template <typename T>
void setStyle(T hist, int styleNumber)
{
  hist->SetLineWidth(3);
  hist->SetLineColor(GetColor(styleNumber));
  hist->SetMarkerStyle(GetMarker(styleNumber));
  hist->SetMarkerColor(GetColor(styleNumber));
}
template <typename T>
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<T> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText)
{
  canvas->cd();
  frame->Draw();
  string drawOption = "same";
  if (setDrawOption != "") {
    drawOption = TString::Format("%s %s", drawOption.c_str(), setDrawOption.c_str()).Data();
  }
  for (unsigned int i = 0; i < histVector.size(); i++) {
    T hist = histVector[i];
    hist->Draw(drawOption.c_str());
  }
  if (legend) { legend->Draw("same"); }
  if (latexText != "") { DrawLatex(0.3, 0.93, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
