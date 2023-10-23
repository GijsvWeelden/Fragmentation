
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

#include "../histUtils.C"

void normaliseHistRowByRow(TH2F* hist);
void normaliseHistColByCol(TH2F* hist);

void plotOneHist(TCanvas* canvas, TH1F* frame, TH2F* hist, TLegend* legend, string saveName, string setDrawOption, string latexText);
template <typename T>
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<T> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);

template <typename T>
T loadHist(string fileName, string histName);
template <typename T>
T loadHist(string fileName, string dirName, string histName);

template <typename T>
void setStyle(T hist, int styleNumber);

void plotRefoldedDetector(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  bool doPt = true;

  string saveName = "testZRefoldedOverDetector";
  string histName = "testZRefoldedOverDetector";
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
    saveName = "testPtRefoldedOverDetector";
    histName = "testPtRefoldedOverDetector";
  }

  // Plotting stuff
  bool setLogY = false;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2;
  if (doPt) { xMaxFrame = 200; }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1D*> histVector;
  std::vector<string> nameVector;
  std::vector<string> legendVector;

  // TODO: Do something smart with nIter here
  nameVector.push_back("./closureTest_binWidth1_pt0-100_ptmin20_z0-1_nIter1.root");
  legendVector.push_back("1");

  for (unsigned int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string legendEntry = legendVector[i];
    TH1D* hist = loadHist<TH1D*>(fileName, histName);
    setStyle(hist, i);
    if (nameVector.size() > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
    histVector.push_back(hist);
  }

  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void plotUnfoldedTruth(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  bool doPt = true;

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
    xTitle = "#it{P}_{T}";
    histName = "testPtUnfoldedOverTruth";
    saveName = "testPtUnfoldedOverTruth";
  }

  // Plotting stuff
  bool setLogY = false;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2;
  if (doPt) { xMaxFrame = 200; }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1D*> histVector;
  std::vector<string> nameVector;
  std::vector<string> legendVector;

  // TODO: Do something smart with nIter here
  nameVector.push_back("./closureTest_binWidth1_pt0-100_ptmin20_z0-1_nIter1.root");
  legendVector.push_back("1");

  for (unsigned int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string legendEntry = legendVector[i];
    TH1D* hist = loadHist<TH1D*>(fileName, histName);
    setStyle(hist, i);
    if (nameVector.size() > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
    histVector.push_back(hist);
  }

  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
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
  if (latexText != "") { DrawLatex(0.3, 0.9, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
