
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
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<TH1F*> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);

template <typename T>
T loadMatchDistHist(string fileName, string dirName = "jet-fragmentation/matching/jets", string histName = "matchPartJetPtEtaPhiMatchDist");

template <typename T>
void setTHnRange(T hist, int axis, double minValue, double maxValue);

template <typename T>
void setStyle(T hist, int styleNumber);

void plotMatchingDistanceComparison(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  bool compareMaxDists = false;
  int normaliseHistograms = 1;
  // double ptMin = 40, ptMax = 60;
  double ptMin = 100, ptMax = 140;
  int maxMatchingDistance = 5;
  int etaTrackSize = 1;
  bool doRatio = true;
  if (doRatio) { normaliseHistograms = 0; }

  string histName = "matchPartJetPtEtaPhiMatchDist";
  string histTitle = "";
  string saveName = "matchDist";
  string saveDir = ".";
  string xTitle = "matching distance";
  string yTitle = "weighted count";
  string legendTitle = "";
  string latexText = "latexText";

  if (normaliseHistograms > 0) { yTitle = "normalised count"; }
  if (normaliseHistograms < 0) { yTitle = "count relative to first bin"; }
  if (doRatio) { yTitle = ""; }
  latexText = TString::Format("#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}", ptMin, ptMax).Data();

  double textSize = 0.03;
  double labelSize = 0.04;
  double titleSize = 0.04;//0.05;
  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 0.4, yMinFrame = 1e-5, yMaxFrame = 2, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int firstBinPt = 1, lastBinPt = 200;
  int rebinNumber = 1;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  if (!compareMaxDists) { xMaxFrame = maxMatchingDistance * 1e-2; }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  std::vector<string> nameVector;
  std::vector<string> legendVector;

  string latexAddition;
  if (compareMaxDists) {
    string etaString = "smallEta";
    latexAddition = "|#eta_{track}^{part}| < 0.9";
    if (etaTrackSize == 2) {
      etaString = "largeEta";
      latexAddition = "|#eta_{track}^{part}| < 0.2";
    }
    latexText = TString::Format("#splitline{%s}{%s}", latexAddition.c_str(), latexText.c_str()).Data();
    nameVector.push_back(TString::Format("AnalysisResults-matchDist05-%s.root", etaString.c_str()).Data());
    legendVector.push_back("Max matching distance 0.05");
    nameVector.push_back(TString::Format("AnalysisResults-matchDist20-%s.root", etaString.c_str()).Data());
    legendVector.push_back("Max matching distance 0.2");
  }
  else {
    latexAddition = "Max matching distance 0.05";
    if (maxMatchingDistance == 20) {
      latexAddition = "Max matching distance 0.2";
    }
    latexText = TString::Format("#splitline{%s}{%s}", latexAddition.c_str(), latexText.c_str()).Data();
    nameVector.push_back(TString::Format("AnalysisResults-matchDist%02d-smallEta.root", maxMatchingDistance).Data());
    legendVector.push_back("|#eta_{track}^{part}| < 0.9");
    nameVector.push_back(TString::Format("AnalysisResults-matchDist%02d-largeEta.root", maxMatchingDistance).Data());
    legendVector.push_back("|#eta_{track}^{part}| < 2.0");
  }


  int ptAxis = 0, etaAxis = 1, phiAxis = 2, matchDistAxis = 3;
  for (unsigned int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string histLegend = legendVector[i];
    THnSparseF* matchPartJetPtEtaPhiMatchDist = loadMatchDistHist<THnSparseF*>(fileName);
    setTHnRange(matchPartJetPtEtaPhiMatchDist, ptAxis, ptMin, ptMax);
    TH1F* hist = (TH1F*)matchPartJetPtEtaPhiMatchDist->Projection(matchDistAxis);
    hist->Rebin(rebinNumber);
    setStyle(hist, i);
    if (normaliseHistograms > 0) { hist->Scale(1./hist->Integral()); }
    else if (normaliseHistograms < 0) { hist->Scale(1./hist->GetBinContent(1)); }
    legend->AddEntry(hist, histLegend.c_str());
    histVector.push_back(hist);
  }
  if (doRatio) {
    TH1F* base = (TH1F*)histVector[0]->Clone("base");
    for (auto hist : histVector) {
      hist->Divide(base);
    }
  }

  if (compareMaxDists) { saveName = TString::Format("%sComparison", saveName.c_str()); }
  else { saveName = TString::Format("%s_trackEtaComparison", saveName.c_str()); }
  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  if (normaliseHistograms > 0) { saveName = TString::Format("%s_normalised", saveName.c_str()); }
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

// Normalise 2D histogram row-by-row
void normaliseHistRowByRow(TH2F* hist)
{
  int firstColBin = 1, lastColBin = hist->GetNbinsX();
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY();
  for (int iRow = 1; iRow <= lastRowBin; iRow++) {
    double integral = hist->Integral(firstColBin, lastColBin, iRow, iRow);
    if (integral < 1) { continue; }
    for (int iCol = 1; iCol <= lastColBin; iCol++) {
      double binContent = hist->GetBinContent(iCol, iRow);
      binContent /= integral;
      hist->SetBinContent(iCol, iRow, binContent);
    }
  }
}
// Normalise 2D histogram col-by-col
void normaliseHistColByCol(TH2F* hist)
{
  int firstColBin = 1, lastColBin = hist->GetNbinsX();
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY();
  for (int iCol = 1; iCol <= lastColBin; iCol++) {
    double integral = hist->Integral(iCol, iCol, firstRowBin, lastRowBin);
    if (integral < 1) { continue; }
    for (int iRow = 1; iRow <= lastRowBin; iRow++) {
      double binContent = hist->GetBinContent(iCol, iRow);
      binContent /= integral;
      hist->SetBinContent(iCol, iRow, binContent);
    }
  }
}

void plotOneHist(TCanvas* canvas, TH1F* frame, TH2F* hist, TLegend* legend, string saveName, string setDrawOption, string latexText)
{
  canvas->cd();
  frame->Draw();
  string drawOption = "same colz";
  if (setDrawOption != "") {
    drawOption = TString::Format("%s %s", drawOption.c_str(), setDrawOption.c_str()).Data();
  }
  hist->Draw(drawOption.c_str());
  if (legend) { legend->Draw("same"); }
  if (latexText != "") { DrawLatex(0.3, 0.9, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<TH1F*> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText)
{
  canvas->cd();
  frame->Draw();
  string drawOption = "same";
  if (setDrawOption != "") {
    drawOption = TString::Format("%s %s", drawOption.c_str(), setDrawOption.c_str()).Data();
  }
  for (unsigned int i = 0; i < histVector.size(); i++) {
    TH1F* hist = histVector[i];
    hist->Draw(drawOption.c_str());
  }
  if (legend) { legend->Draw("same"); }
  if (latexText != "") { DrawLatex(0.3, 0.9, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
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
void setTHnRange(T hist, int axis, double minValue, double maxValue)
{
  int firstBin = 1, lastBin = 1;
  firstBin = hist->GetAxis(axis)->FindBin(minValue + 1e-5);
  lastBin = hist->GetAxis(axis)->FindBin(maxValue - 1e-5);
  hist->GetAxis(axis)->SetRange(firstBin, lastBin);
}

template <typename T>
T loadMatchDistHist(string fileName, string dirName, string histName)
{
  TFile *inFile = TFile::Open(TString::Format("./%s", fileName.c_str()).Data());
  TDirectory* dir = (TDirectory*)inFile->Get(TString::Format("%s", dirName.c_str()).Data());
  T hist = (T)dir->Get(histName.c_str());
  // hist->Sumw2();
  return hist;
}
