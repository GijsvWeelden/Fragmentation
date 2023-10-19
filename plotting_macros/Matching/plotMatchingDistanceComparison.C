
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

TH1F* getMatchDistHist(string fileName, double ptMin, double ptMax)
{
  string histName = "matchPartJetPtMatchDist";
  TFile *inFile = TFile::Open(TString::Format("./%s", fileName.c_str()).Data());
  TDirectory* dir = (TDirectory*)inFile->Get("jet-fragmentation/matching/jets");
  TH2F* matchPartJetPtMatchDist = (TH2F*)dir->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = matchPartJetPtMatchDist->GetNbinsX();
  TH1F* matchingDistance = (TH1F*)matchPartJetPtMatchDist->ProjectionY(TString::Format("%s_pt%.0f-%.0f", fileName.c_str(), ptMin, ptMax).Data(), firstBinPt, lastBinPt);

  return matchingDistance;
}

template <typename T>
void setStyle(T hist, int styleNumber)
{
  hist->SetLineWidth(3);
  hist->SetLineColor(GetColor(styleNumber));
  hist->SetMarkerStyle(GetMarker(styleNumber));
  hist->SetMarkerColor(GetColor(styleNumber));
}

void plotMatchingDistanceComparison(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  bool compareMaxDists = false;
  bool normaliseHistograms = true;
  double ptMin = 40, ptMax = 60;
  // double ptMin = 100, ptMax = 140;

  string histName = "matchPartJetPtMatchDist";
  string histTitle = "";
  string saveName = "matchDist";
  string saveDir = ".";
  string xTitle = "matching distance";
  string yTitle = "count";
  string legendTitle = "";
  string latexText = "latexText";

  if (normaliseHistograms) { yTitle = "normalised count"; }
  if (compareMaxDists) {
    latexText = TString::Format("#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}", ptMin, ptMax).Data();
  }
  else {
    latexText = TString::Format("#splitline{Max matching distance 0.05}{#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}}", ptMin, ptMax).Data();
  }

  double textSize = 0.03;
  double labelSize = 0.04;
  double titleSize = 0.04;//0.05;
  bool setLogY = true;
  bool setLogZ = false;
  double xMinFrame = 0, xMaxFrame = 0.4, yMinFrame = 1e-5, yMaxFrame = 1, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int firstBinPt = 1, lastBinPt = 200;
  int R = 4;
  int rebinNumber = 1;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  if (setLogZ) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  std::vector<string> nameVector;
  std::vector<string> legendVector;

  if (compareMaxDists) {
    nameVector.push_back("AnalysisResults-dist05-noEtaCut.root");
    legendVector.push_back("Max matching distance 0.05");
    nameVector.push_back("AnalysisResults-dist10-noEtaCut.root");
    legendVector.push_back("Max matching distance 0.1");
    nameVector.push_back("AnalysisResults-dist20-noEtaCut.root");
    legendVector.push_back("Max matching distance 0.2");
    nameVector.push_back("AnalysisResults-dist40-noEtaCut.root");
    legendVector.push_back("Max matching distance 0.4");
  }
  else {
    nameVector.push_back("AnalysisResults-dist05-noEtaCut.root");
    legendVector.push_back("|#eta_{track}^{part}| < 0.9");
    nameVector.push_back("AnalysisResults-dist05-noEtaCut-largerEtaTrack.root");
    legendVector.push_back("|#eta_{track}^{part}| < 2.0");
  }

  for (unsigned int i = 0; i < nameVector.size(); i++) {
    string histName = nameVector[i];
    string histLegend = legendVector[i];
    TH1F* hist = getMatchDistHist(histName, ptMin, ptMax);
    hist->Rebin(rebinNumber);
    setStyle(hist, i);
    hist->Scale(1./hist->Integral());
    legend->AddEntry(hist, histLegend.c_str());
    histVector.push_back(hist);
  }

  // saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  if (normaliseHistograms) { saveName = TString::Format("%s_normalised", saveName.c_str()); }
  if (!compareMaxDists) { saveName = TString::Format("%s_trackCut", saveName.c_str()); }
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
