
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

void normaliseHistRowByRow(TH2F* hist);
void normaliseHistColByCol(TH2F* hist);

void plotOneHist(TCanvas* canvas, TH1F* frame, TH2F* hist, TLegend* legend, string saveName, string setDrawOption, string latexText);
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<TH1F*> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);

void plotMatchedJetPt_projected_matchDistComparison(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName05 = "./AnalysisResults-dist05-noEtaCut.root";
  string inName10 = "./AnalysisResults-dist10-noEtaCut.root";
  string inName20 = "./AnalysisResults-dist20-noEtaCut.root";
  string saveDir = ".";
  TFile *inFile05 = TFile::Open(TString::Format("./%s", inName05.c_str()).Data());
  TFile *inFile10 = TFile::Open(TString::Format("./%s", inName10.c_str()).Data());
  TFile *inFile20 = TFile::Open(TString::Format("./%s", inName20.c_str()).Data());
  if(!inFile05){
    std::cout << "File " << inFile05 << " not found. Aborting program." << std::endl;
    return;
  }
  if(!inFile10){
    std::cout << "File " << inFile10 << " not found. Aborting program." << std::endl;
    return;
  }
  if(!inFile20){
    std::cout << "File " << inFile20 << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* dir05 = (TDirectory*)inFile05->Get("jet-fragmentation/matching/jets");
  TDirectory* dir10 = (TDirectory*)inFile10->Get("jet-fragmentation/matching/jets");
  TDirectory* dir20 = (TDirectory*)inFile20->Get("jet-fragmentation/matching/jets");

  // double ptMin = 40, ptMax = 60;
  double ptMin = 100, ptMax = 140;

  string histName = "matchDetJetPtPartJetPt";
  string histTitle = "";
  string saveName = "matchedJetPt_projected_matchDistComparison";
  string xTitle = "#it{p}_{T}^{jet, det}";
  string yTitle = "normalised count";
  string legendTitle = "";
  string latexText = TString::Format("#splitline{Cut after matching}{#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  bool setLogZ = true;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 0, yMaxFrame = .3, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
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
  TH2F* matchedJetPt05 = (TH2F*)dir05->Get(histName.c_str());
  matchedJetPt05->Sumw2();
  TH2F* matchedJetPt10 = (TH2F*)dir10->Get(histName.c_str());
  matchedJetPt10->Sumw2();
  TH2F* matchedJetPt20 = (TH2F*)dir20->Get(histName.c_str());
  matchedJetPt20->Sumw2();

  int firstBinPtTruth = 1, lastBinPtTruth = matchedJetPt05->GetNbinsY() + 1;
  firstBinPtTruth = matchedJetPt05->GetYaxis()->FindBin(ptMin);
  lastBinPtTruth = matchedJetPt05->GetYaxis()->FindBin(ptMax);

  TH1F* matchedJetPt05_projected = (TH1F*)matchedJetPt05->ProjectionX(TString::Format("matchedJetPt05_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  legend->AddEntry(matchedJetPt05_projected, "match dist 0.05");
  histVector.push_back(matchedJetPt05_projected);

  TH1F* matchedJetPt10_projected = (TH1F*)matchedJetPt10->ProjectionX(TString::Format("matchedJetPt10_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  legend->AddEntry(matchedJetPt10_projected, "match dist 0.1");
  histVector.push_back(matchedJetPt10_projected);

  TH1F* matchedJetPt20_projected = (TH1F*)matchedJetPt20->ProjectionX(TString::Format("matchedJetPt20_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  legend->AddEntry(matchedJetPt20_projected, "match dist 0.2");
  histVector.push_back(matchedJetPt20_projected);

  for (int i = 0; i < histVector.size(); i++) {
    auto hist = histVector[i];
    hist->SetLineWidth(3);
    hist->SetMarkerSize(2);
    hist->SetLineColor(GetColor(i));
    hist->SetMarkerColor(GetColor(i));
    hist->SetMarkerStyle(GetMarker(i));
    hist->Scale(1./hist->Integral());
  }

  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
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
  if (latexText != "") { DrawLatex(0.3, 0.8, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
