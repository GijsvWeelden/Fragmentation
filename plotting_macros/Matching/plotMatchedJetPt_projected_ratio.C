
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

void plotMatchedJetPt_projected_ratio(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName1 = "./AnalysisResults-dist20-wEtaCut.root";
  string inName2 = "./AnalysisResults-dist20-noEtaCut.root";
  string saveDir = ".";
  TFile *inFile1 = TFile::Open(TString::Format("./%s", inName1.c_str()).Data());
  TFile *inFile2 = TFile::Open(TString::Format("./%s", inName2.c_str()).Data());
  if(!inFile1){
    std::cout << "File " << inFile1 << " not found. Aborting program." << std::endl;
    return;
  }
  if(!inFile2){
    std::cout << "File " << inFile2 << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* dir1 = (TDirectory*)inFile1->Get("jet-fragmentation/matching/jets");
  TDirectory* dir2 = (TDirectory*)inFile2->Get("jet-fragmentation/matching/jets");

  double ptMin = 40, ptMax = 60;
  // double ptMin = 100, ptMax = 140;
  double xMinFrame = 0, xMaxFrame = 100, yMinFrame = 0.8, yMaxFrame = 1.2, zMinFrame = 1e-5, zMaxFrame = 1;

  string histName = "matchDetJetPtPartJetPt";
  string histTitle = "";
  string saveName = "matchedJetPt_projected";
  string xTitle = "#it{p}_{T}^{jet, det}";
  string yTitle = "normalised count";
  string legendTitle = "";
  string latexText = TString::Format("#splitline{Matching distance: 0.2}{#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  bool setLogZ = false;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  int rebinNumber = 1;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  TH2F* matchedJetPt1 = (TH2F*)dir1->Get(histName.c_str());
  matchedJetPt1->Sumw2();
  TH2F* matchedJetPt2 = (TH2F*)dir2->Get(histName.c_str());
  matchedJetPt2->Sumw2();
  int firstBinPtTruth = 1, lastBinPtTruth = matchedJetPt1->GetNbinsY() + 1;

  firstBinPtTruth = matchedJetPt1->GetYaxis()->FindBin(ptMin);
  lastBinPtTruth = matchedJetPt1->GetYaxis()->FindBin(ptMax);
  TH1F* matchedJetPt1_projected = (TH1F*)matchedJetPt1->ProjectionX(TString::Format("matchedJetPt1_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  matchedJetPt1_projected->SetLineWidth(3);
  matchedJetPt1_projected->SetLineColor(GetColor(0));
  matchedJetPt1_projected->SetMarkerStyle(GetMarker(0));
  matchedJetPt1_projected->SetMarkerColor(GetColor(0));
  TH1F* ratioHist = (TH1F*)matchedJetPt1_projected->Clone("ratio");

  firstBinPtTruth = matchedJetPt2->GetYaxis()->FindBin(ptMin);
  lastBinPtTruth = matchedJetPt2->GetYaxis()->FindBin(ptMax);
  TH1F* matchedJetPt2_projected = (TH1F*)matchedJetPt2->ProjectionX(TString::Format("matchedJetPt2_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);

  ratioHist->Divide(matchedJetPt2_projected);
  frame->GetYaxis()->SetTitle("#frac{cut before matching}{cut after matching}");
  histVector.push_back(ratioHist);

  saveName = TString::Format("%s_ratio", saveName.c_str()).Data();
  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "hist p", latexText);

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
  TLine* line = new TLine(0, 1, 100, 1);
  line->SetLineStyle(9);
  line->Draw("same");
  if (legend) { legend->Draw("same"); }
  if (latexText != "") { DrawLatex(0.3, 0.8, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
