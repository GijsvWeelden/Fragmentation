
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

void plotResponseChargeFrag(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName = "../../data/LHC21k6/train109274.root";
  // string saveDir = "../../Plots/LHC21k6/test109274";
  string saveDir = ".";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* fragDir = (TDirectory*)inFile->Get("jet-fragmentation");
  TDirectory* matchDir = (TDirectory*)fragDir->Get("matching");
  TDirectory* matchJetsDir = (TDirectory*)matchDir->Get("jets");
  TDirectory* matchTracksDir = (TDirectory*)matchDir->Get("tracks");

  string histName = "matchDetJetPtTrackProjPartJetPtTrackProj";
  string histTitle = "";
  string saveName = "responseZproj";
  string xTitle = "#it{z}^{ det}";
  string yTitle = "#it{z}^{ truth}";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.05;
  bool setLogY = false;
  bool setLogZ = true;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 1, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  int rebinNumber = 2;
  double ptMin = 60, ptMax = 80;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 800, 800);
  myCanvas->SetLeftMargin(0.15);
  myCanvas->SetBottomMargin(0.15);
  if (setLogY) { myCanvas->SetLogy(); }
  if (setLogZ) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  THnF* response = (THnF*)matchJetsDir->Get(histName.c_str());
  int firstBinPtTruth = 1, lastBinPtTruth = response->GetAxis(3)->GetNbins();
  int projectionAxisX = 1, projectionAxisY = 3;
  int ptTruthAxis = 2;

  firstBinPtTruth = response->GetAxis(ptTruthAxis)->FindBin(ptMin);
  lastBinPtTruth = response->GetAxis(ptTruthAxis)->FindBin(ptMax);
  response->GetAxis(ptTruthAxis)->SetRange(firstBinPtTruth, lastBinPtTruth);

  TH2F* zTzD = (TH2F*)response->Projection(projectionAxisY, projectionAxisX);
  zTzD->Rebin2D(rebinNumber, rebinNumber);
  normaliseHistRowByRow(zTzD);
  zTzD->GetZaxis()->SetRangeUser(zMinFrame, zMaxFrame);

  saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMin, ptMax).Data();
  saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  latexText = TString::Format("#splitline{PYTHIA, ideal alignment}{#splitline{13.6 TeV pp 500 kHz}{#splitline{anti-kt jets, #it{R} = 0.%d}{#it{p}_{T}^{ch. jet, truth}: %.0f-%.0f GeV/#it{c}}}}", R, ptMin, ptMax).Data();
  plotOneHist(myCanvas, frame, zTzD, legend, saveName, "", latexText);

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
  if (latexText != "") { DrawLatex(0.3, 0.8, latexText.c_str(), legend->GetTextSize()); }
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
