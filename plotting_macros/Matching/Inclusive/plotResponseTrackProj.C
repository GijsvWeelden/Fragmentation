
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

void plotResponseTrackProj(void)
{
  double time = clock();
  string inName = "./AnalysisResults-sample1.root";
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
  string xTitle = "#it{z}_{ch}^{ det.}";
  string yTitle = "#it{z}_{ch}^{ part.}";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.05;
  bool setLogY = false;
  bool setLogZ = true;
  bool centerTitle = false;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 1, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  int rebinNumber = 2;
  // double ptMin = 10, ptMax = 300;
  double ptMin = 60, ptMax = 80;

  // Plotting stuff
  gStyle->SetNdivisions(505, "xy");
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 800, 800);
  myCanvas->SetLeftMargin(0.11);
  myCanvas->SetBottomMargin(0.11);
  myCanvas->SetRightMargin(0.125);
  myCanvas->SetTopMargin(0.02);
  if (setLogY) { myCanvas->SetLogy(); }
  if (setLogZ) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle, 0.1, false);
  frame->GetXaxis()->SetLabelSize(labelSize);
  frame->GetYaxis()->SetLabelSize(labelSize);
  frame->GetXaxis()->SetTitleSize(titleSize);
  frame->GetYaxis()->SetTitleSize(titleSize);
  frame->GetXaxis()->CenterTitle(centerTitle);
  frame->GetYaxis()->CenterTitle(centerTitle);
  frame->GetXaxis()->SetTitleOffset(1.);
  frame->GetYaxis()->SetTitleOffset(1.);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  THnF* response = (THnF*)matchJetsDir->Get(histName.c_str());
  int firstBinPtTruth = 1, lastBinPtTruth = response->GetAxis(3)->GetNbins();
  int projectionAxisX = 1, projectionAxisY = 3;

  int ptDetAxis = 0;
  int zDetAxis = 1;
  int ptTruthAxis = 2;

  firstBinPtTruth = response->GetAxis(ptTruthAxis)->FindBin(ptMin);
  lastBinPtTruth = response->GetAxis(ptTruthAxis)->FindBin(ptMax);
  response->GetAxis(ptTruthAxis)->SetRangeUser(ptMin, ptMax);
  // response->GetAxis(ptTruthAxis)->SetRange(firstBinPtTruth, lastBinPtTruth);

  TH2F* zTzD = (TH2F*)response->Projection(projectionAxisY, projectionAxisX);
  zTzD->Rebin2D(rebinNumber, rebinNumber);
  normaliseHistRowByRow(zTzD);
  zTzD->GetZaxis()->SetRangeUser(zMinFrame, zMaxFrame);

  saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMin, ptMax).Data();
  saveName = TString::Format("%s_binsize%d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s_textsize%d", saveName.c_str(), (int)std::round(100*textSize));
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  latexText = TString::Format("#splitline{ALICE simulation, ideal alignment}{#splitline{pp #sqrt{#it{s} } = 13.6 TeV}{#splitline{Ch-particle jets, anti-#it{k}_{T}, #it{R} = 0.%d}{#splitline{%.0f < #it{p}_{T, ch. jet}^{ part.} < %.0f GeV/#it{c}}{|#it{#eta}_{jet}^{det.}| < 0.5}}}}", R, ptMin, ptMax).Data();
  // latexText = TString::Format("#splitline{ALICE performance, ideal alignment}{#splitline{500 kHz pp #sqrt{#it{s} } = 13.6 TeV}{#splitline{Ch-particle jets, anti-#it{k}_{T}, #it{R} = 0.%d}{#splitline{%.0f < #it{p}_{T, ch. jet}^{ part.} < %.0f GeV/#it{c}}{|#it{#eta}_{jet}^{det.}| < 0.5}}}}", R, ptMin, ptMax).Data();
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
  if (latexText != "") { DrawLatex(0.15, 0.86, latexText.c_str(), legend->GetTextSize()); }
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
