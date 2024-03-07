
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

void plotDetJetPhi(double ptMin = 40, double ptMax = 60, string input = "LHC23d4/train111677")
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName = TString::Format("../../data/%s.root", input.c_str()).Data();
  string saveDir = TString::Format("../../Plots/%s", input.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* fragDir = (TDirectory*)inFile->Get("jet-fragmentation");
  TDirectory* matchDir = (TDirectory*)fragDir->Get("matching");
  TDirectory* matchJetsDir = (TDirectory*)matchDir->Get("jets");
  TDirectory* matchTracksDir = (TDirectory*)matchDir->Get("tracks");

  string histName = "matchDetJetPtEtaPhi";
  string histTitle = "";
  string saveName = "detJetPhi";
  string xTitle = "#phi^{jet, det}";
  string yTitle = "normalised count";
  string legendTitle = "";
  string latexText = "latexText";
  double textSize = 0.03;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  bool setLogZ = false;
  double xMinFrame = 0, xMaxFrame = 2*TMath::Pi(), yMinFrame = 0, yMaxFrame = 2e-1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  int rebinNumber = 8;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  if (setLogZ) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  TH3F* matchedJetPtEtaPhi = (TH3F*)matchJetsDir->Get(histName.c_str());
  matchedJetPtEtaPhi->Sumw2();
  int firstBinPt = 1, lastBinPt = matchedJetPtEtaPhi->GetNbinsX();
  int firstBinEta = 1, lastBinEta = matchedJetPtEtaPhi->GetNbinsY();

  firstBinPt = matchedJetPtEtaPhi->GetXaxis()->FindBin(ptMin);
  lastBinPt = matchedJetPtEtaPhi->GetXaxis()->FindBin(ptMax);
  TH1F* matchedJetPhi = (TH1F*)matchedJetPtEtaPhi->ProjectionZ(TString::Format("eta_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPt, lastBinPt, firstBinEta, lastBinEta);
  matchedJetPhi->Rebin(rebinNumber);
  matchedJetPhi->Scale(1./matchedJetPhi->Integral());
  matchedJetPhi->SetLineWidth(3);
  matchedJetPhi->SetLineColor(GetColor(0));
  matchedJetPhi->SetMarkerStyle(GetMarker(0));
  matchedJetPhi->SetMarkerColor(GetColor(0));
  histVector.push_back(matchedJetPhi);

  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax).Data();
  saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  latexText = TString::Format("#splitline{PYTHIA, LHC23d4}{#splitline{13.6 TeV pp 500 kHz}{#splitline{anti-kt jets, #it{R} = 0.%d}{#it{p}_{T}^{jet}: %.0f-%.0f GeV}}}", R, ptMin, ptMax).Data();
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
