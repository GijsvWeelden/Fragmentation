
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TString.h"
#include "TLegend.h"

#include "/Users/gijsvanweelden/Documents/Fragmentation/plotting_macros/histUtils.C"

int nDim        = 4;
int ptjetAxis   = 0;
int ptratioAxis = 1;
int ptv0Axis    = 2;
int obsAxis     = 3;

string formatHadronName(string hadron);
double getNjets(string inName, double jetptmin, double jetptmax);
void normaliseHistRowByRow(TH2D* hist);
void normaliseHistColByCol(TH2D* hist);
template <typename T>
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<T> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);
template <typename T>
void setStyle(T hist, int styleNumber);

// ----------------------------------------------------------

void plotZ(double ptmin = 10., double ptmax = 100., double v0min = 0., double v0max = 20.)
{
  const int nDim        = 4;
  const int ptjetAxis   = 0;
  const int ptratioAxis = 1;
  const int ptv0Axis    = 2;
  const int zAxis       = 3;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "#it{z}_{V0}";
  yTitle = "#frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{z}}";
  latexText = TString::Format("#splitline{ #it{p}_{T, jet} = %.0f - %.0f }{ #it{p}_{T, V0} = %.0f - %.0f }", ptmin, ptmax, v0min, v0max).Data();

  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 1e-4, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  int firstBinPt, lastBinPt, firstBinV0, lastBinV0;

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Load histograms
  string lfName = "v0studyWithDecays_pthat20-80.root";
  TFile* lfFile = TFile::Open(TString::Format("%s", lfName.c_str()).Data());
  TH1D* jetptUncorrected = (TH1D*)lfFile->Get("hjetpt");
  THnSparseD* zUncorrected = (THnSparseD*)lfFile->Get("hzU");
  TH1D* jetptCorrected = (TH1D*)lfFile->Get("hnewjetpt");
  THnSparseD* zCorrected = (THnSparseD*)lfFile->Get("hzC");

  string hfName = "v0studyWithoutDecays_pthat20-80.root";
  TFile* hfFile = TFile::Open(TString::Format("%s", hfName.c_str()).Data());
  TH1D* jetptNoDecays = (TH1D*)hfFile->Get("hjetpt");
  TH3D* zNoDecays = (TH3D*)hfFile->Get("hz");

  // Get N jets
  firstBinPt = 1;
  lastBinPt = jetptUncorrected->GetNbinsX();
  firstBinPt = jetptUncorrected->GetXaxis()->FindBin(ptmin + 1e-3);
  lastBinPt  = jetptUncorrected->GetXaxis()->FindBin(ptmax - 1e-3);
  double njetsUncorrected = jetptUncorrected->Integral(firstBinPt, lastBinPt);

  firstBinPt = 1;
  lastBinPt = jetptCorrected->GetNbinsX();
  firstBinPt = jetptCorrected->GetXaxis()->FindBin(ptmin + 1e-3);
  lastBinPt  = jetptCorrected->GetXaxis()->FindBin(ptmax - 1e-3);
  double njetsCorrected = jetptCorrected->Integral(firstBinPt, lastBinPt);

  firstBinPt = 1;
  lastBinPt = jetptNoDecays->GetNbinsX();
  firstBinPt = jetptNoDecays->GetXaxis()->FindBin(ptmin + 1e-3);
  lastBinPt  = jetptNoDecays->GetXaxis()->FindBin(ptmax - 1e-3);
  double njetsNoDecays = jetptNoDecays->Integral(firstBinPt, lastBinPt);

  // Project onto z
  firstBinPt = 1;
  lastBinPt = zUncorrected->GetAxis(ptjetAxis)->GetNbins();
  firstBinV0 = 1;
  lastBinV0 = zUncorrected->GetAxis(ptv0Axis)->GetNbins();
  firstBinPt = zUncorrected->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt = zUncorrected->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  firstBinV0 = zUncorrected->GetAxis(ptv0Axis)->FindBin(v0min + 1e-3);
  lastBinV0 = zUncorrected->GetAxis(ptv0Axis)->FindBin(v0max + 1e-3);
  zUncorrected->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  zUncorrected->GetAxis(ptv0Axis)->SetRange(firstBinV0, lastBinV0);
  TH1D* zU = (TH1D*)zUncorrected->Projection(zAxis);
  zU->SetName("zU");
  zU->Rebin(rebinNumber);
  zU->Scale(1./njetsUncorrected);
  setStyle(zU, 0);
  legend->AddEntry(zU, "#it{z}_{V0} (V0 decays, uncorrected #it{p}_{T, jet})");
  histVector.push_back(zU);

  firstBinPt = 1;
  lastBinPt = zCorrected->GetAxis(ptjetAxis)->GetNbins();
  firstBinV0 = 1;
  lastBinV0 = zCorrected->GetAxis(ptv0Axis)->GetNbins();
  firstBinPt = zCorrected->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt = zCorrected->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  firstBinV0 = zCorrected->GetAxis(ptv0Axis)->FindBin(v0min + 1e-3);
  lastBinV0 = zCorrected->GetAxis(ptv0Axis)->FindBin(v0max - 1e-3);
  zCorrected->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  zCorrected->GetAxis(ptv0Axis)->SetRange(firstBinV0, lastBinV0);
  TH1D* zC = (TH1D*)zCorrected->Projection(zAxis);
  zC->SetName("zC");
  zC->Rebin(rebinNumber);
  zC->Scale(1./njetsCorrected);
  setStyle(zC, 1);
  legend->AddEntry(zC, "#it{z}_{V0} (V0 decays, corrected #it{p}_{T, jet})");
  histVector.push_back(zC);

  firstBinPt = 1;
  lastBinPt = zNoDecays->GetNbinsX();
  firstBinV0 = 1;
  lastBinV0 = zNoDecays->GetNbinsY();
  firstBinPt = zNoDecays->GetXaxis()->FindBin(ptmin + 1e-3);
  lastBinPt = zNoDecays->GetXaxis()->FindBin(ptmax - 1e-3);
  firstBinV0 = zNoDecays->GetYaxis()->FindBin(v0min + 1e-3);
  lastBinV0 = zNoDecays->GetYaxis()->FindBin(v0max - 1e-3);
  TH1D* zND = (TH1D*)zNoDecays->ProjectionZ("zND", firstBinPt, lastBinPt, firstBinV0, lastBinV0);
  zND->Rebin(rebinNumber);
  zND->Scale(1./njetsNoDecays);
  setStyle(zND, 2);
  legend->AddEntry(zND, "#it{z}_{V0} (No V0 decays)");
  histVector.push_back(zND);

  saveName = "z_LF-HF";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0min, v0max);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotZRatio(double ptmin = 10., double ptmax = 100., double v0min = 0., double v0max = 20.)
{
  const int nDim        = 4;
  const int ptjetAxis   = 0;
  const int ptratioAxis = 1;
  const int ptv0Axis    = 2;
  const int zAxis       = 3;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "#it{z}_{V0}";
  yTitle = "V0 decays / No decays";
  latexText = TString::Format("#splitline{ #it{p}_{T, jet} = %.0f - %.0f }{ #it{p}_{T, V0} = %.0f - %.0f }", ptmin, ptmax, v0min, v0max).Data();

  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  int firstBinPt, lastBinPt, firstBinV0, lastBinV0;

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Load histograms
  string lfName = "v0studyWithDecays_pthat20-80.root";
  TFile* lfFile = TFile::Open(TString::Format("%s", lfName.c_str()).Data());
  TH1D* jetptUncorrected = (TH1D*)lfFile->Get("hjetpt");
  THnSparseD* zUncorrected = (THnSparseD*)lfFile->Get("hzU");
  TH1D* jetptCorrected = (TH1D*)lfFile->Get("hnewjetpt");
  THnSparseD* zCorrected = (THnSparseD*)lfFile->Get("hzC");

  string hfName = "v0studyWithoutDecays_pthat20-80.root";
  TFile* hfFile = TFile::Open(TString::Format("%s", hfName.c_str()).Data());
  TH1D* jetptNoDecays = (TH1D*)hfFile->Get("hjetpt");
  TH3D* zNoDecays = (TH3D*)hfFile->Get("hz");

  // Get N jets
  firstBinPt = 1;
  lastBinPt = jetptUncorrected->GetNbinsX();
  firstBinPt = jetptUncorrected->GetXaxis()->FindBin(ptmin + 1e-3);
  lastBinPt  = jetptUncorrected->GetXaxis()->FindBin(ptmax - 1e-3);
  double njetsUncorrected = jetptUncorrected->Integral(firstBinPt, lastBinPt);

  firstBinPt = 1;
  lastBinPt = jetptCorrected->GetNbinsX();
  firstBinPt = jetptCorrected->GetXaxis()->FindBin(ptmin + 1e-3);
  lastBinPt  = jetptCorrected->GetXaxis()->FindBin(ptmax - 1e-3);
  double njetsCorrected = jetptCorrected->Integral(firstBinPt, lastBinPt);

  firstBinPt = 1;
  lastBinPt = jetptNoDecays->GetNbinsX();
  firstBinPt = jetptNoDecays->GetXaxis()->FindBin(ptmin + 1e-3);
  lastBinPt  = jetptNoDecays->GetXaxis()->FindBin(ptmax - 1e-3);
  double njetsNoDecays = jetptNoDecays->Integral(firstBinPt, lastBinPt);

  // Project onto z
  firstBinPt = 1;
  lastBinPt = zUncorrected->GetAxis(ptjetAxis)->GetNbins();
  firstBinV0 = 1;
  lastBinV0 = zUncorrected->GetAxis(ptv0Axis)->GetNbins();
  firstBinPt = zUncorrected->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt = zUncorrected->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  firstBinV0 = zUncorrected->GetAxis(ptv0Axis)->FindBin(v0min + 1e-3);
  lastBinV0 = zUncorrected->GetAxis(ptv0Axis)->FindBin(v0max - 1e-3);
  zUncorrected->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  zUncorrected->GetAxis(ptv0Axis)->SetRange(firstBinV0, lastBinV0);
  TH1D* zU = (TH1D*)zUncorrected->Projection(zAxis);
  zU->SetName("zU");
  zU->Rebin(rebinNumber);
  zU->Scale(1./njetsUncorrected);
  setStyle(zU, 0);
  // legend->AddEntry(zU, "#it{z}_{V0} (V0 decays, uncorrected #it{p}_{T, jet})");
  // histVector.push_back(zU);

  firstBinPt = 1;
  lastBinPt = zCorrected->GetAxis(ptjetAxis)->GetNbins();
  firstBinV0 = 1;
  lastBinV0 = zCorrected->GetAxis(ptv0Axis)->GetNbins();
  firstBinPt = zCorrected->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt = zCorrected->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  firstBinV0 = zCorrected->GetAxis(ptv0Axis)->FindBin(v0min + 1e-3);
  lastBinV0 = zCorrected->GetAxis(ptv0Axis)->FindBin(v0max - 1e-3);
  zCorrected->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  zCorrected->GetAxis(ptv0Axis)->SetRange(firstBinV0, lastBinV0);
  TH1D* zC = (TH1D*)zCorrected->Projection(zAxis);
  zC->SetName("zC");
  zC->Rebin(rebinNumber);
  zC->Scale(1./njetsCorrected);
  setStyle(zC, 1);
  // legend->AddEntry(zC, "#it{z}_{V0} (V0 decays, corrected #it{p}_{T, jet})");
  // histVector.push_back(zC);

  firstBinPt = 1;
  lastBinPt = zNoDecays->GetNbinsX();
  firstBinV0 = 1;
  lastBinV0 = zNoDecays->GetNbinsY();
  firstBinPt = zNoDecays->GetXaxis()->FindBin(ptmin + 1e-3);
  lastBinPt = zNoDecays->GetXaxis()->FindBin(ptmax - 1e-3);
  firstBinV0 = zNoDecays->GetYaxis()->FindBin(v0min + 1e-3);
  lastBinV0 = zNoDecays->GetYaxis()->FindBin(v0max - 1e-3);
  TH1D* zND = (TH1D*)zNoDecays->ProjectionZ("zND", firstBinPt, lastBinPt, firstBinV0, lastBinV0);
  zND->Rebin(rebinNumber);
  zND->Scale(1./njetsNoDecays);
  setStyle(zND, 2);
  // legend->AddEntry(zND, "#it{z}_{V0} (No V0 decays)");
  // histVector.push_back(zND);

  // Calculate ratios
  TH1D* rU = (TH1D*)zU->Clone("rU");
  rU->Divide(zND);
  setStyle(rU, 1);
  legend->AddEntry(rU, "Uncorrected jet #it{p}_{T}");
  histVector.push_back(rU);

  TH1D* rC = (TH1D*)zC->Clone("rC");
  rC->Divide(zND);
  setStyle(rC, 2);
  legend->AddEntry(rC, "Corrected jet #it{p}_{T}");
  histVector.push_back(rC);

  TLine* line = new TLine(xMinFrame, 1., xMaxFrame, 1.);
  line->SetLineWidth(3);
  line->SetLineColor(GetColor(0));
  line->SetLineStyle(9);

  saveName = "zRatio_LF-HF";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0min, v0max);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  frame->Draw();
  rU->Draw("same");
  rC->Draw("same");
  line->Draw("same");
  legend->Draw("same");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

// Formats the hadron name to look nice (Greek letters, sub- and superscripts)
string formatHadronName(string hadron)
{
  string had = hadron;
  if (hadron == "pi"){
    had = "#pi^{#pm}";
  }
  else if (hadron == "pi0"){
    had = "#pi^{0}";
  }
  else if (hadron == "K0L"){
    had = "K^{0}_{L}";
  }
  else if (hadron == "K0S"){
    had = "K^{0}_{S}";
  }
  else if (hadron == "K0"){
    had = "#it{K}^{0}";
  }
  else if (hadron == "K"){
    had = "K^{#pm}";
  }
  else if (hadron == "Lambda0"){
    had = "#it{#Lambda}^{0}";
  }
  return had;
}
double getNjets(TFile* file, double jetptmin, double jetptmax)
{
  TH1D* hjetpt = (TH1D*)file->Get("hjetpt");
  int firstBinPt = 1, lastBinPt = hjetpt->GetNbinsX();
  firstBinPt = hjetpt->FindBin(jetptmin);
  lastBinPt = hjetpt->FindBin(jetptmax);
  double integral = hjetpt->Integral(firstBinPt, lastBinPt);
  return integral;
}
// Normalise 2D histogram row-by-row
void normaliseHistRowByRow(TH2D* hist)
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
void normaliseHistColByCol(TH2D* hist)
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
template <typename T>
void setStyle(T hist, int styleNumber)
{
  hist->SetLineWidth(3);
  hist->SetLineColor(GetColor(styleNumber));
  hist->SetMarkerStyle(GetMarker(styleNumber));
  hist->SetMarkerColor(GetColor(styleNumber));
}
