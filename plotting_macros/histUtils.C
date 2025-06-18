//
// This macro contains a set of useful functions for plotting histograms
// and performing standard operations
//

#include <vector>
#include <iostream>
#include <typeinfo>
#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "TString.h"

#include "TLegend.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TROOT.h"

#include "plotUtils.C"

#ifndef HISTUTILS_H
#define HISTUTILS_H

// GLOBAL SETTINGS !!!
const double MassK0S = 0.497611;
const double MassLambda = 1.115683;
const double MassLambda0 = MassLambda;
gStyle->SetNdivisions(505, "xy");

// Helpful functions
string getDataSet(int train)
{
  switch (train) {
    case 210373: return "LHC24b1b";
    case 252064: return "LHC22o_pass6";
    case 271952: return "LHC24b1b";
    case 280432: return "LHC24g4";
    case 282430: return "LHC22o_pass7_small";
    case 350079: return "LHC22o_pass7_small";
    case 349872: return "LHC22o_pass7_small";
    case 349871: return "LHC24g4";
    case 419138: return "LHC22o_pass7_small";
    case 417810: return "LHC22o_pass7_small";
    case 419996: return "LHC22o_pass7_medium";
    case 420215: return "LHC22o_pass7_small";
    default:     return "Could not find dataset";
  }
}

template <typename T>
T loadHist(TFile* inFile, string histName)
{
  T hist = (T)inFile->Get(TString::Format("%s", histName.c_str()).Data());
  return hist;
}

template <typename T>
T loadHist(TDirectory* inDir, string histName)
{
  T hist = (T)inDir->Get(TString::Format("%s", histName.c_str()).Data());
  return hist;
}

// 1D histogram version
void plotOneHist(TH1F* hist, string xTitle, string yTitle, string histTitle, string legendTitle, string saveName,
                 double xMinFrame, double xMaxFrame, double yMinFrame, double yMaxFrame,
                 double xMinLegend, double xMaxLegend, double yMinLegend, double yMaxLegend,
                 bool setLogY, bool drawLegend)
{
  auto nc = new TCanvas("Plot", "Plot", 600, 600);
  nc->SetLeftMargin(0.15);
  nc->cd();

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(TString::Format("%s", histTitle.c_str()).Data());
  frame->Draw();
  hist->SetLineColor(GetColor(0));
  hist->SetMarkerColor(GetColor(0));
  hist->SetMarkerStyle(GetMarker(0));
  hist->SetLineWidth(3);
  hist->Draw("same, hist");

  if (drawLegend) {
    auto legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->AddEntry(hist, TString::Format("%s", histTitle.c_str()).Data());
  }

  if (setLogY) { nc->SetLogy(); }
  nc->SaveAs(TString::Format("./%s.pdf", saveName.c_str()).Data());
  delete nc;
}
// 2D histogram version
void plotOneHist(TH2F* hist, string xTitle, string yTitle, string histTitle, string saveName,
                 double xMinFrame, double xMaxFrame, double yMinFrame, double yMaxFrame,
                 double xMinLegend, double xMaxLegend, double yMinLegend, double yMaxLegend,
                 bool setLogZ, bool drawLegend)
{
  auto nc = new TCanvas("Plot", "Plot", 600, 600);
  nc->SetLeftMargin(0.15);
  nc->cd();

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(TString::Format("%s", histTitle.c_str()).Data());
  frame->GetYaxis()->SetTitleOffset(1.3);
  frame->Draw();
  hist->Draw("same, colz");

  if (drawLegend) {
    auto legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    legend->SetBorderSize(0);
    legend->AddEntry(hist, TString::Format("%s", histTitle.c_str()).Data());
  }

  if (setLogZ) { nc->SetLogz(); }
  nc->SaveAs(TString::Format("./%s.pdf", saveName.c_str()).Data());
  delete nc;
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
  if (latexText != "") { DrawLatex(0.4, 0.8, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
template <typename T>
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<T> histVector, TLegend* legend, TLatex* latex, string saveName, string setDrawOption)
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
  if (latex) { latex->Draw("same"); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
// Plot n histograms, but set all colours, markers, etc. manually
template <typename T>
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<T> histVector, TLegend* legend, vector<TLatex*> latexVector, string setDrawOption = "")
{
  canvas->cd();
  frame->Draw();
  string drawOption = "same";
  if (setDrawOption != "") {
    drawOption += " " + setDrawOption;
  }
  for (unsigned int i = 0; i < histVector.size(); i++) {
    T hist = histVector[i];
    hist->Draw(drawOption.c_str());
  }
  if (legend) { legend->Draw("same"); }
  for (auto latex : latexVector) {
    latex->Draw("same");
  }
  canvas->SaveAs(canvas->GetName());
}
// Plot only histograms
void plotNHists(std::vector<TH1F*> histVector, std::vector<string> histNameVector,
                string xTitle, string yTitle, string histTitle, string legendTitle, string saveName,
                double xMinFrame, double xMaxFrame, double yMinFrame, double yMaxFrame,
                double xMinLegend, double xMaxLegend, double yMinLegend, double yMaxLegend,
                bool setLogY, string setHistDrawOption)
{
  // TODO: how to properly handle exceptions?
  // Check if each histogram/function has an associated name
  if (histVector.size() != histNameVector.size()){
    cout << "plotNHists: hist vector and hist name vector have different sizes (" << histVector.size() << ", " << histNameVector.size() << ")" << endl << "Aborting" << endl;
    return;
  }

  auto nc = new TCanvas("Plot", "Plot", 600, 600);
  nc->SetLeftMargin(0.15);
  nc->cd();

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(TString::Format("%s", histTitle.c_str()).Data());
  frame->Draw();
  string histDrawOption = "same";
  if (setHistDrawOption != "") {
    histDrawOption = TString::Format("%s, %s", histDrawOption.c_str(), setHistDrawOption.c_str()).Data();
  }

  auto legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle);
  legend->SetTextFont(42);
  // legend->SetTextSize(0.04);
  legend->SetTextSize(0.06);
  legend->SetBorderSize(0);

  for (unsigned int i = 0; i < histVector.size(); i++){
    TH1F* hist = histVector[i];
    string name = histNameVector[i];
    hist->SetLineColor(GetColor(i));
    hist->SetMarkerColor(GetColor(i));
    hist->SetMarkerStyle(GetMarker(i));
    hist->SetLineWidth(3);
    legend->AddEntry(hist, TString::Format("%s", name.c_str()).Data());
    hist->Draw(histDrawOption.c_str());
  }
  legend->Draw("same");
  if (setLogY) { nc->SetLogy(); }
  nc->SaveAs(TString::Format("./%s.pdf", saveName.c_str()).Data());
  delete nc;
}
// Plot n histograms and n functions
void plotNHists(std::vector<TH1F*> histVector, std::vector<string> histNameVector,
                std::vector<TF1*> funcVector, std::vector<string> funcNameVector,
                string xTitle, string yTitle, string histTitle, string saveName,
                double xMinFrame, double xMaxFrame, double yMinFrame, double yMaxFrame,
                double xMinLegend, double xMaxLegend, double yMinLegend, double yMaxLegend,
                bool setLogY, string setHistDrawOption, string setFuncDrawOption)
{
  // TODO: how to properly handle exceptions?
  // Check if each histogram/function has an associated name
  if (histVector.size() != histNameVector.size()){
    cout << "plotNHists: hist vector and hist name vector have different sizes (" << histVector.size() << ", " << histNameVector.size() << ")" << endl << "Aborting" << endl;
    return;
  }
  if (funcVector.size() != funcNameVector.size()){
    cout << "plotNHists: func vector and func name vector have different sizes (" << funcVector.size() << ", " << funcNameVector.size() << ")" << endl << "Aborting" << endl;
    return;
  }

  auto nc = new TCanvas("Plot", "Plot", 600, 600);
  nc->SetLeftMargin(0.15);
  nc->cd();

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(TString::Format("%s", histTitle.c_str()).Data());
  frame->Draw();
  string histDrawOption = "same"; string funcDrawOption = "same";
  if (setHistDrawOption != "") {
    histDrawOption = TString::Format("%s, %s", histDrawOption.c_str(), setHistDrawOption.c_str()).Data();
  }
  if (setFuncDrawOption != "") {
    funcDrawOption = TString::Format("%s, %s", funcDrawOption.c_str(), setFuncDrawOption.c_str()).Data();
  }

  auto legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend);
  legend->SetTextFont(42);
  // legend->SetTextSize(0.04);
  legend->SetTextSize(0.06);
  legend->SetBorderSize(0);

  for (unsigned int i = 0; i < histVector.size(); i++){
    TH1F* hist = histVector[i];
    string name = histNameVector[i];
    hist->SetLineColor(GetColor(i));
    hist->SetMarkerColor(GetColor(i));
    hist->SetMarkerStyle(GetMarker(i+4));
    legend->AddEntry(hist, TString::Format("%s", name.c_str()).Data());
    hist->Draw(histDrawOption.c_str());
  }
  for (unsigned int i = 0; i < funcVector.size(); i++){
    TF1* func = funcVector[i];
    string name = funcNameVector[i];
    func->SetLineColor(GetColor(i));
    func->SetMarkerColor(GetColor(i));
    func->SetMarkerStyle(GetMarker(i+4));
    legend->AddEntry(func, TString::Format("%s", name.c_str()).Data());
    func->Draw(funcDrawOption.c_str());
  }
  legend->Draw("same");
  if (setLogY) { nc->SetLogy(); }
  nc->SaveAs(TString::Format("./%s.pdf", saveName.c_str()).Data());
  delete nc;
}

// Project a 2D histogram onto a 1D histogram
TH1F* projectHist(TH2F* inputHist, string projectionAxis, string histName, double axMin, double axMax)
{
  TH1F* outputHist;
  TH2F* inputClone = (TH2F*)inputHist->Clone("CloneOfInput");
  int firstBin = 1; int lastBin = -1;

  if (projectionAxis == "X" || projectionAxis == "x") {
    if (axMin > -900) { firstBin = inputClone->GetYaxis()->FindBin(axMin); }
    if (axMax > -900) { lastBin = inputClone->GetYaxis()->FindBin(axMax); }
    else { lastBin = inputClone->GetNbinsY() + 1; }
    if (firstBin == lastBin) {
      if (firstBin == inputClone->GetNbinsY() + 1 || lastBin == 0) {
        cout << "projectHist: requested bins out of range for Y axis. Aborting" << endl
          << "Requested: (" << axMin << ", " << axMax << ")" << endl
          << "Bins: (" << firstBin << ", " << lastBin << ")" << endl
          << "Max bin: " << inputClone->GetNbinsY() + 1;
        return nullptr;
      }
    }
    outputHist = (TH1F*)inputClone->ProjectionX(TString::Format("%s", histName.c_str()).Data(), firstBin, lastBin);
  }
  else if (projectionAxis == "Y" || projectionAxis == "y") {
    if (axMin > -900) { firstBin = inputClone->GetXaxis()->FindBin(axMin); }
    if (axMax > -900) { lastBin = inputClone->GetXaxis()->FindBin(axMax); }
    else { lastBin = inputClone->GetNbinsX() + 1; }
    if (firstBin == lastBin) {
      if (firstBin == inputClone->GetNbinsX() + 1 || lastBin == 0) {
        cout << "projectHist: requested bins out of range for X axis. Aborting" << endl
          << "Requested: (" << axMin << ", " << axMax << ")" << endl
          << "Bins: (" << firstBin << ", " << lastBin << ")" << endl
          << "Max bin: " << inputClone->GetNbinsX() + 1;
        return nullptr;
      }
    }
    outputHist = (TH1F*)inputClone->ProjectionY(TString::Format("%s", histName.c_str()).Data(), firstBin, lastBin);
  }
  else {
    cout << "make_projection: invalid projection axis " << projectionAxis << ". Aborting." << endl;
    return nullptr;
  }
  delete inputClone;
  return outputHist;
}
// Project a 4D histogram onto a 2D histogram
TH2F* projectHist(THnF* inputHist, int projectionAxisX, int projectionAxisY, string histName,
                  double axis0Min, double axis0Max, double axis1Min, double axis1Max,
                  double axis2Min, double axis2Max, double axis3Min, double axis3Max)
{
  // TH2F* outputHist;
  THnF* inputClone = (THnF*)inputHist->Clone("CloneOfInput");

  int firstBinAxis0 = 1; int lastBinAxis0 = inputClone->GetAxis(0)->GetNbins() + 1;
  int firstBinAxis1 = 1; int lastBinAxis1 = inputClone->GetAxis(1)->GetNbins() + 1;
  int firstBinAxis2 = 1; int lastBinAxis2 = inputClone->GetAxis(2)->GetNbins() + 1;
  int firstBinAxis3 = 1; int lastBinAxis3 = inputClone->GetAxis(3)->GetNbins() + 1;

  if (axis0Min > -900) { firstBinAxis0 = inputClone->GetAxis(0)->FindBin(axis0Min); }
  if (axis0Max > -900) { lastBinAxis0 = inputClone->GetAxis(0)->FindBin(axis0Max); }
  if (axis1Min > -900) { firstBinAxis1 = inputClone->GetAxis(1)->FindBin(axis1Min); }
  if (axis1Max > -900) { lastBinAxis1 = inputClone->GetAxis(1)->FindBin(axis1Max); }
  if (axis2Min > -900) { firstBinAxis2 = inputClone->GetAxis(2)->FindBin(axis2Min); }
  if (axis2Max > -900) { lastBinAxis2 = inputClone->GetAxis(2)->FindBin(axis2Max); }
  if (axis3Min > -900) { firstBinAxis3 = inputClone->GetAxis(3)->FindBin(axis3Min); }
  if (axis3Max > -900) { lastBinAxis3 = inputClone->GetAxis(3)->FindBin(axis3Max); }

  inputClone->GetAxis(0)->SetRange(firstBinAxis0, lastBinAxis0);
  inputClone->GetAxis(1)->SetRange(firstBinAxis1, lastBinAxis1);
  inputClone->GetAxis(2)->SetRange(firstBinAxis2, lastBinAxis2);
  inputClone->GetAxis(3)->SetRange(firstBinAxis3, lastBinAxis3);
  TH2F* outputHist = (TH2F*)inputClone->Projection(projectionAxisY, projectionAxisX);
  outputHist->SetName(TString::Format("%s", histName.c_str()).Data());
  delete inputClone;
  return outputHist;
}

// Get bins corresponding to a range on an axis
// Default behaviour is to include overflow and underflow bins
std::array<int, 2> getProjectionBins(const TAxis* axis, const double min, const double max, double epsilon = 1e-3)
{
  int firstBin = 0, lastBin = axis->GetNbins() + 1;
  if (min > -900) { firstBin = axis->FindBin(min + epsilon); }
  if (max > -900) { lastBin = axis->FindBin(max - epsilon); }
  return std::array{firstBin, lastBin};
}

// Check if a histogram is empty in a given range
bool isHistEmptyInRange(TH1* h, int low, int high, double threshold = 1e-10)
{
  double integral = h->Integral(low, high);
  if (std::isnan(integral))
    return true;
  else
    return (integral < threshold);
}
bool isHistEmptyInRange(TH2* h, int xlow, int xhigh, int ylow, int yhigh, double threshold = 1e-10)
{
  return isHistEmptyInRange(h->ProjectionX("px", ylow, yhigh), xlow, xhigh, threshold);
}

// Returns the upper or lower bound for drawing histogram. Optionally accounts for error
double getHistScale(TH1* h, bool doError, bool doMin = false)
{
  int bin;
  if (!doError) {
    bin = doMin ? h->GetMinimumBin() : h->GetMaximumBin();
    return h->GetBinContent(bin);
  }
  // Initialise with max/min bin content to ensure scale is overwritten when checking min/max
  bin = doMin ? h->GetMaximumBin() : h->GetMinimumBin();
  double scale = h->GetBinContent(bin);
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    if (doMin) scale = min(scale, bc - be);
    else scale = max(scale, bc + be);
  }
  return scale;
}
double getHistScale(vector<TH1*> v, bool doError, bool doSum, bool doMin = false)
{
  double scale = 0;
  if (doSum) {
    TH1* sum = (TH1*)v[0]->Clone("sum");
    sum->Reset();
    for (auto h : v) sum->Add(h);
    scale = getHistScale(sum, doError, doMin);
  }
  else {
    for (auto h : v) {
      double s = getHistScale(h, doError, doMin);
      if (doMin) scale = min(scale, s);
      else scale = max(scale, s);
    }
  }
  return scale;
}
// Returns the lower bound for drawing histogram. Optionally accounts for error
double getHistLowerBound(TH1* h, bool doError)
{
  return getHistScale(h, doError, true);
}
double getHistLowerBound(vector<TH1*> v, bool doError, bool doSum)
{
  return getHistScale(v, doError, doSum, true);
}
// Returns the upper bound for drawing histogram. Optionally accounts for error
double getHistUpperBound(TH1* h, bool doError)
{
  return getHistScale(h, doError, false);
}
double getHistUpperBound(vector<TH1*> v, bool doError, bool doSum)
{
  return getHistScale(v, doError, doSum, false);
}

// Set histogram colours and markers
void setStyle(TH1* hist, int styleNumber, double alpha = -1, int lineWidth = 3)
{
  hist->SetLineWidth(lineWidth);
  hist->SetLineColor(GetColor(styleNumber));
  hist->SetMarkerStyle(GetMarker(styleNumber));
  hist->SetMarkerColor(GetColor(styleNumber));
  if (alpha > 0.) {
    hist->SetFillColorAlpha(GetColor(styleNumber), alpha);
  }
}
void setStyle(TLine* line, int styleNumber, int lineStyle = 9, int lineWidth = 3)
{
  line->SetLineWidth(lineWidth);
  line->SetLineColor(GetColor(styleNumber));
  line->SetLineStyle(lineStyle);
}
void setStyle(TF1* f, int styleNumber, int lineStyle = 1, int lineWidth = 3)
{
  f->SetLineWidth(lineWidth);
  f->SetLineColor(GetColor(styleNumber));
  f->SetLineStyle(lineStyle);
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
  else if (hadron == "K0S"){
    had = "K^{0}_{S}";
  }
  else if (hadron == "Lambda0" || hadron == "Lambda"){
    // had = "#Lambda^{0}";
    had = "#Lambda";
  }
  else if (hadron == "AntiLambda0" || hadron == "AntiLambda"){
    // had = "#bar{#Lambda}^{0}";
    had = "#bar{#Lambda}";
  }
  return had;
}
// Returns decay products given a hadron
string formatHadronDaughters(string hadron)
{
  string daughters = hadron;
  if ("K0S" == hadron) {
    daughters = "#pi^{+}#pi^{-}";
  }
  else if ("Lambda0" == hadron || "Lambda" == hadron) {
    daughters = "p#pi^{-}";
  }
  else if ("AntiLambda0" == hadron || "AntiLambda" == hadron) {
    daughters = "#bar{p}#pi^{+}";
  }
  return daughters;
}
// Normalise 2D histogram row-by-row
void normaliseHistRowByRow(TH2* hist)
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
void normaliseHistColByCol(TH2* hist)
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
#endif // HISTUTILS_H
