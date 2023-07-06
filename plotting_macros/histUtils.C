//
// This macro plots the rejection factor as a function of the trigger patch energy cut
//

#include <vector>
#include <iostream>
#include <typeinfo>
#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TString.h"

#include "TLegend.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TROOT.h"

#include "plotUtils.C"

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
// /* Project a 2D histogram onto a 1D histogram */
// TH1F* projectHist(TH2F* inputHist, string projectionAxis, string histName,
//                   double xMin, double xMax, double yMin, double yMax)
// {
//   TH1F* outputHist;
//   TH2F* inputClone = (TH2F*)inputHist->Clone("CloneOfInput");

//   int firstXBin = 1; int lastXBin = inputHist->GetNbinsX() + 1;
//   int firstYBin = 1; int lastYBin = inputHist->GetNbinsY() + 1;

//   if (xMin > -900) { firstXBin = inputHist->GetXaxis()->FindBin(xMin); }
//   if (xMax > -900) { lastXBin = inputHist->GetXaxis()->FindBin(xMax); }
//   if (yMin > -900) { firstYBin = inputHist->GetYaxis()->FindBin(yMin); }
//   if (yMax > -900) { lastYBin = inputHist->GetYaxis()->FindBin(yMax); }

//   if (projectionAxis == "X" || projectionAxis == "x") {
//     inputClone->GetXaxis()->SetRange(firstXBin, lastXBin);
//     outputHist = (TH1F*)inputHist->ProjectionX(TString::Format("%s", histName.c_str()).Data(), firstYBin, lastYBin);
//   }
//   else if (projectionAxis == "Y" || projectionAxis == "y") {
//     inputClone->GetYaxis()->SetRange(firstYBin, lastYBin);
//     outputHist = (TH1F*)inputHist->ProjectionY(TString::Format("%s", histName.c_str()).Data(), firstXBin, lastXBin);
//   }
//   else {
//     cout << "make_projection: invalid projection axis " << projectionAxis << ". Aborting." << endl;
//   }
//   delete inputClone;
//   return outputHist;
// }
// Project a 3D histogram onto a 1D histogram
// TH1F* projectHist(TH3F* inputHist, string projectionAxis, string histName,
//                   double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
// {
//   TH1F* outputHist;
//   TH3F* inputClone = (TH3F*)inputHist->Clone("CloneOfInput");

//   int firstXBin = 1; int lastXBin = inputClone->GetNbinsX() + 1;
//   int firstYBin = 1; int lastYBin = inputClone->GetNbinsY() + 1;
//   int firstZBin = 1; int lastZBin = inputClone->GetNbinsZ() + 1;

//   if (xMin > -900) { firstXBin = inputClone->GetXaxis()->FindBin(xMin); }
//   if (xMax > -900) { lastXBin = inputClone->GetXaxis()->FindBin(xMax); }
//   if (yMin > -900) { firstYBin = inputClone->GetYaxis()->FindBin(yMin); }
//   if (yMax > -900) { lastYBin = inputClone->GetYaxis()->FindBin(yMax); }
//   if (zMin > -900) { firstZBin = inputClone->GetZaxis()->FindBin(zMin); }
//   if (zMax > -900) { lastZBin = inputClone->GetZaxis()->FindBin(zMax); }

//   if (projectionAxis == "X" || projectionAxis == "x") {
//     outputHist = (TH1F*)inputClone->ProjectionX(TString::Format("%s", histName.c_str()).Data(),
//                                                 firstYBin, lastYBin, firstZBin, lastZBin);
//   }
//   else if (projectionAxis == "Y" || projectionAxis == "y") {
//     outputHist = (TH1F*)inputClone->ProjectionY(TString::Format("%s", histName.c_str()).Data(),
//                                                 firstXBin, lastXBin, firstZBin, lastZBin);
//   }
//   else if (projectionAxis == "Z" || projectionAxis == "z") {
//     outputHist = (TH1F*)inputClone->ProjectionZ(TString::Format("%s", histName.c_str()).Data(),
//                                                 firstXBin, lastXBin, firstYBin, lastYBin);
//   }
//   else {
//     cout << "make_projection: invalid projection axis " << projectionAxis << ". Aborting." << endl;
//   }
//   // delete inputClone;
//   return outputHist;
// }
// Project a 3D histogram onto a 1D or 2D histogram
// TH2F* projectHist(TH3F* inputHist, string projectionAxis,
//                   double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
// {
//   TH2F* outputHist;
//   TH3F* inputClone = (TH3F*)inputHist->Clone("CloneOfInput");

//   int firstXBin = 1; int lastXBin = inputClone->GetNbinsX() + 1;
//   int firstYBin = 1; int lastYBin = inputClone->GetNbinsY() + 1;
//   int firstZBin = 1; int lastZBin = inputClone->GetNbinsZ() + 1;

//   if (xMin > -900) { firstXBin = inputClone->GetXaxis()->FindBin(xMin); }
//   if (xMax > -900) { lastXBin = inputClone->GetXaxis()->FindBin(xMax); }
//   if (yMin > -900) { firstYBin = inputClone->GetYaxis()->FindBin(yMin); }
//   if (yMax > -900) { lastYBin = inputClone->GetYaxis()->FindBin(yMax); }
//   if (zMin > -900) { firstZBin = inputClone->GetZaxis()->FindBin(zMin); }
//   if (zMax > -900) { lastZBin = inputClone->GetZaxis()->FindBin(zMax); }

//   inputClone->GetXaxis()->SetRange(firstXBin, lastXBin);
//   inputClone->GetYaxis()->SetRange(firstYBin, lastYBin);
//   inputClone->GetZaxis()->SetRange(firstZBin, lastZBin);
//   outputHist = (TH2F*)inputClone->Project3D(projectionAxis.c_str());
//   delete inputClone;
//   return outputHist;
// }
// Project a 4D histogram onto a 2D histogram
TH2F* projectHist(THnF* inputHist, int projectionAxisX, int projectionAxisY,
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
  delete inputClone;
  return outputHist;
}
