
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

#include "../histUtils.C"

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;
gStyle->SetNdivisions(505, "xy");

int getLastFilledRowOrCol(TH2* h, bool doRow, double threshold = 1e-5)
{
  if (doRow) {
    for (int i = 0; i < h->GetNbinsY(); i++) {
      int iRow = h->GetNbinsY() - i;
      double integral = h->Integral(1, h->GetNbinsX(), iRow, iRow);
      if (integral != integral) continue; // NaN check
      if (integral > threshold) return iRow;
    }
  } else {
    for (int i = 0; i < h->GetNbinsX(); i++) {
      int iCol = h->GetNbinsX() - i;
      double integral = h->Integral(iCol, iCol, 1, h->GetNbinsY());
      if (integral != integral) continue; // NaN check
      if (integral > threshold) return iCol;
    }
  }
  return -1;
}
int getLastFilledRow(TH2* h, double threshold = 1e-5) { return getLastFilledRowOrCol(h, true, threshold); }
int getLastFilledCol(TH2* h, double threshold = 1e-5) { return getLastFilledRowOrCol(h, false, threshold); }

// Returns a subset of a histogram from minBin to maxBin
template <typename T>
T* makeHistSubset(T* data, int minBin, int maxBin)
{
  T* region = (T*)data->Clone("region");
  region->Reset();
  for (int i = minBin; i <= maxBin; i++) {
    region->SetBinContent(i, data->GetBinContent(i));
  }
  return region;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

void ptAsymmetry(vector<string> inputStrings, bool normalise) // 2D version
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];

  const int nDim       = 4;
  const int v0PtAxis   = 0;
  const int posPtAxis  = 1;
  const int negPtAxis  = 2;
  const int ptDiffAxis = 3;

  string histName = "jet-v0qa/tracks/Pt";
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  TH2D* ptDiff = (TH2D*)thn->Projection(ptDiffAxis, v0PtAxis);

  if (isHistEmptyInRange(ptDiff, 1, ptDiff->GetNbinsX(), 1, ptDiff->GetNbinsY())) {
    cout << "Pt asymmetry: empty hist! Skipping" << endl;
    return;
  }

  string saveName = "daughterPtDiff2D";
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  ptDiff->SetName(saveName.c_str());
  int lastCol = getLastFilledRowOrCol(ptDiff, false);
  double v0max = ptDiff->GetXaxis()->GetBinUpEdge(lastCol);
  if (normalise) {
    normaliseHistColByCol(ptDiff);
    ptDiff->SetMaximum(1.);
    ptDiff->SetMinimum(1e-5);
  }

  string xTitle = "#it{p}_{T, V0} (GeV/#it{c})";
  string yTitle = "#it{p}_{T, pos} - #it{p}_{T, neg} (GeV/#it{c})";
  double xMinFrame = 0., xMaxFrame = v0max;
  double yMinFrame = -1. * v0max, yMaxFrame = v0max;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  canvas->cd();
  frame->Draw();
  ptDiff->Draw("same colz");
  canvas->SaveAs(saveName.c_str());
}

void ptAsymmetry(vector<string> inputStrings, double v0min, double v0max, bool normalise)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];

  const int nDim       = 4;
  const int v0PtAxis   = 0;
  const int posPtAxis  = 1;
  const int negPtAxis  = 2;
  const int ptDiffAxis = 3;

  string histName = "jet-v0qa/tracks/Pt";
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(0), v0min, v0max);
  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);
  TH1D* ptDiff = (TH1D*)thn->Projection(ptDiffAxis);

  if (isHistEmptyInRange(ptDiff, 1, ptDiff->GetNbinsX())) {
    cout << "Pt asymmetry: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << "! Skipping" << endl;
    return;
  }

  string saveName = "daughterPtDiff";
  saveName += TString::Format("_v0pt%.1f-%.1f", lowv0, highv0).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  ptDiff->SetName(saveName.c_str());
  setStyle(ptDiff, 0);
  if (normalise) ptDiff->Scale(1./ptDiff->Integral(), "width");

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string xTitle = "#it{p}_{T, pos} - #it{p}_{T, neg} (GeV/#it{c})";
  string yTitle = normalise ? "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d(#it{p}_{T, pos} - #it{p}_{T, neg})}" : "Counts";
  // double xMinFrame = ptDiff->GetXaxis()->GetXmin(), xMaxFrame = ptDiff->GetXaxis()->GetXmax();
  double xMinFrame = -1. * highv0, xMaxFrame = highv0;
  double yMinFrame = 0., yMaxFrame = 1.5 * getHistScale(ptDiff, true);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  canvas->cd();
  frame->Draw();
  ptDiff->Draw("same");
  canvas->SaveAs(saveName.c_str());
}

// Mass in bins of pt diff
void mass(vector<string> inputStrings, double v0min, double v0max, bool normalise, bool doStack)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  const int nDim            = 5;
  const int v0PtAxis        = 0;
  const int ptDiffAxis      = 1;
  const int K0SAxis         = 2;
  const int Lambda0Axis     = 3;
  const int AntiLambda0Axis = 4;

  string histName = "jet-v0qa/tracks/PtDiffMass";
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(0), v0min, v0max);
  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  TH1D* ptDiffHist = (TH1D*)thn->Projection(ptDiffAxis);
  if (isHistEmptyInRange(ptDiffHist, 1, ptDiffHist->GetNbinsX())) {
    cout << "Pt asymmetry: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << "! Skipping" << endl;
    return;
  }

  string saveName = "daughterPtDiffm";
  saveName += hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowv0, highv0).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  vector<double> ptBins;
  int nBins = 5;
  for (int i = 0; i < nBins; i++) {
    double binWidth = highv0 / nBins;
    ptBins.push_back(highv0 - binWidth * i);
  }
  vector<TH1*> ptDiffs;

  TLegend* legend = CreateLegend(0.25, 0.85, 0.61, 0.9, "", 0.04);
  THnSparseD* thn_copy = (THnSparseD*)thn->Clone();
  int projectionAxis = (hadron == "K0S")*K0SAxis + (hadron == "Lambda0")*Lambda0Axis + (hadron == "AntiLambda0")*AntiLambda0Axis;

  for (int i = 0; i < ptBins.size(); i++) {
    array<int, 2> ptDiffBins = getProjectionBins(thn_copy->GetAxis(ptDiffAxis), -1. * ptBins[i], ptBins[i]);
    thn_copy->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
    thn_copy->GetAxis(ptDiffAxis)->SetRange(ptDiffBins[0], ptDiffBins[1]);

    TH1D* mass = (TH1D*)thn_copy->Projection(projectionAxis);
    mass->SetName((saveName + "_m" + hadron + TString::Format("_dpt%.1f", ptBins[i]).Data()).c_str());
    if (!isHistEmptyInRange(mass, 1, mass->GetNbinsX(), 1.)) {
      ptDiffs.push_back(mass);
      legend->AddEntry(mass, TString::Format("#delta#it{p}_{T} = %.1f - %.1f GeV", ptBins[i+1], ptBins[i]).Data(), "lf");
    } else {
      cout << "Mass: empty histogram for dPt < " << ptBins[i] << " GeV! Skipping" << endl;
    }
  }

  if (ptDiffs.size() == 0) {
    cout << "Pt asymmetry: no histograms filled! Skipping" << endl;
    return;
  }

  // Subtract out the inner range of dPt
  for (int i = 0; i < ptDiffs.size(); i++) {
    if (i < ptDiffs.size() - 1) ptDiffs[i]->Add(ptDiffs[i + 1], -1.);

    double alpha = doStack ? 0.5 : -1.;
    setStyle(ptDiffs[i], i+1, alpha);
    ptDiffs[i]->Rebin(4);
  }

  if (normalise) {
    double scale = getHistScale(ptDiffs, false, false);
    for (auto hist : ptDiffs) {
      hist->Scale(1./scale, "width");
    }
  }

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string xTitle = "M(" + formatHadronDaughters(hadron) + ")";
  string yTitle = "Counts";
  if (normalise) yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
  double xMinFrame = ptDiffs[0]->GetXaxis()->GetXmin(), xMaxFrame = ptDiffs[0]->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 2. * getHistScale(ptDiffs, normalise, doStack);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  double mass = (hadron == "K0S") ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, 0.55*yMaxFrame);
  setStyle(line, 0);

  canvas->cd();
  frame->Draw();
  legend->Draw("same");
  if (doStack) {
    THStack* stack = new THStack("stack", "");
    for (auto h : ptDiffs) stack->Add(h);
    stack->Draw("same");
    line->Draw("same");
  } else {
    line->Draw("same");
    for (auto p : ptDiffs) p->Draw("same");
  }
  canvas->SaveAs(saveName.c_str());
}

void mass2D(vector<string> inputStrings, double v0min, double v0max, bool normalise) // 2D version
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  const int nDim            = 5;
  const int v0PtAxis        = 0;
  const int ptDiffAxis      = 1;
  const int K0SAxis         = 2;
  const int Lambda0Axis     = 3;
  const int AntiLambda0Axis = 4;

  string histName = "jet-v0qa/tracks/PtDiffMass";
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(0), v0min, v0max);
  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  int projectionAxis = (hadron == "K0S")*K0SAxis + (hadron == "Lambda0")*Lambda0Axis + (hadron == "AntiLambda0")*AntiLambda0Axis;
  TH2D* ptDiff = (TH2D*)thn->Projection(ptDiffAxis, projectionAxis);
  if (isHistEmptyInRange(ptDiff, 1, ptDiff->GetNbinsX(), 1, ptDiff->GetNbinsY())) {
    cout << "Pt asymmetry: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << "! Skipping" << endl;
    return;
  }
  ptDiff->Rebin2D(4, 1);

  string saveName = "daughterPtDiffm";
  saveName += hadron;
  saveName += TString::Format("2D_v0pt%.1f-%.1f", lowv0, highv0).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);
  canvas->SetLogz();

  ptDiff->SetName(saveName.c_str());
  int lastRow = getLastFilledRowOrCol(ptDiff, true);
  double diffMax = ptDiff->GetXaxis()->GetBinUpEdge(lastRow);
  if (normalise) {
    normaliseHistRowByRow(ptDiff);
    ptDiff->SetMaximum(1.);
    ptDiff->SetMinimum(1e-5);
  }

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string xTitle = "M(" + formatHadronDaughters(hadron) + ")";
  string yTitle = "#it{p}_{T, pos} - #it{p}_{T, neg} (GeV/#it{c})";
  // double xMinFrame = 0.4, xMaxFrame = 0.6;
  double yMinFrame = -1. * highv0, yMaxFrame = highv0;
  double xMinFrame = ptDiff->GetXaxis()->GetXmin(), xMaxFrame = ptDiff->GetXaxis()->GetXmax();
  // double yMinFrame = -1. * diffMax, yMaxFrame = diffMax;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  double mass = (hadron == "K0S") ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, yMaxFrame);
  setStyle(line, 0);

  canvas->cd();
  frame->Draw();
  ptDiff->Draw("same colz");
  line->Draw("same");
  canvas->SaveAs(saveName.c_str());
}

void trd(vector<string> inputStrings, double v0min, double v0max, bool normalise, bool doStack)
{
  const int nDim            = 6;
  const int v0PtAxis        = 0;
  const int posPtAxis       = 1;
  const int negPtAxis       = 2;
  const int K0SAxis         = 3;
  const int Lambda0Axis     = 4;
  const int AntiLambda0Axis = 5;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  int projectionAxis = (hadron == "K0S")*K0SAxis + (hadron == "Lambda0")*Lambda0Axis + (hadron == "AntiLambda0")*AntiLambda0Axis;

  string histDir = "jet-v0qa/tracks/";
  vector<string> histNames = {"posNoTRD", "posTRD", "negNoTRD", "negTRD"};
  TFile* inFile = TFile::Open(inName.c_str());
  vector<TH1*> hists;

  string saveName = "daughterTRD_m";
  saveName += hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", v0min, v0max).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TLegend* legend = CreateLegend(0.25, 0.85, 0.65, 0.9, "", 0.04);

  for (int i = 0; i < histNames.size(); i++) {
    THnSparseD* thn = (THnSparseD*)inFile->Get((histDir + histNames[i] + "PtMass").c_str());
    array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);

    thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
    double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
    double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

    TH1D* mass = (TH1D*)thn->Projection(projectionAxis);
    if (isHistEmptyInRange(mass, 1, mass->GetNbinsX())) {
      cout << "TRD: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << "! Skipping" << endl;
      continue;
    }
    mass->SetName((saveName + "_" + histNames[i]).c_str());
    legend->AddEntry(mass, histNames[i].c_str(), "lf");
    hists.push_back(mass);
  }
  if (hists.size() == 0) {
    cout << "TRD: no histograms filled! Skipping" << endl;
    return;
  }

  for (int i = 0; i < hists.size(); i++) {
    double alpha = doStack ? 0.5 : -1.;
    setStyle(hists[i], i+1, alpha);
    hists[i]->Rebin(4);
  }
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", v0min, v0max).Data();
  string xTitle = "M(" + formatHadronDaughters(hadron) + ")";
  string yTitle = "Counts";
  if (normalise) yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
  double xMinFrame = hists[0]->GetXaxis()->GetXmin(), xMaxFrame = hists[0]->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.5 * getHistScale(hists, normalise, doStack);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  double mass = (hadron == "K0S") ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, 0.55*yMaxFrame);
  setStyle(line, 0);

  // for (auto h : hists) {
  //   cout << getHistScale(h, false) << endl;
  // }

  canvas->cd();
  frame->Draw();
  legend->Draw("same");
  if (doStack) {
    THStack* stack = new THStack("stack", "");
    for (auto h : hists) stack->Add(h);
    stack->Draw("same");
    line->Draw("same");
  } else {
    line->Draw("same");
    for (auto p : hists) p->Draw("same");
  }
  canvas->SaveAs(saveName.c_str());
}

TH1D* itslayer(vector<string> inputStrings, double v0min, double v0max, bool doPos, int layer)
{
  const int nDim      = 4;
  const int v0PtAxis  = 0;
  const int posPtAxis = 1;
  const int negPtAxis = 2;
  const int massAxis  = 3;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  string histDir = "jet-v0qa/tracks/ITS/";
  TFile* inFile = TFile::Open(inName.c_str());
  vector<TH1*> hists;

  string posNeg = doPos ? "pos" : "neg";
  string histName = histDir + posNeg + "Layer" + to_string(layer) + "Mass" + hadron;
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);
  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  TH1D* mass = (TH1D*)thn->Projection(massAxis);
  if (isHistEmptyInRange(mass, 1, mass->GetNbinsX())) {
    cout << "ITS: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << ", " << posNeg << " layer " << layer << "! Skipping" << endl;
    return nullptr;
  }

  string saveName = TString::Format("_v0pt%.1f-%.1f", lowv0, highv0).Data();
  mass->SetName(saveName.c_str());
  return mass;
}
void itslayer(vector<string> inputStrings, double v0min, double v0max, bool normalise, bool doStack, int layer)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
  vector<TH1*> hists;
  double alpha = doStack ? 0.5 : -1.;

  TLegend* legend = CreateLegend(0.25, 0.85, 0.62, 0.9, ("ITS layer " + to_string(layer)).c_str(), 0.04);
  bool fullHists[2] = {false, false};

  TH1D* posmass = itslayer(inputStrings, v0min, v0max, true, layer);
  if (posmass) {
    fullHists[0] = true;
    posmass->SetName(("posITS" + to_string(layer) + hadron + posmass->GetName()).c_str());
    setStyle(posmass, 1, alpha);
    legend->AddEntry(posmass, "pos", "lf");
    hists.push_back(posmass);
  }
  TH1D* negmass = itslayer(inputStrings, v0min, v0max, false, layer);
  if (negmass) {
    fullHists[1] = true;
    negmass->SetName(("negITS" + to_string(layer) + hadron + negmass->GetName()).c_str());
    setStyle(negmass, 2, alpha);
    legend->AddEntry(negmass, "neg", "lf");
    hists.push_back(negmass);
  }

  if (!fullHists[0] && !fullHists[1]) {
    cout << "ITS: no histograms filled for layer " << layer << "! Skipping" << endl;
    return;
  }

  for (auto h : hists) h->Rebin(4);

  string saveName = "ITS";
  saveName += to_string(layer);
  saveName += "_m";
  saveName += hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", v0min, v0max).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", v0min, v0max).Data();
  string xTitle = "M(" + formatHadronDaughters(hadron) + ")";
  string yTitle = "Counts";
  if (normalise) yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
  double xMinFrame = hists[0]->GetXaxis()->GetXmin(), xMaxFrame = hists[0]->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.5 * getHistScale(hists, normalise, doStack);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  double mass = (hadron == "K0S") ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, 0.6*yMaxFrame);
  setStyle(line, 0);

  canvas->cd();
  frame->Draw();
  legend->Draw("same");
  if (doStack) {
    THStack* stack = new THStack("stack", "");
    for (auto h : hists) stack->Add(h);
    stack->Draw("same");
    line->Draw("same");
  } else {
    line->Draw("same");
    for (auto h : hists) h->Draw("same");
  }
  canvas->SaveAs(saveName.c_str());
}

// FIXME: histograms need to be changed. Currently does not include mass
void itsNCl(vector<string> inputStrings, double v0min, double v0max, bool normalise, bool doStack)
{
  const int nDim            = 5;
  const int v0PtAxis        = 0;
  const int posPtAxis       = 1;
  const int negPtAxis       = 2;
  const int nClAxis         = 3;
  const int massAxis        = 4;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];
}
void itsChi2NCl(){}
void tpc(){}

void dcaXY(vector<string> inputStrings, double v0min, double v0max)
{
  const int nDim       = 6;
  const int v0PtAxis   = 0;
  const int posPtAxis  = 1;
  const int negPtAxis  = 2;
  const int posDCAAxis = 3;
  const int negDCAAxis = 4;
  const int massAxis   = 5;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  string histName = "jet-v0qa/tracks/DCAxyMass";
  histName += hadron;
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);
  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  TH2D* dcaXY = (TH2D*)thn->Projection(posDCAAxis, negDCAAxis);
  dcaXY->Draw("colz");
}
void dcaXYQuadrants(vector<string> inputStrings, double v0min, double v0max, bool normalise, bool doStack)
{
  const int nDim       = 6;
  const int v0PtAxis   = 0;
  const int posPtAxis  = 1;
  const int negPtAxis  = 2;
  const int posDCAAxis = 3;
  const int negDCAAxis = 4;
  const int massAxis   = 5;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  string histName = "jet-v0qa/tracks/DCAxyMass";
  histName += hadron;
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  array<int, 2> posdcabins; array<int, 2> negdcabins;
  double highposdca, highnegdca;
  string posString, negString;
  vector<TH1*> hists;
  TLegend* legend = CreateLegend(0.25, 0.85, 0.65, 0.9, "", 0.04);

  // Four quadrants for DCA: both > 0, both < 0, 2 * (>0, <0)
  // Pos/Pos
  posdcabins = getProjectionBins(thn->GetAxis(posDCAAxis), 0., 0.5);
  negdcabins = getProjectionBins(thn->GetAxis(negDCAAxis), 0., 0.5);
  highposdca = thn->GetAxis(posDCAAxis)->GetBinUpEdge(posdcabins[1]);
  highnegdca = thn->GetAxis(negDCAAxis)->GetBinUpEdge(negdcabins[1]);
  posString = (highposdca > 1e-3) ? "pos DCA > 0" : "pos DCA < 0";
  negString = (highnegdca > 1e-3) ? "neg DCA > 0" : "neg DCA < 0";

  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  thn->GetAxis(posDCAAxis)->SetRange(posdcabins[0], posdcabins[1]);
  thn->GetAxis(negDCAAxis)->SetRange(negdcabins[0], negdcabins[1]);
  TH1D* massPosPos = (TH1D*)thn->Projection(massAxis);
  if (isHistEmptyInRange(massPosPos, 1, massPosPos->GetNbinsX())) {
    cout << "DCA xy: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << ", " << posString << ", " << negString << "! Skipping" << endl;
  } else {
    string name = TString::Format("%s_DCAxyPosPos_v0pt%.1f-%.1f", hadron.c_str(), lowv0, highv0).Data();
    massPosPos->SetName(name.c_str());
    hists.push_back(massPosPos);
    legend->AddEntry(massPosPos, (posString + ", " + negString).c_str(), "lf");
  }

  // Pos/Neg
  posdcabins = getProjectionBins(thn->GetAxis(posDCAAxis), 0., 0.5);
  negdcabins = getProjectionBins(thn->GetAxis(negDCAAxis), -0.5, 0.);
  highposdca = thn->GetAxis(posDCAAxis)->GetBinUpEdge(posdcabins[1]);
  highnegdca = thn->GetAxis(negDCAAxis)->GetBinUpEdge(negdcabins[1]);
  posString = (highposdca > 1e-3) ? "pos DCA > 0" : "pos DCA < 0";
  negString = (highnegdca > 1e-3) ? "neg DCA > 0" : "neg DCA < 0";

  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  thn->GetAxis(posDCAAxis)->SetRange(posdcabins[0], posdcabins[1]);
  thn->GetAxis(negDCAAxis)->SetRange(negdcabins[0], negdcabins[1]);
  TH1D* massPosNeg = (TH1D*)thn->Projection(massAxis);
  if (isHistEmptyInRange(massPosNeg, 1, massPosNeg->GetNbinsX())) {
    cout << "DCA xy: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << ", " << posString << ", " << negString << "! Skipping" << endl;
  } else {
    string name = TString::Format("%s_DCAxyPosNeg_v0pt%.1f-%.1f", hadron.c_str(), lowv0, highv0).Data();
    massPosNeg->SetName(name.c_str());
    hists.push_back(massPosNeg);
    legend->AddEntry(massPosNeg, (posString + ", " + negString).c_str(), "lf");
  }

  // Neg/Pos
  posdcabins = getProjectionBins(thn->GetAxis(posDCAAxis), -0.5, 0.);
  negdcabins = getProjectionBins(thn->GetAxis(negDCAAxis), 0., 0.5);
  highposdca = thn->GetAxis(posDCAAxis)->GetBinUpEdge(posdcabins[1]);
  highnegdca = thn->GetAxis(negDCAAxis)->GetBinUpEdge(negdcabins[1]);
  posString = (highposdca > 1e-3) ? "pos DCA > 0" : "pos DCA < 0";
  negString = (highnegdca > 1e-3) ? "neg DCA > 0" : "neg DCA < 0";

  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  thn->GetAxis(posDCAAxis)->SetRange(posdcabins[0], posdcabins[1]);
  thn->GetAxis(negDCAAxis)->SetRange(negdcabins[0], negdcabins[1]);
  TH1D* massNegPos = (TH1D*)thn->Projection(massAxis);
  if (isHistEmptyInRange(massNegPos, 1, massNegPos->GetNbinsX())) {
    cout << "DCA xy: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << ", " << posString << ", " << negString << "! Skipping" << endl;
  } else {
    string name = TString::Format("%s_DCAxyNegPos_v0pt%.1f-%.1f", hadron.c_str(), lowv0, highv0).Data();
    massNegPos->SetName(name.c_str());
    hists.push_back(massNegPos);
    legend->AddEntry(massNegPos, (posString + ", " + negString).c_str(), "lf");
  }

  // Neg/Neg
  posdcabins = getProjectionBins(thn->GetAxis(posDCAAxis), -0.5, 0.);
  negdcabins = getProjectionBins(thn->GetAxis(negDCAAxis), -0.5, 0.);
  highposdca = thn->GetAxis(posDCAAxis)->GetBinUpEdge(posdcabins[1]);
  highnegdca = thn->GetAxis(negDCAAxis)->GetBinUpEdge(negdcabins[1]);
  posString = (highposdca > 1e-3) ? "pos DCA > 0" : "pos DCA < 0";
  negString = (highnegdca > 1e-3) ? "neg DCA > 0" : "neg DCA < 0";

  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  thn->GetAxis(posDCAAxis)->SetRange(posdcabins[0], posdcabins[1]);
  thn->GetAxis(negDCAAxis)->SetRange(negdcabins[0], negdcabins[1]);
  TH1D* massNegNeg = (TH1D*)thn->Projection(massAxis);
  if (isHistEmptyInRange(massNegNeg, 1, massNegNeg->GetNbinsX())) {
    cout << "DCA xy: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << ", " << posString << ", " << negString << "! Skipping" << endl;
  } else {
    string name = TString::Format("%s_DCAxyNegNeg_v0pt%.1f-%.1f", hadron.c_str(), lowv0, highv0).Data();
    massNegNeg->SetName(name.c_str());
    hists.push_back(massNegNeg);
    legend->AddEntry(massNegNeg, (posString + ", " + negString).c_str(), "lf");
  }

  for (int i = 0; i < hists.size(); i++) {
    double alpha = doStack ? 0.5 : -1.;
    setStyle(hists[i], i+1, alpha);
    hists[i]->Rebin(4);
  }

  if (hists.size() == 0) {
    cout << "DCA xy: no histograms filled! Skipping" << endl;
    return;
  }

  string saveName = "daughterDCAxy_m";
  saveName += hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowv0, highv0).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string xTitle = "M(" + formatHadronDaughters(hadron) + ")";
  string yTitle = "Counts";
  if (normalise) {
    yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
    double scale = getHistScale(hists, normalise, doStack);
    for (auto h : hists) h->Scale(1./scale, "width");
  }
  double xMinFrame = hists[0]->GetXaxis()->GetXmin(), xMaxFrame = hists[0]->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.5 * getHistScale(hists, normalise, doStack);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  double mass = (hadron == "K0S") ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, 0.55*yMaxFrame);
  setStyle(line, 0);

  canvas->cd();
  frame->Draw();
  legend->Draw("same");
  if (doStack) {
    THStack* stack = new THStack("stack", "");
    for (auto h : hists) stack->Add(h);
    stack->Draw("same");
    line->Draw("same");
  } else {
    line->Draw("same");
    for (auto p : hists) p->Draw("same");
  }
  canvas->SaveAs(saveName.c_str());
}
void dcaXY(vector<string> inputStrings, double v0min, double v0max, bool normalise, bool doStack, bool doPos)
{
  const int nDim       = 6;
  const int v0PtAxis   = 0;
  const int posPtAxis  = 1;
  const int negPtAxis  = 2;
  const int posDCAAxis = 3;
  const int negDCAAxis = 4;
  const int massAxis   = 5;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  string histName = "jet-v0qa/tracks/DCAxyMass";
  histName += hadron;
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  string posOrNeg = doPos ? "pos" : "neg";
  string saveName = posOrNeg;
  saveName += "DCAxy_m";
  saveName += hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowv0, highv0).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";

  vector<TH1*> hists;
  TLegend* legend = CreateLegend(0.25, 0.85, 0.63, 0.89, "", 0.04);
  // vector<double> bins = {-0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5};
  vector<double> bins = {-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

  int dcaAxis = doPos ? posDCAAxis : negDCAAxis;
  for (int i = 0; i < bins.size() - 1; i++) {
    array<int, 2> dcabins = getProjectionBins(thn->GetAxis(dcaAxis), bins[i], bins[i+1]);
    double lowdca  = thn->GetAxis(dcaAxis)->GetBinLowEdge(dcabins[0]);
    double highdca = thn->GetAxis(dcaAxis)->GetBinUpEdge(dcabins[1]);

    thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
    thn->GetAxis(dcaAxis)->SetRange(dcabins[0], dcabins[1]);
    TH1D* mass = (TH1D*)thn->Projection(massAxis);

    string name = TString::Format("%s_%sDCA%.2f-%.2f", saveName.c_str(), posOrNeg.c_str(), abs(lowdca), abs(highdca)).Data();
    mass->SetName(name.c_str());
    if (isHistEmptyInRange(mass, 1, mass->GetNbinsX())) {
      cout << "DCA xy: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << ", " << posOrNeg << "DCA " << lowdca << " - " << highdca << "! Skipping" << endl;
    } else {
      string legendEntry = TString::Format("%s DCA = %.2f - %.2f cm", posOrNeg.c_str(), lowdca, highdca).Data();
      if (i == 0) legendEntry = TString::Format("%s DCA < %.2f cm", posOrNeg.c_str(), highdca).Data();
      if (i == bins.size() - 2) legendEntry = TString::Format("%s DCA > %.2f cm", posOrNeg.c_str(), lowdca).Data();

      hists.push_back(mass);
      legend->AddEntry(mass, legendEntry.c_str(), "lf");
    }
  }

  for (int i = 0; i < hists.size(); i++) {
    double alpha = doStack ? 0.5 : -1.;
    setStyle(hists[i], i+1, alpha);
    hists[i]->Rebin(4);
  }

  if (hists.size() == 0) {
    cout << "DCA xy: no histograms filled! Skipping" << endl;
    return;
  }

  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string xTitle = "M(" + formatHadronDaughters(hadron) + ")";
  string yTitle = "Counts";
  if (normalise) {
    yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
    double scale = getHistScale(hists, normalise, doStack);
    for (auto h : hists) h->Scale(1./scale, "width");
  }
  double xMinFrame = hists[0]->GetXaxis()->GetXmin(), xMaxFrame = hists[0]->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.5 * getHistScale(hists, normalise, doStack);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  double mass = (hadron == "K0S") ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, 0.55*yMaxFrame);
  setStyle(line, 0);

  canvas->cd();
  frame->Draw();
  legend->Draw("same");
  if (doStack) {
    THStack* stack = new THStack("stack", "");
    for (auto h : hists) stack->Add(h);
    stack->Draw("same");
    line->Draw("same");
  } else {
    line->Draw("same");
    for (auto p : hists) p->Draw("same");
  }
  canvas->SaveAs(saveName.c_str());
}

void dcaZ(vector<string> inputStrings, double v0min, double v0max)
{
  const int nDim       = 6;
  const int v0PtAxis   = 0;
  const int posPtAxis  = 1;
  const int negPtAxis  = 2;
  const int posDCAAxis = 3;
  const int negDCAAxis = 4;
  const int massAxis   = 5;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  string histName = "jet-v0qa/tracks/DCAzMass";
  histName += hadron;
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);
  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  TH2D* dcaZ = (TH2D*)thn->Projection(posDCAAxis, negDCAAxis);
  dcaZ->Draw("colz");
}

void radiusPtDiff(vector<string> inputStrings, double v0min, double v0max, bool normalise)
{
  const int nDim       = 5;
  const int v0PtAxis   = 0;
  const int posPtAxis  = 1;
  const int negPtAxis  = 2;
  const int v0RAxis    = 3;
  const int ptDiffAxis = 4;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];

  string histName = "jet-v0qa/tracks/V0Radius";
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);
  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  TH2D* v0R = (TH2D*)thn->Projection(ptDiffAxis, v0RAxis);
  if (isHistEmptyInRange(v0R, 1, v0R->GetNbinsX(), 1, v0R->GetNbinsY())) {
    cout << "Radius: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << "! Skipping" << endl;
    return;
  }

  string saveName = "daughterdPtRadius";
  saveName += TString::Format("_v0pt%.1f-%.1f", lowv0, highv0).Data();
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  if (normalise) {
    normaliseHistRowByRow(v0R);
    v0R->SetMaximum(1.);
    v0R->SetMinimum(1e-5);
    canvas->SetLogz();
  }

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string xTitle = "R_{V0} (cm)";
  string yTitle = "#it{p}_{T, pos} - #it{p}_{T, neg} (GeV/#it{c})";
  int lastRow = getLastFilledRow(v0R);
  int lastCol = getLastFilledCol(v0R);
  double dptMax = v0R->GetYaxis()->GetBinUpEdge(lastRow);
  double xMinFrame = v0R->GetXaxis()->GetXmin(), xMaxFrame = v0R->GetXaxis()->GetBinUpEdge(lastCol);
  double yMinFrame = -1. * dptMax, yMaxFrame = dptMax;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  canvas->cd();
  frame->Draw();
  v0R->Draw("same colz");
  canvas->SaveAs(saveName.c_str());
}

void radiusPt(vector<string> inputStrings, bool normalise)
{
  const int nDim       = 5;
  const int v0PtAxis   = 0;
  const int posPtAxis  = 1;
  const int negPtAxis  = 2;
  const int v0RAxis    = 3;
  const int ptDiffAxis = 4;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];

  string histName = "jet-v0qa/tracks/V0Radius";
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  TH2D* v0R = (TH2D*)thn->Projection(v0RAxis, v0PtAxis);
  if (isHistEmptyInRange(v0R, 1, v0R->GetNbinsX(), 1, v0R->GetNbinsY())) {
    cout << "Radius: empty hist! Skipping" << endl;
    return;
  }

  string saveName = "radiusV0Pt";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  if (normalise) {
    normaliseHistColByCol(v0R);
    v0R->SetMaximum(1.);
    v0R->SetMinimum(1e-5);
    canvas->SetLogz();
  }

  string xTitle = "#it{p}_{T, V0} (GeV/#it{c})";
  string yTitle = "R_{V0} (cm)";
  int lastCol = getLastFilledCol(v0R);
  double lastPt = v0R->GetXaxis()->GetBinUpEdge(lastCol);
  double xMinFrame = v0R->GetXaxis()->GetXmin(), xMaxFrame = lastPt;
  double yMinFrame = v0R->GetYaxis()->GetXmin(), yMaxFrame = v0R->GetYaxis()->GetXmax();
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  canvas->cd();
  frame->Draw();
  v0R->Draw("same colz");
  canvas->SaveAs(saveName.c_str());
}

void radiusMass(vector<string> inputStrings, double v0min, double v0max, bool normalise)
{
  const int nDim       = 5;
  const int v0PtAxis   = 0;
  const int posPtAxis  = 1;
  const int negPtAxis  = 2;
  const int v0RAxis    = 3;
  const int massAxis   = 4;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  string histName = "jet-v0qa/tracks/V0RadiusMass";
  histName += hadron;
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);
  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  TH2D* v0R = (TH2D*)thn->Projection(v0RAxis, massAxis);
  if (isHistEmptyInRange(v0R, 1, v0R->GetNbinsX(), 1, v0R->GetNbinsY())) {
    cout << "Radius: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << "! Skipping" << endl;
    return;
  }

  string saveName = "radius_m";
  saveName += hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowv0, highv0).Data();
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  if (normalise) {
    normaliseHistRowByRow(v0R);
    v0R->SetMaximum(1.);
    v0R->SetMinimum(1e-5);
    canvas->SetLogz();
  }

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string xTitle = "#it{M}(" + formatHadronDaughters(hadron) + ")";
  string yTitle = "R_{V0} (cm)";
  double xMinFrame = v0R->GetXaxis()->GetXmin(), xMaxFrame = v0R->GetXaxis()->GetXmax();
  double yMinFrame = v0R->GetYaxis()->GetXmin(), yMaxFrame = v0R->GetYaxis()->GetXmax();
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  canvas->cd();
  frame->Draw();
  v0R->Draw("same colz");
  canvas->SaveAs(saveName.c_str());
}

void radius(vector<string> inputStrings, double v0min, double v0max, bool normalise, bool doStack)
{
  const int nDim      = 5;
  const int v0PtAxis  = 0;
  const int posPtAxis = 1;
  const int negPtAxis = 2;
  const int v0RAxis   = 3;
  const int massAxis  = 4;

  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  string histName = "jet-v0qa/tracks/V0RadiusMass";
  histName += hadron;
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);
  double lowv0 = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);

  vector<double> radiusVals = {0., 10., 20., 30., 40., 50.};
  vector<TH1*> hists;
  TLegend* legend = CreateLegend(0.25, 0.85, 0.63, 0.9, "", 0.04);

  for (int i = 0; i < radiusVals.size() - 1; i++) {
    array<int, 2> rbins = getProjectionBins(thn->GetAxis(v0RAxis), radiusVals[i], radiusVals[i+1]);
    double lowr = thn->GetAxis(v0RAxis)->GetBinLowEdge(rbins[0]);
    double highr = thn->GetAxis(v0RAxis)->GetBinUpEdge(rbins[1]);

    thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
    thn->GetAxis(v0RAxis)->SetRange(rbins[0], rbins[1]);
    TH1D* mass = (TH1D*)thn->Projection(massAxis);
    string name = TString::Format("%s_R%.0f-%.0f_v0pt%.1f-%.1f", hadron.c_str(), lowr, highr, lowv0, highv0).Data();
    mass->SetName(name.c_str());
    if (isHistEmptyInRange(mass, 1, mass->GetNbinsX())) {
      cout << "Radius: empty hist for ptv0 " << lowv0 << " - " <<  highv0 << ", R " << lowr << " - " << highr << "! Skipping" << endl;
    } else {
      hists.push_back(mass);
      legend->AddEntry(mass, TString::Format("%.0f - %.0f cm", lowr, highr).Data(), "lf");
    }
  }
  if (hists.size() == 0) {
    cout << "Radius: no histograms filled for ptv0 " << lowv0 << " - " <<  highv0 << "! Skipping" << endl;
    return;
  }

  for (int i = 0; i < hists.size(); i++) {
    double alpha = doStack ? 0.5 : -1.;
    setStyle(hists[i], i+1, alpha);
    hists[i]->Rebin(4);
  }

  if (normalise) {
    double scale = getHistScale(hists, normalise, doStack);
    for (auto h : hists) h->Scale(1./scale, "width");
  }

  string saveName = "radius_m";
  saveName += hadron;
  saveName += TString::Format("_v0pt%.1f-%.1f", lowv0, highv0).Data();
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string xTitle = "M(" + formatHadronDaughters(hadron) + ")";
  string yTitle = "Counts";
  if (normalise) yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
  double xMinFrame = hists[0]->GetXaxis()->GetXmin(), xMaxFrame = hists[0]->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 2. * getHistScale(hists, normalise, doStack);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  double mass = (hadron == "K0S") ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, 0.55*yMaxFrame);
  setStyle(line, 0);

  canvas->cd();
  frame->Draw();
  legend->Draw("same");
  if (doStack) {
    THStack* stack = new THStack("stack", "");
    for (auto h : hists) stack->Add(h);
    stack->Draw("same");
    line->Draw("same");
  } else {
    line->Draw("same");
    for (auto p : hists) p->Draw("same");
  }
  canvas->SaveAs(saveName.c_str());
}


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


string getDataSet(int train)
{
  if (252064 == train) return "LHC22o_pass6";
  if (282430 == train) return "LHC22o_pass7_small";
  return "Could not find dataset";
}

void plotTrain(int train, string hadron, double v0min, double v0max, int setting)
{
  string inName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string dataSet = getDataSet(train);
  vector<string> inputStrings = {inName, dataSet, hadron};

  switch(setting) {
    case 0:
      ptAsymmetry(inputStrings, v0min, v0max, false);
      break;
    case 1:
      ptAsymmetry(inputStrings, true);
      break;
    case 2:
      mass(inputStrings, v0min, v0max, false, true);
      break;
    case 3:
      mass2D(inputStrings, v0min, v0max, true);
      break;
    case 4:
      trd(inputStrings, v0min, v0max, false, false);
      break;
    case 5:
      {
        vector<int> layers = { 1, 2, 3, 4, 5, 6, 7, 56, 57, 67, 567};
        for (auto layer : layers) {
          itslayer(inputStrings, v0min, v0max, false, false, layer);
        }
        // itslayer(inputStrings, v0min, v0max, false, true, 567);
      }
      break;
    case 6:
      // itsNCl(inputStrings, v0min, v0max, false, true);
      break;
    case 7:
      // itsChi2NCl(inputStrings, v0min, v0max, false, true);
      break;
    case 8:
      dcaXYQuadrants(inputStrings, v0min, v0max, false, true);
      break;
    case 9:
      dcaXY(inputStrings, v0min, v0max); // 2D version
      break;
    case 10:
      dcaXY(inputStrings, v0min, v0max, false, true, true);
      dcaXY(inputStrings, v0min, v0max, false, true, false);
      break;
    case 11:
      dcaZ(inputStrings, v0min, v0max); // 2D version
      break;
    case 12:
      radiusPtDiff(inputStrings, v0min, v0max, true);
      break;
    case 13:
      radiusPt(inputStrings, true);
      break;
    case 14:
      radiusMass(inputStrings, v0min, v0max, true);
      break;
    case 15:
      radius(inputStrings, v0min, v0max, false, true);
      break;
    default:
      cout << "Invalid setting!" << endl;
      break;
  }
}
void plotTrain(int train, string hadron, int setting)
{
  gROOT->SetBatch(kTRUE);
  vector<double> pt   = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0};
  vector<double> bins = pt;
  for (int i = 0; i < bins.size() - 1; i++) {
    plotTrain(train, hadron, bins[i], bins[i+1], setting);
  }
}

void plot282430(string hadron, double v0min, double v0max, int setting)
{
  plotTrain(282430, hadron, v0min, v0max, setting);
}
void plot282430(string hadron, int setting)
{
  plotTrain(282430, hadron, setting);
}
