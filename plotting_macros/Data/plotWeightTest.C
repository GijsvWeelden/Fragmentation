
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

double roundToNextPowerOfTen(double x)
{
  return std::pow(10., std::ceil(std::log10(x)));
}
double getMeanIgnoringFirstBin(TH1* hist)
{
  double sum = 0;
  double n = 0;
  for (int i = 2; i <= hist->GetNbinsX(); ++i) {
    n += hist->GetBinContent(i);
    sum += n * hist->GetBinLowEdge(i);
  }
  return sum / n;
}
// Compare filling inclusive V0 hists with counts vs purity weights
void plotFillComparisonIncl(double p)
{
  int q = std::round(p*10.);
  string fileName = TString::Format("../../inputfiles/weightTest/AnalysisResults-incl-p%.02d.root", q).Data();
  string nAccHistName = "jet-fragmentation/data/V0/V0PtEtaPhi";
  string nWeiHistName = "jet-fragmentation/data/V0/V0PtEtaPhiWeighted";

  TFile* file = TFile::Open(fileName.c_str());
  TH3* countHist = (TH3*)file->Get(nAccHistName.c_str());
  TH1* nAccHist = (TH1*)countHist->ProjectionX("nAccHist");
  setStyle(nAccHist, 0);
  TH3* weightHist = (TH3*)file->Get(nWeiHistName.c_str());
  TH1* nWeiHist = (TH1*)weightHist->ProjectionX("nWeiHist");
  nWeiHist->Sumw2();
  setStyle(nWeiHist, 1);

  string saveName = "weightTest-V0pt-p";
  saveName += TString::Format("%.02d", q).Data();
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 800, 600);
  canvas->SetLogy();

  TLegend* legend = CreateLegend(0.45, 0.9, 0.6, 0.85, "", 0.04);
  legend->AddEntry(nAccHist, "Count V0s", "l");
  legend->AddEntry(nWeiHist, "Weight V0s (weighted)", "l");

  // We should have nAcc * p = nWei for flat p
  double nAcc = nAccHist->Integral(1, nAccHist->GetNbinsX());
  double nWei = nWeiHist->Integral(1, nWeiHist->GetNbinsX());
  double avgWeight = nWei / nAcc;
  string frameTitle = TString::Format("Avg weight %.2f / %.0f = %.2f, expected p = %.1f", nWei, nAcc, avgWeight, p).Data();

  double xMinFrame = 0;
  double xMaxFrame = 20;
  double yMinFrame = 1e-1;
  double yMaxFrame = roundToNextPowerOfTen(nAccHist->GetMaximum());
  TH1* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame);
  frame->SetXTitle("p_{T, V0} (GeV/c)");
  frame->SetYTitle("Counts or weights");
  frame->SetTitle(frameTitle.c_str());

  frame->Draw();
  nAccHist->Draw("same");
  nWeiHist->Draw("same p");
  legend->Draw("same");
  canvas->SaveAs(canvas->GetName());
}
// Compare filling nV0 hists with counts vs purity weights
void plotFillComparisonV0inEvents(double p)
{
  int q = std::round(p*10.);
  string fileName = TString::Format("../../inputfiles/weightTest/AnalysisResults-incl-p%.02d.root", q).Data();
  string nV0sHistName = "jet-fragmentation/data/V0/nV0sEvent";
  string nAccHistName = "jet-fragmentation/data/V0/nV0sEventAcc"; // Filters out rejected candidates
  string nWeiHistName = "jet-fragmentation/data/V0/nV0sEventAccWeighted"; // Purity weight applied

  TFile* file = TFile::Open(fileName.c_str());
  TH1* nV0sHist = (TH1*)file->Get(nV0sHistName.c_str());
  TH1* nAccHist = (TH1*)file->Get(nAccHistName.c_str());
  TH1* nWeiHist = (TH1*)file->Get(nWeiHistName.c_str());
  setStyle(nV0sHist, 0);
  setStyle(nAccHist, 1);
  setStyle(nWeiHist, 2);

  string saveName = "weightTest-nV0inEvents-p";
  saveName += TString::Format("%.02d", q).Data();
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 800, 600);
  canvas->SetLogy();

  TLegend* legend = CreateLegend(0.45, 0.9, 0.6, 0.85, "", 0.04);
  legend->AddEntry(nV0sHist, "All V0s", "l");
  legend->AddEntry(nAccHist, "Acc. V0s", "l");
  legend->AddEntry(nWeiHist, "Acc. V0s (weighted)", "l");

  // We should have mAcc * p = mWei for flat p
  double mAcc = nAccHist->GetMean(1);
  double mWei = nWeiHist->GetMean(1);
  double meanRatio = mWei / mAcc;
  string frameTitle = TString::Format("Ratio of mean %.2f / %.2f = %.2f, expected p = %.1f", mWei, mAcc, meanRatio, p).Data();

  double xMinFrame = 0;
  double xMaxFrame = 10;
  double yMinFrame = 1e-1;
  double yMaxFrame = roundToNextPowerOfTen(nV0sHist->GetMaximum());
  TH1* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame);
  frame->SetXTitle("V0s in event");
  frame->SetYTitle("Counts or weights");
  frame->SetTitle(frameTitle.c_str());

  frame->Draw();
  nV0sHist->Draw("same");
  nAccHist->Draw("same");
  nWeiHist->Draw("same");
  legend->Draw("same");
  canvas->SaveAs(canvas->GetName());
}
// Compare filling jet hists with counts vs state weights
void plotFillComparisonJets(double p)
{
  int q = std::round(p*10.);
  string fileName = TString::Format("../../inputfiles/weightTest/AnalysisResults-frag-p%.02d.root", q).Data();
  string cHistName = "jet-fragmentation/data/jets/jetPtEtaPhi";
  string wHistName = "jet-fragmentation/data/jets/weighted/jetPtEtaPhi";

  TFile* file = TFile::Open(fileName.c_str());
  TH3* countHist = (TH3*)file->Get(cHistName.c_str());
  TH3* weightHist = (TH3*)file->Get(wHistName.c_str());
  TH1* cHist = (TH1*)countHist->ProjectionX("cHist");
  TH1* wHist = (TH1*)weightHist->ProjectionX("wHist");
  setStyle(cHist, 0);
  setStyle(wHist, 1);

  string saveName = "weightTest-jetpt-p";
  saveName += TString::Format("%.02d", q).Data();
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 800, 600);
  canvas->SetLogy();

  TLegend* legend = CreateLegend(0.5, 0.7, 0.6, 0.85, "", 0.04);
  legend->AddEntry(cHist, "Count jets", "l");
  legend->AddEntry(wHist, "Weight jets", "l");

  // We should have counts = int(weights)
  double nJets = cHist->Integral(0, cHist->GetNbinsX()+1);
  double wJets = wHist->Integral(0, wHist->GetNbinsX()+1);
  double wJet0 = wHist->Integral(0, 0);
  string frameTitle = TString::Format("N_{jets} = %.0f, weights = %.2f (%.2f = %.2f%% #it{p}_{T, jet} < 0 GeV/c)", nJets, wJets, wJet0, wJet0/wJets).Data();

  double xMinFrame = 0;
  double xMaxFrame = 50;
  double yMinFrame = 1e-1;
  double yMaxFrame = roundToNextPowerOfTen(cHist->GetMaximum());
  TH1* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame);
  frame->SetXTitle("p_{T, jet} (GeV/c)");
  frame->SetYTitle("Counts or weights");
  frame->SetTitle(frameTitle.c_str());

  frame->Draw();
  cHist->Draw("same");
  wHist->Draw("same");
  legend->Draw("same");
  canvas->SaveAs(canvas->GetName());
}
// Compare filling in-jet V0 hists with counts vs signal weights
void plotFillComparisonV0inJets(double p)
{
  int q = std::round(p*10.);
  string fileName = TString::Format("../../inputfiles/weightTest/AnalysisResults-frag-p%.02d.root", q).Data();
  string cHistName = "jet-fragmentation/data/jets/V0/jetPtnV0nK0SnLambdanAntiLambda";
  string wHistName = "jet-fragmentation/data/jets/weighted/V0/jetPtnV0nK0SnLambdanAntiLambda";

  TFile* file = TFile::Open(fileName.c_str());
  THnSparse* countHist = (THnSparse*)file->Get(cHistName.c_str());
  TH1* cHist = (TH1*)countHist->Projection(1, "cHist");
  setStyle(cHist, 0);
  THnSparse* weightHist = (THnSparse*)file->Get(wHistName.c_str());
  TH1* wHist = (TH1*)weightHist->Projection(1, "wHist");
  setStyle(wHist, 1);

  string saveName = "weightTest-nV0inJets-p";
  saveName += TString::Format("%.02d", q).Data();
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 800, 600);
  canvas->SetLogy();

  TLegend* legend = CreateLegend(0.5, 0.7, 0.6, 0.85, "", 0.04);
  legend->AddEntry(cHist, "Count V0s", "l");
  legend->AddEntry(wHist, "Weight V0s", "l");

  // We should have counts * p = int(weights) for flat p
  double mCou = cHist->GetMean(1);
  double mWei = wHist->GetMean(1);
  double meanRatio = mWei / mCou;
  string frameTitle = TString::Format("Ratio of mean %.2f / %.2f = %.2f, expected %.1f", mWei, mCou, meanRatio, p).Data();

  double xMinFrame = 0;
  double xMaxFrame = 5;
  double yMinFrame = 1e-1;
  double yMaxFrame = roundToNextPowerOfTen(cHist->GetMaximum());
  TH1* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame);
  frame->SetXTitle("V0s in jet");
  frame->SetYTitle("Counts or weights");
  frame->SetTitle(frameTitle.c_str());

  frame->Draw();
  cHist->Draw("same");
  wHist->Draw("same");
  legend->Draw("same");
  canvas->SaveAs(canvas->GetName());
}
