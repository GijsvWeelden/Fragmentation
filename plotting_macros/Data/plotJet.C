
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

double getNjets(TFile* inFile, double jetptmin, double jetptmax)
{
  string histName = "jet-fragmentation/data/jets/jetPtEtaPhi";
  TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
  std::array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
  return jets->Integral(jetptbins[0], jetptbins[1], 1, jets->GetNbinsY(), 1, jets->GetNbinsZ());
}

// -----------------------------------------------------------------------------

// Plot jet pt of all jets (without V0s, with V0s (no corrections))
double plotJetPt(vector<string> inputStrings, double ptmin, double ptmax)
{
  string fileName = inputStrings[0];
  string dataSet  = inputStrings[1];
  double njets = 0;

  TFile* inFile = TFile::Open(fileName.c_str());
  THnSparse* hNchjets = (THnSparse*)inFile->Get("jet-finder-data-charged/hJet");
  hNchjets->Sumw2();
  THnSparse* hNv0jets = (THnSparse*)inFile->Get("jet-finder-v0-data-charged/hJet");
  hNv0jets->Sumw2();

  TH1* chjets = (TH1*)hNchjets->Projection(1);
  TH1* v0jets = (TH1*)hNv0jets->Projection(1);
  array<int, 2> chbins = getProjectionBins(chjets->GetXaxis(), ptmin, ptmax);
  array<int, 2> v0bins = getProjectionBins(v0jets->GetXaxis(), ptmin, ptmax);
  njets = chjets->Integral(chbins[0], chbins[1]);

  string saveName = "chJetPtV0JetPt.pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  string xTitle = "#it{p}_{T, jet} (GeV/#it{c})";
  string yTitle = "counts";
  // double xMinFrame = chjets->GetXaxis()->GetXmin(), xMaxFrame = chjets->GetXaxis()->GetXmax();
  double xMinFrame = chjets->GetXaxis()->GetXmin(), xMaxFrame = 100.;
  double yMinFrame = 0., yMaxFrame = 1.5 * getHistUpperBound(chjets, true);
  yMinFrame = -1. * yMaxFrame;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  TLegend* legend = CreateLegend(0.45, 0.85, 0.62, 0.9, "", 0.04);
  setStyle(chjets, 0); legend->AddEntry(chjets, "Ch jets (incl)");
  setStyle(v0jets, 1); legend->AddEntry(v0jets, "Ch+V0 jets (incl)");

  canvas->cd();
  frame->Draw();
  legend->Draw("same");
  chjets->Draw("same hist");
  v0jets->Draw("same hist");
  canvas->SaveAs(saveName.c_str());
  return njets;
}
void plotJetPtWithV0s(string inName = "AnalysisResults.root", bool doRatio = false)
{
  const int nDim = 5;
  const int jetptAxis = 0;
  const int nv0Axis = 1;
  const int nK0SAxis = 2;
  const int nLambdaAxis = 3;
  const int nAntiLambdaAxis = 4;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 5e-6, yMaxFrame = 2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, ch. jet}";
  yTitle = "normalised count";
  if (doRatio) { yTitle = "jets / all jets"; }
  latexText = "LHC22o_apass6_minBias";

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtnV0nK0SnLambdanAntiLambda";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  TH1D* jetpt = (TH1D*)thn->Projection(jetptAxis);
  TH1D* denominator = (TH1D*)jetpt->Clone("denominator");
  double nJets = jetpt->Integral();
  if (doRatio) { jetpt->Divide(denominator); }
  else { jetpt->Scale(1./nJets); }
  setStyle(jetpt, 0);
  jetpt->SetName("jetpt");
  legend->AddEntry(jetpt, "All jets");
  histVector.push_back(jetpt);

  thn->GetAxis(nv0Axis)->SetRange(1, 1);
  TH1D* jetptNoV0s = (TH1D*)thn->Projection(jetptAxis);
  if (doRatio) { jetptNoV0s->Divide(denominator); }
  else { jetptNoV0s->Scale(1./nJets); }
  setStyle(jetptNoV0s, 1);
  jetptNoV0s->SetName("jetptNoV0s");
  legend->AddEntry(jetptNoV0s, "Jets w.o. V0s");
  histVector.push_back(jetptNoV0s);

  thn->GetAxis(nv0Axis)->SetRange(2, thn->GetAxis(nv0Axis)->GetNbins());
  TH1D* jetPtWithV0s = (TH1D*)thn->Projection(jetptAxis);
  if (doRatio) { jetPtWithV0s->Divide(denominator); }
  else { jetPtWithV0s->Scale(1./nJets); }
  jetPtWithV0s->SetName("jetPtWithV0s");
  setStyle(jetPtWithV0s, 2);
  legend->AddEntry(jetPtWithV0s, "Jets with V0s");
  histVector.push_back(jetPtWithV0s);

  saveName = "jetspectranV0s";
  if (doRatio) { saveName = "jetspectranV0s_ratio"; }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotJetPtWithV0sId(string inName = "AnalysisResults.root", bool doRatio = false)
{
  const int nDim = 5;
  const int jetptAxis = 0;
  const int nv0Axis = 1;
  const int nK0SAxis = 2;
  const int nLambdaAxis = 3;
  const int nAntiLambdaAxis = 4;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 5e-6, yMaxFrame = 2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, ch. jet}";
  yTitle = "normalised count";
  if (doRatio) { yTitle = "jets / all jets"; }
  latexText = "LHC22o_apass6_minBias";

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtnV0nK0SnLambdanAntiLambda";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  TH1D* jetpt = (TH1D*)thn->Projection(jetptAxis);
  TH1D* denominator = (TH1D*)jetpt->Clone("denominator");
  double nJets = jetpt->Integral();
  if (doRatio) { jetpt->Divide(denominator); }
  else { jetpt->Scale(1./nJets); }
  setStyle(jetpt, 0);
  jetpt->SetName("jetpt");
  legend->AddEntry(jetpt, "All jets");
  histVector.push_back(jetpt);

  thn->GetAxis(nv0Axis)->SetRange(1, 1);
  TH1D* jetptNoV0s = (TH1D*)thn->Projection(jetptAxis);
  if (doRatio) { jetptNoV0s->Divide(denominator); }
  else { jetptNoV0s->Scale(1./nJets); }
  setStyle(jetptNoV0s, 1);
  jetptNoV0s->SetName("jetptNoV0s");
  legend->AddEntry(jetptNoV0s, "Jets w.o. V0s");
  histVector.push_back(jetptNoV0s);

  thn->GetAxis(nv0Axis)->SetRange(2, thn->GetAxis(nv0Axis)->GetNbins());
  TH1D* jetPtWithV0s = (TH1D*)thn->Projection(jetptAxis);
  if (doRatio) { jetPtWithV0s->Divide(denominator); }
  else { jetPtWithV0s->Scale(1./nJets); }
  jetPtWithV0s->SetName("jetPtWithV0s");
  setStyle(jetPtWithV0s, 2);
  legend->AddEntry(jetPtWithV0s, "Jets with V0s");
  histVector.push_back(jetPtWithV0s);

  // K0S
  thn->GetAxis(nK0SAxis)->SetRange(2, thn->GetAxis(nK0SAxis)->GetNbins());
  TH1D* jetPtWithK0S = (TH1D*)thn->Projection(jetptAxis);
  if (doRatio) { jetPtWithK0S->Divide(denominator); }
  else { jetPtWithK0S->Scale(1./nJets); }
  jetPtWithK0S->SetName("jetPtWithK0S");
  setStyle(jetPtWithK0S, 3);
  legend->AddEntry(jetPtWithK0S, TString::Format("Jets with %s", formatHadronName("K0S").c_str()).Data());
  histVector.push_back(jetPtWithK0S);

  // Lambda
  thn->GetAxis(nLambdaAxis)->SetRange(2, thn->GetAxis(nLambdaAxis)->GetNbins());
  TH1D* jetPtWithLambda = (TH1D*)thn->Projection(jetptAxis);
  if (doRatio) { jetPtWithLambda->Divide(denominator); }
  else { jetPtWithLambda->Scale(1./nJets); }
  jetPtWithLambda->SetName("jetPtWithLambda");
  setStyle(jetPtWithLambda, 4);
  legend->AddEntry(jetPtWithLambda, TString::Format("Jets with %s", formatHadronName("Lambda0").c_str()).Data());
  histVector.push_back(jetPtWithLambda);

  // AntiLambda
  thn->GetAxis(nAntiLambdaAxis)->SetRange(2, thn->GetAxis(nAntiLambdaAxis)->GetNbins());
  TH1D* jetPtWithAntiLambda = (TH1D*)thn->Projection(jetptAxis);
  if (doRatio) { jetPtWithAntiLambda->Divide(denominator); }
  else { jetPtWithAntiLambda->Scale(1./nJets); }
  jetPtWithAntiLambda->SetName("jetPtWithAntiLambda");
  setStyle(jetPtWithAntiLambda, 5);
  legend->AddEntry(jetPtWithAntiLambda, TString::Format("Jets with %s", formatHadronName("AntiLambda0").c_str()).Data());
  histVector.push_back(jetPtWithAntiLambda);

  saveName = "jetspectranV0s";
  if (doRatio) { saveName = "jetspectranV0s_ratio"; }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}

// -----------------------------------------------------------------------------

void plotJetPt(string inName = "AnalysisResults.root")
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 5e-6, yMaxFrame = 2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, ch. jet}";
  yTitle = "normalised count";
  latexText = "LHC22o_apass6_minBias";

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/jetPtEtaPhi";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());
  TH1D* jetpt = (TH1D*)th3->ProjectionX("jetpt");
  jetpt->Scale(1./jetpt->Integral());
  setStyle(jetpt, 0);
  histVector.push_back(jetpt);

  saveName = "jetpt";
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotJetEta(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = false;
  double xMinFrame = -1., xMaxFrame = 1., yMinFrame = 5e-6, yMaxFrame = 0.2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#eta_{jet}";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{LHC22o_apass6_minBias_small}{#it{p}_{T, ch. jet} = %.0f - %.0f GeV/c}", jetptmin, jetptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/jetPtEtaPhi";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  int firstBinPhi = 1, lastBinPhi = th3->GetNbinsZ();
  firstBinPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
  TH1D* jeteta = (TH1D*)th3->ProjectionY("jeteta", firstBinPt, lastBinPt, firstBinPhi, lastBinPhi);
  jeteta->Scale(1./jeteta->Integral());
  setStyle(jeteta, 0);
  histVector.push_back(jeteta);

  saveName = "jeteta";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotJetPhi(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = false;
  double xMinFrame = 0., xMaxFrame = 6.3, yMinFrame = 0, yMaxFrame = 0.02;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#varphi_{jet}";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{LHC22o_apass6_minBias_small}{#it{p}_{T, ch. jet} = %.0f - %.0f GeV/c}", jetptmin, jetptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/jetPtEtaPhi";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  int firstBinEta = 1, lastBinEta = th3->GetNbinsY();
  firstBinPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinPt, lastBinPt);
  TH1D* jetphi = (TH1D*)th3->ProjectionZ("jetphi", firstBinPt, lastBinPt, firstBinEta, lastBinEta);
  jetphi->Scale(1./jetphi->Integral());
  setStyle(jetphi, 0);
  histVector.push_back(jetphi);

  saveName = "jetphi";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotJetEtaPhi(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  gStyle->SetPalette(77);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = false;
  double xMinFrame = -.6, xMaxFrame = .6, yMinFrame = 0., yMaxFrame = 6.3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#eta_{jet}";
  yTitle = "#varphi_{jet}";
  latexText = TString::Format("LHC22o_apass6_minBias_small, #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c", jetptmin, jetptmax).Data();

  std::vector<TH2D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/jetPtEtaPhi";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  firstBinPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinPt, lastBinPt);
  TH2D* jetetaphi = (TH2D*)th3->Project3D("zy");
  jetetaphi->SetName("jetetaphi");
  setStyle(jetetaphi, 0);
  jetetaphi->Scale(1./jetetaphi->Integral());
  // normaliseHistRowByRow(jetetaphi);
  // jetetaphi->SetMinimum(1e-5);
  // jetetaphi->SetMaximum(1.);

  saveName = "jetetaphi";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  jetetaphi->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

void plotTrain(int train, int setting, double jetmin, double jetmax)
{
  string inName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string dataSet = getDataSet(train);
  vector<string> inputStrings = {inName, dataSet};

  switch(setting) {
    case 0:
      {
        if (train < 350e3) cout << "Warning: unexpected results from nJets < 0!" << endl;
        double njets = plotJetPt(inputStrings, jetmin, jetmax);
        cout << "Njets = " << njets << endl;
      }
      break;
    default:
      cout << "Invalid setting!" << endl;
  }
}

void plot349871(int setting, double jetmin = -1., double jetmax = 1e6)
{ plotTrain(349871, setting, jetmin, jetmax); }
