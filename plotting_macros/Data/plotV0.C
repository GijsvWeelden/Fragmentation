
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

// const double MassK0S = 0.497611;
// const double MassLambda0 = 1.115683;

void plotV0Pt(string inName = "AnalysisResults.root", string dataSet = "dataSet")
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 60, yMinFrame = 1e-8, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, V0}";
  yTitle = "normalised count";
  latexText = dataSet;

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtEtaPhi";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());
  TH1D* v0pt = (TH1D*)th3->ProjectionX("v0pt");
  v0pt->Scale(1./v0pt->Integral());
  setStyle(v0pt, 0);
  histVector.push_back(v0pt);

  saveName = "v0pt";
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0Eta(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = false;
  double xMinFrame = -1., xMaxFrame = 1., yMinFrame = 0, yMaxFrame = 0.2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#eta_{V0}";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtEtaPhi";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  int firstBinPhi = 1, lastBinPhi = th3->GetNbinsZ();
  firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
  TH1D* v0eta = (TH1D*)th3->ProjectionY("v0eta", firstBinPt, lastBinPt, firstBinPhi, lastBinPhi);
  v0eta->Scale(1./v0eta->Integral());
  setStyle(v0eta, 0);
  histVector.push_back(v0eta);

  saveName = "v0eta";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0Phi(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
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
  xTitle = "#varphi_{V0}";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtEtaPhi";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  int firstBinEta = 1, lastBinEta = th3->GetNbinsY();
  firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinPt, lastBinPt);
  TH1D* v0phi = (TH1D*)th3->ProjectionZ("v0phi", firstBinPt, lastBinPt, firstBinEta, lastBinEta);
  v0phi->Scale(1./v0phi->Integral());
  setStyle(v0phi, 0);
  histVector.push_back(v0phi);

  saveName = "v0phi";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0EtaPhi(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  gStyle->SetPalette(77);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = false;
  double xMinFrame = -1., xMaxFrame = 1., yMinFrame = 0., yMaxFrame = 6.3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#eta_{V0}";
  yTitle = "#varphi_{V0}";
  latexText = TString::Format("%s, #it{p}_{T, V0} = %.1f - %.1f GeV/c", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH2D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtEtaPhi";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinPt, lastBinPt);
  TH2D* v0etaphi = (TH2D*)th3->Project3D("zy");
  v0etaphi->SetName("v0etaphi");
  setStyle(v0etaphi, 0);
  v0etaphi->Scale(1./v0etaphi->Integral());
  // normaliseHistRowByRow(v0etaphi);
  // v0etaphi->SetMinimum(1e-5);
  // v0etaphi->SetMaximum(1.);

  saveName = "v0etaphi";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  v0etaphi->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

// ----------------------------------------------------------

void plotV0Radius(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = false;
  double xMinFrame = 0., xMaxFrame = 100., yMinFrame = 0, yMaxFrame = 0.5;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "V0 Radius";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtRadiusCosPA";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  int firstBinCosPA = 1, lastBinCosPA = th3->GetNbinsZ();
  firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
  // th3->GetXaxis()->SetRange(firstBinPt, lastBinPt);
  TH1D* v0Radius = (TH1D*)th3->ProjectionY("v0Radius", firstBinPt, lastBinPt, firstBinCosPA, lastBinCosPA);
  v0Radius->Scale(1./v0Radius->Integral());
  setStyle(v0Radius, 0);
  histVector.push_back(v0Radius);

  saveName = "v0Radius";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0CosPA(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = false;
  double xMinFrame = 0.95, xMaxFrame = 1.01, yMinFrame = 0, yMaxFrame = 0.5;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "V0 cos(PA)";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtRadiusCosPA";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  int firstBinRadius = 1, lastBinRadius = th3->GetNbinsY();
  firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
  TH1D* v0CosPA = (TH1D*)th3->ProjectionZ("v0CosPA", firstBinPt, lastBinPt, firstBinRadius, lastBinRadius);
  v0CosPA->Scale(1./v0CosPA->Integral());
  setStyle(v0CosPA, 0);
  histVector.push_back(v0CosPA);

  saveName = "v0CosPA";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0DCAdaughters(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = false;
  double xMinFrame = 0., xMaxFrame = 1., yMinFrame = 0., yMaxFrame = 0.6;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "DCA daughters";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtDCAd";
  TFile *inFile = TFile::Open(inName.c_str());
  TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th2->GetNbinsX();
  firstBinPt = th2->GetXaxis()->FindBin(v0ptmin + 1e-3);
  lastBinPt  = th2->GetXaxis()->FindBin(v0ptmax - 1e-3);
  th2->GetXaxis()->SetRange(firstBinPt, lastBinPt);
  TH1D* v0DCAdaughters = (TH1D*)th2->ProjectionY("v0DCAdaughters", firstBinPt, lastBinPt);
  v0DCAdaughters->Scale(1./v0DCAdaughters->Integral());
  setStyle(v0DCAdaughters, 0);
  histVector.push_back(v0DCAdaughters);

  saveName = "v0DCAdaughters";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0DCApos(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = -10., xMaxFrame = 10., yMinFrame = 1e-5, yMaxFrame = 0.5;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "DCA pos";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  int firstBinDCA = 1, lastBinDCA = th3->GetNbinsZ();
  firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
  TH1D* v0DCApos = (TH1D*)th3->ProjectionY("v0DCApos", firstBinPt, lastBinPt, firstBinDCA, lastBinDCA);
  v0DCApos->Scale(1./v0DCApos->Integral());
  setStyle(v0DCApos, 0);
  histVector.push_back(v0DCApos);

  saveName = "v0DCApos";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0DCAneg(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = -10., xMaxFrame = 10., yMinFrame = 1e-5, yMaxFrame = 0.5;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "DCA neg";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  int firstBinDCA = 1, lastBinDCA = th3->GetNbinsY();
  firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
  TH1D* v0DCAneg = (TH1D*)th3->ProjectionZ("v0DCAneg", firstBinPt, lastBinPt, firstBinDCA, lastBinDCA);
  v0DCAneg->Scale(1./v0DCAneg->Integral());
  setStyle(v0DCAneg, 0);
  histVector.push_back(v0DCAneg);

  saveName = "v0DCAneg";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0DCAposneg(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = false;
  double xMinFrame = -10., xMaxFrame = 10., yMinFrame = -10., yMaxFrame = 10.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "DCA pos";
  yTitle = "DCA neg";
  latexText = TString::Format("%s, #it{p}_{T, V0} = %.1f - %.1f GeV/c", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
  firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
  lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinPt, lastBinPt);
  TH2D* v0DCAposneg = (TH2D*)th3->Project3D("zy");
  v0DCAposneg->SetName("v0DCAposneg");
  v0DCAposneg->Scale(1./v0DCAposneg->Integral());
  setStyle(v0DCAposneg, 0);

  saveName = "v0DCAposneg";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  v0DCAposneg->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0ctau(string inName = "AnalysisResults.root", string dataSet = "dataSet", int setting = 1, double v0ptmin = 0., double v0ptmax = 100.)
{
  const int nDim           = 4;
  const int ptAxis         = 0;
  const int K0SAxis        = 1;
  const int LambdaAxis     = 2;
  const int AntiLambdaAxis = 3;

  double time = clock();
  gStyle->SetNdivisions(505);
  if (setting < K0SAxis || setting > AntiLambdaAxis) {
    cout << "Invalid setting! 1 for K0S, 2 for Lambda, 3 for antiLambda" << endl;
    return;
  }

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 40., yMinFrame = 1e-5, yMaxFrame = .1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{c}#tau (K^{0}_{S})";
  if (setting == LambdaAxis) { xTitle = "#it{c}#tau (#Lambda)"; }
  if (setting == AntiLambdaAxis) { xTitle = "#it{c}#tau (#bar{#Lambda})"; }
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;


  histName = "jet-fragmentation/data/V0/V0PtCtau";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptAxis)->FindBin(v0ptmax - 1e-3);
  thn->GetAxis(ptAxis)->SetRange(firstBinPt, lastBinPt);
  TH1D* ctau = (TH1D*)thn->Projection(setting);
  ctau->SetName(TString::Format("ctau_%d", setting).Data());
  ctau->Scale(1./ctau->Integral());
  setStyle(ctau, 0);
  histVector.push_back(ctau);

  saveName = "ctauK0S";
  if (setting == LambdaAxis) { saveName = "ctauLambda"; }
  else if (setting == AntiLambdaAxis) { saveName = "ctauAntiLambda"; }
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  TCanvas* myCanvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
// This plot isn't very useful. ctauL = 2*ctauK0S?
void plotV0ctauKL(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  const int nDim           = 4;
  const int ptAxis         = 0;
  const int K0SAxis        = 1;
  const int LambdaAxis     = 2;
  const int AntiLambdaAxis = 3;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = 0., xMaxFrame = 40., yMinFrame = 0., yMaxFrame = 40.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{c}#tau (K^{0}_{S})";
  yTitle = "#it{c}#tau (#Lambda)";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtCtau";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptAxis)->FindBin(v0ptmax - 1e-3);
  thn->GetAxis(ptAxis)->SetRange(firstBinPt, lastBinPt);
  TH2D* ctau = (TH2D*)thn->Projection(LambdaAxis, K0SAxis);
  ctau->SetName("ctau_KL");
  ctau->Scale(1./ctau->Integral());
  setStyle(ctau, 0);

  saveName = "ctauK0SLambda";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  ctau->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0mass(string inName = "AnalysisResults.root", string dataSet = "dataSet", int setting = 1, double v0ptmin = 0., double v0ptmax = 100.)
{
  const int nDim           = 4;
  const int ptAxis         = 0;
  const int K0SAxis        = 1;
  const int LambdaAxis     = 2;
  const int AntiLambdaAxis = 3;

  double time = clock();
  gStyle->SetNdivisions(505);
  if (setting < K0SAxis || setting > AntiLambdaAxis) {
    cout << "Invalid setting! 1 for K0S, 2 for Lambda, 3 for antiLambda" << endl;
    return;
  }

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1e-4, yMaxFrame = .1;
  if (setting == LambdaAxis || setting == AntiLambdaAxis) {
    xMinFrame = 1.015, xMaxFrame = 1.215;
  }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{M} (K^{0}_{S})";
  if (setting == LambdaAxis) { xTitle = "#it{M} (#Lambda)"; }
  if (setting == AntiLambdaAxis) { xTitle = "#it{M} (#bar{#Lambda})"; }
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{%s}{#it{p}_{T, V0} = %.1f - %.1f GeV/c}", dataSet.c_str(), v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptAxis)->FindBin(v0ptmax - 1e-3);
  thn->GetAxis(ptAxis)->SetRange(firstBinPt, lastBinPt);
  TH1D* mass = (TH1D*)thn->Projection(setting);
  mass->SetName(TString::Format("mass_%d", setting).Data());
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);
  histVector.push_back(mass);

  saveName = "massK0S";
  if (setting == LambdaAxis) { saveName = "massLambda"; }
  else if (setting == AntiLambdaAxis) { saveName = "massAntiLambda"; }
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0massKL(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  const int nDim           = 4;
  const int ptAxis         = 0;
  const int K0SAxis        = 1;
  const int LambdaAxis     = 2;
  const int AntiLambdaAxis = 3;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{M}("; xTitle += formatHadronDaughters("K0S").c_str(); xTitle += ") (GeV/#it{c}^{2})";
  yTitle = "#it{M}("; yTitle += formatHadronDaughters("Lambda0").c_str(); yTitle += ") (GeV/#it{c}^{2})";

  histName = "jet-fragmentation/data/V0/V0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  array<int, 2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), v0ptmin, v0ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
  TH2D* mass = (TH2D*)thn->Projection(LambdaAxis, K0SAxis);
  mass->SetName("mass_KL");
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);

  string canvasName = "massK0SLambda";
  canvasName += TString::Format("_v0pt%.1f-%.1f", lowpt, highpt);
  canvasName += ".pdf";
  TCanvas* myCanvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), xCanvas, yCanvas);
  myCanvas->SetLogz();
  myCanvas->cd();

  // double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1.015, yMaxFrame = 1.215;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1.075, yMaxFrame = 1.215;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);

  TLine* lK0S = new TLine(MassK0S, yMinFrame, MassK0S, yMaxFrame);
  lK0S->SetLineWidth(2);
  lK0S->SetLineStyle(9);
  lK0S->SetLineColor(GetColor(0));
  TLine* lLambda = new TLine(xMinFrame, MassLambda0, xMaxFrame, MassLambda0);
  lLambda->SetLineWidth(2);
  lLambda->SetLineStyle(9);
  lLambda->SetLineColor(GetColor(0));

  frame->Draw();
  mass->Draw("same colz");
  lK0S->Draw("same");
  lLambda->Draw("same");
  latexText = TString::Format("%s, #it{p}_{T, V0} = %.1f - %.1f GeV/c", dataSet.c_str(), lowpt, highpt).Data();
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(myCanvas->GetName());
}
void plotV0massLL(string inName = "AnalysisResults.root", string dataSet = "dataSet", double v0ptmin = 0., double v0ptmax = 100.)
{
  const int nDim           = 4;
  const int ptAxis         = 0;
  const int K0SAxis        = 1;
  const int LambdaAxis     = 2;
  const int AntiLambdaAxis = 3;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{M}("; xTitle += formatHadronDaughters("Lambda0").c_str(); xTitle += ") (GeV/#it{c}^{2})";
  yTitle = "#it{M}("; yTitle += formatHadronDaughters("AntiLambda0").c_str(); yTitle += ") (GeV/#it{c}^{2})";

  histName = "jet-fragmentation/data/V0/V0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  array<int, 2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), v0ptmin, v0ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
  TH2D* mass = (TH2D*)thn->Projection(AntiLambdaAxis, LambdaAxis);
  mass->SetName("mass_LL");
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);

  string canvasName = "massLambdaAntiLambda";
  canvasName += TString::Format("_v0pt%.0f-%.0f", v0ptmin, v0ptmax);
  canvasName += ".pdf";
  TCanvas* myCanvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), xCanvas, yCanvas);
  myCanvas->SetLogz();
  myCanvas->cd();

  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 1.015, yMaxFrame = 1.215;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLine* lLambda = new TLine(MassLambda0, yMinFrame, MassLambda0, yMaxFrame);
  lLambda->SetLineWidth(2);
  lLambda->SetLineStyle(9);
  lLambda->SetLineColor(GetColor(0));
  TLine* lAntiLambda = new TLine(xMinFrame, MassLambda0, xMaxFrame, MassLambda0);
  lAntiLambda->SetLineWidth(2);
  lAntiLambda->SetLineStyle(9);
  lAntiLambda->SetLineColor(GetColor(0));

  frame->Draw();
  mass->Draw("same colz");
  lLambda->Draw("same");
  lAntiLambda->Draw("same");
  latexText = TString::Format("%s, #it{p}_{T, V0} = %.1f - %.1f GeV/c", dataSet.c_str(), lowpt, highpt).Data();
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(myCanvas->GetName());
}

// ----------------------------------------------------------

// Compare mass of V0s with R cut to mass of uncut V0s
// Can't do this yet, because I don't have a R vs M hist
void plotV0Rcheck(int train, string hadron, double ptmin, double ptmax)
{
  string inName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string dataSet = getDataSet(train);

  saveName = hadron;
  saveName += "_Rs";
  saveName += "_v0pt%.0f-%.0f", v0ptmin, v0ptmax;
  saveName = TString::Format("%s.pdf", saveName.c_str());
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation";
  if (train == 349871) histName += "_id24580";
  histName += "/data/V0/"
  histName += hadron;
  histName += "PtRadiusCosPA";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  array<int, 2> ptBins = getProjectionBins(th3->GetXaxis(), ptmin, ptmax);
  double lowpt = th3->GetXaxis()->GetBinLowEdge(ptBins[0]);
  double highpt = th3->GetXaxis()->GetBinUpEdge(ptBins[1]);
  double rMin = 0., rMax = 40.;
  array<int, 2> rBins = getProjectionBins(th3->GetYaxis(), rMin, rMax);

  TH1* mass = (TH1*)th3->ProjectionZ("mass", ptBins[0], ptBins[1], 1, th3->GetNbinsY());
  TH1* cut  = (TH1*)th3->ProjectionZ("cut", ptBins[0], ptBins[1], rBins[0], rBins[1]);
  setStyle(mass, 0);
  setStyle(cut, 1);
  mass->Rebin(4); cut->Rebin(4);

  double xMinFrame = mass->GetXaxis()->GetXmin(), xMaxFrame = mass->GetXaxis()->GetXmax();
  double yMinFrame = 0, yMaxFrame = 1.1 * getHistUpperBound(mass);
  // double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  xTitle = formatHadronDaughters(hadron).c_str();
  yTitle = "counts";
  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/c", hadron.c_str(), lowpt, highpt).Data();

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataset + ", " + ptText).c_str());
  frame->Draw();
  mass->Draw("same");
  cut->Draw("same");
  canvas->SaveAs(saveName);
}

// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================

void plotTrain(int train, int setting, double ptmin, double ptmin)
{
  string inName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string dataSet = getDataSet(train);
  vector<string> inputStrings = {inName, dataSet};

    switch(setting) {
    case 0:
      plotV0Radius(inName, dataSet, ptmin, ptmax);
      break;
    case 1:
      plotV0CosPA(inName, dataSet, ptmin, ptmax);
      break;
    case 2:
      plotV0DCAdaughters(inName, dataSet, ptmin, ptmax);
      break;
    case 3:
      plotV0DCApos(inName, dataSet, ptmin, ptmax);
      break;
    case 4:
      plotV0DCAneg(inName, dataSet, ptmin, ptmax);
      return;
    case 5:
      plotV0ctau(inName, dataSet, 1, ptmin, ptmax);
      plotV0ctau(inName, dataSet, 2, ptmin, ptmax);
      plotV0ctau(inName, dataSet, 3, ptmin, ptmax);
      return;
    case 6:
      plotV0massKL(inName, dataSet, ptmin, ptmax);
      return;
    case 7:
      plotV0massLL(inName, dataSet, ptmin, ptmax);
      return;
    default:
      cout << "Invalid setting!" << endl;
      return;
  }
}
void plotTrain(int train, int setting, string hadron, double ptmin, double ptmax)
{
  string inName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string dataSet = getDataSet(train);
  vector<string> inputStrings = {inName, dataSet, hadron};

  switch(setting) {
    case 0:
      plotV0Rcheck(inputStrings, ptmin, ptmax);
      return;
    default:
      cout << "Invalid setting!" << endl;
      return;
  }
}

void plot252064(int setting, double ptmin, double ptmax)
{ plotTrain(252064, setting, ptmin, ptmax); }
void plot252064(int setting, string hadron, double ptmin, double ptmax)
{ plotTrain(252064, setting, hadron, ptmin, ptmax); }

// void plot22o(double ptmin, double ptmax, int setting)
// {
//   string inName = "~/cernbox/TrainOutput/252064/AnalysisResults.root";
//   string dataSet = "LHC22o_pass6";

//   switch(setting) {
//     case 0:
//       plotV0Radius(inName, dataSet, ptmin, ptmax);
//       break;
//     case 1:
//       plotV0CosPA(inName, dataSet, ptmin, ptmax);
//       break;
//     case 2:
//       plotV0DCAdaughters(inName, dataSet, ptmin, ptmax);
//       break;
//     case 3:
//       plotV0DCApos(inName, dataSet, ptmin, ptmax);
//       break;
//     case 4:
//       plotV0DCAneg(inName, dataSet, ptmin, ptmax);
//       return;
//     case 5:
//       plotV0ctau(inName, dataSet, 1, ptmin, ptmax);
//       plotV0ctau(inName, dataSet, 2, ptmin, ptmax);
//       plotV0ctau(inName, dataSet, 3, ptmin, ptmax);
//       return;
//     case 6:
//       plotV0massKL(inName, dataSet, ptmin, ptmax);
//       return;
//     case 7:
//       plotV0massLL(inName, dataSet, ptmin, ptmax);
//       return;
//     default:
//       cout << "Invalid setting!" << endl;
//       return;
//   }
// }

// void plot22o(int setting)
// {
//   vector<double> pt = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0};
//   for (int i = 0; i < pt.size() - 1; i++) {
//     plot22o(pt[i], pt[i + 1], setting);
//   }
// }
