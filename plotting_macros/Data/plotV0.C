
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

void plotV0Pt(string inName = "AnalysisResults.root")
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 100, yMinFrame = 1e-6, yMaxFrame = .1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, V0}";
  yTitle = "normalised count";
  latexText = "LHC23y_pass1";

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtEtaPhi";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());
  TH1D* v0pt = (TH1D*)th3->ProjectionX("v0pt");
  v0pt->Scale(1./v0pt->Integral());
  setStyle(v0pt, 0);
  histVector.push_back(v0pt);

  saveName = "v0pt";
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0Eta(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtEtaPhi";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0Phi(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtEtaPhi";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0EtaPhi(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, V0} = %.0f - %.0f GeV/c", v0ptmin, v0ptmax).Data();

  std::vector<TH2D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtEtaPhi";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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

void plotV0Radius(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtRadiusCosPA";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0CosPA(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtRadiusCosPA";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0DCAdaughters(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtDCAd";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0DCApos(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0DCAneg(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0DCAposneg(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, V0} = %.0f - %.0f GeV/c", v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0ctau(string inName = "AnalysisResults.root", int setting = 1, double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtCtau";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
// This plot isn't very useful. ctauL = 2*ctauK0S?
void plotV0ctauKL(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtCtau";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0mass(string inName = "AnalysisResults.root", int setting = 1, double v0ptmin = 0., double v0ptmax = 100.)
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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, V0} = %.0f - %.0f GeV/c}", v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtMass";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
void plotV0massKL(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1.015, yMaxFrame = 1.215;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{M} (K^{0}_{S})";
  yTitle = "#it{M} (#Lambda)";
  latexText = TString::Format("LHC23y_pass1, #it{p}_{T, V0} = %.0f - %.0f GeV/c", v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtMass";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptAxis)->FindBin(v0ptmax - 1e-3);
  thn->GetAxis(ptAxis)->SetRange(firstBinPt, lastBinPt);
  TH2D* mass = (TH2D*)thn->Projection(LambdaAxis, K0SAxis);
  mass->SetName("mass_KL");
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);

  saveName = "massK0SLambda";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0massLL(string inName = "AnalysisResults.root", double v0ptmin = 0., double v0ptmax = 100.)
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
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 1.015, yMaxFrame = 1.215;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{M} (#Lambda)";
  yTitle = "#it{M} (#bar{#Lambda})";
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, V0} = %.0f - %.0f GeV/c", v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/V0/V0PtMass";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptAxis)->FindBin(v0ptmax - 1e-3);
  thn->GetAxis(ptAxis)->SetRange(firstBinPt, lastBinPt);
  TH2D* mass = (TH2D*)thn->Projection(AntiLambdaAxis, LambdaAxis);
  mass->SetName("mass_LL");
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);

  saveName = "massLambdaAntiLambda";
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
