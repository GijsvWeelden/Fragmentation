
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
  xTitle = "#it{p}_{T, jet}";
  yTitle = "normalised count";
  latexText = "LHC23y_pass1";

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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, jet} = %.0f - %.0f GeV/c}", jetptmin, jetptmax).Data();

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
  latexText = TString::Format("#splitline{LHC23y_pass1_small}{#it{p}_{T, jet} = %.0f - %.0f GeV/c}", jetptmin, jetptmax).Data();

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
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, jet} = %.0f - %.0f GeV/c", jetptmin, jetptmax).Data();

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
