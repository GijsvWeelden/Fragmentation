
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

string formatHadronName(string hadron);
template <typename T>
T loadHist(string fileName, string histName);
void normaliseHistRowByRow(TH2D* hist);
void normaliseHistColByCol(TH2D* hist);
template <typename T>
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<T> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);
template <typename T>
void setStyle(T hist, int styleNumber);

// ----------------------------------------------------------

void plotRadius(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = 0., xMaxFrame = 1., yMinFrame = 0., yMaxFrame = 100.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}";
  yTitle = "R";
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, jet} = %.0f - %.0f GeV/c", jetptmin, jetptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtAntiLambdaTrackProjRadius";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = th3->GetNbinsX();
  firstBinJetPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinJetPt, lastBinJetPt);

  TH2D* zRadius = (TH2D*)th3->Project3D("zy");
  zRadius->SetName(TString::Format("zRadius_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  zRadius->Scale(1./zRadius->Integral());
  setStyle(zRadius, 0);
  zRadius->SetMinimum(1e-7);
  zRadius->SetMaximum(1.);

  saveName = "zRadius";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  zRadius->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.95, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotCosPA(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = 0., xMaxFrame = 1., yMinFrame = 0.95, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}";
  yTitle = "cos(PA)";
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, jet} = %.0f - %.0f GeV/c", jetptmin, jetptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtAntiLambdaTrackProjCosPA";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = th3->GetNbinsX();
  firstBinJetPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinJetPt, lastBinJetPt);

  TH2D* zCosPA = (TH2D*)th3->Project3D("zy");
  zCosPA->SetName(TString::Format("zCosPA_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  zCosPA->Scale(1./zCosPA->Integral());
  setStyle(zCosPA, 0);
  zCosPA->SetMinimum(1e-7);
  zCosPA->SetMaximum(1.);

  saveName = "zCosPA";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  zCosPA->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.95, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotDCAdaughters(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = 0., xMaxFrame = 1., yMinFrame = 0., yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}";
  yTitle = "DCA daughters";
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, jet} = %.0f - %.0f GeV/c", jetptmin, jetptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtAntiLambdaTrackProjDCAd";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = th3->GetNbinsX();
  firstBinJetPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinJetPt, lastBinJetPt);

  TH2D* zDCAdaughters = (TH2D*)th3->Project3D("zy");
  zDCAdaughters->SetName(TString::Format("zDCAdaughters_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  zDCAdaughters->Scale(1./zDCAdaughters->Integral());
  setStyle(zDCAdaughters, 0);
  zDCAdaughters->SetMinimum(1e-7);
  zDCAdaughters->SetMaximum(1.);

  saveName = "zDCAdaughters";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  zDCAdaughters->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.95, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotDCApos(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0zAxis    = 1;
  const int dcaposAxis = 2;
  const int dcanegAxis = 3;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = -10., yMaxFrame = 10.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}";
  yTitle = "DCA pos";
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, jet} = %.0f - %.0f GeV/c", jetptmin, jetptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtAntiLambdaTrackProjDCAposneg";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  TH2D* zDCApos = (TH2D*)thn->Projection(dcaposAxis, v0zAxis);
  zDCApos->SetName(TString::Format("zDCApos_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  zDCApos->Scale(1./zDCApos->Integral());
  setStyle(zDCApos, 0);

  saveName = "zDCApos";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  zDCApos->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.95, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotDCAneg(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0zAxis    = 1;
  const int dcaposAxis = 2;
  const int dcanegAxis = 3;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = -10., yMaxFrame = 10.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}";
  yTitle = "DCA neg";
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, jet} = %.0f - %.0f GeV/c", jetptmin, jetptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtAntiLambdaTrackProjDCAposneg";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  TH2D* zDCAneg = (TH2D*)thn->Projection(dcaposAxis, v0zAxis);
  zDCAneg->SetName(TString::Format("zDCAneg_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  zDCAneg->Scale(1./zDCAneg->Integral());
  setStyle(zDCAneg, 0);

  saveName = "zDCAneg";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  zDCAneg->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.95, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotctau(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = 0., xMaxFrame = 1., yMinFrame = 0., yMaxFrame = 40.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}";
  yTitle = "#it{c}#tau (#bar{#Lambda})";
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, jet} = %.0f - %.0f GeV/c", jetptmin, jetptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtAntiLambdaTrackProjCtau";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = th3->GetNbinsX();
  firstBinJetPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinJetPt, lastBinJetPt);

  TH2D* zCtau = (TH2D*)th3->Project3D("zy");
  zCtau->SetName(TString::Format("zCtau_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  zCtau->Scale(1./zCtau->Integral());
  setStyle(zCtau, 0);
  zCtau->SetMinimum(1e-7);
  zCtau->SetMaximum(1.);

  saveName = "zCtau";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  zCtau->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.95, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotMass(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = 0., xMaxFrame = 1., yMinFrame = 1.015, yMaxFrame = 1.215;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}";
  yTitle = "#it{M} (#bar{#Lambda})";
  latexText = TString::Format("LHC23y_pass1_small, #it{p}_{T, jet} = %.0f - %.0f GeV/c", jetptmin, jetptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtAntiLambdaTrackProjMass";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = th3->GetNbinsX();
  firstBinJetPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
  th3->GetXaxis()->SetRange(firstBinJetPt, lastBinJetPt);

  TH2D* zMass = (TH2D*)th3->Project3D("zy");
  zMass->SetName(TString::Format("zMass_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  zMass->Scale(1./zMass->Integral());
  setStyle(zMass, 0);
  zMass->SetMinimum(1e-7);
  zMass->SetMaximum(1.);

  saveName = "zMass";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  zMass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.95, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

// ----------------------------------------------------------
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
  else if (hadron == "Lambda0"){
    had = "#it{#Lambda}^{0}";
  }
  return had;
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
  if (latexText != "") { DrawLatex(0.4, 0.8, latexText.c_str(), legend->GetTextSize()); }
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
