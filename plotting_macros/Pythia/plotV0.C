
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
template <typename T>
T loadHist(string fileName, string histName);
void normaliseHistRowByRow(TH2D* hist);
void normaliseHistColByCol(TH2D* hist);
template <typename T>
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<T> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);
template <typename T>
void setStyle(T hist, int styleNumber);

// ----------------------------------------------------------

void plotNdaughtersV0in(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;

  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 20, yMinFrame = 1e-5, yMaxFrame = 2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, V0}";
  yTitle = "Fraction of V0s";
  latexText = TString::Format("#it{p}_{T, jet} = %.0f - %.0f, V0 #in jet", ptmin, ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hNdauV0in";
  THnSparseD* thn = loadHist<THnSparseD*>(inName, histName);

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptjetAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  thn->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  TH2D* v0ndau = (TH2D*)thn->Projection(obsAxis, ptv0Axis);
  v0ndau->SetName("hNdauV0in");
  v0ndau->Rebin2D(rebinNumber, 1);
  normaliseHistColByCol(v0ndau);

  // Fraction of jets that do not contain n daughters
  for (int i = 0; i <= 2; i++) {
    string name = TString::Format("h%ddau", i).Data();
    string legendEntry = TString::Format("%d daughter(s) in jet", i).Data();
    TH1D* hist = (TH1D*)v0ndau->ProjectionX(name.c_str(), i+1, i+1);
    hist->SetTitle("V0s in jets");
    setStyle(hist, i);
    histVector.push_back(hist);
    legend->AddEntry(hist, legendEntry.c_str());
  }
  saveName = "ndaughters-v0in";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotNdaughtersV0out(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;

  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 20, yMinFrame = 1e-3, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, V0}";
  yTitle = "Fraction of V0s";
  latexText = TString::Format("#it{p}_{T, jet} = %.0f - %.0f, V0 #notin jet", ptmin, ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hNdauV0out";
  THnSparseD* thn = loadHist<THnSparseD*>(inName, histName);

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptjetAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  thn->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  TH2D* v0ndau = (TH2D*)thn->Projection(obsAxis, ptv0Axis);
  v0ndau->SetName("hNdauV0out");
  v0ndau->Rebin2D(rebinNumber, 1);
  normaliseHistColByCol(v0ndau);

  // Fraction of jets that do not contain n daughters
  for (int i = 0; i <= 2; i++) {
    string name = TString::Format("h%ddau", i).Data();
    string legendEntry = TString::Format("%d daughter(s) in jet", i).Data();
    TH1D* hist = (TH1D*)v0ndau->ProjectionX(name.c_str(), i+1, i+1);
    hist->SetTitle("V0s not in jets");
    setStyle(hist, i);
    histVector.push_back(hist);
    legend->AddEntry(hist, legendEntry.c_str());
  }
  saveName = "ndaughters-v0out";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotNdaughtersV0in2D(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  gStyle->SetPalette(77);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;

  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 20, yMinFrame = 0, yMaxFrame = 100;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 2000, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, V0}";
  yTitle = "#it{p}_{T, jet}";

  histName = "hNdauV0in";
  saveName = "ndaughters-v0in2D";
  THnSparseD* thn = loadHist<THnSparseD*>(inName, histName);
  TCanvas* canvas = new TCanvas("c", "c", xCanvas, yCanvas);
  canvas->Divide(3, 1);
  // TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);

  for (int i = 0; i <= 2; i++) {
    canvas->cd(i+1);
    // DrawFrame in loop because it uses gPad, which can confuse the pads,
    //   leading to us only having one plot in the end
    TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
    thn->GetAxis(obsAxis)->SetRange(i+1, i+1);
    TH2D* hist = (TH2D*)thn->Projection(ptjetAxis, ptv0Axis);
    hist->SetName(TString::Format("h%ddau", i).Data());
    hist->SetStats(0);
    hist->SetTitle("");
    normaliseHistRowByRow(hist);
    hist->SetMinimum(1e-5);
    hist->SetMaximum(1);
    auto* pad = canvas->cd(i+1);
    pad->SetLogz();
    frame->Draw();
    hist->Draw("colz same");
    latexText = TString::Format("V0 #in jet, %d dau #in jet, norm. per #it{p}_{T, jet}", i).Data();
    DrawLatex(0.3, 0.93, latexText.c_str(), textSize);
  }
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->cd(0);
  canvas->SaveAs(saveName.c_str());
}
void plotDeltaPtV0in(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "#it{p}_{T, jet}^{orig.} - #it{p}_{T, jet}^{corrected}";
  yTitle = "";
  latexText = TString::Format("#it{p}_{T, jet} = %.0f - %.0f, V0 #in jet", ptmin, ptmax).Data();

  bool setLogY = false;
  double xMinFrame = -50, xMaxFrame = 50, yMinFrame = 0, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;

  std::vector<TH1D*> histVector;
  std::vector<double> v0CutVector = {0., 1., 3., 5., 10., 20.};

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hDeltaPtV0in";
  THnSparse* thn = loadHist<THnSparse*>(inName, histName);

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptjetAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  thn->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  TH2D* v0dpt = (TH2D*)thn->Projection(obsAxis, ptv0Axis);
  v0dpt->SetName("v0dpt");
  v0dpt->Rebin2D(rebinNumber, 1);
  normaliseHistColByCol(v0dpt);
  // v0dpt->Draw("colz");

  // Fraction of jets that do not contain n daughters
  for (int i = 0; i < v0CutVector.size() - 1; i++) {
    double v0min = v0CutVector[i];
    double v0max = v0CutVector[i+1];
    double binv0min = v0dpt->GetXaxis()->FindBin(v0min + 1e-3);
    double binv0max = v0dpt->GetXaxis()->FindBin(v0max - 1e-3);
    string name = TString::Format("dpt_v0%.0f_%.0f", v0min, v0max).Data();
    string legendEntry = TString::Format("#it{p}_{T, V0} = %.0f - %.0f", v0min, v0max).Data();
    TH1D* hist = (TH1D*)v0dpt->ProjectionY(name.c_str(), i+1, i+1);
    setStyle(hist, i);
    histVector.push_back(hist);
    legend->AddEntry(hist, legendEntry.c_str());
  }
  saveName = "dptin";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotDeltaPtV0out(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "#it{p}_{T, jet}^{orig.} - #it{p}_{T, jet}^{corrected}";
  yTitle = "";
  latexText = TString::Format("#it{p}_{T, jet} = %.0f - %.0f, V0 #notin jet", ptmin, ptmax).Data();

  bool setLogY = false;
  double xMinFrame = -50, xMaxFrame = 50, yMinFrame = 0, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;

  std::vector<TH1D*> histVector;
  std::vector<double> v0CutVector = {0., 1., 3., 5., 10., 20.};

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hDeltaPtV0out";
  THnSparseD* thn = loadHist<THnSparseD*>(inName, histName);

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptjetAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  thn->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  TH2D* v0dpt = (TH2D*)thn->Projection(obsAxis, ptv0Axis);
  v0dpt->SetName("v0dpt");
  v0dpt->Rebin2D(rebinNumber, 1);
  normaliseHistColByCol(v0dpt);
  // v0dpt->Draw("colz");

  // Fraction of jets that do not contain n daughters
  for (int i = 0; i < v0CutVector.size() - 1; i++) {
    double v0min = v0CutVector[i];
    double v0max = v0CutVector[i+1];
    double binv0min = v0dpt->GetXaxis()->FindBin(v0min + 1e-3);
    double binv0max = v0dpt->GetXaxis()->FindBin(v0max - 1e-3);
    string name = TString::Format("dpt_v0%.0f_%.0f", v0min, v0max).Data();
    string legendEntry = TString::Format("#it{p}_{T, V0} = %.0f - %.0f", v0min, v0max).Data();
    TH1D* hist = (TH1D*)v0dpt->ProjectionY(name.c_str(), i+1, i+1);
    setStyle(hist, i);
    histVector.push_back(hist);
    legend->AddEntry(hist, legendEntry.c_str());
  }
  saveName = "dptout";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotAvgNV0in(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "#it{p}_{T, jet}";
  yTitle = "Average #it{N}(V0) per jet";
  latexText = TString::Format("#it{p}_{T, jet} = %.0f - %.0f", ptmin, ptmax).Data();

  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 0, yMaxFrame = 3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;

  std::vector<TH1D*> histVector;
  std::vector<double> v0CutVector = {0., 1., 3., 5., 10., 20.};

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hNv0";
  TH3D* th3 = loadHist<TH3D*>(inName, histName);
  TH2D* jetptnv0 = (TH2D*)th3->Project3D("zx");
  jetptnv0->SetName("jetptnv0");
  // jetptnv0->Rebin2D(rebinNumber, 1);
  TH1D* nv0 = (TH1D*)th3->ProjectionX("nv0");
  // TH1D* nv0 = (TH1D*)jetptnv0->ProjectionX("nv0");
  nv0->Reset();
  nv0->SetLineWidth(3);
  // Get average number of V0s per jet pt bin
  for (int iCol = 1; iCol <= jetptnv0->GetNbinsX(); iCol++) {
    double njets = 0;
    double avg = 0;
    for (int iRow = 1; iRow <= jetptnv0->GetNbinsY(); iRow++) {
      double binContent = jetptnv0->GetBinContent(iCol, iRow);
      double binValue = jetptnv0->GetYaxis()->GetBinCenter(iRow);
      njets += binContent;
      avg += binContent * binValue;
    }
    if (njets > 0) { avg /= njets; }
    nv0->Fill(iCol, avg);
  }
  histVector.push_back(nv0);

  saveName = "avgNv0injet";
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "hist", latexText);
}
void plotNV0in(string inName = "withdecays_pthat20_80.root")
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "#it{N}(V0) #in jet";
  yTitle = "Fraction";

  bool setLogY = false;
  double xMinFrame = -0.5, xMaxFrame = 20.5, yMinFrame = 0, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;

  std::vector<TH1D*> histVector;
  std::vector<double> jetCutVector = {10., 20., 40., 60., 80., 100.};
  // std::vector<double> v0CutVector = {0., 1., 3., 5., 10., 20.};

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hNv0";
  TH3D* th3 = loadHist<TH3D*>(inName, histName);

  for (int i = 0; i < jetCutVector.size() - 1; i++) {
    double jetptmin = jetCutVector[i];
    double jetptmax = jetCutVector[i+1];
    int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
    int firstBinRatio = 1, lastBinRatio = th3->GetNbinsY();
    firstBinPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
    lastBinPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
    TH1D* nv0 = (TH1D*)th3->ProjectionZ("nv0", firstBinPt, lastBinPt, firstBinRatio, lastBinRatio);
    nv0->SetName(TString::Format("nv0_%d", i).Data());
    nv0->Scale(1./nv0->Integral());
    setStyle(nv0, i);
    legend->AddEntry(nv0, TString::Format("#it{p}_{T, jet} = %.0f - %.0f", jetptmin, jetptmax).Data());
    histVector.push_back(nv0);
  }

  saveName = "nv0injet";
  // saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  // plotNHists(myCanvas, frame, histVector, legend, saveName, "hist", latexText);
  // We want to draw a line and the markers, so we draw "manually" here
  myCanvas->cd();
  frame->Draw();
  for (unsigned int i = 0; i < histVector.size(); i++) {
    TH1D* hist = histVector[i];
    hist->Draw("hist same");
    hist->Draw("hist p same");
  }
  if (legend) { legend->Draw("same"); }
  if (latexText != "") { DrawLatex(0.3, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotdRV0in(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "d#it{R} (V0, jet)";
  yTitle = "normalised count";
  latexText = TString::Format("#it{p}_{T, jet} = %.0f - %.0f", ptmin, ptmax).Data();

  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 0.5;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;

  std::vector<TH1D*> histVector;
  std::vector<double> v0CutVector = {0., 1., 3., 5., 10., 20.};

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hdRV0in";
  THnSparseD* thn = loadHist<THnSparseD*>(inName, histName);

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptjetAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  thn->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  TH2D* v0dR = (TH2D*)thn->Projection(obsAxis, ptv0Axis);
  v0dR->SetName("v0dR");
  v0dR->Rebin2D(rebinNumber, 5);

  // Fraction of jets that do not contain n daughters
  for (int i = 0; i < v0CutVector.size() - 1; i++) {
    double v0min = v0CutVector[i];
    double v0max = v0CutVector[i+1];
    double binv0min = v0dR->GetXaxis()->FindBin(v0min + 1e-3);
    double binv0max = v0dR->GetXaxis()->FindBin(v0max - 1e-3);
    string name = TString::Format("dR_v0%.0f_%.0f", v0min, v0max).Data();
    string legendEntry = TString::Format("#it{p}_{T, V0} = %.0f - %.0f", v0min, v0max).Data();
    TH1D* hist = (TH1D*)v0dR->ProjectionY(name.c_str(), i+1, i+1);
    hist->Scale(1./hist->Integral());
    setStyle(hist, i);
    histVector.push_back(hist);
    legend->AddEntry(hist, legendEntry.c_str());
  }
  saveName = "dRv0in";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotdRV0out(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "d#it{R} (V0, jet)";
  yTitle = "normalised count";
  latexText = TString::Format("#it{p}_{T, jet} = %.0f - %.0f", ptmin, ptmax).Data();

  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;

  std::vector<TH1D*> histVector;
  std::vector<double> v0CutVector = {0., 1., 3., 5., 10., 20.};

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hdRV0out";
  THnSparseD* thn = loadHist<THnSparseD*>(inName, histName);

  int firstBinPt = 1, lastBinPt = thn->GetAxis(ptjetAxis)->GetNbins();
  firstBinPt = thn->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt  = thn->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  thn->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  TH2D* v0dR = (TH2D*)thn->Projection(obsAxis, ptv0Axis);
  v0dR->SetName("v0dR");
  v0dR->Rebin2D(rebinNumber, 5);

  // Fraction of jets that do not contain n daughters
  for (int i = 0; i < v0CutVector.size() - 1; i++) {
    double v0min = v0CutVector[i];
    double v0max = v0CutVector[i+1];
    double binv0min = v0dR->GetXaxis()->FindBin(v0min + 1e-3);
    double binv0max = v0dR->GetXaxis()->FindBin(v0max - 1e-3);
    string name = TString::Format("dR_v0%.0f_%.0f", v0min, v0max).Data();
    string legendEntry = TString::Format("#it{p}_{T, V0} = %.0f - %.0f", v0min, v0max).Data();
    TH1D* hist = (TH1D*)v0dR->ProjectionY(name.c_str(), i+1, i+1);
    hist->Scale(1./hist->Integral());
    setStyle(hist, i);
    histVector.push_back(hist);
    legend->AddEntry(hist, legendEntry.c_str());
  }
  saveName = "dRv0out";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotdRdauout(string inName = "withdecays_pthat20_80.root", double jetmin = 10., double jetmax = 100., double v0min = 0., double v0max = 20.)
{
  double time = clock();
  gStyle->SetNdivisions(505, "xy");
  gStyle->SetPalette(77);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "d#it{R}(jet, daughter 1)";
  yTitle = "d#it{R}(jet, daughter 2)";
  latexText  = TString::Format("#it{p}_{T, jet} = %.0f - %.0f, #it{p}_{T, V0} = %.0f - %.0f", jetmin, jetmax, v0min, v0max).Data();
  string latexText2 = "V0 #in jet, dau #notin jet";

  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  std::vector<TH2D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  myCanvas->SetLogz();
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  TLine* lhori = new TLine(0., 0.4, 1., 0.4);
  lhori->SetLineWidth(3);
  lhori->SetLineColor(GetColor(0));
  lhori->SetLineStyle(9);
  TLine* lvert = new TLine(0.4, 0., 0.4, 1.);
  lvert->SetLineWidth(3);
  lvert->SetLineColor(GetColor(0));
  lvert->SetLineStyle(9);

  histName = "hMissDaudR";
  THnSparseD* thn = loadHist<THnSparseD*>(inName, histName);
  int _ptjetAxis = 0, _ptv0Axis = 1, _dRv0axis = 2, _dRd0axis = 3, _dRd1axis = 4; // Special axes

  int firstBinPt = 1, lastBinPt = thn->GetAxis(_ptjetAxis)->GetNbins();
  int firstBinV0 = 1, lastBinV0 = thn->GetAxis(_ptv0Axis)->GetNbins();
  firstBinPt = thn->GetAxis(_ptjetAxis)->FindBin(jetmin + 1e-3);
  lastBinPt  = thn->GetAxis(_ptjetAxis)->FindBin(jetmax - 1e-3);
  firstBinV0 = thn->GetAxis(_ptv0Axis)->FindBin(v0min + 1e-3);
  lastBinV0  = thn->GetAxis(_ptv0Axis)->FindBin(v0max - 1e-3);
  thn->GetAxis(_ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  thn->GetAxis(_ptv0Axis)->SetRange(firstBinV0, lastBinV0);
  TH2D* daudr = (TH2D*)thn->Projection(_dRd0axis, _dRd1axis);
  daudr->SetName("daudr");
  // daudr->Rebin2D(rebinNumber, 5);
  daudr->Scale(1./daudr->Integral());
  daudr->SetMinimum(1e-5);
  daudr->SetMaximum(1);
  histVector.push_back(daudr);
  saveName = "dRmissdau";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetmin, jetmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0min, v0max);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  // plotNHists(myCanvas, frame, histVector, legend, saveName, "colz", latexText);

  frame->Draw();
  daudr->Draw("same, colz");
  lhori->Draw("same");
  lvert->Draw("same");
  if (latexText != "") { DrawLatex(0.3, 0.93, latexText.c_str(), textSize); }
  if (latexText2 != "") { DrawLatex(0.5, 0.83, latexText2.c_str(), textSize); }
  myCanvas->SaveAs(saveName.c_str());
  // myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotZ(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100., double v0min = 0., double v0max = 20.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "#it{z}_{V0}";
  yTitle = "#frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{z}}";
  latexText = TString::Format("#it{p}_{T, jet} = %.0f - %.0f, #it{p}_{T, V0} = %.0f - %.0f", ptmin, ptmax, v0min, v0max).Data();

  bool setLogY = true;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 1e-4, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;

  std::vector<TH1D*> histVector;
  // std::vector<double> v0CutVector = {0., 1., 3., 5., 10., 20.};

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hzU";
  THnSparseD* uncorrected = loadHist<THnSparseD*>(inName, histName);
  histName = "hzC";
  THnSparseD* corrected = loadHist<THnSparseD*>(inName, histName);

  double njets = getNjets(inName, ptmin, ptmax);

  int firstBinPt = 1, lastBinPt = uncorrected->GetAxis(ptjetAxis)->GetNbins();
  int firstBinV0 = 1, lastBinV0 = uncorrected->GetAxis(ptv0Axis)->GetNbins();

  firstBinPt = uncorrected->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt  = uncorrected->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  firstBinV0 = uncorrected->GetAxis(ptv0Axis)->FindBin(v0min + 1e-3);
  lastBinV0  = uncorrected->GetAxis(ptv0Axis)->FindBin(v0max - 1e-3);

  uncorrected->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  uncorrected->GetAxis(ptv0Axis)->SetRange(firstBinV0, lastBinV0);
  TH1D* zU = (TH1D*)uncorrected->Projection(obsAxis);
  zU->SetName("zU");
  zU->Rebin(rebinNumber);
  zU->Scale(1./njets);
  setStyle(zU, 0);
  legend->AddEntry(zU, "#it{z}_{V0} (uncorrected #it{p}_{T, jet})");
  histVector.push_back(zU);

  corrected->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  corrected->GetAxis(ptv0Axis)->SetRange(firstBinV0, lastBinV0);
  TH1D* zC = (TH1D*)corrected->Projection(obsAxis);
  zC->SetName("zC");
  zC->Rebin(rebinNumber);
  zC->Scale(1./njets);
  setStyle(zC, 1);
  legend->AddEntry(zC, "#it{z}_{V0} (corrected #it{p}_{T, jet})");
  histVector.push_back(zC);

  saveName = "z";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotZRatio(string inName = "withdecays_pthat20_80.root", double ptmin = 10., double ptmax = 100., double v0min = 0., double v0max = 20.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  xTitle = "#it{z}_{V0}";
  yTitle = "corrected / uncorrected";
  latexText = TString::Format("#it{p}_{T, jet} = %.0f - %.0f, #it{p}_{T, V0} = %.0f - %.0f", ptmin, ptmax, v0min, v0max).Data();

  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;

  std::vector<TH1D*> histVector;
  std::vector<double> v0CutVector = {0., 1., 3., 5., 10., 20.};

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "hzU";
  THnSparseD* uncorrected = loadHist<THnSparseD*>(inName, histName);
  histName = "hzC";
  THnSparseD* corrected = loadHist<THnSparseD*>(inName, histName);

  // double njets = getNjets(inName, ptmin, ptmax);

  int firstBinPt = 1, lastBinPt = uncorrected->GetAxis(ptjetAxis)->GetNbins();
  int firstBinV0 = 1, lastBinV0 = uncorrected->GetAxis(ptv0Axis)->GetNbins();

  firstBinPt = uncorrected->GetAxis(ptjetAxis)->FindBin(ptmin + 1e-3);
  lastBinPt  = uncorrected->GetAxis(ptjetAxis)->FindBin(ptmax - 1e-3);
  firstBinV0 = uncorrected->GetAxis(ptv0Axis)->FindBin(v0min + 1e-3);
  lastBinV0  = uncorrected->GetAxis(ptv0Axis)->FindBin(v0max - 1e-3);

  uncorrected->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  uncorrected->GetAxis(ptv0Axis)->SetRange(firstBinV0, lastBinV0);
  TH1D* zU = (TH1D*)uncorrected->Projection(obsAxis);
  zU->SetName("zU");
  zU->Rebin(rebinNumber);
  // zU->Scale(1./njets);

  corrected->GetAxis(ptjetAxis)->SetRange(firstBinPt, lastBinPt);
  corrected->GetAxis(ptv0Axis)->SetRange(firstBinV0, lastBinV0);
  TH1D* zC = (TH1D*)corrected->Projection(obsAxis);
  zC->SetName("zC");
  zC->Rebin(rebinNumber);
  // zC->Scale(1./njets);

  TH1D* zR = (TH1D*)zC->Clone("zR");
  zR->Divide(zU);
  setStyle(zR, 0);
  histVector.push_back(zR);

  saveName = "zRatio";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  frame->Draw();
  zR->Draw("same");
  if (latexText != "") { DrawLatex(0.3, 0.93, latexText.c_str(), textSize); }
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
double getNjets(string inName, double jetptmin, double jetptmax)
{
  TH1D* hjetpt = loadHist<TH1D*>(inName, "hjetpt");
  int firstBinPt = 1, lastBinPt = hjetpt->GetNbinsX();
  firstBinPt = hjetpt->FindBin(jetptmin + 1e-3);
  lastBinPt  = hjetpt->FindBin(jetptmax - 1e-3);
  double integral = hjetpt->Integral(firstBinPt, lastBinPt);
  return integral;
}
template <typename T>
T loadHist(string fileName, string histName)
{
  TFile *inFile = TFile::Open(TString::Format("./%s", fileName.c_str()).Data());
  T hist = (T)inFile->Get(histName.c_str());
  return hist;
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
