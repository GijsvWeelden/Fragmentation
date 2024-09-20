
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

void plotV0Radius(string inName, string dataSet)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  xTitle = "#it{p}_{T, V0} (GeV/#it{c})";
  yTitle = "V0 Radius";

  histName = "jet-fragmentation/data/V0/V0PtRadiusCosPA";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  TH2D* v0Radius = (TH2D*)th3->Project3D("yx");
  v0Radius->SetName("v0Radius");
  normaliseHistColByCol(v0Radius);
  v0Radius->SetMinimum(1e-5);
  v0Radius->SetMaximum(1.);

  string canvasName = "v0PtRadius";
  canvasName += ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);
  canvas->SetLogz();

  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0, yMaxFrame = 100.;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  vector<TH2D*> histVector = {v0Radius};
  plotNHists(canvas, frame, histVector, nullptr, {}, "colz");
}
void plotV0CosPA(string inName, string dataSet)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  xTitle = "#it{p}_{T, V0} (GeV/#it{c})";
  yTitle = "V0 cos(PA)";

  histName = "jet-fragmentation/data/V0/V0PtRadiusCosPA";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  TH2D* v0CosPA = (TH2D*)th3->Project3D("zx");
  v0CosPA->SetName("v0CosPA");
  normaliseHistColByCol(v0CosPA);
  v0CosPA->SetMinimum(1e-5);
  v0CosPA->SetMaximum(1.);

  string canvasName = "v0PtCosPA";
  canvasName += ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);
  canvas->SetLogz();

  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0.95, yMaxFrame = 1.0;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  vector<TH2D*> histVector = {v0CosPA};
  plotNHists(canvas, frame, histVector, nullptr, {}, "colz");
}
void plotV0DCAdaughters(string inName, string dataSet)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  xTitle = "#it{p}_{T, V0} (GeV/#it{c})";
  yTitle = "V0 DCA daughters";

  histName = "jet-fragmentation/data/V0/V0PtDCAd";
  TFile *inFile = TFile::Open(inName.c_str());
  TH2D* v0DCAd = (TH2D*)inFile->Get(histName.c_str());
  normaliseHistColByCol(v0DCAd);
  v0DCAd->SetMinimum(1e-5);
  v0DCAd->SetMaximum(1.);

  string canvasName = "v0PtDCAd";
  canvasName += ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);
  canvas->SetLogz();

  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0., yMaxFrame = 1.0;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  vector<TH2D*> histVector = {v0DCAd};
  plotNHists(canvas, frame, histVector, nullptr, {}, "colz");
}
void plotV0DCApos(string inName, string dataSet)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  xTitle = "#it{p}_{T, V0} (GeV/#it{c})";
  yTitle = "V0 DCA pos";

  histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  TH2D* v0DCAp = (TH2D*)th3->Project3D("yx");
  v0DCAp->SetName("v0DCAp");
  normaliseHistColByCol(v0DCAp);
  v0DCAp->SetMinimum(1e-5);
  v0DCAp->SetMaximum(1.);

  string canvasName = "v0PtDCAp";
  canvasName += ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);
  canvas->SetLogz();

  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = -10., yMaxFrame = 10.;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  vector<TH2D*> histVector = {v0DCAp};
  plotNHists(canvas, frame, histVector, nullptr, {}, "colz");
}
void plotV0DCAneg(string inName, string dataSet)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  xTitle = "#it{p}_{T, V0} (GeV/#it{c})";
  yTitle = "V0 DCA neg";

  histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  TH2D* v0DCAn = (TH2D*)th3->Project3D("zx");
  v0DCAn->SetName("v0DCAn");
  normaliseHistColByCol(v0DCAn);
  v0DCAn->SetMinimum(1e-5);
  v0DCAn->SetMaximum(1.);

  string canvasName = "v0PtDCAn";
  canvasName += ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);
  canvas->SetLogz();

  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = -10., yMaxFrame = 10.;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  vector<TH2D*> histVector = {v0DCAn};
  plotNHists(canvas, frame, histVector, nullptr, {}, "colz");
}
void plotV0ctau(string inName, string dataSet, int setting)
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
  string hadron;
  if (setting == K0SAxis) { hadron = "K0S"; }
  if (setting == LambdaAxis) { hadron = "Lambda"; }
  if (setting == AntiLambdaAxis) { hadron = "AntiLambda"; }

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  xTitle = "#it{p}_{T, V0} (GeV/#it{c})";
  yTitle = "V0 #it{c}#tau (" + formatHadronName(hadron) + ")";

  histName = "jet-fragmentation/data/V0/V0PtCtau";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  TH2D* v0ctau = (TH2D*)thn->Projection(setting, ptAxis);
  v0ctau->SetName("v0ctau");
  normaliseHistColByCol(v0ctau);
  v0ctau->SetMinimum(1e-5);
  v0ctau->SetMaximum(1.);

  string canvasName = "v0Ptctau" + hadron;
  canvasName += ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);
  canvas->SetLogz();

  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0., yMaxFrame = 40.;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  vector<TH2D*> histVector = {v0ctau};
  plotNHists(canvas, frame, histVector, nullptr, {}, "colz");
}
// void plotV0DCAdaughters(string inName = "AnalysisResults.root", string dataSet = "dataSet")
// {
//   double time = clock();
//   gStyle->SetNdivisions(505);

//   string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
//   double textSize = 0.04;
//   double labelSize = 0.04;
//   double titleSize = 0.04;

//   bool setLogY = false;
//   double xMinFrame = 0., xMaxFrame = 1., yMinFrame = 0., yMaxFrame = 0.6;
//   double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
//   int xCanvas = 900, yCanvas = 900;
//   int rebinNumber = 5;
//   xTitle = "DCA daughters";
//   yTitle = "normalised count";
//   latexText = dataSet;

//   std::vector<TH1D*> histVector;

//   TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
//   if (setLogY) { myCanvas->SetLogy(); }
//   TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
//   TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

//   histName = "jet-fragmentation/data/V0/V0PtDCAd";
//   TFile *inFile = TFile::Open(inName.c_str());
//   TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());

//   int firstBinPt = 1, lastBinPt = th2->GetNbinsX();
//   firstBinPt = th2->GetXaxis()->FindBin(v0ptmin + 1e-3);
//   lastBinPt  = th2->GetXaxis()->FindBin(v0ptmax - 1e-3);
//   th2->GetXaxis()->SetRange(firstBinPt, lastBinPt);
//   TH1D* v0DCAdaughters = (TH1D*)th2->ProjectionY("v0DCAdaughters", firstBinPt, lastBinPt);
//   v0DCAdaughters->Scale(1./v0DCAdaughters->Integral());
//   setStyle(v0DCAdaughters, 0);
//   histVector.push_back(v0DCAdaughters);

//   saveName = "v0DCAdaughters";
//   saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
//   saveName = TString::Format("%s.pdf", saveName.c_str());
//   plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
// }
// void plotV0DCApos(string inName = "AnalysisResults.root", string dataSet = "dataSet")
// {
//   double time = clock();
//   gStyle->SetNdivisions(505);

//   string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
//   double textSize = 0.04;
//   double labelSize = 0.04;
//   double titleSize = 0.04;

//   bool setLogY = true;
//   double xMinFrame = -10., xMaxFrame = 10., yMinFrame = 1e-5, yMaxFrame = 0.5;
//   double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
//   int xCanvas = 900, yCanvas = 900;
//   int rebinNumber = 5;
//   xTitle = "DCA pos";
//   yTitle = "normalised count";
//   latexText = dataSet;

//   std::vector<TH1D*> histVector;

//   TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
//   if (setLogY) { myCanvas->SetLogy(); }
//   TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
//   TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

//   histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
//   TFile *inFile = TFile::Open(inName.c_str());
//   TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

//   int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
//   int firstBinDCA = 1, lastBinDCA = th3->GetNbinsZ();
//   firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
//   lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
//   TH1D* v0DCApos = (TH1D*)th3->ProjectionY("v0DCApos", firstBinPt, lastBinPt, firstBinDCA, lastBinDCA);
//   v0DCApos->Scale(1./v0DCApos->Integral());
//   setStyle(v0DCApos, 0);
//   histVector.push_back(v0DCApos);

//   saveName = "v0DCApos";
//   saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
//   saveName = TString::Format("%s.pdf", saveName.c_str());
//   plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
// }
// void plotV0DCAneg(string inName = "AnalysisResults.root", string dataSet = "dataSet")
// {
//   double time = clock();
//   gStyle->SetNdivisions(505);

//   string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
//   double textSize = 0.04;
//   double labelSize = 0.04;
//   double titleSize = 0.04;

//   bool setLogY = true;
//   double xMinFrame = -10., xMaxFrame = 10., yMinFrame = 1e-5, yMaxFrame = 0.5;
//   double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
//   int xCanvas = 900, yCanvas = 900;
//   int rebinNumber = 5;
//   xTitle = "DCA neg";
//   yTitle = "normalised count";
//   latexText = dataSet;

//   std::vector<TH1D*> histVector;

//   TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
//   if (setLogY) { myCanvas->SetLogy(); }
//   TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
//   TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

//   histName = "jet-fragmentation/data/V0/V0PtDCAposneg";
//   TFile *inFile = TFile::Open(inName.c_str());
//   TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

//   int firstBinPt = 1, lastBinPt = th3->GetNbinsX();
//   int firstBinDCA = 1, lastBinDCA = th3->GetNbinsY();
//   firstBinPt = th3->GetXaxis()->FindBin(v0ptmin + 1e-3);
//   lastBinPt  = th3->GetXaxis()->FindBin(v0ptmax - 1e-3);
//   TH1D* v0DCAneg = (TH1D*)th3->ProjectionZ("v0DCAneg", firstBinPt, lastBinPt, firstBinDCA, lastBinDCA);
//   v0DCAneg->Scale(1./v0DCAneg->Integral());
//   setStyle(v0DCAneg, 0);
//   histVector.push_back(v0DCAneg);

//   saveName = "v0DCAneg";
//   saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
//   saveName = TString::Format("%s.pdf", saveName.c_str());
//   plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
// }
// void plotV0ctau(string inName = "AnalysisResults.root", string dataSet = "dataSet", int setting = 1)
// {
//   const int nDim           = 4;
//   const int ptAxis         = 0;
//   const int K0SAxis        = 1;
//   const int LambdaAxis     = 2;
//   const int AntiLambdaAxis = 3;

//   double time = clock();
//   gStyle->SetNdivisions(505);
//   if (setting < K0SAxis || setting > AntiLambdaAxis) {
//     cout << "Invalid setting! 1 for K0S, 2 for Lambda, 3 for antiLambda" << endl;
//     return;
//   }

//   string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
//   double textSize = 0.04;
//   double labelSize = 0.04;
//   double titleSize = 0.04;

//   bool setLogY = true;
//   double xMinFrame = 0., xMaxFrame = 40., yMinFrame = 1e-5, yMaxFrame = .1;
//   double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
//   int xCanvas = 900, yCanvas = 900;
//   xTitle = "#it{c}#tau (K^{0}_{S})";
//   if (setting == LambdaAxis) { xTitle = "#it{c}#tau (#Lambda)"; }
//   if (setting == AntiLambdaAxis) { xTitle = "#it{c}#tau (#bar{#Lambda})"; }
//   yTitle = "normalised count";
//   latexText = dataSet;

//   std::vector<TH1D*> histVector;


//   histName = "jet-fragmentation/data/V0/V0PtCtau";
//   TFile *inFile = TFile::Open(inName.c_str());
//   THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
//   thn->Sumw2();

//   int firstBinPt = 1, lastBinPt = thn->GetAxis(ptAxis)->GetNbins();
//   firstBinPt = thn->GetAxis(ptAxis)->FindBin(v0ptmin + 1e-3);
//   lastBinPt  = thn->GetAxis(ptAxis)->FindBin(v0ptmax - 1e-3);
//   thn->GetAxis(ptAxis)->SetRange(firstBinPt, lastBinPt);
//   TH1D* ctau = (TH1D*)thn->Projection(setting);
//   ctau->SetName(TString::Format("ctau_%d", setting).Data());
//   ctau->Scale(1./ctau->Integral());
//   setStyle(ctau, 0);
//   histVector.push_back(ctau);

//   saveName = "ctauK0S";
//   if (setting == LambdaAxis) { saveName = "ctauLambda"; }
//   else if (setting == AntiLambdaAxis) { saveName = "ctauAntiLambda"; }
//   saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
//   saveName = TString::Format("%s.pdf", saveName.c_str());

//   TCanvas* myCanvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);
//   if (setLogY) { myCanvas->SetLogy(); }
//   TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
//   TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
//   plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
// }
// void plotV0mass(string inName = "AnalysisResults.root", string dataSet = "dataSet", int setting = 1)
// {
//   const int nDim           = 4;
//   const int ptAxis         = 0;
//   const int K0SAxis        = 1;
//   const int LambdaAxis     = 2;
//   const int AntiLambdaAxis = 3;

//   double time = clock();
//   gStyle->SetNdivisions(505);
//   if (setting < K0SAxis || setting > AntiLambdaAxis) {
//     cout << "Invalid setting! 1 for K0S, 2 for Lambda, 3 for antiLambda" << endl;
//     return;
//   }

//   string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
//   double textSize = 0.04;
//   double labelSize = 0.04;
//   double titleSize = 0.04;

//   bool setLogY = true;
//   double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1e-4, yMaxFrame = .1;
//   if (setting == LambdaAxis || setting == AntiLambdaAxis) {
//     xMinFrame = 1.015, xMaxFrame = 1.215;
//   }
//   double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
//   int xCanvas = 900, yCanvas = 900;
//   xTitle = "#it{M} (K^{0}_{S})";
//   if (setting == LambdaAxis) { xTitle = "#it{M} (#Lambda)"; }
//   if (setting == AntiLambdaAxis) { xTitle = "#it{M} (#bar{#Lambda})"; }
//   yTitle = "normalised count";
//   latexText = dataSet;

//   std::vector<TH1D*> histVector;

//   TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
//   if (setLogY) { myCanvas->SetLogy(); }
//   TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
//   TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

//   histName = "jet-fragmentation/data/V0/V0PtMass";
//   TFile *inFile = TFile::Open(inName.c_str());
//   THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
//   thn->Sumw2();

//   int firstBinPt = 1, lastBinPt = thn->GetAxis(ptAxis)->GetNbins();
//   firstBinPt = thn->GetAxis(ptAxis)->FindBin(v0ptmin + 1e-3);
//   lastBinPt  = thn->GetAxis(ptAxis)->FindBin(v0ptmax - 1e-3);
//   thn->GetAxis(ptAxis)->SetRange(firstBinPt, lastBinPt);
//   TH1D* mass = (TH1D*)thn->Projection(setting);
//   mass->SetName(TString::Format("mass_%d", setting).Data());
//   mass->Scale(1./mass->Integral());
//   setStyle(mass, 0);
//   histVector.push_back(mass);

//   saveName = "massK0S";
//   if (setting == LambdaAxis) { saveName = "massLambda"; }
//   else if (setting == AntiLambdaAxis) { saveName = "massAntiLambda"; }
//   saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
//   saveName = TString::Format("%s.pdf", saveName.c_str());
//   plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
// }

void plot22o(int setting)
{
  string inName = "~/cernbox/TrainOutput/252064/AnalysisResults.root";
  string dataSet = "LHC22o_pass6";

  switch(setting) {
    case 0:
      plotV0Radius(inName, dataSet);
      break;
    case 1:
      plotV0CosPA(inName, dataSet);
      break;
    case 2:
      plotV0DCAdaughters(inName, dataSet);
      break;
    case 3:
      plotV0DCApos(inName, dataSet);
      break;
    case 4:
      plotV0DCAneg(inName, dataSet);
      break;
    case 5:
      plotV0ctau(inName, dataSet, 1);
      plotV0ctau(inName, dataSet, 2);
      plotV0ctau(inName, dataSet, 3);
      return;
    default:
      cout << "Invalid setting!" << endl;
      return;
  }
}
