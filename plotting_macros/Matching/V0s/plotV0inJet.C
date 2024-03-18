
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

// ----------------------------------------------------------

void matchedJetPt(string inName = "AnalysisResults.root")
{
  double time = clock();
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 200., yMinFrame = 0., yMaxFrame = 200.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.1, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, V0}^{det.} (GeV/#it{c})";
  yTitle = "#it{p}_{T, V0}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/matchDetJetPtPartJetPt";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH2D* matchedjetpt = (TH2D*)inFile->Get(histName.c_str());
  normaliseHistRowByRow(matchedjetpt);
  matchedjetpt->SetMinimum(1e-5);
  matchedjetpt->SetMaximum(1.);
  matchedjetpt->SetName("matchedjetpt");
  histVector.push_back(matchedjetpt);

  latexText = TString::Format("%s, geometrically matched jets", dataSet.c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedJetPt";
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void matchedV0Pt(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{
  const int nDim          = 4;
  const int partJetPtAxis = 0;
  const int partV0PtAxis  = 1;
  const int detJetPtAxis  = 2;
  const int detV0PtAxis   = 3;

  double time = clock();
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0., yMaxFrame = 60.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, V0}^{det.} (GeV/#it{c})";
  yTitle = "#it{p}_{T, V0}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/V0/partJetPtV0PtDetJetPtV0Pt";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0pt = (TH2D*)thn->Projection(partV0PtAxis, detV0PtAxis);
  normaliseHistRowByRow(v0pt);
  v0pt->SetMinimum(1e-5);
  v0pt->SetMaximum(1.);
  v0pt->SetName("v0pt");
  histVector.push_back(v0pt);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched V0s in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedV0Pt";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

void matchedV0Z(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{
  const int nDim          = 4;
  const int detJetPtAxis  = 0;
  const int detV0ZAxis    = 1;
  const int partJetPtAxis = 2;
  const int partV0ZAxis   = 3;

  double time = clock();
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = false;
  double xMinFrame = 1e-3, xMaxFrame = 1. + 1e-3, yMinFrame = 1e-3, yMaxFrame = 1. + 1e-3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}_{V0}^{det.} (GeV/#it{c})";
  yTitle = "#it{z}_{V0}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/V0/matchDetJetPtV0TrackProjPartJetPtV0TrackProj";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0z = (TH2D*)thn->Projection(partV0ZAxis, detV0ZAxis);
  normaliseHistRowByRow(v0z);
  v0z->SetMinimum(1e-5);
  v0z->SetMaximum(1.);
  v0z->SetName("v0z");
  histVector.push_back(v0z);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched V0s in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedV0Z";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

// ----------------------------------------------------------
// Lambda0
// ----------------------------------------------------------

void matchedLambda0Pt(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{
  const int nDim          = 4;
  const int partJetPtAxis = 0;
  const int partV0PtAxis  = 1;
  const int detJetPtAxis  = 2;
  const int detV0PtAxis   = 3;

  double time = clock();
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = false;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0., yMaxFrame = 60.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, #Lambda^{0}}^{det.} (GeV/#it{c})";
  yTitle = "#it{p}_{T, #Lambda^{0}}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/V0/partJetPtLambda0PtDetJetPtLambda0Pt";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0pt = (TH2D*)thn->Projection(partV0PtAxis, detV0PtAxis);
  normaliseHistRowByRow(v0pt);
  v0pt->SetName("Lambda0Pt");
  histVector.push_back(v0pt);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched #Lambda^{0} in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedLambda0Pt";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

void matchedLambda0Z(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{
  const int nDim          = 4;
  const int detJetPtAxis  = 0;
  const int detV0ZAxis    = 1;
  const int partJetPtAxis = 2;
  const int partV0ZAxis   = 3;

  double time = clock();
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = false;
  double xMinFrame = 1e-3, xMaxFrame = 1. + 1e-3, yMinFrame = 1e-3, yMaxFrame = 1. + 1e-3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}_{#Lambda^{0}}^{det.} (GeV/#it{c})";
  yTitle = "#it{z}_{#Lambda^{0}}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/V0/matchDetJetPtLambda0TrackProjPartJetPtLambda0TrackProj";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0z = (TH2D*)thn->Projection(partV0ZAxis, detV0ZAxis);
  normaliseHistRowByRow(v0z);
  v0z->SetName("Lambda0z");
  histVector.push_back(v0z);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched #Lambda^{0} in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedLambda0Z";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

// ----------------------------------------------------------
// AntiLambda0
// ----------------------------------------------------------

void matchedAntiLambda0Pt(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{
  const int nDim          = 4;
  const int partJetPtAxis = 0;
  const int partV0PtAxis  = 1;
  const int detJetPtAxis  = 2;
  const int detV0PtAxis   = 3;

  double time = clock();
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = false;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0., yMaxFrame = 60.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, #bar{Lambda}^{0}}^{det.} (GeV/#it{c})";
  yTitle = "#it{p}_{T, #bar{Lambda}^{0}}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/V0/partJetPtAntiLambda0PtDetJetPtAntiLambda0Pt";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0pt = (TH2D*)thn->Projection(partV0PtAxis, detV0PtAxis);
  normaliseHistRowByRow(v0pt);
  v0pt->SetName("AntiLambda0Pt");
  histVector.push_back(v0pt);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched #bar{Lambda}^{0} in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedAntiLambda0Pt";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

void matchedAntiLambda0Z(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{
  const int nDim          = 4;
  const int detJetPtAxis  = 0;
  const int detV0ZAxis    = 1;
  const int partJetPtAxis = 2;
  const int partV0ZAxis   = 3;

  double time = clock();
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = false;
  double xMinFrame = 1e-3, xMaxFrame = 1. + 1e-3, yMinFrame = 1e-3, yMaxFrame = 1. + 1e-3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}_{#bar{#Lambda}^{0}}^{det.} (GeV/#it{c})";
  yTitle = "#it{z}_{#bar{#Lambda}^{0}}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/V0/matchDetJetPtAntiLambda0TrackProjPartJetPtAntiLambda0TrackProj";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0z = (TH2D*)thn->Projection(partV0ZAxis, detV0ZAxis);
  normaliseHistRowByRow(v0z);
  v0z->SetName("AntiLambda0z");
  histVector.push_back(v0z);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched #bar{#Lambda}^{0} in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedAntiLambda0Z";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

// ----------------------------------------------------------
// K0S
// ----------------------------------------------------------

void matchedK0SPt(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{
  const int nDim          = 4;
  const int partJetPtAxis = 0;
  const int partV0PtAxis  = 1;
  const int detJetPtAxis  = 2;
  const int detV0PtAxis   = 3;

  double time = clock();
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = false;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0., yMaxFrame = 60.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, K^{0}_{S}}^{det.} (GeV/#it{c})";
  yTitle = "#it{p}_{T, K^{0}_{S}}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/V0/partJetPtK0SPtDetJetPtK0SPt";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0pt = (TH2D*)thn->Projection(partV0PtAxis, detV0PtAxis);
  normaliseHistRowByRow(v0pt);
  v0pt->SetName("K0SPt");
  histVector.push_back(v0pt);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched K^{0}_{S} in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedK0SPt";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

void matchedK0SZ(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{
  const int nDim          = 4;
  const int detJetPtAxis  = 0;
  const int detV0ZAxis    = 1;
  const int partJetPtAxis = 2;
  const int partV0ZAxis   = 3;

  double time = clock();
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = false;
  double xMinFrame = 1e-3, xMaxFrame = 1. + 1e-3, yMinFrame = 1e-3, yMaxFrame = 1. + 1e-3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}_{K^{0}_{S}}^{det.} (GeV/#it{c})";
  yTitle = "#it{z}_{K^{0}_{S}}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/V0/matchDetJetPtK0STrackProjPartJetPtK0STrackProj";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0z = (TH2D*)thn->Projection(partV0ZAxis, detV0ZAxis);
  normaliseHistRowByRow(v0z);
  v0z->SetName("K0Sz");
  histVector.push_back(v0z);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched K^{0}_{S} in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedK0SZ";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}