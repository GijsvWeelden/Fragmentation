
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

double getNjets(TFile* inFile, double jetptmin, double jetptmax)
{
  string histName = "jet-fragmentation/matching/jets/matchDetJetPtEtaPhi"; // Should this be particle level?
  TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
  std::array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
  return jets->Integral(jetptbins[0], jetptbins[1], 1, jets->GetNbinsY(), 1, jets->GetNbinsZ());
}

void matchedJetPt(string inName = "AnalysisResults.root")
{
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

// ----------------------------------------------------------

void matchedPt(string inName = "AnalysisResults.root", string hadron = "", double partjetptmin = 10., double partjetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim          = 4;
  const int partJetPtAxis = 0;
  const int partV0PtAxis  = 1;
  const int detJetPtAxis  = 2;
  const int detV0PtAxis   = 3;

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
  xTitle = TString::Format("#it{p}_{T, %s}^{det.} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  yTitle = TString::Format("#it{p}_{T, %s}^{part.} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("partJetPt%sPtDetJetPt%sPt", hadron.c_str(), hadron.c_str()).Data();
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
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
  v0pt->SetName(TString::Format("%sPt", hadron.c_str()).Data());
  histVector.push_back(v0pt);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched %s in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("matched%sPt", hadron.c_str()).Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void matchedZ(string inName = "AnalysisResults.root", string hadron = "", double partjetptmin = 10., double partjetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim          = 4;
  const int detJetPtAxis  = 0;
  const int detV0ZAxis    = 1;
  const int partJetPtAxis = 2;
  const int partV0ZAxis   = 3;

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 1e-3, xMaxFrame = 1. + 1e-3, yMinFrame = 1e-3, yMaxFrame = 1. + 1e-3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{z}_{%s}^{det.} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  yTitle = TString::Format("#it{z}_{%s}^{part.} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("matchDetJetPt%sTrackProjPartJetPt%sTrackProj", hadron.c_str(), hadron.c_str()).Data();
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
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
  v0z->SetName(TString::Format("%sZ", hadron.c_str()).Data());
  histVector.push_back(v0z);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/c }{ Matched V0s in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("matched%sZ", hadron.c_str()).Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void matchedCtau(string inName = "", string hadron = "", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1 ) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim          = 5;
  const int partJetPtAxis = 0;
  const int partV0PtAxis  = 1;
  const int detJetPtAxis  = 2;
  const int detV0PtAxis   = 3;
  const int detV0CtauAxis = 4;

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
  xTitle = TString::Format("#it{c}#tau_{%s} (cm)", formatHadronName(hypothesis).c_str()).Data();
  yTitle = TString::Format("#it{p}_{T, %s}^{part.} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("partJetPt%sPtDetJetPt%sPtCtau%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data();
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0ctau = (TH2D*)thn->Projection(partV0PtAxis, detV0CtauAxis);
  normaliseHistRowByRow(v0ctau);
  v0ctau->SetMinimum(1e-5);
  v0ctau->SetMaximum(1.);
  v0ctau->SetName(TString::Format("%sPtCtau%s", hadron.c_str(), hypothesis.c_str()).Data());
  histVector.push_back(v0ctau);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ Matched %s in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("matched%sPtCtau%s", hadron.c_str(), hypothesis.c_str()).Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void matchedMass(string inName = "", string hadron = "", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1 ) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim          = 5;
  const int partJetPtAxis = 0;
  const int partV0PtAxis  = 1;
  const int detJetPtAxis  = 2;
  const int detV0PtAxis   = 3;
  const int detV0MassAxis = 4;

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 60.;
  if ("K0S" == hypothesis) { xMinFrame = 0.4; xMaxFrame = 0.6; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronName(hypothesis).c_str()).Data();
  yTitle = TString::Format("#it{p}_{T, %s}^{part.} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("partJetPt%sPtDetJetPt%sPtMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data();
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH2D* v0mass = (TH2D*)thn->Projection(partV0PtAxis, detV0MassAxis);
  normaliseHistRowByRow(v0mass);
  v0mass->SetMinimum(1e-5);
  v0mass->SetMaximum(1.);
  v0mass->SetName(TString::Format("%sPtMass%s", hadron.c_str(), hypothesis.c_str()).Data());
  histVector.push_back(v0mass);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ Matched %s in matched jets }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("matched%sPtMass%s", hadron.c_str(), hypothesis.c_str()).Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

// ----------------------------------------------------------
// Lambda0
// ----------------------------------------------------------

void matchedLambda0Pt(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedPt(inName, "Lambda0", partjetptmin, partjetptmax); }

void matchedLambda0Z(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedZ(inName, "Lambda0", partjetptmin, partjetptmax); }

void matchedLambda0Ctau(string inName = "AnalysisResults.root", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedCtau(inName, "Lambda0", hypothesis, partjetptmin, partjetptmax); }

void matchedLambda0Mass(string inName = "AnalysisResults.root", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedMass(inName, "Lambda0", hypothesis, partjetptmin, partjetptmax); }

// ----------------------------------------------------------
// AntiLambda0
// ----------------------------------------------------------

void matchedAntiLambda0Pt(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedPt(inName, "AntiLambda0", partjetptmin, partjetptmax); }

void matchedAntiLambda0Z(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedZ(inName, "AntiLambda0", partjetptmin, partjetptmax); }

void matchedAntiLambda0Ctau(string inName = "AnalysisResults.root", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedCtau(inName, "AntiLambda0", hypothesis, partjetptmin, partjetptmax); }

void matchedAntiLambda0Mass(string inName = "AnalysisResults.root", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedMass(inName, "AntiLambda0", hypothesis, partjetptmin, partjetptmax); }

// ----------------------------------------------------------
// K0S
// ----------------------------------------------------------

void matchedK0SPt(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedPt(inName, "K0S", partjetptmin, partjetptmax); }

void matchedK0SZ(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedZ(inName, "K0S", partjetptmin, partjetptmax); }

void matchedK0SCtau(string inName = "AnalysisResults.root", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedCtau(inName, "K0S", hypothesis, partjetptmin, partjetptmax); }

void matchedK0SMass(string inName = "AnalysisResults.root", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200.)
{ matchedMass(inName, "K0S", hypothesis, partjetptmin, partjetptmax); }
