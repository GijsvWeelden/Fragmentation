
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

// ----------------------------------------------------------

void fakePt(string inName = "AnalysisResults.root", string hadron = "", double jetptmin = 10., double jetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim      = 4;
  const int jetPtAxis = 0;
  const int v0PtAxis  = 1;
  const int v0EtaAxis = 2;
  const int v0PhiAxis = 3;

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0., yMaxFrame = 200.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  yTitle = "#it{p}_{T, jet} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("fakeJetPt%sPtEtaPhi", hadron.c_str()).Data();
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetPtAxis), jetptmin, jetptmax);
  thn->GetAxis(jetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetPtAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0pt = (TH2D*)thn->Projection(jetPtAxis, v0PtAxis);
  normaliseHistRowByRow(v0pt);
  v0pt->SetMinimum(1e-5);
  v0pt->SetMaximum(1.);
  v0pt->SetName(TString::Format("%sPt", hadron.c_str()).Data());
  histVector.push_back(v0pt);

  latexText = TString::Format("#splitline{ %s }{ Fake %s in jets }", dataSet.c_str(), formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("fake%sPt", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void fakeZ(string inName = "AnalysisResults.root", string hadron = "", double jetptmin = 10., double jetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 200., yMinFrame = 1e-3, yMaxFrame = 1. + 1e-3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet} (GeV/#it{c})";
  yTitle = TString::Format("#it{z}_{%s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("fakeJetPt%sTrackProj", hadron.c_str()).Data();
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(th2->GetXaxis(), jetptmin, jetptmax);
  double lowjetpt = th2->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th2->GetXaxis()->GetBinUpEdge(jetptbins[1]);

  // normaliseHistRowByRow(th2);
  normaliseHistColByCol(th2);
  th2->SetMinimum(1e-5);
  th2->SetMaximum(1.);
  th2->SetName(TString::Format("%sZ", hadron.c_str()).Data());
  histVector.push_back(th2);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ Fake %s in jets }}", dataSet.c_str(), lowjetpt, highjetpt, hadron.c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("fake%sZ", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void fakeCtau(string inName = "", string hadron = "", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1 ) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim                = 5;
  const int jetPtAxis           = 0;
  const int v0PtAxis            = 1;
  const int K0SctauAxis         = 2;
  const int Lambda0ctauAxis     = 3;
  const int antiLambda0ctauAxis = 4;
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0., yMaxFrame = 60.;
  if (doZ) { yMinFrame = 1e-3; yMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{c}#tau_{%s} (cm)", formatHadronName(hypothesis).c_str()).Data();
  yTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { yTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("fakeJetPt%sPtCtau", hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("fakeJetPt%sTrackProjCtau", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetPtAxis), jetptmin, jetptmax);
  thn->GetAxis(jetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetPtAxis)->GetBinUpEdge(jetptbins[1]);

  // Ugly but compact way of projecting the correct axis
  TH2D* v0ctau = (TH2D*)thn->Projection(v0PtAxis, K0SctauAxis * ("K0S" == hypothesis) + Lambda0ctauAxis * ("Lambda0" == hypothesis) + antiLambda0ctauAxis * ("AntiLambda0" == hypothesis));
  normaliseHistRowByRow(v0ctau);
  v0ctau->SetMinimum(1e-5);
  v0ctau->SetMaximum(1.);
  v0ctau->SetName(TString::Format("%sPtCtau%s", hadron.c_str(), hypothesis.c_str()).Data());
  if (doZ) { v0ctau->SetName(TString::Format("%sZCtau%s", hadron.c_str(), hypothesis.c_str()).Data()); }
  histVector.push_back(v0ctau);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/#it{c} }{ Fake %s in jets }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("fake%sPtCtau%s", hadron.c_str(), hypothesis.c_str()).Data();
  if (doZ) { saveName = TString::Format("fake%sZCtau%s", hadron.c_str(), hypothesis.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void fakeMass(string inName = "", string hadron = "", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1 ) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim                = 5;
  const int jetPtAxis           = 0;
  const int v0PtAxis            = 1;
  const int K0SmassAxis         = 2;
  const int Lambda0massAxis     = 3;
  const int antiLambda0massAxis = 4;
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 60.;
  if ("K0S" == hypothesis) { xMinFrame = 0.4; xMaxFrame = 0.6; }
  if (doZ) { yMinFrame = 1e-3; yMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronName(hypothesis).c_str()).Data();
  yTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { yTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("fakeJetPt%sPtMass", hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("fakeJetPt%sTrackProjMass", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetPtAxis), jetptmin, jetptmax);
  thn->GetAxis(jetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetPtAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0mass = (TH2D*)thn->Projection(v0PtAxis, K0SmassAxis * ("K0S" == hypothesis) + Lambda0massAxis * ("Lambda0" == hypothesis) + antiLambda0massAxis * ("AntiLambda0" == hypothesis));
  normaliseHistRowByRow(v0mass);
  v0mass->SetMinimum(1e-5);
  v0mass->SetMaximum(1.);
  v0mass->SetName(TString::Format("%sPtMass%s", hadron.c_str(), hypothesis.c_str()).Data());
  if (doZ) { v0mass->SetName(TString::Format("%sTrackProjMass%s", hadron.c_str(), hypothesis.c_str()).Data()); }
  histVector.push_back(v0mass);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/#it{c} }{ Fake %s in jets }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("fake%sPtMass%s", hadron.c_str(), hypothesis.c_str()).Data();
  if (doZ) { saveName = TString::Format("fake%sZMass%s", hadron.c_str(), hypothesis.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void fakeRadius(string inName = "", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim       = 4;
  const int jetPtAxis  = 0;
  const int v0PtAxis   = 1;
  const int radiusAxis = 2;
  const int cosPAAxis  = 3;

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 100., yMinFrame = 0., yMaxFrame = 60.;
  if (doZ) { yMinFrame = 1e-3; yMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{R} (cm)";
  yTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { yTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("fakeJetPt%sPtRadiusCosPA", hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("fakeJetPt%sTrackProjRadiusCosPA", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetPtAxis), jetptmin, jetptmax);
  thn->GetAxis(jetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetPtAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0radius = (TH2D*)thn->Projection(v0PtAxis, radiusAxis);
  normaliseHistRowByRow(v0radius);
  v0radius->SetMinimum(1e-5);
  v0radius->SetMaximum(1.);
  v0radius->SetName(TString::Format("%sPtRadius", hadron.c_str()).Data());
  if (doZ) { v0radius->SetName(TString::Format("%sTrackProjRadius", hadron.c_str()).Data()); }
  histVector.push_back(v0radius);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/#it{c} }{ Fake %s in jets }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("fake%sPtRadius", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("fake%sZRadius", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void fakeCosPA(string inName = "", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim       = 4;
  const int jetPtAxis  = 0;
  const int v0PtAxis   = 1;
  const int radiusAxis = 2;
  const int cosPAAxis  = 3;

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0.95, xMaxFrame = 1., yMinFrame = 0., yMaxFrame = 60.;
  if (doZ) { yMinFrame = 1e-3; yMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "cos(PA)";
  yTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { yTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("fakeJetPt%sPtRadiusCosPA", hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("fakeJetPt%sTrackProjRadiusCosPA", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetPtAxis), jetptmin, jetptmax);
  thn->GetAxis(jetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetPtAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0cospa = (TH2D*)thn->Projection(v0PtAxis, cosPAAxis);
  normaliseHistRowByRow(v0cospa);
  v0cospa->SetMinimum(1e-5);
  v0cospa->SetMaximum(1.);
  v0cospa->SetName(TString::Format("%sPtCosPA", hadron.c_str()).Data());
  if (doZ) { v0cospa->SetName(TString::Format("%sZCosPA", hadron.c_str()).Data()); }
  histVector.push_back(v0cospa);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/#it{c} }{ Fake %s in jets }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("fake%sPtCosPA", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("fake%sZCosPA", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void fakeDCAposneg(string inName = "", string hadron = "", bool doNeg = false, double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim       = 4;
  const int jetPtAxis  = 0;
  const int v0PtAxis   = 1;
  const int DCAPosAxis = 2;
  const int DCANegAxis = 3;

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = -10., xMaxFrame = 10., yMinFrame = 0., yMaxFrame = 60.;
  if (doZ) { yMinFrame = 1e-3; yMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = doNeg ? "DCA neg. (cm)" : "DCA pos. (cm)";
  yTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { yTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("fakeJetPt%sPtDCAposneg", hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("fakeJetPt%sTrackProjDCAposneg", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetPtAxis), jetptmin, jetptmax);
  thn->GetAxis(jetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetPtAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0dca = (TH2D*)thn->Projection(v0PtAxis, doNeg ? DCANegAxis : DCAPosAxis);
  normaliseHistRowByRow(v0dca);
  v0dca->SetMinimum(1e-5);
  v0dca->SetMaximum(1.);
  v0dca->SetName(TString::Format("%sPtDCA", hadron.c_str()).Data());
  if (doZ) { v0dca->SetName(TString::Format("%sZDCA", hadron.c_str()).Data()); }
  histVector.push_back(v0dca);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/#it{c} }{ Fake %s in jets }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("fake%sPtDCA", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("fake%sZDCA", hadron.c_str()).Data(); }
  saveName = TString::Format(doNeg ? "%sneg" : "%spos", saveName.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void fakeDCAd(string inName = "", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 10., yMinFrame = 0., yMaxFrame = 60.;
  if (doZ) { yMinFrame = 1e-3; yMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "DCA daughters (cm^{2})";
  yTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { yTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("fakeJetPt%sPtDCAd", hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("fakeJetPt%sTrackProjDCAd", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(th3->GetXaxis(), jetptmin, jetptmax);
  th3->GetXaxis()->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = th3->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th3->GetXaxis()->GetBinUpEdge(jetptbins[1]);

  TH2D* v0dca = (TH2D*)th3->Project3D("yz");
  normaliseHistRowByRow(v0dca);
  v0dca->SetMinimum(1e-5);
  v0dca->SetMaximum(1.);
  v0dca->SetName(TString::Format("%sPtDCAd", hadron.c_str()).Data());
  if (doZ) { v0dca->SetName(TString::Format("%sZDCAd", hadron.c_str()).Data()); }
  histVector.push_back(v0dca);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/#it{c} }{ Fake %s in jets }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("fake%sPtDCAd", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("fake%sZDCAd", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

// ----------------------------------------------------------

void missPt(string inName = "AnalysisResults.root", string hadron = "", double jetptmin = 10., double jetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim      = 4;
  const int jetPtAxis = 0;
  const int v0PtAxis  = 1;
  const int v0EtaAxis = 2;
  const int v0PhiAxis = 3;

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 0., yMaxFrame = 200.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  yTitle = "#it{p}_{T, jet} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("missJetPt%sPtEtaPhi", hadron.c_str()).Data();
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetPtAxis), jetptmin, jetptmax);
  thn->GetAxis(jetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetPtAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0pt = (TH2D*)thn->Projection(jetPtAxis, v0PtAxis);
  normaliseHistRowByRow(v0pt);
  v0pt->SetMinimum(1e-5);
  v0pt->SetMaximum(1.);
  v0pt->SetName(TString::Format("%sPt", hadron.c_str()).Data());
  histVector.push_back(v0pt);

  latexText = TString::Format("#splitline{ %s }{ Missed %s in jets }", dataSet.c_str(), formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("miss%sPt", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void missZ(string inName = "AnalysisResults.root", string hadron = "", double jetptmin = 10., double jetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("V0" == hadron) + ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 200., yMinFrame = 1e-3, yMaxFrame = 1. + 1e-3;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet} (GeV/#it{c})";
  yTitle = TString::Format("#it{z}_{%s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("missJetPt%sTrackProj", hadron.c_str()).Data();
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(th2->GetXaxis(), jetptmin, jetptmax);
  double lowjetpt = th2->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th2->GetXaxis()->GetBinUpEdge(jetptbins[1]);

  // normaliseHistRowByRow(th2);
  normaliseHistColByCol(th2);
  th2->SetMinimum(1e-5);
  th2->SetMaximum(1.);
  th2->SetName(TString::Format("%sZ", hadron.c_str()).Data());
  histVector.push_back(th2);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ Missed %s in jets }}", dataSet.c_str(), lowjetpt, highjetpt, hadron.c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("fake%sZ", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

// ----------------------------------------------------------
// Lambda0
// ----------------------------------------------------------

void fakeLambda0Pt(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ fakePt(inName, "Lambda0", jetptmin, jetptmax); }

void fakeLambda0Z(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ fakeZ(inName, "Lambda0", jetptmin, jetptmax); }

void fakeLambda0Ctau(string inName = "AnalysisResults.root", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeCtau(inName, "Lambda0", hypothesis, jetptmin, jetptmax, doZ); }

void fakeLambda0Mass(string inName = "AnalysisResults.root", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeMass(inName, "Lambda0", hypothesis, jetptmin, jetptmax, doZ); }

void fakeLambda0Radius(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeRadius(inName, "Lambda0", jetptmin, jetptmax, doZ); }

void fakeLambda0CosPA(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeCosPA(inName, "Lambda0", jetptmin, jetptmax, doZ); }

void fakeLambda0DCAposneg(string inName = "AnalysisResults.root", bool doNeg = false, double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeDCAposneg(inName, "Lambda0", doNeg, jetptmin, jetptmax, doZ); }

void fakeLambda0DCAd(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeDCAd(inName, "Lambda0", jetptmin, jetptmax, doZ); }

void missLambda0Pt(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ missPt(inName, "Lambda0", jetptmin, jetptmax); }

void missLambda0Z(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ missZ(inName, "Lambda0", jetptmin, jetptmax); }

// ----------------------------------------------------------
// AntiLambda0
// ----------------------------------------------------------

void fakeAntiLambda0Pt(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ fakePt(inName, "AntiLambda0", jetptmin, jetptmax); }

void fakeAntiLambda0Z(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ fakeZ(inName, "AntiLambda0", jetptmin, jetptmax); }

void fakeAntiLambda0Ctau(string inName = "AnalysisResults.root", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeCtau(inName, "AntiLambda0", hypothesis, jetptmin, jetptmax, doZ); }

void fakeAntiLambda0Mass(string inName = "AnalysisResults.root", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeMass(inName, "AntiLambda0", hypothesis, jetptmin, jetptmax, doZ); }

void fakeAntiLambda0Radius(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeRadius(inName, "AntiLambda0", jetptmin, jetptmax, doZ); }

void fakeAntiLambda0CosPA(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeCosPA(inName, "AntiLambda0", jetptmin, jetptmax, doZ); }

void fakeAntiLambda0DCAposneg(string inName = "AnalysisResults.root", bool doNeg = false, double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeDCAposneg(inName, "AntiLambda0", doNeg, jetptmin, jetptmax, doZ); }

void fakeAntiLambda0DCAd(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeDCAd(inName, "AntiLambda0", jetptmin, jetptmax, doZ); }

void missAntiLambda0Pt(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ missPt(inName, "AntiLambda0", jetptmin, jetptmax); }

void missAntiLambda0Z(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ missZ(inName, "AntiLambda0", jetptmin, jetptmax); }

// ----------------------------------------------------------
// K0S
// ----------------------------------------------------------

void fakeK0SPt(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ fakePt(inName, "K0S", jetptmin, jetptmax); }

void fakeK0SZ(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ fakeZ(inName, "K0S", jetptmin, jetptmax); }

void fakeK0SCtau(string inName = "AnalysisResults.root", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeCtau(inName, "K0S", hypothesis, jetptmin, jetptmax, doZ); }

void fakeK0SMass(string inName = "AnalysisResults.root", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeMass(inName, "K0S", hypothesis, jetptmin, jetptmax, doZ); }

void fakeK0SRadius(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeRadius(inName, "K0S", jetptmin, jetptmax, doZ); }

void fakeK0SCosPA(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeCosPA(inName, "K0S", jetptmin, jetptmax, doZ); }

void fakeK0SDCAposneg(string inName = "AnalysisResults.root", bool doNeg = false, double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeDCAposneg(inName, "K0S", doNeg, jetptmin, jetptmax, doZ); }

void fakeK0SDCAd(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{ fakeDCAd(inName, "K0S", jetptmin, jetptmax, doZ); }

void missK0SPt(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ missPt(inName, "K0S", jetptmin, jetptmax); }

void missK0SZ(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{ missZ(inName, "K0S", jetptmin, jetptmax); }
