
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

double getNjets(TFile* inFile, double jetptmin, double jetptmax)
{
  string histName = "jet-fragmentation/data/jets/jetPtEtaPhi";
  TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
  std::array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
  return jets->Integral(jetptbins[0], jetptbins[1], 1, jets->GetNbinsY(), 1, jets->GetNbinsZ());
}

// ----------------------------------------------------------

void plotV0Radius(string inName = "", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  string hadron = "V0";
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int radiusAxis = 2;
  const int cospaAxis  = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 0., yMaxFrame = 100.;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.45, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("%s Radius (cm)", formatHadronName(hadron).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  // For data, histograms do not contain "0" in the name
  histName = TString::Format("jetPt%sPtRadiusCosPA", hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjRadiusCosPA", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0Radius = (TH2D*)thn->Projection(radiusAxis, v0ptAxis);
  v0Radius->SetName(TString::Format("v0Radius_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  // v0Radius->Rebin2D(doZ ? 1 : 5, rebinNumber);
  normaliseHistColByCol(v0Radius);
  v0Radius->SetMinimum(1e-5);
  v0Radius->SetMaximum(1.);
  histVector.push_back(v0Radius);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c",
                              dataSet.c_str(), lowjetpt, highjetpt).Data();
  // latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtRadius", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZRadius", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotV0CosPA(string inName = "", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  string hadron = "V0";
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int radiusAxis = 2;
  const int cospaAxis  = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 0.95, yMaxFrame = 1.;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.2, yLatex = 0.95;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("%s cos PA", formatHadronName(hadron).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("jetPt%sPtRadiusCosPA", hadron.c_str()).Data();
  if(doZ) { histName = TString::Format("jetPt%sTrackProjRadiusCosPA", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0CosPA = (TH2D*)thn->Projection(cospaAxis, v0ptAxis);
  v0CosPA->SetName(TString::Format("v0CosPA_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  // // v0CosPA->Rebin2D(doZ ? 1 : 5, rebinNumber);
  normaliseHistColByCol(v0CosPA);
  v0CosPA->SetMinimum(1e-5);
  v0CosPA->SetMaximum(1.);
  histVector.push_back(v0CosPA);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c", dataSet.c_str(), lowjetpt, highjetpt).Data();
  // latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtCosPA", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZCosPA", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotV0ctau(string inName = "", string dataSet = "dataSet", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1 ) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  string hadron = "V0";
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 0., yMaxFrame = 40.;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.45, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("#it{c}#tau (%s)", formatHadronName(hypothesis).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jetPtV0PtCtau";
  if (doZ) { histName = "jetPtV0TrackProjCtau"; }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* ctau = (TH2D*)thn->Projection(K0SAxis * ("K0S" == hypothesis) + LambdaAxis * ("Lambda0" == hypothesis) + AntiLambdaAxis * ("AntiLambda0" == hypothesis), v0ptAxis);
  // ctau->Rebin2D(doZ ? 1 : 5, rebinNumber);
  normaliseHistColByCol(ctau);
  ctau->SetName(TString::Format("ctau_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  histVector.push_back(ctau);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c", dataSet.c_str(), lowjetpt, highjetpt).Data();
  // TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  TLatex* latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("V0Ptctau%s", hypothesis.c_str());
  if (doZ) { saveName = TString::Format("V0Zctau%s", hypothesis.c_str()); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotV0Mass(string inName = "", string dataSet = "dataSet", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1 ) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  string hadron = "V0";
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 1.05, yMaxFrame = 1.215;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  if ("K0S" == hypothesis) { yMinFrame = 0.4, yMaxFrame = 0.6; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.45, yLatex = 0.25;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("#it{M} (%s)", formatHadronName(hypothesis).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = TString::Format("jetPt%sPtMass", hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjMass", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* mass = (TH2D*)thn->Projection(K0SAxis * ("K0S" == hypothesis) + LambdaAxis * ("Lambda0" == hypothesis) + AntiLambdaAxis * ("AntiLambda0" == hypothesis), v0ptAxis);
  mass->SetName("v0mass");
  // mass->Rebin2D(doZ ? 1 : 5, rebinNumber);
  normaliseHistColByCol(mass);
  mass->SetMinimum(1e-5);
  mass->SetMaximum(1.);
  histVector.push_back(mass);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c",
                              dataSet.c_str(), jetptmin, jetptmax).Data();
  // TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  TLatex* latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtmass%s", hadron.c_str(), hypothesis.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZmass%s", hadron.c_str(), hypothesis.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotV0Masses(string inName = "", string dataSet = "dataSet", std::array<string, 2> hypotheses = {"", ""}, double jetptmin = 10., double jetptmax = 200., double v0min = 0., double v0max = 100., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if (hypotheses[0] == hypotheses[1]) {
    cout << "Error: hypotheses must be different" << endl;
    return;
  }
  if ( ("K0S" == hypotheses[0]) + ("Lambda0" == hypotheses[0]) + ("AntiLambda0" == hypotheses[0]) != 1 ) {
    cout << "Error: hypothesis 0 must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if ( ("K0S" == hypotheses[1]) + ("Lambda0" == hypotheses[1]) + ("AntiLambda0" == hypotheses[1]) != 1 ) {
    cout << "Error: hypothesis 1 must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  string hadron = "V0";
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 1.05, xMaxFrame = 1.215, yMinFrame = 1.05, yMaxFrame = 1.215;
  if ("K0S" == hypotheses[0]) { xMinFrame = 0.4, xMaxFrame = 0.6; }
  if ("K0S" == hypotheses[1]) { yMinFrame = 0.4, yMaxFrame = 0.6; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.45, yLatex = 0.25;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M} (%s)", formatHadronName(hypotheses[0]).c_str()).Data();
  yTitle = TString::Format("#it{M} (%s)", formatHadronName(hypotheses[1]).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = TString::Format("jetPt%sPtMass", hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjMass", hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0ptAxis), v0min, v0max);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(v0ptAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0 = thn->GetAxis(v0ptAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0ptAxis)->GetBinUpEdge(v0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; }

  int axisX = K0SAxis * ("K0S" == hypotheses[0]) + LambdaAxis * ("Lambda0" == hypotheses[0]) + AntiLambdaAxis * ("AntiLambda0" == hypotheses[0]);
  int axisY = K0SAxis * ("K0S" == hypotheses[1]) + LambdaAxis * ("Lambda0" == hypotheses[1]) + AntiLambdaAxis * ("AntiLambda0" == hypotheses[1]);
  TH2D* mass = (TH2D*)thn->Projection(axisY, axisX);
  mass->SetName("v0mass");
  mass->Scale(1./mass->Integral());
  mass->SetMinimum(1e-5);
  mass->SetMaximum(1.);
  histVector.push_back(mass);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/#it{c} }}",
                              dataSet.c_str(), jetptmin, jetptmax, hadron.c_str(), lowv0, highv0).Data();
  if (doZ) {
    latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c }{ #it{z}_{%s} = %.0f - %.0f }}",
                                dataSet.c_str(), jetptmin, jetptmax, hadron.c_str(), lowv0, highv0).Data();
  }
  // TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  TLatex* latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%smass%s-%s", hadron.c_str(), hypotheses[0].c_str(), hypotheses[1].c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  if (doZ) { saveName = TString::Format("%s_v0z%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  else { saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

// ----------------------------------------------------------

void plotRadius(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( "V0" == hadron ) {
    cout << "V0 selected. Redirecting to plotV0Radius" << endl;
    plotV0Radius(inName, dataSet, jetptmin, jetptmax, doZ);
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 0., yMaxFrame = 50.;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.45, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("%s Radius (cm)", formatHadronName(hadron).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  // For data, histograms do not contain "0" in the name
  string hadronForHist = hadron;
  if ("Lambda0" == hadronForHist) { hadronForHist = "Lambda"; }
  if ("AntiLambda0" == hadronForHist) { hadronForHist = "AntiLambda"; }
  histName = TString::Format("jetPt%sPtRadius", hadronForHist.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjRadius", hadronForHist.c_str()).Data(); }

  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(th3->GetXaxis(), jetptmin, jetptmax);
  th3->GetXaxis()->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = th3->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th3->GetXaxis()->GetBinUpEdge(jetptbins[1]);

  TH2D* v0Radius = (TH2D*)th3->Project3D("zy");
  v0Radius->SetName(TString::Format("v0Radius_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  // v0Radius->Rebin2D(doZ ? 1 : 5, rebinNumber);
  normaliseHistColByCol(v0Radius);
  v0Radius->SetMinimum(1e-5);
  v0Radius->SetMaximum(1.);
  histVector.push_back(v0Radius);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c",
                              dataSet.c_str(), lowjetpt, highjetpt).Data();
  // latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtRadius", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZRadius", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotCosPA(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( "V0" == hadron ) {
    cout << "V0 selected. Redirecting to plotV0CosPA" << endl;
    plotV0CosPA(inName, dataSet, jetptmin, jetptmax, doZ);
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int radiusAxis = 2;
  const int cospaAxis  = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 0.995, yMaxFrame = 1.;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.5, yLatex = 0.25;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("%s cos PA", formatHadronName(hadron).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  // For data, histograms do not contain "0" in the name
  string hadronForHist = hadron;
  if ("Lambda0" == hadronForHist) { hadronForHist = "Lambda"; }
  if ("AntiLambda0" == hadronForHist) { hadronForHist = "AntiLambda"; }
  histName = TString::Format("jetPt%sPtCosPA", hadronForHist.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjCosPA", hadronForHist.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(th3->GetXaxis(), jetptmin, jetptmax);
  th3->GetXaxis()->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = th3->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th3->GetXaxis()->GetBinUpEdge(jetptbins[1]);

  TH2D* v0CosPA = (TH2D*)th3->Project3D("zy");
  v0CosPA->SetName(TString::Format("v0CosPA_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  // v0CosPA->Rebin2D(doZ ? 1 : 5, rebinNumber);
  normaliseHistColByCol(v0CosPA);
  v0CosPA->SetMinimum(1e-5);
  v0CosPA->SetMaximum(1.);
  histVector.push_back(v0CosPA);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c", dataSet.c_str(), lowjetpt, highjetpt).Data();
  // latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtCosPA", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZCosPA", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotDCAdaughters(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
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
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 0., yMaxFrame = 1.;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();;
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("%s DCA daughters (cm)", formatHadronName(hadron).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  // For data, histograms do not contain "0" in the name
  string hadronForHist = hadron;
  if ("Lambda0" == hadronForHist) { hadronForHist = "Lambda"; }
  if ("AntiLambda0" == hadronForHist) { hadronForHist = "AntiLambda"; }
  histName = TString::Format("jetPt%sPtDCAd", hadronForHist.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjDCAd", hadronForHist.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(th3->GetXaxis(), jetptmin, jetptmax);
  th3->GetXaxis()->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = th3->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th3->GetXaxis()->GetBinUpEdge(jetptbins[1]);

  TH2D* v0DCAdaughters = (TH2D*)th3->Project3D("zy");
  v0DCAdaughters->SetName(TString::Format("v0DCAdaughters_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  normaliseHistColByCol(v0DCAdaughters);
  v0DCAdaughters->SetMinimum(1e-5);
  v0DCAdaughters->SetMaximum(1.);
  histVector.push_back(v0DCAdaughters);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c",
                              dataSet.c_str(), jetptmin, jetptmax).Data();
  // latex = CreateLatex(0.45, 0.8, latexText, textSize);
  latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtDCAdaughters", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZDCAdaughters", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotDCApos(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
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
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int dcaposAxis = 2;
  const int dcanegAxis = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = -10., yMaxFrame = 10.;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("%s DCA pos", formatHadronName(hadron).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  // For data, histograms do not contain "0" in the name
  string hadronForHist = hadron;
  if ("Lambda0" == hadronForHist) { hadronForHist = "Lambda"; }
  if ("AntiLambda0" == hadronForHist) { hadronForHist = "AntiLambda"; }
  histName = TString::Format("jetPt%sPtDCAposneg", hadronForHist.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjDCAposneg", hadronForHist.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0DCApos = (TH2D*)thn->Projection(dcaposAxis, v0ptAxis);
  v0DCApos->SetName(TString::Format("v0DCApos_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  normaliseHistColByCol(v0DCApos);
  v0DCApos->SetMinimum(1e-5);
  v0DCApos->SetMaximum(1.);
  histVector.push_back(v0DCApos);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c",
                              dataSet.c_str(), lowjetpt, highjetpt).Data();
  // latex = CreateLatex(0.4, 0.8, latexText, textSize);
  latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtDCApos", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZDCApos", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotDCAneg(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
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
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int dcaposAxis = 2;
  const int dcanegAxis = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = -10., yMaxFrame = 10.;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("%s DCA neg", formatHadronName(hadron).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // For data, histograms do not contain "0" in the name
  string hadronForHist = hadron;
  if ("Lambda0" == hadronForHist) { hadronForHist = "Lambda"; }
  if ("AntiLambda0" == hadronForHist) { hadronForHist = "AntiLambda"; }
  histName = TString::Format("jetPt%sPtDCAposneg", hadronForHist.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjDCAposneg", hadronForHist.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0DCAneg = (TH2D*)thn->Projection(dcanegAxis, v0ptAxis);
  v0DCAneg->SetName(TString::Format("v0DCAneg_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  normaliseHistColByCol(v0DCAneg);
  v0DCAneg->SetMinimum(1e-5);
  v0DCAneg->SetMaximum(1.);
  histVector.push_back(v0DCAneg);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c",
                              dataSet.c_str(), lowjetpt, highjetpt).Data();
  // TLatex* latex = CreateLatex(0.4, 0.8, latexText, textSize);
  TLatex* latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtDCAneg", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZDCAneg", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotctau(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ("V0" == hadron) {
    cout << "V0 selected. Use plotV0ctau instead" << endl;
    // plotV0ctau(inName, dataSet, hadron, jetptmin, jetptmax, doZ);
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 0., yMaxFrame = 40.;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.45, yLatex = 0.7;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("#it{c}#tau (%s)", formatHadronName(hadron).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // For data, histograms do not contain "0" in the name
  string hadronForHist = hadron;
  if ("Lambda0" == hadronForHist) { hadronForHist = "Lambda"; }
  if ("AntiLambda0" == hadronForHist) { hadronForHist = "AntiLambda"; }
  histName = TString::Format("jetPt%sPtCtau", hadronForHist.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjCtau", hadronForHist.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());
  th3->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(th3->GetXaxis(), jetptmin, jetptmax);
  th3->GetXaxis()->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = th3->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th3->GetXaxis()->GetBinUpEdge(jetptbins[1]);

  TH2D* ctau = (TH2D*)th3->Project3D("zy");
  // ctau->Rebin2D(doZ ? 1 : 5, rebinNumber);
  normaliseHistColByCol(ctau);
  ctau->SetMinimum(1e-5);
  ctau->SetMaximum(1.);
  ctau->SetName(TString::Format("ctau_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  histVector.push_back(ctau);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c",
                              dataSet.c_str(), lowjetpt, highjetpt).Data();
  // TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  TLatex* latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtctau", hadron.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZctau", hadron.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotMass(string inName = "", string dataSet = "dataSet", string hadron = "", string hypothesis = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ("V0" == hadron) { // Need different plotting function for V0s
    cout << "V0 selected. Redirecting to plotV0Mass" << endl;
    plotV0Mass(inName, dataSet, hypothesis, jetptmin, jetptmax, doZ);
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
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 1.05, yMaxFrame = 1.215;
  if (doZ) { xMinFrame = 1e-3; xMaxFrame = 1. + 1e-3; }
  if ("K0S" == hypothesis) { yMinFrame = 0.4, yMaxFrame = 0.6; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.45, yLatex = 0.7;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("#it{M} (%s)", formatHadronName(hypothesis).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // For data, histograms do not contain "0" in the name
  string hadronForHist = hadron;
  if ("Lambda0" == hadronForHist) { hadronForHist = "Lambda"; }
  if ("AntiLambda0" == hadronForHist) { hadronForHist = "AntiLambda"; }
  histName = TString::Format("jetPt%sPtAllMasses", hadronForHist.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjAllMasses", hadronForHist.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* mass = (TH2D*)thn->Projection(K0SAxis * ("K0S" == hypothesis) + LambdaAxis * ("Lambda0" == hypothesis) + AntiLambdaAxis * ("AntiLambda0" == hypothesis), v0ptAxis);
  mass->SetName("v0mass");
  // mass->Rebin2D(doZ ? 1 : 5, rebinNumber);
  normaliseHistColByCol(mass);
  mass->SetMinimum(1e-5);
  mass->SetMaximum(1.);
  histVector.push_back(mass);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c}{ %s tagged V0s } }",
                              dataSet.c_str(), jetptmin, jetptmax, hadron.c_str()).Data();
  // TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  TLatex* latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%sPtmass%s", hadron.c_str(), hypothesis.c_str()).Data();
  if (doZ) { saveName = TString::Format("%sZmass%s", hadron.c_str(), hypothesis.c_str()).Data(); }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotMasses(string inName = "", string dataSet = "dataSet", string hadron = "", std::array<string, 2> hypotheses = {"", ""}, double jetptmin = 10., double jetptmax = 200., double v0min = 0., double v0max = 100., bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ("V0" == hadron) {
    cout << "V0 selected. Use plotV0Masses instead" << endl;
    plotV0Masses(inName, dataSet, hypotheses, jetptmin, jetptmax, v0min, v0max, doZ);
    return;
  }
  if (("K0S" != hadron) + ("Lambda0" != hadron) + ("AntiLambda0" != hadron) != 1) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if (hypotheses[0] == hypotheses[1]) {
    cout << "Error: hypotheses must be different" << endl;
    return;
  }
  if ( ("K0S" == hypotheses[0]) + ("Lambda0" == hypotheses[0]) + ("AntiLambda0" == hypotheses[0]) != 1 ) {
    cout << "Error: hypothesis 0 must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if ( ("K0S" == hypotheses[1]) + ("Lambda0" == hypotheses[1]) + ("AntiLambda0" == hypotheses[1]) != 1 ) {
    cout << "Error: hypothesis 1 must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogZ = true;
  double xMinFrame = 1.05, xMaxFrame = 1.215, yMinFrame = 1.05, yMaxFrame = 1.215;
  if ("K0S" == hypotheses[0]) { xMinFrame = 0.4, xMaxFrame = 0.6; }
  if ("K0S" == hypotheses[1]) { yMinFrame = 0.4, yMaxFrame = 0.6; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.45, yLatex = 0.25;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M} (%s)", formatHadronName(hypotheses[0]).c_str()).Data();
  yTitle = TString::Format("#it{M} (%s)", formatHadronName(hypotheses[1]).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = TString::Format("partJetPt%sPtDetJetPt%sPtAllMasses", hadron.c_str(), hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProjAllMasses", hadron.c_str(), hadron.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0bins = getProjectionBins(thn->GetAxis(v0ptAxis), v0min, v0max);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(v0ptAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0 = thn->GetAxis(v0ptAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0ptAxis)->GetBinUpEdge(v0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; }

  int axisX = K0SAxis * ("K0S" == hypotheses[0]) + LambdaAxis * ("Lambda0" == hypotheses[0]) + AntiLambdaAxis * ("AntiLambda0" == hypotheses[0]);
  int axisY = K0SAxis * ("K0S" == hypotheses[1]) + LambdaAxis * ("Lambda0" == hypotheses[1]) + AntiLambdaAxis * ("AntiLambda0" == hypotheses[1]);
  TH2D* mass = (TH2D*)thn->Projection(axisY, axisX);
  mass->SetName("v0mass");
  mass->Scale(1./mass->Integral());
  mass->SetMinimum(1e-5);
  mass->SetMaximum(1.);
  histVector.push_back(mass);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/#it{c} }}",
                              dataSet.c_str(), jetptmin, jetptmax, hadron.c_str(), lowv0, highv0).Data();
  if (doZ) {
    latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c }{ #it{z}_{%s} = %.0f - %.0f }}",
                                dataSet.c_str(), jetptmin, jetptmax, hadron.c_str(), lowv0, highv0).Data();
  }
  // TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  TLatex* latex = nullptr;
  frame->SetTitle(latexText.c_str());

  saveName = TString::Format("%smass%s-%s", hadron.c_str(), hypotheses[0].c_str(), hypotheses[1].c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  if (doZ) { saveName = TString::Format("%s_v0z%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  else { saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

/*
void plotDCAposneg(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
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
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int dcaposAxis = 2;
  const int dcanegAxis = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  // string dataSet = "LHC23y_pass1";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogz = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = -10, yMaxFrame = 10;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.2, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("%s DCA pos", formatHadronName(hadron).c_str()).Data();
  yTitle = TString::Format("%s DCA neg", formatHadronName(hadron).c_str()).Data();

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = TString::Format("jet-fragmentation/data/jets/V0/jetPt%sPtDCAposneg", hadron.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);

  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);

  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);



  TH2D* v0DCAposneg = (TH2D*)thn->Projection(dcanegAxis, dcaposAxis);
  v0DCAposneg->SetName(TString::Format("v0DCAposneg_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  v0DCAposneg->Scale(1./v0DCAposneg->Integral());
  v0DCAposneg->SetMinimum(1e-5);
  v0DCAposneg->SetMaximum(1.);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c", dataSet.c_str(), lowjetpt, highjetpt).Data();
  // TLatex* latex = CreateLatex(0.2, 0.93, latexText, textSize);
  frame->SetTitle

  saveName = TString::Format("%sDCAposneg", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);

  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->cd();
  frame->Draw();
  v0DCAposneg->Draw("same colz");
  latex->Draw("same");
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}


// ----------------------------------------------------------
// Unidentified V0s
// ----------------------------------------------------------
void plotRadius(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., bool doZ = false)
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
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int radiusAxis = 2;
  const int cospaAxis  = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 0., yMaxFrame = 50.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  if (doZ) { xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = TString::Format("%s Radius (cm)", formatHadronName(hadron).c_str()).Data();
  // string dataSet = "LHC23y_pass1";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  // For data, histograms do not contain "0" in the name
  string hadronForHist = hadron;
  if ("Lambda0" == hadronForHist) { hadronForHist = "Lambda"; }
  if ("AntiLambda0" == hadronForHist) { hadronForHist = "AntiLambda"; }
  histName = TString::Format("jetPt%sPtRadiusCosPA", hadronForHist.c_str()).Data();
  if (doZ) { histName = TString::Format("jetPt%sTrackProjRadiusCosPA", hadronForHist.c_str()).Data(); }

  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH2D* v0Radius = (TH2D*)thn->Projection(radiusAxis, v0ptAxis);
  v0Radius->SetName(TString::Format("v0Radius_jetpt%.0f-%.0f", lowjetpt, highjetpt).Data());
  // v0Radius->Rebin2D(rebinNumber, rebinNumber);
  normaliseHistColByCol(v0Radius);
  histVector.push_back(v0Radius);

  latexText = TString::Format("%s, #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c",
                              dataSet.c_str(), lowjetpt, highjetpt).Data();
  // latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  frame->SetTitle

  saveName = TString::Format("%sRadius", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

// This plot isn't very useful. ctauL = 2*ctauK0S
void plotV0ctauKL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;


  gStyle->SetNdivisions(505, "xy");

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
  // string dataSet = "LHC23y_pass1";

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtCtau";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);

  TH2D* ctau = (TH2D*)thn->Projection(LambdaAxis, K0SAxis);
  ctau->Scale(1./ctau->Integral());

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax).Data();

  saveName = "ctauK0SLambda";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);

  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->cd();
  frame->Draw();
  ctau->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0Mass(string inName = "AnalysisResults.root", int setting = 2, double jetptmin = 10., double jetptmax = 200.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;


  gStyle->SetNdivisions(505, "xy");
  if (setting < K0SAxis || setting > AntiLambdaAxis) {
    cout << "Invalid setting! 2 for K0S, 3 for Lambda, 4 for antiLambda" << endl;
    return;
  }

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1e-4, yMaxFrame = .1;
  if (setting == LambdaAxis || setting == AntiLambdaAxis) {
    xMinFrame = 1.05, xMaxFrame = 1.215;
    yMinFrame = 1e-3, yMaxFrame = 0.1;
  }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{M} (K^{0}_{S})";
  if (setting == LambdaAxis) { xTitle = "#it{M} (#Lambda)"; }
  if (setting == AntiLambdaAxis) { xTitle = "#it{M} (#bar{#Lambda})"; }
  xTitle = "normalised count";
  // string dataSet = "LHC23y_pass1";
  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax).Data();

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);

  TH2D* mass = (TH2D*)thn->Projection(setting);
  mass->Scale(1./mass->Integral());
  mass->SetName(TString::Format("mass_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  histVector.push_back(mass);

  saveName = "massK0S";
  if (setting == LambdaAxis) { saveName = "massLambda"; }
  else if (setting == AntiLambdaAxis) { saveName = "massAntiLambda"; }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);

  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0MassKL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;


  gStyle->SetNdivisions(505, "xy");

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
  // string dataSet = "LHC23y_pass1";
  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax).Data();

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);

  TH2D* mass = (TH2D*)thn->Projection(LambdaAxis, K0SAxis);
  mass->SetName(TString::Format("massKL_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  mass->Rebin2D(5,5);
  mass->Scale(1./mass->Integral());

  saveName = "massK0SLambda";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);

  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0MassKaL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;


  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1.015, yMaxFrame = 1.215;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{M} (K^{0}_{S})";
  yTitle = "#it{M} (#bar{#Lambda})";
  // string dataSet = "LHC23y_pass1";
  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax).Data();

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);

  TH2D* mass = (TH2D*)thn->Projection(AntiLambdaAxis, K0SAxis);
  mass->SetName(TString::Format("massKaL_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  mass->Rebin2D(5,5);
  mass->Scale(1./mass->Integral());

  saveName = "massK0SAntiLambda";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);

  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0MassLL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;


  gStyle->SetNdivisions(505, "xy");

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
  // string dataSet = "LHC23y_pass1";
  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax).Data();

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);

  TH2D* mass = (TH2D*)thn->Projection(AntiLambdaAxis, LambdaAxis);
  mass->SetName(TString::Format("massLL_jetpt%.0f-%.0f", jetptmin, jetptmax).Data());
  mass->Rebin2D(5,5);
  mass->Scale(1./mass->Integral());

  saveName = "massLambdaAntiLambda";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);

  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotmassWithCuts(string inName = "AnalysisResults.root", int setting = 3, double jetptmin = 10., double jetptmax = 200.,
                      double K0SMin = 0.4, double K0SMax = 0.6, double LambdaMin = 1.015, double LambdaMax = 1.215, double AntiLambdaMin = 1.015, double AntiLambdaMax = 1.215)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;


  gStyle->SetNdivisions(505, "xy");
  if (setting < K0SAxis || setting > AntiLambdaAxis) {
    cout << "Invalid setting! 2 for K0S, 3 for Lambda, 4 for antiLambda" << endl;
    return;
  }

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1e-4, yMaxFrame = .1;
  if (setting == LambdaAxis || setting == AntiLambdaAxis) {
    xMinFrame = 1.05, xMaxFrame = 1.215;
    yMinFrame = 1e-3, yMaxFrame = 0.1;
  }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  xTitle = "#it{M} (K^{0}_{S})";
  if (setting == LambdaAxis) { xTitle = "#it{M} (#Lambda)"; }
  if (setting == AntiLambdaAxis) { xTitle = "#it{M} (#bar{#Lambda})"; }
  xTitle = "normalised count";
  // string dataSet = "LHC23y_pass1";
  latexText = "";

  std::vector<TH2D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  int firstBinK0SMass = 1, lastBinK0SMass = thn->GetAxis(K0SAxis)->GetNbins();
  int firstBinLambdaMass = 1, lastBinLambdaMass = thn->GetAxis(LambdaAxis)->GetNbins();
  int firstBinAntiLambdaMass = 1, lastBinAntiLambdaMass = thn->GetAxis(AntiLambdaAxis)->GetNbins();

  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);
  firstBinK0SMass = thn->GetAxis(K0SAxis)->FindBin(K0SMin + 1e-3);
  lastBinK0SMass  = thn->GetAxis(K0SAxis)->FindBin(K0SMax - 1e-3);
  firstBinLambdaMass = thn->GetAxis(LambdaAxis)->FindBin(LambdaMin + 1e-3);
  lastBinLambdaMass  = thn->GetAxis(LambdaAxis)->FindBin(LambdaMax - 1e-3);
  firstBinAntiLambdaMass = thn->GetAxis(AntiLambdaAxis)->FindBin(AntiLambdaMin + 1e-3);
  lastBinAntiLambdaMass  = thn->GetAxis(AntiLambdaAxis)->FindBin(AntiLambdaMax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  TH2D* onlyJet = (TH2D*)thn->Projection(setting);
  onlyJet->SetName("onlyJet");
  onlyJet->Scale(1./onlyJet->Integral());
  onlyJet->Rebin(5);
  setStyle(onlyJet, 0);
  legend->AddEntry(onlyJet, TString::Format("#it{p}_{T, ch+V0 jet} = %.0f - %.0f GeV/#it{c}", jetptmin, jetptmax).Data());
  histVector.push_back(onlyJet);

  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);
  TH2D* jetAndV0 = (TH2D*)thn->Projection(setting);
  jetAndV0->SetName("jetAndV0");
  jetAndV0->Scale(1./jetAndV0->Integral());
  jetAndV0->Rebin(5);
  setStyle(jetAndV0, 1);
  legend->AddEntry(jetAndV0, TString::Format("& #it{p}_{T, V0} = %.0f - %.0f GeV/#it{c}").Data());
  histVector.push_back(jetAndV0);

  thn->GetAxis(K0SAxis)->SetRange(firstBinK0SMass, lastBinK0SMass);
  thn->GetAxis(LambdaAxis)->SetRange(firstBinLambdaMass, lastBinLambdaMass);
  thn->GetAxis(AntiLambdaAxis)->SetRange(firstBinAntiLambdaMass, lastBinAntiLambdaMass);
  TH2D* mass = (TH2D*)thn->Projection(setting);
  mass->SetName("mass");
  mass->Scale(1./mass->Integral());
  mass->Rebin(5);
  setStyle(mass, 2);
  legend->AddEntry(mass, TString::Format("& #it{M}_{K^{0}_{S}} = %.0f - %.0f, #it{M}_{#Lambda} = %.0f - %.0f, #it{M}_{#bar{#Lambda}} = %.0f - %.0f", K0SMin, K0SMax, LambdaMin, LambdaMax, AntiLambdaMin, AntiLambdaMax).Data());
  histVector.push_back(mass);

  saveName = "massK0S";
  if (setting == LambdaAxis) { saveName = "massLambda"; }
  else if (setting == AntiLambdaAxis) { saveName = "massAntiLambda"; }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);

  saveName = TString::Format("%s_mK0S%.0f-%.0f", saveName.c_str(), K0SMin, K0SMax);
  saveName = TString::Format("%s_mL%.0f-%.0f", saveName.c_str(), LambdaMin, LambdaMax);
  saveName = TString::Format("%s_mAL%.0f-%.0f", saveName.c_str(), AntiLambdaMin, AntiLambdaMax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, saveName, "", latexText);
}
*/
// ----------------------------------------------------------
// Lambda0
// ----------------------------------------------------------
void plotLambda0Radius(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotRadius(inName, dataSet, "Lambda0", jetptmin, jetptmax); }
void plotLambda0CosPA(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotCosPA(inName, dataSet, "Lambda0", jetptmin, jetptmax); }
void plotLambda0DCAdaughters(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotDCAdaughters(inName, dataSet, "Lambda0", jetptmin, jetptmax); }
void plotLambda0DCApos(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotDCApos(inName, dataSet, "Lambda0", jetptmin, jetptmax); }
void plotLambda0DCAneg(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotDCAneg(inName, dataSet, "Lambda0", jetptmin, jetptmax); }
// void plotLambda0DCAposneg(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
// {, dataSet plotDCAposneg(inName, "Lambda0", jetptmin, jetptmax); }
void plotLambda0ctau(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotctau(inName, dataSet, "Lambda0", jetptmin, jetptmax); }
void plotLambda0Mass(string inName = "AnalysisResults.root", string dataSet = "dataSet", string hypothesis = "", double jetptmin = 10., double jetptmax = 200.)
{ plotMass(inName, dataSet, "Lambda0", hypothesis, jetptmin, jetptmax); }

// ----------------------------------------------------------
// AntiLambda0
// ----------------------------------------------------------
void plotAntiLambda0Radius(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotRadius(inName, dataSet, "AntiLambda0", jetptmin, jetptmax); }
void plotAntiLambda0CosPA(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotCosPA(inName, dataSet, "AntiLambda0", jetptmin, jetptmax); }
void plotAntiLambda0DCAdaughters(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotDCAdaughters(inName, dataSet, "AntiLambda0", jetptmin, jetptmax); }
void plotAntiLambda0DCApos(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotDCApos(inName, dataSet, "AntiLambda0", jetptmin, jetptmax); }
void plotAntiLambda0DCAneg(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotDCAneg(inName, dataSet, "AntiLambda0", jetptmin, jetptmax); }
// void plotAntiLambda0DCAposneg(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
// {, dataSet plotDCAposneg(inName, "AntiLambda0", jetptmin, jetptmax); }
void plotAntiLambda0ctau(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotctau(inName, dataSet, "AntiLambda0", jetptmin, jetptmax); }
void plotAntiLambda0Mass(string inName = "AnalysisResults.root", string dataSet = "dataSet", string hypothesis = "", double jetptmin = 10., double jetptmax = 200.)
{ plotMass(inName, dataSet, "AntiLambda0", hypothesis, jetptmin, jetptmax); }

// ----------------------------------------------------------
// K0S
// ----------------------------------------------------------
void plotK0SRadius(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotRadius(inName, dataSet, "K0S", jetptmin, jetptmax); }
void plotK0SCosPA(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotCosPA(inName, dataSet, "K0S", jetptmin, jetptmax); }
void plotK0SDCAdaughters(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotDCAdaughters(inName, dataSet, "K0S", jetptmin, jetptmax); }
void plotK0SDCApos(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotDCApos(inName, dataSet, "K0S", jetptmin, jetptmax); }
void plotK0SDCAneg(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotDCAneg(inName, dataSet, "K0S", jetptmin, jetptmax); }
// void plotK0SDCAposneg(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
// {, dataSet plotDCAposneg(inName, "K0S", jetptmin, jetptmax); }
void plotK0Sctau(string inName = "AnalysisResults.root", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{ plotctau(inName, dataSet, "K0S", jetptmin, jetptmax); }
void plotK0SMass(string inName = "AnalysisResults.root", string dataSet = "dataSet", string hypothesis = "", double jetptmin = 10., double jetptmax = 200.)
{ plotMass(inName, dataSet, "K0S", hypothesis, jetptmin, jetptmax); }

// ----------------------------------------------------------
// ----------------------------------------------------------
// ----------------------------------------------------------

void plot22o(double jetptmin, double jetptmax, string hadron, int setting, bool doZ)
{
  string inName = "~/cernbox/TrainOutput/252064/AnalysisResults.root";
  string dataSet = "LHC22o_pass6";

  switch(setting) {
    case 0:
      plotRadius(inName, dataSet, hadron, jetptmin, jetptmax, doZ);
      break;
    case 1:
      plotCosPA(inName, dataSet, hadron, jetptmin, jetptmax, doZ);
      break;
    case 2:
      plotDCAdaughters(inName, dataSet, hadron, jetptmin, jetptmax, doZ);
      break;
    case 3:
      plotDCApos(inName, dataSet, hadron, jetptmin, jetptmax, doZ);
      break;
    case 4:
      plotDCAneg(inName, dataSet, hadron, jetptmin, jetptmax, doZ);
      break;
    case 5:
      plotV0ctau(inName, dataSet, hadron, jetptmin, jetptmax, doZ);
      break;
    default:
      cout << "Invalid setting! Pick a number from [0, 4]." << endl;
      break;
  }
}
