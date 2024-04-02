
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

double getNjets(TFile* inFile, double jetptmin, double jetptmax)
{
  string histName = "jet-fragmentation/data/jets/jetPtEtaPhi";
  TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
  std::array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
  return jets->Integral(jetptbins[0], jetptbins[1], 1, jets->GetNbinsY(), 1, jets->GetNbinsZ());
}

// ----------------------------------------------------------

void plotPt(string inName = "AnalysisResults.root", string hadron = "", double jetptmin = 10., double jetptmax = 200.)
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
  const int v0etaAxis  = 2;
  const int v0phiAxis  = 3;

  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogy = true;
  double xMinFrame = 0., xMaxFrame = 100., yMinFrame = 1e-6, yMaxFrame = 0.1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";
  dataSet = "LHC23y_pass1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("jetPt%sPtEtaPhi", hadron.c_str()).Data();
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH1D* v0pt = (TH1D*)thn->Projection(v0ptAxis);
  v0pt->Scale(1./v0pt->Integral());
  // normaliseHistRowByRow(v0pt);
  // v0pt->SetMinimum(1e-5);
  // v0pt->SetMaximum(1.);
  v0pt->SetName(TString::Format("%sPt", hadron.c_str()).Data());
  setStyle(v0pt, 0);
  histVector.push_back(v0pt);

  latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%sPt", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotZ(string inName = "AnalysisResults.root", string hadron = "", double jetptmin = 10., double jetptmax = 200.)
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
  const int v0zAxis    = 1;
  const int dcaposAxis = 2;
  const int dcanegAxis = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogy = true;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 1e-4, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data();
  yTitle = TString::Format("#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}_{%s}}", formatHadronName(hadron).c_str()).Data();
  dataSet = "LHC23y_pass1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("jetPt%sTrackProjDCAposneg", hadron.c_str()).Data();
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH1D* v0z = (TH1D*)thn->Projection(v0zAxis);
  v0z->Scale(1./getNjets(inFile, jetptmin, jetptmax), "width");
  v0z->SetName(TString::Format("%sZ", hadron.c_str()).Data());
  setStyle(v0z, 0);
  histVector.push_back(v0z);

  latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%sZ", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void plotV0Z(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  string dataSet = "LHC23y_pass1";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 1e-4, yMaxFrame = 0.2;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.4, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}_{V0}";
  yTitle = "#frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{z}_{V0}}";

  std::vector<TH1D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/data/jets/V0/jetPtV0TrackProj";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(th2->GetXaxis(), jetptmin, jetptmax);
  TH1D* v0z = (TH1D*)th2->ProjectionY("v0z", jetptbins[0], jetptbins[1]);
  v0z->Scale(1./getNjets(inFile, jetptmin, jetptmax), "width");
  setStyle(v0z, 0);
  histVector.push_back(v0z);

  double lowjetpt = th2->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th2->GetXaxis()->GetBinUpEdge(jetptbins[1]);
  latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "v0z";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}

// ----------------------------------------------------------

void plotV0Radius(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int radiusAxis = 2;
  const int cospaAxis  = 3;

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
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtRadiusCosPA";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);

  TH1D* v0Radius = (TH1D*)thn->Projection(radiusAxis);
  v0Radius->SetName(TString::Format("v0Radius_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  v0Radius->Scale(1./v0Radius->Integral());
  setStyle(v0Radius, 0);
  histVector.push_back(v0Radius);

  saveName = "v0Radius";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0CosPA(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int radiusAxis = 2;
  const int cospaAxis  = 3;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0.95, xMaxFrame = 1., yMinFrame = 1e-3, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "V0 Radius";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtRadiusCosPA";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);

  TH1D* v0CosPA = (TH1D*)thn->Projection(cospaAxis);
  v0CosPA->SetName(TString::Format("v0CosPA_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  v0CosPA->Scale(1./v0CosPA->Integral());
  setStyle(v0CosPA, 0);
  histVector.push_back(v0CosPA);

  saveName = "v0CosPA";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0DCAdaughters(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 1., yMinFrame = 5e-3, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "DCA daughters";
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtDCAd";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = th3->GetNbinsX();
  int firstBinV0Pt = 1, lastBinV0Pt = th3->GetNbinsY();
  firstBinJetPt = th3->GetXaxis()->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = th3->GetXaxis()->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = th3->GetYaxis()->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = th3->GetYaxis()->FindBin(v0ptmax - 1e-3);
  // th3->GetXaxis()->SetRange(firstBinJetPt, lastBinJetPt);

  TH1D* v0DCAdaughters = (TH1D*)th3->ProjectionZ(TString::Format("v0DCAdaughters_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data(), firstBinJetPt, lastBinJetPt, firstBinV0Pt, lastBinV0Pt);
  v0DCAdaughters->Scale(1./v0DCAdaughters->Integral());
  setStyle(v0DCAdaughters, 0);
  histVector.push_back(v0DCAdaughters);

  saveName = "v0DCAdaughters";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0DCApos(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int dcaposAxis = 2;
  const int dcanegAxis = 3;

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
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtDCAposneg";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);
  TH1D* v0DCApos = (TH1D*)thn->Projection(dcaposAxis);
  v0DCApos->SetName(TString::Format("v0DCApos_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  v0DCApos->Scale(1./v0DCApos->Integral());
  setStyle(v0DCApos, 0);
  histVector.push_back(v0DCApos);

  saveName = "v0DCApos";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0DCAneg(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int dcaposAxis = 2;
  const int dcanegAxis = 3;

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
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtDCAposneg";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);
  TH1D* v0DCAneg = (TH1D*)thn->Projection(dcanegAxis);
  v0DCAneg->SetName(TString::Format("v0DCAneg_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  v0DCAneg->Scale(1./v0DCAneg->Integral());
  setStyle(v0DCAneg, 0);
  histVector.push_back(v0DCAneg);

  saveName = "v0DCAneg";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0DCAposneg(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int v0ptAxis   = 1;
  const int dcaposAxis = 2;
  const int dcanegAxis = 3;

  double time = clock();
  gStyle->SetNdivisions(505);

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogz = true;
  double xMinFrame = -10., xMaxFrame = 10., yMinFrame = -10, yMaxFrame = 10;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "DCA pos";
  yTitle = "DCA neg";
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtDCAposneg";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinV0Pt = 1, lastBinV0Pt = thn->GetAxis(v0ptAxis)->GetNbins();
  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);
  firstBinV0Pt = thn->GetAxis(v0ptAxis)->FindBin(v0ptmin + 1e-3);
  lastBinV0Pt  = thn->GetAxis(v0ptAxis)->FindBin(v0ptmax - 1e-3);

  thn->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);
  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);
  TH2D* v0DCAposneg = (TH2D*)thn->Projection(dcanegAxis, dcaposAxis);
  v0DCAposneg->SetName(TString::Format("v0DCAposneg_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  v0DCAposneg->Scale(1./v0DCAposneg->Integral());
  setStyle(v0DCAposneg, 0);
  v0DCAposneg->SetMinimum(1e-5);
  v0DCAposneg->SetMaximum(1.);

  saveName = "v0DCAposneg";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  v0DCAposneg->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0ctau(string inName = "AnalysisResults.root", int setting = 2, double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

  double time = clock();
  gStyle->SetNdivisions(505);
  if (setting < K0SAxis || setting > AntiLambdaAxis) {
    cout << "Invalid setting! 2 for K0S, 3 for Lambda, 4 for antiLambda" << endl;
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
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtCtau";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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

  TH1D* ctau = (TH1D*)thn->Projection(setting);
  ctau->Scale(1./ctau->Integral());
  ctau->SetName(TString::Format("ctau_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  setStyle(ctau, 0);
  histVector.push_back(ctau);

  saveName = "ctauK0S";
  if (setting == LambdaAxis) { saveName = "ctauLambda"; }
  else if (setting == AntiLambdaAxis) { saveName = "ctauAntiLambda"; }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
// This plot isn't very useful. ctauL = 2*ctauK0S?
void plotV0ctauKL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

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
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtCtau";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
  setStyle(ctau, 0);

  saveName = "ctauK0SLambda";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  ctau->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0mass(string inName = "AnalysisResults.root", int setting = 2, double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

  double time = clock();
  gStyle->SetNdivisions(505);
  if (setting < K0SAxis || setting > AntiLambdaAxis) {
    cout << "Invalid setting! 2 for K0S, 3 for Lambda, 4 for antiLambda" << endl;
    return;
  }

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
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
  yTitle = "normalised count";
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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

  TH1D* mass = (TH1D*)thn->Projection(setting);
  mass->Scale(1./mass->Integral());
  mass->SetName(TString::Format("mass_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  setStyle(mass, 0);
  histVector.push_back(mass);

  saveName = "massK0S";
  if (setting == LambdaAxis) { saveName = "massLambda"; }
  else if (setting == AntiLambdaAxis) { saveName = "massAntiLambda"; }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0massKL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

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
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
  mass->SetName(TString::Format("massKL_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  mass->Rebin2D(5,5);
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);

  saveName = "massK0SLambda";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0massKaL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

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
  yTitle = "#it{M} (#bar{#Lambda})";
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
  mass->SetName(TString::Format("massKaL_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  mass->Rebin2D(5,5);
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);

  saveName = "massK0SAntiLambda";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0massLL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

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
  latexText = TString::Format("#splitline{ LHC23y_pass1 }{ #splitline{ #it{p}_{T, jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
  mass->SetName(TString::Format("massLL_jetpt%.0f-%.0f_v0pt%.0f-%.0f", jetptmin, jetptmax, v0ptmin, v0ptmax).Data());
  mass->Rebin2D(5,5);
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);

  saveName = "massLambdaAntiLambda";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  myCanvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  myCanvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

// ----------------------------------------------------------
void plotmassWithCuts(string inName = "AnalysisResults.root", int setting = 3, double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 10.,
                      double K0SMin = 0.4, double K0SMax = 0.6, double LambdaMin = 1.015, double LambdaMax = 1.215, double AntiLambdaMin = 1.015, double AntiLambdaMax = 1.215)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0ptAxis       = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

  double time = clock();
  gStyle->SetNdivisions(505);
  if (setting < K0SAxis || setting > AntiLambdaAxis) {
    cout << "Invalid setting! 2 for K0S, 3 for Lambda, 4 for antiLambda" << endl;
    return;
  }

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
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
  yTitle = "normalised count";
  latexText = "";

  std::vector<TH1D*> histVector;

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtMass";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
  TH1D* onlyJet = (TH1D*)thn->Projection(setting);
  onlyJet->SetName("onlyJet");
  onlyJet->Scale(1./onlyJet->Integral());
  onlyJet->Rebin(5);
  setStyle(onlyJet, 0);
  legend->AddEntry(onlyJet, TString::Format("#it{p}_{T, jet} = %.0f - %.0f GeV/#it{c}", jetptmin, jetptmax).Data());
  histVector.push_back(onlyJet);

  thn->GetAxis(v0ptAxis)->SetRange(firstBinV0Pt, lastBinV0Pt);
  TH1D* jetAndV0 = (TH1D*)thn->Projection(setting);
  jetAndV0->SetName("jetAndV0");
  jetAndV0->Scale(1./jetAndV0->Integral());
  jetAndV0->Rebin(5);
  setStyle(jetAndV0, 1);
  legend->AddEntry(jetAndV0, TString::Format("& #it{p}_{T, V0} = %.0f - %.0f GeV/#it{c}", v0ptmin, v0ptmax).Data());
  histVector.push_back(jetAndV0);

  thn->GetAxis(K0SAxis)->SetRange(firstBinK0SMass, lastBinK0SMass);
  thn->GetAxis(LambdaAxis)->SetRange(firstBinLambdaMass, lastBinLambdaMass);
  thn->GetAxis(AntiLambdaAxis)->SetRange(firstBinAntiLambdaMass, lastBinAntiLambdaMass);
  TH1D* mass = (TH1D*)thn->Projection(setting);
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
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s_mK0S%.0f-%.0f", saveName.c_str(), K0SMin, K0SMax);
  saveName = TString::Format("%s_mL%.0f-%.0f", saveName.c_str(), LambdaMin, LambdaMax);
  saveName = TString::Format("%s_mAL%.0f-%.0f", saveName.c_str(), AntiLambdaMin, AntiLambdaMax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
}
