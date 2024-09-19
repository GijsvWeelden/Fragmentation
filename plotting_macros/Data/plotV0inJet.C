
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

void plotPt(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200.)
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
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = true;
  double xMinFrame = 0., xMaxFrame = 100., yMinFrame = 1e-7, yMaxFrame = 0.1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.35, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";
  // dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  // histName = TString::Format("jetPt%sPtEtaPhi", hadron.c_str()).Data();
  string hadronForHistname = hadron;
  if ("Lambda0" == hadron) { hadronForHistname = "Lambda"; }
  if ("AntiLambda0" == hadron) { hadronForHistname = "AntiLambda"; }
  histName = TString::Format("jetPt%sPtDCAposneg", hadronForHistname.c_str()).Data();
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

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ d#it{R}(jet, %s) #leq #it{R}} }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  // latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%sPt", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
void plotV0Z(string inName = "", string dataSet = "dataSet", double jetptmin = 10., double jetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  string hadron = "V0";
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
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogy = true;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 1e-4, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.4, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data();
  yTitle = TString::Format("#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}_{%s}}", formatHadronName(hadron).c_str()).Data();
  // dataSet = "LHC22o_pass6_minBias";

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

  // latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ d#it{R}(jet, %s) #leq #it{R}} }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%sZ", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
void plotZ(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ("V0" == hadron) {
    cout << "V0 selected. Redirecting to plotV0Z()" << endl;
    plotV0Z(inName, jetptmin, jetptmax);
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim                = 5;
  const int jetptAxis           = 0;
  const int v0zAxis             = 1;
  const int K0SMassAxis         = 2;
  const int Lambda0MassAxis     = 3;
  const int AntiLambda0MassAxis = 4;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogy = true;
  // bool SetLogy = false;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 1e-7, yMaxFrame = 1e-3;
  // double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 0, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.4, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 4;
  xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data();
  yTitle = TString::Format("#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}_{%s}}", formatHadronName(hadron).c_str()).Data();
  // dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  string hadronForHistname = hadron;
  if ("Lambda0" == hadron) { hadronForHistname = "Lambda"; }
  if ("AntiLambda0" == hadron) { hadronForHistname = "AntiLambda"; }
  histName = TString::Format("jetPt%sTrackProjAllMasses", hadronForHistname.c_str()).Data();
  histName = TString::Format("jet-fragmentation/data/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH1D* v0z = (TH1D*)thn->Projection(v0zAxis);
  v0z->Rebin(rebinNumber);
  v0z->Scale(1./getNjets(inFile, jetptmin, jetptmax), "width");
  v0z->SetName(TString::Format("%sZ", hadron.c_str()).Data());
  setStyle(v0z, 0);
  histVector.push_back(v0z);

  // latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ d#it{R}(jet, %s) #leq #it{R}} }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%sZ", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}

// ----------------------------------------------------------

void plotRadius(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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

  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 100., yMinFrame = 1e-6, yMaxFrame = 0.5;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("%s Radius (cm)", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";
  // string dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("jet-fragmentation/data/jets/V0/jetPt%sPtRadiusCosPA", hadron.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0ptbins = getProjectionBins(thn->GetAxis(v0ptAxis), v0ptmin, v0ptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(v0ptAxis)->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0pt = thn->GetAxis(v0ptAxis)->GetBinLowEdge(v0ptbins[0]);
  double highv0pt = thn->GetAxis(v0ptAxis)->GetBinUpEdge(v0ptbins[1]);

  TH1D* v0Radius = (TH1D*)thn->Projection(radiusAxis);
  v0Radius->SetName(TString::Format("v0Radius_jetpt%.0f-%.0f_v0pt%.0f-%.0f", lowjetpt, highjetpt, lowv0pt, highv0pt).Data());
  v0Radius->Scale(1./v0Radius->Integral());
  setStyle(v0Radius, 0);
  histVector.push_back(v0Radius);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/c } }",
                              dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0pt, highv0pt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%sRadius", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), lowv0pt, highv0pt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");

}
void plotCosPA(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  bool setLogY = true;
  double xMinFrame = 0.95, xMaxFrame = 1., yMinFrame = 1e-3, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("%s cos PA", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";
  // string dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("jet-fragmentation/data/jets/V0/jetPt%sPtRadiusCosPA", hadron.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0ptbins = getProjectionBins(thn->GetAxis(v0ptAxis), v0ptmin, v0ptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(v0ptAxis)->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0pt = thn->GetAxis(v0ptAxis)->GetBinLowEdge(v0ptbins[0]);
  double highv0pt = thn->GetAxis(v0ptAxis)->GetBinUpEdge(v0ptbins[1]);

  TH1D* v0CosPA = (TH1D*)thn->Projection(cospaAxis);
  v0CosPA->SetName(TString::Format("v0CosPA_jetpt%.0f-%.0f_v0pt%.0f-%.0f", lowjetpt, highjetpt, lowv0pt, highv0pt).Data());
  v0CosPA->Scale(1./v0CosPA->Integral());
  setStyle(v0CosPA, 0);
  histVector.push_back(v0CosPA);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/c } }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0pt, highv0pt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%sCosPA", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), lowv0pt, highv0pt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
void plotDCAdaughters(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 1., yMinFrame = 5e-3, yMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("%s DCA daughters (cm)", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";
  // string dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("jet-fragmentation/data/jets/V0/jetPt%sPtDCAd", hadron.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(th3->GetXaxis(), jetptmin, jetptmax);
  std::array<int, 2> v0ptbins = getProjectionBins(th3->GetYaxis(), v0ptmin, v0ptmax);
  th3->GetXaxis()->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt = th3->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th3->GetXaxis()->GetBinUpEdge(jetptbins[1]);
  th3->GetYaxis()->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowv0pt = th3->GetYaxis()->GetBinLowEdge(v0ptbins[0]);
  double highv0pt = th3->GetYaxis()->GetBinUpEdge(v0ptbins[1]);

  TH1D* v0DCAdaughters = (TH1D*)th3->ProjectionZ(TString::Format("v0DCAdaughters_jetpt%.0f-%.0f_v0pt%.0f-%.0f", lowjetpt, highjetpt, lowv0pt, highv0pt).Data(), jetptbins[0], jetptbins[1], v0ptbins[0], v0ptbins[1]);
  v0DCAdaughters->Scale(1./v0DCAdaughters->Integral());
  setStyle(v0DCAdaughters, 0);
  histVector.push_back(v0DCAdaughters);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax, formatHadronName(hadron).c_str(), v0ptmin, v0ptmax).Data();
  latex = CreateLatex(0.3, 0.8, latexText, textSize);

  saveName = TString::Format("%sDCAdaughters", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
void plotDCApos(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  bool setLogY = true;
  double xMinFrame = -10., xMaxFrame = 10., yMinFrame = 1e-5, yMaxFrame = 0.5;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("%s DCA pos", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";
  // string dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("jet-fragmentation/data/jets/V0/jetPt%sPtDCAposneg", hadron.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0ptbins = getProjectionBins(thn->GetAxis(v0ptAxis), v0ptmin, v0ptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(v0ptAxis)->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0pt = thn->GetAxis(v0ptAxis)->GetBinLowEdge(v0ptbins[0]);
  double highv0pt = thn->GetAxis(v0ptAxis)->GetBinUpEdge(v0ptbins[1]);

  TH1D* v0DCApos = (TH1D*)thn->Projection(dcaposAxis);
  v0DCApos->SetName(TString::Format("v0DCApos_jetpt%.0f-%.0f_v0pt%.0f-%.0f", lowjetpt, highjetpt, lowv0pt, highv0pt).Data());
  v0DCApos->Scale(1./v0DCApos->Integral());
  setStyle(v0DCApos, 0);
  histVector.push_back(v0DCApos);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/c } }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0pt, highv0pt).Data();
  latex = CreateLatex(0.3, 0.8, latexText, textSize);

  saveName = TString::Format("%sDCApos", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), lowv0pt, highv0pt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
void plotDCAneg(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  bool setLogY = true;
  double xMinFrame = -10., xMaxFrame = 10., yMinFrame = 1e-5, yMaxFrame = 0.5;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("%s DCA neg", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";
  // string dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = TString::Format("jet-fragmentation/data/jets/V0/jetPt%sPtDCAposneg", hadron.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0ptbins = getProjectionBins(thn->GetAxis(v0ptAxis), v0ptmin, v0ptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(v0ptAxis)->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0pt = thn->GetAxis(v0ptAxis)->GetBinLowEdge(v0ptbins[0]);
  double highv0pt = thn->GetAxis(v0ptAxis)->GetBinUpEdge(v0ptbins[1]);

  TH1D* v0DCAneg = (TH1D*)thn->Projection(dcanegAxis);
  v0DCAneg->SetName(TString::Format("v0DCAneg_jetpt%.0f-%.0f_v0pt%.0f-%.0f", lowjetpt, highjetpt, lowv0pt, highv0pt).Data());
  v0DCAneg->Scale(1./v0DCAneg->Integral());
  setStyle(v0DCAneg, 0);
  histVector.push_back(v0DCAneg);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/c } }",
                              dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0pt, highv0pt).Data();
  TLatex* latex = CreateLatex(0.3, 0.8, latexText, textSize);

  saveName = TString::Format("%sDCAneg", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), lowv0pt, highv0pt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
void plotDCAposneg(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  // string dataSet = "LHC22o_pass6_minBias";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogz = true;
  double xMinFrame = -10., xMaxFrame = 10., yMinFrame = -10, yMaxFrame = 10;
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
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0ptbins = getProjectionBins(thn->GetAxis(v0ptAxis), v0ptmin, v0ptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(v0ptAxis)->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0pt = thn->GetAxis(v0ptAxis)->GetBinLowEdge(v0ptbins[0]);
  double highv0pt = thn->GetAxis(v0ptAxis)->GetBinUpEdge(v0ptbins[1]);

  TH2D* v0DCAposneg = (TH2D*)thn->Projection(dcanegAxis, dcaposAxis);
  v0DCAposneg->SetName(TString::Format("v0DCAposneg_jetpt%.0f-%.0f_v0pt%.0f-%.0f", lowjetpt, highjetpt, lowv0pt, highv0pt).Data());
  v0DCAposneg->Scale(1./v0DCAposneg->Integral());
  setStyle(v0DCAposneg, 0);
  v0DCAposneg->SetMinimum(1e-5);
  v0DCAposneg->SetMaximum(1.);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/c } }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0pt, highv0pt).Data();
  TLatex* latex = CreateLatex(0.2, 0.93, latexText, textSize);

  saveName = TString::Format("%sDCAposneg", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), lowv0pt, highv0pt);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->cd();
  frame->Draw();
  v0DCAposneg->Draw("same colz");
  latex->Draw("same");
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotctau(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if ("V0" == hadron) {
    cout << "Error: hadron must not be V0. To plot V0 ctau, use plotV0ctau" << endl;
    return;
  }
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 40., yMinFrame = 1e-5, yMaxFrame = .1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.2, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  xTitle = TString::Format("#it{c}#tau (%s)", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";
  // string dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = TString::Format("jet-fragmentation/data/jets/V0/jetPt%sPtCtau", formatHadronName(hadron).c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());
  th3->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(th3->GetXaxis(), jetptmin, jetptmax);
  std::array<int, 2> v0ptbins = getProjectionBins(th3->GetYaxis(), v0ptmin, v0ptmax);
  th3->GetXaxis()->SetRange(jetptbins[0], jetptbins[1]);
  th3->GetYaxis()->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowjetpt = th3->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th3->GetXaxis()->GetBinUpEdge(jetptbins[1]);
  double lowv0pt = th3->GetYaxis()->GetBinLowEdge(v0ptbins[0]);
  double highv0pt = th3->GetYaxis()->GetBinUpEdge(v0ptbins[1]);

  TH1D* ctau = (TH1D*)th3->ProjectionZ("v0ctau", jetptbins[0], jetptbins[1], v0ptbins[0], v0ptbins[1]);
  ctau->Scale(1./ctau->Integral());
  ctau->SetName(TString::Format("ctau_jetpt%.0f-%.0f_v0pt%.0f-%.0f", lowjetpt, highjetpt, lowv0pt, highv0pt).Data());
  setStyle(ctau, 0);
  histVector.push_back(ctau);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/c } }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0pt, highv0pt).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("ctau%s", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), lowv0pt, highv0pt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
void plotMass(string inName = "", string dataSet = "dataSet", string hadron = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if ("V0" == hadron) {
    cout << "Error: hadron must not be V0. To plot V0 ctau, use plotV0ctau" << endl;
    return;
  }
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1e-4, yMaxFrame = .1;
  if ("Lambda0" == hadron || "AntiLambda0" == hadron) {
    xMinFrame = 1.05, xMaxFrame = 1.215;
    yMinFrame = 1e-3, yMaxFrame = 0.1;
  }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.2, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  xTitle = TString::Format("#it{M} (%s)", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";
  // string dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = TString::Format("jet-fragmentation/data/jets/V0/jetPt%sPtMass", hadron.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());
  th3->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(th3->GetXaxis(), jetptmin, jetptmax);
  std::array<int, 2> v0ptbins = getProjectionBins(th3->GetYaxis(), v0ptmin, v0ptmax);
  th3->GetXaxis()->SetRange(jetptbins[0], jetptbins[1]);
  th3->GetYaxis()->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowjetpt = th3->GetXaxis()->GetBinLowEdge(jetptbins[0]);
  double highjetpt = th3->GetXaxis()->GetBinUpEdge(jetptbins[1]);
  double lowv0pt = th3->GetYaxis()->GetBinLowEdge(v0ptbins[0]);
  double highv0pt = th3->GetYaxis()->GetBinUpEdge(v0ptbins[1]);

  TH1D* mass = (TH1D*)th3->ProjectionZ("v0mass", jetptbins[0], jetptbins[1], v0ptbins[0], v0ptbins[1]);
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);
  histVector.push_back(mass);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, %s} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax, formatHadronName(hadron).c_str(), v0ptmin, v0ptmax).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%smass", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}

// ----------------------------------------------------------
// Unidentified V0s
// ----------------------------------------------------------
void plotV0ctau(string inName = "", string dataSet = "dataSet", int setting = 2, double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
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
  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 40., yMinFrame = 1e-5, yMaxFrame = .1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.2, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  xTitle = TString::Format("#it{c}#tau (%s)", formatHadronName("K0S").c_str()).Data();
  if (setting == LambdaAxis) { xTitle = TString::Format("#it{c}#tau (%s)", formatHadronName("Lambda0").c_str()).Data(); }
  if (setting == AntiLambdaAxis) { xTitle = TString::Format("#it{c}#tau (%s)", formatHadronName("AntiLambda0").c_str()).Data(); }
  yTitle = "normalised count";
  // string dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtV0PtCtau";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0ptbins = getProjectionBins(thn->GetAxis(v0ptAxis), v0ptmin, v0ptmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(v0ptAxis)->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0pt = thn->GetAxis(v0ptAxis)->GetBinLowEdge(v0ptbins[0]);
  double highv0pt = thn->GetAxis(v0ptAxis)->GetBinUpEdge(v0ptbins[1]);

  TH1D* ctau = (TH1D*)thn->Projection(setting);
  ctau->Scale(1./ctau->Integral());
  ctau->SetName(TString::Format("ctau_jetpt%.0f-%.0f_v0pt%.0f-%.0f", lowjetpt, highjetpt, lowv0pt, highv0pt).Data());
  setStyle(ctau, 0);
  histVector.push_back(ctau);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), lowjetpt, highjetpt, lowv0pt, highv0pt).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "ctauK0S";
  if (setting == LambdaAxis) { saveName = "ctauLambda"; }
  else if (setting == AntiLambdaAxis) { saveName = "ctauAntiLambda"; }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), lowv0pt, highv0pt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
// This plot isn't very useful. ctauL = 2*ctauK0S
void plotV0ctauKL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  // string dataSet = "LHC22o_pass6_minBias";

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { canvas->SetLogz(); }
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

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  saveName = "ctauK0SLambda";
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s_v0pt%.0f-%.0f", saveName.c_str(), v0ptmin, v0ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->cd();
  frame->Draw();
  ctau->Draw("same colz");
  if (latexText != "") { DrawLatex(0.2, 0.93, latexText.c_str(), textSize); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0mass(string inName = "AnalysisResults.root", int setting = 2, double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  // string dataSet = "LHC22o_pass6_minBias";
  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  std::vector<TH1D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
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
  plotNHists(canvas, frame, histVector, legend, saveName, "", latexText);
}
void plotV0massKL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  // string dataSet = "LHC22o_pass6_minBias";
  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { canvas->SetLogz(); }
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

  canvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0massKaL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  // string dataSet = "LHC22o_pass6_minBias";
  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { canvas->SetLogz(); }
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

  canvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotV0massLL(string inName = "AnalysisResults.root", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
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
  // string dataSet = "LHC22o_pass6_minBias";
  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{p}_{T, V0} = %.0f - %.0f GeV/c } }", dataSet.c_str(), jetptmin, jetptmax, v0ptmin, v0ptmax).Data();

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogz) { canvas->SetLogz(); }
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

  canvas->cd();
  frame->Draw();
  mass->Draw("same colz");
  if (latexText != "") { DrawLatex(0.3, 0.3, latexText.c_str(), textSize); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotmassWithCuts(string inName = "AnalysisResults.root", int setting = 3, double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.,
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
  // string dataSet = "LHC22o_pass6_minBias";
  latexText = "";

  std::vector<TH1D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
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
  legend->AddEntry(onlyJet, TString::Format("#it{p}_{T, ch. jet} = %.0f - %.0f GeV/#it{c}", jetptmin, jetptmax).Data());
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
  plotNHists(canvas, frame, histVector, legend, saveName, "", latexText);
}

// ----------------------------------------------------------
// Number of taggedV0s in jets
// ----------------------------------------------------------
void plotNV0sInJet(string inName, string hadron, double jetptmin = 10., double jetptmax = 200.)
{
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0Axis         = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = -0.5, xMaxFrame = 5, yMinFrame = 1e-1, yMaxFrame = 1e7;
  double xMinLegend = 0.4, xMaxLegend = 0.8, yMinLegend = 0.6, yMaxLegend = 0.85;

  int xCanvas = 900, yCanvas = 900;
  xTitle = TString::Format("#it{N} (%s) #in jet", formatHadronName(hadron).c_str()).Data();
  yTitle = "#it{N}_{jets}";
  // string dataSet = "LHC22o_pass6_minBias_small";
  latexText = "";

  std::vector<TH1D*> histVector;

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0/jetPtnV0nK0SnLambdanAntiLambda";
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  int firstBinJetPt = 1, lastBinJetPt = thn->GetAxis(jetptAxis)->GetNbins();
  int firstBinK0SMass = 1, lastBinK0SMass = thn->GetAxis(K0SAxis)->GetNbins();
  int firstBinLambdaMass = 1, lastBinLambdaMass = thn->GetAxis(LambdaAxis)->GetNbins();
  int firstBinAntiLambdaMass = 1, lastBinAntiLambdaMass = thn->GetAxis(AntiLambdaAxis)->GetNbins();

  firstBinJetPt = thn->GetAxis(jetptAxis)->FindBin(jetptmin + 1e-3);
  lastBinJetPt  = thn->GetAxis(jetptAxis)->FindBin(jetptmax - 1e-3);

  TH1D* hJets = (TH1D*) thn->Projection(jetptAxis);
  double nJets = hJets->Integral();

  int projectionAxis;
  array<int, 3> otherAxes = {0, 0, 0};
  array<string, 4> labels = {"No sel", "", "", ""};

  if ("K0S" == hadron) {
    projectionAxis = K0SAxis;
    otherAxes[1] = LambdaAxis;
    otherAxes[2] =  AntiLambdaAxis;
    labels[1] = "#it{N} (#Lambda) > 0";
    labels[2] = "#it{N} (#bar{#Lambda}) > 0";
  }
  if ("Lambda0" == hadron) {
    projectionAxis = LambdaAxis;
    otherAxes[1] = K0SAxis;
    otherAxes[2] =  AntiLambdaAxis;
    labels[1] = "#it{N} (K^{0}_{S}) > 0";
    labels[2] = "#it{N} (#bar{#Lambda}) > 0";
  }
  if ("AntiLambda0" == hadron) {
    projectionAxis = AntiLambdaAxis;
    otherAxes[1] = LambdaAxis;
    otherAxes[2] =  K0SAxis;
    labels[1] = "#it{N} (K^{0}_{S}) > 0";
    labels[2] = "#it{N} (#Lambda) > 0";
  }
  labels[3] = TString::Format("%s, %s", labels[1].c_str(), labels[2].c_str());

  for (int i = 0; i < 4; i++) {
    THnSparseD* tmp = (THnSparseD*) thn->Clone("temporary");
    tmp->GetAxis(jetptAxis)->SetRange(firstBinJetPt, lastBinJetPt);

    if (otherAxes.size() == i) {
      tmp->GetAxis(otherAxes[1])->SetRange(2, thn->GetAxis(otherAxes[1])->GetNbins());
      tmp->GetAxis(otherAxes[2])->SetRange(2, thn->GetAxis(otherAxes[2])->GetNbins());
    }
    else if (0 < i) {
      tmp->GetAxis(otherAxes[i])->SetRange(2, thn->GetAxis(otherAxes[i])->GetNbins());
    }

    TH1D* h = tmp->Projection(projectionAxis);
    h->SetName(TString::Format("h%d", i));
    // h->Scale(1./nJets);
    // h->Scale(1./h->Integral());
    setStyle(h, i);
    legend->AddEntry(h, labels[i].c_str());
    histVector.push_back(h);
  }

  latexText = TString::Format("%s, %.0f < #it{p}_{T, jet} < %.0f", dataSet.c_str(), jetptmin, jetptmax);
  TLatex* latex = CreateLatex(0.2, 0.94, latexText, textSize);

  frame->Draw();
  legend->Draw("same");
  latex->Draw("same");
  for (auto hist : histVector) {
    hist->Draw("same");
    hist->Print();
    for (int i = 1; i < 4; i++) {
      cout << hist->GetBinContent(i) << endl;
    }
  }

  saveName = "nK";
  if ("Lambda0" == hadron) { saveName = "nL"; }
  if ("AntiLambda0" == hadron) { saveName = "nAL"; }
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), jetptmin, jetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  canvas->SaveAs(saveName.c_str());
}

// ----------------------------------------------------------
// Lambda0
// ----------------------------------------------------------
void plotLambda0Pt(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200.)
{ plotPt(inName, "Lambda0", jetptmin, jetptmax); }
void plotLambda0Z(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200.)
{ plotZ(inName, "Lambda0", jetptmin, jetptmax); }

void plotLambda0Radius(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotRadius(inName, "Lambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotLambda0CosPA(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotCosPA(inName, "Lambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotLambda0DCAdaughters(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCAdaughters(inName, "Lambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotLambda0DCApos(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCApos(inName, "Lambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotLambda0DCAneg(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCAneg(inName, "Lambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotLambda0DCAposneg(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCAposneg(inName, "Lambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotLambda0ctau(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotctau(inName, "Lambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotLambda0Mass(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotMass(inName, "Lambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }

// ----------------------------------------------------------
// AntiLambda0
// ----------------------------------------------------------
void plotAntiLambda0Pt(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200.)
{ plotPt(inName, "AntiLambda0", jetptmin, jetptmax); }
void plotAntiLambda0Z(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200.)
{ plotZ(inName, "AntiLambda0", jetptmin, jetptmax); }

void plotAntiLambda0Radius(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotRadius(inName, "AntiLambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotAntiLambda0CosPA(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotCosPA(inName, "AntiLambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotAntiLambda0DCAdaughters(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCAdaughters(inName, "AntiLambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotAntiLambda0DCApos(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCApos(inName, "AntiLambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotAntiLambda0DCAneg(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCAneg(inName, "AntiLambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotAntiLambda0DCAposneg(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCAposneg(inName, "AntiLambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotAntiLambda0ctau(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotctau(inName, "AntiLambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotAntiLambda0Mass(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotMass(inName, "AntiLambda0", jetptmin, jetptmax, v0ptmin, v0ptmax); }

// ----------------------------------------------------------
// K0S
// ----------------------------------------------------------
void plotK0SPt(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200.)
{ plotPt(inName, "K0S", jetptmin, jetptmax); }
void plotK0SZ(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200.)
{ plotZ(inName, "K0S", jetptmin, jetptmax); }

void plotK0SRadius(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotRadius(inName, "K0S", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotK0SCosPA(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotCosPA(inName, "K0S", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotK0SDCAdaughters(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCAdaughters(inName, "K0S", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotK0SDCApos(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCApos(inName, "K0S", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotK0SDCAneg(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCAneg(inName, "K0S", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotK0SDCAposneg(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotDCAposneg(inName, "K0S", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotK0Sctau(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotctau(inName, "K0S", jetptmin, jetptmax, v0ptmin, v0ptmax); }
void plotK0SMass(string inName = "AnalysisResults.root", string dataSet = "", double jetptmin = 10., double jetptmax = 200., double v0ptmin = 0., double v0ptmax = 100.)
{ plotMass(inName, "K0S", jetptmin, jetptmax, v0ptmin, v0ptmax); }

// ----------------------------------------------------------
// ----------------------------------------------------------
// ----------------------------------------------------------

void plot22o(double jetptmin, double jetptmax, double v0ptmin, double v0ptmax, string hadron, int setting)
{
  string inName = "~/cernbox/TrainOutput/252064/AnalysisResults.root";
  string dataSet = "LHC22o_pass6";

  switch(setting) {
    case 0:
      plotRadius(inName, dataSet, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax);
      break;
    case 1:
      plotCosPA(inName, dataSet, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax);
      break;
    case 2:
      plotDCAdaughters(inName, dataSet, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax);
      break;
    case 3:
      plotDCApos(inName, dataSet, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax);
      break;
    case 4:
      plotDCAneg(inName, dataSet, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax);
      break;
    case 5:
      plotctau(inName, dataSet, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax);
      break;
    default:
      cout << "Invalid setting!" << endl;
      break;
  }
}

void plotAll22o(double jetptmin, double jetptmax, string hadron, int setting)
{
  vector<double> pt = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0};
  for (int i = 0; i < pt.size() - 1; i++) {
    if pt[i] > jetptmax { break; }
    plot22o(jetptmin, jetptmax, pt[i], pt[i+1], hadron, setting);
  }
}