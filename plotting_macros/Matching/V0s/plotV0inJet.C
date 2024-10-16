
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

#include "../../histUtils.C"

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

double getHistRMS(TH1* hInput)
{
  TH1* h = (TH1*)hInput->Clone("h");
  h->Scale(1./h->Integral());
  double rms = 0.;
  for (int i = 1; i < h->GetNbinsX(); i++) {
    double binContent = h->GetBinContent(i);
    double binCenter = h->GetBinCenter(i);
    rms += binContent * binCenter * binCenter;
  }
  rms = TMath::Sqrt(rms);
  return rms;
}

double getNjets(TFile* inFile, double jetptmin, double jetptmax)
{
  string histName = "jet-fragmentation/matching/jets/matchPartJetPtEtaPhi";
  TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
  std::array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
  return jets->Integral(jetptbins[0], jetptbins[1], 0, jets->GetNbinsY() + 1, 0, jets->GetNbinsZ() + 1);
}

bool isHistEmptyInRange(TH1* h, int low, int high, double threshold = 1e-10)
{
  double integral = h->Integral(low, high);
  if (integral != integral) // NaN check
    return true;
  else
    return (integral < threshold);
}

// Returns the scale for drawing histogram. Accounts for bin content and error
double getHistScale(TH1* h, bool doError)
{
  if (!doError) { return h->GetBinContent(h->GetMaximumBin()); }

  double scale = -900.;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    double s = be + bc;
    if (s > scale) { scale = s; }
  }
  return scale;
}

// -------------------------------------------------------------------------------------------------

void ptResolution(string inName, string dataSet, double partjetptmin, double partjetptmax, double partv0min, double partv0max)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  const int nDim = 5;
  const int partJetPtAxis = 0;
  const int detJetPtAxis = 1;
  const int partV0PtAxis = 2;
  const int ptRatioAxis = 3;
  const int ptRelDiffAxis = 4;

  string hadron = "V0";
  gStyle->SetNdivisions(505, "xy");
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  string histName = "jet-fragmentation/matching/jets/V0/";
  histName += "partJetPtDetJetPtPartV0PtRatioPtRelDiffPt";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt  = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(jetptbins[1]);

  std::array<int, 2> v0ptbins  = getProjectionBins(thn->GetAxis(partV0PtAxis), partv0min, partv0max);
  thn->GetAxis(partV0PtAxis)->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowv0  = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(v0ptbins[0]);
  double highv0 = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(v0ptbins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  TH1D* v0res = (TH1D*)thn->Projection(ptRelDiffAxis);
  v0res->SetName("v0resolution");
  v0res->Scale(1./v0res->Integral());
  setStyle(v0res, 0);

  double rms = getHistRMS(v0res);
  string rmsText   = TString::Format("RMS: %.2f", rms).Data();
  string ptV0Text  = TString::Format("#it{p}_{T, V0}^{part.} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string ptJetText = TString::Format("#it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c}", lowjetpt, highjetpt).Data();

  double xLatex = 0.4, yLatex = 0.8;
  string latexText = TString::Format("#splitline{%s}{#splitline{%s}{#splitline{%s}{%s}}}", dataSet.c_str(), ptJetText.c_str(), ptV0Text.c_str(), rmsText.c_str()).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  string saveName = "V0PtResolution";
  saveName += TString::Format("_jetpt%.0f-%.0f", lowjetpt, highjetpt);
  saveName += TString::Format("_v0pt%.1f-%.1f", lowv0, highv0);
  saveName += ".pdf";
  int xCanvas = 900, yCanvas = 900;
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);
  canvas->SetLogy();

  double xMinFrame = -1., xMaxFrame = 10., yMinFrame = 1e-7, yMaxFrame = 2.;
  string xTitle = "(#it{p}_{T, V0}^{det.} - #it{p}_{T, V0}^{part.})/#it{p}_{T, V0}^{part.}";
  string yTitle = "normalised count";
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);

  canvas->cd();
  frame->Draw();
  v0res->Draw("same");
  latex->Draw("same");
  canvas->SaveAs(canvas->GetName());
}
void dauPtResolution(string inName, string dataSet, double partjetptmin, double partjetptmax, double partmin, double partmax, bool doPos)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  const int nDim          = 6;
  const int partJetPtAxis = 0;
  const int detJetPtAxis  = 1;
  const int partV0PtAxis  = 2;
  const int partDauPtAxis = 3;
  const int ptRatioAxis   = 4;
  const int ptRelDiffAxis = 5;

  gStyle->SetNdivisions(505, "xy");
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  string histName = "jet-fragmentation/matching/jets/V0/";
  if (doPos) histName += "partJetPtDetJetPtPartV0PtPosPtRatioPtRelDiffPt";
  else histName += "partJetPtDetJetPtPartV0PtNegPtRatioPtRelDiffPt";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt  = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(jetptbins[1]);

  std::array<int, 2> v0ptbins  = getProjectionBins(thn->GetAxis(partDauPtAxis), partmin, partmax);
  thn->GetAxis(partV0PtAxis)->SetRange(v0ptbins[0], v0ptbins[1]);
  double lowv0  = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(v0ptbins[0]);
  double highv0 = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(v0ptbins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  TH1D* daures = (TH1D*)thn->Projection(ptRelDiffAxis);
  daures->SetName("dauresolution");
  daures->Scale(1./daures->Integral());
  setStyle(daures, 0);

  double rms = getHistRMS(daures);
  string rmsText   = TString::Format("RMS: %.2f", rms).Data();
  string ptV0Text  = TString::Format("#it{p}_{T, %s dau}^{part.} = %.1f - %.1f GeV/#it{c}", (doPos ? "pos" : "neg"), lowv0, highv0).Data();
  string ptJetText = TString::Format("#it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c}", lowjetpt, highjetpt).Data();

  double xLatex = 0.4, yLatex = 0.8;
  string latexText = TString::Format("#splitline{%s}{#splitline{%s}{#splitline{%s}{%s}}}", dataSet.c_str(), ptJetText.c_str(), ptV0Text.c_str(), rmsText.c_str()).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  string saveName = "V0";
  saveName += (doPos ? "Pos" : "Neg");
  saveName += "PtResolution";
  saveName += TString::Format("_jetpt%.0f-%.0f", lowjetpt, highjetpt);
  saveName += TString::Format("_pt%.1f-%.1f", lowv0, highv0);
  saveName += ".pdf";
  int xCanvas = 900, yCanvas = 900;
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);
  canvas->SetLogy();

  double xMinFrame = -1., xMaxFrame = 10., yMinFrame = 1e-7, yMaxFrame = 2.;
  string xTitle = "(#it{p}_{T}^{det.} - #it{p}_{T}^{part.})/#it{p}_{T}^{part.}";
  string yTitle = "normalised count";
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);

  canvas->cd();
  frame->Draw();
  daures->Draw("same");
  latex->Draw("same");
  canvas->SaveAs(canvas->GetName());
}

void matchedPtZ(string inName = "", double partjetptmin = 10., double partjetptmax = 1e6, bool detector = false, bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
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

  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 1e-8, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.35, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, V0}^{%s.} (GeV/#it{c})", detector ? "det" : "part").Data();
  if (doZ) {
    xTitle = TString::Format("#it{z}_{V0}^{%s.}", detector ? "det" : "part").Data();
    xMinFrame = 1e-3; xMaxFrame = 1.+1e-3;
  }
  yTitle = "normalised count";
  dataSet = "LHC24b1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "partJetPtV0PtDetJetPtV0Pt";
  if (doZ) { histName = "matchDetJetPtV0TrackProjPartJetPtV0TrackProj"; }
  histName = TString::Format("jet-fragmentation_id10235/matching/jets/V0/%s", histName.c_str());
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH1D* matchedjetpt = (TH1D*)thn->Projection(detector ? detV0PtAxis : partV0PtAxis);
  matchedjetpt->SetName("matchedjetpt");
  matchedjetpt->Scale(1./ matchedjetpt->Integral());
  setStyle(matchedjetpt, 0);
  histVector.push_back(matchedjetpt);

  latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f (GeV/#it{c}) }", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("matched%sV0%s", detector ? "Det" : "Part" , doZ ? "Z" : "Pt").Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}

void matchedCtau(string inName = "", string hadron = "", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200., double partv0min = -1., double partv0max = 1e6, bool doZ = false)
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

  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 1e-6, yMaxFrame = 1.;
  // if (doZ) { yMinFrame = 1e-3; yMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{c}#tau_{%s} (cm)", formatHadronName(hypothesis).c_str()).Data();
  // yTitle = TString::Format("#it{p}_{T, %s}^{part.} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  // if (doZ) { yTitle = TString::Format("#it{z}_{%s}^{part.}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = "normalised count";
  dataSet = "LHC24b1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("partJetPt%sPtDetJetPt%sPtCtau%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data();
  if (doZ) { histName = TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProjCtau%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data(); }
  // histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  histName = TString::Format("jet-fragmentation_id10235/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  std::array<int, 2> partv0bins = getProjectionBins(thn->GetAxis(partV0PtAxis), partv0min, partv0max);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  thn->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);
  double lowv0 = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(partv0bins[0]);
  double highv0 = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(partv0bins[1]);

  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  TH1D* v0ctau = (TH1D*)thn->Projection(detV0CtauAxis);
  v0ctau->Scale(1./v0ctau->Integral());
  // normaliseHistRowByRow(v0ctau);
  // v0ctau->SetMinimum(1e-5);
  // v0ctau->SetMaximum(1.);
  v0ctau->SetName(TString::Format("%sCtau%s", hadron.c_str(), hypothesis.c_str()).Data());
  // if (doZ) { v0ctau->SetName(TString::Format("%sZCtau%s", hadron.c_str(), hypothesis.c_str()).Data()); }
  histVector.push_back(v0ctau);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) {
    latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.0f - %.0f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  }
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("matched%sCtau%s", hadron.c_str(), hypothesis.c_str()).Data();
  // if (doZ) { saveName = TString::Format("matched%sZCtau%s", hadron.c_str(), hypothesis.c_str()).Data(); }
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  if (doZ) { saveName = TString::Format("%s_partv0z%.1f-%.1f", saveName.c_str(), lowv0, highv0); }
  else { saveName = TString::Format("%s_partv0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void matchedMass(string inName = "", string hadron = "", string hypothesis = "", double partjetptmin = 10., double partjetptmax = 200., double partv0min = -1., double partv0max = 1e6, bool doZ = false)
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

  bool setLogY = true;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 1e-7, yMaxFrame = 0.1;
  if ( "K0S" == hypothesis ) {
    xMinFrame = 0.45; xMaxFrame = 0.55;
  }
  // else if ("K0S" == hadron) { yMaxFrame = 15e-3; }
  // if (doZ) { yMinFrame = 1e-3; yMaxFrame = 1. + 1e-3; }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronName(hypothesis).c_str()).Data();
  // yTitle = TString::Format("#it{p}_{T, %s}^{part.} (GeV/#it{c})", formatHadronName(hadron).c_str()).Data();
  // if (doZ) { yTitle = TString::Format("#it{z}_{%s}^{part.}", formatHadronName(hadron).c_str()).Data(); }
  yTitle = "normalised count";
  dataSet = "LHC24b1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("partJetPt%sPtDetJetPt%sPtMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data();
  if (doZ) { histName = TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProjMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation_id10235/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  std::array<int, 2> partv0bins = getProjectionBins(thn->GetAxis(partV0PtAxis), partv0min, partv0max);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  thn->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);
  double lowv0 = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(partv0bins[0]);
  double highv0 = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(partv0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  TH1D* v0mass = (TH1D*)thn->Projection(detV0MassAxis);
  // normaliseHistRowByRow(v0mass);
  // v0mass->SetMinimum(1e-5);
  // v0mass->SetMaximum(1.);
  v0mass->Scale(1./v0mass->Integral());
  setStyle(v0mass, 0);
  v0mass->SetName(TString::Format("%sMass%s", hadron.c_str(), hypothesis.c_str()).Data());
  // if (doZ) { v0mass->SetName(TString::Format("%sTrackProjMass%s", hadron.c_str(), hypothesis.c_str()).Data()); }
  histVector.push_back(v0mass);

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) { latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.0f - %.0f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data(); }
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("matched%sMass%s", hadron.c_str(), hypothesis.c_str()).Data();
  // if (doZ) { saveName = TString::Format("matched%sZMass%s", hadron.c_str(), hypothesis.c_str()).Data(); }
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  if (doZ) {saveName = TString::Format("%s_partv0z%.1f-%.1f", saveName.c_str(), lowv0, highv0); }
  else { saveName = TString::Format("%s_partv0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

void mass(vector<string> inputStrings, double partjetptmin, double partjetptmax, double v0min, double v0max, bool doZ, bool normalise)
{
  gStyle->SetNdivisions(505, "xy");
  string inName     = inputStrings[0];
  string dataSet    = inputStrings[1];
  string hadron     = inputStrings[2];
  string hypothesis = inputStrings[3];

  const int nDim          = 5;
  const int partJetPtAxis = 0;
  const int partV0PtAxis  = 1;
  const int detJetPtAxis  = 2;
  const int detV0PtAxis   = 3;
  const int massAxis      = 4;

  double textSize = 0.04;

  string histName = "jet-fragmentation/matching/jets/V0";
  if (doZ) histName += TString::Format("/partJetPt%sTrackProjDetJetPt%sTrackProjMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str());
  else histName += TString::Format("/partJetPt%sPtDetJetPt%sPtMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str());
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt  = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);
  array<int, 2> partv0bins = getProjectionBins(thn->GetAxis(partV0PtAxis), v0min, v0max, 2e-4);
  thn->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);
  double lowv0  = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(partv0bins[0]);
  double highv0 = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(partv0bins[1]);

  string saveName = hadron;
  saveName += "_mass";
  saveName += hypothesis;
  saveName += TString::Format("_partjetpt%.0f-%.0f", lowjetpt, highjetpt).Data();
  if (doZ) saveName += TString::Format("_partv0z%.2f-%.2f", lowv0, highv0).Data();
  else saveName += TString::Format("_partv0pt%.1f-%.1f", lowv0, highv0).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  TH1D* mass = (TH1D*)thn->Projection(massAxis);
  if (isHistEmptyInRange(mass, 1, mass->GetNbinsX())) {
    cout << "Error: Empty mass histogram for jetpt " << lowjetpt << " - " << highjetpt << ", " << (doZ ? "zv0 " : "ptv0 ") << lowv0 << " - " << highv0 << endl;
    return;
  }
  // mass->Rebin(4);
  if (normalise) mass->Scale(1./mass->Integral(), "width");
  setStyle(mass, 0);
  mass->SetName(saveName.c_str());

  string xTitle = "M(" + formatHadronDaughters(hypothesis) + ") (GeV/#it{c}^{2})";
  string yTitle = "counts";
  if (normalise) yTitle = "#frac{1}{#it{N}_{" + formatHadronName(hadron) + "}} #frac{d#it{N_{" + formatHadronName(hadron) + "}}}{d#it{M}}";
  double xMinFrame = mass->GetXaxis()->GetXmin(), xMaxFrame = mass->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 2.0 * getHistScale(mass, true);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  vector<TObject*> objList;
  string jetText = TString::Format("#it{p}_{T, ch+V0 jet}^{part.} = %.0f-%.0f GeV/#it{c}", lowjetpt, highjetpt).Data();
  string v0Text  = TString::Format("#it{p}_{T, %s} = %.1f-%.1f GeV/#it{c}", hadron.c_str(), lowv0, highv0).Data();
  if (doZ) { v0Text = TString::Format("#it{z}_{%s} = %.2f-%.2f", hadron.c_str(), lowv0, highv0).Data(); }
  string latexText = TString::Format("#splitline{%s}{%s}", jetText.c_str(), v0Text.c_str()).Data();
  double xLatex = 0.3, yLatex = 0.8;
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  objList.push_back(latex);

  double v0mass = ("K0S" == hypothesis) ? MassK0S : MassLambda0;
  TLine* line = new TLine(v0mass, yMinFrame, v0mass, 0.75 * yMaxFrame);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(1));
  line->SetLineWidth(3);
  objList.push_back(line);

  canvas->cd();
  frame->Draw();
  for (auto obj : objList) { obj->Draw("same"); }
  mass->Draw("same");
  canvas->SaveAs(canvas->GetName());
}

// ----------------------------------------------------------
// Number of taggedV0s in jets
// ----------------------------------------------------------
void plotNV0sInJet(vector<string> inputStrings, double jetptmin, double jetptmax, bool doRatio)
{
  string inName  = inputStrings[0];
  string hadron  = inputStrings[1];
  string dataSet = inputStrings[2];
  const int nDim           = 5;
  const int jetptAxis      = 0;
  const int v0Axis         = 1;
  const int K0SAxis        = 2;
  const int LambdaAxis     = 3;
  const int AntiLambdaAxis = 4;
  gStyle->SetNdivisions(500, "x");
  gStyle->SetNdivisions(505, "y");
  double textSize  = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  string histName = "jet-fragmentation/matching/jets/V0/jetPtnV0MatchednK0SnLambdanAntiLambda";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

  TH1D* hJets = (TH1D*) thn->Projection(jetptAxis);
  double nJets = hJets->Integral(jetptbins[0], jetptbins[1]);

  int projectionAxis;
  array<int, 3> otherAxes = {0, 0, 0};
  array<string, 4> labels = {"No sel", "", "", ""};

  if ("K0S" == hadron) {
    projectionAxis = K0SAxis;
    otherAxes[1] = LambdaAxis;
    otherAxes[2] = AntiLambdaAxis;
    labels[1] = "#it{N}(#Lambda)>0";
    labels[2] = "#it{N}(#bar{#Lambda})>0";
  }
  if ("Lambda0" == hadron) {
    projectionAxis = LambdaAxis;
    otherAxes[1] = K0SAxis;
    otherAxes[2] = AntiLambdaAxis;
    labels[1] = "#it{N}(K^{0}_{S})>0";
    labels[2] = "#it{N}(#bar{#Lambda})>0";
  }
  if ("AntiLambda0" == hadron) {
    projectionAxis = AntiLambdaAxis;
    otherAxes[1] = K0SAxis;
    otherAxes[2] = LambdaAxis;
    labels[1] = "#it{N}(K^{0}_{S})>0";
    labels[2] = "#it{N}(#Lambda)>0";
  }
  labels[3] = TString::Format("%s, %s", labels[1].c_str(), labels[2].c_str());

  double xMinLegend = 0.4, xMaxLegend = 0.79, yMinLegend = 0.6, yMaxLegend = 0.85;
  string legendTitle;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  std::vector<TH1D*> histVector;

  for (int i = 0; i < 4; i++) {
    THnSparseD* tmp = (THnSparseD*) thn->Clone("temporary");
    tmp->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);

    if (otherAxes.size() == i) {
      // Require 1 of both other hadrons
      tmp->GetAxis(otherAxes[1])->SetRange(2, thn->GetAxis(otherAxes[1])->GetNbins());
      tmp->GetAxis(otherAxes[2])->SetRange(2, thn->GetAxis(otherAxes[2])->GetNbins());
    }
    else if (0 < i) {
      // Require 1 of the other hadron
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

  string xTitle = TString::Format("#it{N}(%s)#in jet", formatHadronName(hadron).c_str()).Data();
  string yTitle = "#it{N}_{jets}";
  double xMinFrame = -0.5, xMaxFrame = 9.5, yMinFrame = 1e-1, yMaxFrame = 1e7;

  // Scale hists with content of first bin, to get some measure of relative probability (?)
  if (doRatio) {
    double firstBin = 0.;
    for (auto hist : histVector) {
      firstBin = hist->GetBinContent(1);
      hist->Scale(1./firstBin);
    }
    yTitle = "#it{N}/#it{N}_{0}";
    yMinFrame = 1e-3, yMaxFrame = 10.;
  }

  int xCanvas = 900, yCanvas = 900;
  string saveName = "nK";
  if ("Lambda0" == hadron) { saveName = "nL"; }
  if ("AntiLambda0" == hadron) { saveName = "nAL"; }
  saveName += TString::Format("_jetpt%.0f-%.0f", jetptmin, jetptmax);
  if (doRatio) { saveName += "_ratio"; }
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);
  canvas->SetLogy();

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  string ptJetText = TString::Format("%.0f#leq#it{p}_{T, ch+V0 jet} < %.0f GeV/#it{c}", lowjetpt, highjetpt).Data();
  string latexText = dataSet + ", " + ptJetText;
  frame->SetTitle(latexText.c_str());
  frame->Draw();
  legend->Draw("same");
  for (auto hist : histVector) {
    hist->Draw("same");
    hist->Print();
    // Print hist content
    for (int i = 1; i < hist->GetNbinsX(); i++) {
      double binContent = hist->GetBinContent(i);
      if (binContent < 1e-15) {
        break;
      }
      cout << hist->GetBinContent(i) << endl;
    }
  }
  canvas->SaveAs(canvas->GetName());
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

void plotTrain(string train, string dataSet, string hadron, double jetptmin, double jetptmax, double v0ptmin, double v0ptmax, bool doZ, int setting)
{
  string inputName = "~/cernbox/TrainOutput/" + train + "/AnalysisResults.root";
  switch (setting) {
    case 0:
      ptResolution(inputName, dataSet, jetptmin, jetptmax, v0ptmin, v0ptmax);
      break;
    case 1:
      {
        dauPtResolution(inputName, dataSet, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ);
      }
      break;
    case 2:
      {
        vector<string> inputStrings = {inputName, hadron, dataSet};
        plotNV0sInJet(inputStrings, jetptmin, jetptmax, false /* doRatio */);
        plotNV0sInJet(inputStrings, jetptmin, jetptmax, true);
      }
      break;
    case 3:
      {
        vector<string> inputStrings = {inputName, dataSet, hadron, "K0S"};
        mass(inputStrings, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, false /* normalise */);
        mass(inputStrings, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, true);
      }
      break;
    case 4:
      {
        vector<string> inputStrings = {inputName, dataSet, hadron, "Lambda0"};
        mass(inputStrings, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, false /* normalise */);
        mass(inputStrings, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, true);
      }
      break;
    case 5:
      {
        vector<string> inputStrings = {inputName, dataSet, hadron, "AntiLambda0"};
        mass(inputStrings, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, false /* normalise */);
        mass(inputStrings, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, true);
      }
      break;
    default:
      cout << "Error: invalid setting" << endl;
      return;
  }
}
void plot271952(string hadron, double jetptmin, double jetptmax, double v0ptmin, double v0ptmax, bool doZ, int setting)
{
  string train = "271952";
  string dataSet = "LHC24b1b";
  plotTrain(train, dataSet, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, setting);
}
void plot271952(string hadron, double jetptmin, double jetptmax, bool doZ, int setting)
{
  gROOT->SetBatch();
  vector<double> pt = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};
  vector<double> z  = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  vector<double> bins = doZ ? z : pt;
  for (int i = 0; i < bins.size()-1; i++) {
    if (bins[i] > jetptmin) break;
    plot271952(hadron, jetptmin, jetptmax, bins[i], bins[i+1], doZ, setting);
  }
}

void plot210373(string hadron, double jetptmin, double jetptmax, double v0ptmin, double v0ptmax, bool doZ, int setting)
{
  string train = "210373";
  string dataSet = "LHC24b1b";
  plotTrain(train, dataSet, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, setting);
}
