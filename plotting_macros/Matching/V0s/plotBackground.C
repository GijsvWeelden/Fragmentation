
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

/*
 * Plots background from mistagged V0s
 * E.g. K0SMass("Lambda0") plots and fits the shape of mK0S for true Lambdas
 */

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

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
bool isHistEmptyInRange(TH1* h, double low, double high, double threshold = 1e-10)
{
  array<int, 2> bins = getProjectionBins(h->GetXaxis(), low, high);
  return isHistEmptyInRange(h, bins[0], bins[1], threshold);
}

// -------------------------------------------------------------------------------------------------

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
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) {
    latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.0f - %.0f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
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
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) { latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.0f - %.0f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data(); }
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("matched%sMass%s", hadron.c_str(), hypothesis.c_str()).Data();
  // if (doZ) { saveName = TString::Format("matched%sZMass%s", hadron.c_str(), hypothesis.c_str()).Data(); }
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  if (doZ) {saveName = TString::Format("%s_partv0z%.1f-%.1f", saveName.c_str(), lowv0, highv0); }
  else { saveName = TString::Format("%s_partv0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

// Mass under K0S hypothesis from true Lambda/antiLambda
void K0SMass(string inName = "", string hadron = "", double partjetptmin = 10., double partjetptmax = 200., double partv0min = -1., double partv0max = 1e6, bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either Lambda0, or AntiLambda0" << endl;
    return;
  }
  string hypothesis = "K0S";
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

  bool setLogY = false;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 0., yMaxFrame = 0.1;
  double xMinLegend = 0.6, xMaxLegend = 0.9, yMinLegend = 0.4, yMaxLegend = 0.6;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 10;
  xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronName(hypothesis).c_str()).Data();
  yTitle = "normalised count";
  dataSet = "LHC24b1";

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("partJetPt%sPtDetJetPt%sPtMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data();
  if (doZ) { histName = TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProjMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  // histName = TString::Format("jet-fragmentation_id10235/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
  v0mass->Rebin(rebinNumber);
  v0mass->Scale(1./v0mass->Integral());
  setStyle(v0mass, 0);
  v0mass->SetName(TString::Format("%sMass%s", hadron.c_str(), hypothesis.c_str()).Data());
  legend->AddEntry(v0mass, "data");

  canvas->cd();
  frame->Draw();
  v0mass->Draw("same");

  double parameters[3];
  double fitRegion[2] = {0.42, 0.55};
  TF1* background = new TF1("background", "pol2", fitRegion[0], fitRegion[1]);
  TF1* background_extrapolated = new TF1("background_extrapolated", "pol2", 0.4, 0.6);

  setStyle(background, 1);
  TFitResultPtr bkgPtr = v0mass->Fit(background, "S R");
  background->GetParameters(parameters);
  legend->AddEntry(background, "pol2 fit");

  setStyle(background_extrapolated, 1);
  background_extrapolated->SetLineStyle(9);
  background_extrapolated->SetParameters(parameters);
  background_extrapolated->Draw("same");

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) { latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.2f - %.2f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data(); }
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  TLatex* latexChi2 = CreateLatex(0.3, 0.6, TString::Format("#chi^{2}/NDF = %.2f / %d = %.2f", bkgPtr->Chi2(), bkgPtr->Ndf(), bkgPtr->Chi2() / bkgPtr->Ndf()).Data(), textSize);

  TLine* line = new TLine(MassK0S, yMinFrame, MassK0S, 0.05);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(0));
  line->SetLineWidth(2);

  legend->Draw("same");
  latex->Draw("same");
  latexChi2->Draw("same");
  line->Draw("same");

  saveName = TString::Format("%sMass%s", hadron.c_str(), hypothesis.c_str()).Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  if (doZ) {saveName = TString::Format("%s_partv0z%.1f-%.1f", saveName.c_str(), lowv0, highv0); }
  else { saveName = TString::Format("%s_partv0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->SaveAs(saveName.c_str());
}
// Mass under Lambda hypothesis from true K0S/antiLambda
void Lambda0Mass(string inName = "", string hadron = "", double partjetptmin = 10., double partjetptmax = 200., double partv0min = -1., double partv0max = 1e6, bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, or AntiLambda0" << endl;
    return;
  }
  string hypothesis = "Lambda0";
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

  bool setLogY = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 0.07;
  double xMinLegend = 0.25, xMaxLegend = 0.9, yMinLegend = 0.2, yMaxLegend = 0.35;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronName(hypothesis).c_str()).Data();
  yTitle = "normalised count";
  dataSet = "LHC24b1";

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("partJetPt%sPtDetJetPt%sPtMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data();
  if (doZ) { histName = TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProjMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  // histName = TString::Format("jet-fragmentation_id10235/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
  // v0mass->Rebin(rebinNumber);
  v0mass->Scale(1./v0mass->Integral());
  setStyle(v0mass, 0);
  v0mass->SetName(TString::Format("%sMass%s", hadron.c_str(), hypothesis.c_str()).Data());
  legend->AddEntry(v0mass, "data");

  canvas->cd();
  frame->Draw();
  v0mass->Draw("same");

  double parameters_pol[2];
  double parameters_pol2[3];
  double parameters_exp[2];
  double fitRegion[4] = {1.1, 1.15, 1.17, 1.215};
  TF1* pol = new TF1("pol", "pol1", fitRegion[0], fitRegion[1]);
  TF1* pol_extrapolated = new TF1("pol_extrapolated", "pol1", 1.05, 1.215);
  TF1* pol2 = new TF1("pol2", "pol2", fitRegion[0], fitRegion[1]);
  TF1* pol2_extrapolated = new TF1("pol2_extrapolated", "pol2", 1.05, 1.215);
  TF1* exp1 = new TF1("exp1", "expo", fitRegion[0], fitRegion[1]);
  TF1* exp1_extrapolated = new TF1("exp1_extrapolated", "expo", 1.05, 1.215);
  TF1* exp2 = new TF1("exp2", "expo", fitRegion[2], fitRegion[3]);
  TF1* total = new TF1("total", "pol1(0) + expo(2)", fitRegion[0], fitRegion[3]);

  setStyle(pol, 1);
  TFitResultPtr polptr = v0mass->Fit(pol, "S R");
  pol->GetParameters(&parameters_pol[0]);
  legend->AddEntry(pol, "pol1 fit");

  setStyle(pol_extrapolated, 1);
  pol_extrapolated->SetLineStyle(9);
  pol_extrapolated->SetParameters(parameters_pol);
  pol_extrapolated->Draw("same");

  setStyle(exp1, 2);
  TFitResultPtr expptr = v0mass->Fit(exp1, "S R+");
  exp1->GetParameters(&parameters_exp[0]);
  legend->AddEntry(exp1, "expo fit");

  setStyle(exp1_extrapolated, 2);
  exp1_extrapolated->SetLineStyle(9);
  exp1_extrapolated->SetParameters(parameters_exp);
  exp1_extrapolated->Draw("same");

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) { latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.2f - %.2f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data(); }
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  TLatex* latexPolChi2 = CreateLatex(0.25, 0.6, TString::Format("Pol: #chi^{2}/NDF = %.2f / %d = %.2f", polptr->Chi2(), polptr->Ndf(), polptr->Chi2() / polptr->Ndf()).Data(), textSize);
  TLatex* latexPol = CreateLatex(0.25, 0.55, TString::Format("Pol = %.2f + %.2f M", parameters_pol[0], parameters_pol[1]).Data(), textSize);
  TLatex* latexExpChi2 = CreateLatex(0.25, 0.5, TString::Format("Exp: #chi^{2}/NDF = %.2f / %d = %.2f", expptr->Chi2(), expptr->Ndf(), expptr->Chi2() / expptr->Ndf()).Data(), textSize);
  TLatex* latexExp = CreateLatex(0.25, 0.45, TString::Format("Exp = exp(%.2f + %.2f M)", parameters_exp[0], parameters_exp[1]).Data(), textSize);

  TLine* line = new TLine(MassLambda0, yMinFrame, MassLambda0, 0.05);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(0));
  line->SetLineWidth(2);

  legend->Draw("same");
  latex->Draw("same");
  latexPolChi2->Draw("same");
  latexPol->Draw("same");
  latexExpChi2->Draw("same");
  latexExp->Draw("same");
  if ("K0S" == hadron) {
    setStyle(pol2, 3);
    TFitResultPtr pol2ptr = v0mass->Fit(pol2, "S R+");
    pol2->GetParameters(&parameters_pol2[0]);
    legend->AddEntry(pol2, "pol2 fit");

    setStyle(pol2_extrapolated, 3);
    pol2_extrapolated->SetLineStyle(9);
    pol2_extrapolated->SetParameters(parameters_pol2);
    pol2_extrapolated->Draw("same");

    TLatex* latexPol2Chi2 = CreateLatex(0.25, 0.4, TString::Format("Pol2: #chi^{2}/NDF = %.2f / %d = %.2f", pol2ptr->Chi2(), pol2ptr->Ndf(), pol2ptr->Chi2() / pol2ptr->Ndf()).Data(), textSize);
    TLatex* latexPol2 = CreateLatex(0.25, 0.35, TString::Format("Pol2 = %.2f + %.2f M + %.2f M^{2}", parameters_pol2[0], parameters_pol2[1], parameters_pol2[2]).Data(), textSize);

    latexPol2Chi2->Draw("same");
    latexPol2->Draw("same");
  }
  line->Draw("same");

  saveName = TString::Format("%sMass%s", hadron.c_str(), hypothesis.c_str()).Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  if (doZ) {saveName = TString::Format("%s_partv0z%.1f-%.1f", saveName.c_str(), lowv0, highv0); }
  else { saveName = TString::Format("%s_partv0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->SaveAs(saveName.c_str());
}
// Mass under antiLambda hypothesis from true Lambda/K0S
void AntiLambda0Mass(string inName = "", string hadron = "", double partjetptmin = 10., double partjetptmax = 200., double partv0min = -1., double partv0max = 1e6, bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, or Lambda0" << endl;
    return;
  }
  string hypothesis = "AntiLambda0";
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

  bool setLogY = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 0.1;
  double xMinLegend = 0.25, xMaxLegend = 0.9, yMinLegend = 0.2, yMaxLegend = 0.35;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronName(hypothesis).c_str()).Data();
  yTitle = "normalised count";
  dataSet = "LHC24b1";

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("partJetPt%sPtDetJetPt%sPtMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data();
  if (doZ) { histName = TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProjMass%s", hadron.c_str(), hadron.c_str(), hypothesis.c_str()).Data(); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  // histName = TString::Format("jet-fragmentation_id10235/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
  v0mass->Scale(1./v0mass->Integral());
  setStyle(v0mass, 0);
  v0mass->SetName(TString::Format("%sMass%s", hadron.c_str(), hypothesis.c_str()).Data());
  legend->AddEntry(v0mass, "data");

  canvas->cd();
  frame->Draw();
  v0mass->Draw("same");

  double parameters_pol[2];
  double parameters_pol2[3];
  double parameters_exp[2];
  double fitRegion[4] = {1.1, 1.15, 1.17, 1.215};
  TF1* pol = new TF1("pol", "pol1", fitRegion[0], fitRegion[1]);
  TF1* pol_extrapolated = new TF1("pol_extrapolated", "pol1", 1.05, 1.215);
  TF1* pol2 = new TF1("pol2", "pol2", fitRegion[0], fitRegion[1]);
  TF1* pol2_extrapolated = new TF1("pol2_extrapolated", "pol2", 1.05, 1.215);
  TF1* exp1 = new TF1("exp1", "expo", fitRegion[0], fitRegion[1]);
  TF1* exp1_extrapolated = new TF1("exp1_extrapolated", "expo", 1.05, 1.215);
  TF1* exp2 = new TF1("exp2", "expo", fitRegion[2], fitRegion[3]);
  TF1* total = new TF1("total", "pol1(0) + expo(2)", fitRegion[0], fitRegion[3]);

  setStyle(pol, 1);
  TFitResultPtr polptr = v0mass->Fit(pol, "S R");
  pol->GetParameters(&parameters_pol[0]);
  legend->AddEntry(pol, "pol1 fit");

  setStyle(pol_extrapolated, 1);
  pol_extrapolated->SetLineStyle(9);
  pol_extrapolated->SetParameters(parameters_pol);
  pol_extrapolated->Draw("same");

  setStyle(exp1, 2);
  TFitResultPtr expptr = v0mass->Fit(exp1, "S R+");
  exp1->GetParameters(&parameters_exp[0]);
  legend->AddEntry(exp1, "expo fit");

  setStyle(exp1_extrapolated, 2);
  exp1_extrapolated->SetLineStyle(9);
  exp1_extrapolated->SetParameters(parameters_exp);
  exp1_extrapolated->Draw("same");

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) { latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.2f - %.2f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data(); }
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  TLatex* latexPolChi2 = CreateLatex(0.25, 0.6, TString::Format("Pol: #chi^{2}/NDF = %.2f / %d = %.2f", polptr->Chi2(), polptr->Ndf(), polptr->Chi2() / polptr->Ndf()).Data(), textSize);
  TLatex* latexPol = CreateLatex(0.25, 0.55, TString::Format("Pol = %.2f + %.2f M", parameters_pol[0], parameters_pol[1]).Data(), textSize);
  TLatex* latexExpChi2 = CreateLatex(0.25, 0.5, TString::Format("Exp: #chi^{2}/NDF = %.2f / %d = %.2f", expptr->Chi2(), expptr->Ndf(), expptr->Chi2() / expptr->Ndf()).Data(), textSize);
  TLatex* latexExp = CreateLatex(0.25, 0.45, TString::Format("Exp = exp(%.2f + %.2f M)", parameters_exp[0], parameters_exp[1]).Data(), textSize);

  TLine* line = new TLine(MassLambda0, yMinFrame, MassLambda0, 0.02);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(0));
  line->SetLineWidth(2);

  legend->Draw("same");
  latex->Draw("same");
  latexPolChi2->Draw("same");
  latexPol->Draw("same");
  latexExpChi2->Draw("same");
  latexExp->Draw("same");
  if ("K0S" == hadron) {
    setStyle(pol2, 3);
    TFitResultPtr pol2ptr = v0mass->Fit(pol2, "S R+");
    pol2->GetParameters(&parameters_pol2[0]);
    legend->AddEntry(pol2, "pol2 fit");

    setStyle(pol2_extrapolated, 3);
    pol2_extrapolated->SetLineStyle(9);
    pol2_extrapolated->SetParameters(parameters_pol2);
    pol2_extrapolated->Draw("same");

    TLatex* latexPol2Chi2 = CreateLatex(0.25, 0.4, TString::Format("Pol2: #chi^{2}/NDF = %.2f / %d = %.2f", pol2ptr->Chi2(), pol2ptr->Ndf(), pol2ptr->Chi2() / pol2ptr->Ndf()).Data(), textSize);
    TLatex* latexPol2 = CreateLatex(0.25, 0.35, TString::Format("Pol2 = %.2f + %.2f M + %.2f M^{2}", parameters_pol2[0], parameters_pol2[1], parameters_pol2[2]).Data(), textSize);

    latexPol2Chi2->Draw("same");
    latexPol2->Draw("same");
  }
  line->Draw("same");

  saveName = TString::Format("%sMass%s", hadron.c_str(), hypothesis.c_str()).Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  if (doZ) {saveName = TString::Format("%s_partv0z%.1f-%.1f", saveName.c_str(), lowv0, highv0); }
  else { saveName = TString::Format("%s_partv0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  saveName = TString::Format("%s.pdf", saveName.c_str());

  canvas->SaveAs(saveName.c_str());
}
// -------------------------------------------------------------------------------------------------

// Inv mass of combinatorial V0s
void combMass(vector<string> inputStrings, double v0min, double v0max, double dM /* in MeV */)
{
  gStyle->SetNdivisions(505, "xy");
  string inName  = inputStrings[0];
  string hadron  = inputStrings[1];
  string dataSet = inputStrings[2];

  const int nDim            = 4;
  const int ptAxis          = 0;
  const int K0SAxis         = 1;
  const int Lambda0Axis     = 2;
  const int AntiLambda0Axis = 3;

  string histName = "jet-fragmentation/matching/V0/";
  histName += "fakeV0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> ptbins = getProjectionBins(thn->GetAxis(ptAxis), v0min, v0max);
  thn->GetAxis(ptAxis)->SetRange(ptbins[0], ptbins[1]);
  double lowpt  = thn->GetAxis(ptAxis)->GetBinLowEdge(ptbins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptbins[1]);
  if (lowpt < 0.) { lowpt = 0.; } // Avoids ugly pt>-0 in latextext

  // Apply competing mass cut
  string saveSuffix, latexSuffix;
  if (dM > 1e-3) {
    saveSuffix = TString::Format("_dM%.0f", dM).Data();
    latexSuffix = TString::Format(", #DeltaM = %.0f MeV/#it{c}^{2}", dM).Data();
    double massWindow = dM * 1e-3; // Convert to GeV
    int* coord = new int[nDim]; //Carries the bin coordinates

    array<int, 2> mKBins  = getProjectionBins(thn->GetAxis(K0SAxis), MassK0S - massWindow, MassK0S + massWindow);
    array<int, 2> mLBins  = getProjectionBins(thn->GetAxis(Lambda0Axis), MassLambda0 - massWindow, MassLambda0 + massWindow);
    array<int, 2> mALBins = getProjectionBins(thn->GetAxis(AntiLambda0Axis), MassLambda0 - massWindow, MassLambda0 + massWindow);

    for (int i = 0; i < thn->GetNbins(); i++) {
      double w = thn->GetBinContent(i, coord);

      if ("K0S" == hadron) {
        int mb = coord[Lambda0Axis];
        if (mb >= mLBins[0] && mb <= mLBins[1]) {
          thn->SetBinContent(i, 0.);
        }
        mb = coord[AntiLambda0Axis];
        if (mb >= mALBins[0] && mb <= mALBins[1]) {
          thn->SetBinContent(i, 0.);
        }
      } else {
        int mb = coord[K0SAxis];
        if (mb >= mKBins[0] && mb <= mKBins[1]) {
          thn->SetBinContent(i, 0.);
        }
      }
    } // bin loop
  } // if dM > 1e-3

  int projectionAxis = (hadron == "K0S")*K0SAxis + (hadron == "Lambda0")*Lambda0Axis + (hadron == "AntiLambda0")*AntiLambda0Axis;
  TH1D* combMass = (TH1D*)thn->Projection(projectionAxis);
  combMass->Scale(1./combMass->Integral(), "width");
  setStyle(combMass, 0);
  combMass->SetName(TString::Format("combMass%spt%.1f-%.1f", hadron.c_str(), lowpt, highpt).Data());
  if (isHistEmptyInRange(combMass, 1, combMass->GetNbinsX())) {
    cout << "combMass: histogram empty in non-over/underflow for " << hadron << " in pt range " << lowpt << " - " << highpt << endl;
    return;
  }
  vector<TObject*> objList = {combMass};

  string saveName;
  saveName += "combMass";
  saveName += hadron;
  saveName += TString::Format("_v0pt%.0f-%.0f", lowpt, highpt);
  saveName += saveSuffix;
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  string xTitle = "M(" + formatHadronDaughters(hadron) + ") (GeV/#it{c}^{2})";
  string yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", lowpt, highpt).Data();
  double xMinFrame = combMass->GetXaxis()->GetXmin(), xMaxFrame = combMass->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * combMass->GetBinContent(combMass->GetMaximumBin());
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText + latexSuffix).c_str());

  double mass = ("K0S" == hadron) ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, yMaxFrame);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(0));
  objList.push_back(line);

  canvas->cd();
  frame->Draw();
  for (auto obj : objList) { obj->Draw("same"); }
  canvas->SaveAs(canvas->GetName());
}
void combMass(vector<string> inputStrings, double jetmin, double jetmax, double v0min, double v0max, double dM /* in MeV */, bool doZ = false)
{
  gStyle->SetNdivisions(505, "xy");
  string inName  = inputStrings[0];
  string hadron  = inputStrings[1];
  string dataSet = inputStrings[2];

  double textSize = 0.04;

  const int nDim            = 5;
  const int jetPtAxis       = 0;
  const int v0PtAxis        = 1;
  const int K0SAxis         = 2;
  const int Lambda0Axis     = 3;
  const int AntiLambda0Axis = 4;

  string histName = "jet-fragmentation/matching/jets/V0/";
  if (doZ) histName += "fakeJetPtV0TrackProjMass";
  else histName += "fakeJetPtV0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetPtAxis), jetmin, jetmax);
  thn->GetAxis(jetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  double lowjetpt  = thn->GetAxis(jetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetPtAxis)->GetBinUpEdge(jetptbins[1]);
  array<int, 2> v0bins    = getProjectionBins(thn->GetAxis(v0PtAxis), v0min, v0max);
  thn->GetAxis(v0PtAxis)->SetRange(v0bins[0], v0bins[1]);
  double lowv0  = thn->GetAxis(v0PtAxis)->GetBinLowEdge(v0bins[0]);
  double highv0 = thn->GetAxis(v0PtAxis)->GetBinUpEdge(v0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  // Apply competing mass cut
  string saveSuffix, latexSuffix;
  if (dM > 1e-3) {
    saveSuffix  = TString::Format("_dM%.0f", dM).Data();
    latexSuffix = TString::Format(", #DeltaM = %.0f MeV/#it{c}^{2}", dM).Data();
    double massWindow = dM * 1e-3; // Convert to GeV
    int* coord = new int[nDim]; //Carries the bin coordinates

    array<int, 2> mKBins  = getProjectionBins(thn->GetAxis(K0SAxis), MassK0S - massWindow, MassK0S + massWindow);
    array<int, 2> mLBins  = getProjectionBins(thn->GetAxis(Lambda0Axis), MassLambda0 - massWindow, MassLambda0 + massWindow);
    array<int, 2> mALBins = getProjectionBins(thn->GetAxis(AntiLambda0Axis), MassLambda0 - massWindow, MassLambda0 + massWindow);

    for (int i = 0; i < thn->GetNbins(); i++) {
      double w = thn->GetBinContent(i, coord);

      if ("K0S" == hadron) {
        int mb = coord[Lambda0Axis];
        if (mb >= mLBins[0] && mb <= mLBins[1]) {
          thn->SetBinContent(i, 0.);
        }
        mb = coord[AntiLambda0Axis];
        if (mb >= mALBins[0] && mb <= mALBins[1]) {
          thn->SetBinContent(i, 0.);
        }
      } else {
        int mb = coord[K0SAxis];
        if (mb >= mKBins[0] && mb <= mKBins[1]) {
          thn->SetBinContent(i, 0.);
        }
      }
    } // bin loop
  } // if dM > 1e-3

  int projectionAxis = (hadron == "K0S")*K0SAxis + (hadron == "Lambda0")*Lambda0Axis + (hadron == "AntiLambda0")*AntiLambda0Axis;
  TH1D* combMass = (TH1D*)thn->Projection(projectionAxis);
  combMass->Scale(1./combMass->Integral(), "width");
  setStyle(combMass, 0);
  combMass->SetName(TString::Format("combMass%spt%.1f-%.1f", hadron.c_str(), lowv0, highv0).Data());
  if (isHistEmptyInRange(combMass, 1, combMass->GetNbinsX())) {
    cout << "combMass: histogram empty in non-over/underflow for ptjet " << lowjetpt << " - " << highjetpt << ", " << (doZ ? "z" : "pt") << "v0 " << lowv0 << " - " << highv0 << ". Skipping..." << endl;
    return;
  }
  vector<TObject*> objList = {combMass};

  string saveName;
  saveName += "combMass";
  saveName += hadron;
  saveName += TString::Format("_v0pt%.0f-%.0f", lowv0, highv0);
  saveName += saveSuffix;
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  string xTitle = "M(" + formatHadronDaughters(hadron) + ") (GeV/#it{c}^{2})";
  string yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
  double xMinFrame = combMass->GetXaxis()->GetXmin(), xMaxFrame = combMass->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * combMass->GetBinContent(combMass->GetMaximumBin());
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + latexSuffix).c_str());

  string jetText = TString::Format("#it{p}_{T, ch+V0 jet} = %.0f-%.0f GeV/#it{c}", lowjetpt, highjetpt).Data();
  string v0Text  = TString::Format("#it{p}_{T, V0} = %.1f-%.1f GeV/#it{c}", lowv0, highv0).Data();
  if (doZ) { v0Text = TString::Format("#it{z}_{V0} = %.2f-%.2f", lowv0, highv0).Data(); }
  string latexText = TString::Format("#splitline{%s}{%s}", jetText.c_str(), v0Text.c_str()).Data();
  double xLatex = 0.3, yLatex = 0.8;
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  objList.push_back(latex);

  double mass = ("K0S" == hadron) ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, yMaxFrame);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(0));
  objList.push_back(line);

  canvas->cd();
  frame->Draw();
  for (auto obj : objList) { obj->Draw("same"); }
  canvas->SaveAs(canvas->GetName());
}


// -------------------------------------------------------------------------------------------------

void Reflection(vector<string> inputStrings, double partv0min, double partv0max, double dM /* in MeV */)
{
  gStyle->SetNdivisions(505, "xy");
  string inName  = inputStrings[0];
  string hadron  = inputStrings[1];
  string dataSet = inputStrings[2];

  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim                = 6;
  const int partV0PtAxis        = 0;
  const int detV0PtAxis         = 1;
  const int K0SMassAxis         = 2;
  const int Lambda0MassAxis     = 3;
  const int AntiLambda0MassAxis = 4;
  const int ReflectionAxis      = 5;

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = false;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "Reflected mass (GeV/#it{c}^{2})";
  yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";

  histName = "jet-fragmentation/matching/V0/";
  histName += TString::Format("%sReflection", ("Lambda0" == hadron) ? "Lambda0" : "antiLambda0").Data();
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partv0bins = getProjectionBins(thn->GetAxis(partV0PtAxis), partv0min, partv0max);
  thn->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);
  double lowv0  = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(partv0bins[0]);
  double highv0 = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(partv0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  // Apply competing mass cut to K0S mass
  string saveSuffix;
  string latexSuffix;
  if (dM > 1e-3) {
    saveSuffix = TString::Format("_dM%.0f", dM).Data();
    latexSuffix = TString::Format(", #DeltaM = %.0f MeV/#it{c}^{2}", dM).Data();
    double massWindow = dM * 1e-3; // Convert to GeV
    array<int, 2> mBins = getProjectionBins(thn->GetAxis(K0SMassAxis), MassK0S - massWindow, MassK0S + massWindow);
    int* coord = new int[nDim]; //Carries the bin coordinates
    for (int ib = 0; ib <= thn->GetNbins(); ib++) {
      double w = thn->GetBinContent(ib, coord);
      int mb = coord[K0SMassAxis];
      if (mb >= mBins[0] && mb <= mBins[1]) {
        thn->SetBinContent(ib, 0.);
      }
    }
  }

  TH1D* v0mass = (TH1D*)thn->Projection(ReflectionAxis);
  v0mass->Scale(1./v0mass->Integral(), "width");
  setStyle(v0mass, 0);
  v0mass->SetName(TString::Format("%sReflection", hadron.c_str()).Data());
  if (isHistEmptyInRange(v0mass, 1, v0mass->GetNbinsX())) {
    cout << "Reflection: histogram empty in non-over/underflow for ptv0 " << lowv0 << " - " << highv0 << ". Skipping..." << endl;
    return;
  }

  double xMinLegend = 0.3, xMaxLegend = 0.9, yMinLegend = 0.3, yMaxLegend = 0.5;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  array<double, 2> fitRange = {1.09, 1.15};
  TF1* pol = new TF1("pol", "pol1", fitRange[0], fitRange[1]);
  for (int i = 0; i < pol->GetNpar(); i++) {
    pol->SetParameter(i, 0.);
    pol->SetParLimits(i, -1., 1.);
  }
  TF1* exp1 = new TF1("exp1", "expo", fitRange[0], fitRange[1]);
  for (int i = 0; i < exp1->GetNpar(); i++) {
    exp1->SetParameter(i, 0.);
    exp1->SetParLimits(i, -1., 1.);
  }

  vector<TF1*> fits = {pol, exp1};
  vector<TLatex*> latexFits;
  double xLatex = 0.3, yLatex = 0.8;

  bool emptyFitRange = isHistEmptyInRange(v0mass, fitRange[0], fitRange[1]);
  for (int i = 0; i < fits.size(); i++) {
    if (emptyFitRange) {
      cout << "Reflection: empty fit range for ptv0 " << lowv0 << " - " << highv0 << ". Skipping fits ..." << endl;
      break;
    }
    fits[i]->SetLineColor(GetColor(i+1));
    fits[i]->SetLineWidth(3);
    fits[i]->SetLineStyle(9);
    TFitResultPtr fitptr = v0mass->Fit(fits[i], "RSBQ0");
    string fitName = fits[i]->GetName();
    double chiSq = fitptr->Chi2();
    int ndf = fitptr->Ndf();
    double chiSqNdf = chiSq / ndf;
    TLatex* latexFit = CreateLatex(xLatex, yLatex - 0.05 * i, TString::Format("%s: #chi^{2}/NDF = %.2f / %d = %.2f", fitName.c_str(), chiSq, ndf, chiSqNdf).Data(), textSize);
    legend->AddEntry(fits[i], fitName.c_str(), "l");
    latexFits.push_back(latexFit);
    fits[i]->SetRange(v0mass->GetXaxis()->GetXmin(), v0mass->GetXaxis()->GetXmax());
  }

  vector<TObject*> objList = {v0mass};
  for (auto f : fits) {
    objList.push_back(f);
  }
  for (auto l : latexFits) {
    objList.push_back(l);
  }

  saveName  = hadron;
  saveName += "Reflection";
  saveName += TString::Format("_partv0pt%.0f-%.0f", lowv0, highv0);
  saveName += saveSuffix;
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);

  latexText = dataSet;
  latexText += TString::Format(", #it{p}_{T, %s}^{part.} = %.1f - %.1f GeV/#it{c}", formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  latexText += latexSuffix;

  double xMinFrame = v0mass->GetXaxis()->GetXmin(), xMaxFrame = v0mass->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * v0mass->GetBinContent(v0mass->GetMaximumBin());
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(latexText.c_str());

  TLine* line = new TLine(MassLambda0, yMinFrame, MassLambda0, yMaxFrame);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(0));
  objList.push_back(line);

  if (fits.size() > 0) { objList.push_back(legend); } // Only add legend when plotting fits

  canvas->cd();
  frame->Draw();
  for (auto obj : objList) {
    obj->Draw("same");
  }
  canvas->SaveAs(canvas->GetName());
}
void Reflection(vector<string> inputStrings, double partjetptmin, double partjetptmax, double partv0min, double partv0max, double dM /* in MeV */, bool doZ = false)
{
  string inName  = inputStrings[0];
  string hadron  = inputStrings[1];
  string dataSet = inputStrings[2];
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either Lambda0, or AntiLambda0" << endl;
    return;
  }

  const int nDim                = 8;
  const int partJetPtAxis       = 0;
  const int partV0PtAxis        = 1;
  const int detJetPtAxis        = 2;
  const int detV0PtAxis         = 3;
  const int K0SMassAxis         = 4;
  const int Lambda0MassAxis     = 5;
  const int AntiLambda0MassAxis = 6;
  const int ReflectionAxis      = 7;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  int rebinNumber = 5;
  xTitle = "Reflected mass (GeV/#it{c}^{2})";
  yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";

  histName = TString::Format("partJetPt%sPtDetJetPt%sPt%sReflection", hadron.c_str(), hadron.c_str(), hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProj%sReflection", hadron.c_str(), hadron.c_str(), hadron.c_str()).Data();}
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetbins[0], partjetbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetbins[1]);

  std::array<int, 2> partv0bins = getProjectionBins(thn->GetAxis(partV0PtAxis), partv0min, partv0max);
  thn->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);
  double lowv0 = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(partv0bins[0]);
  double highv0 = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(partv0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  // Competing mass cut
  string saveSuffix;
  string latexSuffix;
  if (dM > 1e-3) {
    saveSuffix = TString::Format("_dM%.0f", dM).Data();
    latexSuffix = TString::Format(", #DeltaM = %.0f MeV/#it{c}^{2}", dM).Data();
    double massWindow = dM * 1e-3; // Convert to GeV
    array<int, 2> mBins = getProjectionBins(thn->GetAxis(K0SMassAxis), MassK0S - massWindow, MassK0S + massWindow);
    int* coord = new int[nDim]; //Carries the bin coordinates
    for (int ib = 0; ib <= thn->GetNbins(); ib++) {
      double w = thn->GetBinContent(ib, coord);
      int mb = coord[K0SMassAxis];
      if (mb >= mBins[0] && mb <= mBins[1]) {
        thn->SetBinContent(ib, 0.);
      }
    }
  }

  TH1D* v0mass = (TH1D*)thn->Projection(ReflectionAxis);
  v0mass->Scale(1./v0mass->Integral(), "width");
  v0mass->SetName(TString::Format("%sReflection", hadron.c_str()).Data());
  setStyle(v0mass, 0);
  if (isHistEmptyInRange(v0mass, 1, v0mass->GetNbinsX())) {
    cout << "Reflection: histogram empty in non-over/underflow for ptjet " << lowjetpt << " - " << highjetpt << ", " << (doZ ? "z" : "pt") << "v0 " << lowv0 << " - " << highv0 << ". Skipping..." << endl;
    return;
  }

  double xMinLegend = 0.3, xMaxLegend = 0.9, yMinLegend = 0.3, yMaxLegend = 0.5;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  array<double, 2> fitRange = {1.09, 1.15};
  TF1* pol = new TF1("pol", "pol1", fitRange[0], fitRange[1]);
  for (int i = 0; i < pol->GetNpar(); i++) {
    pol->SetParameter(i, 0.);
    pol->SetParLimits(i, -1., 1.);
  }
  TF1* exp1 = new TF1("exp1", "expo", fitRange[0], fitRange[1]);
  for (int i = 0; i < exp1->GetNpar(); i++) {
    exp1->SetParameter(i, 0.);
    exp1->SetParLimits(i, -1., 1.);
  }
  vector<TF1*> fits = {pol, exp1};
  vector<TLatex*> latexFits;
  double xLatex = 0.25, yLatex = 0.65;

  bool emptyFitRange = isHistEmptyInRange(v0mass, fitRange[0], fitRange[1]);
  for (int i = 0; i < fits.size(); i++) {
    if (emptyFitRange) {
      cout << "Reflection: empty fit range for ptjet " << lowjetpt << " - " << highjetpt << ", " << (doZ ? "z" : "pt") << "v0 " << lowv0 << " - " << highv0 << ". Skipping fitting..." << endl;
      break;
    }
    fits[i]->SetLineColor(GetColor(i+1));
    fits[i]->SetLineWidth(3);
    fits[i]->SetLineStyle(9);
    TFitResultPtr fitptr = v0mass->Fit(fits[i], "RSBQ0");
    string fitName = fits[i]->GetName();
    double chiSq = fitptr->Chi2();
    int ndf = fitptr->Ndf();
    double chiSqNdf = chiSq / ndf;
    TLatex* latexFit = CreateLatex(xLatex, yLatex - 0.05 * i, TString::Format("%s: #chi^{2}/NDF = %.2f / %d = %.2f", fitName.c_str(), chiSq, ndf, chiSqNdf).Data(), textSize);
    legend->AddEntry(fits[i], fitName.c_str(), "l");
    latexFits.push_back(latexFit);
    fits[i]->SetRange(v0mass->GetXaxis()->GetXmin(), v0mass->GetXaxis()->GetXmax());
  }

  saveName = hadron;
  saveName += "Reflection";
  saveName += TString::Format("_partjetpt%.0f-%.0f", lowjetpt, highjetpt);
  if (doZ) { saveName += TString::Format("_partv0z%.1f-%.1f", lowv0, highv0); }
  else { saveName += TString::Format("_partv0pt%.0f-%.0f", lowv0, highv0); }
  saveName += saveSuffix;
  saveName += ".pdf";
  int xCanvas = 900, yCanvas = 900;
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);

  double xMinFrame = v0mass->GetXaxis()->GetXmin(), xMaxFrame = v0mass->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = max(0.1, 1.1 * v0mass->GetBinContent(v0mass->GetMaximumBin()));
  if (yMaxFrame < 1e-3) { yMaxFrame = 1e-3; } // Change this number if necessary
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + latexSuffix).c_str());

  TLine* line = new TLine(MassLambda0, yMinFrame, MassLambda0, 0.1*yMaxFrame);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(0));
  line->SetLineWidth(3);

  string jetText = TString::Format("#it{p}_{T, ch+V0 jet}^{part.} = %.0f-%.0f GeV/#it{c}", lowjetpt, highjetpt).Data();
  string v0Text  = TString::Format("#it{p}_{T, %s}^{part.} = %.1f-%.1f GeV/#it{c}", formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) { v0Text = TString::Format("#it{z}_{%s}^{part.} = %.2f-%.2f", formatHadronName(hadron).c_str(), lowv0, highv0).Data(); }
  // latexText = TString::Format("#splitline{%s}{#splitline{%s}{%s}}",
  //                             dataSet.c_str(), jetText.c_str(), v0Text.c_str()).Data();
  latexText = TString::Format("#splitline{%s}{%s}",
                              jetText.c_str(), v0Text.c_str()).Data();
  xLatex = 0.35, yLatex = 0.8;
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  vector<TObject*> objList = {v0mass, line};
  for (auto f : fits) { objList.push_back(f); }
  for (auto l : latexFits) { objList.push_back(l); }
  objList.push_back(latex);
  if (fits.size() > 0) { objList.push_back(legend); } // Only add legend when plotting fits

  canvas->cd();
  frame->Draw();
  for (auto obj : objList) { obj->Draw("same"); }
  canvas->SaveAs(canvas->GetName());
}

// -------------------------------------------------------------------------------------------------

void plotTrain(string train, string dataSet, string hadron, double partv0min, double partv0max, double dM /* in MeV */, int setting)
{
  string inName = "~/cernbox/TrainOutput/" + train + "/AnalysisResults.root";
  vector<string> inputStrings = {inName, hadron, dataSet};

  switch (setting) {
    case 0:
      combMass(inputStrings, partv0min, partv0max, dM);
      break;
    case 1:
      Reflection(inputStrings, partv0min, partv0max, dM);
      break;
    default:
      cout << "Error: invalid setting" << endl;
      break;
  }
}
void plotTrain(string train, string dataSet, string hadron, double partjetptmin, double partjetptmax, double partv0min, double partv0max, double dM, bool doZ, int setting)
{
  string inName = "~/cernbox/TrainOutput/" + train + "/AnalysisResults.root";
  vector<string> inputStrings = {inName, hadron, dataSet};
  switch (setting) {
    case 0:
      combMass(inputStrings, partjetptmin, partjetptmax, partv0min, partv0max, dM, doZ);
      break;
    case 1:
      Reflection(inputStrings, partjetptmin, partjetptmax, partv0min, partv0max, dM, doZ);
      break;
    default:
      cout << "Error: invalid setting" << endl;
      break;
  }
}

// For inclusive V0s
void plot271952(string hadron, double v0ptmin, double v0ptmax, double dM, int setting)
{
  string train = "271952";
  string dataSet = "LHC24b1b";
  plotTrain(train, dataSet, hadron, v0ptmin, v0ptmax, dM, setting);
}
// For V0s in jets
void plot271952(string hadron, double jetptmin, double jetptmax, double v0ptmin, double v0ptmax, double dM, bool doZ, int setting)
{
  string train = "271952";
  string dataSet = "LHC24b1b";
  plotTrain(train, dataSet, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax, dM, doZ, setting);
}

// Plot all pt bins, inclusive V0s
void plot271952(string hadron, double dM, int setting)
{
  vector<double> pt = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};
  for (int i = 0; i < pt.size() - 1; i++) {
    plot271952(hadron, pt[i], pt[i+1], dM, setting);
  }
}
// Plot all bins, V0s in jets
void plot271952(string hadron, double partjetptmin, double partjetptmax, double dM, bool doZ, int setting)
{
  vector<double> pt = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};
  vector<double> z = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  vector<double> bins = (doZ) ? z : pt;
  for (int i = 0; i < bins.size() - 1; i++) {
    if (bins[i] > partjetptmax) { break; }
    plot271952(hadron, partjetptmin, partjetptmax, bins[i], bins[i+1], dM, doZ, setting);
  }
}

// Plot all hadrons, jet pt bins, v0 pt bins, v0 z bins, and dM bins
void plotAll271952(int setting)
{
  gROOT->SetBatch(kTRUE);
  vector<string> hadrons = {"K0S", "Lambda0", "AntiLambda0"};
  vector<double> jetpts  = {10.0, 20.0, 40.0};
  vector<double> dMs     = {0.0, 10.0};

  for (auto hadron : hadrons) {
    for (auto dM : dMs) {
      plot271952(hadron, dM, setting);
      for (int i = 0; i < jetpts.size() - 1; i++) {
        plot271952(hadron, jetpts[i], jetpts[i+1], dM, false, setting);
        plot271952(hadron, jetpts[i], jetpts[i+1], dM, true, setting);
      }
    }
  }
}
