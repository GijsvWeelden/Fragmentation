
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

double pol1bkg(double *x, double *par)
{
  if (x[0] > 0.46221947 && x[0] < 0.53135842) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0];
}
// Is there a smart way to implement these separately for Lambda and antiLambda?
// * Lambda0     = 1.10192774, 1.12922828
// * AntiLambda0 = 1.10181581, 1.12952655
double pol2bkg(double *x, double *par)
{
  if (x[0] > 1.10192774 && x[0] < 1.12922828) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}
double pol3bkg(double *x, double *par)
{
  if (x[0] > par[0] && x[0] < par[1]) {
    TF1::RejectPoint();
    return 0;
  }
  return par[2] + par[3]*x[0] + par[4]*x[0]*x[0] + par[5]*x[0]*x[0]*x[0];
}

// -------------------------------------------------------------------------------------------------

double getNjets(TFile* inFile, double jetptmin, double jetptmax)
{
  string histName = "jet-fragmentation/matching/jets/matchPartJetPtEtaPhi";
  TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
  std::array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
  return jets->Integral(jetptbins[0], jetptbins[1], 0, jets->GetNbinsY() + 1, 0, jets->GetNbinsZ() + 1);
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

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) { latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.2f - %.2f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data(); }
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

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) { latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.2f - %.2f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data(); }
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

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  if (doZ) { latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{z}_{%s}^{part.} = %.2f - %.2f }}", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data(); }
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

void ReflectionNoJets(string inName, string hadron, double partv0min = -1., double partv0max = 1e6)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either Lambda0, or AntiLambda0" << endl;
    return;
  }
  string hypothesis = ("Lambda0" == hadron) ? "AntiLambda0" : "Lambda0";
  const int nDim                = 6;
  const int partV0PtAxis        = 0;
  const int detV0PtAxis         = 1;
  const int K0SMassAxis         = 2;
  const int Lambda0MassAxis     = 3;
  const int AntiLambda0MassAxis = 4;
  const int ReflectionAxis      = 5;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 0.08;
  double xMinLegend = 0.3, xMaxLegend = 0.9, yMinLegend = 0.3, yMaxLegend = 0.5;
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

  histName = TString::Format("%sReflection", ("Lambda0" == hadron) ? "Lambda0" : "antiLambda0").Data();
  histName = TString::Format("jet-fragmentation/matching/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partv0bins = getProjectionBins(thn->GetAxis(partV0PtAxis), partv0min, partv0max);
  thn->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);
  double lowv0 = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(partv0bins[0]);
  double highv0 = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(partv0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  TH1D* v0mass = (TH1D*)thn->Projection(ReflectionAxis);
  v0mass->Scale(1./v0mass->Integral());
  setStyle(v0mass, 0);
  v0mass->SetName(TString::Format("%sReflection", hadron.c_str()).Data());
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
  TF1* exp1 = new TF1("exp1", "expo", fitRegion[0], fitRegion[1]);
  TF1* exp1_extrapolated = new TF1("exp1_extrapolated", "expo", 1.05, 1.215);

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

  latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }", dataSet.c_str(), formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  TLatex* latexPolChi2 = CreateLatex(0.26, 0.7, TString::Format("Pol: #chi^{2}/NDF = %.2f / %d = %.2f", polptr->Chi2(), polptr->Ndf(), polptr->Chi2() / polptr->Ndf()).Data(), textSize);
  TLatex* latexPol = CreateLatex(0.26, 0.65, TString::Format("Pol = %.2f + %.2f M", parameters_pol[0], parameters_pol[1]).Data(), textSize);
  TLatex* latexExpChi2 = CreateLatex(0.26, 0.6, TString::Format("Exp: #chi^{2}/NDF = %.2f / %d = %.2f", expptr->Chi2(), expptr->Ndf(), expptr->Chi2() / expptr->Ndf()).Data(), textSize);
  TLatex* latexExp = CreateLatex(0.26, 0.55, TString::Format("Exp = exp(%.2f + %.2f M)", parameters_exp[0], parameters_exp[1]).Data(), textSize);

  TLine* line = new TLine(MassLambda0, yMinFrame, MassLambda0, yMaxFrame);

  legend->Draw("same");
  latex->Draw("same");
  latexPolChi2->Draw("same");
  latexPol->Draw("same");
  latexExpChi2->Draw("same");
  latexExp->Draw("same");

  saveName = TString::Format("%sReflection", hadron.c_str()).Data();
  saveName = TString::Format("%s_partv0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  canvas->SaveAs(saveName.c_str());
}
void Reflection(string inName, string hadron, double partjetptmin = 10., double partjetptmax = 200., double partv0min = -1., double partv0max = 1e6, bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either Lambda0, or AntiLambda0" << endl;
    return;
  }
  string hypothesis = ("Lambda0" == hadron) ? "AntiLambda0" : "Lambda0";
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
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 0.08;
  double xMinLegend = 0.3, xMaxLegend = 0.9, yMinLegend = 0.3, yMaxLegend = 0.5;
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

  histName = TString::Format("partJetPt%sPtDetJetPt%sPt%sReflection", hadron.c_str(), hadron.c_str(), hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProj%sReflection", hadron.c_str(), hadron.c_str(), hadron.c_str()).Data();}
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  std::array<int, 2> partv0bins = getProjectionBins(thn->GetAxis(partV0PtAxis), partv0min, partv0max);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetbins[0], partjetbins[1]);
  thn->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetbins[1]);
  double lowv0 = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(partv0bins[0]);
  double highv0 = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(partv0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  TH1D* v0mass = (TH1D*)thn->Projection(ReflectionAxis);
  v0mass->Scale(1./v0mass->Integral());
  setStyle(v0mass, 0);
  v0mass->SetName(TString::Format("%sReflection", hadron.c_str()).Data());
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
  TF1* exp1 = new TF1("exp1", "expo", fitRegion[0], fitRegion[1]);
  TF1* exp1_extrapolated = new TF1("exp1_extrapolated", "expo", 1.05, 1.215);

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

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch. jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} } }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0, highv0).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  TLatex* latexPolChi2 = CreateLatex(0.26, 0.7, TString::Format("Pol: #chi^{2}/NDF = %.2f / %d = %.2f", polptr->Chi2(), polptr->Ndf(), polptr->Chi2() / polptr->Ndf()).Data(), textSize);
  TLatex* latexPol = CreateLatex(0.26, 0.65, TString::Format("Pol = %.2f + %.2f M", parameters_pol[0], parameters_pol[1]).Data(), textSize);
  TLatex* latexExpChi2 = CreateLatex(0.26, 0.6, TString::Format("Exp: #chi^{2}/NDF = %.2f / %d = %.2f", expptr->Chi2(), expptr->Ndf(), expptr->Chi2() / expptr->Ndf()).Data(), textSize);
  TLatex* latexExp = CreateLatex(0.26, 0.55, TString::Format("Exp = exp(%.2f + %.2f M)", parameters_exp[0], parameters_exp[1]).Data(), textSize);

  TLine* line = new TLine(MassLambda0, yMinFrame, MassLambda0, yMaxFrame);

  legend->Draw("same");
  latex->Draw("same");
  latexPolChi2->Draw("same");
  latexPol->Draw("same");
  latexExpChi2->Draw("same");
  latexExp->Draw("same");

  saveName = TString::Format("%sReflection", hadron.c_str()).Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  if (doZ) {saveName = TString::Format("%s_partv0z%.1f-%.1f", saveName.c_str(), lowv0, highv0); }
  else { saveName = TString::Format("%s_partv0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0); }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  canvas->SaveAs(saveName.c_str());
}
