
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TString.h"
#include "TLegend.h"

#include "../histUtils.C"

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

// -------------------------------------------------------------------------------------------------
// Functions for background fitting, rejects peak region
// Parameter [0] -> Hi LeftBg Boundary
// Parameter [1] -> Lo RightBg Boundary
double pol1bkg(double *x, double *par)
{
  if (x[0] > par[0] && x[0] < par[1]) {
    TF1::RejectPoint();
    return 0;
  }
  return par[2] + par[3]*x[0];
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
double pol2bkgVar(double *x, double *par)
{
  if (x[0] > par[0] && x[0] < par[1]) {
    TF1::RejectPoint();
    return 0;
  }
  return par[2] + par[3]*x[0] + par[4]*x[0]*x[0];
}

string formatLatexText(int cutAxis, vector<string> axisNames, double binEdge)
{
  const int ptAxis              = 0;
  const int K0SmassAxis         = 1;
  const int Lambda0massAxis     = 2;
  const int AntiLambda0massAxis = 3;
  const int RAxis               = 4;
  const int ctauAxis            = 5;
  const int cosPAAxis           = 6;
  const int DCApAxis            = 7;
  const int DCAnAxis            = 8;
  const int DCAdAxis            = 9;

  string latexText;
  switch (cutAxis) {
    case RAxis:
      latexText = TString::Format("%s > %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case ctauAxis:
      latexText = TString::Format("%s < %.0f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case cosPAAxis:
      latexText = TString::Format("%s > %.3f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case DCApAxis:
      latexText = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case DCAnAxis:
      latexText = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
    case DCAdAxis:
      latexText = TString::Format("%s < %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
      break;
  } // switch (cutAxis)
  return latexText;
}

// -------------------------------------------------------------------------------------------------

array<TFitResultPtr, 3> fitK0S(TH1D* mass, double* sidebandRegion, double* signalRegion, double* fitRegion)
{
  double parameters[5];

  TF1* bkg = new TF1("bkg", pol1bkg, fitRegion[0], fitRegion[1], 4);
  TF1* bkgE = new TF1("bkgE", "pol1", fitRegion[0], fitRegion[1]);
  TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
  TF1* total = new TF1("total", "pol1(0) + gaus(2)", fitRegion[0], fitRegion[1]);

  bkg->SetParameter(0, sidebandRegion[0]); bkg->SetParameter(1, sidebandRegion[1]);
  bkg->FixParameter(0, sidebandRegion[0]); bkg->FixParameter(1, sidebandRegion[1]);
  setStyle(bkg, 1, 2);
  TFitResultPtr bkgPtr = mass->Fit("bkg", "S R");
  parameters[0] = bkg->GetParameter(2); parameters[1] = bkg->GetParameter(3);

  setStyle(bkgE, 1, 2);
  bkgE->SetParameters(parameters);
  bkgE->SetLineStyle(2);
  bkgE->Draw("same");

  setStyle(signal, 2, 2);
  TFitResultPtr sigPtr = mass->Fit("signal", "S R+");
  signal->GetParameters(&parameters[2]);

  setStyle(total, 3, 2);
  total->SetParameters(parameters);
  TFitResultPtr totPtr = mass->Fit("total", "S R+");

  return {bkgPtr, sigPtr, totPtr};
}
array<TFitResultPtr, 3> fitLambda(TH1D* mass, double* sidebandRegion, double* signalRegion, double* fitRegion)
{
  double parameters[6];
  TF1* bkg = new TF1("bkg", pol2bkgVar, fitRegion[0], fitRegion[1], 5);
  // TF1* bkg = new TF1("bkg", "pol2", fitRegion[0], fitRegion[1]);
  TF1* sig = new TF1("sig", "gaus", signalRegion[0], signalRegion[1]);
  TF1* tot = new TF1("tot", "pol2(0) + gaus(3)", fitRegion[0], fitRegion[1]);
  bkg->SetParameter(0, sidebandRegion[0]); bkg->SetParameter(1, sidebandRegion[1]);
  bkg->FixParameter(0, sidebandRegion[0]); bkg->FixParameter(1, sidebandRegion[1]);

  setStyle(bkg, 1);
  TFitResultPtr bkgPtr = mass->Fit("bkg", "S R");
  parameters[0] = bkg->GetParameter(2);
  parameters[1] = bkg->GetParameter(3);
  parameters[2] = bkg->GetParameter(4);
  // return {nullptr, nullptr, nullptr};

  TF1* bkgExtr = new TF1("bkgExtr", "pol2", fitRegion[0], fitRegion[1]);
  bkgExtr->SetParameter(0, bkg->GetParameter(2));
  bkgExtr->SetParameter(1, bkg->GetParameter(3));
  bkgExtr->SetParameter(2, bkg->GetParameter(4));
  setStyle(bkgExtr, 1);
  bkgExtr->SetLineStyle(9);
  bkgExtr->Draw("same");

  setStyle(sig, 2);
  TFitResultPtr sigPtr = mass->Fit("sig", "S R+");
  sig->GetParameters(&parameters[3]);

  setStyle(tot, 3);
  tot->SetParameters(parameters);
  TFitResultPtr totPtr = mass->Fit("tot", "S R+");

  return {bkgPtr, sigPtr, totPtr};
}
void cutVarK0SPurity(string inName = "", int cutAxis = -1, double ptmin = 0., double ptmax = 100.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  string hypothesis = "K0S";
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if (cutAxis < 4 || cutAxis > 9) {
    cout << "Error: cutAxis must be specified (4-9)" << endl;
    return;
  }
  const int nDim                = 10;
  const int ptAxis              = 0;
  const int K0SmassAxis         = 1;
  const int Lambda0massAxis     = 2;
  const int AntiLambda0massAxis = 3;
  const int RAxis               = 4;
  const int ctauAxis            = 5;
  const int cosPAAxis           = 6;
  const int DCApAxis            = 7;
  const int DCAnAxis            = 8;
  const int DCAdAxis            = 9;

  vector<string> axisNames = {"pt", "K0Smass", "Lambda0mass", "AntiLambda0mass", "R", "ctau", "cosPA", "DCAp", "DCAn", "DCAd"};
  int projectionAxis = ("K0S" == hypothesis)*K0SmassAxis + ("Lambda0" == hypothesis)*Lambda0massAxis + ("AntiLambda0" == hypothesis)*AntiLambda0massAxis;
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.06;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 0.08;
  if ("K0S" == hypothesis) { xMinFrame = 0.4, xMaxFrame = 0.6, yMaxFrame = 0.1; }
  double xMinLegend = 0.1, xMaxLegend = 0.7, yMinLegend = 0.1, yMaxLegend = 0.45;
  double xLatex = 0.1, yLatex = 0.9;
  int xCanvas = 1800, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}_{%s} (GeV/#it{c}^{2})", formatHadronName(hypothesis).c_str()).Data();
  yTitle = "#frac{N}{N(uncut)}";
  dataSet = "LHC22o_pass6_minBias_small";

  std::vector<TH1D*> histVector; std::vector<TF1*> funcVector;
  histName = "V0CutVariation";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  // Apply pt selection for all histograms
  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);

  double normalisationFactor = 1.;
  double uncutSignal = 1.;

  double sidebandRegion[2] = {0.46221947, 0.53135842}; // mu ± 7 sigma
  double signalRegion[2] = {0.47209646, 0.52148143}; // mu ± 5 sigma
  double fitRegion[2] = {0.4, 0.6};

  int underflowBin = 0, overflowBin = 1 + thn->GetAxis(cutAxis)->GetNbins();
  for (int iBin = underflowBin; iBin <= overflowBin; iBin++) {
    THnSparseD* thn_copy = (THnSparseD*)thn->Clone("thn_copy");
    thn_copy->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]); // Just to be sure

    // Cut the axis off at the bottom (var > binLowEdge)
    thn_copy->GetAxis(cutAxis)->SetRange(iBin, overflowBin);
    if ( (ctauAxis == cutAxis) || (DCAdAxis == cutAxis) ) {
      // Cut the axis off at the top (var < binUpEdge)
      thn_copy->GetAxis(cutAxis)->SetRange(underflowBin, overflowBin - iBin);
    }

    TH1D* mass = (TH1D*)thn_copy->Projection(projectionAxis);
    mass->SetName(TString::Format("mass_axis%d_bin%d", cutAxis, iBin).Data());
    setStyle(mass, 0);

    // Latex stuff
    double binEdge = thn_copy->GetAxis(cutAxis)->GetBinLowEdge(iBin);
    if ( (ctauAxis == cutAxis) || (DCAdAxis == cutAxis) ) {
      binEdge = thn_copy->GetAxis(cutAxis)->GetBinUpEdge(iBin);
    }
    latexText = formatLatexText(cutAxis, axisNames, binEdge);
    if (0 == iBin) {
      normalisationFactor = mass->Integral();
      latexText = "No cut";
    }
    mass->Scale(1./normalisationFactor);

    TCanvas* canvas = new TCanvas(TString::Format("Plot_a%d_b%d", cutAxis, iBin).Data(), TString::Format("Plot_a%d_b%d", cutAxis, iBin).Data(), xCanvas, yCanvas);
    TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
    canvas->Divide(2, 1);
    canvas->cd(1);
    TH1F* tmpframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
    gPad->SetRightMargin(0.05);
    tmpframe->Draw();
    mass->Draw("same");

    // Purity extraction
    array<TFitResultPtr, 3> fitresults = fitK0S(mass, sidebandRegion, signalRegion, fitRegion);
    TFitResultPtr bkgPtr = fitresults[0];

    TF1* bkgE = new TF1("bkgE", "pol1", fitRegion[0], fitRegion[1]);
    TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
    TF1* total = new TF1("total", "pol1(0) + gaus(2)", fitRegion[0], fitRegion[1]);

    setStyle(bkgE, 1);
    setStyle(signal, 2);
    setStyle(total, 3);

    bkgE->SetParameter(0, *(bkgPtr->GetParams() + 2));
    bkgE->SetParameter(1, *(bkgPtr->GetParams() + 3));

    signal->SetParameter(0, *(fitresults[1]->GetParams()));
    signal->SetParameter(1, *(fitresults[1]->GetParams() + 1));

    total->SetParameter(0, *(fitresults[2]->GetParams()));
    total->SetParameter(1, *(fitresults[2]->GetParams() + 1));
    total->SetParameter(2, *(fitresults[2]->GetParams() + 2));
    total->SetParameter(3, *(fitresults[2]->GetParams() + 3));
    total->SetParameter(4, *(fitresults[2]->GetParams() + 4));

    legend->AddEntry(mass, "data");
    legend->AddEntry(bkgE, "background");
    legend->AddEntry(signal, "signal");
    legend->AddEntry(total, "combined");

    // Draw latex with purity in here.
    double bkgEstimate = bkgE->Integral(signalRegion[0], signalRegion[1]);
    double bkgErr = bkgE->IntegralError(signalRegion[0], signalRegion[1], bkgE->GetParameters(), bkgPtr->GetCovarianceMatrix().GetMatrixArray());
    double sigPlusBkgErr;
    double sigPlusBkg = mass->IntegralAndError(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]), sigPlusBkgErr, "width");
    double purity = 1 - bkgEstimate / sigPlusBkg;
    double sigEst = sigPlusBkg * purity;

    if (0 == iBin) { uncutSignal = sigEst; }
    double efficiency = sigEst / uncutSignal;

    // if (cosPAAxis != cutAxis) {
    //   latexText = TString::Format("%s, purity %.2f %%", latexText.c_str(), 1e2*purity).Data();
    // }
    canvas->cd(2);
    gPad->SetLeftMargin(0.);
    // TLatex* tmplatex = CreateLatex(xLatex, yLatex, latexText.c_str(), textSize);
    // tmplatex->Draw("same");

    TLatex* latexData       = CreateLatex(xLatex, yLatex, dataSet.c_str(), textSize);
    TLatex* latexCut        = CreateLatex(xLatex, yLatex - 0.1, latexText.c_str(), textSize);
    TLatex* latexPurity     = CreateLatex(xLatex, yLatex - 0.2, TString::Format("Purity: %.2f %%", 1e2*purity).Data(), textSize);
    TLatex* latexEfficiency = CreateLatex(xLatex, yLatex - 0.3, TString::Format("Rel. Efficiency: %.2f %%", 1e2*efficiency).Data(), textSize);
    TLatex* latexPt         = CreateLatex(xLatex, yLatex - 0.4, TString::Format("#it{p}_{T, V0} = %.0f - %.0f GeV/#it{c}", ptmin, ptmax).Data(), textSize);

    // if (cosPAAxis == cutAxis) {
    //   string sublatextext = TString::Format("%.3f: p %.2f%%, #eta %.2f%%", binEdge, 1e2*purity, 1e2*efficiency).Data();
    //   if (0 == iBin) { sublatextext = TString::Format("No cut: p %.2f%%", 1e2*purity).Data(); }
    //   TLatex* sublatex = CreateLatex(0.05, 0.8 - iBin*0.15, sublatextext.c_str(), textSize);
    //   canvas->cd(overflowBin + 3);
    //   sublatex->Draw("same");
    // }
    // else if (0 != iBin){
    //   string sublatextext = TString::Format("#eta %.2f%%", 1e2*efficiency).Data();
    //   TLatex* sublatex = CreateLatex(xLatex + 0.1, yLatex - 0.1, sublatextext.c_str(), textSize);
    //   sublatex->Draw("same");
    // }

    // latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, V0} = %.0f - %.0f GeV/#it{c} }", dataSet.c_str(), ptmin, ptmax).Data();
    // TLatex* latex = CreateLatex(0.2, 0.7, latexText, textSize);
    // latex->Draw();
    latexData->Draw();
    latexCut->Draw("same");
    latexPurity->Draw("same");
    latexEfficiency->Draw("same");
    latexPt->Draw("same");
    legend->Draw("same");

    saveName = TString::Format("purity%s", hypothesis.c_str()).Data();
    saveName = TString::Format("%s_%s", saveName.c_str(), axisNames[cutAxis].c_str()).Data();
    saveName = TString::Format("%s_bin%d", saveName.c_str(), iBin).Data();
    saveName = TString::Format("%s_v0pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax).Data();
    saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
    canvas->SaveAs(saveName.c_str());
  }


}
void cutVarLambdaPurity(string inName = "", int cutAxis = -1, double ptmin = 0., double ptmax = 100.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  string hypothesis = "Lambda0";
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if (cutAxis < 4 || cutAxis > 9) {
    cout << "Error: cutAxis must be specified (4-9)" << endl;
    return;
  }
  const int nDim                = 10;
  const int ptAxis              = 0;
  const int K0SmassAxis         = 1;
  const int Lambda0massAxis     = 2;
  const int AntiLambda0massAxis = 3;
  const int RAxis               = 4;
  const int ctauAxis            = 5;
  const int cosPAAxis           = 6;
  const int DCApAxis            = 7;
  const int DCAnAxis            = 8;
  const int DCAdAxis            = 9;

  std::array<string, nDim> axisNames = {"pt", "K0Smass", "Lambda0mass", "AntiLambda0mass", "R", "ctau", "cosPA", "DCAp", "DCAn", "DCAd"};
  int projectionAxis = ("K0S" == hypothesis)*K0SmassAxis + ("Lambda0" == hypothesis)*Lambda0massAxis + ("AntiLambda0" == hypothesis)*AntiLambda0massAxis;
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.08; //0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 0.08;
  if ("K0S" == hypothesis) { xMinFrame = 0.4, xMaxFrame = 0.6, yMaxFrame = 0.1; }
  double xMinLegend = 0.2, xMaxLegend = 0.8, yMinLegend = 0.1, yMaxLegend = 0.5;
  double xLatex = 0.23, yLatex = 0.93;
  int xCanvas = 1800, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}_{%s} (GeV/#it{c}^{2})", formatHadronName(hypothesis).c_str()).Data();
  yTitle = "#frac{N}{N(uncut)}";
  dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector; std::vector<TF1*> funcVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  if (RAxis == cutAxis || cosPAAxis == cutAxis) { canvas->Divide(4, 2); }
  else { canvas->Divide(3, 2);}

  histName = "V0CutVariation";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  // Apply pt selection for all histograms
  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);

  double normalisationFactor = 1.;
  double uncutSignal = 1.;

  double sidebandRegion[2] = {1.105, 1.125}; // Hand picked, improve in future
  double signalRegion[2] = {1.1125, 1.12}; // Hand picked, improve in future
  double fitRegion[2] = {1.09, 1.215};

  int underflowBin = 0, overflowBin = 1 + thn->GetAxis(cutAxis)->GetNbins();
  for (int iBin = underflowBin; iBin <= overflowBin; iBin++) {
    THnSparseD* thn_copy = (THnSparseD*)thn->Clone("thn_copy");
    thn_copy->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]); // Just to be sure

    // Cut the axis off at the bottom (var > binLowEdge)
    thn_copy->GetAxis(cutAxis)->SetRange(iBin, overflowBin);
    if ( (ctauAxis == cutAxis) || (DCAdAxis == cutAxis) ) {
      // Cut the axis off at the top (var < binUpEdge)
      thn_copy->GetAxis(cutAxis)->SetRange(underflowBin, overflowBin - iBin);
    }

    TH1D* mass = (TH1D*)thn_copy->Projection(projectionAxis);
    mass->SetName(TString::Format("mass_axis%d_bin%d", cutAxis, iBin).Data());
    setStyle(mass, 0);

    // Latex stuff
    double binEdge = thn_copy->GetAxis(cutAxis)->GetBinLowEdge(iBin);
    if ( (ctauAxis == cutAxis) || (DCAdAxis == cutAxis) ) {
      binEdge = thn_copy->GetAxis(cutAxis)->GetBinUpEdge(iBin);
    }
    switch (cutAxis) {
      case RAxis:
        latexText = TString::Format("%s > %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case ctauAxis:
        latexText = TString::Format("%s < %.0f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case cosPAAxis:
        latexText = TString::Format("%s > %.3f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCApAxis:
        latexText = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCAnAxis:
        latexText = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCAdAxis:
        latexText = TString::Format("%s < %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
    } // switch (cutAxis)

    if (0 == iBin) {
      normalisationFactor = mass->Integral();
      latexText = "No cut";
    }
    mass->Scale(1./normalisationFactor);

    canvas->cd(iBin+1);
    TH1F* tmpframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
    tmpframe->Draw();
    mass->Draw("same");

    // Purity extraction
    array<TFitResultPtr, 3> fitresults = fitLambda(mass, sidebandRegion, signalRegion, fitRegion);
    TFitResultPtr bkgPtr = fitresults[0];

    TF1* bkgE = new TF1("bkgE", "pol2", fitRegion[0], fitRegion[1]);
    TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
    TF1* total = new TF1("total", "pol2(0) + gaus(3)", fitRegion[0], fitRegion[1]);

    setStyle(bkgE, 1);
    setStyle(signal, 2);
    setStyle(total, 3);

    bkgE->SetParameter(0, *(fitresults[0]->GetParams() + 2));
    bkgE->SetParameter(1, *(fitresults[0]->GetParams() + 3));
    bkgE->SetParameter(2, *(fitresults[0]->GetParams() + 4));

    signal->SetParameter(0, *(fitresults[1]->GetParams()));
    signal->SetParameter(1, *(fitresults[1]->GetParams() + 1));

    total->SetParameter(0, *(fitresults[2]->GetParams()));
    total->SetParameter(1, *(fitresults[2]->GetParams() + 1));
    total->SetParameter(2, *(fitresults[2]->GetParams() + 2));
    total->SetParameter(3, *(fitresults[2]->GetParams() + 3));
    total->SetParameter(4, *(fitresults[2]->GetParams() + 4));

    if (0 == iBin) {
      legend->AddEntry(mass, "data");
      legend->AddEntry(bkgE, "background");
      legend->AddEntry(signal, "signal");
      legend->AddEntry(total, "combined");
    }

    // Draw latex with purity in here.
    double bkgEstimate = bkgE->Integral(signalRegion[0], signalRegion[1]);
    double bkgErr = bkgE->IntegralError(signalRegion[0], signalRegion[1], bkgE->GetParameters(), bkgPtr->GetCovarianceMatrix().GetMatrixArray());
    double sigPlusBkgErr;
    double sigPlusBkg = mass->IntegralAndError(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]), sigPlusBkgErr, "width");
    double purity = 1 - bkgEstimate / sigPlusBkg;
    double sigEst = sigPlusBkg * purity;

    if (0 == iBin) { uncutSignal = sigEst; }
    double efficiency = sigEst / uncutSignal;

    if (cosPAAxis != cutAxis) {
      latexText = TString::Format("%s, purity %.2f %%", latexText.c_str(), 1e2*purity).Data();
    }
    TLatex* tmplatex = CreateLatex(xLatex, yLatex, latexText.c_str(), textSize);
    tmplatex->Draw("same");

    if (cosPAAxis == cutAxis) {
      string sublatextext = TString::Format("%.3f: p %.2f%%, #eta %.2f%%", binEdge, 1e2*purity, 1e2*efficiency).Data();
      if (0 == iBin) { sublatextext = TString::Format("No cut: p %.2f%%", 1e2*purity).Data(); }
      TLatex* sublatex = CreateLatex(0.05, 0.8 - iBin*0.15, sublatextext.c_str(), textSize);
      canvas->cd(overflowBin + 3);
      sublatex->Draw("same");
    }
    else if (0 != iBin){
      string sublatextext = TString::Format("#eta %.2f%%", 1e2*efficiency).Data();
      TLatex* sublatex = CreateLatex(xLatex + 0.1, yLatex - 0.1, sublatextext.c_str(), textSize);
      sublatex->Draw("same");
    }
  }

  canvas->cd(overflowBin + 2);
  latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, V0} = %.0f - %.0f GeV/#it{c} }", dataSet.c_str(), ptmin, ptmax).Data();
  TLatex* latex = CreateLatex(0.2, 0.7, latexText, textSize);
  latex->Draw();
  legend->Draw("same");

  saveName = TString::Format("purity%s", hypothesis.c_str()).Data();
  saveName = TString::Format("%s_%s", saveName.c_str(), axisNames[cutAxis].c_str()).Data();
  saveName = TString::Format("%s_v0pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax).Data();
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(saveName.c_str());
}
void cutVarAntiLambdaPurity(string inName = "", int cutAxis = -1, double ptmin = 0., double ptmax = 100.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  string hypothesis = "AntiLambda0";
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if (cutAxis < 4 || cutAxis > 9) {
    cout << "Error: cutAxis must be specified (4-9)" << endl;
    return;
  }
  const int nDim                = 10;
  const int ptAxis              = 0;
  const int K0SmassAxis         = 1;
  const int Lambda0massAxis     = 2;
  const int AntiLambda0massAxis = 3;
  const int RAxis               = 4;
  const int ctauAxis            = 5;
  const int cosPAAxis           = 6;
  const int DCApAxis            = 7;
  const int DCAnAxis            = 8;
  const int DCAdAxis            = 9;

  std::array<string, nDim> axisNames = {"pt", "K0Smass", "Lambda0mass", "AntiLambda0mass", "R", "ctau", "cosPA", "DCAp", "DCAn", "DCAd"};
  int projectionAxis = ("K0S" == hypothesis)*K0SmassAxis + ("Lambda0" == hypothesis)*Lambda0massAxis + ("AntiLambda0" == hypothesis)*AntiLambda0massAxis;
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.08; //0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 0.08;
  if ("K0S" == hypothesis) { xMinFrame = 0.4, xMaxFrame = 0.6, yMaxFrame = 0.1; }
  double xMinLegend = 0.2, xMaxLegend = 0.8, yMinLegend = 0.1, yMaxLegend = 0.5;
  double xLatex = 0.23, yLatex = 0.93;
  int xCanvas = 1800, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}_{%s} (GeV/#it{c}^{2})", formatHadronName(hypothesis).c_str()).Data();
  yTitle = "#frac{N}{N(uncut)}";
  dataSet = "LHC22o_pass6_minBias";

  std::vector<TH1D*> histVector; std::vector<TF1*> funcVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  if (RAxis == cutAxis || cosPAAxis == cutAxis) { canvas->Divide(4, 2); }
  else { canvas->Divide(3, 2);}

  histName = "V0CutVariation";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  // Apply pt selection for all histograms
  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);

  double normalisationFactor = 1.;
  double uncutSignal = 1.;

  double sidebandRegion[2] = {1.105, 1.125}; // Hand picked, improve in future
  double signalRegion[2] = {1.1125, 1.12}; // Hand picked, improve in future
  double fitRegion[2] = {1.09, 1.215};

  int underflowBin = 0, overflowBin = 1 + thn->GetAxis(cutAxis)->GetNbins();
  for (int iBin = underflowBin; iBin <= overflowBin; iBin++) {
    THnSparseD* thn_copy = (THnSparseD*)thn->Clone("thn_copy");
    thn_copy->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]); // Just to be sure

    // Cut the axis off at the bottom (var > binLowEdge)
    thn_copy->GetAxis(cutAxis)->SetRange(iBin, overflowBin);
    if ( (ctauAxis == cutAxis) || (DCAdAxis == cutAxis) ) {
      // Cut the axis off at the top (var < binUpEdge)
      thn_copy->GetAxis(cutAxis)->SetRange(underflowBin, overflowBin - iBin);
    }

    TH1D* mass = (TH1D*)thn_copy->Projection(projectionAxis);
    mass->SetName(TString::Format("mass_axis%d_bin%d", cutAxis, iBin).Data());
    setStyle(mass, 0);

    // Latex stuff
    double binEdge = thn_copy->GetAxis(cutAxis)->GetBinLowEdge(iBin);
    if ( (ctauAxis == cutAxis) || (DCAdAxis == cutAxis) ) {
      binEdge = thn_copy->GetAxis(cutAxis)->GetBinUpEdge(iBin);
    }
    switch (cutAxis) {
      case RAxis:
        latexText = TString::Format("%s > %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case ctauAxis:
        latexText = TString::Format("%s < %.0f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case cosPAAxis:
        latexText = TString::Format("%s > %.3f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCApAxis:
        latexText = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCAnAxis:
        latexText = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCAdAxis:
        latexText = TString::Format("%s < %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
    } // switch (cutAxis)

    if (0 == iBin) {
      normalisationFactor = mass->Integral();
      latexText = "No cut";
    }
    mass->Scale(1./normalisationFactor);

    canvas->cd(iBin+1);
    TH1F* tmpframe = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
    tmpframe->Draw();
    mass->Draw("same");

    // Purity extraction
    array<TFitResultPtr, 3> fitresults = fitLambda(mass, sidebandRegion, signalRegion, fitRegion);
    TFitResultPtr bkgPtr = fitresults[0];

    TF1* bkgE = new TF1("bkgE", "pol2", fitRegion[0], fitRegion[1]);
    TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
    TF1* total = new TF1("total", "pol2(0) + gaus(3)", fitRegion[0], fitRegion[1]);

    setStyle(bkgE, 1);
    setStyle(signal, 2);
    setStyle(total, 3);

    bkgE->SetParameter(0, *(fitresults[0]->GetParams() + 2));
    bkgE->SetParameter(1, *(fitresults[0]->GetParams() + 3));
    bkgE->SetParameter(2, *(fitresults[0]->GetParams() + 4));

    signal->SetParameter(0, *(fitresults[1]->GetParams()));
    signal->SetParameter(1, *(fitresults[1]->GetParams() + 1));

    total->SetParameter(0, *(fitresults[2]->GetParams()));
    total->SetParameter(1, *(fitresults[2]->GetParams() + 1));
    total->SetParameter(2, *(fitresults[2]->GetParams() + 2));
    total->SetParameter(3, *(fitresults[2]->GetParams() + 3));
    total->SetParameter(4, *(fitresults[2]->GetParams() + 4));

    if (0 == iBin) {
      legend->AddEntry(mass, "data");
      legend->AddEntry(bkgE, "background");
      legend->AddEntry(signal, "signal");
      legend->AddEntry(total, "combined");
    }

    // Draw latex with purity in here.
    double bkgEstimate = bkgE->Integral(signalRegion[0], signalRegion[1]);
    double bkgErr = bkgE->IntegralError(signalRegion[0], signalRegion[1], bkgE->GetParameters(), bkgPtr->GetCovarianceMatrix().GetMatrixArray());
    double sigPlusBkgErr;
    double sigPlusBkg = mass->IntegralAndError(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]), sigPlusBkgErr, "width");
    double purity = 1 - bkgEstimate / sigPlusBkg;
    double sigEst = sigPlusBkg * purity;

    if (0 == iBin) { uncutSignal = sigEst; }
    double efficiency = sigEst / uncutSignal;

    if (cosPAAxis != cutAxis) {
      latexText = TString::Format("%s, purity %.2f %%", latexText.c_str(), 1e2*purity).Data();
    }
    TLatex* tmplatex = CreateLatex(xLatex, yLatex, latexText.c_str(), textSize);
    tmplatex->Draw("same");

    if (cosPAAxis == cutAxis) {
      string sublatextext = TString::Format("%.3f: p %.2f%%, #eta %.2f%%", binEdge, 1e2*purity, 1e2*efficiency).Data();
      if (0 == iBin) { sublatextext = TString::Format("No cut: p %.2f%%", 1e2*purity).Data(); }
      TLatex* sublatex = CreateLatex(0.05, 0.8 - iBin*0.15, sublatextext.c_str(), textSize);
      canvas->cd(overflowBin + 3);
      sublatex->Draw("same");
    }
    else if (0 != iBin){
      string sublatextext = TString::Format("#eta %.2f%%", 1e2*efficiency).Data();
      TLatex* sublatex = CreateLatex(xLatex + 0.1, yLatex - 0.1, sublatextext.c_str(), textSize);
      sublatex->Draw("same");
    }
  }

  canvas->cd(overflowBin + 2);
  latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, V0} = %.0f - %.0f GeV/#it{c} }", dataSet.c_str(), ptmin, ptmax).Data();
  TLatex* latex = CreateLatex(0.2, 0.7, latexText, textSize);
  latex->Draw();
  legend->Draw("same");

  saveName = TString::Format("purity%s", hypothesis.c_str()).Data();
  saveName = TString::Format("%s_%s", saveName.c_str(), axisNames[cutAxis].c_str()).Data();
  saveName = TString::Format("%s_v0pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax).Data();
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(saveName.c_str());
}
