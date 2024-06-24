
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

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// Makes histogram to be used as a shaded region
template <typename T>
T* fillShadedRegion(T* data, int minBin, int maxBin)
{
  T* region = (T*)data->Clone("region");
  region->Reset();
  for (int i = minBin; i <= maxBin; i++) {
    region->SetBinContent(i, data->GetBinContent(i));
  }
  return region;
}

// -------------------------------------------------------------------------------------------------
// Initialise parameters and set bounds
void parsK0S(TF1* f, bool doubleGauss, int bkg)
{
  double gausPar[3] = {1., MassK0S, 0.01};
  double gausVar[6] = {0.4, 1., MassK0S - 5e-2, MassK0S + 5e-2, 1e-3, 0.02};
  double doubleGaussPar[3] = {0.1, MassK0S, 0.05};
  double doubleGaussVar[6] = {0.1, 0.6, MassK0S - 0.01, MassK0S + 0.01, 1e-3, 0.1};
  double bkgPar[3] = {0., 0., 0.};
  double bkgVar[6] = {-0.1, 0.1, -0.1, 0.1, -0.1, 0.1};

  f->SetParameter(0, gausPar[0]); f->SetParLimits(0, gausVar[0], gausVar[1]);
  f->SetParameter(1, gausPar[1]); f->SetParLimits(1, gausVar[2], gausVar[3]);
  f->SetParameter(2, gausPar[2]); f->SetParLimits(2, gausVar[4], gausVar[5]);
  if (doubleGauss) {
    f->SetParameter(3, doubleGaussPar[0]); f->SetParLimits(3, doubleGaussVar[0], doubleGaussVar[1]);
    f->SetParameter(4, doubleGaussPar[1]); f->SetParLimits(4, doubleGaussVar[2], doubleGaussVar[3]);
    f->SetParameter(5, doubleGaussPar[2]); f->SetParLimits(5, doubleGaussVar[4], doubleGaussVar[5]);
    f->SetParameter(6, bkgPar[0]);         f->SetParLimits(6, bkgVar[0], bkgVar[1]);
    f->SetParameter(7, bkgPar[1]);         f->SetParLimits(7, bkgVar[2], bkgVar[3]);
    if (2 >= bkg) {
      f->SetParameter(8, bkgPar[2]);       f->SetParLimits(8, bkgVar[4], bkgVar[5]);
    }
  }
  else {
    f->SetParameter(3, bkgPar[0]);   f->SetParLimits(3, bkgVar[0], bkgVar[1]);
    f->SetParameter(4, bkgPar[1]);   f->SetParLimits(4, bkgVar[2], bkgVar[3]);
    if (2 >= bkg) {
      f->SetParameter(5, bkgPar[2]); f->SetParLimits(5, bkgVar[4], bkgVar[5]);
    }
  }
}
void parsL0(TF1* f, bool doubleGauss, int bkg)
{
  double gausPar[3] = {1., MassLambda0, 0.01};
  double gausVar[6] = {0.4, 1., MassLambda0 - 0.01, MassLambda0 + 0.01, 1e-4, 0.02};
  double doubleGaussPar[3] = {0.1, MassLambda0, 0.05};
  double doubleGaussVar[6] = {0.1, 0.6, MassLambda0 - 0.03, MassLambda0 + 0.03, 1e-3, 0.1};
  double bkgPar[3] = {0., 0., 0.};
  double bkgVar[6] = {-3., 3., -3., 3., -3., 3.};

  f->SetParameter(0, gausPar[0]); f->SetParLimits(0, gausVar[0], gausVar[1]);
  f->SetParameter(1, gausPar[1]); f->SetParLimits(1, gausVar[2], gausVar[3]);
  f->SetParameter(2, gausPar[2]); f->SetParLimits(2, gausVar[4], gausVar[5]);
  if (doubleGauss) {
    f->SetParameter(3, doubleGaussPar[0]); f->SetParLimits(3, doubleGaussVar[0], doubleGaussVar[1]);
    f->SetParameter(4, doubleGaussPar[1]); f->SetParLimits(4, doubleGaussVar[2], doubleGaussVar[3]);
    f->SetParameter(5, doubleGaussPar[2]); f->SetParLimits(5, doubleGaussVar[4], doubleGaussVar[5]);
    f->SetParameter(6, bkgPar[0]);         f->SetParLimits(6, bkgVar[0], bkgVar[1]);
    f->SetParameter(7, bkgPar[1]);         f->SetParLimits(7, bkgVar[2], bkgVar[3]);
    if (2 >= bkg) {
      f->SetParameter(8, bkgPar[2]);       f->SetParLimits(8, bkgVar[4], bkgVar[5]);
    }
  }
  else {
    f->SetParameter(3, bkgPar[0]);   f->SetParLimits(3, bkgVar[0], bkgVar[1]);
    f->SetParameter(4, bkgPar[1]);   f->SetParLimits(4, bkgVar[2], bkgVar[3]);
    if (2 >= bkg) {
      f->SetParameter(5, bkgPar[2]); f->SetParLimits(5, bkgVar[4], bkgVar[5]);
    }
  }
}
// Get fit function for signal and background
TF1* getSigBkgFit(string hadron, int bkg, bool doubleGauss)
{
  double fitRegion[2] = {1.08, 1.215};
  if ("K0S" == hadron) { fitRegion[0] = 0.44, fitRegion[1] = 0.555; }

  int arg = 3;
  string sFit = "gaus(0)";
  if (doubleGauss) {
    sFit += " + gaus(" + to_string(arg) + ")";
    arg += 3;
  }
  sFit += " + pol" + to_string(bkg) + "(" + to_string(arg) + ")";
  TF1* fit = new TF1("fit", sFit.c_str(), fitRegion[0], fitRegion[1]);

  if ("K0S" == hadron) {
    parsK0S(fit, doubleGauss, bkg);
  }
  else {
    parsL0(fit, doubleGauss, bkg);
  }
  return fit;
}
// Get shaded regions for signal and background
array<TH1*, 2> getSigBkgRegions(TH1* mass, double mean, double sigma, double nSigmaSignal, double nSigmaBkgMin, double nSigmaBkgMax)
{
  array<int, 2> sigBins      = getProjectionBins(mass->GetXaxis(), mean - nSigmaSignal * sigma, mean + nSigmaSignal * sigma);
  array<int, 2> leftBkgBins  = getProjectionBins(mass->GetXaxis(), mean - nSigmaBkgMax * sigma, mean - nSigmaBkgMin * sigma);
  array<int, 2> rightBkgBins = getProjectionBins(mass->GetXaxis(), mean + nSigmaBkgMin * sigma, mean + nSigmaBkgMax * sigma);
  // Prevent overlap of signal and background regions
  if (leftBkgBins[1] >= sigBins[0]) {
    int shift = leftBkgBins[1] - sigBins[0] + 1;
    leftBkgBins[0] -= shift;
    leftBkgBins[1] -= shift;
    cout << "Warning: Overlap of signal and background regions. Shifting left background region by " << shift << " bins." << endl;
  }
  if (rightBkgBins[0] <= sigBins[1]) {
    int shift = sigBins[1] - rightBkgBins[0] + 1;
    rightBkgBins[0] += shift;
    rightBkgBins[1] += shift;
    cout << "Warning: Overlap of signal and background regions. Shifting right background region by " << shift << " bins." << endl;
  }
  TH1D* sigRegion = (TH1D*)fillShadedRegion(mass, sigBins[0], sigBins[1]);
  TH1D* bkgRegion = (TH1D*)fillShadedRegion(mass, leftBkgBins[0], leftBkgBins[1]);
  bkgRegion->Add(fillShadedRegion(mass, rightBkgBins[0], rightBkgBins[1]));
  return {sigRegion, bkgRegion};
}

void cutVarPurity(string inName, string dataSet, string hadron, int cutAxis, double ptmin, double ptmax, int bkg, bool doubleGauss, bool flipGaussians)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
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

  int nSigmaSignal = 3, nSigmaBkgMin = 5, nSigmaBkgMax = 8;
  gStyle->SetNdivisions(505, "xy");
  int projectionAxis = ("K0S" == hadron)*K0SmassAxis + ("Lambda0" == hadron)*Lambda0massAxis + ("AntiLambda0" == hadron)*AntiLambda0massAxis;
  bool reverseCuts = (ctauAxis == cutAxis) || (DCAdAxis == cutAxis);
  double refValue = 1.;
  int refBin = 0;

  vector<string> axisNames = {"pt", "K0Smass", "Lambda0mass", "AntiLambda0mass", "R", "ctau", "cosPA", "DCAp", "DCAn", "DCAd"};
  string histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 1.07, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 1.1;
  if ("K0S" == hadron) { xMinFrame = 0.4, xMaxFrame = 0.6; }
  double xMinLegend = 0.1, xMaxLegend = 0.7, yMinLegend = 0.1, yMaxLegend = 0.45;
  double xLatex = 0.23, yLatex = 0.93;
  int xCanvas = 1800, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  yTitle = "arb. units";

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Line for hadron mass
  TLine* lineMass = new TLine(("K0S" == hadron) ? MassK0S : MassLambda0, yMinFrame, ("K0S" == hadron) ? MassK0S : MassLambda0, yMaxFrame);
  lineMass->SetLineColor(GetColor(0));
  lineMass->SetLineWidth(2);
  lineMass->SetLineStyle(9);

  std::vector<TH1D*> histVector;
  histName = "V0CutVariation";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
  string sPt = TString::Format("%.1f < #it{p}_{T, V0} < %.1f GeV/#it{c}", lowpt, highpt).Data();

  int underflowBin = 0, overflowBin = 1 + thn->GetAxis(cutAxis)->GetNbins();
  TH1D* hEfficiency = new TH1D("hEfficiency", "", overflowBin+1, thn->GetAxis(cutAxis)->GetBinLowEdge(0), thn->GetAxis(cutAxis)->GetBinUpEdge(overflowBin));
  TH1D* hPurity = (TH1D*)hEfficiency->Clone("hPurity");

  for (int iBin = underflowBin; iBin <= overflowBin; iBin++) {
    THnSparseD* thn_copy = (THnSparseD*)thn->Clone("thn_copy");
    thn_copy->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);

    // Cut the axis off at the bottom (var > binLowEdge)
    thn_copy->GetAxis(cutAxis)->SetRange(iBin, overflowBin);
    if (reverseCuts) { // Cut the axis off at the top (var < binUpEdge)
      thn_copy->GetAxis(cutAxis)->SetRange(underflowBin, overflowBin - iBin);
    }
    TH1D* mass = (TH1D*)thn_copy->Projection(projectionAxis);
    mass->SetName(TString::Format("mass_axis%d_bin%d", cutAxis, iBin).Data());
    setStyle(mass, 0);

    // Normalise to max value (usually peak)
    double normalisationFactor = mass->GetBinContent(mass->GetMaximumBin());
    mass->Scale(1./normalisationFactor);
    if (0 == iBin) {
      latexText = "No cut";
    }
    else {
      double binEdge = thn_copy->GetAxis(cutAxis)->GetBinLowEdge(iBin);
      if (reverseCuts) {
        binEdge = thn_copy->GetAxis(cutAxis)->GetBinUpEdge(iBin);
      }
      latexText = formatLatexText(cutAxis, axisNames, binEdge);
    }

    string canvasName = TString::Format("canvas_%s_%s_%d", hadron.c_str(), axisNames[cutAxis].c_str(), iBin).Data();
    TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), xCanvas, yCanvas);
    TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
    canvas->Divide(2,1);
    canvas->cd(1);
    TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
    gPad->SetRightMargin(0.05);
    frame->Draw();
    mass->Draw("same");
    legend->AddEntry(mass, "Data");

    TF1* fit = getSigBkgFit(hadron, bkg, doubleGauss);
    TFitResultPtr fitResult = mass->Fit(fit, "RBSQ");
    legend->AddEntry(fit, "Fit", "l");
    double mean  = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    if (doubleGauss && flipGaussians) {
      mean  = fit->GetParameter(4);
      sigma = fit->GetParameter(5);
    }

    array<TH1*, 2> regions = getSigBkgRegions(mass, mean, sigma, nSigmaSignal, nSigmaBkgMin, nSigmaBkgMax);
    TH1D* sigRegion = (TH1D*)regions[0];
    TH1D* bkgRegion = (TH1D*)regions[1];
    sigRegion->SetFillColorAlpha(kGreen, 0.3);
    bkgRegion->SetFillColorAlpha(kRed, 0.3);
    legend->AddEntry(sigRegion, "Signal+Background", "f");
    legend->AddEntry(bkgRegion, "Background", "f");
    sigRegion->Draw("bars same");
    bkgRegion->Draw("bars same");
    lineMass->Draw("same");

    double sigPlusBkg = sigRegion->Integral();
    double background = bkgRegion->Integral();
    double purity = 1. - background/sigPlusBkg;
    if (0 == iBin) {
      refValue = normalisationFactor * sigPlusBkg * purity;
    }
    double efficiency = normalisationFactor * sigPlusBkg * purity;
    efficiency /= refValue;

    hEfficiency->SetBinContent(iBin + 1, efficiency);
    hPurity->SetBinContent(iBin + 1, purity);

    canvas->cd(2);
    gPad->SetLeftMargin(0.);
    TLatex* lData   = CreateLatex(0.1, 0.85, dataSet.c_str(), textSize);
    TLatex* lPt     = CreateLatex(0.1, 0.8, sPt.c_str(), textSize);
    TLatex* lCut    = CreateLatex(0.1, 0.7, latexText.c_str(), textSize);
    TLatex* lPurity = CreateLatex(0.1, 0.65, TString::Format("Purity: %.1f%%", purity*1e2), textSize);
    TLatex* lEff    = CreateLatex(0.1, 0.6, TString::Format("Rel. Efficiency: %.2f %%", 1e2*efficiency).Data(), textSize);
    if (legend) legend->Draw("same");
    lData->Draw("same");
    lPt->Draw("same");
    lCut->Draw("same");
    lPurity->Draw("same");
    lEff->Draw("same");

    string saveName = TString::Format("purity%s", hadron.c_str()).Data();
    saveName = TString::Format("%s_%s", saveName.c_str(), axisNames[cutAxis].c_str()).Data();
    saveName = TString::Format("%s_bin%d", saveName.c_str(), iBin).Data();
    saveName = TString::Format("%s_v0pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax).Data();
    saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
    canvas->SaveAs(saveName.c_str());
  } // for iBin

  TCanvas* cEP  = new TCanvas("cEP", "cEP", xCanvas, yCanvas);
  TLatex* lData = CreateLatex(0.25, 0.35, dataSet.c_str(), textSize);
  TLatex* lPt   = CreateLatex(0.25, 0.3, sPt.c_str(), textSize);
  TLegend* lEP  = CreateLegend(0.25, 0.5, 0.15, 0.25, legendTitle, textSize);
  setStyle(hEfficiency, 0);
  setStyle(hPurity, 1);
  lEP->AddEntry(hEfficiency, "Relative Efficiency");
  lEP->AddEntry(hPurity, "Purity");

  for (int i = underflowBin; i <= overflowBin; i++) {
    string latexText;
    if (0 == i) {
      latexText = "No cut";
    }
    else {
      double binEdge = thn->GetAxis(cutAxis)->GetBinLowEdge(i);
      if (reverseCuts) {
        binEdge = thn->GetAxis(cutAxis)->GetBinUpEdge(i);
      }
      latexText = formatLatexText(cutAxis, axisNames, binEdge);
    }
    hEfficiency->GetXaxis()->SetBinLabel(i + 1, latexText.c_str());
  }
  hEfficiency->SetStats(0);

  hEfficiency->Draw();
  lEP->Draw("same");
  hPurity->Draw("same");

  string saveName = "eff-pur";
  saveName = TString::Format("%s_%s", saveName.c_str(), axisNames[cutAxis].c_str()).Data();
  saveName = TString::Format("%s_v0pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax).Data();
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  cEP->SaveAs(saveName.c_str());
}
