
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
#include "RooFit.h"

#include "../histUtils.C"

const double MassK0S     = 0.497611;
const double MassLambda0 = 1.115683;
gStyle->SetNdivisions(505, "xy");

string getDataSet(int train)
{
  if (252064 == train) return "LHC22o_pass6";
  if (282430 == train) return "LHC22o_pass7_small";
  return "Could not find dataset";
}

// Prints the parameter limits of a function
void printParLimits(TF1* f)
{
  for (int i = 0; i < f->GetNpar(); i++) {
    double min, max;
    f->GetParLimits(i, min, max);
    cout << f->GetName() << "(" << i << " = " << f->GetParName(i) << ") " << f->GetParameter(i) << " (" << min << ", " << max << ")" << endl;
  }
}

// Returns a subset of a histogram from minBin to maxBin
template <typename T>
T* makeHistSubset(T* data, int minBin, int maxBin)
{
  T* region = (T*)data->Clone("region");
  region->Reset();
  for (int i = minBin; i <= maxBin; i++) {
    region->SetBinContent(i, data->GetBinContent(i));
  }
  return region;
}

// Add RooFit to legend, ensures that line width in legend is appropriate
void AddToLegend(TLegend* legend, TH1* h, string name, int styleNumber)
{
  TH1* k = (TH1*)h->Clone(name.c_str());
  setStyle(k, styleNumber);
  legend->AddEntry(k, name.c_str(), "l");
}

// -------------------------------------------------------------------------------------------------
//
// Functions to use with setup
//
// -------------------------------------------------------------------------------------------------

// Get mass hist for given train, hadron, pt
TH1* getHist(string inName, string hadron, double& ptmin, double& ptmax)
{
  // string inName = "~/cernbox/TrainOutput/252064/AnalysisResults.root";
  // string inName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string histName = "jet-fragmentation/data/V0/V0CutVariation";

  int ptAxis = 0;
  int mAxis;
  if ("K0S" == hadron) {
    mAxis = 1;
  }
  else if ("Lambda0" == hadron) {
    mAxis = 2;
  }
  else if ("AntiLambda0" == hadron) {
    mAxis = 3;
  }
  else {
    cout << "Hadron " << hadron << " not recognized" << endl;
    return nullptr;
  }

  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*) inFile->Get(histName.c_str());

  array<int, 2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  ptmin = ptBins[0];
  ptmax = ptBins[1];
  return (thn->Projection(mAxis));
}

// Plots multiple fits for signal+bkg
void plotBkgs(vector<string> inputStrings, double ptmin, double ptmax)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  const int ptAxis = 0;
  const int K0SAxis = 1;
  const int Lambda0Axis = 2;
  const int AntiLambda0Axis = 3;
  int projectionAxis = ("K0S" == hadron)*K0SAxis + ("Lambda0" == hadron)*Lambda0Axis + ("AntiLambda0" == hadron)*AntiLambda0Axis;

  double lowpt = ptmin, highpt = ptmax;
  TH1D* hist = (TH1D*)getHist(inName, hadron, lowpt, highpt);
  setStyle(hist, 0);
  hist->SetName("data");

  string saveName = hadron;
  saveName += TString::Format("_pt%.1f-%.1f", ptmin, ptmax).Data();
  saveName += "_roofit";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 1800, 900);

  array<double, 2> fitRange = {1.08, 1.2};
  double mass = MassLambda0;
  double dM = 1e-2;
  double signalWidth = 1e-3;
  if ("K0S" == hadron) {
    fitRange = {0.45, 0.55};
    mass = MassK0S;
    dM = 5e-2;
    signalWidth = 1e-2;
  }

  string ptText = TString::Format("#it{p}_{T, V0} = %.1f - %.1f GeV/#it{c}", ptmin, ptmax).Data();
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(hadron).c_str()).Data();
  string yTitle = "counts";

  // RooRealVar x = RooRealVar("x", "x", 1.08, 1.215);
  // x.setRange("FULL", 1.08, 1.215);
  // x.setRange("SIGNAL", 1.1, 1.13);
  // x.setRange("LHS", 1.08, 1.1);
  // x.setRange("RHS", 0.54, 1.215);
  // K0S values
  RooRealVar x = RooRealVar("x", "x", 0.4, 0.6);
  x.setRange("FULL", 0.4, 0.6);
  x.setRange("SIGNAL", 0.48, 0.52);
  x.setRange("LHS", 0.4, 0.46);
  x.setRange("RHS", 0.54, 0.6);
  x.setRange("RTAIL", 0.5, 0.6);
  double yMaxFrame = 1.5 * getHistScale(hist, false);
  RooPlot* frame = x.frame();
  frame->SetTitle((dataSet + ", " + ptText).c_str());
  frame->GetXaxis()->SetTitle(xTitle.c_str());
  frame->GetYaxis()->SetTitle(yTitle.c_str());
  frame->GetYaxis()->SetRangeUser(0., yMaxFrame);

  RooPlot* frame2 = x.frame();
  frame2->GetXaxis()->SetTitle(xTitle.c_str());
  frame2->GetYaxis()->SetTitle("Residual");

  RooPlot* frame3 = x.frame();
  frame3->GetXaxis()->SetTitle(xTitle.c_str());
  frame3->GetYaxis()->SetTitle("Pull");

  RooDataHist dh("dh", "title", RooArgList(x), hist);

  dh.plotOn(frame, RooFit::Name("dh"), RooFit::LineWidth(3), RooFit::MarkerColor(GetColor(0)), RooFit::LineColor(GetColor(0)), RooFit::MarkerStyle(GetMarker(0)));
  TLegend* legend = CreateLegend(0.15, 0.8, 0.7, 0.9, "", 0.04);
  legend->AddEntry(hist, hist->GetName(), "l");

  // Gaussian
  RooRealVar sigmean("sigmean", "K mass", mass, mass - 1e-2, mass + 1e-2);
  RooRealVar sigwidth("sigwidth", "K width", signalWidth, 1e-6, 5 * signalWidth);
  RooGaussian sigGaus("sigGaus", "sigGaus", x, sigmean, sigwidth);
  // sigGaus.fitTo(dh, RooFit::Range("SIGNAL"));

  // Gaussian
  RooRealVar bkgMean("bkgMean", "bkgMean", mass, mass - 1e-2, mass + 1e-2);
  RooRealVar bkgSigma("bkgSigma", "bkgSigma", 1e-2, 0.2);
  RooGaussian bkgGaus("bkgGaus", "bkgGaus", x, bkgMean, bkgSigma);

  // Linear
  RooRealVar linA("linA", "linA", 0., 4e4);
  RooRealVar linB("linB", "linB", 0., 4e6);
  RooGenericPdf linear("linear","linA+linB*(x)", RooArgList(linA, linB, x));
  // linear.fitTo(dh, RooFit::Range("RHS"));

  // Quadratic
  RooRealVar quadA("quadA", "c", 0., 4e4);
  RooRealVar quadB("quadB", "b", 0., 4e6);
  RooRealVar quadC("quadC", "a", 0., 4e6);
  RooGenericPdf quadratic("quadratic","quadA+quadB*x+quadC*x*x", RooArgList(quadA, quadB, quadC, x));
  // quadratic.fitTo(dh, RooFit::Range("RHS"));

  // Power law
  RooRealVar powA("powA", "powA", 3e6);
  RooRealVar powB("powB", "powB", -4e6, 4e6);
  RooRealVar powC("powC", "powC", -4e6, 4e6);
  RooRealVar powD("powD", "powD", -4e6, 4e6);
  RooRealVar powE("powE", "powE", -4e6, 4e6);
  RooRealVar powF("powF", "powF", -4e6, 4e6);
  RooGenericPdf power("power","powA*pow((x-powB),powC)", RooArgList(powA, powB, powC, x));
  power.fitTo(dh, RooFit::Range("RTAIL"));
  power.plotOn(frame, RooFit::Range("FULL"), RooFit::LineColor(GetColor(5)));
  AddToLegend(legend, hist, "(x-a)^b", 5);

  RooGenericPdf power2("power2","powA*pow((x-powB),-1*(x-powD))", RooArgList(powA, powB, powD, x));
  power2.fitTo(dh, RooFit::Range("RTAIL"));
  power2.plotOn(frame, RooFit::Range("FULL"), RooFit::LineColor(GetColor(6)));
  AddToLegend(legend, hist, "(x-a)^-(x-b)", 6);

  RooGenericPdf power3("power3","powA*pow((x-powB),-1*(x-powE)*(x-powE))", RooArgList(powA, powB, powE, x));
  power3.fitTo(dh, RooFit::Range("RTAIL"));
  power3.plotOn(frame, RooFit::Range("FULL"), RooFit::LineColor(GetColor(7)));
  AddToLegend(legend, hist, "(x-a)^-(x-b)^2", 7);


  // Gaussian + Linear
  RooRealVar fSignalLin("fSignalLin", "fSignalLin", 0.5, 1.);
  RooAddPdf sigPlusLin("sigPlusLin", "sigPlusLin", RooArgList(sigGaus, linear), fSignalLin);
  sigPlusLin.fitTo(dh, RooFit::Range("SIGNAL"));
  // sigPlusLin.plotOn(frame, RooFit::Range("FULL"), RooFit::Components(sigGaus), RooFit::LineColor(GetColor(1)));
  // sigPlusLin.plotOn(frame, RooFit::Range("FULL"), RooFit::Components(linear), RooFit::LineColor(GetColor(2)));
  sigPlusLin.plotOn(frame, RooFit::Range("FULL"), RooFit::LineColor(GetColor(1)));
  AddToLegend(legend, hist, "sigPlusLin", 1);

  // Gaussian + Gaussian + Linear
  RooRealVar fSignalGausLin("fSignalGausLin", "fSignalGausLin", 0.5, 1.);
  RooRealVar fBkgGausLin("fBkgGausLin", "fBkgGausLin", 0., 0.3);
  RooAddPdf sigPlusGausLin("sigPlusGausLin", "sigPlusGausLin", RooArgList(sigGaus, bkgGaus, linear), RooArgList(fSignalGausLin, fBkgGausLin));
  // sigPlusGausLin.fitTo(dh, RooFit::Range("SIGNAL"));
  // sigPlusGausLin.plotOn(frame, RooFit::Range("FULL"), RooFit::LineColor(GetColor(5)));
  // AddToLegend(legend, hist, "sigPlusGausLin", 5);

  // Crystal ball
  RooRealVar c_x0("c_x0", "c_x0", mass, mass - 3*dM, mass + 3*dM);
  RooRealVar c_sigmaL("c_sigmaL", "c_sigmaL", signalWidth, 1e-6, 5 * signalWidth);
  RooRealVar c_alphaL("c_alphaL", "c_alphaL", 2., 1., 10.);
  RooRealVar c_nL("c_nL", "c_nL", 1., 1e-5, 10.);
  RooRealVar c_sigmaR("c_sigmaR", "c_sigmaR", signalWidth, 1e-6, 5 * signalWidth);
  RooRealVar c_alphaR("c_alphaR", "c_alphaR", 2., 1., 10.);
  RooRealVar c_nR("c_nR", "c_nR", 5., 1e-5, 10.);
  // RooCrystalBall c_cb("c_cb", "c_cb", x, c_x0, c_sigmaL, c_sigmaR, c_alphaL, c_nL, c_alphaR, c_nR);
  RooCrystalBall c_cb("c_cb", "c_cb", x, c_x0, c_sigmaL, c_sigmaL, c_alphaL, c_nL, c_alphaR, c_nR);
  // c_cb.fitTo(dh, RooFit::Range("SIGNAL"));
  c_cb.fitTo(dh, RooFit::Range("FULL"));
  c_cb.plotOn(frame, RooFit::Name("c_cb"), RooFit::Range("FULL"), RooFit::NormRange("SIGNAL"), RooFit::LineColor(GetColor(3)));
  double c_chiSq = frame->chiSquare();
  string c_name = TString::Format("c_cb: #chi^{2}/NDF=%.1f, #sigma=%.1e, #alpha_{L}=%.1f, n_{L}=%.1f, #alpha_{R}=%.1f, n_{R}=%.1f", c_chiSq, c_sigmaL.getValV(), c_alphaL.getValV(), c_nL.getValV(), c_alphaR.getValV(), c_nR.getValV()).Data();
  // string c_name = TString::Format("c_cb: #splitline{#chi^{2}/NDF=%.1f}{#splitline{#sigma_{L}=%.1e, #alpha_{L}=%.1f, n_{L}=%.1f}{#sigma_{R}=%.1e, #alpha_{R}=%.1f, n_{R}=%.1f}}", c_chiSq, c_sigmaL.getValV(), c_alphaL.getValV(), c_nL.getValV(), c_sigmaR.getValV(), c_alphaR.getValV(), c_nR.getValV()).Data();
  AddToLegend(legend, hist, c_name, 3);
  TLine* c_lineL = new TLine(c_x0.getValV() - c_alphaL.getValV() * c_sigmaL.getValV(), 0., c_x0.getValV() - c_alphaL.getValV() * c_sigmaL.getValV(), 0.5 * yMaxFrame);
  setStyle(c_lineL, 3);
  TLine* c_lineR = new TLine(c_x0.getValV() + c_alphaL.getValV() * c_sigmaL.getValV(), 0., c_x0.getValV() + c_alphaL.getValV() * c_sigmaL.getValV(), 0.5 * yMaxFrame);
  // TLine* c_lineR = new TLine(c_x0.getValV() + c_alphaR.getValV() * c_sigmaR.getValV(), 0., c_x0.getValV() + c_alphaR.getValV() * c_sigmaR.getValV(), 0.5 * yMaxFrame);
  setStyle(c_lineR, 3);

  RooCrystalBall d_cb("d_cb", "d_cb", x, c_x0, c_sigmaL, c_sigmaR, c_alphaL, c_nL, c_alphaR, c_nR);
  d_cb.fitTo(dh, RooFit::Range("FULL"));
  d_cb.plotOn(frame, RooFit::Name("d_cb"), RooFit::Range("FULL"), RooFit::NormRange("SIGNAL"), RooFit::LineColor(GetColor(4)));
  double d_chiSq = frame->chiSquare();
  string d_name = TString::Format("d_cb: #chi^{2}/NDF=%.1f, #sigma_{L}=%.1e, #alpha_{L}=%.1f, n_{L}=%.1f, #sigma_{R}=%.1e, #alpha_{R}=%.1f, n_{R}=%.1f", d_chiSq,  c_sigmaL.getValV(),  c_alphaL.getValV(),  c_nL.getValV(),  c_sigmaR.getValV(),  c_alphaR.getValV(),  c_nR.getValV()).Data();
  AddToLegend(legend, hist, d_name, 4);

  // vector<RooAbsPDF> pdfs = {sigPlusLin, sigPlusGausLin, c_cb, d_cb};

  frame->Draw();
  legend->Draw("same");
  // a_line->Draw("same");
  // b_line->Draw("same");
  // c_lineL->Draw("same");
  // c_lineR->Draw("same");
  canvas->SaveAs(saveName.c_str());

  // RooHist* pull = frame->residHist("data", "c_cb", true, false);
  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 1800, 900);
  canvas2->cd();
  RooHist* c_resid = frame->residHist(dh.GetName(), c_cb.GetName());
  c_resid->SetLineColor(GetColor(3));
  c_resid->SetMarkerColor(GetColor(3));
  frame2->addPlotable(c_resid, "P");
  RooHist* d_resid = frame->residHist(dh.GetName(), d_cb.GetName());
  d_resid->SetLineColor(GetColor(4));
  d_resid->SetMarkerColor(GetColor(4));
  frame2->addPlotable(d_resid, "P");
  frame2->Draw();

  TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 1800, 900);
  canvas3->cd();
  RooHist* c_pull = frame->pullHist(dh.GetName(), c_cb.GetName());
  c_pull->SetLineColor(GetColor(3));
  c_pull->SetMarkerColor(GetColor(3));
  frame3->addPlotable(c_pull, "P");
  RooHist* d_pull = frame->pullHist(dh.GetName(), d_cb.GetName());
  d_pull->SetLineColor(GetColor(4));
  d_pull->SetMarkerColor(GetColor(4));
  frame3->addPlotable(d_pull, "P");
  frame3->Draw();
}

// -------------------------------------------------------------------------------------------------
//
// Interface
//
// -------------------------------------------------------------------------------------------------

void plotTrain(int train, string hadron, double v0min, double v0max, int setting)
{
  string inName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string dataSet = getDataSet(train);
  vector<string> inputStrings = {inName, dataSet, hadron};

  switch(setting) {
    case 0:
      plotBkgs(inputStrings, v0min, v0max);
      break;
    // case 1:
    //   plotBkgParts(inputStrings, v0min, v0max);
    //   break;
    default:
      cout << "Setting " << setting << " not recognized" << endl;
  }
}

void plot252064(string hadron, double v0min, double v0max, int setting) { plotTrain(252064, hadron, v0min, v0max, setting); }
void plot282430(string hadron, double v0min, double v0max, int setting) { plotTrain(282430, hadron, v0min, v0max, setting); }

void test(){plot252064("K0S", 5., 10., 0);}