
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

// -------------------------------------------------------------------------------------------------

void K0SPurity(string inName = "AnalysisResults.root", double ptmin = 0., double ptmax = 100.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  const int nDim = 4;
  const int ptAxis = 0;
  const int K0SmassAxis = 1;
  const int Lambda0massAxis = 2;
  const int AntiLambda0massAxis = 3;
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 0., yMaxFrame = 60e6;
  double xMinLegend = 0.25, xMaxLegend = 0.5, yMinLegend = 0.65, yMaxLegend = 0.85;
  double xLatex = 0.23, yLatex = 0.93;
  int xCanvas = 1800, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}_{%s} (GeV/#it{c}^{2})", formatHadronName("K0S").c_str()).Data();
  yTitle = "count";
  dataSet = "LHC23y_pass1";

  std::vector<TH1D*> histVector; std::vector<TF1*> funcVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "V0PtMass";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
  TH1D* mass = (TH1D*)thn->Projection(K0SmassAxis);
  // mass->Scale(1./mass->Integral()); // TODO: Normalising shows S+B=0 in TLatex. How to deal with this?
  setStyle(mass, 0);
  legend->AddEntry(mass, "data");

  canvas->cd();
  frame->Draw();
  mass->Draw("same");

  double parameters[5];
  double signalRegion[2] = {0.47209646, 0.52148143}; // mu ± 5 sigma
  double fitRegion[2] = {0.4, 0.6};
  TF1* background_extrapolated = new TF1("background_extrapolated", "pol1", fitRegion[0], fitRegion[1]);
  // TF1* background = new TF1("background", "pol1", fitRegion[0], fitRegion[1]);
  TF1* background = new TF1("background", pol1bkg, fitRegion[0], fitRegion[1], 2);
  TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
  TF1* total = new TF1("total", "pol1(0) + gaus(2)", fitRegion[0], fitRegion[1]);

  setStyle(background, 1);
  TFitResultPtr bkgPtr = mass->Fit("background", "S");
  background->GetParameters(&parameters[0]);
  legend->AddEntry(background, "background");

  setStyle(background_extrapolated, 1);
  background_extrapolated->SetParameters(parameters);
  background_extrapolated->SetLineStyle(9);

  setStyle(signal, 2);
  mass->Fit("signal", "R+");
  signal->GetParameters(&parameters[2]);
  legend->AddEntry(signal, "signal");

  setStyle(total, 3);
  total->SetParameters(parameters);
  mass->Fit("total", "R+");
  legend->AddEntry(total, "combined");

  if (legend) legend->Draw("same");
  latexText = TString::Format("%s, %.0f < #it{p}_{T, V0} < %.0f GeV/#it{c}", dataSet.c_str(), ptmin, ptmax).Data();
  TLatex* latex = CreateLatex(0.3, 0.95, latexText, textSize);
  latex->Draw("same");

  double bkgOffset = total->GetParameter(0); double bkgOffsetErr = total->GetParError(0);
  double bkgSlope = total->GetParameter(1);  double bkgSlopeErr = total->GetParError(1);
  double norm = total->GetParameter(2);      double normErr = total->GetParError(2);
  double mu = total->GetParameter(3);        double muErr = total->GetParError(3);
  double sigma = total->GetParameter(4);     double sigmaErr = total->GetParError(4);

  TLatex* latexChi = CreateLatex(0.6, 0.8, TString::Format("#chi^{2}/NDF = %.2f / %d = %.2f", total->GetChisquare(), total->GetNDF(), total->GetChisquare() / total->GetNDF()).Data(), textSize);
  TLatex* latexMu = CreateLatex(0.6, 0.7, TString::Format("#mu = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * mu, 1e3 * muErr).Data(), textSize);
  TLatex* latexSigma = CreateLatex(0.6, 0.6, TString::Format("#sigma = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * sigma, 1e3 * sigmaErr).Data(), textSize);

  // double bkgEstimate = background->Integral(signalRegion[0], signalRegion[1]);
  double bkgEstimate = background_extrapolated->Integral(signalRegion[0], signalRegion[1]);
  double bkgErr = background->IntegralError(signalRegion[0], signalRegion[1], background->GetParameters(), bkgPtr->GetCovarianceMatrix().GetMatrixArray());
  double sigPlusBkgErr;
  double sigPlusBkg = mass->IntegralAndError(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]), sigPlusBkgErr, "width");
  double purity = 1 - bkgEstimate / sigPlusBkg;
  TLatex* latexS = CreateLatex(0.6, 0.5, TString::Format("S+B = %.1f #pm %.1f", sigPlusBkg, sigPlusBkgErr).Data(), textSize);
  TLatex* latexB = CreateLatex(0.6, 0.4, TString::Format("B = %.1f #pm %.1f", bkgEstimate, bkgErr).Data(), textSize);
  TLatex* latexPurity = CreateLatex(0.6, 0.3, TString::Format("Purity = %.2f %%", purity * 1e2).Data(), textSize);

  latexMu->Draw("same");
  latexSigma->Draw("same");
  latexS->Draw("same");
  latexB->Draw("same");
  latexChi->Draw("same");
  latexPurity->Draw("same");
  background_extrapolated->Draw("same");

  saveName = "K0S_purity";
  saveName = TString::Format("%s_pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  canvas->SaveAs(saveName.c_str());

  cout << TString::Format("#mu - 12 #sigma = %.8f GeV/#it{c}^{2}", (mu - 12 * sigma)) << endl;
  cout << TString::Format("#mu - 7 #sigma = %.8f GeV/#it{c}^{2}", (mu - 7 * sigma)) << endl;
  cout << TString::Format("#mu - 5 #sigma = %.8f GeV/#it{c}^{2}", (mu - 5 * sigma)) << endl;
  cout << TString::Format("#mu + 5 #sigma = %.8f GeV/#it{c}^{2}", (mu + 5 * sigma)) << endl;
  cout << TString::Format("#mu + 7 #sigma = %.8f GeV/#it{c}^{2}", (mu + 7 * sigma)) << endl;
  cout << TString::Format("#mu + 12 #sigma = %.8f GeV/#it{c}^{2}", (mu + 12 * sigma)) << endl;

  cout << "S+B count in signal region: " << sigPlusBkg << endl;
  cout << "B integral in signal region: " << bkgEstimate << endl;
  cout << "S estimate: " << sigPlusBkg - bkgEstimate << endl;
  cout << "Purity: " << purity << " = " << 1e2*purity << "%" << endl;
}
void Lambda0Purity(string inName = "AnalysisResults.root", double ptmin = 0., double ptmax = 100.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  const int nDim = 4;
  const int ptAxis = 0;
  const int K0SmassAxis = 1;
  const int Lambda0massAxis = 2;
  const int AntiLambda0massAxis = 3;
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 80e6;
  double xMinLegend = 0.25, xMaxLegend = 0.5, yMinLegend = 0.65, yMaxLegend = 0.85;
  double xLatex = 0.23, yLatex = 0.93;
  int xCanvas = 1800, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}_{%s} (GeV/#it{c}^{2})", formatHadronName("Lambda0").c_str()).Data();
  yTitle = "count";
  dataSet = "LHC23y_pass1";

  std::vector<TH1D*> histVector; std::vector<TF1*> funcVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "V0PtMass";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
  TH1D* mass = (TH1D*)thn->Projection(Lambda0massAxis);
  setStyle(mass, 0);
  legend->AddEntry(mass, "data");

  canvas->cd();
  frame->Draw();
  mass->Draw("same");

  double parameters[6];
  double signalRegion[2] = {1.10582782, 1.12532820}; // mu ± 5 sigma
  double fitRegion[2] = {1.08, 1.215};
  TF1* background_extrapolated = new TF1("background_extrapolated", "pol2", fitRegion[0], fitRegion[1]);
  // TF1* background = new TF1("background", "pol2", fitRegion[0], fitRegion[1]);
  TF1* background = new TF1("background", pol2bkg, fitRegion[0], fitRegion[1], 3);
  TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
  TF1* total = new TF1("total", "pol2(0) + gaus(3)", fitRegion[0], fitRegion[1]);

  setStyle(background, 1);
  TFitResultPtr bkgPtr = mass->Fit("background", "S");
  background->GetParameters(&parameters[0]);
  legend->AddEntry(background, "background");

  setStyle(background_extrapolated, 1);
  background_extrapolated->SetParameters(parameters);
  background_extrapolated->SetLineStyle(9);

  setStyle(signal, 2);
  mass->Fit("signal", "R+");
  signal->GetParameters(&parameters[3]);
  legend->AddEntry(signal, "signal");

  setStyle(total, 3);
  total->SetParameters(parameters);
  mass->Fit("total", "R+");
  legend->AddEntry(total, "combined");

  if (legend) legend->Draw("same");
  latexText = TString::Format("%s, %.0f < #it{p}_{T, V0} < %.0f GeV/#it{c}", dataSet.c_str(), ptmin, ptmax).Data();
  TLatex* latex = CreateLatex(0.3, 0.95, latexText, textSize);
  latex->Draw("same");

  double bkgOffset = total->GetParameter(0); double bkgOffsetErr = total->GetParError(0);
  double bkgSlope = total->GetParameter(1);  double bkgSlopeErr = total->GetParError(1);
  double bkgQuad = total->GetParameter(2);   double bkgQuadErr = total->GetParError(2);
  double norm = total->GetParameter(3);      double normErr = total->GetParError(3);
  double mu = total->GetParameter(4);        double muErr = total->GetParError(4);
  double sigma = total->GetParameter(5);     double sigmaErr = total->GetParError(5);

  TLatex* latexChi = CreateLatex(0.6, 0.8, TString::Format("#chi^{2}/NDF = %.2f / %d = %.2f", total->GetChisquare(), total->GetNDF(), total->GetChisquare() / total->GetNDF()).Data(), textSize);
  TLatex* latexMu = CreateLatex(0.6, 0.7, TString::Format("#mu = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * mu, 1e3 * muErr).Data(), textSize);
  TLatex* latexSigma = CreateLatex(0.6, 0.6, TString::Format("#sigma = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * sigma, 1e3 * sigmaErr).Data(), textSize);

  // double bkgEstimate = background->Integral(signalRegion[0], signalRegion[1]);
  double bkgEstimate = background_extrapolated->Integral(signalRegion[0], signalRegion[1]);
  double bkgErr = background->IntegralError(signalRegion[0], signalRegion[1], background->GetParameters(), bkgPtr->GetCovarianceMatrix().GetMatrixArray());
  double sigPlusBkgErr;
  double sigPlusBkg = mass->IntegralAndError(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]), sigPlusBkgErr, "width");
  double purity = 1 - bkgEstimate / sigPlusBkg;
  TLatex* latexS = CreateLatex(0.6, 0.5, TString::Format("S+B = %.1f #pm %.1f", sigPlusBkg, sigPlusBkgErr).Data(), textSize);
  TLatex* latexB = CreateLatex(0.6, 0.4, TString::Format("B = %.1f #pm %.1f", bkgEstimate, bkgErr).Data(), textSize);
  TLatex* latexPurity = CreateLatex(0.6, 0.3, TString::Format("Purity = %.2f %%", purity * 1e2).Data(), textSize);

  latexMu->Draw("same");
  latexSigma->Draw("same");
  latexS->Draw("same");
  latexB->Draw("same");
  latexChi->Draw("same");
  latexPurity->Draw("same");
  background_extrapolated->Draw("same");

  saveName = "Lambda0_purity";
  saveName = TString::Format("%s_pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  canvas->SaveAs(saveName.c_str());

  cout << TString::Format("#mu - 12 #sigma = %.8f GeV/#it{c}^{2}", (mu - 12 * sigma)) << endl;
  cout << TString::Format("#mu - 7 #sigma = %.8f GeV/#it{c}^{2}", (mu - 7 * sigma)) << endl;
  cout << TString::Format("#mu - 5 #sigma = %.8f GeV/#it{c}^{2}", (mu - 5 * sigma)) << endl;
  cout << TString::Format("#mu + 5 #sigma = %.8f GeV/#it{c}^{2}", (mu + 5 * sigma)) << endl;
  cout << TString::Format("#mu + 7 #sigma = %.8f GeV/#it{c}^{2}", (mu + 7 * sigma)) << endl;
  cout << TString::Format("#mu + 12 #sigma = %.8f GeV/#it{c}^{2}", (mu + 12 * sigma)) << endl;

  cout << "S+B count in signal region: " << sigPlusBkg << endl;
  cout << "B integral in signal region: " << bkgEstimate << endl;
  cout << "S estimate: " << sigPlusBkg - bkgEstimate << endl;
  cout << "Purity: " << purity << " = " << 1e2 * purity << "%" << endl;
}
void AntiLambda0Purity(string inName = "AnalysisResults.root", double ptmin = 0., double ptmax = 100.)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  const int nDim = 4;
  const int ptAxis = 0;
  const int K0SmassAxis = 1;
  const int Lambda0massAxis = 2;
  const int AntiLambda0massAxis = 3;
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 80e6;
  double xMinLegend = 0.25, xMaxLegend = 0.5, yMinLegend = 0.65, yMaxLegend = 0.85;
  double xLatex = 0.23, yLatex = 0.93;
  int xCanvas = 1800, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}_{%s} (cm)", formatHadronName("AntiLambda0").c_str()).Data();
  yTitle = "count";
  dataSet = "LHC23y_pass1";

  std::vector<TH1D*> histVector; std::vector<TF1*> funcVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "V0PtMass";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
  TH1D* mass = (TH1D*)thn->Projection(AntiLambda0massAxis);
  setStyle(mass, 0);
  legend->AddEntry(mass, "data");

  canvas->cd();
  frame->Draw();
  mass->Draw("same");

  double parameters[6];
  double signalRegion[2] = {1.10582782, 1.12532820}; // mu ± 5 sigma
  double fitRegion[2] = {1.08, 1.215};
  TF1* background_extrapolated = new TF1("background_extrapolated", "pol2", fitRegion[0], fitRegion[1]);
  // TF1* background = new TF1("background", "pol2", fitRegion[0], fitRegion[1]);
  TF1* background = new TF1("background", pol2bkg, fitRegion[0], fitRegion[1], 3);
  TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
  TF1* total = new TF1("total", "pol2(0) + gaus(3)", fitRegion[0], fitRegion[1]);

  setStyle(background, 1);
  TFitResultPtr bkgPtr = mass->Fit("background", "S");
  background->GetParameters(&parameters[0]);
  legend->AddEntry(background, "background");

  setStyle(background_extrapolated, 1);
  background_extrapolated->SetParameters(parameters);
  background_extrapolated->SetLineStyle(9);

  setStyle(signal, 2);
  mass->Fit("signal", "R+");
  signal->GetParameters(&parameters[3]);
  legend->AddEntry(signal, "signal");

  setStyle(total, 3);
  total->SetParameters(parameters);
  mass->Fit("total", "R+");
  legend->AddEntry(total, "combined");

  if (legend) legend->Draw("same");
  latexText = TString::Format("%s, %.0f < #it{p}_{T, V0} < %.0f GeV/#it{c}", dataSet.c_str(), ptmin, ptmax).Data();
  TLatex* latex = CreateLatex(0.3, 0.95, latexText, textSize);
  latex->Draw("same");

  double bkgOffset = total->GetParameter(0); double bkgOffsetErr = total->GetParError(0);
  double bkgSlope = total->GetParameter(1);  double bkgSlopeErr = total->GetParError(1);
  double bkgQuad = total->GetParameter(2);   double bkgQuadErr = total->GetParError(2);
  double norm = total->GetParameter(3);      double normErr = total->GetParError(3);
  double mu = total->GetParameter(4);        double muErr = total->GetParError(4);
  double sigma = total->GetParameter(5);     double sigmaErr = total->GetParError(5);

  TLatex* latexChi = CreateLatex(0.6, 0.8, TString::Format("#chi^{2}/NDF = %.2f / %d = %.2f", total->GetChisquare(), total->GetNDF(), total->GetChisquare() / total->GetNDF()).Data(), textSize);
  TLatex* latexMu = CreateLatex(0.6, 0.7, TString::Format("#mu = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * mu, 1e3 * muErr).Data(), textSize);
  TLatex* latexSigma = CreateLatex(0.6, 0.6, TString::Format("#sigma = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * sigma, 1e3 * sigmaErr).Data(), textSize);

  // double bkgEstimate = background->Integral(signalRegion[0], signalRegion[1]);
  double bkgEstimate = background_extrapolated->Integral(signalRegion[0], signalRegion[1]);
  double bkgErr = background->IntegralError(signalRegion[0], signalRegion[1], background->GetParameters(), bkgPtr->GetCovarianceMatrix().GetMatrixArray());
  double sigPlusBkgErr;
  double sigPlusBkg = mass->IntegralAndError(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]), sigPlusBkgErr, "width");
  double purity = 1 - bkgEstimate / sigPlusBkg;
  TLatex* latexS = CreateLatex(0.6, 0.5, TString::Format("S+B = %.1f #pm %.1f", sigPlusBkg, sigPlusBkgErr).Data(), textSize);
  TLatex* latexB = CreateLatex(0.6, 0.4, TString::Format("B = %.1f #pm %.1f", bkgEstimate, bkgErr).Data(), textSize);
  TLatex* latexPurity = CreateLatex(0.6, 0.3, TString::Format("Purity = %.2f %%", purity * 1e2).Data(), textSize);

  latexMu->Draw("same");
  latexSigma->Draw("same");
  latexS->Draw("same");
  latexB->Draw("same");
  latexChi->Draw("same");
  latexPurity->Draw("same");
  background_extrapolated->Draw("same");

  saveName = "AntiLambda0_purity";
  saveName = TString::Format("%s_pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  canvas->SaveAs(saveName.c_str());

  cout << TString::Format("#mu - 12 #sigma = %.8f GeV/#it{c}^{2}", (mu - 12 * sigma)) << endl;
  cout << TString::Format("#mu - 7 #sigma = %.8f GeV/#it{c}^{2}", (mu - 7 * sigma)) << endl;
  cout << TString::Format("#mu - 5 #sigma = %.8f GeV/#it{c}^{2}", (mu - 5 * sigma)) << endl;
  cout << TString::Format("#mu + 5 #sigma = %.8f GeV/#it{c}^{2}", (mu + 5 * sigma)) << endl;
  cout << TString::Format("#mu + 7 #sigma = %.8f GeV/#it{c}^{2}", (mu + 7 * sigma)) << endl;
  cout << TString::Format("#mu + 12 #sigma = %.8f GeV/#it{c}^{2}", (mu + 12 * sigma)) << endl;

  cout << "S+B count in signal region: " << sigPlusBkg << endl;
  cout << "B integral in signal region: " << bkgEstimate << endl;
  cout << "S estimate: " << sigPlusBkg - bkgEstimate << endl;
  cout << "Purity: " << purity << " = " << 1e2 * purity << "%" << endl;
}

// -------------------------------------------------------------------------------------------------

void cutVarPurityWrtNoCut(string inName = "", string hypothesis = "", int cutAxis = -1, double ptmin = 0., double ptmax = 100., bool doRatio = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if (cutAxis < 0) {
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
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 0.08;
  if ("K0S" == hypothesis) { xMinFrame = 0.4, xMaxFrame = 0.6, yMaxFrame = 0.1; }
  if (doRatio) { yMinFrame = 0.5, yMaxFrame = 1.2; }
  double xMinLegend = 0.25, xMaxLegend = 0.5, yMinLegend = 0.65, yMaxLegend = 0.85;
  double xLatex = 0.23, yLatex = 0.93;
  int xCanvas = 1800, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}_{%s} (cm)", formatHadronName(hypothesis).c_str()).Data();
  yTitle = "#frac{N}{N(uncut)}";
  dataSet = "LHC23y_pass1";

  std::vector<TH1D*> histVector; std::vector<TF1*> funcVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "V0CutVariation";
  histName = TString::Format("jet-fragmentation/data/V0/%s", histName.c_str()).Data();
  TFile* inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  // Apply pt selection for all histograms
  std::array<int,2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);

  TH1D* stdmass = (TH1D*)thn->Projection(projectionAxis);
  double stdmassIntegral = stdmass->Integral();
  stdmass->Scale(1./stdmassIntegral);
  setStyle(stdmass, 0);
  legend->AddEntry(stdmass, "No cut");
  latexText = TString::Format("%s, %.0f < #it{p}_{T, V0} < %.0f GeV/#it{c}", dataSet.c_str(), ptmin, ptmax).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  int underflowBin = 0, overflowBin = 1 + thn->GetAxis(cutAxis)->GetNbins();
  for (int iBin = underflowBin + 1; iBin <= overflowBin; iBin++) {
    THnSparseD* thn_copy = (THnSparseD*)thn->Clone("thn_copy");
    thn_copy->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]); // Just to be sure

    // Cut the axis off at the bottom (var > binLowEdge)
    thn_copy->GetAxis(cutAxis)->SetRange(iBin, overflowBin);
    if ( (ctauAxis == cutAxis) + (DCAdAxis == cutAxis) ) {
      // Cut the axis off at the top (var < binUpEdge)
      thn_copy->GetAxis(cutAxis)->SetRange(underflowBin, overflowBin - iBin);
    }
    TH1D* mass = (TH1D*)thn_copy->Projection(projectionAxis);
    mass->SetName(TString::Format("mass_axis%d_bin%d", cutAxis, iBin).Data());
    mass->Scale(1./stdmassIntegral);
    setStyle(mass, iBin);

    if (doRatio) { mass->Divide(stdmass); }
    histVector.push_back(mass);

    // Legend stuff
    string legendEntry;
    double binEdge = thn_copy->GetAxis(cutAxis)->GetBinLowEdge(iBin);
    if ( (ctauAxis == cutAxis) + (DCAdAxis == cutAxis) ) {
      binEdge = thn_copy->GetAxis(cutAxis)->GetBinUpEdge(iBin);
    }
    switch (cutAxis) {
      case RAxis:
        legendEntry = TString::Format("%s > %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case ctauAxis:
        legendEntry = TString::Format("%s < %.0f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case cosPAAxis:
        legendEntry = TString::Format("%s > %.3f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCApAxis:
        legendEntry = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCAnAxis:
        legendEntry = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCAdAxis:
        legendEntry = TString::Format("%s < %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
    }
    legend->AddEntry(mass, legendEntry.c_str());
  }
  if (doRatio) { stdmass->Divide(stdmass); }
  frame->Draw();
  stdmass->Draw("same hist");
  for (auto hist : histVector) {
    hist->Draw("same hist");
  }
  legend->Draw("same");
  latex->Draw("same");

  saveName = TString::Format("m%s_wrtNoCut", hypothesis.c_str()).Data();
  saveName = TString::Format("%s_%s", saveName.c_str(), axisNames[cutAxis].c_str()).Data();
  saveName = TString::Format("%s_%s", saveName.c_str(), doRatio ? "ratio" : "spectrum").Data();
  saveName = TString::Format("%s_v0pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax).Data();
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(saveName.c_str());
}

void cutVarPurityWrtStdCut(string inName = "", string hypothesis = "", int cutAxis = -1, double ptmin = 0., double ptmax = 100., bool doRatio = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hypothesis) + ("Lambda0" == hypothesis) + ("AntiLambda0" == hypothesis) != 1) {
    cout << "Error: hypothesis must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  if (cutAxis < 0) {
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

  // Bins corresponding to standard cuts. We don't cut on the first 4 axes, so -1 to ignore
  std::array<int, nDim> stdBins = {-1, -1, -1, -1, 3, 2 - 1*("K0S" == hypothesis), 3, 2, 2, 1};
  std::array<string, nDim> axisNames = {"pt", "K0Smass", "Lambda0mass", "AntiLambda0mass", "R", "ctau", "cosPA", "DCAp", "DCAn", "DCAd"};
  int projectionAxis = ("K0S" == hypothesis)*K0SmassAxis + ("Lambda0" == hypothesis)*Lambda0massAxis + ("AntiLambda0" == hypothesis)*AntiLambda0massAxis;
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool SetLogy = false;
  double xMinFrame = 1.015, xMaxFrame = 1.215, yMinFrame = 0., yMaxFrame = 0.08;
  if ("K0S" == hypothesis) {
    xMinFrame = 0.4, xMaxFrame = 0.6, yMaxFrame = 0.11;
    if (doRatio) { yMinFrame = 0.8, yMaxFrame = 1.3; }
  }
  else { if (doRatio) { yMinFrame = 0.9, yMaxFrame = 1.1; } }
  double xMinLegend = 0.25, xMaxLegend = 0.5, yMinLegend = 0.65, yMaxLegend = 0.85;
  double xLatex = 0.23, yLatex = 0.93;
  int xCanvas = 1800, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{M}_{%s} (cm)", formatHadronName(hypothesis).c_str()).Data();
  yTitle = "#frac{N}{N(std)}";
  dataSet = "LHC23y_pass1";

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

  thn->GetAxis(cutAxis)->SetRange(1, 1 + thn->GetAxis(cutAxis)->GetNbins());
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  latexText = TString::Format("%s, %.0f < #it{p}_{T, V0} < %.0f GeV/#it{c}", dataSet.c_str(), ptmin, ptmax).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  int underflowBin = 0, overflowBin = 1 + thn->GetAxis(cutAxis)->GetNbins();
  for (int iBin = underflowBin + 1; iBin <= overflowBin; iBin++) {
    THnSparseD* thn_copy = (THnSparseD*)thn->Clone("thn_copy");
    thn_copy->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]); // Just to be sure

    // Cut the axis off at the bottom (var > binLowEdge)
    thn_copy->GetAxis(cutAxis)->SetRange(iBin, overflowBin);
    if ( (ctauAxis == cutAxis) + (DCAdAxis == cutAxis) ) {
      // Cut the axis off at the top (var < binUpEdge)
      thn_copy->GetAxis(cutAxis)->SetRange(underflowBin, overflowBin - iBin);
    }
    TH1D* mass = (TH1D*)thn_copy->Projection(projectionAxis);
    mass->SetName(TString::Format("mass_axis%d_bin%d", cutAxis, iBin).Data());
    // mass->Scale(1./stdmassIntegral);
    setStyle(mass, iBin);

    // if (doRatio) { mass->Divide(stdmass); }
    histVector.push_back(mass);

    // Legend stuff
    string legendEntry;
    double binEdge = thn_copy->GetAxis(cutAxis)->GetBinLowEdge(iBin);
    if ( (ctauAxis == cutAxis) + (DCAdAxis == cutAxis) ) {
      binEdge = thn_copy->GetAxis(cutAxis)->GetBinUpEdge(iBin);
    }
    switch (cutAxis) {
      case RAxis:
        legendEntry = TString::Format("%s > %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case ctauAxis:
        legendEntry = TString::Format("%s < %.0f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case cosPAAxis:
        legendEntry = TString::Format("%s > %.3f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCApAxis:
        legendEntry = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCAnAxis:
        legendEntry = TString::Format("%s > %.2f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
      case DCAdAxis:
        legendEntry = TString::Format("%s < %.1f", axisNames[cutAxis].c_str(), binEdge).Data();
        break;
    }
    legend->AddEntry(mass, legendEntry.c_str());
    if (iBin == stdBins[cutAxis]) {
      yTitle = TString::Format("#frac{N}{N(%s)}", legendEntry.c_str()).Data();
      setStyle(mass, 0);
    }
  }
  TH1D* stdmass = (TH1D*)histVector[stdBins[cutAxis]-1]->Clone("stdmass");
  double stdmassIntegral = stdmass->Integral();
  for (auto hist : histVector) {
  if (doRatio) { hist->Divide(stdmass); }
  else { hist->Scale(1./stdmassIntegral); }
  }

  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->Draw();
  for (auto hist : histVector) {
    hist->Draw("same hist");
  }
  legend->Draw("same");
  latex->Draw("same");

  saveName = TString::Format("m%s_wrtStdCut", hypothesis.c_str()).Data();
  saveName = TString::Format("%s_%s", saveName.c_str(), axisNames[cutAxis].c_str()).Data();
  saveName = TString::Format("%s_%s", saveName.c_str(), doRatio ? "ratio" : "spectrum").Data();
  saveName = TString::Format("%s_v0pt%.1f-%.1f", saveName.c_str(), ptmin, ptmax).Data();
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(saveName.c_str());
}

