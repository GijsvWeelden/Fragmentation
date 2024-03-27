
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

#include "/Users/gijsvanweelden/Documents/Fragmentation/plotting_macros/histUtils.C"

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
double pol2bkg(double *x, double *par)
{
  if (x[0] > par[0] && x[0] < par[1]) {
    TF1::RejectPoint();
    return 0;
  }
  return par[2] + par[3]*x[0] + par[4]*x[0]*x[0];
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
  xTitle = TString::Format("#it{M}_{%s} (cm)", formatHadronName("K0S").c_str()).Data();
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
  setStyle(mass, 0);
  legend->AddEntry(mass, "data");

  canvas->cd();
  frame->Draw();
  mass->Draw("same");

  double parameters[5];
  double signalRegion[2] = {0.47209646, 0.52148143}; // mu ± 5 sigma
  TF1* background = new TF1("background", "pol1", 0.4, 0.6);
  TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
  TF1* total = new TF1("total", "pol1(0) + gaus(2)", 0.4, 0.6);

  setStyle(background, 1);
  TFitResultPtr bkgPtr = mass->Fit("background", "S");
  background->GetParameters(&parameters[0]);
  legend->AddEntry(background, "background");

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

  TLatex* latexChi = CreateLatex(0.6, 0.4, TString::Format("#chi^{2}/NDF = %.2f / %d = %.2f", total->GetChisquare(), total->GetNDF(), total->GetChisquare() / total->GetNDF()).Data(), textSize);
  TLatex* latexMu = CreateLatex(0.6, 0.5, TString::Format("#mu = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * mu, 1e3 * muErr).Data(), textSize);
  TLatex* latexSigma = CreateLatex(0.6, 0.6, TString::Format("#sigma = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * sigma, 1e3 * sigmaErr).Data(), textSize);

  double bkgEstimate = background->Integral(signalRegion[0], signalRegion[1]);
  double bkgErr = background->IntegralError(signalRegion[0], signalRegion[1], background->GetParameters(), bkgPtr->GetCovarianceMatrix().GetMatrixArray());
  double sigPlusBkgErr;
  double sigPlusBkg = mass->IntegralAndError(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]), sigPlusBkgErr, "width");
  double purity = 1 - bkgEstimate / sigPlusBkg;
  TLatex* latexS = CreateLatex(0.6, 0.8, TString::Format("S+B = %.1f #pm %.1f", sigPlusBkg, sigPlusBkgErr).Data(), textSize);
  TLatex* latexB = CreateLatex(0.6, 0.7, TString::Format("B = %.1f #pm %.1f", bkgEstimate, bkgErr).Data(), textSize);

  latexMu->Draw("same");
  latexSigma->Draw("same");
  latexS->Draw("same");
  latexB->Draw("same");
  latexChi->Draw("same");

  saveName = "K0S_purity.pdf";
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
  xTitle = TString::Format("#it{M}_{%s} (cm)", formatHadronName("Lambda0").c_str()).Data();
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
  TF1* background = new TF1("background", "pol2", 1.08, 1.215);
  TF1* fBKG = new TF1("fBKG", pol2bkg, 1.08, 1.215, 3);
  TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
  TF1* total = new TF1("total", "pol2(0) + gaus(3)", 1.08, 1.215);

  setStyle(background, 1);
  TFitResultPtr bkgPtr = mass->Fit("background", "S");
  background->GetParameters(&parameters[0]);
  legend->AddEntry(background, "background");

  // setStyle(fBKG, 2);
  // mass->Fit("fBKG", "R+");
  // legend->AddEntry(fBKG, "fBKG");
  // return;

  setStyle(signal, 2);
  mass->Fit("signal", "R+");
  signal->GetParameters(&parameters[3]);
  legend->AddEntry(signal, "signal");
  // return;

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

  TLatex* latexChi = CreateLatex(0.6, 0.4, TString::Format("#chi^{2}/NDF = %.2f / %d = %.2f", total->GetChisquare(), total->GetNDF(), total->GetChisquare() / total->GetNDF()).Data(), textSize);
  TLatex* latexMu = CreateLatex(0.6, 0.5, TString::Format("#mu = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * mu, 1e3 * muErr).Data(), textSize);
  TLatex* latexSigma = CreateLatex(0.6, 0.6, TString::Format("#sigma = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * sigma, 1e3 * sigmaErr).Data(), textSize);

  double bkgEstimate = background->Integral(signalRegion[0], signalRegion[1]);
  double bkgErr = background->IntegralError(signalRegion[0], signalRegion[1], background->GetParameters(), bkgPtr->GetCovarianceMatrix().GetMatrixArray());
  double sigPlusBkgErr;
  double sigPlusBkg = mass->IntegralAndError(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]), sigPlusBkgErr, "width");
  double purity = 1 - bkgEstimate / sigPlusBkg;
  TLatex* latexS = CreateLatex(0.6, 0.8, TString::Format("S+B = %.1f #pm %.1f", sigPlusBkg, sigPlusBkgErr).Data(), textSize);
  TLatex* latexB = CreateLatex(0.6, 0.7, TString::Format("B = %.1f #pm %.1f", bkgEstimate, bkgErr).Data(), textSize);

  latexMu->Draw("same");
  latexSigma->Draw("same");
  latexS->Draw("same");
  latexB->Draw("same");
  latexChi->Draw("same");

  saveName = "Lambda0_purity.pdf";
  canvas->SaveAs(saveName.c_str());

  cout << TString::Format("#mu - 12 #sigma = %.8f GeV/#it{c}^{2}", (mu - 12 * sigma)) << endl;
  cout << TString::Format("#mu - 7 #sigma = %.8f GeV/#it{c}^{2}", (mu - 7 * sigma)) << endl;
  cout << TString::Format("#mu - 5 #sigma = %.8f GeV/#it{c}^{2}", (mu - 5 * sigma)) << endl;
  cout << TString::Format("#mu + 5 #sigma = %.8f GeV/#it{c}^{2}", (mu + 5 * sigma)) << endl;
  cout << TString::Format("#mu + 7 #sigma = %.8f GeV/#it{c}^{2}", (mu + 7 * sigma)) << endl;
  cout << TString::Format("#mu + 12 #sigma = %.8f GeV/#it{c}^{2}", (mu + 12 * sigma)) << endl;

  // double bkgEstimate = background->Integral(signalRegion[0], signalRegion[1]);
  // double sigPlusBkg = mass->Integral(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]));
  // double purity = 1 - bkgEstimate / sigPlusBkg;

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
  TF1* background = new TF1("background", "pol2", 1.08, 1.215);
  TF1* fBKG = new TF1("fBKG", pol2bkg, 1.08, 1.215, 3);
  TF1* signal = new TF1("signal", "gaus", signalRegion[0], signalRegion[1]);
  TF1* total = new TF1("total", "pol2(0) + gaus(3)", 1.08, 1.215);

  setStyle(background, 1);
  TFitResultPtr bkgPtr = mass->Fit("background", "S");
  background->GetParameters(&parameters[0]);
  legend->AddEntry(background, "background");

  // setStyle(fBKG, 2);
  // mass->Fit("fBKG", "R+");
  // legend->AddEntry(fBKG, "fBKG");
  // return;

  setStyle(signal, 2);
  mass->Fit("signal", "R+");
  signal->GetParameters(&parameters[3]);
  legend->AddEntry(signal, "signal");
  // return;

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

  TLatex* latexChi = CreateLatex(0.6, 0.4, TString::Format("#chi^{2}/NDF = %.2f / %d = %.2f", total->GetChisquare(), total->GetNDF(), total->GetChisquare() / total->GetNDF()).Data(), textSize);
  TLatex* latexMu = CreateLatex(0.6, 0.5, TString::Format("#mu = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * mu, 1e3 * muErr).Data(), textSize);
  TLatex* latexSigma = CreateLatex(0.6, 0.6, TString::Format("#sigma = %.3f #pm %.3f MeV/#it{c}^{2}", 1e3 * sigma, 1e3 * sigmaErr).Data(), textSize);

  double bkgEstimate = background->Integral(signalRegion[0], signalRegion[1]);
  double bkgErr = background->IntegralError(signalRegion[0], signalRegion[1], background->GetParameters(), bkgPtr->GetCovarianceMatrix().GetMatrixArray());
  double sigPlusBkgErr;
  double sigPlusBkg = mass->IntegralAndError(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]), sigPlusBkgErr, "width");
  double purity = 1 - bkgEstimate / sigPlusBkg;
  TLatex* latexS = CreateLatex(0.6, 0.8, TString::Format("S+B = %.1f #pm %.1f", sigPlusBkg, sigPlusBkgErr).Data(), textSize);
  TLatex* latexB = CreateLatex(0.6, 0.7, TString::Format("B = %.1f #pm %.1f", bkgEstimate, bkgErr).Data(), textSize);

  latexMu->Draw("same");
  latexSigma->Draw("same");
  latexS->Draw("same");
  latexB->Draw("same");
  latexChi->Draw("same");

  saveName = "AntiLambda0_purity.pdf";
  canvas->SaveAs(saveName.c_str());

  cout << TString::Format("#mu - 12 #sigma = %.8f GeV/#it{c}^{2}", (mu - 12 * sigma)) << endl;
  cout << TString::Format("#mu - 7 #sigma = %.8f GeV/#it{c}^{2}", (mu - 7 * sigma)) << endl;
  cout << TString::Format("#mu - 5 #sigma = %.8f GeV/#it{c}^{2}", (mu - 5 * sigma)) << endl;
  cout << TString::Format("#mu + 5 #sigma = %.8f GeV/#it{c}^{2}", (mu + 5 * sigma)) << endl;
  cout << TString::Format("#mu + 7 #sigma = %.8f GeV/#it{c}^{2}", (mu + 7 * sigma)) << endl;
  cout << TString::Format("#mu + 12 #sigma = %.8f GeV/#it{c}^{2}", (mu + 12 * sigma)) << endl;

  // double bkgEstimate = background->Integral(signalRegion[0], signalRegion[1]);
  // double sigPlusBkg = mass->Integral(mass->FindBin(signalRegion[0]), mass->FindBin(signalRegion[1]));
  // double purity = 1 - bkgEstimate / sigPlusBkg;

  cout << "S+B count in signal region: " << sigPlusBkg << endl;
  cout << "B integral in signal region: " << bkgEstimate << endl;
  cout << "S estimate: " << sigPlusBkg - bkgEstimate << endl;
  cout << "Purity: " << purity << " = " << 1e2 * purity << "%" << endl;
}

