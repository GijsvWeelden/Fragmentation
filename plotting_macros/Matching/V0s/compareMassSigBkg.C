
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
// -------------------------------------------------------------------------------------------------
// -----------------          Comparing mass of matched vs unmatched V0s           -----------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

THnSparse* getMatchedMass(vector<string> inputStrings, bool doZ)
{
  string inName = inputStrings[0];
  string hadron = inputStrings[2];
  string histName = "jet-fragmentation/matching/jets/V0/";
  if (doZ) histName += TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProjMass%s", hadron.c_str(), hadron.c_str(), hadron.c_str());
  else histName += TString::Format("partJetPt%sPtDetJetPt%sPtMass%s", hadron.c_str(), hadron.c_str(), hadron.c_str());
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  return thn;
}
THnSparse* getMatchedMass(vector<string> inputStrings)
{
  string inName = inputStrings[0];
  string hadron = inputStrings[2];
  if (hadron == "Lambda0") hadron = "Lambda";
  if (hadron == "AntiLambda0") hadron = "antiLambda";
  string histName = "jet-fragmentation/matching/V0/";
  histName += hadron;
  histName += "PtCtauMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  return thn;
}
THnSparse* getFakeMass(vector<string> inputStrings, bool doZ)
{
  string inName = inputStrings[0];
  string histName = "jet-fragmentation/matching/jets/V0/";
  if (doZ) histName += "fakeJetPtV0TrackProjMass";
  else histName += "fakeJetPtV0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  return thn;
}
THnSparse* getFakeMass(vector<string> inputStrings)
{
  string inName = inputStrings[0];
  string histName = "jet-fragmentation/matching/V0/";
  histName += "fakeV0PtMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  return thn;
}

void comparemass(vector<string> inputStrings, double v0min, double v0max, bool normalise)
{
  double textSize = 0.04;
  gStyle->SetNdivisions(505, "xy");
  string inName     = inputStrings[0];
  string dataSet    = inputStrings[1];
  string hadron     = inputStrings[2];

  const int partV0PtAxis    = 0;
  const int detV0PtAxis     = 1;
  const int ctauAxis        = 2;
  const int massAxis        = 3;

  const int K0SAxis         = 1;
  const int Lambda0Axis     = 2;
  const int AntiLambda0Axis = 3;

  THnSparseD* thn_matchedMass = (THnSparseD*) getMatchedMass(inputStrings);
  THnSparseD* thn_fakeMass    = (THnSparseD*) getFakeMass(inputStrings);

  array<int, 2> partv0bins = getProjectionBins(thn_matchedMass->GetAxis(partV0PtAxis), v0min, v0max, 2e-4);
  thn_matchedMass->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);
  double lowv0  = thn_matchedMass->GetAxis(partV0PtAxis)->GetBinLowEdge(partv0bins[0]);
  double highv0 = thn_matchedMass->GetAxis(partV0PtAxis)->GetBinUpEdge(partv0bins[1]);

  partv0bins = getProjectionBins(thn_fakeMass->GetAxis(partV0PtAxis), v0min, v0max, 2e-4);
  thn_fakeMass->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);

  int projectionAxis = (hadron == "K0S")*K0SAxis + (hadron == "Lambda0")*Lambda0Axis + (hadron == "AntiLambda0")*AntiLambda0Axis;
  TH1D* matchedMass = (TH1D*)thn_matchedMass->Projection(massAxis);
  TH1D* fakeMass = (TH1D*)thn_fakeMass->Projection(projectionAxis);
  if (isHistEmptyInRange(matchedMass, 1, matchedMass->GetNbinsX())
      || isHistEmptyInRange(fakeMass, 1, fakeMass->GetNbinsX())) {
    cout << "compareMass: Empty mass histogram for ptv0 " << lowv0 << " - " << highv0 << endl;
    return;
  }
  setStyle(matchedMass, 0);
  setStyle(fakeMass, 1);

  TLegend* legend = CreateLegend(0.6, 0.9, 0.55, 0.75, "", textSize);
  legend->AddEntry(matchedMass, "Signal");
  legend->AddEntry(fakeMass, "Bkg");

  string saveName = hadron;
  saveName += "_massComparison";
  saveName += TString::Format("_partv0pt%.1f-%.1f", lowv0, highv0).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  matchedMass->Rebin(4);
  fakeMass->Rebin(4);
  if (normalise) {
    double scale = TMath::Max(getHistScale(matchedMass, false), getHistScale(fakeMass, false));
    matchedMass->Scale(1./scale, "width");
    fakeMass->Scale(1./scale, "width");
  }
  matchedMass->SetName(("matched" + saveName).c_str());
  fakeMass->SetName(("fake" + saveName).c_str());

  string xTitle = "M(" + formatHadronDaughters(hadron) + ") (GeV/#it{c}^{2})";
  string yTitle = "Counts";
  if (normalise) yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
  double xMinFrame = matchedMass->GetXaxis()->GetXmin(), xMaxFrame = matchedMass->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.5 * TMath::Max(getHistScale(matchedMass, true), getHistScale(fakeMass, true));
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  string v0Text  = TString::Format("#it{p}_{T, V0} = %.1f-%.1f GeV/#it{c}", lowv0, highv0).Data();
  frame->SetTitle((dataSet + ", " + v0Text).c_str());

  double mass = ("K0S" == hadron) ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, 0.75 * yMaxFrame);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(2));
  line->SetLineWidth(3);

  canvas->cd();
  frame->Draw();
  line->Draw("same");
  legend->Draw("same");
  matchedMass->Draw("same");
  fakeMass->Draw("same");
  canvas->SaveAs(saveName.c_str());
}
void comparemass(vector<string> inputStrings, double partjetptmin, double partjetptmax, double v0min, double v0max, bool doZ, bool normalise)
{
  double textSize = 0.04;
  gStyle->SetNdivisions(505, "xy");
  string inName     = inputStrings[0];
  string dataSet    = inputStrings[1];
  string hadron     = inputStrings[2];

  const int partJetPtAxis   = 0;
  const int partV0PtAxis    = 1;
  const int detJetPtAxis    = 2;
  const int detV0PtAxis     = 3;
  const int massAxis        = 4;

  const int K0SAxis         = 2;
  const int Lambda0Axis     = 3;
  const int AntiLambda0Axis = 4;

  THnSparseD* thn_matchedMass = (THnSparseD*) getMatchedMass(inputStrings, doZ);
  THnSparseD* thn_fakeMass    = (THnSparseD*) getFakeMass(inputStrings, doZ);

  array<int, 2> jetptbins = getProjectionBins(thn_matchedMass->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  array<int, 2> partv0bins = getProjectionBins(thn_matchedMass->GetAxis(partV0PtAxis), v0min, v0max, 2e-4);
  thn_matchedMass->GetAxis(partJetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn_matchedMass->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);
  double lowjetpt  = thn_matchedMass->GetAxis(partJetPtAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn_matchedMass->GetAxis(partJetPtAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0     = thn_matchedMass->GetAxis(partV0PtAxis)->GetBinLowEdge(partv0bins[0]);
  double highv0    = thn_matchedMass->GetAxis(partV0PtAxis)->GetBinUpEdge(partv0bins[1]);

  jetptbins  = getProjectionBins(thn_fakeMass->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  partv0bins = getProjectionBins(thn_fakeMass->GetAxis(partV0PtAxis), v0min, v0max, 2e-4);
  thn_fakeMass->GetAxis(partJetPtAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn_fakeMass->GetAxis(partV0PtAxis)->SetRange(partv0bins[0], partv0bins[1]);

  int projectionAxis = (hadron == "K0S")*K0SAxis + (hadron == "Lambda0")*Lambda0Axis + (hadron == "AntiLambda0")*AntiLambda0Axis;
  TH1D* matchedMass = (TH1D*)thn_matchedMass->Projection(massAxis);
  TH1D* fakeMass = (TH1D*)thn_fakeMass->Projection(projectionAxis);
  if (isHistEmptyInRange(matchedMass, 1, matchedMass->GetNbinsX())
      || isHistEmptyInRange(fakeMass, 1, fakeMass->GetNbinsX())) {
    cout << "compareMass: Empty mass histogram for jetpt " << lowjetpt << " - " << highjetpt << ", " << (doZ ? "zv0 " : "ptv0 ") << lowv0 << " - " << highv0 << endl;
    return;
  }
  setStyle(matchedMass, 0);
  setStyle(fakeMass, 1);

  TLegend* legend = CreateLegend(0.6, 0.9, 0.55, 0.75, "", textSize);
  legend->AddEntry(matchedMass, "Signal");
  legend->AddEntry(fakeMass, "Bkg");

  string saveName = hadron;
  saveName += "_massComparison";
  saveName += TString::Format("_partjetpt%.0f-%.0f", lowjetpt, highjetpt).Data();
  if (doZ) saveName += TString::Format("_partv0z%.2f-%.2f", lowv0, highv0).Data();
  else saveName += TString::Format("_partv0pt%.1f-%.1f", lowv0, highv0).Data();
  if (!normalise) saveName += "_counts";
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  matchedMass->Rebin(4);
  fakeMass->Rebin(4);
  if (normalise) {
    double scale = TMath::Max(getHistScale(matchedMass, false), getHistScale(fakeMass, false));
    matchedMass->Scale(1./scale, "width");
    fakeMass->Scale(1./scale, "width");
  }
  matchedMass->SetName(("matched" + saveName).c_str());
  fakeMass->SetName(("fake" + saveName).c_str());

  string xTitle = "M(" + formatHadronDaughters(hadron) + ") (GeV/#it{c}^{2})";
  string yTitle = "Counts";
  if (normalise) yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
  double xMinFrame = matchedMass->GetXaxis()->GetXmin(), xMaxFrame = matchedMass->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.5 * TMath::Max(getHistScale(matchedMass, true), getHistScale(fakeMass, true));
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  string jetText = TString::Format("#it{p}_{T, ch+V0 jet} = %.0f-%.0f GeV/#it{c}", lowjetpt, highjetpt).Data();
  string v0Text  = TString::Format("#it{p}_{T, V0} = %.1f-%.1f GeV/#it{c}", lowv0, highv0).Data();
  if (doZ) { v0Text = TString::Format("#it{z}_{V0} = %.2f-%.2f", lowv0, highv0).Data(); }
  string latexText = TString::Format("#splitline{%s}{%s}", jetText.c_str(), v0Text.c_str()).Data();
  double xLatex = 0.3, yLatex = 0.8;
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  double mass = ("K0S" == hadron) ? MassK0S : MassLambda0;
  TLine* line = new TLine(mass, yMinFrame, mass, 0.75 * yMaxFrame);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(2));
  line->SetLineWidth(3);

  canvas->cd();
  frame->Draw();
  line->Draw("same");
  latex->Draw("same");
  legend->Draw("same");
  matchedMass->Draw("same");
  fakeMass->Draw("same");
  canvas->SaveAs(saveName.c_str());
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

void plotTrain(string train, string dataSet, string hadron, double jetptmin, double jetptmax, double v0ptmin, double v0ptmax, bool doZ, int setting)
{
  string inputName = "~/cernbox/TrainOutput/" + train + "/AnalysisResults.root";
  vector<string> inputStrings = {inputName, dataSet, hadron};
  switch (setting) {
    case 0:
      comparemass(inputStrings, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, false);
      comparemass(inputStrings, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, true);
      break;
    case 1:
      comparemass(inputStrings, v0ptmin, v0ptmax, false);
      comparemass(inputStrings, v0ptmin, v0ptmax, true);
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
