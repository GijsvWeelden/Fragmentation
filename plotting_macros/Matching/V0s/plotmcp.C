
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
  mass->Rebin(4);
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

string getDataSet(int train)
{
  if (210373 == train) return "LHC24b1b";
  if (271952 == train) return "LHC24b1b";
  if (280432 == train) return "LHC24g4";

  return "Could not find dataset";
}

void plotTrain(int train, string hadron, double jetptmin, double jetptmax, double v0min, double v0max, bool doZ, int setting)
{
  string inputName = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  string dataSet   = getDataSet(train);

  switch (setting) {
    case 0:
      ptResolution(inputName, dataSet, jetptmin, jetptmax, v0min, v0max);
      break;
    case 1:
      {
        dauPtResolution(inputName, dataSet, jetptmin, jetptmax, v0min, v0max, doZ);
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
        string hypothesis = "K0S";
        vector<string> inputStrings = {inputName, dataSet, hadron, hypothesis};
        mass(inputStrings, jetptmin, jetptmax, v0min, v0max, doZ, false /* normalise */);
        mass(inputStrings, jetptmin, jetptmax, v0min, v0max, doZ, true);
      }
      break;
    case 4:
      {
        string hypothesis = "Lambda0";
        vector<string> inputStrings = {inputName, dataSet, hadron, hypothesis};
        mass(inputStrings, jetptmin, jetptmax, v0min, v0max, doZ, false /* normalise */);
        mass(inputStrings, jetptmin, jetptmax, v0min, v0max, doZ, true);
      }
      break;
    case 5:
      {
        string hypothesis = "AntiLambda0";
        vector<string> inputStrings = {inputName, dataSet, hadron, hypothesis};
        mass(inputStrings, jetptmin, jetptmax, v0min, v0max, doZ, false /* normalise */);
        mass(inputStrings, jetptmin, jetptmax, v0min, v0max, doZ, true);
      }
      break;
    default:
      cout << "Error: invalid setting" << endl;
      return;
  }
}
void plotTrain(int train, string hadron, double jetptmin, double jetptmax, bool doZ, int setting)
{
  gROOT->SetBatch();
  vector<double> pt   = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0};
  vector<double> z    = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  vector<double> bins = (doZ) ? z : pt;
  for (int i = 0; i < bins.size() - 1; i++) {
    if (bins[i] > jetptmax) break;
    plotTrain(train, hadron, jetptmin, jetptmax, bins[i], bins[i+1], doZ, setting);
  }
}

void plot280432(string hadron, double jetptmin, double jetptmax, double v0ptmin, double v0ptmax, bool doZ, int setting)
{
  plotTrain(280432, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, setting);
}
void plot280432(string hadron, double jetptmin, double jetptmax, bool doZ, int setting)
{
  plotTrain(280432, hadron, jetptmin, jetptmax, doZ, setting);
}


void plot271952(string hadron, double jetptmin, double jetptmax, double v0ptmin, double v0ptmax, bool doZ, int setting)
{
  plotTrain(271952, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, setting);
}
void plot271952(string hadron, double jetptmin, double jetptmax, bool doZ, int setting)
{
  plotTrain(271952, hadron, jetptmin, jetptmax, doZ, setting);
}

void plot210373(string hadron, double jetptmin, double jetptmax, double v0ptmin, double v0ptmax, bool doZ, int setting)
{
  plotTrain(210373, hadron, jetptmin, jetptmax, v0ptmin, v0ptmax, doZ, setting);
}
