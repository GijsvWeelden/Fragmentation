
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
 * Plots for inclusive V0s (not necessarily in jets)
 */

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
// -------------------------------------------------------------------------------------------------

void ptResolution(string inName, string dataSet, double partv0min, double partv0max)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }

  string hadron = "V0";
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  histName  = "jet-fragmentation/matching/V0/";
  histName += "V0PartPtRatioPtRelDiffPt";
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  std::array<int, 2> partv0bins = getProjectionBins(th3->GetXaxis(), partv0min, partv0max);
  TH1D* v0res = (TH1D*)th3->ProjectionZ("v0resolution", partv0bins[0], partv0bins[1], 0, th3->GetNbinsY() + 1);
  v0res->Scale(1./v0res->Integral());
  setStyle(v0res, 0);

  double lowv0 = th3->GetXaxis()->GetBinLowEdge(partv0bins[0]);
  double highv0 = th3->GetXaxis()->GetBinUpEdge(partv0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  // Get RMS = V0 resolution for given pt range
  double rms = 0.;
  for (int i = 1; i < v0res->GetNbinsX(); i++) {
    double binContent = v0res->GetBinContent(i);
    double binCenter = v0res->GetBinCenter(i);
    // Hist is self-normalised, so no need to divide by integral
    rms += binContent * binCenter * binCenter;
  }
  rms = TMath::Sqrt(rms);
  // v0res->Print("all");

  double xLatex = 0.4, yLatex = 0.8;
  string ptText  = TString::Format("#it{p}_{T, V0}^{part.} = %.1f - %.1f GeV/#it{c}", lowv0, highv0).Data();
  string rmsText = TString::Format("RMS: %.2f", rms).Data();
  latexText = TString::Format("#splitline{%s}{#splitline{%s}{%s}}", dataSet.c_str(), ptText.c_str(), rmsText.c_str()).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "V0PtResolution";
  saveName += TString::Format("_partv0pt%.1f-%.1f", lowv0, highv0);
  saveName += ".pdf";
  int xCanvas = 900, yCanvas = 900;
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);
  canvas->SetLogy();

  double xMinFrame = -1., xMaxFrame = 10., yMinFrame = 1e-8, yMaxFrame = 2.;
  xTitle = "(#it{p}_{T, V0}^{det.} - #it{p}_{T, V0}^{part.})/#it{p}_{T, V0}^{part.}";
  yTitle = "normalised count";
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);

  canvas->cd();
  frame->Draw();
  v0res->Draw("same");
  latex->Draw("same");
  canvas->SaveAs(canvas->GetName());
}
void dauPtResolution(string inName, string dataSet, double partmin, double partmax, bool doPos)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  gStyle->SetNdivisions(505, "xy");
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  string histName  = "jet-fragmentation/matching/V0/";
  if (doPos) histName += "V0PosPartPtRatioPtRelDiffPt";
  else histName += "V0NegPartPtRatioPtRelDiffPt";
  TFile *inFile = TFile::Open(inName.c_str());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  std::array<int, 2> ptbins = getProjectionBins(th3->GetXaxis(), partmin, partmax);
  TH1D* daures = (TH1D*)th3->ProjectionZ("resolution", ptbins[0], ptbins[1], 0, th3->GetNbinsY() + 1);
  daures->Scale(1./daures->Integral());
  setStyle(daures, 0);

  double lowpt = th3->GetXaxis()->GetBinLowEdge(ptbins[0]);
  double highpt = th3->GetXaxis()->GetBinUpEdge(ptbins[1]);
  if (lowpt < 0.) { lowpt = 0.; } // Avoids ugly pt>-0 in latextext

  // Get RMS = daughter resolution for given pt range
  double rms = getHistRMS(daures);

  double xLatex = 0.4, yLatex = 0.8;
  string ptText  = TString::Format("#it{p}_{T, %s dau} = %.1f - %.1f GeV/#it{c}", (doPos ? "pos" : "neg"), lowpt, highpt).Data();
  string rmsText = TString::Format("RMS: %.2f", rms).Data();
  string latexText = TString::Format("#splitline{%s}{#splitline{%s}{%s}}", dataSet.c_str(), ptText.c_str(), rmsText.c_str()).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  string saveName = "V0";
  saveName += (doPos ? "Pos" : "Neg");
  saveName += "PtResolution";
  saveName += TString::Format("_pt%.1f-%.1f", lowpt, highpt);
  saveName += ".pdf";

  int xCanvas = 900, yCanvas = 900;
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);
  canvas->SetLogy();

  double xMinFrame = -1., xMaxFrame = 10., yMinFrame = 1e-8, yMaxFrame = 2.;
  string xTitle = "(#it{p}_{T}^{det.} - #it{p}_{T}^{part.})/#it{p}_{T}^{part.}";
  string yTitle = "normalised count";
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);

  canvas->cd();
  frame->Draw();
  daures->Draw("same");
  latex->Draw("same");
  canvas->SaveAs(canvas->GetName());
}

void matchedPt(string inName = "", bool detector = false)
{
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 60., yMinFrame = 1e-8, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, V0}^{part.} (GeV/#it{c})";
  if (detector) { xTitle = "#it{p}_{T, V0}^{det.} (GeV/#it{c})"; }
  yTitle = "normalised count";
  dataSet = "LHC24b1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/V0/V0PartPtDetPt";
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());
  TH1D* matchedjetpt = (TH1D*)th2->ProjectionX("matchedjetpt", 0, th2->GetNbinsY()+1);
  if (detector) { matchedjetpt = (TH1D*)th2->ProjectionY("matchedjetpt", 0, th2->GetNbinsX()+1);}
  matchedjetpt->Scale(1./ matchedjetpt->Integral());
  setStyle(matchedjetpt, 0);
  histVector.push_back(matchedjetpt);

  latexText = TString::Format("%s, matched V0", dataSet.c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedPartV0Pt";
  if (detector) { saveName = "matchedDetV0Pt"; }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}

void mass(vector<string> inputStrings, double v0min, double v0max, bool normalise)
{
  gStyle->SetNdivisions(505, "xy");
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  string hadron  = inputStrings[2];

  const int nDim         = 4;
  const int partV0PtAxis = 0;
  const int detV0PtAxis  = 1;
  const int ctauAxis     = 2;
  const int massAxis     = 3;

  string histName = "jet-fragmentation/matching/V0/";
  histName += hadron;
  histName += "PtCtauMass";
  TFile *inFile = TFile::Open(inName.c_str());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  // To work well with file name and formatting conventions
  if (hadron == "Lambda") hadron = "Lambda0";
  if (hadron == "antiLambda") hadron = "AntiLambda0";

  array<int, 2> ptbins = getProjectionBins(thn->GetAxis(partV0PtAxis), v0min, v0max);
  thn->GetAxis(partV0PtAxis)->SetRange(ptbins[0], ptbins[1]);
  double lowpt = thn->GetAxis(partV0PtAxis)->GetBinLowEdge(ptbins[0]);
  double highpt = thn->GetAxis(partV0PtAxis)->GetBinUpEdge(ptbins[1]);
  if (lowpt < 0.) { lowpt = 0.; } // Avoids ugly pt>-0 in latextext

  string saveName = hadron;
  saveName += "Mass";
  saveName += TString::Format("_pt%.1f-%.1f", lowpt, highpt);
  saveName += ".pdf";
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);

  TH1D* mass = (TH1D*)thn->Projection(massAxis);
  if (normalise) { mass->Scale(1./mass->Integral(), "width"); }
  // mass->Rebin(4);
  setStyle(mass, 0);
  mass->SetName(saveName.c_str());

  if (isHistEmptyInRange(mass, 1, mass->GetNbinsX())) {
    cout << "Mass: empty histogram" << endl;
    return;
  }

  string xTitle = "M(" + formatHadronDaughters(hadron) + ") (GeV/#it{c}^{2})";
  string yTitle = "counts";
  if (normalise) yTitle = "#frac{1}{#it{N}_{V0}} #frac{d#it{N}}{d#it{M}}";
  string ptText = TString::Format("#it{p}_{T, %s} = %.1f - %.1f GeV/#it{c}", formatHadronName(hadron).c_str(), lowpt, highpt).Data();
  double xMinFrame = mass->GetXaxis()->GetXmin(), xMaxFrame = mass->GetXaxis()->GetXmax();
  double yMinFrame = 0., yMaxFrame = 1.1 * getHistScale(mass, true);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle((dataSet + ", " + ptText).c_str());

  double v0mass = ("K0S" == hadron) ? MassK0S : MassLambda0;
  TLine* line = new TLine(v0mass, yMinFrame, v0mass, yMaxFrame);
  line->SetLineStyle(9);
  line->SetLineColor(GetColor(1));
  line->SetLineWidth(3);

  canvas->cd();
  frame->Draw();
  line->Draw("same");
  mass->Draw("same");
  canvas->SaveAs(canvas->GetName());
}

void plotTrain(string inName, string dataSet, double partv0min, double partv0max, int setting)
{
  switch (setting) {
    case 0:
      ptResolution(inName, dataSet, partv0min, partv0max);
      break;
    case 1:
      dauPtResolution(inName, dataSet, partv0min, partv0max, true);
      dauPtResolution(inName, dataSet, partv0min, partv0max, false);
      break;
    case 2:
      {
        vector<string> inputStrings = {inName, dataSet, "K0S"};
        mass(inputStrings, partv0min, partv0max, false);
        inputStrings[2] = "Lambda";
        mass(inputStrings, partv0min, partv0max, false);
        inputStrings[2] = "antiLambda";
        mass(inputStrings, partv0min, partv0max, false);
      }
      break;
    default:
      cout << "Error: invalid setting" << endl;
      break;
  }
}
void plot271952(double partv0min, double partv0max, int setting)
{
  string inName = "~/cernbox/TrainOutput/271952/AnalysisResults.root";
  string dataSet = "LHC24b1b";
  plotTrain(inName, dataSet, partv0min, partv0max, setting);
}
void plot271952(int setting)
{
  gROOT->SetBatch();
  vector<double> pt = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};
  for (int i = 0; i < pt.size() - 1; i++) {
    plot271952(pt[i], pt[i+1], setting);
  }
}
