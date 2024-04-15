
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

double getNjets(TFile* inFile, double jetptmin, double jetptmax)
{
  string histName = "jet-fragmentation/matching/jets/matchPartJetPtEtaPhi";
  TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
  std::array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
  return jets->Integral(jetptbins[0], jetptbins[1], 0, jets->GetNbinsY() + 1, 0, jets->GetNbinsZ() + 1);
}

void matchedJetPt(string inName = "AnalysisResults.root", bool detector = false)
{
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 200., yMinFrame = 1e-6, yMaxFrame = 2.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.1, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet}^{part.} (GeV/#it{c})";
  if (detector) { xTitle = "#it{p}_{T, jet}^{det.} (GeV/#it{c})"; }
  yTitle = "normalised count";
  dataSet = "LHC23k4b_pass1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation_id10235/matching/jets/matchDetJetPtPartJetPt";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());
  TH1D* matchedjetpt = (TH1D*)th2->ProjectionY("matchedjetpt", 0, th2->GetNbinsX()+1);
  if (detector) { matchedjetpt = (TH1D*)th2->ProjectionX("matchedjetpt", 0, th2->GetNbinsY()+1);}
  matchedjetpt->Scale(1./ matchedjetpt->Integral());
  setStyle(matchedjetpt, 0);
  histVector.push_back(matchedjetpt);

  latexText = TString::Format("%s, geometrically matched jets", dataSet.c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedPartJetPt";
  if (detector) { saveName = "matchedDetJetPt"; }
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
void matchedJetPt2D(string inName = "AnalysisResults.root")
{
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool setLogZ = true;
  double xMinFrame = 0., xMaxFrame = 200., yMinFrame = 0., yMaxFrame = 200.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.1, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, V0}^{det.} (GeV/#it{c})";
  yTitle = "#it{p}_{T, V0}^{part.} (GeV/#it{c})";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH2D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogZ) { canvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/jets/matchDetJetPtPartJetPt";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH2D* matchedjetpt = (TH2D*)inFile->Get(histName.c_str());
  normaliseHistRowByRow(matchedjetpt);
  matchedjetpt->SetMinimum(1e-5);
  matchedjetpt->SetMaximum(1.);
  matchedjetpt->SetName("matchedjetpt");
  histVector.push_back(matchedjetpt);

  latexText = TString::Format("%s, geometrically matched jets", dataSet.c_str()).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "matchedJetPt";
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}
void jetEnergyScale(string inName = "AnalysisResults.root", double partjetptmin = 10., double partjetptmax = 200.)
{
  gStyle->SetNdivisions(505, "xy");

  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;

  bool SetLogy = true;
  double xMinFrame = -1., xMaxFrame = 10., yMinFrame = 1e-5, yMaxFrame = 2.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.4, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "(#it{p}_{T, jet}^{det.} - #it{p}_{T, jet}^{part.})/#it{p}_{T, jet}^{part.}";
  // yTitle = "#it{p}_{T, jet}^{part.} (GeV/#it{c})";
  yTitle = "normalised count";
  dataSet = "LHC23k4b_pass1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (SetLogy) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  // histName = "jet-fragmentation/matching/jets/matchPartJetptRelDiffPt";
  histName = "jet-fragmentation_id10235/matching/jets/matchPartJetPtRelDiffPt";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(th2->GetXaxis(), partjetptmin, partjetptmax);
  TH1D* jes = th2->ProjectionY("jetEnergyScale", partjetptbins[0], partjetptbins[1]);
  jes->Scale(1./ jes->Integral());
  setStyle(jes, 0);
  histVector.push_back(jes);

  // Get RMS = jet resolution for given pt range
  double rms = 0.;
  for (int i = 1; i < jes->GetNbinsX(); i++) {
    double binContent = jes->GetBinContent(i);
    double binCenter = jes->GetBinCenter(i);
    // Hist is self-normalised, so no need to divide by integral
    rms += binContent * binCenter * binCenter;
  }
  rms = TMath::Sqrt(rms);
  jes->Print("all");

  double lowjetpt = th2->GetXaxis()->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = th2->GetXaxis()->GetBinUpEdge(partjetptbins[1]);
  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ RMS: %.2f } }", dataSet.c_str(), lowjetpt, highjetpt, rms).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "jetEnergyScale";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), partjetptmin, partjetptmax);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "colz");
}

