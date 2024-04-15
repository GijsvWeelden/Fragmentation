
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

void ptResolution(string inName = "", double partv0min = -1., double partv0max = 1e6)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }

  string hadron = "V0";
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = false;
  double xMinFrame = -1., xMaxFrame = 10., yMinFrame = 0., yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.3, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("(#it{p}_{T, %s}^{det.} - #it{p}_{T, %s}^{part.})/#it{p}_{T, %s}^{part.}", hadron.c_str(), hadron.c_str(), hadron.c_str()).Data();
  yTitle = "normalised count";
  dataSet = "LHC23k4b_pass1_small";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "V0PartPtRatioPtRelDiffPt";
  histName = TString::Format("jet-fragmentation/matching/V0/%s", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());

  std::array<int, 2> partv0bins = getProjectionBins(th3->GetXaxis(), partv0min, partv0max);
  th3->GetXaxis()->SetRange(partv0bins[0], partv0bins[1]);
  double lowv0 = th3->GetXaxis()->GetBinLowEdge(partv0bins[0]);
  double highv0 = th3->GetXaxis()->GetBinUpEdge(partv0bins[1]);
  if (lowv0 < 0.) { lowv0 = 0.; } // Avoids ugly pt>-0 in latextext

  TH1D* v0res = (TH1D*)th3->ProjectionZ("v0resolution", partv0bins[0], partv0bins[1], 0, th3->GetNbinsY() + 1);
  v0res->Scale(1./v0res->Integral());
  setStyle(v0res, 0);
  histVector.push_back(v0res);

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

  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, %s}^{part.} = %.0f - %.0f GeV/#it{c} }{ RMS: %.2f } }", dataSet.c_str(), formatHadronName(hadron).c_str(), lowv0, highv0, rms).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%sPtResolution", hadron.c_str()).Data();
  saveName = TString::Format("%s_partv0pt%.0f-%.0f", saveName.c_str(), lowv0, highv0);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
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
  dataSet = "LHC23k4b_pass1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = "jet-fragmentation/matching/V0/V0PartPtDetPt";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
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
