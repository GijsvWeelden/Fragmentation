
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
  histVector.push_back(v0res);

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
  saveName += TString::Format("_partv0pt%.1f-%.1f", saveName.c_str(), lowv0, highv0);
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

  histName  = "jet-fragmentation/matching/V0/";
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
