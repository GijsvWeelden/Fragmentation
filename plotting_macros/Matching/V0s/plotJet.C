
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

void matchedJetPt(string inName, bool detector = false)
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
  xTitle = "#it{p}_{T, ch+V0 jet}^{part.} (GeV/#it{c})";
  if (detector) { xTitle = "#it{p}_{T, ch+V0 jet}^{det.} (GeV/#it{c})"; }
  yTitle = "normalised count";
  dataSet = "LHC24b1";

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
void matchedJetPt2D(string inName)
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
  dataSet = "LHC24b1";

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
void jetEnergyScale(vector<string> inputStrings, double partjetptmin, double partjetptmax)
{
  gStyle->SetNdivisions(505, "xy");

  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];
  histName = "jet-fragmentation/matching/jets/matchPartJetPtRelDiffPt";
  TFile *inFile = TFile::Open(inName.c_str());
  TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(th2->GetXaxis(), partjetptmin, partjetptmax);
  TH1D* jes = th2->ProjectionY(TString::Format("jetEnergyScale_%s_pt%.0f-%.0f").Data(), inName.c_str(), partjetptbins[0], partjetptbins[1]);
  jes->Scale(1./ jes->Integral());
  setStyle(jes, 0);

  // Get RMS = jet resolution for given pt range
  double rms = 0.;
  for (int i = 1; i < jes->GetNbinsX(); i++) {
    double binContent = jes->GetBinContent(i);
    double binCenter = jes->GetBinCenter(i);
    // Hist is self-normalised, so no need to divide by integral
    rms += binContent * binCenter * binCenter;
  }
  rms = TMath::Sqrt(rms);
  // jes->Print("all");

  double xLatex = 0.4, yLatex = 0.8;
  double lowjetpt = th2->GetXaxis()->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = th2->GetXaxis()->GetBinUpEdge(partjetptbins[1]);
  latexText = TString::Format("#splitline{ %s }{#splitline{ #it{p}_{T, ch+V0 jet}^{part.} = %.0f - %.0f GeV/#it{c} }{ RMS: %.2f } }", dataSet.c_str(), lowjetpt, highjetpt, rms).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = "jes";
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  // double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  // TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  int xCanvas = 900, yCanvas = 900;
  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), xCanvas, yCanvas);
  canvas->SetLogy();
  double xMinFrame = -1., xMaxFrame = 10., yMinFrame = 1e-5, yMaxFrame = 2.;
  xTitle = "(#it{p}_{T, ch+V0 jet}^{det.} - #it{p}_{T, ch+V0 jet}^{part.})/#it{p}_{T, ch+V0 jet}^{part.}";
  yTitle = "normalised count";
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);

  frame->Draw();
  jes->Draw("same");
  latex->Draw("same");
  canvas->SaveAs(canvas->GetName());
}

void jetEnergyScale(vector<string> inputStrings, vector<vector<double> > partjetptranges)
{
  string inName  = inputStrings[0];
  string dataSet = inputStrings[1];

  string saveName = "jes";
  TLegend* legend = CreateLegend(0.5, 0.9, 0.6, 0.8, "", 0.04);

  string histName = "jet-fragmentation/matching/jets/matchPartJetPtRelDiffPt";
  TFile *inFile = TFile::Open(inName.c_str());
  TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());
  vector<TH1D*> jesVector;
  for (int i = 0; i < partjetptranges.size(); i++) {
    array<int, 2> partjetptbins = getProjectionBins(th2->GetXaxis(), partjetptranges[i][0], partjetptranges[i][1]);
    double lowpt = th2->GetXaxis()->GetBinLowEdge(partjetptbins[0]);
    double highpt = th2->GetXaxis()->GetBinUpEdge(partjetptbins[1]);
    saveName += TString::Format("_ptjet%.0f-%.0f", lowpt, highpt);
    TH1D* jes = th2->ProjectionY(TString::Format("jetEnergyScale_jetpt%.0f-%.0f").Data(), lowpt, highpt);
    jes->Scale(1./ jes->Integral());

    double rms = 0;
    for (int i = 1; i < jes->GetNbinsX(); i++) {
      double binContent = jes->GetBinContent(i);
      double binCenter  = jes->GetBinCenter(i);
      rms += binContent * binCenter * binCenter; // Hist is self-normalised, so no need to divide by integral
    }
    rms = TMath::Sqrt(rms);

    setStyle(jes, i);
    jesVector.push_back(jes);
    legend->AddEntry(jes, TString::Format("%.0f < #it{p}_{T, ch+V0 jet}^{part.} < %.0f GeV/#it{c} (RMS: %.2f)", lowpt, highpt, rms).Data());
  }

  TCanvas* canvas = new TCanvas(saveName.c_str(), saveName.c_str(), 900, 900);
  canvas->SetLogy();

  double xMinFrame = -1., xMaxFrame = 10., yMinFrame = 1e-5, yMaxFrame = 2.;
  string xTitle = "(#it{p}_{T, ch+V0 jet}^{det.} - #it{p}_{T, ch+V0 jet}^{part.})/#it{p}_{T, ch+V0 jet}^{part.}";
  string yTitle = "normalised count";
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle(dataSet.c_str());

  frame->Draw();
  for (auto jes : jesVector) {
    jes->Draw("same");
  }
  legend->Draw("same");
  canvas->SaveAs(canvas->GetName());
}

