
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

#include "./histUtils.C"

////////////////////////////////////////////////////////////////////////
//
// This macro checks that the MC event weights are correctly implemented in O2Physics
//
////////////////////////////////////////////////////////////////////////

void normaliseHistRowByRow(TH2F* hist);
void normaliseHistColByCol(TH2F* hist);

template <typename T>
void setHistStyle(T hist, int styleNumber);

void plotOneHist(TCanvas* canvas, TH1F* frame, TH2F* hist, TLegend* legend, string saveName, string setDrawOption, string latexText);
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<TH1F*> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);

void plotCheckEventWeights(string input = "AnalysisResults.root", int setting = 1)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName = TString::Format("./%s", input.c_str()).Data();
  string saveDir = TString::Format("../Plots/%s", input.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* qaDir = (TDirectory*)inFile->Get("jet-finder-qa");
  TDirectory* mcpDir = (TDirectory*)inFile->Get("jet-fragmentation/particle-level/jets");
  TDirectory* mcdDir = (TDirectory*)inFile->Get("jet-fragmentation/detector-level/jets");

  double textSize = 0.03;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.75, yMaxLegend = 0.85;

  // Particle level
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetLogy();
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, "", textSize);
  std::vector<TH1F*> histVector;

  if (setting < 0) // Event weights
  {
    TH1F* eventWeights = (TH1F*)qaDir->Get("h_collision_eventweight_part");
    setHistStyle(eventWeights, 0);
    eventWeights->Scale(1./eventWeights->Integral());
    histVector.push_back(eventWeights);
    TH1F* frame = DrawFrame(1e-9, 1e1, 1e-4, 1, "weight", "normalised count");
    frame->SetTitle("Event weights");
    canvas->SetLogx();
    plotNHists(canvas, frame, histVector, legend, "event-weights.pdf", "", "");
  }
  else if (setting < 3) // Particle level
  {
    TH3F* mcp3D = (TH3F*)mcpDir->Get("partJetPtEtaPhi");
    TH1F* mcp = (TH1F*)mcp3D->ProjectionX("mcp", 1, mcp3D->GetNbinsY()+1, 1, mcp3D->GetNbinsZ()+1);
    setHistStyle(mcp, 0);
    mcp->Scale(1./mcp->Integral()); // Self-normalise
    TH1F* mcpWeighted = (TH1F*)qaDir->Get("h_jet_pt_part");
    setHistStyle(mcpWeighted, 1);
    mcpWeighted->Rebin(5);
    mcpWeighted->Scale(1./mcpWeighted->Integral()); // Self-normalise
    TH1F* mcpRatio = (TH1F*)mcpWeighted->Clone("mcpRatio");
    mcpRatio->Divide(mcp);

    if (setting == 1) {
      histVector.push_back(mcp);
      legend->AddEntry(mcp, "mcp unweighted");
      histVector.push_back(mcpWeighted);
      legend->AddEntry(mcpWeighted, "mcp weighted");
      TH1F* frame = DrawFrame(0, 200, 1e-8, 2, "#it{p}_{T}", "normalised count");
      frame->SetTitle("mcp");
      plotNHists(canvas, frame, histVector, legend, "mcp_comparison.pdf", "", "");
    }
    else {
      histVector.push_back(mcpRatio);
      TH1F* frame = DrawFrame(0, 200, 1e-5, 1e1, "#it{p}_{T}", "weighted/unweighted");
      frame->SetTitle("mcp");
      plotNHists(canvas, frame, histVector, legend, "mcp_ratio.pdf", "", "");
    }
  }    // Particle level
  else // Detector level
  {
    TH3F* mcd3D = (TH3F*)mcdDir->Get("detJetPtEtaPhi");
    TH1F* mcd = (TH1F*)mcd3D->ProjectionX("mcd", 1, mcd3D->GetNbinsY()+1, 1, mcd3D->GetNbinsZ()+1);
    setHistStyle(mcd, 0);
    mcd->Scale(1./mcd->Integral()); // Self-normalise
    TH1F* mcdWeighted = (TH1F*)qaDir->Get("h_jet_pt");
    setHistStyle(mcdWeighted, 1);
    mcdWeighted->Rebin(5);
    mcdWeighted->Scale(1./mcdWeighted->Integral()); // Self-normalise
    TH1F* mcdRatio = (TH1F*)mcdWeighted->Clone("mcdRatio");
    mcdRatio->Divide(mcd);

    if (setting == 3) {
      histVector.push_back(mcd);
      legend->AddEntry(mcd, "mcd unweighted");
      histVector.push_back(mcdWeighted);
      legend->AddEntry(mcdWeighted, "mcd weighted");
      TH1F* frame = DrawFrame(0, 200, 1e-8, 2, "#it{p}_{T}", "normalised count");
      frame->SetTitle("mcd");
      plotNHists(canvas, frame, histVector, legend, "mcd_comparison.pdf", "", "");
    }
    else {
      histVector.push_back(mcdRatio);
      TH1F* frame = DrawFrame(0, 200, 1e-5, 1e1, "#it{p}_{T}", "weighted/unweighted");
      frame->SetTitle("mcd");
      plotNHists(canvas, frame, histVector, legend, "mcd_ratio.pdf", "", "");
    }
  } // Detector level

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

template <typename T>
void setHistStyle(T hist, int styleNumber)
{
  hist->SetLineWidth(3);
  hist->SetLineColor(GetColor(styleNumber));
  hist->SetMarkerStyle(GetMarker(styleNumber));
  hist->SetMarkerColor(GetColor(styleNumber));
}

// Normalise 2D histogram row-by-row
void normaliseHistRowByRow(TH2F* hist)
{
  int firstColBin = 1, lastColBin = hist->GetNbinsX();
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY();
  for (int iRow = 1; iRow <= lastRowBin; iRow++) {
    double integral = hist->Integral(firstColBin, lastColBin, iRow, iRow);
    if (integral < 1) { continue; }
    for (int iCol = 1; iCol <= lastColBin; iCol++) {
      double binContent = hist->GetBinContent(iCol, iRow);
      binContent /= integral;
      hist->SetBinContent(iCol, iRow, binContent);
    }
  }
}
// Normalise 2D histogram col-by-col
void normaliseHistColByCol(TH2F* hist)
{
  int firstColBin = 1, lastColBin = hist->GetNbinsX();
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY();
  for (int iCol = 1; iCol <= lastColBin; iCol++) {
    double integral = hist->Integral(iCol, iCol, firstRowBin, lastRowBin);
    if (integral < 1) { continue; }
    for (int iRow = 1; iRow <= lastRowBin; iRow++) {
      double binContent = hist->GetBinContent(iCol, iRow);
      binContent /= integral;
      hist->SetBinContent(iCol, iRow, binContent);
    }
  }
}

void plotOneHist(TCanvas* canvas, TH1F* frame, TH2F* hist, TLegend* legend, string saveName, string setDrawOption, string latexText)
{
  canvas->cd();
  frame->Draw();
  string drawOption = "same colz";
  if (setDrawOption != "") {
    drawOption = TString::Format("%s %s", drawOption.c_str(), setDrawOption.c_str()).Data();
  }
  hist->Draw(drawOption.c_str());
  if (legend) { legend->Draw("same"); }
  if (latexText != "") { DrawLatex(0.3, 0.8, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<TH1F*> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText)
{
  canvas->cd();
  frame->Draw();
  string drawOption = "same";
  if (setDrawOption != "") {
    drawOption = TString::Format("%s %s", drawOption.c_str(), setDrawOption.c_str()).Data();
  }
  for (unsigned int i = 0; i < histVector.size(); i++) {
    TH1F* hist = histVector[i];
    hist->Draw(drawOption.c_str());
  }

  // auto fa = new TF1("fa","3 * TMath::Power(5, 2)/(TMath::Power(x,2))", 5, 200);
  // fa->SetLineWidth(5);
  // fa->SetLineColor(GetColor(0));
  // legend->AddEntry(fa, "1/x^{2}");
  // fa->Draw("same");

  // auto fb = new TF1("fb","3 * TMath::Power(5, 3)/(TMath::Power(x,3))", 5, 200);
  // fb->SetLineWidth(5);
  // fb->SetLineColor(GetColor(histVector.size() + 2));
  // legend->AddEntry(fb, "1/x^{3}");
  // fb->Draw("same");

  // auto fc = new TF1("fc","3 * TMath::Power(5, 4)/(TMath::Power(x,4))", 5, 200);
  // fc->SetLineWidth(5);
  // fc->SetLineColor(GetColor(histVector.size() + 1));
  // legend->AddEntry(fc, "1/x^{4}");
  // fc->Draw("same");

  if (legend) { legend->Draw("same"); }
  if (latexText != "") { DrawLatex(0.3, 0.8, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
