
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

#include "/Users/gijsvanweelden/Documents/Fragmentation/plotting_macros/histUtils.C"

void comparePt(string dataSet, vector<string> inNames, vector<string> legendEntries, bool doPart = false, bool doRatio = false)
{
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 200., yMinFrame = 1e-7, yMaxFrame = 1.;
  double xMinLegend = 0.3, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.4, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet}";
  yTitle = "normalised count";
  drawoption = "same";
  saveName = "ptComparison";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  latexText = TString::Format("%s", dataSet.c_str()).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  for (int i = 0; i < inNames.size(); i++) {
    TFile *inFile = TFile::Open(TString::Format("./%s", inNames[i].c_str()).Data());
    histName = "jet-finder-charged-qa/h_jet_pt";
    if (2 == i) { histName = "jet-finder-charged-v0-qa/h_jet_pt"; }
    if (doPart) { histName = TString::Format("%s_part", histName.c_str()).Data(); }

    TH1D* th1 = (TH1D*)inFile->Get(histName.c_str());
    th1->Print();
    setStyle(th1, i);
    legend->AddEntry(th1, legendEntries[i].c_str(), "l");
    histVector.push_back(th1);
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.5; yMaxFrame = 1.1;
    drawoption = "same hist p";
    saveName = "ptRatio";
    yTitle = TString::Format("ratio w.r.t. %s", legendEntries[0].c_str()).Data();
    for (int i = 1; i < histVector.size(); i++) {
      histVector[i]->Divide(histVector[0]);
    }
    histVector[0]->Divide(histVector[0]);
  }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->Draw();
  for (auto hist : histVector) {
    if (!doRatio) {
      for (auto hist : histVector) {
        hist->Scale(1./hist->Integral());
      }
    }
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");
  latex->Draw("same");
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void compareEta(string dataSet, vector<string> inNames, vector<string> legendEntries, double jetptmin = 0., double jetptmax = 1e6, bool doPart = false, bool doRatio = false)
{
  gStyle->SetNdivisions(505, "xy");
  string saveName, histDir, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = false;
  double xMinFrame = -1., xMaxFrame = +1., yMinFrame = 0., yMaxFrame = .05;
  double xMinLegend = 0.3, xMaxLegend = 0.9, yMinLegend = 0.58, yMaxLegend = 0.78;
  double xLatex = 0.25, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#eta_{jet}";
  yTitle = "normalised count";
  drawoption = "same";
  saveName = "etaComparison";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  for (int i = 0; i < inNames.size(); i++) {
    TFile *inFile = TFile::Open(TString::Format("./%s", inNames[i].c_str()).Data());
    histDir = "jet-finder-charged-qa";
    histName = "jet-finder-charged-qa/h3_jet_r_jet_pt_jet_eta";
    if (2 == i) { histDir = "jet-finder-charged-v0-qa"; }
    if (doPart) { histName = "h3_jet_r_part_jet_pt_part_jet_eta_part"; }

    TH3D* th3 = (TH3D*)inFile->Get(TString::Format("%s/%s", histDir.c_str(), histName.c_str()).Data());
    th3->Print();
    std::array<int, 2> jetptbins = getProjectionBins(th3->GetYaxis(), jetptmin, jetptmax);
    TH1D* th1 = (TH1D*)th3->ProjectionZ(TString::Format("jeteta_%.0f-%.0f_%s",
                                          jetptmin, jetptmax, inNames[i].c_str()).Data(),
                                        0, th3->GetNbinsX()+1, jetptbins[0], jetptbins[1]);
    th1->Print();
    setStyle(th1, i);
    legend->AddEntry(th1, legendEntries[i].c_str(), "l");
    histVector.push_back(th1);

    double lowjetpt = th3->GetYaxis()->GetBinLowEdge(jetptbins[0]);
    double highjetpt = th3->GetYaxis()->GetBinUpEdge(jetptbins[1]);
    latexText = TString::Format("%s, #it{p}_{T, jet}: %.0f - %.0f GeV/#it{c}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.; yMaxFrame = 1.1;
    drawoption = "same hist p";
    saveName = "etaRatio";
    yTitle = TString::Format("ratio w.r.t. %s", legendEntries[0].c_str()).Data();
    for (int i = 1; i < histVector.size(); i++) {
      histVector[i]->Divide(histVector[0]);
    }
    histVector[0]->Divide(histVector[0]);
  }
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->Draw();
  for (auto hist : histVector) {
    if (!doRatio) {
      for (auto hist : histVector) {
        hist->Scale(1./hist->Integral());
      }
    }
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");
  latex->Draw("same");
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
void comparePhi(string dataSet, vector<string> inNames, vector<string> legendEntries, double jetptmin = 0., double jetptmax = 1e6, bool doPart = false, bool doRatio = false)
{
  gStyle->SetNdivisions(505, "xy");
  string saveName, histDir, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = false;
  double xMinFrame = 0., xMaxFrame = 6.3, yMinFrame = 0., yMaxFrame = .05;
  double xMinLegend = 0.3, xMaxLegend = 0.9, yMinLegend = 0.58, yMaxLegend = 0.78;
  double xLatex = 0.25, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#varphi_{jet}";
  yTitle = "normalised count";
  drawoption = "same";
  saveName = "phiComparison";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  for (int i = 0; i < inNames.size(); i++) {
    TFile *inFile = TFile::Open(TString::Format("./%s", inNames[i].c_str()).Data());
    histDir = "jet-finder-charged-qa";
    histName = "h3_jet_r_jet_pt_jet_phi";
    if (2 == i) { histDir = "jet-finder-charged-v0-qa"; }
    if (doPart) { histName = "h3_jet_r_part_jet_pt_part_jet_phi_part"; }
    histDir = "jet-finder-charged-qa";

    TH3D* th3 = (TH3D*)inFile->Get(TString::Format("%s/%s", histDir.c_str(), histName.c_str()).Data());
    th3->Print();
    std::array<int, 2> jetptbins = getProjectionBins(th3->GetYaxis(), jetptmin, jetptmax);
    TH1D* th1 = (TH1D*)th3->ProjectionZ(TString::Format("jetphi_%.0f-%.0f_%s",
                                          jetptmin, jetptmax, inNames[i].c_str()).Data(),
                                        0, th3->GetNbinsX()+1, jetptbins[0], jetptbins[1]);
    th1->Print();
    setStyle(th1, i);
    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);

    double lowjetpt = th3->GetYaxis()->GetBinLowEdge(jetptbins[0]);
    double highjetpt = th3->GetYaxis()->GetBinUpEdge(jetptbins[1]);
    latexText = TString::Format("%s, #it{p}_{T, jet}: %.0f - %.0f GeV/#it{c}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.; yMaxFrame = 1.1;
    drawoption = "same hist p";
    saveName = "phiRatio";
    yTitle = TString::Format("ratio w.r.t. %s", legendEntries[0].c_str()).Data();
    for (int i = 1; i < histVector.size(); i++) {
      histVector[i]->Divide(histVector[0]);
    }
    histVector[0]->Divide(histVector[0]);
  }
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->Draw();
  for (auto hist : histVector) {
    if (!doRatio) {
      for (auto hist : histVector) {
        hist->Scale(1./hist->Integral());
      }
    }
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");
  latex->Draw("same");
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

void compareptdata(bool doRatio = false)
{ comparePt("Data: LHC23y", {"ChJets:ChJets.root", "ChJets:V0Jets.root", "V0Jets:V0Jets.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)", "Ch.+V0 jets"}, false, doRatio); }
void compareetadata(double jetptmin = 0., double jetptmax = 1e6, bool doRatio = false)
{ compareEta("Data: LHC23y", {"ChJets:ChJets.root", "ChJets:V0Jets.root", "V0Jets:V0Jets.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)", "Ch.+V0 jets"}, jetptmin, jetptmax, false, doRatio); }
void comparephidata(double jetptmin = 0., double jetptmax = 1e6, bool doRatio = false)
{ comparePhi("Data: LHC23y", {"ChJets:ChJets.root", "ChJets:V0Jets.root", "V0Jets:V0Jets.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)", "Ch.+V0 jets"}, jetptmin, jetptmax, false, doRatio); }

void compareptmcd(bool doRatio = false)
{ comparePt("MC: LHC23d4", {"26-04-ChJets:ChJets_MC.root", "26-04-ChJets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)"}, false, doRatio); }
// { comparePt("MC: LHC23d4", {"ChJets:ChJets_MC.root", "ChJets:V0Jets_MC.root", "V0Jets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)", "Ch.+V0 jets"}, false, doRatio); }
void compareetamcd(double jetptmin = 0., double jetptmax = 1e6, bool doRatio = false)
{ compareEta("MC: LHC23d4", {"26-04-ChJets:ChJets_MC.root", "26-04-ChJets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)"}, jetptmin, jetptmax, false, doRatio); }
// { compareEta("MC: LHC23d4", {"ChJets:ChJets_MC.root", "ChJets:V0Jets_MC.root", "V0Jets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)", "Ch.+V0 jets"}, jetptmin, jetptmax, false, doRatio); }
void comparephimcd(double jetptmin = 0., double jetptmax = 1e6, bool doRatio = false)
{ comparePhi("MC: LHC23d4", {"26-04-ChJets:ChJets_MC.root", "26-04-ChJets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)"}, jetptmin, jetptmax, false, doRatio); }
// { comparePhi("MC: LHC23d4", {"ChJets:ChJets_MC.root", "ChJets:V0Jets_MC.root", "V0Jets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)", "Ch.+V0 jets"}, jetptmin, jetptmax, false, doRatio); }

void compareptmcp(bool doRatio = false)
{ comparePt("MC: LHC23d4", {"26-04-ChJets:ChJets_MC.root", "26-04-ChJets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)"}, true, doRatio); }
// { comparePt("MC: LHC23d4", {"ChJets:ChJets_MC.root", "ChJets:V0Jets_MC.root", "V0Jets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)", "Ch.+V0 jets"}, true, doRatio); }
void compareetamcp(double jetptmin = 0., double jetptmax = 1e6, bool doRatio = false)
{ compareEta("MC: LHC23d4", {"26-04-ChJets:ChJets_MC.root", "26-04-ChJets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)"}, jetptmin, jetptmax, true, doRatio); }
// { compareEta("MC: LHC23d4", {"ChJets:ChJets_MC.root", "ChJets:V0Jets_MC.root", "V0Jets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)", "Ch.+V0 jets"}, jetptmin, jetptmax, true, doRatio); }
void comparephimcp(double jetptmin = 0., double jetptmax = 1e6, bool doRatio = false)
{ comparePhi("MC: LHC23d4", {"26-04-ChJets:ChJets_MC.root", "26-04-ChJets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)"}, jetptmin, jetptmax, true, doRatio); }
// { comparePhi("MC: LHC23d4", {"ChJets:ChJets_MC.root", "ChJets:V0Jets_MC.root", "V0Jets:V0Jets_MC.root"}, {"Ch. jets (master branch)", "Ch. jets (V0 jets branch)", "Ch.+V0 jets"}, jetptmin, jetptmax, true, doRatio); }
