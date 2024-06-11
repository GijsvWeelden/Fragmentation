
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

// This script compares the jet spectra of Ch. jets and Ch+V0 jets

double getNevts(TFile* f, bool withV0s = false)
{
  TH1D* hNEvents = (TH1D*)f->Get("hNEvts");
  return hNEvents->GetBinContent(1 + (int)withV0s);
}

void ptSpectraAndRatio()
{
  double minEta = -0.5, maxEta = 0.5;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 5., xMaxFrame = 50., yMinFrame = 5e-8, yMaxFrame = 5e-3;
  double xMinLegend = 0.4, xMaxLegend = 0.8, yMinLegend = 0.5, yMaxLegend = 0.65;
  double xLatex = 0.3, yLatex = 0.81;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet}";
  yTitle = "normalised count";
  yTitle = "#it{N}_{jets} / #it{N}_{evts}";
  drawoption = "same";
  saveName = "Pythia-jetptComparison";
  vector<string> histNames = { "hChJet", "hChJetCorrected", "hV0Jet"};
  vector<string> legendEntries = { "Ch. jet", "Ch. jet (corrected)", "Ch+V0 jet"};

  vector<TH1D*> histVector; vector<TH1D*> ratioVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  TPad* padSpectrum = new TPad("padSpectrum", "padSpectrum", 0, 0.3, 1, 1);
  TPad* padRatio = new TPad("padRatio", "padRatio", 0, 0, 1, 0.3);
  padRatio->SetTopMargin(0.);
  if (setLogY) { padSpectrum->SetLogy(); }
  // if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  // latexText = "PYTHIA simulation";
  latexText = TString::Format("#splitline{ %s }{#splitline{ %s }{ %s }}",
                                "Pythia simulation pp",
                                "#sqrt{#it{s}} = 13.6 TeV",
                                TString::Format("Anti-k_{T} jets, #it{R} = 0.4, |#it{#eta}_{jet}| < %.1f", maxEta).Data()
                              );
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  string inName = "./Pythia-chjets-v0jets.root";
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  for (int i = 0; i < histNames.size(); i++) {
    string histName = histNames[i];

    TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
    array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), minEta, maxEta);
    TH1D* th1 = jets->ProjectionX(TString::Format("jetpt%s", histName.c_str()).Data(), etabins[0], etabins[1]);
    setStyle(th1, i);
    th1->Rebin(5);

    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);
  }
  // Make ratios
  for (int i = 0; i < histVector.size(); i++) {
    TH1D* ratio = (TH1D*)histVector[i]->Clone();
    ratio->Divide(histVector[0]);
    ratioVector.push_back(ratio);
  }

  padSpectrum->Draw();
  padSpectrum->cd();
  // TH1F* fSpectrum = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, "", yTitle);
  gPad->SetLeftMargin(0.22);
  gPad->SetBottomMargin(0.1);
  gPad->SetRightMargin(0.15);//0.05);
  gPad->SetTopMargin(0.);//0.05);
  TH1F *fSpectrum = gPad->DrawFrame(xMinFrame, yMinFrame, xMaxFrame, yMaxFrame);
  fSpectrum->SetXTitle("");
  fSpectrum->SetYTitle(yTitle.c_str());
  fSpectrum->GetXaxis()->SetLabelSize(0.05);
  fSpectrum->GetYaxis()->SetLabelSize(0.05);
  fSpectrum->GetXaxis()->SetTitleSize(0.06);
  fSpectrum->GetYaxis()->SetTitleSize(0.06);
  fSpectrum->GetXaxis()->SetTitleOffset(1.0);
  fSpectrum->GetYaxis()->SetTitleOffset(1.6);//1.3);
  fSpectrum->GetXaxis()->CenterTitle(true);
  fSpectrum->GetYaxis()->CenterTitle(true);
  gPad->SetTicks(1,1);
  fSpectrum->Draw();
  for (int i = 0; i < histVector.size(); i++) {
    TH1D* hist = histVector[i];
    hist->Scale(1./getNevts(inFile));
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");
  latex->Draw("same");

  canvas->cd();
  padRatio->Draw();
  padRatio->cd();
  // TH1F* fRatio = DrawFrame(xMinFrame, xMaxFrame, 0.8, 3, xTitle, "ratio");
  gPad->SetLeftMargin(0.22);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.15);//0.05);
  gPad->SetTopMargin(0.);//0.05);
  TH1F *fRatio = gPad->DrawFrame(xMinFrame, 0.5, xMaxFrame, 3.5);
  fRatio->SetXTitle(xTitle.c_str());
  fRatio->SetYTitle("Ratio");
  fRatio->GetXaxis()->SetLabelSize(0.05);
  fRatio->GetYaxis()->SetLabelSize(0.05);
  fRatio->GetXaxis()->SetTitleSize(0.06);
  fRatio->GetYaxis()->SetTitleSize(0.06);
  fRatio->GetXaxis()->SetTitleOffset(1.0);
  fRatio->GetYaxis()->SetTitleOffset(1.6);//1.3);
  fRatio->GetXaxis()->CenterTitle(true);
  fRatio->GetYaxis()->CenterTitle(true);
  gPad->SetTicks(1,1);

  fRatio->Draw();
  for (auto hist : ratioVector) {
    hist->Draw("same");
  }

  canvas->cd();
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

void comparePt(bool doRatio = false)
{
  double minEta = -0.5, maxEta = 0.5;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 150., yMinFrame = 1e-10, yMaxFrame = 5e-3;
  double xMinLegend = 0.4, xMaxLegend = 0.8, yMinLegend = 0.5, yMaxLegend = 0.65;
  double xLatex = 0.3, yLatex = 0.81;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet}";
  yTitle = "normalised count";
  yTitle = "#it{N}_{jets} / #it{N}_{evts}";
  drawoption = "same";
  saveName = "Pythia-jetptComparison";
  vector<string> histNames = { "hChJet", "hChJetCorrected", "hV0Jet"};
  vector<string> legendEntries = { "Ch. jet", "Ch. jet (corrected)", "Ch+V0 jet"};

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  // latexText = "PYTHIA simulation";
  latexText = TString::Format("#splitline{ %s }{#splitline{ %s }{ %s }}",
                                "Pythia simulation pp",
                                "#sqrt{#it{s}} = 13.6 TeV",
                                TString::Format("Anti-k_{T} jets, #it{R} = 0.4, |#it{#eta}_{jet}| < %.1f", maxEta).Data()
                              );
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  string inName = "./Pythia-chjets-v0jets.root";
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  for (int i = 0; i < histNames.size(); i++) {
    string histName = histNames[i];

    TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
    array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), minEta, maxEta);
    TH1D* th1 = jets->ProjectionX(TString::Format("jetpt%s", histName.c_str()).Data(), etabins[0], etabins[1]);
    setStyle(th1, i);
    th1->Rebin(5);

    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.5; yMaxFrame = 1.1;
    // drawoption = "same p";
    saveName = "Pythia-jetptRatio";
    yTitle = TString::Format("ratio w.r.t. %s", legendEntries[0].c_str()).Data();
    for (int i = 1; i < histVector.size(); i++) {
      histVector[i]->Divide(histVector[0]);
    }
    histVector[0]->Divide(histVector[0]);
  }

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->Draw();
  for (int i = 0; i < histVector.size(); i++) {
    TH1D* hist = histVector[i];
    if (!doRatio) {
      hist->Scale(1./getNevts(inFile));
    }
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");
  latex->Draw("same");

  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

void compareZ(double jetptmin = 10., double jetptmax = 1e3, bool doRatio = false)
{
  double minEta = -0.5, maxEta = 0.5;
  cout << "jetptmin: " << jetptmin << ", jetptmax: " << jetptmax << endl;

  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int jetetaAxis = 1;
  const int jetphiAxis = 2;
  const int zAxis      = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double lowjetpt, highjetpt;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 1e-7, yMaxFrame = 1e-3;
  double xMinLegend = 0.25, xMaxLegend = 0.55, yMinLegend = 0.2, yMaxLegend = 0.35;
  // double xMinLegend = 0.25, xMaxLegend = 0.8, yMinLegend = 0.175, yMaxLegend = 0.325;
  double xLatex = 0.25, yLatex = 0.82;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 10;
  xTitle = "#it{z}_{V0}";
  yTitle = "#frac{1}{#it{N}_{evts}} #frac{d #it{N}}{d #it{z}_{V0}}";
  drawoption = "same";
  saveName = "Pythia-zComparison";

  vector<string> histNames = { "hzCh", "hzChCorrected", "hzV0"};
  vector<string> legendEntries = { "Ch. jet", "Ch. jet (corrected)", "Ch+V0 jet"};

  vector<double> integrals; vector<double> nevents;
  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  string inName = "./Pythia-chjets-v0jets.root";
  TFile *inFile = TFile::Open(inName.c_str());
  for (int i = 0; i < histNames.size(); i++) {
    string histName = histNames[i];
    THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

    std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
    thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
    lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
    highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
    array<int, 2> etabins = getProjectionBins(thn->GetAxis(jetetaAxis), minEta, maxEta);
    thn->GetAxis(jetetaAxis)->SetRange(etabins[0], etabins[1]);

    TH1D* th1 = (TH1D*)thn->Projection(zAxis);
    th1->SetName(TString::Format("z_%s", histName.c_str()).Data());
    th1->Rebin(rebinNumber);
    setStyle(th1, i);
    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.; yMaxFrame = 1.1;
    // drawoption = "same hist p";
    saveName = "Pythia-zRatio";
    yTitle = TString::Format("ratio w.r.t. %s", legendEntries[0].c_str()).Data();
    for (int i = 1; i < histVector.size(); i++) {
      histVector[i]->Divide(histVector[0]);
    }
    histVector[0]->Divide(histVector[0]);
  }

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->Draw();
  for (int i = 0; i < histVector.size(); i++) {
    TH1D* hist = histVector[i];
    if (!doRatio) {
      hist->Scale(1./getNevts(inFile), "width");
    }
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");

  // latexText = TString::Format("PYTHIA simulation, #it{p}_{T, jet} = %.0f - %.0f GeV/#it{c}", lowjetpt, highjetpt).Data();
  latexText = TString::Format("#splitline{ %s }{#splitline{ %s }{ %s }}",
                                "Pythia simulation pp, #sqrt{#it{s}} = 13.6 TeV",
                                "Anti-k_{T} jets, #it{R} = 0.4, |#it{#eta}_{jet}| < 0.5",
                                TString::Format("%.0f < #it{p}_{T, ch.+V0 jet} < %.0f GeV/#it{c}", lowjetpt, highjetpt).Data()
                              );
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex->Draw("same");

  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

