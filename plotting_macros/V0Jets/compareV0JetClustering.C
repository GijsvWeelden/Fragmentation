
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

#include "../histUtils.C"

// This script compares jets when we assign Lambdas the K0S mass or the Lambda mass:
// * Jet pt spectrum
// * K0S z spectrum
// * Lambda z spectrum

double getNevts(TFile* f, bool withV0s = false)
{
  TH1D* hNEvents = (TH1D*)f->Get("hNEvts");
  return hNEvents->GetBinContent(1 + (int)withV0s);
}

double getNjets(TH1* h, int minBin, int maxBin)
{
  return h->Integral(minBin, maxBin);
}

void comparePt(vector<string> inNames, vector<string> histNames, vector<string> legendEntries, TLatex* additionalLatex = nullptr, bool doRatio = false)
{
  double minEta = -0.5, maxEta = 0.5;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 150., yMinFrame = 1e-1, yMaxFrame = 2.;
  double xMinLegend = 0.25, xMaxLegend = 0.65, yMinLegend = 0.25, yMaxLegend = 0.4;
  double xLatex = 0.25, yLatex = 0.81;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet}";
  yTitle = "normalised count";
  yTitle = "#it{N}_{jets}";
  drawoption = "same";
  saveName = "Pythia-V0JetClustering-jetpt";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  latexText = TString::Format("#splitline{ %s }{#splitline{ %s }{#splitline{ %s }{ %s }}}",
                                "Pythia simulation pp",
                                "#sqrt{#it{s}} = 13.6 TeV",
                                "Anti-k_{T} jets, #it{R} = 0.4",
                                TString::Format("|#it{#eta}_{jet}| < %.1f", maxEta).Data()
                              );
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  double scale = 1.;
  for (int i = 0; i < histNames.size(); i++) {
    string inName = inNames[i];
    TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
    string histName = histNames[i];
    TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
    array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), minEta, maxEta);
    TH1D* th1 = jets->ProjectionX("jetpt", etabins[0], etabins[1]);
    th1->SetName(TString::Format("jetpt_%s_%s", histName.c_str(), inName.c_str()).Data());
    setStyle(th1, i);
    th1->Rebin(rebinNumber);
    th1->Sumw2();
    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);
    if (th1->GetBinContent(th1->GetMaximumBin()) > scale) {
      scale = th1->GetBinContent(th1->GetMaximumBin());
    }
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.9; yMaxFrame = 1.1;
    scale = 1.;
    drawoption = "same hist";
    yTitle = TString::Format("ratio w.r.t. %s", legendEntries[0].c_str()).Data();
    for (int i = 1; i < histVector.size(); i++) {
      histVector[i]->Divide(histVector[0]);
    }
    histVector[0]->Divide(histVector[0]);
  }

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame * scale, xTitle, yTitle);
  frame->Draw();
  for (int i = 0; i < histVector.size(); i++) {
    TH1D* hist = histVector[i];
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");
  latex->Draw("same");
  if (additionalLatex) additionalLatex->Draw("same");

  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += ".pdf";
  canvas->SaveAs(saveName.c_str());
}

// /*
void compareZ(vector<string> inNames, vector<string> histNames, vector<string> njetsNames, vector<string> legendEntries, string hadron, double jetptmin = 10., double jetptmax = 1e3, TLatex* additionalLatex = nullptr, bool doRatio = false)
{
  double minEta = -0.5, maxEta = 0.5;
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("V0" == hadron) != 1) {
    cout << "Hadron " << hadron << " not recognised. Should be K0S, Lambda0 or V0" << endl;
    return;
  }

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
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 1e-3, yMaxFrame = 1e3;
  double xMinLegend = 0.25, xMaxLegend = 0.55, yMinLegend = 0.2, yMaxLegend = 0.35;
  double xLatex = 0.25, yLatex = 0.82;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 10;
  xTitle = "#it{z}_{" + formatHadronName(hadron) + "}";
  yTitle = TString::Format("#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}_{%s}}", formatHadronName(hadron).c_str()).Data();
  drawoption = "same";
  saveName = "Pythia-z";
  saveName += hadron;

  vector<double> njets = { 0., 0. };
  vector<double> integrals; vector<double> nevents;
  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  for (int i = 0; i < histNames.size(); i++) {
    string inName = inNames[i];
    string histName = histNames[i];
    TFile *inFile = TFile::Open(inName.c_str());
    THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

    std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
    thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
    lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
    highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
    array<int, 2> etabins = getProjectionBins(thn->GetAxis(jetetaAxis), minEta, maxEta);
    thn->GetAxis(jetetaAxis)->SetRange(etabins[0], etabins[1]);
    TH1D* th1 = (TH1D*)thn->Projection(zAxis);
    th1->SetName(TString::Format("z_%s_%s", histName.c_str(), inName.c_str()).Data());
    th1->Rebin(rebinNumber);
    th1->Sumw2();
    setStyle(th1, i);
    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);

    TH3D* jetPtEtaPhi = (TH3D*)inFile->Get(njetsNames[i].c_str());
    jetptbins = getProjectionBins(jetPtEtaPhi->GetXaxis(), jetptmin, jetptmax);
    etabins = getProjectionBins(jetPtEtaPhi->GetYaxis(), minEta, maxEta);
    TH1D* jetPt = jetPtEtaPhi->ProjectionX("jetpt", etabins[0], etabins[1]);
    jetPt->SetName(TString::Format("jetpt_%s", histName.c_str()).Data());
    njets[i] = getNjets(jetPt, jetptbins[0], jetptbins[1]);
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.9; yMaxFrame = 1.1;
    drawoption = "same hist";
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
      // hist->Scale(1./getNevts(inFile), "width");
      hist->Scale(1./njets[i], "width");
    }
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");

  latexText = TString::Format("#splitline{ %s }{#splitline{ %s }{#splitline{ %s }{ %s }}}",
                                "Pythia simulation pp",
                                "#sqrt{#it{s}} = 13.6 TeV",
                                TString::Format("Anti-k_{T} jets, #it{R} = 0.4, |#it{#eta}_{jet}| < %.1f", maxEta).Data(),
                                TString::Format("%.0f < #it{p}_{T, jet} < %.0f GeV/#it{c}", lowjetpt, highjetpt).Data()
                              );
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex->Draw("same");
  if (additionalLatex) additionalLatex->Draw("same");

  saveName += (doRatio ? "_Ratio" : "_Comparison");
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName += ".pdf";
  canvas->SaveAs(saveName.c_str());
}
// */

void LasK_Escheme(string hadron, bool doRatio = false, double jetptmin = 10., double jetptmax = 1e3)
{
  string inName = "../../inputfiles/pythia/V0Study/v0jetclustering-Escheme.root";
  vector<string> jHistNames = { "hV0Jet", "hK0Jet"};
  vector<string> zHistNames = { "hzV0_" + hadron, "hzK0_" + hadron};
  vector<string> legendEntries = { "#Lambda mass = #it{m}_{#Lambda}", "#Lambda mass = #it{m}_{K^{0}_{S}}"};

  TLatex* additionalLatex = CreateLatex(0.65, 0.81, "#it{E} scheme", 0.04);

  if ("" == hadron) {
    comparePt({inName, inName}, jHistNames, legendEntries, additionalLatex, doRatio);
  }
  else {
    compareZ({inName, inName}, zHistNames, jHistNames, legendEntries, hadron, jetptmin, jetptmax, additionalLatex, doRatio);
  }
}

void LasK_Ptscheme(string hadron, bool doRatio = false, double jetptmin = 10., double jetptmax = 1e3)
{
  string inName = "../../inputfiles/pythia/V0Study/v0jetclustering-ptscheme.root";
  vector<string> jHistNames = { "hV0Jet", "hK0Jet"};
  vector<string> zHistNames = { "hzV0_" + hadron, "hzK0_" + hadron};
  vector<string> legendEntries = { "#Lambda mass = #it{m}_{#Lambda}", "#Lambda mass = #it{m}_{K^{0}_{S}}"};

  TLatex* additionalLatex = CreateLatex(0.65, 0.81, "#it{p}_{T} scheme", 0.04);

  if ("" == hadron) {
    comparePt({inName, inName}, jHistNames, legendEntries, additionalLatex, doRatio);
  }
  else {
    compareZ({inName, inName}, zHistNames, jHistNames, legendEntries, hadron, jetptmin, jetptmax, additionalLatex, doRatio);
  }
}

void EvsPtScheme(string hadron, bool V0Jets, bool doRatio = false, double jetptmin = 10., double jetptmax = 1e3)
{
  string EName = "../../inputfiles/pythia/V0Study/v0jetclustering-Escheme.root";
  string ptName = "../../inputfiles/pythia/V0Study/v0jetclustering-ptscheme.root";

  string jHistName = "hK0Jet";
  string zHistName = "hzK0_" + hadron;
  string legendEntry = "#Lambda mass = #it{m}_{K^{0}_{S}}";
  if (V0Jets) {
    jHistName = "hV0Jet";
    zHistName = "hzV0_" + hadron;
    legendEntry = "#Lambda mass = #it{m}_{#Lambda}";
  }

  vector<string> jHistNames = {jHistName, jHistName};
  vector<string> zHistNames = {zHistName, zHistName};
  vector<string> legendEntries = {"#it{E} scheme", "#it{p}_{T} scheme"};

  TLatex* additionalLatex = CreateLatex(0.6, 0.6, legendEntry.c_str(), 0.04);

  if ("" == hadron) {
    comparePt({EName, ptName}, jHistNames, legendEntries, additionalLatex, doRatio);
  }
  else {
    compareZ({EName, ptName}, zHistNames, jHistNames, legendEntries, hadron, jetptmin, jetptmax, additionalLatex, doRatio);
  }
}
