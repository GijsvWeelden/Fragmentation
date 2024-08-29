
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

// This script compares ch jets with ch+V0 jets and ch+V0 jets where we have subtracted the V0

double getNevts(TFile* f, bool withV0s = false)
{
  TH1D* hNEvents = (TH1D*)f->Get("hNEvts");
  return hNEvents->GetBinContent(1 + (int)withV0s);
}

double getNjets(TH1* h, int minBin, int maxBin)
{
  return h->Integral(minBin, maxBin);
}

void comparePt(vector<vector<string> > inputStrings, TLatex* additionalLatex = nullptr, bool doRatio = false)
{
  double minEta = -0.5, maxEta = 0.5;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  // double xMinFrame = 0., xMaxFrame = 150., yMinFrame = 1e-1, yMaxFrame = 2.;
  double xMinFrame = 0., xMaxFrame = 100., yMinFrame = 1e-1, yMaxFrame = 2.;
  double xMinLegend = 0.35, xMaxLegend = 0.75, yMinLegend = 0.7, yMaxLegend = 0.85;
  double xLatex = 0.25, yLatex = 0.3;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet}";
  yTitle = "#it{N}_{jets}";
  drawoption = "same";

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
  // size - 1 because the last entry is the saveName
  for (int i = 0; i < inputStrings.size() - 1; i++) {
    string inName      = inputStrings[i][0];
    string histName    = inputStrings[i][1];
    string legendEntry = inputStrings[i][2];

    TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
    TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
    jets->Sumw2();
    array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), minEta, maxEta);

    TH1D* jetpt = jets->ProjectionX("jetpt", etabins[0], etabins[1]);
    jetpt->SetName(TString::Format("jetpt_%s_%s", histName.c_str(), inName.c_str()).Data());
    setStyle(jetpt, i);
    jetpt->Rebin(rebinNumber);
    if (jetpt->GetBinContent(jetpt->GetMaximumBin()) > scale) {
      scale = jetpt->GetBinContent(jetpt->GetMaximumBin());
    }

    legend->AddEntry(jetpt, legendEntry.c_str());
    histVector.push_back(jetpt);
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.9; yMaxFrame = 1.1;
    scale = 1.;
    drawoption = "same hist";
    yTitle = TString::Format("ratio w.r.t. %s", inputStrings[0][2].c_str()).Data();
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

  saveName = inputStrings[inputStrings.size() - 1][0];
  canvas->SaveAs(saveName.c_str());
}
void compareZ(vector<vector<string> > inputStrings, double jetptmin = 10., double jetptmax = 1e3, TLatex* additionalLatex = nullptr, bool doRatio = false)
{
  double minEta = -0.5, maxEta = 0.5;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 1e-3, xMaxFrame = 1. + 1e-3, yMinFrame = 1e-3, yMaxFrame = 1e3;
  double xMinLegend = 0.25, xMaxLegend = 0.75, yMinLegend = 0.7, yMaxLegend = 0.85;
  // double xMinLegend = 0.35, xMaxLegend = 0.75, yMinLegend = 0.7, yMaxLegend = 0.85;
  double xLatex = 0.25, yLatex = 0.3;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{z}_{V0}";
  yTitle = "#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}_{V0}}";
  drawoption = "same";

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

  double scale = 0.;
  double lowjetpt = -1., highjetpt = -1.;
  vector<double> njets;
  // size - 1 because the last entry is the saveName
  for (int i = 0; i < inputStrings.size() - 1; i++) {
    string inName      = inputStrings[i][0];
    string histName    = inputStrings[i][1];
    string legendEntry = inputStrings[i][2];
    string jetHistName = inputStrings[i][3];

    TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
    TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
    jets->Sumw2();
    array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
    array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), minEta, maxEta);

    TH1D* zv0 = jets->ProjectionZ("zv0", jetptbins[0], jetptbins[1], etabins[0], etabins[1]);
    zv0->SetName(TString::Format("zv0_%s_%s", histName.c_str(), inName.c_str()).Data());
    setStyle(zv0, i);
    zv0->Rebin(rebinNumber);
    lowjetpt = jetptbins[0]; highjetpt = jetptbins[1];

    legend->AddEntry(zv0, legendEntry.c_str());

    TH3D* jetPtEtaPhi = (TH3D*)inFile->Get(jetHistName.c_str());
    jetptbins = getProjectionBins(jetPtEtaPhi->GetXaxis(), jetptmin, jetptmax);
    etabins = getProjectionBins(jetPtEtaPhi->GetYaxis(), minEta, maxEta);
    TH1D* jetPt = jetPtEtaPhi->ProjectionX("jetpt", etabins[0], etabins[1]);
    jetPt->SetName(TString::Format("jetpt_%s", histName.c_str()).Data());
    njets.push_back(getNjets(jetPt, jetptbins[0], jetptbins[1]));

    if (zv0->GetBinContent(zv0->GetMaximumBin()) / njets[i] > scale) {
      scale = zv0->GetBinContent(zv0->GetMaximumBin()) / njets[i];
    }
    zv0->Scale(1./njets[i], "width");
    histVector.push_back(zv0);
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.9; yMaxFrame = 1.1;
    scale = 1.;
    drawoption = "same hist";
    yTitle = TString::Format("ratio w.r.t. %s", inputStrings[0][2].c_str()).Data();
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

  // saveName = "Pythia-WeightedJetFinder";
  // saveName += "-zv0";
  // saveName += (doRatio ? "Ratio" : "Comparison");
  // saveName += TString::Format("_jetpt%.0f-%.0f", jetptmin, jetptmax).Data();
  // saveName += ".pdf";
  saveName = inputStrings[inputStrings.size() - 1][0];
  canvas->SaveAs(saveName.c_str());
}

// -------------------------------------------------------------------------------------------------

void ptjet(bool doRatio = false, bool ptscheme = true)
{
  vector<vector<string> > inputStrings;
  string inName = "../../inputfiles/pythia/V0Study/weightedjetfinder.root";
  string saveName = "WeightedJetFinder";
  saveName += "-jetpt";
  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += "-all";
  saveName += ".pdf";

  vector<string> v0Jets = { "hV0Jet", "Ch+V0 jets" };
  vector<string> w0Jets = { "hW0Jet", "Ch+V0 jets (weighted)" };
  vector<string> eJets  = { "hEJet",  "Ch+V0 jets - #it{p}_{V0} (weighted)" };
  vector<string> ptJets = { "hPtJet", "Ch+V0 jets - #it{p}_{T, V0} (weighted)" };

  string latexText = ptscheme ? "#it{p}_{T} scheme" : "#it{E} scheme";
  TLatex* additionalLatex = CreateLatex(0.5, 0.6, latexText.c_str(), 0.04);

  // Spectra
  inputStrings.push_back({inName, v0Jets[0], v0Jets[1]});
  inputStrings.push_back({inName, w0Jets[0], w0Jets[1]});
  inputStrings.push_back({inName, eJets[0],  eJets[1]});
  inputStrings.push_back({inName, ptJets[0], ptJets[1]});

  inputStrings.push_back({saveName});
  comparePt(inputStrings, additionalLatex, doRatio);
}
void ptjetmatched(bool doRatio = false, bool ptscheme = true)
{
  vector<vector<string> > inputStrings;
  string inName = "../../inputfiles/pythia/V0Study/weightedjetfinder.root";
  string saveName = "WeightedJetFinder";
  saveName += "-jetpt";
  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += "-matched-all";
  saveName += ".pdf";

  vector<string> v0JetsWV0 = { "hV0JetMatched", "Ch+V0 jets" };
  vector<string> w0JetsWV0 = { "hW0JetMatched", "Ch+V0 jets (weighted)" };
  vector<string> eJets     = { "hEJetMatched",  "Ch+V0 jets - #it{p}_{V0} (weighted)" };
  vector<string> ptJets    = { "hPtJetMatched", "Ch+V0 jets - #it{p}_{T, V0} (weighted)" };

  string latexText = ptscheme ? "#it{p}_{T} scheme" : "#it{E} scheme";
  TLatex* additionalLatex = CreateLatex(0.65, 0.81, latexText.c_str(), 0.04);

  // Matched spectra
  inputStrings.push_back({inName, v0JetsWV0[0], v0JetsWV0[1]});
  inputStrings.push_back({inName, w0JetsWV0[0], w0JetsWV0[1]});
  inputStrings.push_back({inName, eJets[0],  eJets[1]});
  inputStrings.push_back({inName, ptJets[0], ptJets[1]});

  inputStrings.push_back({saveName});
  comparePt(inputStrings, additionalLatex, doRatio);
}
void ptjetweighted(bool doRatio = false, bool ptscheme = true)
{
  vector<vector<string> > inputStrings;
  string inName = "../../inputfiles/pythia/V0Study/weightedjetfinder.root";
  string saveName = "WeightedJetFinder";
  saveName += "-jetpt";
  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += "-weighted";
  saveName += ".pdf";

  vector<string> w0JetsWV0 = { "hW0Jet", "Ch+V0 jets (weighted)" };
  vector<string> eJets     = { "hEJet",  "Ch+V0 jets - #it{p}_{V0} (weighted)" };
  vector<string> ptJets    = { "hPtJet", "Ch+V0 jets - #it{p}_{T, V0} (weighted)" };

  string latexText = ptscheme ? "#it{p}_{T} scheme" : "#it{E} scheme";
  TLatex* additionalLatex = CreateLatex(0.65, 0.81, latexText.c_str(), 0.04);

  // Matched spectra
  inputStrings.push_back({inName, w0JetsWV0[0], w0JetsWV0[1]});
  inputStrings.push_back({inName, eJets[0],  eJets[1]});
  inputStrings.push_back({inName, ptJets[0], ptJets[1]});

  inputStrings.push_back({saveName});
  comparePt(inputStrings, additionalLatex, doRatio);
}
void ptjetmatchedweighted(bool doRatio = false, bool ptscheme = true)
{
  vector<vector<string> > inputStrings;
  string inName = "../../inputfiles/pythia/V0Study/weightedjetfinder.root";
  string saveName = "WeightedJetFinder";
  saveName += "-jetpt";
  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += "-matched-weighted";
  saveName += ".pdf";

  vector<string> w0JetsWV0 = { "hW0JetMatched", "Ch+V0 jets (weighted)" };
  vector<string> eJets     = { "hEJetMatched",  "Ch+V0 jets - #it{p}_{V0} (weighted)" };
  vector<string> ptJets    = { "hPtJetMatched", "Ch+V0 jets - #it{p}_{T, V0} (weighted)" };

  string latexText = ptscheme ? "#it{p}_{T} scheme" : "#it{E} scheme";
  TLatex* additionalLatex = CreateLatex(0.65, 0.81, latexText.c_str(), 0.04);

  // Matched spectra
  inputStrings.push_back({inName, w0JetsWV0[0], w0JetsWV0[1]});
  inputStrings.push_back({inName, eJets[0],  eJets[1]});
  inputStrings.push_back({inName, ptJets[0], ptJets[1]});

  inputStrings.push_back({saveName});
  comparePt(inputStrings, additionalLatex, doRatio);
}

void v0z(double jetptmin, double jetptmax, bool doRatio = false)
{
  vector<vector<string> > inputStrings;
  string inName = "../../inputfiles/pythia/V0Study/weightedjetfinder.root";
  bool ptscheme = true;
  string saveName = "WeightedJetFinder";
  saveName += "-zv0";
  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += TString::Format("_jetpt%.0f-%.0f", jetptmin, jetptmax).Data();
  saveName += "-all";
  saveName += ".pdf";

  vector<string> v0JetsZ = { "hV0JetZ", "Unweighted", "hV0Jet" };
  vector<string> w0JetsZ = { "hW0JetZ", "Weighted", "hW0Jet" };
  vector<string> eJetsZ  = { "hEJetZ",  "Ch+V0 jets - #it{p}_{V0} (weighted)", "hEJet" };
  vector<string> ptJetsZ = { "hPtJetZ", "Ch+V0 jets - #it{p}_{T, V0} (weighted)", "hPtJet" };

  string latexText = TString::Format("#it{p}_{T} scheme, %.0f #leq #it{p}_{T, jet} < %.0f GeV/#it{c}", jetptmin, jetptmax).Data();
  TLatex* additionalLatex = CreateLatex(0.25, 0.93, latexText.c_str(), 0.04);

  inputStrings.push_back({inName, v0JetsZ[0], v0JetsZ[1], v0JetsZ[2]});
  inputStrings.push_back({inName, w0JetsZ[0], w0JetsZ[1], w0JetsZ[2]});
  inputStrings.push_back({inName, eJetsZ[0],  eJetsZ[1],  eJetsZ[2]});
  inputStrings.push_back({inName, ptJetsZ[0], ptJetsZ[1], ptJetsZ[2]});

  inputStrings.push_back({saveName});
  compareZ(inputStrings, jetptmin, jetptmax, additionalLatex, doRatio);
}
void v0zmatched(double jetptmin, double jetptmax, bool doRatio = false)
{
  vector<vector<string> > inputStrings;
  string inName = "../../inputfiles/pythia/V0Study/weightedjetfinder.root";
  bool ptscheme = true;
  string saveName = "WeightedJetFinder";
  saveName += "-zv0";
  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += TString::Format("_jetpt%.0f-%.0f", jetptmin, jetptmax).Data();
  saveName += "-matched-all";
  saveName += ".pdf";

  vector<string> v0JetsZ = { "hV0JetZ", "Unweighted", "hV0JetMatched" };
  vector<string> w0JetsZ = { "hW0JetZ", "Weighted", "hW0JetMatched" };
  vector<string> eJetsZ  = { "hEJetZ",  "Ch+V0 jets - #it{p}_{V0} (weighted)", "hEJetMatched" };
  vector<string> ptJetsZ = { "hPtJetZ", "Ch+V0 jets - #it{p}_{T, V0} (weighted)", "hPtJetMatched" };

  string latexText = TString::Format("#it{p}_{T} scheme, %.0f #leq #it{p}_{T, jet} < %.0f GeV/#it{c}", jetptmin, jetptmax).Data();
  TLatex* additionalLatex = CreateLatex(0.25, 0.93, latexText.c_str(), 0.04);

  inputStrings.push_back({inName, v0JetsZ[0], v0JetsZ[1], v0JetsZ[2]});
  inputStrings.push_back({inName, w0JetsZ[0], w0JetsZ[1], w0JetsZ[2]});
  inputStrings.push_back({inName, eJetsZ[0],  eJetsZ[1],  eJetsZ[2]});
  inputStrings.push_back({inName, ptJetsZ[0], ptJetsZ[1], ptJetsZ[2]});

  inputStrings.push_back({saveName});
  compareZ(inputStrings, jetptmin, jetptmax, additionalLatex, doRatio);
}

void v0zweighted(double jetptmin, double jetptmax, bool doRatio = false)
{
  vector<vector<string> > inputStrings;
  string inName = "../../inputfiles/pythia/V0Study/weightedjetfinder.root";
  bool ptscheme = true;
  string saveName = "WeightedJetFinder";
  saveName += "-zv0";
  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += TString::Format("_jetpt%.0f-%.0f", jetptmin, jetptmax).Data();
  saveName += "-weighted";
  saveName += ".pdf";

  vector<string> w0JetsZ = { "hW0JetZ", "Weighted", "hW0Jet" };
  vector<string> eJetsZ  = { "hEJetZ",  "Ch+V0 jets - #it{p}_{V0} (weighted)", "hEJet" };
  vector<string> ptJetsZ = { "hPtJetZ", "Ch+V0 jets - #it{p}_{T, V0} (weighted)", "hPtJet" };

  string latexText = TString::Format("#it{p}_{T} scheme, %.0f #leq #it{p}_{T, jet} < %.0f GeV/#it{c}", jetptmin, jetptmax).Data();
  TLatex* additionalLatex = CreateLatex(0.25, 0.93, latexText.c_str(), 0.04);

  inputStrings.push_back({inName, w0JetsZ[0], w0JetsZ[1], w0JetsZ[2]});
  inputStrings.push_back({inName, eJetsZ[0],  eJetsZ[1],  eJetsZ[2]});
  inputStrings.push_back({inName, ptJetsZ[0], ptJetsZ[1], ptJetsZ[2]});

  inputStrings.push_back({saveName});
  compareZ(inputStrings, jetptmin, jetptmax, additionalLatex, doRatio);
}
void v0zweightedmatched(double jetptmin, double jetptmax, bool doRatio = false)
{
  vector<vector<string> > inputStrings;
  string inName = "../../inputfiles/pythia/V0Study/weightedjetfinder.root";
  bool ptscheme = true;
  string saveName = "WeightedJetFinder";
  saveName += "-zv0";
  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += TString::Format("_jetpt%.0f-%.0f", jetptmin, jetptmax).Data();
  saveName += "-matched-weighted";
  saveName += ".pdf";

  vector<string> w0JetsZ = { "hW0JetZ", "Weighted", "hW0JetMatched" };
  vector<string> eJetsZ  = { "hEJetZ",  "Ch+V0 jets - #it{p}_{V0} (weighted)", "hEJetMatched" };
  vector<string> ptJetsZ = { "hPtJetZ", "Ch+V0 jets - #it{p}_{T, V0} (weighted)", "hPtJetMatched" };

  string latexText = TString::Format("#it{p}_{T} scheme, %.0f #leq #it{p}_{T, jet} < %.0f GeV/#it{c}", jetptmin, jetptmax).Data();
  TLatex* additionalLatex = CreateLatex(0.25, 0.93, latexText.c_str(), 0.04);

  inputStrings.push_back({inName, w0JetsZ[0], w0JetsZ[1], w0JetsZ[2]});
  inputStrings.push_back({inName, eJetsZ[0],  eJetsZ[1],  eJetsZ[2]});
  inputStrings.push_back({inName, ptJetsZ[0], ptJetsZ[1], ptJetsZ[2]});

  inputStrings.push_back({saveName});
  compareZ(inputStrings, jetptmin, jetptmax, additionalLatex, doRatio);
}
