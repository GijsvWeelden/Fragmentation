
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

double getNevts(string inName)
{
  TFile* inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  string histName = "jet-fragmentation/data/V0/nV0sEvent";
  TH1D* h0 = (TH1D*)inFile->Get(histName.c_str());
  TH1D* hNEvents = (TH1D*)h0->Clone();
  if ("~/Documents/TrainOutput/211540/AnalysisResults.root" == inName) {
    TFile* f1 = new TFile("~/Documents/TrainOutput/211540/449460.root");
    TH1D* h1 = (TH1D*)f1->Get(histName.c_str());
    TFile* f2 = new TFile("~/Documents/TrainOutput/211540/449462.root");
    TH1D* h2 = (TH1D*)f2->Get(histName.c_str());
    TFile* f3 = new TFile("~/Documents/TrainOutput/211540/449466.root");
    TH1D* h3 = (TH1D*)f3->Get(histName.c_str());

    hNEvents->Add(h1);
    hNEvents->Add(h2);
    hNEvents->Add(h3);
  }
  return hNEvents->Integral();
}

void comparePtData(string dataSet, vector<string> inNames, vector<string> legendEntries, bool doRatio = false)
{
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 200., yMinFrame = 1e-11, yMaxFrame = 1e-3;
  double xMinLegend = 0.3, xMaxLegend = 0.9, yMinLegend = 0.4, yMaxLegend = 0.55;
  double xLatex = 0.3, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet}";
  yTitle = "normalised count";
  yTitle = "#it{N}_{jets} / #it{N}_{evts}";
  drawoption = "same";
  saveName = "ptComparison";
  histName = "jet-fragmentation/data/jets/jetPtEtaPhi";

  vector<double> integrals; vector<double> nevents;

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  latexText = TString::Format("%s", dataSet.c_str()).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  for (int i = 0; i < inNames.size(); i++) {
    TH1D* th1;
    string inName = inNames[i];
    TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());

    if (inName.find("203281") != std::string::npos) {
      histName = "jet-fragmentation/data/jets/V0/jetPtnV0";
      TH2D* th2 = (TH2D*)inFile->Get(histName.c_str());
      th1 = th2->ProjectionX(TString::Format("pt_%s", inName.c_str()).Data(), 2, th2->GetYaxis()->GetNbins());
    }
    else {
      histName = "jet-fragmentation/data/jets/jetPtEtaPhi";
      // TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());
      TH3D* h0 = (TH3D*)inFile->Get(histName.c_str());

      TFile* f1 = new TFile("~/Documents/TrainOutput/211540/449460.root");
      TH3D* h1 = (TH3D*)f1->Get(histName.c_str());
      TFile* f2 = new TFile("~/Documents/TrainOutput/211540/449462.root");
      TH3D* h2 = (TH3D*)f2->Get(histName.c_str());
      TFile* f3 = new TFile("~/Documents/TrainOutput/211540/449466.root");
      TH3D* h3 = (TH3D*)f3->Get(histName.c_str());

      TH3D* th3 = (TH3D*)h0->Clone();
      th3->Add(h1);
      th3->Add(h2);
      th3->Add(h3);

      th1 = th3->ProjectionX(TString::Format("pt_%s", inName.c_str()).Data());
    }

    setStyle(th1, i);
    legend->AddEntry(th1, legendEntries[i].c_str());
    th1->Scale(1./getNevts(inName));

    histVector.push_back(th1);
    integrals.push_back(th1->Integral());
    nevents.push_back(getNevts(inName));
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.5; yMaxFrame = 1.1;
    // drawoption = "same hist p";
    saveName = "ptRatio";
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
      // hist->Scale(1./nevents[i]);
    }
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");
  latex->Draw("same");

  TLatex* integralsLatex0 = CreateLatex(0.4, 0.8, TString::Format("#splitline{%.2g ch. jets}{ %.2g events}", integrals[0], nevents[0]).Data(), textSize);
  integralsLatex0->Draw("same");
  TLatex* integralsLatex1 = CreateLatex(0.4, 0.7, TString::Format("#splitline{%.2g ch.+V0 jets}{ %.2g events}", integrals[1], nevents[1]).Data(), textSize);
  integralsLatex1->Draw("same");

  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

void compareZData(string dataSet, vector<string> inNames, vector<string> legendEntries, string hadron, double jetptmin = 10., double jetptmax = 1e3, bool doRatio = false, bool doCounts = false)
{
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }

  const int nDim                = 5;
  const int jetptAxis           = 0;
  const int v0zAxis             = 1;
  const int K0SMassAxis         = 2;
  const int Lambda0MassAxis     = 3;
  const int AntiLambda0MassAxis = 4;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double lowjetpt, highjetpt;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 1e-10, yMaxFrame = 1e-5;
  double xMinLegend = 0.5, xMaxLegend = 0.8, yMinLegend = 0.74, yMaxLegend = 0.89;
  // double xMinLegend = 0.25, xMaxLegend = 0.8, yMinLegend = 0.175, yMaxLegend = 0.325;
  double xLatex = 0.15, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 4;
  xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data();
  // yTitle = TString::Format("#it{N}_{%s}/#it{N}_{evts}", formatHadronName(hadron).c_str()).Data();
  yTitle = TString::Format("#frac{1}{#it{N}_{evts}} #frac{d #it{N}}{d #it{z}_{%s}}", formatHadronName(hadron).c_str()).Data();
  if (doCounts) { yTitle = "Counts"; }
  drawoption = "same";
  saveName = TString::Format("z%sComparison", hadron.c_str()).Data();

  string hadronForHistname = hadron;
  if ("Lambda0" == hadron) { hadronForHistname = "Lambda"; }
  if ("AntiLambda0" == hadron) { hadronForHistname = "AntiLambda"; }
  histName = "jet-fragmentation/data/jets/V0";
  histName = TString::Format("%s/jetPt%sTrackProjAllMasses", histName.c_str(), hadronForHistname.c_str()).Data();

  vector<double> integrals; vector<double> nevents;
  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  for (int i = 0; i < inNames.size(); i++) {
    TFile *inFile = TFile::Open(TString::Format("%s", inNames[i].c_str()).Data());
    THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

    if ("~/Documents/TrainOutput/211540/AnalysisResults.root" == inNames[i]) {
      THnSparseD* h0 = (THnSparseD*)inFile->Get(histName.c_str());
      TFile* f1 = new TFile("~/Documents/TrainOutput/211540/449460.root");
      THnSparseD* h1 = (THnSparseD*)f1->Get(histName.c_str());
      TFile* f2 = new TFile("~/Documents/TrainOutput/211540/449462.root");
      THnSparseD* h2 = (THnSparseD*)f2->Get(histName.c_str());
      TFile* f3 = new TFile("~/Documents/TrainOutput/211540/449466.root");
      THnSparseD* h3 = (THnSparseD*)f3->Get(histName.c_str());

      thn = (THnSparseD*)h0->Clone();
      thn->Add(h1);
      thn->Add(h2);
      thn->Add(h3);

      h0->Print();
      h1->Print();
      h2->Print();
      h3->Print();
      thn->Print();
    }

    std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
    thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
    lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
    highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);

    TH1D* th1 = thn->Projection(v0zAxis);
    th1->SetName(TString::Format("z%s_%s", hadron.c_str(), inNames[i].c_str()).Data());
    th1->Rebin(rebinNumber);
    setStyle(th1, i);
    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);
    integrals.push_back(th1->Integral());
    nevents.push_back(getNevts(inNames[i]));
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.; yMaxFrame = 1.1;
    // drawoption = "same hist p";
    saveName = TString::Format("z%sRatio", hadron.c_str()).Data();
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
    if (!doRatio && !doCounts) {
      hist->Scale(1./nevents[i], "width");
    }
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");

  latexText = TString::Format("%s, #it{p}_{T, jet} = %.0f - %.0f GeV/#it{c}", dataSet.c_str(), lowjetpt, highjetpt).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex->Draw("same");

  TLatex* integralsLatex0 = CreateLatex(0.25, 0.3, TString::Format("#splitline{%g %s #in ch. jets}{ %g events}", integrals[0], formatHadronName(hadron).c_str(), nevents[0]).Data(), textSize);
  integralsLatex0->Draw("same");
  TLatex* integralsLatex1 = CreateLatex(0.25, 0.2, TString::Format("#splitline{%g %s #in ch.+V0 jets}{ %g events}", integrals[1], formatHadronName(hadron).c_str(), nevents[1]).Data(), textSize);
  integralsLatex1->Draw("same");

  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

void compare22o(bool doRatio = false, string hadron = "", double jetptmin = 10., double jetptmax = 1e3, bool doCounts = false)
{
  vector<string> inNames; vector<string> legendEntries;
  string dataSet = "LHC22o_pass6_minBias";

  // array<string, 2> file203281 = {"~/Documents/TrainOutput/203281/AnalysisResults.root", "Ch. jets"};
  array<string, 2> file203281 = {"~/Documents/TrainOutput/203281/AnalysisResults.root", "Ch. jets w. V0"};
  inNames.push_back(file203281[0]);
  legendEntries.push_back(file203281[1]);

  // array<string, 2> file210823 = {"~/Documents/TrainOutput/210823/AnalysisResults.root", "Ch+V0 jets (small dataset)"};
  // inNames.push_back(file210823[0]);
  // legendEntries.push_back(file210823[1]);

  array<string, 2> file211540 = {"~/Documents/TrainOutput/211540/AnalysisResults.root", "Ch+V0 jets"};
  inNames.push_back(file211540[0]);
  legendEntries.push_back(file211540[1]);

  if ("" == hadron) {
    comparePtData(dataSet, inNames, legendEntries, doRatio);
  }
  else {
    compareZData(dataSet, inNames, legendEntries, hadron, jetptmin, jetptmax, doRatio, doCounts);
  }
}

