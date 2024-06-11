
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

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

double getNjets(TFile* inFile, double jetptmin, double jetptmax)
{
  string histName = "jet-fragmentation/data/jets/jetPtEtaPhi";
  TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
  std::array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
  return jets->Integral(jetptbins[0], jetptbins[1], 1, jets->GetNbinsY(), 1, jets->GetNbinsZ());
}

// ----------------------------------------------------------

void plotK0SMass(double jetptmin = 10., double jetptmax = 1e3, double zmin = 0., double zmax = 1.+1e-3)
{
  string inName = "~/Documents/TrainOutput/211540/AnalysisResults.root";
  string dataSet = "LHC22o_pass6_minBias";

  string hadron = "K0S";
  const int nDim = 5;
  const int jetptAxis = 0;
  const int zAxis = 1;
  const int K0SMassAxis = 2;
  const int Lambda0MassAxis = 3;
  const int AntiLambda0MassAxis = 4;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  int rebinNumber = 5;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1e-4, yMaxFrame = .1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.2, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  xTitle = TString::Format("#it{M} (%s)", formatHadronName(hadron).c_str()).Data();
  yTitle = "normalised count";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  histName = "jet-fragmentation/data/jets/V0";
  histName = TString::Format("%s/jetPt%sTrackProjAllMasses", histName.c_str(), hadron.c_str()).Data();
  cout << "histName: " << histName << endl;
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  inFile->Print();
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0zbins = getProjectionBins(thn->GetAxis(zAxis), zmin, zmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(zAxis)->SetRange(v0zbins[0], v0zbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0z = thn->GetAxis(zAxis)->GetBinLowEdge(v0zbins[0]);
  double highv0z = thn->GetAxis(zAxis)->GetBinUpEdge(v0zbins[1]);

  TH1D* mass = (TH1D*)thn->Projection(K0SMassAxis);
  mass->SetName(TString::Format("K0Smass_jetpt%.0f-%.0f_z%.1f-%.1f", lowjetpt, highjetpt, lowv0z, highv0z).Data());
  mass->Rebin(rebinNumber);
  mass->Scale(1./mass->Integral());
  setStyle(mass, 0);
  histVector.push_back(mass);

  latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{z}_{%s} = %.1f - %.1f } }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0z, highv0z).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%smass", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s_v0z%.1f-%.1f", saveName.c_str(), lowv0z, highv0z);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}

void plotV0Mass(string hadron, double jetptmin = 10., double jetptmax = 1e3, double zmin = 0., double zmax = 1.+1e-3)
{
  string inName = "~/Documents/TrainOutput/211540/AnalysisResults.root";
  string dataSet = "LHC22o_pass6_minBias";

  if (hadron != "K0S" && hadron != "Lambda0" && hadron != "AntiLambda0") {
    cout << "Error: hadron should be K0S, Lambda0 or AntiLambda0" << endl;
    return;
  }
  const int nDim = 5;
  const int jetptAxis = 0;
  const int zAxis = 1;
  const int K0SMassAxis = 2;
  const int Lambda0MassAxis = 3;
  const int AntiLambda0MassAxis = 4;

  int projectionAxis = ("K0S" == hadron)*K0SMassAxis + ("Lambda0" == hadron)*Lambda0MassAxis + ("AntiLambda0" == hadron)*AntiLambda0MassAxis;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = false;
  int rebinNumber = 5;
  double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 0.1, yMaxFrame = 80;
  // double xMinFrame = 0.4, xMaxFrame = 0.6, yMinFrame = 1e-4, yMaxFrame = .1;
  if ("Lambda0" == hadron || "AntiLambda0" == hadron) {
    // xMinFrame = 1.05, xMaxFrame = 1.215;
    xMinFrame = 1.1, xMaxFrame = 1.15;
    yMinFrame = 0.1, yMaxFrame = 110.;
    rebinNumber = 2;
  }
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.25, yLatex = 0.82;
  int xCanvas = 900, yCanvas = 900;
  xTitle = TString::Format("#it{M} (%s)", formatHadronName(hadron).c_str()).Data();
  // yTitle = TString::Format("#frac{1}{#it{N}_{%s}} #frac{d#it{N}}{d#it{M}}" , formatHadronName(hadron).c_str()).Data();
  // yTitle = "#frac{1}{#it{N}} #frac{d#it{N}}{d#it{M}}";
  yTitle = "probability density";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  TLine* line = new TLine(("K0S" == hadron) ? MassK0S : MassLambda0, yMinFrame, ("K0S" == hadron) ? MassK0S : MassLambda0, 0.5*yMaxFrame);
  line->SetLineColor(GetColor(0));
  line->SetLineWidth(2);
  line->SetLineStyle(9);

  histName = "jet-fragmentation/data/jets/V0";
  histName = TString::Format("%s/jetPtV0TrackProjMass", histName.c_str()).Data();
  TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());
  thn->Sumw2();

  std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
  std::array<int, 2> v0zbins = getProjectionBins(thn->GetAxis(zAxis), zmin, zmax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(zAxis)->SetRange(v0zbins[0], v0zbins[1]);
  double lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
  double highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
  double lowv0z = thn->GetAxis(zAxis)->GetBinLowEdge(v0zbins[0]);
  double highv0z = thn->GetAxis(zAxis)->GetBinUpEdge(v0zbins[1]);
  if (lowv0z < 0) { lowv0z = 0.; } // Prevents z = -0 in latex

  TH1D* mass = (TH1D*)thn->Projection(projectionAxis);
  mass->SetName(TString::Format("%smass_jetpt%.0f-%.0f_z%.1f-%.1f", hadron.c_str(), lowjetpt, highjetpt, lowv0z, highv0z).Data());
  mass->Rebin(rebinNumber);
  mass->Scale(1./mass->Integral(), "width");
  setStyle(mass, 0);
  histVector.push_back(mass);

  // latexText = TString::Format("#splitline{ %s }{ #splitline{ #it{p}_{T, ch. jet} = %.0f - %.0f GeV/c }{ #it{z}_{%s} = %.1f - %.1f } }", dataSet.c_str(), lowjetpt, highjetpt, formatHadronName(hadron).c_str(), lowv0z, highv0z).Data();
  latexText = TString::Format("#splitline{ %s }{#splitline{ %s }{#splitline{ %s }{ %s }}}",
                                "ALICE 2022 data pp, #sqrt{#it{s}} = 13.6 TeV",
                                "Anti-k_{T} jets, |#it{#eta}_{jet}| < 0.75 - #it{R}",
                                TString::Format("%.0f < #it{p}_{T, ch.+V0 jet} < %.0f GeV/#it{c}", lowjetpt, highjetpt).Data(),
                                ("K0S" == hadron) ? TString::Format("%.1f < #it{z}_{V0} < %.1f", lowv0z, highv0z).Data() : TString::Format("                    %.1f < #it{z}_{V0} < %.1f", lowv0z, highv0z).Data() // 20 spaces
                              );
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("%smass", hadron.c_str()).Data();
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s_v0z%.1f-%.1f", saveName.c_str(), lowv0z, highv0z);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  // plotNHists(canvas, frame, histVector, legend, latex, saveName, "");

  frame->Draw();
  for (auto& hist : histVector) {
    hist->Draw("same");
  }
  line->Draw("same");
  latex->Draw("same");
  canvas->SaveAs(saveName.c_str());
}
