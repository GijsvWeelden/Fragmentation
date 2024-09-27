
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
 * Plots for V0s in the perpendicular cones
 */

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

void myPlotter(TCanvas* c, TH1F* frame, vector<TObject*> objList)
{
  c->cd();
  frame->Draw();
  for (auto obj : objList) {
    obj->Draw("same");
  }
  c->SaveAs(c->GetName());
}

void plotMass(array<string, 4> inputStrings, double ptmin, double ptmax, bool isTrueV0)
{
  gStyle->SetNdivisions(505);
  string fileName = inputStrings[0];
  string histName = inputStrings[1];
  string dataSet  = inputStrings[2];
  string hadron   = inputStrings[3];
  if (("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1) {
    cout << "Invalid input! Hadron (" << hadron << ") should be K0S, Lambda0 or AntiLambda0" << endl;
    return;
  }

  const int nDim           = 4;
  const int ptAxis         = 0;
  const int K0SAxis        = 1;
  const int LambdaAxis     = 2;
  const int AntiLambdaAxis = 3;

  int projectionAxis = K0SAxis;
  if ("Lambda0" == hadron) projectionAxis = LambdaAxis;
  if ("AntiLambda0" == hadron) projectionAxis = AntiLambdaAxis;

  TFile* f = TFile::Open(fileName.c_str());
  THnSparse* thn = (THnSparse*)f->Get(histName.c_str());

  array<int, 2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
  thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
  double lowpt  = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
  double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
  TH1D* mass = thn->Projection(projectionAxis);
  mass->SetName(Form("mass_%s_%s_pt%.1f-%.1f", hadron.c_str(), histName.c_str(), lowpt, highpt));
  setStyle(mass, 0);

  TLine* massLine = new TLine(MassK0S, 0., MassK0S, 1.1*mass->GetMaximum());
  massLine->SetLineColor(GetColor(1));
  massLine->SetLineStyle(9);
  massLine->SetLineWidth(3);

  string canvasName = "PCV0" + (isTrueV0 ? "True" : "Fake");
  canvasName += TString::Format("_m%s_pt%.1f-%.1f", hadron.c_str(), lowpt, highpt).Data();
  canvasName += ".pdf";
  TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);

  string xTitle = "#it{M}("; xTitle += formatHadronDaughters(hadron).c_str(); xTitle += ") (GeV/#it{c}^{2})";
  string yTitle = "Counts";

  double xMinFrame = mass->GetXaxis()->GetXmin(), xMaxFrame = mass->GetXaxis()->GetXmax(),
         yMinFrame = 0., yMaxFrame = 1.1 * mass->GetBinContent(mass->GetMaximumBin());
  string text = dataSet;
  text += TString::Format(", #it{p}_{T} = %.1f-%.1f GeV/#it{c}", lowpt, highpt).Data();
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle.c_str(), yTitle.c_str());
  frame->SetTitle(text.c_str());

  vector<TObject*> objList = {mass, massLine};
  for (auto l : latex) { objList.push_back(l); }
  myPlotter(c, frame, objList);
}

void plot24b1bMatched(string hadron, double ptmin, double ptmax)
{
  string fileName = "~/cernbox/TrainOutput/267466/AnalysisResults.root";
  string histName = "jet-fragmentation/matching/PC/matchedV0PtMass";
  string dataSet  = "LHC24b1b";
  array<string, 4> inputStrings = {fileName, histName, dataSet, hadron};
  plotMass(inputStrings, ptmin, ptmax, true);
}
void plot24b1bFake(string hadron, double ptmin, double ptmax)
{
  string fileName = "~/cernbox/TrainOutput/267466/AnalysisResults.root";
  string histName = "jet-fragmentation/matching/PC/fakeV0PtMass";
  string dataSet  = "LHC24b1b";
  array<string, 4> inputStrings = {fileName, histName, dataSet, hadron};
  TLatex* latex = new TLatex(0.15, 0.85, "Combinatorial V0");
  plotMass(inputStrings, ptmin, ptmax, false);
}