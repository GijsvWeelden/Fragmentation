
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
 * Plots V0 reconstruction efficiency inside and outside of jets
 */

const double MassK0S = 0.497611;
const double MassLambda = 1.115683;

void myPlotter(TCanvas* c, TH1F* frame, vector<TObject*> objList)
{
  c->cd();
  frame->Draw();
  for (auto obj : objList) {
    obj->Draw("same");
  }
  c->SaveAs(c->GetName());
}
// Returns rebinned version of a V0 pT histogram
// Necessary due to non-constant bin width for the output
TH1D* rebinV0Pt(TH1* h)
{
  string name = h->GetName();
  name += "_Rebinned";
  double ptBinEdges[12] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0};
  TH1D* hnew = new TH1D(name.c_str(), name.c_str(), 11, ptBinEdges);

  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double pt = h->GetBinCenter(i);
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    hnew->Fill(pt, bc);

    double w = hnew->GetBinError(hnew->FindBin(pt));
    hnew->SetBinError(hnew->FindBin(pt), sqrt(w*w + be*be)); // Sum of squares of weights
  }
  return hnew;
}
// Returns the scale for drawing histogram. Accounts for bin content and error
double getHistScale(TH1* h, bool doError)
{
  if (!doError) { return h->GetBinContent(h->GetMaximumBin()); }

  double scale = -900.;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    double s = be + bc;
    if (s > scale) { scale = s; }
  }
  return scale;
}

// Reconstruction efficiency of inclusive V0s
// Only plots one efficiency curve. Currently obsolete, but may be useful in the future
void plotV0Efficiency(array<string, 5> inputStrings)
{
  gStyle->SetNdivisions(505);
  string fileName        = inputStrings[0];
  string numeratorName   = inputStrings[1];
  string denominatorName = inputStrings[2];
  string dataSet         = inputStrings[3];
  string hadron          = inputStrings[4];
  if (("K0S" == hadron) + ("Lambda" == hadron) + ("AntiLambda" == hadron) != 1) {
    cout << "Invalid input! Hadron (" << hadron << ") should be K0S, Lambda or AntiLambda" << endl;
    return;
  }

  TFile* f = TFile::Open(fileName.c_str());
  TH3D* th3_det  = (TH3D*)f->Get(numeratorName.c_str());
  TH3D* th3_part = (TH3D*)f->Get(denominatorName.c_str());

  TH1D* numerator   = th3_det->ProjectionX();
  TH1D* denominator = th3_part->ProjectionX();

  TH1D* numRebinned = rebinV0Pt(numerator);
  TH1D* denRebinned = rebinV0Pt(denominator);
  TH1D* efficiency  = (TH1D*)numRebinned->Clone("efficiency");
  efficiency->SetName(TString::Format("efficiency%s", hadron.c_str()).Data());
  efficiency->Divide(denRebinned);
  setStyle(efficiency, 0);

  string canvasName = hadron;
  canvasName += "Efficiency";
  canvasName += ".pdf";
  TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);

  string xTitle = TString::Format("#it{p}_{T, %s} (GeV/#it{c})", hadron.c_str()).Data();
  string yTitle = "Efficiency";

  double xMinFrame = efficiency->GetXaxis()->GetXmin(), xMaxFrame = efficiency->GetXaxis()->GetXmax(),
         yMinFrame = 0., yMaxFrame = 1.1 * efficiency->GetBinContent(efficiency->GetMaximumBin());
  string text = dataSet;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle.c_str(), yTitle.c_str());
  frame->SetTitle(text.c_str());

  vector<TObject*> objList = {efficiency};
  myPlotter(c, frame, objList);
}
// Reconstruction efficiency of inclusive V0s
// Plots multiple efficiency curves
void plotV0Efficiency(vector<vector<string> > inputStrings, vector<TLatex*> textList, bool drawErr)
{
  gStyle->SetNdivisions(505);
  double textSize = 0.04;
  TLegend* legend = CreateLegend(0.5, 0.8, 0.6, 0.85, "", textSize);
  vector<TObject*> objList;

  double xMinFrame = 900., xMaxFrame = -900., yMinFrame = 0., yMaxFrame = -900.;
  for (int i = 0; i < inputStrings.size(); i++) {
    string fileName        = inputStrings[i][0];
    string numeratorName   = inputStrings[i][1];
    string denominatorName = inputStrings[i][2];
    string legendEntry     = inputStrings[i][3];

    TFile* f = TFile::Open(fileName.c_str());
    TH3D* th3_det  = (TH3D*)f->Get(numeratorName.c_str());
    TH3D* th3_part = (TH3D*)f->Get(denominatorName.c_str());

    TH1D* numerator   = th3_det->ProjectionX();
    TH1D* denominator = th3_part->ProjectionX();

    TH1D* numRebinned = rebinV0Pt(numerator);
    TH1D* denRebinned = rebinV0Pt(denominator);
    TH1D* efficiency  = (TH1D*)numRebinned->Clone("efficiency");
    efficiency->Divide(denRebinned);
    efficiency->SetName(TString::Format("efficiency%s", legendEntry.c_str()).Data());
    setStyle(efficiency, i);

    double xmin = efficiency->GetXaxis()->GetXmin();
    double xmax = efficiency->GetXaxis()->GetXmax();
    double ymax = getHistScale(efficiency, drawErr);

    if (xmin < xMinFrame) { xMinFrame = xmin; }
    if (xmax > xMaxFrame) { xMaxFrame = xmax; }
    if (ymax > yMaxFrame) { yMaxFrame = ymax; }
    objList.push_back(efficiency);
    legend->AddEntry(efficiency, legendEntry.c_str());
  }
  objList.push_back(legend);
  for (auto l : textList) { objList.push_back(l); }

  string canvasName = "V0Efficiency";
  canvasName += ".pdf";
  TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);
  string xTitle = "#it{p}_{T} (GeV/#it{c})";
  string yTitle = "Efficiency";
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, 1.1*yMaxFrame, xTitle.c_str(), yTitle.c_str());

  c->cd();
  frame->Draw();
  string drawOption = "same";
  if (!drawErr) { drawOption += " hist"; }
  for (auto obj : objList) {
    obj->Draw(drawOption.c_str());
  }
  c->SaveAs(c->GetName());
}
// Reconstruction efficiency of V0s in jets
void plotV0Efficiency(vector<vector<string> > inputStrings, double jetptmin, double jetptmax, vector<TLatex*> textList, bool drawErr)
{
  gStyle->SetNdivisions(505);
  double textSize = 0.04;
  TLegend* legend = CreateLegend(0.5, 0.8, 0.6, 0.85, "", textSize);
  vector<TObject*> objList;
  string lowptText, highptText;

  double xMinFrame = 0., xMaxFrame = jetptmax, yMinFrame = 0., yMaxFrame = -900.;
  for (int i = 0; i < inputStrings.size(); i++) {
    string fileName        = inputStrings[i][0];
    string numeratorName   = inputStrings[i][1];
    string denominatorName = inputStrings[i][2];
    string legendEntry     = inputStrings[i][3];

    array<int, 2> jetPtBins;
    TFile* f = TFile::Open(fileName.c_str());
    TH3D* th3_det  = (TH3D*)f->Get(numeratorName.c_str());
    TH3D* th3_part = (TH3D*)f->Get(denominatorName.c_str());

    jetPtBins = getProjectionBins(th3_det->GetXaxis(), jetptmin, jetptmax);
    TH1D* numerator   = th3_det->ProjectionY("det", jetPtBins[0], jetPtBins[1], 0, th3_det->GetNbinsZ()+1);
    TH1D* numRebinned = rebinV0Pt(numerator);

    jetPtBins = getProjectionBins(th3_part->GetXaxis(), jetptmin, jetptmax);
    TH1D* denominator = th3_part->ProjectionZ("part", jetPtBins[0], jetPtBins[1], 0, th3_part->GetNbinsY()+1);
    TH1D* denRebinned = rebinV0Pt(denominator);

    double lowpt  = th3_part->GetXaxis()->GetBinLowEdge(jetPtBins[0]);
    double highpt = th3_part->GetXaxis()->GetBinUpEdge(jetPtBins[1]);
    lowptText  = TString::Format("%.0f", lowpt).Data();
    highptText = TString::Format("%.0f", highpt).Data();

    TH1D* efficiency  = (TH1D*)numRebinned->Clone("efficiency");
    efficiency->SetName(TString::Format("efficiency%s_jetpt%s-%s", legendEntry.c_str(), lowptText.c_str(), highptText.c_str()).Data());
    efficiency->Divide(denRebinned);
    setStyle(efficiency, i);

    if (highpt > xMaxFrame) { xMaxFrame = highpt; }
    double ymax = getHistScale(efficiency, drawErr);
    if (ymax > yMaxFrame)   { yMaxFrame = ymax; }
    objList.push_back(efficiency);
    legend->AddEntry(efficiency, legendEntry.c_str());
  }
  objList.push_back(legend);
  // Add jet pT range to the first latex
  if (textList.size() > 0) {
    string title = textList[0]->GetTitle();
    title += TString::Format(", #it{p}_{T} = %s-%s GeV/#it{c}", lowptText.c_str(), highptText.c_str()).Data();
    textList[0]->SetTitle(title.c_str());
  }
  for (auto l : textList) { objList.push_back(l); }

  string canvasName = "V0Efficiency";
  canvasName += TString::Format("_jetpt%s-%s", lowptText.c_str(), highptText.c_str()).Data();
  canvasName += ".pdf";
  TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);

  string xTitle = "#it{p}_{T} (GeV/#it{c})";
  string yTitle = "Efficiency";
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, 1.1*yMaxFrame, xTitle.c_str(), yTitle.c_str());

  c->cd();
  frame->Draw();
  string drawOption = "same";
  if (!drawErr) { drawOption += " hist"; }
  for (auto obj : objList) {
    obj->Draw(drawOption.c_str());
  }
  c->SaveAs(c->GetName());
}

void plot24b1bIncl(string hadron, bool drawErr = true)
{
  string fileName = "~/cernbox/TrainOutput/267466/AnalysisResults.root";
  string detName  = TString::Format("jet-v0qa_id13556/InvMass%sTrue", hadron.c_str()).Data();
  string partName = TString::Format("jet-v0qa_id13556/Generated%s", hadron.c_str()).Data();
  string dataSet  = "LHC24b1b";
  vector<string> inputStrings = {fileName, detName, partName, formatHadronName(hadron).c_str()};
  vector<TLatex*> textList = { CreateLatex(0.2, 0.95, dataSet.c_str()) };
  plotV0Efficiency({inputStrings}, textList, drawErr);
}
void plot24b1bIncl(bool drawErr = true)
{
  double textSize = 0.04;
  string fileName       = "~/cernbox/TrainOutput/267466/AnalysisResults.root";
  string dataSet        = "LHC24b1b";
  string detK0S         = "jet-v0qa_id13556/InvMassK0STrue";
  string partK0S        = "jet-v0qa_id13556/GeneratedK0S";
  string detLambda      = "jet-v0qa_id13556/InvMassLambdaTrue";
  string partLambda     = "jet-v0qa_id13556/GeneratedLambda";
  string detAntiLambda  = "jet-v0qa_id13556/InvMassAntiLambdaTrue";
  string partAntiLambda = "jet-v0qa_id13556/GeneratedAntiLambda";

  vector<vector<string> > inputStrings = {
    { fileName, detK0S, partK0S, formatHadronName("K0S")} ,
    { fileName, detLambda, partLambda, formatHadronName("Lambda")} ,
    { fileName, detAntiLambda, partAntiLambda, formatHadronName("AntiLambda") }
  };
  vector<TLatex*> textList = { CreateLatex(0.2, 0.95, dataSet.c_str(), textSize) };
  plotV0Efficiency(inputStrings, textList, drawErr);
}
void plot24b1bJet(string hadron, double jetptmin, double jetptmax, bool drawErr = true)
{
  double textSize = 0.04;
  string fileName = "~/cernbox/TrainOutput/267466/AnalysisResults.root";
  string dataSet  = "LHC24b1b";
  string detK0S  = "jet-v0qa_id17174/InvMassJetK0STrue";
  string partK0S = "jet-v0qa_id17174/GeneratedJetK0S";
  string detLambda  = "jet-v0qa_id17174/InvMassJetLambdaTrue";
  string partLambda = "jet-v0qa_id17174/GeneratedJetLambda";
  string detAntiLambda  = "jet-v0qa_id17174/InvMassJetAntiLambdaTrue";
  string partAntiLambda = "jet-v0qa_id17174/GeneratedJetAntiLambda";

  vector<vector<string> > inputStrings = {
    { fileName, detK0S, partK0S, formatHadronName("K0S")} ,
    { fileName, detLambda, partLambda, formatHadronName("Lambda")} ,
    { fileName, detAntiLambda, partAntiLambda, formatHadronName("AntiLambda") }
  };
  vector<TLatex*> textList = { CreateLatex(0.2, 0.95, dataSet.c_str(), textSize) };

  plotV0Efficiency(inputStrings, jetptmin, jetptmax, textList, drawErr);
}
