
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
double getNevts(string inName, bool isData = true)
{
  TFile* inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  TH1D* hNEvents = (TH1D*)inFile->Get(TString::Format("jet-fragmentation/%s/V0/nV0sEvent", isData ? "data" : "matching").Data());
  return hNEvents->Integral();
}

double getNevtsData(string inName)
{
  return getNevts(inName, true);
  // TFile* inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  // TH1D* hNEvents = (TH1D*)inFile->Get("jet-fragmentation/data/V0/nV0sEvent");
  // return hNEvents->Integral();
}

double getNevtsMc(string inName)
{
  return getNevts(inName, false);
  // TFile* inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  // TH1D* hNEvents = (TH1D*)inFile->Get("jet-fragmentation/matching/V0/nV0sEvent");
  // return hNEvents->Integral();
}

double getNevtsEfficiency(string inName)
{
  TFile* inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
  TH1D* hNEvents = (TH1D*)inFile->Get("track-efficiency/hMcCollCutsCounts");
  return hNEvents->GetBinContent(1);
}

void compareDataPt(string dataSet, vector<string> inNames, vector<string> legendEntries, bool doRatio = false)
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
      TH3D* th3 = (TH3D*)inFile->Get(histName.c_str());
      th1 = th3->ProjectionX(TString::Format("pt_%s", inName.c_str()).Data());
    }

    setStyle(th1, i);
    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);
    integrals.push_back(th1->Integral());
    nevents.push_back(getNevtsData(inName));
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
      hist->Scale(1./nevents[i]);
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

void compareDataZ(string dataSet, vector<string> inNames, vector<string> legendEntries, string hadron, double jetptmin = 10., double jetptmax = 1e3, bool doRatio = false, bool doCounts = false)
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
    nevents.push_back(getNevtsData(inNames[i]));
    // legend->AddEntry(th1, TString::Format("#splitline{%g %s #in %s}{  in %g events}", th1->Integral(), formatHadronName(hadron).c_str(), legendEntries[i].c_str(), getNevtsData(inNames[i])).Data());
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

void compareMcZ(string dataSet, vector<string> inNames, vector<string> legendEntries, string hadron, double jetptmin = 10., double jetptmax = 1e3, bool doParticle = false, bool doRatio = false, bool doCounts = false)
{
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }

  const int nDim                = 7;
  const int partjetptAxis       = 0;
  const int partv0zAxis         = 1;
  const int detjetptAxis        = 2;
  const int detv0zAxis          = 3;
  const int K0SMassAxis         = 4;
  const int Lambda0MassAxis     = 5;
  const int AntiLambda0MassAxis = 6;

  int jetptAxis = doParticle ? partjetptAxis : detjetptAxis;
  int v0zAxis   = doParticle ? partv0zAxis : detv0zAxis;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double lowjetpt, highjetpt;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 1e-9, yMaxFrame = 1e-2;
  double xMinLegend = 0.5, xMaxLegend = 0.8, yMinLegend = 0.74, yMaxLegend = 0.89;
  // double xMinLegend = 0.25, xMaxLegend = 0.8, yMinLegend = 0.175, yMaxLegend = 0.325;
  double xLatex = 0.3, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 4;
  xTitle = TString::Format("#it{z}_{%s}", formatHadronName(hadron).c_str()).Data();
  // yTitle = TString::Format("#it{N}_{%s}/#it{N}_{evts}", formatHadronName(hadron).c_str()).Data();
  yTitle = TString::Format("#frac{1}{#it{N}_{evts}} #frac{d #it{N}}{d #it{z}_{%s}}", formatHadronName(hadron).c_str()).Data();
  if (doCounts) { yTitle = "Counts"; }
  drawoption = "same";
  saveName = TString::Format("z%sComparison", hadron.c_str()).Data();

  string hadronForHistname = hadron;
  // if ("Lambda0" == hadron) { hadronForHistname = "Lambda"; }
  // if ("AntiLambda0" == hadron) { hadronForHistname = "AntiLambda"; }
  histName = "jet-fragmentation/matching/jets/V0";
  histName = TString::Format("%s/partJetPt%sTrackProjDetJetPt%sTrackProjAllMasses", histName.c_str(), hadronForHistname.c_str(), hadronForHistname.c_str()).Data();

  vector<double> integrals; vector<double> nevents;
  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  for (int i = 0; i < inNames.size(); i++) {
    TFile *inFile = TFile::Open(TString::Format("%s", inNames[i].c_str()).Data());
    THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

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
    nevents.push_back(getNevtsMc(inNames[i]));
    // legend->AddEntry(th1, TString::Format("#splitline{%g %s #in %s}{  in %g events}", th1->Integral(), formatHadronName(hadron).c_str(), legendEntries[i].c_str(), getNevts(inNames[i])).Data());
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

  latexText = TString::Format("%s, #it{p}_{T, jet} = %.0f - %.0f GeV/#it{c}", doParticle ? "particle-level" : "detector-level", lowjetpt, highjetpt).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex->Draw("same");

  TLatex* integralsLatex0 = CreateLatex(0.25, 0.3, TString::Format("#splitline{%s: %g %s in jets}{ %g events}", legendEntries[0].c_str(), integrals[0], formatHadronName(hadron).c_str(), nevents[0]).Data(), textSize);
  integralsLatex0->Draw("same");
  TLatex* integralsLatex1 = CreateLatex(0.25, 0.2, TString::Format("#splitline{%s: %g %s in jets}{ %g events}", legendEntries[1].c_str(), integrals[1], formatHadronName(hadron).c_str(), nevents[1]).Data(), textSize);
  integralsLatex1->Draw("same");

  saveName = TString::Format("%s_%s", saveName.c_str(), doParticle ? "part" : "det");
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

void compareMcTrackingEfficiency(vector<string> inNames, vector<string> legendEntries)
{
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double lowjetpt, highjetpt;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = false;
  double xMinFrame = 0., xMaxFrame = 10., yMinFrame = 0., yMaxFrame = 1.1;
  double xMinLegend = 0.5, xMaxLegend = 0.8, yMinLegend = 0.74, yMaxLegend = 0.89;
  double xLatex = 0.15, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 4;
  xTitle = "#it{p}_{T}";
  yTitle = "detector-level / particle-level";
  // yTitle = "#frac{1}{#it{N}_{evts}} #frac{d #it{N}}{d #it{p}_{T}}";
  drawoption = "same";

  string histDir = "track-efficiency";
  string histNameDet = TString::Format("%s/h3_track_pt_track_eta_track_phi_associatedtrack_primary", histDir.c_str()).Data();
  // string histNamePart = TString::Format("%s/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary", histDir.c_str()).Data();
  string histNamePart = TString::Format("%s/h3_particle_pt_particle_eta_particle_phi_mcpartofinterest", histDir.c_str()).Data();

  vector<double> integrals; vector<double> nevents;
  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  for (int i = 0; i < inNames.size(); i++) {
    TFile *inFile = TFile::Open(TString::Format("%s", inNames[i].c_str()).Data());
    TH3D* th3_det = (TH3D*)inFile->Get(histNameDet.c_str());
    TH3D* th3_part = (TH3D*)inFile->Get(histNamePart.c_str());

    TH1D* th1_det = th3_det->ProjectionX();
    TH1D* th1_part = th3_part->ProjectionX();

    th1_det->Rebin(rebinNumber);
    th1_part->Rebin(rebinNumber);

    TH1D* efficiency = (TH1D*)th1_det->Clone(TString::Format("efficiency_%s", inNames[i].c_str()).Data());
    efficiency->Divide(th1_part);

    setStyle(efficiency, i);
    legend->AddEntry(efficiency, TString::Format("%s: %g events", legendEntries[i].c_str(), getNevtsEfficiency(inNames[i])).Data());
    histVector.push_back(efficiency);
    integrals.push_back(efficiency->Integral());
    nevents.push_back(getNevtsEfficiency(inNames[i]));
    // legend->AddEntry(th1, TString::Format("#splitline{%g %s #in %s}{  in %g events}", th1->Integral(), formatHadronName(hadron).c_str(), legendEntries[i].c_str(), getNevts(inNames[i])).Data());
  }

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->Draw();
  for (int i = 0; i < histVector.size(); i++) {
    TH1D* hist = histVector[i];
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");

  latexText = "Tracking efficiency";
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex->Draw("same");

  // TLatex* integralsLatex0 = CreateLatex(0.25, 0.3, TString::Format("#splitline{%g %s #in ch. jets}{ %g events}", integrals[0], formatHadronName(hadron).c_str(), nevents[0]).Data(), textSize);
  // integralsLatex0->Draw("same");
  // TLatex* integralsLatex1 = CreateLatex(0.25, 0.2, TString::Format("#splitline{%g %s #in ch.+V0 jets}{ %g events}", integrals[1], formatHadronName(hadron).c_str(), nevents[1]).Data(), textSize);
  // integralsLatex1->Draw("same");

  saveName = "trackingEfficiencyComparison";
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

void compareMcV0ReconstructionEfficiency(vector<string> inNames, vector<string> legendEntries, string hadron)
{
  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double lowjetpt, highjetpt;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = false;
  double xMinFrame = 0., xMaxFrame = 10., yMinFrame = 0., yMaxFrame = 0.5;
  double xMinLegend = 0.25, xMaxLegend = 0.6, yMinLegend = 0.74, yMaxLegend = 0.89;
  double xLatex = 0.15, yLatex = 0.93;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 4;
  xTitle = "#it{p}_{T}";
  yTitle = "detector-level / particle-level";
  // yTitle = "#frac{1}{#it{N}_{evts}} #frac{d #it{N}}{d #it{p}_{T}}";
  drawoption = "same";

  string hadronForHistname = hadron;
  if ("Lambda0" == hadron) { hadronForHistname = "Lambda"; }
  if ("AntiLambda0" == hadron) { hadronForHistname = "AntiLambda"; }
  string histDir = "v0cascades-q-a";
  string histNameDet = TString::Format("%s/histos-V0/InvMass%sTrue", histDir.c_str(), hadronForHistname.c_str()).Data();
  // string histNamePart = TString::Format("%s/h3_particle_pt_particle_eta_particle_phi_associatedtrack_primary", histDir.c_str()).Data();
  string histNamePart = TString::Format("%s/histos-eve/GeneratedParticles", histDir.c_str()).Data();

  vector<double> integrals; vector<double> nevents;
  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  for (int i = 0; i < inNames.size(); i++) {
    TFile *inFile = TFile::Open(TString::Format("%s", inNames[i].c_str()).Data());
    TH3D* th3_det = (TH3D*)inFile->Get(histNameDet.c_str());
    TH3D* th3_part = (TH3D*)inFile->Get(histNamePart.c_str());
    int projectionBin = 1 + 2 * ("Lambda0" == hadron) + 4 * ("AntiLambda0" == hadron);
    th3_part->GetXaxis()->SetRange(projectionBin, projectionBin);

    TH1D* th1_det = th3_det->ProjectionX();
    TH1D* th1_part = th3_part->ProjectionY();

    th1_det->Rebin(rebinNumber);
    th1_part->Rebin(rebinNumber);

    TH1D* efficiency = (TH1D*)th1_det->Clone(TString::Format("efficiency_%s", inNames[i].c_str()).Data());
    efficiency->Divide(th1_part);

    setStyle(efficiency, i);
    legend->AddEntry(efficiency, TString::Format("%s: %g events", legendEntries[i].c_str(), getNevtsEfficiency(inNames[i])).Data());
    histVector.push_back(efficiency);
    integrals.push_back(efficiency->Integral());
    nevents.push_back(getNevtsEfficiency(inNames[i]));
    // legend->AddEntry(th1, TString::Format("#splitline{%g %s #in %s}{  in %g events}", th1->Integral(), formatHadronName(hadron).c_str(), legendEntries[i].c_str(), getNevts(inNames[i])).Data());
  }

  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->Draw();
  for (int i = 0; i < histVector.size(); i++) {
    TH1D* hist = histVector[i];
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");

  latexText = TString::Format("%s reconstruction efficiency", formatHadronName(hadron).c_str()).Data();
  TLatex* latex = CreateLatex(xLatex, yLatex, latexText, textSize);
  latex->Draw("same");

  // TLatex* integralsLatex0 = CreateLatex(0.25, 0.3, TString::Format("#splitline{%g %s #in ch. jets}{ %g events}", integrals[0], formatHadronName(hadron).c_str(), nevents[0]).Data(), textSize);
  // integralsLatex0->Draw("same");
  // TLatex* integralsLatex1 = CreateLatex(0.25, 0.2, TString::Format("#splitline{%g %s #in ch.+V0 jets}{ %g events}", integrals[1], formatHadronName(hadron).c_str(), nevents[1]).Data(), textSize);
  // integralsLatex1->Draw("same");

  saveName = TString::Format("%sReconstructionEfficiencyComparison", hadron.c_str()).Data();
  saveName = TString::Format("%s.pdf", saveName.c_str()).Data();
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

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
    compareDataPt(dataSet, inNames, legendEntries, doRatio);
  }
  else {
    compareDataZ(dataSet, inNames, legendEntries, hadron, jetptmin, jetptmax, doRatio, doCounts);
  }
}

void compareMassCuts(bool doRatio = false, string hadron = "", double jetptmin = 10., double jetptmax = 1e3, bool doCounts = false)
{
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }

  vector<string> inNames; vector<string> legendEntries;
  string dataSet = "LHC22o_pass6_minBias";

  string sK = "|#it{M}(K^{0}_{S}) - #it{m}_{K^{0}_{S}}| <";
  string sL = "|#it{M}(#Lambda) - #it{m}_{#Lambda}| <";
  string sAL = "|#it{M}(#bar{#Lambda}) - #it{m}_{#bar{#Lambda}}| <";

  // Train with tight mass acceptance
  string sCut = sK;
  double cutValue = 0.03;
  if ("Lambda0" == hadron) {
    sCut = sL;
    cutValue = 0.01;
  }
  if ("AntiLambda0" == hadron) {
    sCut = sAL;
    cutValue = 0.01;
  }

  // array<string, 2> file211540 = {"~/Documents/TrainOutput/211540/AnalysisResults.root", TString::Format("%s %.2f", sCut.c_str(), cutValue).Data()};
  // inNames.push_back(file211540[0]);
  // legendEntries.push_back(file211540[1]);
  array<string, 2> file225406 = {"~/Documents/TrainOutput/225406/AnalysisResults.root", TString::Format("%s %.2f", sCut.c_str(), cutValue).Data()};
  inNames.push_back(file225406[0]);
  legendEntries.push_back(file225406[1]);

  // Train with wide mass acceptance
  cutValue = 0.06;
  if ("Lambda0" == hadron) { cutValue = 0.02; }
  if ("AntiLambda0" == hadron) { cutValue = 0.02; }

  array<string, 2> file224894 = {"~/Documents/TrainOutput/224894/AnalysisResults.root", TString::Format("%s %.2f", sCut.c_str(), cutValue).Data()};
  inNames.push_back(file224894[0]);
  legendEntries.push_back(file224894[1]);

  compareDataZ(dataSet, inNames, legendEntries, hadron, jetptmin, jetptmax, doRatio, doCounts);
}

void compareZ_23d4_24b1(string hadron = "", double jetptmin = 10., double jetptmax = 1e3, bool doParticle = false, bool doRatio = false, bool doCounts = false)
{
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  vector<string> inNames; vector<string> legendEntries;
  string dataSet = "";

  array<string, 2> file210373 = {"~/Documents/TrainOutput/210373/AnalysisResults.root", "24b1"};
  inNames.push_back(file210373[0]);
  legendEntries.push_back(file210373[1]);

  array<string, 2> file210376 = {"~/Documents/TrainOutput/210376/AnalysisResults.root", "23d4"};
  inNames.push_back(file210376[0]);
  legendEntries.push_back(file210376[1]);

  compareMcZ(dataSet, inNames, legendEntries, hadron, jetptmin, jetptmax, doParticle, doRatio, doCounts);
}

void compareTrackingEfficiency_23d4_24b1()
{
  vector<string> inNames; vector<string> legendEntries;

  array<string, 2> file224905 = {"~/Documents/TrainOutput/224905/AnalysisResults.root", "24b1"};
  inNames.push_back(file224905[0]);
  legendEntries.push_back(file224905[1]);

  array<string, 2> file224906 = {"~/Documents/TrainOutput/224906/AnalysisResults.root", "23d4"};
  inNames.push_back(file224906[0]);
  legendEntries.push_back(file224906[1]);

  compareMcTrackingEfficiency(inNames, legendEntries);
}

void compareV0ReconstructionEfficiency_23d4_24b1(string hadron)
{
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  vector<string> inNames; vector<string> legendEntries;

  array<string, 2> file224905 = {"~/Documents/TrainOutput/224905/AnalysisResults.root", "24b1"};
  inNames.push_back(file224905[0]);
  legendEntries.push_back(file224905[1]);

  array<string, 2> file224906 = {"~/Documents/TrainOutput/224906/AnalysisResults.root", "23d4"};
  inNames.push_back(file224906[0]);
  legendEntries.push_back(file224906[1]);

  compareMcV0ReconstructionEfficiency(inNames, legendEntries, hadron);
}
