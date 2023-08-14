
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

enum { mQuark = 0,
       mGluon,
       mInclusive,
       mQuarkGluon,
       mQuarkGluonInclusive
       };
enum { kData = 0,
       kMCParticle,
       kMCDetector
       };
enum { kProj = 0, // Definition of z
       kPt
       };

void plotOneHist(TCanvas* canvas, TH1F* frame, TH2F* hist, TLegend* legend, string saveName, string setDrawOption, string latexText);

void comparePtBins(TDirectory* inDir, int dataOrMC, int projOrPt, string saveDir);
void compareZDefinitions(TDirectory* inDir, int dataOrMC, double ptMin, double ptMax, string saveDir);

string formatHadronName(string hadron);

int getNJets(TFile* inFile, double ptMin, double ptMax, int quarkOrGluon);

TH1F* loadFragmentation(TDirectory* inDir, int dataOrMC, int projOrPt, double ptMin, double ptMax);
void normaliseHistRowByRow(TH2F* hist);
void normaliseHistColByCol(TH2F* hist);

void plotJetMatchingQA(TDirectory* inDir, double ptMinTruth, double ptMaxTruth, string saveDir);
void plotJetQA(TDirectory* inDir, int dataOrMC, double ptMin, double ptMax, string saveDir);
void plotResponse(TDirectory* inDir, double ptMinTruth, double ptMaxTruth, string saveDir);
void plotTrackQA(TDirectory* inDir, int dataOrMC, double ptMin, double ptMax, string saveDir);

// TH1F* projectHist(TH2F* inputHist, string projectionAxis, string histName, double axMin, double axMax);

// TODO: Add theory (selection from various sets; give a fragmentation scale)

void plotResponseChargeFrag_old(void)
{
  double time = clock();
  // gROOT->SetBatch();
  gStyle->SetNdivisions(505);
  string inName = "../data/LHC21k6/train107435.root";
  string saveDir = "../Plots/LHC21k6/train107435";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* fragDir = (TDirectory*)inFile->Get("jet-fragmentation");
  TDirectory* matchJetsDir = (TDirectory*)fragDir->Get("jets");
  // TDirectory* matchDir = (TDirectory*)fragDir->Get("matching");
  // TDirectory* matchJetsDir = (TDirectory*)matchDir->Get("jets");
  // TDirectory* matchTracksDir = (TDirectory*)matchDir->Get("tracks");

  // Matching plots
  // plotResponse(matchJetsDir, -999, -999, saveDir);
  // plotResponse(matchJetsDir, 5, 20, saveDir);
  // plotResponse(matchJetsDir, 20, 30, saveDir);
  // plotResponse(matchJetsDir, 30, 40, saveDir);
  // plotResponse(matchJetsDir, 40, 60, saveDir);
  // plotResponse(matchJetsDir, 60, 80, saveDir);
  // plotResponse(matchJetsDir, 80, 100, saveDir);
  // plotResponse(matchJetsDir, 100, 120, saveDir);
  // plotResponse(matchJetsDir, 120, 160, saveDir);
  // plotResponse(matchJetsDir, 160, 200, saveDir);

  // return;

  string histName = "matchDetJetPtTrackProjPartJetPtTrackProj";
  string histTitle = "";
  string saveName = "responseZproj";
  string xTitle = "#it{z}^{det}";
  string yTitle = "#it{z}^{part}";
  string legendTitle = "";
  string latexText = "";
    // "#splitline{ALICE 2022 data}{#splitline{13.6 TeV pp 500 kHz}{#splitline{anti-kt full jets}{#splitline{Trigger: #it{R} = 0.2, #it{p}_{T}^{jet} > 8 GeV}{ }}}}";
  double textSize = 0.03;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  bool setLogZ = true;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 1, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  double ptMin = 40, ptMax = 60;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  if (setLogZ) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  THnF* response = (THnF*)matchJetsDir->Get(histName.c_str());
  int firstBinPtTruth = 1, lastBinPtTruth = response->GetAxis(3)->GetNbins();
  int projectionAxisX = 1, projectionAxisY = 3;
  int ptTruthAxis = 2;

  firstBinPtTruth = response->GetAxis(ptTruthAxis)->FindBin(ptMin);
  lastBinPtTruth = response->GetAxis(ptTruthAxis)->FindBin(ptMax);

  response->GetAxis(ptTruthAxis)->SetRange(firstBinPtTruth, lastBinPtTruth);
  TH2F* zTzD = (TH2F*)response->Projection(projectionAxisY, projectionAxisX);
  // TH2F* zTzD = projectHist(response, 1, 3, "trackProj", -999, -999, -999, -999, ptMin, ptMax, -999, -999);
  // zTzD->Rebin2D(5, 5);
  normaliseHistRowByRow(zTzD);
  zTzD->GetZaxis()->SetRangeUser(zMinFrame, zMaxFrame);
  // zTzD->Draw("colz");
  // return;
  saveName = TString::Format("%s_ptTruth%.0f-%.0f_old.pdf", saveName.c_str(), ptMin, ptMax).Data();
  latexText = TString::Format("#splitline{PYTHIA, ideal alignment}{#splitline{13.6 TeV pp 500 kHz}{#splitline{anti-kt jets, #it{R} = 0.%d}{#it{p}_{T}^{jet, truth}: %.0f-%.0f GeV}}}", R, ptMin, ptMax).Data();
  plotOneHist(myCanvas, frame, zTzD, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
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
// Response matrix plots
void plotResponse(TDirectory* inDir, double ptMinTruth, double ptMaxTruth, string saveDir)
{
  string xTitle, yTitle;
  double xMinFrame, xMaxFrame, yMinFrame, yMaxFrame;
  double ptMin = 0, ptMax = 200, zMin = 0, zMax = 1;
  int ptMinTruthBin = 1, ptMaxTruthBin = -1;
  string histTitle; string saveName; string histName;
  double xMinLegend = 0.6, xMaxLegend = 0.8, yMinLegend = 0.6, yMaxLegend = 0.8;
  bool setLogZ = true, drawLegend = false;

  // Track projection
  xTitle = "#it{z}^{det}", yTitle = "#it{z}^{part}";
  xMinFrame = zMin, xMaxFrame = zMax, yMinFrame = zMin, yMaxFrame = zMax;
  // histTitle = "Matched #it{z}";
  histTitle = "";
  // saveName = TString::Format("%s/ResponseZproj", saveDir.c_str());
  saveName = "responseZproj";
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    // histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
  }

  histName = "matchDetJetPtTrackProjPartJetPtTrackProj";
  THnF* responseTrackProj = (THnF*)inDir->Get(histName.c_str());
  TH2F* zTzD_TrackProj = projectHist(responseTrackProj, 1, 3, "trackProj", -999, -999, -999, -999, ptMinTruth, ptMaxTruth, -999, -999);
  normaliseHistRowByRow(zTzD_TrackProj);
  zTzD_TrackProj->GetZaxis()->SetRangeUser(1e-5, 1);
  // histTitle = TString::Format("%s (normalised row-by-row)", histTitle.c_str());
  plotOneHist(zTzD_TrackProj, xTitle, yTitle, histTitle, TString::Format("%s", saveName.c_str()).Data(),
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogZ, drawLegend);

  /*
  // Frag: pT ratio
  xTitle = "#it{z}_{#it{p}_{T}}^{det}", yTitle = "#it{z}_{#it{p}_{T}}^{part}";
  xMinFrame = zMin, xMaxFrame = zMax, yMinFrame = zMin, yMaxFrame = zMax;
  histTitle = "Matched #it{z}_{#it{p}_{T}}";
  saveName = TString::Format("%s/ResponseZpt", saveDir.c_str());
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
  }

  histName = "matchDetJetPtFragPartJetPtFrag";
  THnF* responseFrag = (THnF*)inDir->Get(histName.c_str());
  TH2F* zTzD_Frag = projectHist(responseFrag, 1, 3, "trackProj", -999, -999, -999, -999, ptMinTruth, ptMaxTruth, -999, -999);
  normaliseHistRowByRow(zTzD_Frag);
  zTzD_Frag->GetZaxis()->SetRangeUser(1e-5, 1);
  histTitle = TString::Format("%s (normalised row-by-row)", histTitle.c_str());
  plotOneHist(zTzD_Frag, xTitle, yTitle, histTitle, TString::Format("%s", saveName.c_str()).Data(),
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogZ, drawLegend);
  // */
}

// Get number of jets in pt range
int getNJets(TDirectory* inDir, double ptMin, double ptMax, int dataOrMC)
{
  string projectionAxis = "X", dataType = "";
  if ( dataOrMC == kData ) { dataType = "jet"; }
  else if ( dataOrMC == kMCDetector ) { dataType = "detJet"; }
  else if ( dataOrMC == kMCParticle ) { dataType = "partJet"; }
  string histName = TString::Format("%sPtEtaPhi", dataType.c_str()).Data();

  TH3F* h3nJets = (TH3F*)inDir->Get(histName.c_str());
  // TH1F* hNJets = projectHist(hNJetTypes, projectionAxis, "nJets", ptMin, ptMax);
  int firstXBin = h3nJets->GetXaxis()->FindBin(ptMin);
  int lastXBin = h3nJets->GetXaxis()->FindBin(ptMax);
  int firstZBin = 1;
  int lastZBin = h3nJets->GetNbinsZ() + 1;
  TH1F* h1nJets = (TH1F*)h3nJets->ProjectionY("nJetsProjected", firstXBin, lastXBin, firstZBin, lastZBin);

  int nJets = h1nJets->Integral();
  // cout << dataType << " pt " << ptMin << "-" << ptMax << " nJets: " << nJets << endl;

  return nJets;
}

// Returns the fragmentation histogram into a hadron for a given pt range
TH1F* loadFragmentation(TDirectory* inDir, int dataOrMC, int projOrPt, double ptMin, double ptMax)
{
  int nJets = -1;
  string projectionAxis = "Y", histName = "";
  string zDefinition = "";
  if ( projOrPt == kProj ) { zDefinition = "TrackProj"; }
  else if ( projOrPt == kPt ) { zDefinition = "Frag"; }
  else {
    cout << "loadFragmentation: invalid input. projOrPt = " << projOrPt << std::endl
      << "Allowed inputs: projOrPt = " << kProj << ", " << kPt << std::endl;
    return nullptr;
  }

  switch (dataOrMC) {
    case kData:
      histName = TString::Format("jetPt%s", zDefinition.c_str()).Data();
      break;
    case kMCParticle:
      histName = TString::Format("partJetPt%s", zDefinition.c_str()).Data();
      break;
    case kMCDetector:
      histName = TString::Format("detJetPt%s", zDefinition.c_str()).Data();
      break;
    default:
      std::cout << "loadFragmentation: invalid input. dataOrMC = " << dataOrMC << std::endl
        << "Allowed inputs: dataOrMC = " << kData << ", " << kMCParticle << ", " << kMCDetector << std::endl;
      return nullptr;
  }

  TH2F* jets = loadHist<TH2F*>(inDir, histName);
  TH1F* fragmentation = projectHist(jets, projectionAxis, TString::Format("%s_projected_pt%.0f-%.0f", histName.c_str(), ptMin, ptMax).Data(), ptMin, ptMax);
  fragmentation->Rebin(5);
  nJets = getNJets(inDir, ptMin, ptMax, dataOrMC); // Get number of jets of appropriate type
  fragmentation->Scale(1./nJets, "width"); // Normalise by number of jets and bin width
  return fragmentation;
}

// Formats the hadron name to look nice (Greek letters, sub- and superscripts)
string formatHadronName(string hadron)
{
  string had = hadron;
  if (hadron == "pi"){
    had = "#pi^{#pm}";
  }
  else if (hadron == "pi0"){
    had = "#pi^{0}";
  }
  else if (hadron == "K0L"){
    had = "K^{0}_{L}";
  }
  else if (hadron == "K0S"){
    had = "K^{0}_{S}";
  }
  else if (hadron == "K0"){
    had = "K^{0}";
  }
  else if (hadron == "K"){
    had = "K^{#pm}";
  }
  else if (hadron == "Lambda0"){
    had = "#Lambda^{0}";
  }
  return had;
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
  if (latexText != "") { DrawLatex(0.3, 0.9, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
