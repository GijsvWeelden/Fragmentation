
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

void plotFragsO2_old(void)
{
  double time = clock();
  gROOT->SetBatch();
  gStyle->SetNdivisions(505);
  // string inName = "../data/LHC21k6/test109192.root";
  // string saveDir = "../Plots/LHC21k6/test109192";
  string inName = "../data/LHC21k6/train107435.root";
  string saveDir = "../Plots/LHC21k6/train107435";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* fragDir = (TDirectory*)inFile->Get("jet-fragmentation");
  TDirectory* jetsDir = (TDirectory*)fragDir->Get("jets");
  // TDirectory* matchJetsDir = (TDirectory*)matchDir->Get("jets");
  // TDirectory* matchTracksDir = (TDirectory*)matchDir->Get("tracks");

  // comparePtBins(jetsDir, kData, kPt, saveDir);
  // comparePtBins(jetsDir, kData, kProj, saveDir);
  // compareZDefinitions(jetsDir, kData, 20, 40, saveDir);
  // compareZDefinitions(jetsDir, kData, 40, 60, saveDir);
  // compareZDefinitions(jetsDir, kData, 60, 80, saveDir);
  // plotJetQA(jetsDir, kData, -999, -999, saveDir);
  // plotJetQA(jetsDir, kData, 20, 40, saveDir);
  // plotJetQA(jetsDir, kData, 40, 60, saveDir);
  // plotJetQA(jetsDir, kData, 60, 80, saveDir);
  // plotTrackQA(tracksDir, kData, -999, -999, saveDir);

  // Matching plots
  // plotJetMatchingQA(matchJetsDir, -999, -999, saveDir);
  plotJetMatchingQA(jetsDir, 40, 60, saveDir);
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

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

// Comparison of jet pt bins
void comparePtBins(TDirectory* inDir, int dataOrMC, int projOrPt, string saveDir)
{
  std::vector<TH1F*> histVector; std::vector<string> nameVector;
  string xTitle = "xTitle", yTitle = "#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}}";
  string histTitle = "Jet fragmentation", saveName = "saveName", legendTitle = "";
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 1e-1, yMaxFrame = 1e2;
  double xMinLegend = 0.5, xMaxLegend = 0.75, yMinLegend = 0.6, yMaxLegend = 0.9;
  double ptMin = -1, ptMax = -1;
  bool setLogY = true;

  switch (projOrPt) {
    case kProj:
      xTitle = "#it{p}^{proj} / #it{p}^{jet}";
      saveName = "ptProj";
      break;
    case kPt:
      xTitle = "#it{p}_{T}^{track} / #it{p}_{T}^{jet}";
      saveName = "ptRatio";
      break;
    default:
      std::cout << "comparePtBins: invalid input. projOrPt = " << projOrPt << std::endl
        << "Allowed inputs: projOrPt = " << kProj << ", " << kPt << std::endl;
      return;
  }

  // Pt bins to compare
  // ptMin = 5, ptMax = 20;
  // TH1F* frag_5_20 = loadFragmentation(inDir, dataOrMC, projOrPt, ptMin, ptMax);
  // histVector.push_back(frag_5_20);
  // nameVector.push_back("5-20 GeV");
  // saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);

  ptMin = 20, ptMax = 40;
  TH1F* frag_20_40 = loadFragmentation(inDir, dataOrMC, projOrPt, ptMin, ptMax);
  histVector.push_back(frag_20_40);
  nameVector.push_back("20-40 GeV");
  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);

  ptMin = 40, ptMax = 60;
  TH1F* frag_40_60 = loadFragmentation(inDir, dataOrMC, projOrPt, ptMin, ptMax);
  histVector.push_back(frag_40_60);
  nameVector.push_back("40-60 GeV");
  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);

  ptMin = 60, ptMax = 80;
  TH1F* frag_60_80 = loadFragmentation(inDir, dataOrMC, projOrPt, ptMin, ptMax);
  histVector.push_back(frag_60_80);
  nameVector.push_back("60-80 GeV");
  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);

  // ptMin = 80, ptMax = 100;
  // TH1F* frag_80_100 = loadFragmentation(inDir, dataOrMC, projOrPt, ptMin, ptMax);
  // histVector.push_back(frag_80_100);
  // nameVector.push_back("80-100 GeV");
  // saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);

  // ptMin = 100, ptMax = 120;
  // TH1F* frag_100_120 = loadFragmentation(inDir, dataOrMC, projOrPt, ptMin, ptMax);
  // histVector.push_back(frag_100_120);
  // nameVector.push_back("100-120 GeV");
  // saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);

  // ptMin = 120, ptMax = 140;
  // TH1F* frag_120_140 = loadFragmentation(inDir, dataOrMC, projOrPt, ptMin, ptMax);
  // histVector.push_back(frag_120_140);
  // nameVector.push_back("120-140 GeV");
  // saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);

  // ptMin = 140, ptMax = 160;
  // TH1F* frag_140_160 = loadFragmentation(inDir, dataOrMC, projOrPt, ptMin, ptMax);
  // histVector.push_back(frag_140_160);
  // nameVector.push_back("140-160 GeV");
  // saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);

  // ptMin = 160, ptMax = 200;
  // TH1F* frag_160_200 = loadFragmentation(inDir, dataOrMC, projOrPt, ptMin, ptMax);
  // histVector.push_back(frag_160_200);
  // nameVector.push_back("160-200 GeV");
  // saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);

  saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str());
  plotNHists(histVector, nameVector,
             xTitle, yTitle, histTitle, legendTitle, saveName,
             xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
             xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
             setLogY, "");
}
// Comparison of z definitions
void compareZDefinitions(TDirectory* inDir, int dataOrMC, double ptMin, double ptMax, string saveDir)
{
  std::vector<TH1F*> histVector; std::vector<string> histNameVector;
  string xTitle = "#it{z}", yTitle = "#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}}";
  string histTitle = "Jet fragmentation: #it{p}^{proj} / #it{p}^{jet} vs. #it{p}_{T}^{track} / #it{p}_{T}^{jet}";
  string legendTitle = TString::Format("Jet #it{p}_{T}: %.0f-%.0f GeV", ptMin, ptMax).Data();
  string saveName = TString::Format("zComparison_pt%.0f-%.0f", ptMin, ptMax).Data();
  string setHistDrawOption = "", setFuncDrawOption = "";

  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 1e-1, yMaxFrame = 90;
  double xMinLegend = 0.5, xMaxLegend = 0.8, yMinLegend = 0.6, yMaxLegend = 0.9;
  bool setLogY = true;

  TH1F* frag = loadFragmentation(inDir, dataOrMC, kPt, ptMin, ptMax);
  histVector.push_back(frag);
  histNameVector.push_back("#it{p}_{T}^{track} / #it{p}_{T}^{jet}");
  TH1F* ptProj = loadFragmentation(inDir, dataOrMC, kProj, ptMin, ptMax);
  histVector.push_back(ptProj);
  histNameVector.push_back("#it{p}^{proj} / #it{p}^{jet}");

  saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str());
  plotNHists(histVector, histNameVector,
             xTitle, yTitle, histTitle, legendTitle, saveName,
             xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
             xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
             setLogY, setHistDrawOption);

  // Make ratio plot
  saveName = TString::Format("zRatio_pt%.0f-%.0f", ptMin, ptMax).Data();
  histTitle = TString::Format("Jet fragmentation ratio: #it{z}_{proj} / #it{z}_{#it{p}_{T}} (%.0f-%.0f GeV)", ptMin, ptMax).Data();
  xTitle = "#it{z}_{proj}", yTitle = "N(#it{z}_{proj}) / N(#it{z}_{#it{p}_{T}})";
  xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0.8, yMaxFrame = 1.2;
  setLogY = false;

  TH1F* zRatio = (TH1F*)ptProj->Clone("zRatio");
  zRatio->Divide(frag);
  saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str());
  plotOneHist(zRatio, xTitle, yTitle, histTitle, legendTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, false);

  // Make ratio of difference plot
  saveName = TString::Format("zDiffRatio_pt%.0f-%.0f", ptMin, ptMax).Data();
  histTitle = TString::Format("Jet fragmentation ratio: (#it{z}_{proj} - #it{z}_{#it{p}_{T}}) / #it{z}_{proj} (%.0f-%.0f GeV)", ptMin, ptMax).Data();
  xTitle = "#it{z}_{proj}", yTitle = "#frac{N(#it{z}_{proj}) - N(#it{z}_{#it{p}_{T}})}{N(#it{z}_{#it{p}_{T}})}";
  xMinFrame = 0, xMaxFrame = 1, yMinFrame = -1, yMaxFrame = 1;

  TH1F* zDiffRatio = (TH1F*)ptProj->Clone("zDiffRatio");
  zDiffRatio->Add(frag, -1);
  zDiffRatio->Divide(ptProj);
  saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str());
  plotOneHist(zDiffRatio, xTitle, yTitle, histTitle, legendTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, false);
}
// Normalise 2D histogram row-by-row
void normaliseHistRowByRow(TH2F* hist)
{
  int firstColBin = 1, lastColBin = hist->GetNbinsX()+1;
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY()+1;
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
  int firstColBin = 1, lastColBin = hist->GetNbinsX()+1;
  int firstRowBin = 1, lastRowBin = hist->GetNbinsY()+1;
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
// Jet QA plots
void plotJetQA(TDirectory* inDir, int dataOrMC, double ptMin, double ptMax, string saveDir)
{
  string xTitle = "#it{p}_{T}^{jet}", yTitle = "#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{p}_{T}}";
  string histTitle; string legendTitle; string saveName;
  double xMinFrame = 0, xMaxFrame = 100, yMinFrame = 1e-7, yMaxFrame = 1e0;
  double xMinLegend = 0.6, xMaxLegend = 0.8, yMinLegend = 0.6, yMaxLegend = 0.8;
  bool setLogY = true, drawLegend = false;

  string projectionAxis = "X", dataType = "";
  if ( dataOrMC == kData ) { dataType = "jet"; }
  else if ( dataOrMC == kMCDetector ) { dataType = "detJet"; }
  else if ( dataOrMC == kMCParticle ) { dataType = "partJet"; }
  string histName = TString::Format("%sPtEtaPhi", dataType.c_str()).Data();
  // TH3F* h3Jets_1 = (TH3F*)inDir->Get(histName.c_str());
  TH3F* h3Jets_2 = (TH3F*)inDir->Get(histName.c_str());

  int firstXBin = 1;
  if ( ptMin > -900 ) { h3Jets_2->GetXaxis()->FindBin(ptMin); }
  int lastXBin = h3Jets_2->GetNbinsX() + 1;
  if ( ptMax > -900 ) { h3Jets_2->GetXaxis()->FindBin(ptMax); }
  int firstYBin = 1;
  int lastYBin = h3Jets_2->GetNbinsY() + 1;
  int firstZBin = 1;
  int lastZBin = h3Jets_2->GetNbinsZ() + 1;

  // Pt spectrum
  TH1F* hJetPt = (TH1F*)h3Jets_2->ProjectionX("hJetPt", firstYBin, lastYBin, firstZBin, lastZBin);
  histTitle = "Jet #it{p}_{T} spectrum";
  hJetPt->Scale(1./hJetPt->Integral(), "width");
  saveName = TString::Format("%s/jetPtSpectrum", saveDir.c_str());
  plotOneHist(hJetPt, xTitle, yTitle, histTitle, legendTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, drawLegend);

  // Eta-Phi plot
  xTitle = "#eta^{jet}", yTitle = "#phi^{jet}";
  histTitle = TString::Format("Jet #eta, #phi (%.0f - %.0f GeV)", ptMin, ptMax);
  xMinFrame = -1, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2*TMath::Pi();
  setLogY = false;
  h3Jets_2->GetXaxis()->SetRangeUser(ptMin, ptMax); // SetRangeUser is risky!
  TH2F* hJetEtaPhi = (TH2F*)h3Jets_2->Project3D("zy");
  hJetEtaPhi->Scale(1./hJetEtaPhi->Integral());
  hJetEtaPhi->SetName(TString::Format("etaphi%.0f_%.0f", ptMin, ptMax));
  saveName = TString::Format("%s/jetEtaPhi_pt%.0f-%.0f", saveDir.c_str(), ptMin, ptMax);
  plotOneHist(hJetEtaPhi, xTitle, yTitle, histTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, drawLegend);
}
// Jet matching QA plots
void plotJetMatchingQA(TDirectory* inDir, double ptMinTruth, double ptMaxTruth, string saveDir)
{
  std::vector<TH1F*> histVector; std::vector<string> nameVector;
  string xTitle, yTitle;
  double xMinFrame, xMaxFrame, yMinFrame, yMaxFrame;
  double ptMin = 0, ptMax = 200, etaMin = -1, etaMax = 1, phiMin = 0, phiMax = 2*TMath::Pi();
  int binPtMin, binPtMax;
  string histTitle; string legendTitle; string saveName; string histName;
  double xMinLegend = 0.6, xMaxLegend = 0.8, yMinLegend = 0.6, yMaxLegend = 0.8;
  bool setLogY = false, setLogZ = true, drawLegend = false;

  histName = "matchDetJetPtPartJetPt";
  TH2F* matchDetJetPtPartJetPt = (TH2F*)inDir->Get(histName.c_str());
  // Projected version
  xTitle = "#it{p}_{T}^{jet, det}", yTitle = "normalised count";
  histTitle = "Matched #it{p}_{T}^{jet, det}";
  string projectedHistName = "matchDetJetPtPartJetPt_projected";
  saveName = TString::Format("%s/matchedJetPt_projected", saveDir.c_str());
  binPtMin = 1, binPtMax = matchDetJetPtPartJetPt->GetNbinsX() + 1;
  xMinFrame = 0, xMaxFrame = 200, yMinFrame = 0, yMaxFrame = 0.3;
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    double min = ptMinTruth * (ptMinTruth > 0), max = 200;
    if (ptMaxTruth > 0 && ptMaxTruth < 200) { max = ptMaxTruth; }
    histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
  //   // yMinFrame = min; yMaxFrame = max;
    xMaxFrame = max * 1.5;
  //   projectedHistName = TString::Format("matchDetJetPtPartJetPt_projected_ptTruth%.0f-%.0f", ptMinTruth, ptMaxTruth).Data();
  }
  TH1F* matchDetJetPtPartJetPt_projected =
        (TH1F*)matchDetJetPtPartJetPt->ProjectionX(projectedHistName.c_str(), 9, 13);
  matchDetJetPtPartJetPt_projected->Scale(1./matchDetJetPtPartJetPt_projected->Integral());
  plotOneHist(matchDetJetPtPartJetPt_projected, xTitle, yTitle, histTitle, "", saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, drawLegend);

  // 2D version
  normaliseHistRowByRow(matchDetJetPtPartJetPt);
  matchDetJetPtPartJetPt->GetZaxis()->SetRangeUser(1e-5, 1);
  xTitle = "#it{p}_{T}^{jet, det}", yTitle = "#it{p}_{T}^{jet, part}";
  histTitle = "Matched #it{p}_{T} (normalised row-by-row)";
  saveName = TString::Format("%s/matchedJetPt", saveDir.c_str());
  xMinFrame = 0, xMaxFrame = 200, yMinFrame = 0, yMaxFrame = 200;
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    double min = ptMinTruth * (ptMinTruth > 0), max = 200;
    if (ptMaxTruth > 0 && ptMaxTruth < 200) { max = ptMaxTruth; }
    histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
    yMinFrame = min; yMaxFrame = max;
  }
  plotOneHist(matchDetJetPtPartJetPt, xTitle, yTitle, histTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogZ, drawLegend);

  // Matched eta
  // histName = "matchPartJetPtDetJetEtaPartJetEta";
  // TH3F* matchPartJetPtDetJetEtaPartJetEta = (TH3F*)inDir->Get(histName.c_str());
  // if (ptMinTruth > -900 || ptMaxTruth > -900) {
  //   double min = ptMinTruth * (ptMinTruth > 0), max = 200;
  //   if (ptMaxTruth > 0 && ptMaxTruth < 200) { max = ptMaxTruth; }
  //   histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
  //   saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
  //   matchPartJetPtDetJetEtaPartJetEta->GetXaxis()->SetRangeUser(min, max); // SetRange is better, but is bugged for project3D
  // }
  histName = "matchDetJetEtaPartJetEta";
  TH2F* matchDetJetEtaPartJetEta = (TH2F*)inDir->Get(histName.c_str());
  // 2D version
  normaliseHistRowByRow(matchDetJetEtaPartJetEta);
  matchDetJetEtaPartJetEta->GetZaxis()->SetRangeUser(1e-5, 1);
  xTitle = "#eta^{jet, det}", yTitle = "#eta^{jet, part}";
  xMinFrame = etaMin, xMaxFrame = etaMax, yMinFrame = etaMin, yMaxFrame = etaMax;
  histTitle = "Matched #eta (normalised row-by-row)";
  saveName = TString::Format("%s/matchedJetEta", saveDir.c_str());
  plotOneHist(matchDetJetEtaPartJetEta, xTitle, yTitle, histTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogZ, drawLegend);

  // histName = "matchPartJetPtDetJetPhiPartJetPhi";
  histName = "matchDetJetPhiPartJetPhi";
  xTitle = "#phi^{jet, det}", yTitle = "#phi^{jet, part}";
  xMinFrame = phiMin, xMaxFrame = phiMax, yMinFrame = phiMin, yMaxFrame = phiMax;
  // TH3F* matchPartJetPtDetJetPhiPartJetPhi = (TH3F*)inDir->Get(histName.c_str());
  TH2F* matchDetJetPhiPartJetPhi = (TH2F*)inDir->Get(histName.c_str());
  // if (ptMinTruth > -900 || ptMaxTruth > -900) {
  //   double min = ptMinTruth * (ptMinTruth > 0), max = 200;
  //   if (ptMaxTruth > 0 && ptMaxTruth < 200) { max = ptMaxTruth; }
  //   histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
  //   saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
  //   matchPartJetPtDetJetPhiPartJetPhi->GetXaxis()->SetRangeUser(min, max); // SetRange is better, but is bugged for project3D
  // }
  // TH2F* matchDetJetPhiPartJetPhi = (TH2F*)matchPartJetPtDetJetPhiPartJetPhi->Project3D("zy");
  normaliseHistRowByRow(matchDetJetPhiPartJetPhi);
  matchDetJetPhiPartJetPhi->GetZaxis()->SetRangeUser(1e-5, 1);
  histTitle = "Matched #phi (normalised row-by-row)";
  saveName = TString::Format("%s/matchedJetPhi", saveDir.c_str());
  plotOneHist(matchDetJetPhiPartJetPhi, xTitle, yTitle, histTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogZ, drawLegend);

  histName = "matchPartJetPtResolutionPt";
  TH2F* matchPartJetPtResolutionPt = (TH2F*)inDir->Get(histName.c_str());
  xTitle = "#it{p}_{T}^{jet, part} - #it{p}_{T}^{jet, det}", yTitle = "";
  xMinFrame = -2, xMaxFrame = 2, yMinFrame = 0, yMaxFrame = 2;
  histTitle = "Matched #it{p}_{T} resolution";
  saveName = TString::Format("%s/matchedJetResolutionPt", saveDir.c_str());
  binPtMin = 1, binPtMax = matchPartJetPtResolutionPt->GetNbinsX() + 1;
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    double min = ptMinTruth * (ptMinTruth > 0), max = 200;
    if (ptMaxTruth > 0 && ptMaxTruth < 200) { max = ptMaxTruth; }
    histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), min, max);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), min, max);
    binPtMin = matchPartJetPtResolutionPt->GetXaxis()->FindBin(min);
    binPtMax = matchPartJetPtResolutionPt->GetXaxis()->FindBin(max);
  }
  TH1F* jetResolutionPt = (TH1F*)matchPartJetPtResolutionPt->ProjectionY("resolutionPt", binPtMin, binPtMax);
  jetResolutionPt->Scale(1./jetResolutionPt->Integral());
  plotOneHist(jetResolutionPt, xTitle, yTitle, histTitle, legendTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, drawLegend);

  // FIXME: eta phi QA plots
  histName = "matchPartJetPtResolutionEta";
  TH3F* jetPtResolutionEta = (TH3F*)inDir->Get(histName.c_str());
  xTitle = "#eta^{jet, part}", yTitle = "(#eta^{jet, part} - #eta^{jet, det})";
  xMinFrame = etaMin, xMaxFrame = etaMax, yMinFrame = -5, yMaxFrame = 5;
  histTitle = "Matched #eta resolution (normalised column-by-column)";
  saveName = TString::Format("%s/jetResolutionEta", saveDir.c_str());
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    double min = ptMinTruth * (ptMinTruth > 0), max = 200;
    if (ptMaxTruth > 0 && ptMaxTruth < 200) { max = ptMaxTruth; }
    histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
    jetPtResolutionEta->GetXaxis()->SetRangeUser(min, max); // SetRange is better, but is bugged for project3D
  }
  TH2F* jetResolutionEta = (TH2F*)jetPtResolutionEta->Project3D("zy");
  normaliseHistColByCol(jetResolutionEta);
  // matchDetJetPhiPartJetPhi->GetZaxis()->SetRangeUser(1e-5, 1);
  plotOneHist(jetResolutionEta, xTitle, yTitle, histTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogZ, drawLegend);

  histName = "matchPartJetPtResolutionPhi";
  TH3F* jetPtResolutionPhi = (TH3F*)inDir->Get(histName.c_str());
  xTitle = "#phi^{jet, part}", yTitle = "(#phi^{jet, part} - #phi^{jet, det})";
  xMinFrame = phiMin, xMaxFrame = phiMax, yMinFrame = -5, yMaxFrame = 5;
  histTitle = "Matched #phi resolution (normalised column-by-column)";
  saveName = TString::Format("%s/jetResolutionPhi", saveDir.c_str());
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    double min = ptMinTruth * (ptMinTruth > 0), max = 200;
    if (ptMaxTruth > 0 && ptMaxTruth < 200) { max = ptMaxTruth; }
    histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
    jetPtResolutionPhi->GetXaxis()->SetRangeUser(min, max); // SetRange is better, but is bugged for project3D
  }
  TH2F* jetResolutionPhi = (TH2F*)jetPtResolutionPhi->Project3D("zy");
  normaliseHistColByCol(jetResolutionPhi);
  plotOneHist(jetResolutionPhi, xTitle, yTitle, histTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogZ, drawLegend);

  histName = "matchPartJetPtMatchDist";
  TH2F* jetPtMatchDist = (TH2F*)inDir->Get(histName.c_str());
  xTitle = "#Delta^{jet}", yTitle = "normalised count";
  xMinFrame = 0, xMaxFrame = 0.2, yMinFrame = 1e-2, yMaxFrame = 1;
  histTitle = "Distance between matched jets";
  saveName = TString::Format("%s/jetMatchDist", saveDir.c_str());
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    double min = ptMinTruth * (ptMinTruth > 0), max = 200;
    if (ptMaxTruth > 0 && ptMaxTruth < 200) { max = ptMaxTruth; }
    histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
  }
  TH1F* jetMatchDist = projectHist(jetPtMatchDist, "y", "matchDist", ptMinTruth, ptMaxTruth);
  jetMatchDist->Scale(1./jetMatchDist->Integral());
  legendTitle = "";
  setLogY = true;
  plotOneHist(jetMatchDist, xTitle, yTitle, histTitle, legendTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, drawLegend);

  histName = "matchPartJetPtEnergyScale";
  TH2F* jetPtEnergyScale = (TH2F*)inDir->Get(histName.c_str());
  // Projected version
  xTitle = "#it{p}_{T}^{ jet, det} / #it{p}_{T}^{ jet, truth}", yTitle = "normalised count";
  histTitle = "Jet Energy Scale";
  xMinFrame = 0, xMaxFrame = 1.6, yMinFrame = 1e-5, yMaxFrame = .8;
  saveName = TString::Format("%s/jetEnergyScale_projected", saveDir.c_str());
  // binPtMin = 1, binPtMax = jetPtEnergyScale->GetNbinsX() + 1;
  binPtMin = 9, binPtMax = 13;
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    double min = ptMinTruth * (ptMinTruth > 0), max = 200;
    if (ptMaxTruth > 0 && ptMaxTruth < 200) { max = ptMaxTruth; }
    histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
    binPtMin = jetPtEnergyScale->GetXaxis()->FindBin(min);
    binPtMax = jetPtEnergyScale->GetXaxis()->FindBin(max);
  }
  TH1F* jetPtEnergyScale_projected = (TH1F*)jetPtEnergyScale->ProjectionY("jetEnergyScaleProjected", binPtMin, binPtMax);
  jetPtEnergyScale_projected->Scale(1./jetPtEnergyScale_projected->Integral());
  setLogY = true;
  plotOneHist(jetPtEnergyScale_projected, xTitle, yTitle, histTitle, legendTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, drawLegend);

  // 2D version
  xTitle = "#it{p}_{T}^{ jet, truth}", yTitle = "#it{p}_{T}^{ jet, det} / #it{p}_{T}^{ jet, truth}";
  xMinFrame = ptMinTruth, xMaxFrame = ptMaxTruth, yMinFrame = 0, yMaxFrame = 2;
  normaliseHistColByCol(jetPtEnergyScale);
  jetPtEnergyScale->GetZaxis()->SetRangeUser(1e-5, 1);
  histTitle = "Jet energy scale (normalised column-by-column)";
  saveName = TString::Format("%s/jetEnergyScale", saveDir.c_str());
  setLogZ = true;
  plotOneHist(jetPtEnergyScale, xTitle, yTitle, histTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogZ, drawLegend);

  // TODO:
  //  * fakes
  //  * misses
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
  xTitle = "#it{z}_{proj}^{det}", yTitle = "#it{z}_{proj}^{part}";
  xMinFrame = zMin, xMaxFrame = zMax, yMinFrame = zMin, yMaxFrame = zMax;
  histTitle = "Matched #it{z}_{proj}";
  saveName = TString::Format("%s/ResponseZproj", saveDir.c_str());
  if (ptMinTruth > -900 || ptMaxTruth > -900) {
    histTitle = TString::Format("%s (#it{p}_{T}^{ jet, truth}: %.0f - %.0f)", histTitle.c_str(), ptMinTruth, ptMaxTruth);
    saveName = TString::Format("%s_ptTruth%.0f-%.0f", saveName.c_str(), ptMinTruth, ptMaxTruth);
  }

  histName = "matchDetJetPtTrackProjPartJetPtTrackProj";
  THnF* responseTrackProj = (THnF*)inDir->Get(histName.c_str());
  TH2F* zTzD_TrackProj = projectHist(responseTrackProj, 1, 3, "trackProj", -999, -999, -999, -999, ptMinTruth, ptMaxTruth, -999, -999);
  normaliseHistRowByRow(zTzD_TrackProj);
  zTzD_TrackProj->GetZaxis()->SetRangeUser(1e-5, 1);
  histTitle = TString::Format("%s (normalised row-by-row)", histTitle.c_str());
  plotOneHist(zTzD_TrackProj, xTitle, yTitle, histTitle, TString::Format("%s", saveName.c_str()).Data(),
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogZ, drawLegend);

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
}
// Track QA plots
void plotTrackQA(TDirectory* inDir, int dataOrMC, double ptMin, double ptMax, string saveDir)
{
  string xTitle = "#it{p}_{T}^{track}", yTitle = "#frac{1}{#it{N}_{tracks}} #frac{d #it{N}}{d #it{p}_{T}}";
  string histTitle; string legendTitle; string saveName;
  double xMinFrame = 0, xMaxFrame = 100, yMinFrame = 1e-7, yMaxFrame = 1e0;
  double xMinLegend = 0.6, xMaxLegend = 0.8, yMinLegend = 0.6, yMaxLegend = 0.8;
  bool setLogY = true, drawLegend = false;

  string projectionAxis = "X", dataType = "";
  if ( dataOrMC == kData ) { dataType = "track"; }
  else if ( dataOrMC == kMCDetector ) { dataType = "detTrack"; }
  else if ( dataOrMC == kMCParticle ) { dataType = "partTrack"; }
  string histName = TString::Format("%sPtEtaPhi", dataType.c_str()).Data();

  TH3F* h3Tracks = (TH3F*)inDir->Get(histName.c_str());

  int firstXBin = 1;
  if ( ptMin > -900 ) { h3Tracks->GetXaxis()->FindBin(ptMin); }
  int lastXBin = h3Tracks->GetNbinsX() + 1;
  if ( ptMax > -900 ) { h3Tracks->GetXaxis()->FindBin(ptMax); }
  int firstYBin = 1;
  int lastYBin = h3Tracks->GetNbinsY() + 1;
  int firstZBin = 1;
  int lastZBin = h3Tracks->GetNbinsZ() + 1;

  // Pt spectrum
  TH1F* hTrackPt = (TH1F*)h3Tracks->ProjectionX("hTrackPt", firstYBin, lastYBin, firstZBin, lastZBin);
  histTitle = "Track #it{p}_{T} spectrum";
  hTrackPt->Scale(1./hTrackPt->Integral(), "width");
  saveName = TString::Format("%s/trackPtSpectrum", saveDir.c_str());
  plotOneHist(hTrackPt, xTitle, yTitle, histTitle, legendTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, drawLegend);

  // Eta-Phi plot
  xTitle = "#eta^{track}", yTitle = "#phi^{track}";
  histTitle = TString::Format("Track #eta, #phi (%.0f - %.0f GeV)", ptMin, ptMax);
  xMinFrame = -1, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2*TMath::Pi();
  setLogY = false;

  h3Tracks->GetXaxis()->SetRangeUser(ptMin, ptMax); // This doesn't work with SetRange()
  TH2F* hTrackEtaPhi = (TH2F*)h3Tracks->Project3D("zy");
  hTrackEtaPhi->Scale(1./hTrackEtaPhi->Integral());
  hTrackEtaPhi->SetName(TString::Format("etaphi%.0f_%.0f", ptMin, ptMax));
  saveName = TString::Format("%s/trackEtaPhi_pt%.0f-%.0f", saveDir.c_str(), ptMin, ptMax);
  plotOneHist(hTrackEtaPhi, xTitle, yTitle, histTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, drawLegend);
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
