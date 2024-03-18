
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

#include "/Users/gijsvanweelden/Documents/Fragmentation/plotting_macros/histUtils.C"

void normaliseHistRowByRow(TH2F* hist);
void normaliseHistColByCol(TH2F* hist);

void plotOneHist(TCanvas* canvas, TH1F* frame, TH2F* hist, TLegend* legend, string saveName, string setDrawOption, string latexText);
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<TH1F*> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);
template <typename T>
T loadMatchedPtHist(string fileName, string dirName = "jet-fragmentation/matching/jets", string histName = "matchDetJetPtPartJetPt");
template <typename T>
void setStyle(T hist, int styleNumber);

// Compare the matched jet pt for a narrow or wide particle eta cut
void plotMatchedJetPt_ParticleEtaCutNarrowOrWide(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  double ptMin = 40, ptMax = 60;
  // double ptMin = 100, ptMax = 140;

  int normaliseHistograms = 1;
  bool doRatio = true;
  if (doRatio) { normaliseHistograms = 0; }

  string histName = "matchDetJetPtPartJetPt";
  string histTitle = "";
  string saveName = "matchedJetPt_projected_trackCutComparison";
  string xTitle = "#it{p}_{T}^{jet, det}";
  string yTitle = "normalised count";
  string legendTitle = "";
  string latexText = TString::Format("#splitline{Matching distance: 0.05}{#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  bool setLogZ = true;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 0, yMaxFrame = .3, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  int rebinNumber = 1;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  if (setLogZ) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  std::vector<string> nameVector;
  std::vector<string> legendVector;

  nameVector.push_back("./AnalysisResults-matchDist05-smallEta.root");
  legendVector.push_back("|#eta_{track}^{part}| < 0.9");
  nameVector.push_back("./AnalysisResults-matchDist05-largeEta.root");
  legendVector.push_back("|#eta_{track}^{part}| < 2.0");


  TH2F* matchedJetPt1 = loadMatchedPtHist<TH2F*>(inName1);
  matchedJetPt1->Sumw2();

  TH2F* matchedJetPt2 = loadMatchedPtHist<TH2F*>(inName2);
  matchedJetPt2->Sumw2();

  int firstBinPtTruth = 1, lastBinPtTruth = 999;
  for (unsigned int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string histLegend = legendVector[i];
    TH2F* hMatchedPt = loadMatchedPtHist<TH2F*>(fileName);
    firstBinPtTruth = hMatchedPt->GetYaxis()->FindBin(ptMin);
    lastBinPtTruth = hMatchedPt->GetYaxis()->FindBin(ptMax);
    TH1F* hist = (TH1F*)hMatchedPt->ProjectionX(firstBinPtTruth, lastBinPtTruth);
    hist->Rebin(rebinNumber);
    setStyle(hist, i);
    if (normaliseHistograms > 0) { hist->Scale(1./hist->Integral()); }
    else if (normaliseHistograms < 0) { hist->Scale(1./hist->GetBinContent(1)); }
    legend->AddEntry(hist, histLegend.c_str());
    histVector.push_back(hist);
  }
  if (doRatio) {
    TH1F* base = (TH1F*)histVector[0]->Clone("base");
    for (auto hist : histVector) {
      hist->Divide(base);
    }
    saveName = TString::Format("%s_ratio", saveName.c_str());
  }

  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  // saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

// Compare the matched jet pt for different matching distances
void plotMatchedJetPt_MatchingDistance(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName05 = "./AnalysisResults-dist05-noEtaCut.root";
  string inName10 = "./AnalysisResults-dist10-noEtaCut.root";
  string inName20 = "./AnalysisResults-dist20-noEtaCut.root";
  string saveDir = ".";
  TFile *inFile05 = TFile::Open(TString::Format("./%s", inName05.c_str()).Data());
  TFile *inFile10 = TFile::Open(TString::Format("./%s", inName10.c_str()).Data());
  TFile *inFile20 = TFile::Open(TString::Format("./%s", inName20.c_str()).Data());
  if(!inFile05){
    std::cout << "File " << inFile05 << " not found. Aborting program." << std::endl;
    return;
  }
  if(!inFile10){
    std::cout << "File " << inFile10 << " not found. Aborting program." << std::endl;
    return;
  }
  if(!inFile20){
    std::cout << "File " << inFile20 << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* dir05 = (TDirectory*)inFile05->Get("jet-fragmentation/matching/jets");
  TDirectory* dir10 = (TDirectory*)inFile10->Get("jet-fragmentation/matching/jets");
  TDirectory* dir20 = (TDirectory*)inFile20->Get("jet-fragmentation/matching/jets");

  // double ptMin = 40, ptMax = 60;
  double ptMin = 100, ptMax = 140;

  string histName = "matchDetJetPtPartJetPt";
  string histTitle = "";
  string saveName = "matchedJetPt_projected_matchDistComparison";
  string xTitle = "#it{p}_{T}^{jet, det}";
  string yTitle = "normalised count";
  string legendTitle = "";
  string latexText = TString::Format("#splitline{Cut after matching}{#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  bool setLogZ = true;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 0, yMaxFrame = .3, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  int rebinNumber = 1;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  if (setLogZ) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  TH2F* matchedJetPt05 = (TH2F*)dir05->Get(histName.c_str());
  matchedJetPt05->Sumw2();
  TH2F* matchedJetPt10 = (TH2F*)dir10->Get(histName.c_str());
  matchedJetPt10->Sumw2();
  TH2F* matchedJetPt20 = (TH2F*)dir20->Get(histName.c_str());
  matchedJetPt20->Sumw2();

  int firstBinPtTruth = 1, lastBinPtTruth = matchedJetPt05->GetNbinsY() + 1;
  firstBinPtTruth = matchedJetPt05->GetYaxis()->FindBin(ptMin);
  lastBinPtTruth = matchedJetPt05->GetYaxis()->FindBin(ptMax);

  TH1F* matchedJetPt05_projected = (TH1F*)matchedJetPt05->ProjectionX(TString::Format("matchedJetPt05_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  legend->AddEntry(matchedJetPt05_projected, "match dist 0.05");
  histVector.push_back(matchedJetPt05_projected);

  TH1F* matchedJetPt10_projected = (TH1F*)matchedJetPt10->ProjectionX(TString::Format("matchedJetPt10_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  legend->AddEntry(matchedJetPt10_projected, "match dist 0.1");
  histVector.push_back(matchedJetPt10_projected);

  TH1F* matchedJetPt20_projected = (TH1F*)matchedJetPt20->ProjectionX(TString::Format("matchedJetPt20_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  legend->AddEntry(matchedJetPt20_projected, "match dist 0.2");
  histVector.push_back(matchedJetPt20_projected);

  for (int i = 0; i < histVector.size(); i++) {
    auto hist = histVector[i];
    hist->SetLineWidth(3);
    hist->SetMarkerSize(2);
    hist->SetLineColor(GetColor(i));
    hist->SetMarkerColor(GetColor(i));
    hist->SetMarkerStyle(GetMarker(i));
    hist->Scale(1./hist->Integral());
  }

  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}
void plotMatchedJetPt_MatchingDistance_Ratio(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName05 = "./AnalysisResults-dist05-noEtaCut.root";
  string inName10 = "./AnalysisResults-dist10-noEtaCut.root";
  string inName20 = "./AnalysisResults-dist20-noEtaCut.root";
  string saveDir = ".";
  TFile *inFile05 = TFile::Open(TString::Format("./%s", inName05.c_str()).Data());
  TFile *inFile10 = TFile::Open(TString::Format("./%s", inName10.c_str()).Data());
  TFile *inFile20 = TFile::Open(TString::Format("./%s", inName20.c_str()).Data());
  if(!inFile05){
    std::cout << "File " << inFile05 << " not found. Aborting program." << std::endl;
    return;
  }
  if(!inFile10){
    std::cout << "File " << inFile10 << " not found. Aborting program." << std::endl;
    return;
  }
  if(!inFile20){
    std::cout << "File " << inFile20 << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* dir05 = (TDirectory*)inFile05->Get("jet-fragmentation/matching/jets");
  TDirectory* dir10 = (TDirectory*)inFile10->Get("jet-fragmentation/matching/jets");
  TDirectory* dir20 = (TDirectory*)inFile20->Get("jet-fragmentation/matching/jets");

  // double ptMin = 40, ptMax = 60;
  double ptMin = 100, ptMax = 140;
  int baseId = 2;
  int base = 5 * (baseId == 0) + 10 * (baseId == 1) + 20 * (baseId == 2);

  string histName = "matchDetJetPtPartJetPt";
  string histTitle = "";
  string saveName = "matchedJetPt_projected_matchDistRatio";
  string xTitle = "#it{p}_{T}^{jet, det}";
  string yTitle = "ratio";
  string legendTitle = "";
  string latexText = TString::Format("#splitline{Cut after matching}{#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 0, yMaxFrame = 5, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  int rebinNumber = 1;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  TH2F* matchedJetPt05 = (TH2F*)dir05->Get(histName.c_str());
  matchedJetPt05->Sumw2();
  TH2F* matchedJetPt10 = (TH2F*)dir10->Get(histName.c_str());
  matchedJetPt10->Sumw2();
  TH2F* matchedJetPt20 = (TH2F*)dir20->Get(histName.c_str());
  matchedJetPt20->Sumw2();

  int firstBinPtTruth = 1, lastBinPtTruth = matchedJetPt05->GetNbinsY() + 1;
  firstBinPtTruth = matchedJetPt05->GetYaxis()->FindBin(ptMin);
  lastBinPtTruth = matchedJetPt05->GetYaxis()->FindBin(ptMax);

  TH1F* matchedJetPt05_projected = (TH1F*)matchedJetPt05->ProjectionX(TString::Format("matchedJetPt05_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  legend->AddEntry(matchedJetPt05_projected, "match dist 0.05");
  histVector.push_back(matchedJetPt05_projected);

  TH1F* matchedJetPt10_projected = (TH1F*)matchedJetPt10->ProjectionX(TString::Format("matchedJetPt10_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  legend->AddEntry(matchedJetPt10_projected, "match dist 0.1");
  histVector.push_back(matchedJetPt10_projected);

  TH1F* matchedJetPt20_projected = (TH1F*)matchedJetPt20->ProjectionX(TString::Format("matchedJetPt20_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  legend->AddEntry(matchedJetPt20_projected, "match dist 0.2");
  histVector.push_back(matchedJetPt20_projected);

  TH1F* denominator = (TH1F*)histVector[baseId]->Clone("denominator");

  for (int i = 0; i < histVector.size(); i++) {
    auto hist = histVector[i];
    hist->SetLineWidth(3);
    hist->SetMarkerSize(2);
    hist->SetLineColor(GetColor(i));
    hist->SetMarkerColor(GetColor(i));
    hist->SetMarkerStyle(GetMarker(i));
    hist->Divide(denominator);
  }

  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s_base%02d", saveName.c_str(), base);
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

// Compare the jet pt when the jet eta cut is applied before or after matching
void plotMatchedJetPt_JetEtaCutBeforeOrAfterMatching(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName1 = "./AnalysisResults-dist40-wEtaCut.root";
  string inName2 = "./AnalysisResults-dist40-noEtaCut.root";
  string saveDir = ".";
  TFile *inFile1 = TFile::Open(TString::Format("./%s", inName1.c_str()).Data());
  TFile *inFile2 = TFile::Open(TString::Format("./%s", inName2.c_str()).Data());
  if(!inFile1){
    std::cout << "File " << inFile1 << " not found. Aborting program." << std::endl;
    return;
  }
  if(!inFile2){
    std::cout << "File " << inFile2 << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* dir1 = (TDirectory*)inFile1->Get("jet-fragmentation/matching/jets");
  TDirectory* dir2 = (TDirectory*)inFile2->Get("jet-fragmentation/matching/jets");

  double ptMin = 40, ptMax = 60;
  // double ptMin = 100, ptMax = 140;

  string histName = "matchDetJetPtPartJetPt";
  string histTitle = "";
  string saveName = "matchedJetPt_projected";
  string xTitle = "#it{p}_{T}^{jet, det}";
  string yTitle = "normalised count";
  string legendTitle = "";
  string latexText = TString::Format("#splitline{Matching distance: 0.4}{#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  bool setLogZ = true;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 0, yMaxFrame = .3, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  int rebinNumber = 1;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  if (setLogZ) { myCanvas->SetLogz(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  TH2F* matchedJetPt1 = (TH2F*)dir1->Get(histName.c_str());
  matchedJetPt1->Sumw2();
  TH2F* matchedJetPt2 = (TH2F*)dir2->Get(histName.c_str());
  matchedJetPt2->Sumw2();
  int firstBinPtTruth = 1, lastBinPtTruth = matchedJetPt1->GetNbinsY() + 1;
  firstBinPtTruth = matchedJetPt1->GetYaxis()->FindBin(ptMin);
  lastBinPtTruth = matchedJetPt1->GetYaxis()->FindBin(ptMax);

  TH1F* matchedJetPt1_projected = (TH1F*)matchedJetPt1->ProjectionX(TString::Format("matchedJetPt1_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  matchedJetPt1_projected->SetLineWidth(3);
  matchedJetPt1_projected->SetLineColor(GetColor(0));
  matchedJetPt1_projected->SetMarkerStyle(GetMarker(0));
  matchedJetPt1_projected->SetMarkerColor(GetColor(0));
  TH1F* ratioHist = (TH1F*)matchedJetPt1->Clone("ratio");

  TH1F* matchedJetPt2_projected = (TH1F*)matchedJetPt2->ProjectionX(TString::Format("matchedJetPt2_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  matchedJetPt2_projected->SetLineWidth(3);
  matchedJetPt2_projected->SetLineColor(GetColor(1));
  matchedJetPt2_projected->SetMarkerStyle(GetMarker(1));
  matchedJetPt2_projected->SetMarkerColor(GetColor(1));

  matchedJetPt1_projected->Scale(1./matchedJetPt1_projected->Integral());
  matchedJetPt2_projected->Scale(1./matchedJetPt2_projected->Integral());
  legend->AddEntry(matchedJetPt1_projected, "#eta cut before matching");
  histVector.push_back(matchedJetPt1_projected);
  legend->AddEntry(matchedJetPt2_projected, "#eta cut after matching");
  histVector.push_back(matchedJetPt2_projected);

  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}
void plotMatchedJetPt_JetEtaCutBeforeOrAfterMatching_Ratio(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName1 = "./AnalysisResults-dist20-wEtaCut.root";
  string inName2 = "./AnalysisResults-dist20-noEtaCut.root";
  string saveDir = ".";
  TFile *inFile1 = TFile::Open(TString::Format("./%s", inName1.c_str()).Data());
  TFile *inFile2 = TFile::Open(TString::Format("./%s", inName2.c_str()).Data());
  if(!inFile1){
    std::cout << "File " << inFile1 << " not found. Aborting program." << std::endl;
    return;
  }
  if(!inFile2){
    std::cout << "File " << inFile2 << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* dir1 = (TDirectory*)inFile1->Get("jet-fragmentation/matching/jets");
  TDirectory* dir2 = (TDirectory*)inFile2->Get("jet-fragmentation/matching/jets");

  double ptMin = 40, ptMax = 60;
  // double ptMin = 100, ptMax = 140;
  double xMinFrame = 0, xMaxFrame = 100, yMinFrame = 0.8, yMaxFrame = 1.2, zMinFrame = 1e-5, zMaxFrame = 1;

  string histName = "matchDetJetPtPartJetPt";
  string histTitle = "";
  string saveName = "matchedJetPt_projected";
  string xTitle = "#it{p}_{T}^{jet, det}";
  string yTitle = "normalised count";
  string legendTitle = "";
  string latexText = TString::Format("#splitline{Matching distance: 0.2}{#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  bool setLogZ = false;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int R = 4;
  int rebinNumber = 1;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  TH2F* matchedJetPt1 = (TH2F*)dir1->Get(histName.c_str());
  matchedJetPt1->Sumw2();
  TH2F* matchedJetPt2 = (TH2F*)dir2->Get(histName.c_str());
  matchedJetPt2->Sumw2();
  int firstBinPtTruth = 1, lastBinPtTruth = matchedJetPt1->GetNbinsY() + 1;

  firstBinPtTruth = matchedJetPt1->GetYaxis()->FindBin(ptMin);
  lastBinPtTruth = matchedJetPt1->GetYaxis()->FindBin(ptMax);
  TH1F* matchedJetPt1_projected = (TH1F*)matchedJetPt1->ProjectionX(TString::Format("matchedJetPt1_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  matchedJetPt1_projected->SetLineWidth(3);
  matchedJetPt1_projected->SetLineColor(GetColor(0));
  matchedJetPt1_projected->SetMarkerStyle(GetMarker(0));
  matchedJetPt1_projected->SetMarkerColor(GetColor(0));
  TH1F* ratioHist = (TH1F*)matchedJetPt1_projected->Clone("ratio");

  firstBinPtTruth = matchedJetPt2->GetYaxis()->FindBin(ptMin);
  lastBinPtTruth = matchedJetPt2->GetYaxis()->FindBin(ptMax);
  TH1F* matchedJetPt2_projected = (TH1F*)matchedJetPt2->ProjectionX(TString::Format("matchedJetPt2_projected_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);

  ratioHist->Divide(matchedJetPt2_projected);
  frame->GetYaxis()->SetTitle("#frac{cut before matching}{cut after matching}");
  histVector.push_back(ratioHist);

  saveName = TString::Format("%s_ratio", saveName.c_str()).Data();
  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s/%s.pdf", saveDir.c_str(), saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "hist p", latexText);

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
void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<TH1F*> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText)
{
  canvas->cd();
  frame->Draw();
  string drawOption = "same";
  if (setDrawOption != "") {
    drawOption = TString::Format("%s %s", drawOption.c_str(), setDrawOption.c_str()).Data();
  }
  for (unsigned int i = 0; i < histVector.size(); i++) {
    TH1F* hist = histVector[i];
    hist->Draw(drawOption.c_str());
  }
  if (legend) { legend->Draw("same"); }
  if (latexText != "") { DrawLatex(0.3, 0.8, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

template <typename T>
T loadMatchedPtHist(string fileName, string dirName, string histName)
{
  TFile *inFile = TFile::Open(TString::Format("./%s", fileName.c_str()).Data());
  TDirectory* dir = (TDirectory*)inFile->Get(TString::Format("%s", dirName.c_str()).Data());
  T hist = (T)dir->Get(histName.c_str());
  // hist->Sumw2();
  return hist;
}
template <typename T>
void setStyle(T hist, int styleNumber)
{
  hist->SetLineWidth(3);
  hist->SetLineColor(GetColor(styleNumber));
  hist->SetMarkerStyle(GetMarker(styleNumber));
  hist->SetMarkerColor(GetColor(styleNumber));
}