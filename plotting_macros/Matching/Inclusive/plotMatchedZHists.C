
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
T loadHist(string fileName, string histName);
template <typename T>
T loadHist(string fileName, string dirName, string histName);

template <typename T>
void setStyle(T hist, int styleNumber);

void plotMatchedJetPt(double ptMin = 40, double ptMax = 60, bool doRatio = false)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  int normaliseHistograms = 1;
  if (doRatio) { normaliseHistograms = 0; }

  string dirName  = "jet-fragmentation/matching/jets";
  string histName = "matchDetJetPtPartJetPt";
  string histTitle = "";
  string saveName = "matchedJetPt_projected";
  string xTitle = "#it{p}_{T}^{jet, det}";
  string yTitle = "weighted count";
  if (normaliseHistograms > 0) { yTitle = "normalised count"; }
  if (doRatio) { yTitle = "ratio"; }
  string legendTitle = "";
  string latexText = TString::Format("#splitline{Max matching distance: 0.05}{#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  bool setLogY = false;
  double xMinFrame = 0, xMaxFrame = 100, yMinFrame = 0, yMaxFrame = .3, zMinFrame = 1e-5, zMaxFrame = 1;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  int rebinNumber = 1;

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  std::vector<string> nameVector;
  std::vector<string> legendVector;

  nameVector.push_back("./AnalysisResults-matchDist05-smallEta.root");
  legendVector.push_back("|#eta_{track}^{part}| < 0.9");
  // nameVector.push_back("./AnalysisResults-matchDist05-largeEta.root");
  // legendVector.push_back("|#eta_{track}^{part}| < 2.0");

  int firstBinPtTruth = 1, lastBinPtTruth = 999;
  for (unsigned int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string legendEntry = legendVector[i];

    TH2F* h2D = loadHist<TH2F*>(fileName, dirName, histName);
    h2D->Sumw2();
    firstBinPtTruth = h2D->GetXaxis()->FindBin(ptMin + 1e-3);
    lastBinPtTruth = h2D->GetXaxis()->FindBin(ptMax - 1e-3);
    TH1F* hist = (TH1F*)h2D->ProjectionY(TString::Format("jetEnergyScale_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
    hist->Rebin(rebinNumber);
    setStyle(hist, i);
    if (normaliseHistograms > 0) { hist->Scale(1./hist->Integral()); }
    else if (normaliseHistograms < 0) { hist->Scale(1./hist->GetBinContent(1)); }
    if (nameVector.size() > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
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
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void plotMatchedJetPt2D(double ptMin = 40, double ptMax = 60, bool doRatio = false)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string fileName = "./AnalysisResults-matchDist05-largeEta.root";
  string dirName  = "jet-fragmentation/matching/jets";
  string histName = "matchDetJetPtPartJetPt";
  string histTitle = "";
  string saveName = "matchedJetPt";
  string xTitle = "#it{p}_{T}^{jet, det}";
  string yTitle = "#it{p}_{T}^{jet, truth}";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  bool setLogY = false;
  bool setLogZ = true;
  double xMinFrame = 0, xMaxFrame = 200, yMinFrame = 0, yMaxFrame = 200, zMinFrame = 1e-5, zMaxFrame = 1;
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
  TH2F* matchedJetPt = loadHist<TH2F*>(fileName, dirName, histName);
  matchedJetPt->Sumw2();
  matchedJetPt->Rebin2D(rebinNumber, rebinNumber);
  normaliseHistRowByRow(matchedJetPt);
  matchedJetPt->GetZaxis()->SetRangeUser(zMinFrame, zMaxFrame);

  // saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotOneHist(myCanvas, frame, matchedJetPt, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void plotMatchedJetEnergyScale(double ptMin = 40, double ptMax = 60, bool doRatio = false)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  int normaliseHistograms = 1;
  if (doRatio) { normaliseHistograms = 0; }

  string dirName  = "jet-fragmentation/matching/jets";
  string histName = "matchPartJetPtEnergyScale";
  string histTitle = "";
  string saveName = "matchedJetPtResolution";
  string xTitle = "#it{p}_{T}^{jet, det} / #it{p}_{T}^{jet, truth}";
  string yTitle = "weighted count";
  if (normaliseHistograms > 0) { yTitle = "normalised count"; }
  if (doRatio) { yTitle = "ratio"; }
  string legendTitle = "";
  string latexText = TString::Format("#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;
  int rebinNumber = 1;

  // Plotting stuff
  bool setLogY = false;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = 0, xMaxFrame = 10, yMinFrame = 0, yMaxFrame = 1;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  std::vector<string> nameVector;
  std::vector<string> legendVector;

  nameVector.push_back("./AnalysisResults-sample1.root");
  legendVector.push_back("");
  // nameVector.push_back("./AnalysisResults-matchDist05-smallEta.root");
  // legendVector.push_back("|#eta_{track}^{part}| < 0.9");
  // nameVector.push_back("./AnalysisResults-matchDist05-largeEta.root");
  // legendVector.push_back("|#eta_{track}^{part}| < 2.0");
  // latexText = TString::Format("#splitline{Max matching distance: 0.05}{%s}", latexText.c_str()).Data();

  int firstBinPtTruth = 1, lastBinPtTruth = 999;
  for (unsigned int i = 0; i < nameVector.size(); i++) {
    string fileName = nameVector[i];
    string legendEntry = legendVector[i];
    TH2F* h2D = loadHist<TH2F*>(fileName, dirName, histName);
    h2D->Sumw2();
    firstBinPtTruth = h2D->GetXaxis()->FindBin(ptMin);
    lastBinPtTruth = h2D->GetXaxis()->FindBin(ptMax);
    TH1F* hist = (TH1F*)h2D->ProjectionY(TString::Format("jetEnergyScale_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
    hist->Rebin(rebinNumber);
    setStyle(hist, i);
    if (normaliseHistograms > 0) { hist->Scale(1./hist->Integral()); }
    else if (normaliseHistograms < 0) { hist->Scale(1./hist->GetBinContent(1)); }
    if (nameVector.size() > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
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
  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void plotResolutionTrackProj(double ptMin = 40, double ptMax = 60, bool doRatio = false)
{
  double time = clock();
  gStyle->SetNdivisions(505);

  int normaliseHistograms = 1;
  if (doRatio) { normaliseHistograms = 0; }

  string dirName  = "";
  string histName = "jet-fragmentation/matching/jets/matchPartJetPtResolutionTrackProj";
  string histTitle = "";
  string saveName = "matchedJetResolutionTrackProj";
  string xTitle = "#it{z}_{ch}^{part} - #it{z}_{ch}^{det}";
  string yTitle = "";
  string legendTitle = "";
  string latexText = TString::Format("#it{p}_{T, jet}^{part}: %.0f-%.0f GeV/#it{c}", ptMin, ptMax).Data();
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.03;//0.05;
  int R = 4;
  int rebinNumber = 2;

  // Plotting stuff
  bool setLogY = true;
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  double xMinFrame = -1, xMaxFrame = 1, yMinFrame = 1e-9, yMaxFrame = 1;
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Histogram stuff
  std::vector<TH1F*> histVector;
  std::vector<string> nameVector;
  std::vector<string> legendVector;

  nameVector.push_back("./AnalysisResults-sample1.root");
  legendVector.push_back("");

  TH3F* h3 = loadHist<TH3F*>("AnalysisResults-sample1.root", histName);
  // TH1F* hist = new TH1F("dZ", "dZ; #it{z}_{ch}^{truth} - #it{z}_{ch}^{det}",
  //                       h3->GetNbinsZ(), h3->GetZaxis()->GetXmin(), h3->GetZaxis()->GetZmax()
  //                      );
  int firstBinPtTruth = 1, lastBinPtTruth = h3->GetNbinsX();
  int firstBinZ = 1, lastBinZ = h3->GetNbinsY();
  firstBinPtTruth = h3->GetXaxis()->FindBin(ptMin);
  lastBinPtTruth  = h3->GetXaxis()->FindBin(ptMax);
  TH1F* hist = (TH1F*)h3->ProjectionZ("dz", firstBinPtTruth, lastBinPtTruth, firstBinZ, lastBinZ);
  histVector.push_back(hist);

  // for (unsigned int i = 0; i < nameVector.size(); i++) {
  //   string fileName = nameVector[i];
  //   string legendEntry = legendVector[i];
  //   TH3F* h3D = loadHist<TH3F*>(fileName, dirName, histName);
  //   h3D->Sumw2();
  //   firstBinPtTruth = h3D->GetXaxis()->FindBin(ptMin);
  //   lastBinPtTruth = h3D->GetXaxis()->FindBin(ptMax);
  //   TH1F* hist = (TH1F*)h3D->ProjectionZ(TString::Format("deltaPt_pt%.0f-%.0f", ptMin, ptMax).Data(), firstBinPtTruth, lastBinPtTruth);
  //   hist->Rebin(rebinNumber);
  //   setStyle(hist, i);
  //   if (normaliseHistograms > 0) { hist->Scale(1./hist->Integral()); }
  //   else if (normaliseHistograms < 0) { hist->Scale(1./hist->GetBinContent(1)); }
  //   if (nameVector.size() > 1) { legend->AddEntry(hist, legendEntry.c_str()); }
  //   histVector.push_back(hist);
  // }
  // if (doRatio) {
  //   TH1F* base = (TH1F*)histVector[0]->Clone("base");
  //   for (auto hist : histVector) {
  //     hist->Divide(base);
  //   }
  //   saveName = TString::Format("%s_ratio", saveName.c_str());
  // }

  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  // saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("./%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void plotPtJetLevelVsTrackLevel(string inName,
                                double ptMin = 40, double ptMax = 60,
                                int rebinNumber = 1,
                                int compareTruth = 1, bool doRatio = false,
                                string drawOption = "")
{
  double time = clock();
  gStyle->SetNdivisions(505);

  string histTitle = "";
  string saveName = "matchedJetPtRmPt";
  string xTitle = "#it{p}_{T}^{jet}";
  string yTitle = "count/N_{jet-pairs}";
  string legendTitle = "";
  string latexText = "";
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;//0.05;
  int R = 4;
  // int rebinNumber = 1;

  int normaliseHistograms = 1;
  if (doRatio) {
    normaliseHistograms = 0;
    yTitle = "Ratio";
    saveName = TString::Format("%s_ratio", saveName.c_str());
  }

  // Plotting stuff
  bool setLogY = true;
  if (doRatio) { setLogY = false; }
  double xMinFrame = 0, xMaxFrame = 300, yMinFrame = 1e-8, yMaxFrame = 10;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.7, yMaxLegend = 0.8;

  string dirName    = "jet-fragmentation/matching/jets";
  string jetsName   = TString::Format("%s/matchDetJetPtPartJetPt", dirName.c_str()).Data();
  string tracksName = TString::Format("%s/matchDetJetPtTrackProjPartJetPtTrackProj", dirName.c_str()).Data();

  TH2F* hJets = loadHist<TH2F*>(inName, jetsName);
  hJets->Sumw2();
  THnSparseF* RM = loadHist<THnSparseF*>(inName, tracksName);
  RM->Sumw2();
  TH2F* hRM = (TH2F*)RM->Projection(2,0);

  TH1F* detJets   = (TH1F*)hJets->ProjectionX();
  TH1F* truthJets = (TH1F*)hJets->ProjectionY();
  TH1F* detRM     = (TH1F*)hRM->ProjectionX();
  TH1F* truthRM   = (TH1F*)hRM->ProjectionY();

  // Histogram stuff
  std::vector<TH1F*>  histVector; std::vector<TH1F*>  tmpVector;
  std::vector<TH1F*>  detHists      = { detJets, detRM };
  std::vector<TH1F*>  truthHists    = { truthJets, truthRM };
  std::vector<string> legendEntries = { "Matched jets", "From RM" };

  if (compareTruth > 0) {
    tmpVector = truthHists;
    xTitle = "#it{p}_{T, jet}^{part}";
    saveName = TString::Format("%s_truth", saveName.c_str());
  }
  else {
    tmpVector = detHists;
    xTitle = "#it{p}_{T, jet}^{det}";
    saveName = TString::Format("%s_det", saveName.c_str());
  }

  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 900, 900);
  if (setLogY) { myCanvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  double integral = tmpVector[0]->Integral();
  for (unsigned int i = 0; i < tmpVector.size(); i++) {
    string legendEntry = legendEntries[i];
    TH1F* hist = tmpVector[i];
    int firstBinPt = 1, lastBinPt = hist->GetNbinsX();
    firstBinPt = hist->GetXaxis()->FindBin(ptMin + 0.1);
    lastBinPt = hist->GetXaxis()->FindBin(ptMax + 0.1);
    if (ptMax > hist->GetXaxis()->GetBinLowEdge(lastBinPt)+ 0.1) {
      lastBinPt++;
    }
    hist->Rebin(rebinNumber);
    setStyle(hist, i);
    if (normaliseHistograms > 0) { hist->Scale(1./integral); }
    // if (normaliseHistograms > 0) { hist->Scale(1./hist->Integral()); }
    // else if (normaliseHistograms < 0) { hist->Scale(1./hist->GetBinContent(1)); }
    legend->AddEntry(hist, legendEntry.c_str());
    histVector.push_back(hist);
  }
  if (doRatio) {
    TH1F* base = (TH1F*)tmpVector[0]->Clone("base");
    for (auto hist : histVector) {
      hist->Divide(base);
    }
  }

  saveName = TString::Format("%s_pt%.0f-%.0f", saveName.c_str(), ptMin, ptMax);
  // saveName = TString::Format("%s_binsize%.d", saveName.c_str(), rebinNumber);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(myCanvas, frame, histVector, legend, saveName, drawOption, latexText);

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
  if (latexText != "") { DrawLatex(0.3, 0.9, latexText.c_str(), legend->GetTextSize()); }
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}
template <typename T>
T loadHist(string fileName, string histName)
{
  TFile *inFile = TFile::Open(TString::Format("./%s", fileName.c_str()).Data());
  T hist = (T)inFile->Get(histName.c_str());
  return hist;
}
template <typename T>
T loadHist(string fileName, string dirName, string histName)
{
  TFile *inFile = TFile::Open(TString::Format("./%s", fileName.c_str()).Data());
  TDirectory* dir = (TDirectory*)inFile->Get(TString::Format("%s", dirName.c_str()).Data());
  T hist = (T)dir->Get(histName.c_str());
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
