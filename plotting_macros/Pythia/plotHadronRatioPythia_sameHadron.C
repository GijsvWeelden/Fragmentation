
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

#include "../histUtils.C"

enum { mQuark = 0,
       mGluon,
       mInclusive,
       mQuarkGluon,
       mQuarkGluonInclusive
       };

string formatHadronName(string hadron);

int getNJets(TFile* inFile, double ptMin, double ptMax, int quarkOrGluon);

TF1* loadKKP();
TH1F* loadFragmentation(TFile* inFile, string hadron, int quarkOrGluon, double ptMin, double ptMax);
TF1* loadTheory(string hadron, int quarkOrGluon, double QSquared);

TH1F* makeRatio(TFile* inFile, string nominatorHadron, int nominatorType, string denominatorHadron, int denominatorType, double ptMin, double ptMax);

void plotHadronRatios(TFile* inFile, std::vector<string> nominators, std::vector<string> denominators,
                      std::vector<int> nominatorTypes, std::vector<int> denominatorTypes,
                      string theory = "",
                      double ptMin = -999, double ptMax = -999);
void plotMultipleHadrons(TFile* inFile, std::vector<string> hadrons, int quarkOrGluon, string theory = "",
                         double ptMin = -999, double ptMax = -999);
void plotSingleHadron(TFile* inFile, string hadron, int quarkOrGluon, string theory = "",
                      double ptMin = -999, double ptMax = -999);

void plotGluonQuarkJetsFraction(TFile* inFile);

void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<TH1F*> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText);

void plotHadronRatioPythia_sameHadron(void)
{
  double time = clock();
  gStyle->SetNdivisions(505);
  string inName = "PythiaResult_pthat20-80";
  TFile *inFile = TFile::Open(TString::Format("../../data/pythia/%s.root", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  std::vector<double> jetPtBins = {40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};
	std::vector<string> partons = {"q", "g"};

  int R = 4;
  double jetR = R * 1e-1;
  // string saveDir = "Plots/67485/JetSubstructure";
  string saveDir = ".";
  double ptMin = 40, ptMax = 60;
  string saveName = "";
  string hadron = "";

  std::vector<TH1F*> histVector;
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0, yMaxFrame = 2.5;
  double xMinLegend = 0.6, xMaxLegend = 0.7, yMinLegend = 0.45, yMaxLegend = 0.55;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.05;
  int lineWidth = 3; int markerSize = 1;
  bool setLogY = false;
  string setDrawOption = "";
  string xTitle = "#it{z}";
  string yTitle = "#frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{z}} (gluon) / #frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{z}} (quark)";
  string histTitle = "";
  string legendTitle = "";
  string latexText =
    TString::Format("#splitline{PYTHIA Tune 4c}{#splitline{13.6 TeV pp}{#splitline{anti-kt jets, #it{R} = 0.%d}{#it{p}_{T}^{ jet}: %.0f - %.0f GeV/#it{c}}}}", R, ptMin, ptMax).Data();
  // string latexText =
  //   TString::Format("#splitline{PYTHIA Tune 4c, 13.6 TeV pp}{anti-kt jets, #it{R} = 0.%d, #it{p}_{T}^{ jet}: %.0f - %.0f GeV/#it{c}}", R, ptMin, ptMax).Data();
  string obsName = "#it{z}";

  // Plotting stuff
  TCanvas* myCanvas = new TCanvas("Plot", "Plot", 800, 800);
  myCanvas->SetLeftMargin(0.15);
  myCanvas->SetBottomMargin(0.15);
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->GetXaxis()->SetLabelSize(labelSize);
  frame->GetYaxis()->SetLabelSize(labelSize);
  frame->GetXaxis()->SetTitleSize(titleSize);
  frame->GetYaxis()->SetTitleSize(titleSize);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  // Make histograms
  hadron = "K0";
  TH1F* gluon0 = loadFragmentation(inFile, hadron, mGluon, ptMin, ptMax);
  TH1F* quark0 = loadFragmentation(inFile, hadron, mQuark, ptMin, ptMax);
  TH1F* ratio0 = (TH1F*)gluon0->Clone("ratio0_g-q");
  ratio0->Divide(quark0);
  ratio0->SetLineWidth(lineWidth);
  ratio0->SetMarkerSize(markerSize);
  ratio0->SetLineColor(GetColor(0));
  ratio0->SetMarkerStyle(GetMarker(0));
  ratio0->SetMarkerColor(GetColor(0));
  histVector.push_back(ratio0);
  legend->AddEntry(ratio0, formatHadronName(hadron).c_str());
  saveName = TString::Format("%s_g%s-q%s", saveName.c_str(), hadron.c_str(), hadron.c_str()).Data();

  hadron = "Lambda0";
  TH1F* gluon1 = loadFragmentation(inFile, hadron, mGluon, ptMin, ptMax);
  TH1F* quark1 = loadFragmentation(inFile, hadron, mQuark, ptMin, ptMax);
  TH1F* ratio1 = (TH1F*)gluon1->Clone("ratio1_g-q");
  ratio1->Divide(quark1);
  ratio1->SetLineWidth(lineWidth);
  ratio1->SetMarkerSize(markerSize);
  ratio1->SetLineColor(GetColor(1));
  ratio1->SetMarkerStyle(GetMarker(1));
  ratio1->SetMarkerColor(GetColor(1));
  histVector.push_back(ratio1);
  legend->AddEntry(ratio1, formatHadronName(hadron).c_str());
  saveName = TString::Format("%s_g%s-q%s", saveName.c_str(), hadron.c_str(), hadron.c_str()).Data();

  saveName = TString::Format("%s_pt%.0f-%.0f.pdf", saveName.c_str(), ptMin, ptMax).Data();
  plotNHists(myCanvas, frame, histVector, legend, saveName, "", latexText);
  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void plotNHists(TCanvas* canvas, TH1F* frame, std::vector<TH1F*> histVector, TLegend* legend, string saveName, string setDrawOption, string latexText)
{
  canvas->cd();
  frame->Draw();
  string drawOption = "same";
  if (setDrawOption != "") {
    drawOption = TString::Format("%s %s", drawOption.c_str(), setDrawOption.c_str()).Data();
  }
  for (auto hist : histVector) {
    hist->Draw(drawOption.c_str());
  }
  if (legend) { legend->Draw("same"); }
  if (latexText != "") { DrawLatex(0.375, 0.7, latexText.c_str(), legend->GetTextSize()); }
  TLine* line = new TLine(0, 1, 1, 1);
  line->SetLineStyle(9);
  line->Draw("same");
  canvas->SaveAs(TString::Format("./%s", saveName.c_str()).Data());
}

int getNJets(TFile* inFile, double ptMin, double ptMax, int quarkOrGluon)
{
  string projectionAxis = "X";
  TH2F* hNJetTypes = (TH2F*) inFile->Get("hNJetTypes");
  TH1F* hNJets = projectHist(hNJetTypes, projectionAxis, "nJets", ptMin, ptMax);

  int nGluons = (int) hNJets->GetBinContent(1);
  int nQuarks = (int) hNJets->GetBinContent(2);
  int nJets = nQuarks + nGluons;

  if (quarkOrGluon == mQuark) { return nQuarks; }
  else if (quarkOrGluon == mGluon) { return nGluons; }
  else if (quarkOrGluon == mInclusive) { return nJets; }
  else { return -1; }
}

double kkpParameter(string hadron, int quarkOrGluon, int parameter, double sHat)
{
  double outPar = -999.;
  if (parameter > 4) { cout << "kkpParameter: invalid input!" << endl; }
  auto f3 = [](double a, double b, double c, double d, double x)
  {
    return a + b*x + c*x*x + d*x*x*x;
  };

  if (hadron == "pi") {
    if (quarkOrGluon == mQuark) { // Currently only up quark
      switch (parameter) {
        case 0: // N
          outPar = f3(0.44809, -0.13828, -0.06951, 0.01354, sHat);
          break;
        case 1: // alpha
          outPar = f3(-1.47598, -0.30498, -0.01863, -0.12529, sHat);
          break;
        case 2: // beta
          outPar = f3(91338, 0.64145, 0.07270, -0.16989, sHat);
          break;
        case 3: // gamma
          outPar = f3(0, 0.07396, 0.07757, 0, sHat);
          break;
        default:
          break;
      }
    }
    else if (quarkOrGluon == mGluon) {
      switch (parameter) {
        case 0: // N
          outPar = f3(3.73331, -3.16946, -0.47683, 0.70270 , sHat);
          break;
        case 1: // alpha
          outPar = f3(-0.74159, -0.51377, -0.19705, -0.17917, sHat);
          break;
        case 2: // beta
          outPar = f3(2.33092, 2.03394, -0.50764, -0.08565, sHat);
          break;
        case 3: // gamma
          outPar = f3(0, 0.09466, -0.10222, 0, sHat);
          break;
        default:
          break;
      }
    }
  }
  else { cout << "kkpParameter: this hadron type is not implemented" << endl; }
  return outPar;
}

// Returns the fragmentation histogram into a hadron for a given pt range
TH1F* loadFragmentation(TFile* inFile, string hadron, int quarkOrGluon, double ptMin, double ptMax)
{
  string projectionAxis = "X", histName = "";
  int nQuarks = -1, nGluons = -1, nJets = -1;

  switch (quarkOrGluon) {
    case mQuark:
      histName = TString::Format("hqFrags_%s", hadron.c_str()).Data();
      break;
    case mGluon:
      histName = TString::Format("hgFrags_%s", hadron.c_str()).Data();
      break;
    case mInclusive:
      histName = TString::Format("hJetFrag_%s", hadron.c_str()).Data();
      break;
    default:
      std::cout << "loadFragmentation: invalid input. quarkOrGluon = " << quarkOrGluon << std::endl
        << "Allowed inputs: quarkOrGluon = " << mQuark << ", " << mGluon << ", " << mInclusive << std::endl;
      return nullptr;
  }

  TH2F* jets = loadHist<TH2F*>(inFile, histName);
  TH1F* fragmentation = projectHist(jets, projectionAxis, TString::Format("%s_projected", histName.c_str()).Data(), ptMin, ptMax);
  // fragmentation->Rebin(10);
  fragmentation->Rebin(5);
  nJets = getNJets(inFile, ptMin, ptMax, quarkOrGluon); // Get number of jets of appropriate type
  fragmentation->Scale(1./nJets, "width"); // Normalise by number of jets and bin width
  return fragmentation;
}

// Load a theory prediction for fragmentation in TF1* format
TF1* loadTheory(string hadron, int quarkOrGluon, double QSquared)
{
  // KKP: https://arxiv.org/pdf/hep-ph/0011155.pdf
  double Lambda = 213e-3; // In GeV, at NLO
  double LambdaSquared = Lambda * Lambda;
  double Q0Squared = TMath::Sqrt(2); // In GeV
  double N = 22.2815, alpha = 0.12732, beta = 6.13697, gamma = 0;
  double sHat = log( log(QSquared / LambdaSquared) / log (Q0Squared / LambdaSquared) );

  TF1* theory = new TF1("kkp", "[0]*TMath::Power(x, [1])*TMath::Power(1 - x, [2])*(1 + [3]/x)", 0, 1);
  double kkpParameters[4] = { 0, 0, 0, 0 };
  for (int iPar = 0; iPar <= 4; iPar++) {
    kkpParameters[iPar] = kkpParameter(hadron, quarkOrGluon, iPar, sHat);
  }
  theory->SetParameters(kkpParameters[0], kkpParameters[1], kkpParameters[2], kkpParameters[3]);
  cout << "parameters: " << kkpParameters[0] << ", " << kkpParameters[1] << ", " << kkpParameters[2] << ", " << kkpParameters[3] << endl;
  theory->Draw();
  return theory;
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

TH1F* makeRatio(TFile* inFile, string nominatorHadron, int nominatorType, string denominatorHadron, int denominatorType, double ptMin, double ptMax)
{
  string histName = "";

  // Format histogram name such that it will be unique
  switch (nominatorType) {
    case mQuark:
      histName = TString::Format("_q%s", nominatorHadron.c_str()).Data();
      break;
    case mGluon:
      histName = TString::Format("_g%s", nominatorHadron.c_str()).Data();
      break;
    case mInclusive:
      histName = TString::Format("_%s", nominatorHadron.c_str()).Data();
      break;
    default:
      std::cout << "makeRatio: invalid input. nominatorType = " << nominatorType << std::endl
        << "Allowed inputs: " << mQuark << ", " << mGluon << ", " << mInclusive << std::endl;
      return nullptr;
  }
  switch (denominatorType) {
    case mQuark:
      histName = TString::Format("_q%s", denominatorHadron.c_str()).Data();
      break;
    case mGluon:
      histName = TString::Format("_g%s", denominatorHadron.c_str()).Data();
      break;
    case mInclusive:
      histName = TString::Format("_%s", denominatorHadron.c_str()).Data();
      break;
    default:
      std::cout << "makeRatio: invalid input. denominatorType = " << denominatorType << std::endl
        << "Allowed inputs: " << mQuark << ", " << mGluon << ", " << mInclusive << std::endl;
      return nullptr;
  }
  if (ptMin > -900) { histName = TString::Format("_pt%.0f-%.0f", ptMin, ptMax).Data(); }

  // Make ratio histogram
  TH1F* nominator = loadFragmentation(inFile, nominatorHadron, nominatorType, ptMin, ptMax);
  TH1F* denominator = loadFragmentation(inFile, denominatorHadron, denominatorType, ptMin, ptMax);
  TH1F* fragmentation = (TH1F*)nominator->Clone(histName.c_str());
  fragmentation->Divide(denominator);
  return fragmentation;
}

// Plots the fragmentation into a single hadron species for quark, gluon, and/or inclusive sample
void plotSingleHadron(TFile* inFile, string hadron, int quarkOrGluon, string theory,
                      double ptMin, double ptMax)
{
  TH2F* inclusiveJets; TH2F* gluonJets; TH2F* quarkJets;
  TH1F* inclusiveFragmentation; TH1F* gluonFragmentation; TH1F* quarkFragmentation;
  std::vector<TH1F*> histVector; std::vector<string> histNameVector;
  std::vector<TF1*> funcVector; std::vector<string> funcNameVector;

  string saveDir = "../Plots/pthat20_80/SingleHadrons";
  string saveName = TString::Format("%s", hadron.c_str()).Data();
  string saveSuffix = "";

  // Histogram settings
  if (ptMin > -900) { saveSuffix = TString::Format("%s_pt%.0f_%.0f", saveSuffix.c_str(), ptMin, ptMax).Data(); }

  // Names for legend
  string formattedHadron = formatHadronName(hadron);
  string quarkHistName = TString::Format("q #rightarrow %s", formattedHadron.c_str()).Data();
  string gluonHistName = TString::Format("g #rightarrow %s", formattedHadron.c_str()).Data();
  string inclusiveHistName = TString::Format("incl. #rightarrow %s", formattedHadron.c_str()).Data();

  // Plot settings
  string xTitle = "#it{z}", yTitle = "#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}}";
  string histTitle = TString::Format("Fragmentation into %s", formattedHadron.c_str()).Data();
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 1e-4, yMaxFrame = 1e4;
  double xMinLegend = 0.6, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.9;
  bool setLogY = true;
  string setHistDrawOption = "hist", setFuncDrawOption = "";

  switch(quarkOrGluon) {
    case mQuark:
      quarkFragmentation = loadFragmentation(inFile, hadron, mQuark, ptMin, ptMax);
      histVector.push_back(quarkFragmentation);
      histNameVector.push_back(quarkHistName);
      saveSuffix = TString::Format("_q%s", saveSuffix.c_str()).Data();
      break;
    case mGluon:
      gluonFragmentation = loadFragmentation(inFile, hadron, mGluon, ptMin, ptMax);
      histVector.push_back(gluonFragmentation);
      histNameVector.push_back(gluonHistName);
      saveSuffix = TString::Format("_g%s", saveSuffix.c_str()).Data();
      break;
    case mQuarkGluon:
      quarkFragmentation = loadFragmentation(inFile, hadron, mQuark, ptMin, ptMax);
      histVector.push_back(quarkFragmentation);
      histNameVector.push_back(quarkHistName);
      gluonFragmentation = loadFragmentation(inFile, hadron, mGluon, ptMin, ptMax);
      histVector.push_back(gluonFragmentation);
      histNameVector.push_back(gluonHistName);
      saveSuffix = TString::Format("_qg%s", saveSuffix.c_str()).Data();
      break;
    default: // mQuarkGluonInclusive
      inclusiveFragmentation = loadFragmentation(inFile, hadron, mInclusive, ptMin, ptMax);
      histVector.push_back(inclusiveFragmentation);
      histNameVector.push_back(inclusiveHistName);
      quarkFragmentation = loadFragmentation(inFile, hadron, mQuark, ptMin, ptMax);
      histVector.push_back(quarkFragmentation);
      histNameVector.push_back(quarkHistName);
      gluonFragmentation = loadFragmentation(inFile, hadron, mGluon, ptMin, ptMax);
      histVector.push_back(gluonFragmentation);
      histNameVector.push_back(gluonHistName);
  }

  if (saveSuffix != "") { saveName = TString::Format("%s%s", saveName.c_str(), saveSuffix.c_str()).Data(); }
  saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str()).Data();
  plotNHists(histVector, histNameVector, funcVector, funcNameVector,
             xTitle, yTitle, histTitle, saveName,
             xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
             xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
             setLogY, setHistDrawOption, setFuncDrawOption);
}

void plotMultipleHadrons(TFile* inFile, std::vector<string> hadrons, int quarkOrGluon, string theory,
                         double ptMin, double ptMax)
{
  std::vector<TH1F*> histVector; std::vector<string> histNameVector;
  std::vector<TF1*> funcVector; std::vector<string> funcNameVector;

  string saveDir = "../Plots/pthat20_80/MultipleHadrons";
  string saveName = "";
  string saveSuffix = "";

  // Plot settings
  string xTitle = "#it{z}", yTitle = "#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}}";
  string histTitle = "";
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 1e-5, yMaxFrame = 1e3;
  double xMinLegend = 0.6, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.9;
  bool setLogY = true;
  string setHistDrawOption = "hist", setFuncDrawOption = "";

  if (quarkOrGluon == mQuark) {
    histTitle = "Quark jet fragmentation";
    saveSuffix = TString::Format("_q%s", saveSuffix.c_str()).Data();
  }
  else if (quarkOrGluon == mGluon) {
    histTitle = "Gluon jet fragmentation";
    saveSuffix = TString::Format("_g%s", saveSuffix.c_str()).Data();
  }
  else {
    histTitle = "Jet fragmentation";
  }

  if (ptMin > -900) { saveSuffix = TString::Format("%s_pt%.0f_%.0f", saveSuffix.c_str(), ptMin, ptMax).Data(); }

  for (auto hadron : hadrons) {
    saveName = TString::Format("%s_%s", saveName.c_str(), hadron.c_str()).Data();
    string histName = formatHadronName(hadron);
    TH1F* fragmentation = loadFragmentation(inFile, hadron, quarkOrGluon, ptMin, ptMax);
    histVector.push_back(fragmentation);
    histNameVector.push_back(histName);
  }

  if (saveSuffix != "") { saveName = TString::Format("%s%s", saveName.c_str(), saveSuffix.c_str()).Data(); }
  saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str()).Data();
  plotNHists(histVector, histNameVector, funcVector, funcNameVector,
             xTitle, yTitle, histTitle, saveName,
             xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
             xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
             setLogY, setHistDrawOption, setFuncDrawOption);
}

void plotHadronRatios(TFile* inFile, std::vector<string> nominators, std::vector<string> denominators,
                      std::vector<int> nominatorTypes, std::vector<int> denominatorTypes,
                      string theory = "",
                      double ptMin, double ptMax)
{
  // TODO: Add check to see if all input vectors are the same length
  std::vector<TH1F*> histVector; std::vector<string> histNameVector;
  std::vector<TF1*> funcVector; std::vector<string> funcNameVector;

  // string saveDir = "../Plots/Pythia/pthat20_80/HadronRatios";
  string saveDir = ".";
  string saveName = "";
  string saveSuffix = "";

  // Plot settings
  string xTitle = "#it{z}";
  string yTitle = "";
  string histTitle = "Jet fragmentation ratio";
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0e-2, yMaxFrame = 1e0;
  double xMinLegend = 0.6, xMaxLegend = 0.8, yMinLegend = 0.7, yMaxLegend = 0.9;
  bool setLogY = false;
  string setHistDrawOption = "", setFuncDrawOption = "";

  if (ptMin > -900) {
    saveSuffix = TString::Format("%s_pt%.0f_%.0f", saveSuffix.c_str(), ptMin, ptMax).Data();
    histTitle = TString::Format("%s (%.0f - %.0f GeV)", histTitle.c_str(), ptMin, ptMax).Data();
  }

  for (unsigned int i = 0; i < nominators.size(); i++) {
    string nomName = nominators[i], denomName = denominators[i]; string histName;
    int nomType = nominatorTypes[i], denomType = denominatorTypes[i];
    string nomTypeName, denomTypeName;
    if (nomType == mQuark) { nomTypeName = "q"; }
    else { nomTypeName = "g"; }
    if (denomTypeName == mQuark) { denomTypeName = "q"; }
    else { denomTypeName = "g"; }

    // Divide the histograms
    TH1F* nominator = loadFragmentation(inFile, nomName, nomType, ptMin, ptMax);
    TH1F* denominator = loadFragmentation(inFile, denomName, denomType, ptMin, ptMax);
    TH1F* fragmentation = (TH1F*)nominator->Clone("fragmentation");
    fragmentation->Divide(denominator);
    fragmentation->SetLineWidth(3);

    if (nomType == denomType) { // Same parton
      saveName = TString::Format("%s_%s-%s", saveName.c_str(), nomName.c_str(), denomName.c_str()).Data();
      histName = TString::Format("%s: %s/%s", nomTypeName.c_str(), formatHadronName(nomName).c_str(), formatHadronName(denomName).c_str()).Data();
    } // Same parton
    else { // Different parton
      saveName = TString::Format("%s_%s%s-%s%s", saveName.c_str(), nomTypeName.c_str(), nomName.c_str(), denomTypeName.c_str(), denomName.c_str()).Data();
      if (nomName == denomName) { // Same hadron
        histName = TString::Format("%s: %s/%s", formatHadronName(nomName).c_str(), nomTypeName.c_str(), denomTypeName.c_str()).Data();
      }
      else { // Different hadron
        histName = TString::Format("(%s/%s) / (%s/%s)", formatHadronName(nomName).c_str(), nomTypeName.c_str(), formatHadronName(denomName).c_str(), denomTypeName.c_str()).Data();
      }
    } // Different parton

    histVector.push_back(fragmentation);
    // histName = TString::Format("%s", formatHadronName(nomName).c_str());
    if (nomType == mQuark) { histName = "quark"; }
    else { histName = "gluon"; }
    yTitle = TString::Format("#it{N}(%s) / #it{N}(%s)", formatHadronName(nomName).c_str(), formatHadronName(denomName).c_str()).Data();
    histNameVector.push_back(histName);
    // delete nominator; delete denominator; delete fragmentation;
  } // for i < nominators.size()

  // if (quarkOrGluon == mQuarkGluon) { saveName = TString::Format("%s_q-g", saveName.c_str()).Data(); }
  // else if (quarkOrGluon < 0 ) { saveName = TString::Format("%s_g-q", saveName.c_str()).Data(); }

  if (saveSuffix != "") {
    saveName = TString::Format("%s%s", saveName.c_str(), saveSuffix.c_str()).Data();
  }
  saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str()).Data();
  plotNHists(histVector, histNameVector, funcVector, funcNameVector,
             xTitle, yTitle, histTitle, saveName,
             xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
             xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
             setLogY, setHistDrawOption, setFuncDrawOption);
}

void plotGluonQuarkJetsFraction(TFile* inFile)
{
  int nPtBins = 20;
  double ptBins[21] = { 0., 10, 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200. };
  TH1F* gluonFraction = new TH1F("gluonFraction", "Gluon Fraction", nPtBins, ptBins);
  TH1F* quarkFraction = new TH1F("quarkFraction", "Quark Fraction", nPtBins, ptBins);
  int nGluons = -1, nQuarks = -1, nJets = -1;
  double fGluons = -1, fQuarks = -1;
  double ptMin = -999, ptMax = -999;
  double ptMean = -999;

  std::vector<TH1F*> histVector; std::vector<string> histNameVector;
  std::vector<TF1*> funcVector; std::vector<string> funcNameVector;
  string saveDir = "../Plots/";
  string saveName = "qgJetsFraction";
  string saveSuffix = "";
  // Plot settings
  string xTitle = "#it{p}_{T}", yTitle = "";
  string histTitle = "Relative quark and gluon fraction";
  double xMinFrame = 0, xMaxFrame = ptBins[nPtBins], yMinFrame = 0e-2, yMaxFrame = 1e0;
  double xMinLegend = 0.6, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.9;
  bool setLogY = false;
  string setHistDrawOption = "hist", setFuncDrawOption = "";

  for (auto iPt = 0; iPt < nPtBins; iPt++) {
    ptMin = ptBins[iPt];
    ptMax = ptBins[iPt+1];
    ptMean = (ptMin + ptMax)/2;
    nGluons = getNJets(inFile, ptMin, ptMax, mGluon);
    nQuarks = getNJets(inFile, ptMin, ptMax, mQuark);
    nJets = getNJets(inFile, ptMin, ptMax, mInclusive);
    if (nJets != 0) {
      fGluons = ( 1. * nGluons) / (1. * nJets);
      fQuarks = ( 1. * nQuarks) / (1. * nJets);
    }
    else {
      fGluons = 0;
      fQuarks = 0;
    }
    cout << "pT = (" << ptMin << ", " << ptMax << ")" << endl
      << "nJets = " << nJets << " = " << nGluons << " + " << nQuarks << endl
      << fGluons + fQuarks << " = " << fGluons << " + " << fQuarks << endl;
    gluonFraction->Fill(ptMean, fGluons);
    quarkFraction->Fill(ptMean, fQuarks);
  }

  histVector.push_back(gluonFraction);
  histNameVector.push_back("g");
  histVector.push_back(quarkFraction);
  histNameVector.push_back("q");

  saveName = TString::Format("%s%s", saveDir.c_str(), saveName.c_str()).Data();

  plotNHists(histVector, histNameVector, funcVector, funcNameVector,
             xTitle, yTitle, histTitle, saveName,
             xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
             xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
             setLogY, setHistDrawOption, setFuncDrawOption);
}
