
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

string formatHadronName(string hadron);

int getNJets(TFile* inFile, double ptMin, double ptMax, int quarkOrGluon);

TF1* loadKKP();
TH1F* loadFragmentation(TFile* inFile, string hadron, int quarkOrGluon, double ptMin, double ptMax);

TH1F* makeRatio(TFile* inFile, string nominatorHadron, int nominatorType, string denominatorHadron, int denominatorType, double ptMin, double ptMax);

void plotHadronRatios(TFile* inFile, std::vector<string> nominators, std::vector<string> denominators,
                      std::vector<int> nominatorTypes, std::vector<int> denominatorTypes,
                      string theory = "",
                      double ptMin = -999, double ptMax = -999);
void plotMultipleHadrons(TFile* inFile, std::vector<string> hadrons, int quarkOrGluon, string theory = "",
                         double ptMin = -999, double ptMax = -999);
void plotSingleHadron(TFile* inFile, string hadron, int quarkOrGluon, string theory = "",
                      double ptMin = -999, double ptMax = -999);

TH1F* projectHist(TH2F* inputHist, string projectionAxis, string histName, double axMin, double axMax);

// TODO: Add theory (selection from various sets; give a fragmentation scale)

void plotFrags(void)
{
  double time = clock();
  gROOT->SetBatch();
  gStyle->SetNdivisions(505);
  string inName = "PythiaResultJob12128342_66_pthat80_200";
  TFile *inFile = TFile::Open(TString::Format("../data/pythia/%s.root", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  std::vector<double> jetPtBins = {40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};
	std::vector<string> partons = {"q", "g"};

  // Quark vs gluon (with and without inclusive)
  // for (unsigned int i = 0; i < hadrons.size(); i++) {
  //   plotSingleHadron(inFile, hadrons[i], mQuarkGluon);
  //   plotSingleHadron(inFile, hadrons[i], mQuarkGluonInclusive);
  // }

  // Multiple hadrons in one plot
  // std::vector<string> hadronsVector = { "pi", "p", "K0", "Lambda0" };
  // plotMultipleHadrons(inFile, hadronsVector, mQuark);
  // plotMultipleHadrons(inFile, hadronsVector, mGluon);
  // plotMultipleHadrons(inFile, hadronsVector, mInclusive);

  // std::vector<string> nominators = { hadrons[2], hadrons[7], hadrons[0], hadrons[1], hadrons[5], hadrons[5], hadrons[7] };
  // std::vector<string> denominators = { hadrons[0], hadrons[5], hadrons[0], hadrons[6], hadrons[1], hadrons[6], hadrons[6] };
  // std::vector<string> nominators = { hadrons[2], hadrons[7] };
  // std::vector<string> denominators = { hadrons[0], hadrons[6] };
  // plotHadronRatios(inFile, nominators, denominators, mQuark);
  // plotHadronRatios(inFile, nominators, denominators, mGluon);
  // plotHadronRatios(inFile, nominators, denominators, mInclusive);


  // Gluon over quark
  // std::vector<string> gluonNominators = { "pi", "Lambda0" };
  // std::vector<int> gluonNominatorTypes = { mGluon, mGluon };
  // std::vector<string> quarkDenominators = { "pi", "Lambda0" };
  // std::vector<int> quarkDenominatorTypes = { mQuark, mQuark };
  // plotHadronRatios(inFile, gluonNominators, quarkDenominators, gluonNominatorTypes, quarkDenominatorTypes); // No pt selection
  // plotHadronRatios(inFile, gluonNominators, quarkDenominators, gluonNominatorTypes, quarkDenominatorTypes, "", 5, 10);
  // plotHadronRatios(inFile, gluonNominators, quarkDenominators, gluonNominatorTypes, quarkDenominatorTypes, "", 10, 15);
  // plotHadronRatios(inFile, gluonNominators, quarkDenominators, gluonNominatorTypes, quarkDenominatorTypes, "", 100, 150);
  // plotHadronRatios(inFile, gluonNominators, quarkDenominators, gluonNominatorTypes, quarkDenominatorTypes, "", 500, 700);

  // Lambda / K0
  std::vector<string> nominators = { hadrons[7], hadrons[7] };
  std::vector<int> nominatorTypes = { mQuark, mGluon };
  std::vector<string> denominators = { hadrons[6], hadrons[6] };
  std::vector<int> denominatorTypes = { mQuark, mGluon };
  plotHadronRatios(inFile, nominators, denominators, nominatorTypes, denominatorTypes); // No pt selection
  plotHadronRatios(inFile, nominators, denominators, nominatorTypes, denominatorTypes, "", 5, 10);
  plotHadronRatios(inFile, nominators, denominators, nominatorTypes, denominatorTypes, "", 10, 15);
  plotHadronRatios(inFile, nominators, denominators, nominatorTypes, denominatorTypes, "", 50, 60);
  plotHadronRatios(inFile, nominators, denominators, nominatorTypes, denominatorTypes, "", 60, 70);
  plotHadronRatios(inFile, nominators, denominators, nominatorTypes, denominatorTypes, "", 100, 150);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

// TH1F* projectHist(TH2F* inputHist, string projectionAxis, string histName, double axMin, double axMax)
// {
//   TH1F* outputHist;
//   TH2F* inputClone = (TH2F*)inputHist->Clone("CloneOfInput");
//   int firstBin = 1; int lastBin = -1;

//   if (projectionAxis == "X" || projectionAxis == "x") {
//     if (axMin > -900) { firstBin = inputClone->GetYaxis()->FindBin(axMin); }
//     if (axMax > -900) { lastBin = inputClone->GetYaxis()->FindBin(axMax); }
//     else { lastBin = inputClone->GetNbinsY() + 1; }
//     if (firstBin == lastBin) {
//       if (firstBin == inputClone->GetNbinsY() + 1 || lastBin == 0) {
//         cout << "projectHist: requested bins out of range for Y axis. Aborting" << endl
//           << "Requested: (" << axMin << ", " << axMax << ")" << endl
//           << "Bins: (" << firstBin << ", " << lastBin << ")" << endl
//           << "Max bin: " << inputClone->GetNbinsY() + 1;
//         return nullptr;
//       }
//     }
//     outputHist = (TH1F*)inputClone->ProjectionX(TString::Format("%s", histName.c_str()).Data(), firstBin, lastBin);
//   }
//   else if (projectionAxis == "Y" || projectionAxis == "y") {
//     if (axMin > -900) { firstBin = inputClone->GetXaxis()->FindBin(axMin); }
//     if (axMax > -900) { lastBin = inputClone->GetXaxis()->FindBin(axMax); }
//     else { lastBin = inputClone->GetNbinsX() + 1; }
//     if (firstBin == lastBin) {
//       if (firstBin == inputClone->GetNbinsX() + 1 || lastBin == 0) {
//         cout << "projectHist: requested bins out of range for X axis. Aborting" << endl
//           << "Requested: (" << axMin << ", " << axMax << ")" << endl
//           << "Bins: (" << firstBin << ", " << lastBin << ")" << endl
//           << "Max bin: " << inputClone->GetNbinsX() + 1;
//         return nullptr;
//       }
//     }
//     outputHist = (TH1F*)inputClone->ProjectionY(TString::Format("%s", histName.c_str()).Data(), firstBin, lastBin);
//   }
//   else {
//     cout << "make_projection: invalid projection axis " << projectionAxis << ". Aborting." << endl;
//     return nullptr;
//   }
//   delete inputClone;
//   return outputHist;
// }

int getNJets(TFile* inFile, double ptMin, double ptMax, int quarkOrGluon)
{
  string projectionAxis = "X";
  // double zMin = -999., zMax = -999;

  TH2F* hNJetTypes = (TH2F*) inFile->Get("hNJetTypes");
  // TH1F* hNJets = projectHist(hNJetTypes, projectionAxis, "nJets", zMin, zMax, ptMin, ptMax);
  TH1F* hNJets = projectHist(hNJetTypes, projectionAxis, "nJets", ptMin, ptMax);

  int nGluons = (int) hNJets->GetBinContent(1);
  int nQuarks = (int) hNJets->GetBinContent(2);
  int nJets = nQuarks + nGluons;

  if (quarkOrGluon == mQuark) { return nQuarks; }
  else if (quarkOrGluon == mGluon) { return nGluons; }
  else if (quarkOrGluon == mInclusive) { return nJets; }
  else { return -1; }
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
  nJets = getNJets(inFile, ptMin, ptMax, quarkOrGluon); // Get number of jets of appropriate type
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

  string saveDir = "../Plots/SingleHadrons";
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

  string saveDir = "../Plots/MultipleHadrons";
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

  string saveDir = "../Plots/HadronRatios";
  string saveName = "";
  string saveSuffix = "";

  // Plot settings
  string xTitle = "#it{z}", yTitle = "";
  string histTitle = "Jet fragmentation ratio";
  double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0e-2, yMaxFrame = 1e0;
  double xMinLegend = 0.6, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.9;
  bool setLogY = false;
  string setHistDrawOption = "hist", setFuncDrawOption = "";

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
