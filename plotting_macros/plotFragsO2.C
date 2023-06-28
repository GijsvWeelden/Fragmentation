
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

TF1* loadKKP();
TH1F* loadFragmentation(TDirectory* inDir, int dataOrMC, int projOrPt, double ptMin, double ptMax);
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

void plotJetQA(TDirectory* inDir, int dataOrMC, double ptMin, double ptMax, string saveDir);
void plotTrackQA(TDirectory* inDir, int dataOrMC, double ptMin, double ptMax, string saveDir);

// TH1F* projectHist(TH2F* inputHist, string projectionAxis, string histName, double axMin, double axMax);

// TODO: Add theory (selection from various sets; give a fragmentation scale)

void plotFragsData(void)
{
  double time = clock();
  gROOT->SetBatch();
  gStyle->SetNdivisions(505);
  string inName = "LHC22o_pass4_small-train90120.root";
  string saveDir = "../Plots/LHC22o_pass4_small";
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TDirectory* dir = (TDirectory*)inFile->Get("jet-fragmentation");
  TDirectory* jetsDir = (TDirectory*)dir->Get("jets");
  TDirectory* tracksDir = (TDirectory*)dir->Get("tracks");

  // comparePtBins(jetsDir, kData, kPt, saveDir);
  // comparePtBins(jetsDir, kData, kProj, saveDir);

  // compareZDefinitions(jetsDir, kData, 20, 40, saveDir);
  // compareZDefinitions(jetsDir, kData, 40, 60, saveDir);
  // compareZDefinitions(jetsDir, kData, 60, 80, saveDir);

  // plotJetQA(jetsDir, kData, -999, -999, saveDir);
  // plotJetQA(jetsDir, kData, 20, 40, saveDir);
  // plotJetQA(jetsDir, kData, 40, 60, saveDir);
  // plotJetQA(jetsDir, kData, 60, 80, saveDir);

  plotTrackQA(tracksDir, kData, -999, -999, saveDir);

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

// Comparison of jet pt bins
void comparePtBins(TDirectory* inDir, int dataOrMC, int projOrPt, string saveDir)
{
  std::vector<TH1F*> histVector; std::vector<string> nameVector;
  std::vector<TF1*> funcVector; std::vector<string> funcNameVector;
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
  std::vector<TF1*> funcVector; std::vector<string> funcNameVector;
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
  // h3Jets_2->GetXaxis()->SetRange(1, h3Jets_2->GetNbinsX()); // FIXME: This doesn't work!!
  h3Jets_2->GetXaxis()->SetRangeUser(ptMin, ptMax); // FIXME: This doesn't work!!
  // TH2F* hJetEtaPhi = dynamic_cast<TH2F*>(h3Jets_2->Project3D("zy"));
  TH2F* hJetEtaPhi = (TH2F*)h3Jets_2->Project3D("zy");
  hJetEtaPhi->Scale(1./hJetEtaPhi->Integral());
  hJetEtaPhi->SetName(TString::Format("etaphi%.0f_%.0f", ptMin, ptMax));
  saveName = TString::Format("%s/jetEtaPhi_pt%.0f-%.0f", saveDir.c_str(), ptMin, ptMax);
  plotOneHist(hJetEtaPhi, xTitle, yTitle, histTitle, saveName,
              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
              setLogY, drawLegend);
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

// TH1F* makeRatio(TFile* inFile, string nominatorHadron, int nominatorType, string denominatorHadron, int denominatorType, double ptMin, double ptMax)
// {
//   string histName = "";

//   // Format histogram name such that it will be unique
//   switch (nominatorType) {
//     case mQuark:
//       histName = TString::Format("_q%s", nominatorHadron.c_str()).Data();
//       break;
//     case mGluon:
//       histName = TString::Format("_g%s", nominatorHadron.c_str()).Data();
//       break;
//     case mInclusive:
//       histName = TString::Format("_%s", nominatorHadron.c_str()).Data();
//       break;
//     default:
//       std::cout << "makeRatio: invalid input. nominatorType = " << nominatorType << std::endl
//         << "Allowed inputs: " << mQuark << ", " << mGluon << ", " << mInclusive << std::endl;
//       return nullptr;
//   }
//   switch (denominatorType) {
//     case mQuark:
//       histName = TString::Format("_q%s", denominatorHadron.c_str()).Data();
//       break;
//     case mGluon:
//       histName = TString::Format("_g%s", denominatorHadron.c_str()).Data();
//       break;
//     case mInclusive:
//       histName = TString::Format("_%s", denominatorHadron.c_str()).Data();
//       break;
//     default:
//       std::cout << "makeRatio: invalid input. denominatorType = " << denominatorType << std::endl
//         << "Allowed inputs: " << mQuark << ", " << mGluon << ", " << mInclusive << std::endl;
//       return nullptr;
//   }
//   if (ptMin > -900) { histName = TString::Format("_pt%.0f-%.0f", ptMin, ptMax).Data(); }

//   // Make ratio histogram
//   TH1F* nominator = loadFragmentation(inFile, nominatorHadron, nominatorType, ptMin, ptMax);
//   TH1F* denominator = loadFragmentation(inFile, denominatorHadron, denominatorType, ptMin, ptMax);
//   TH1F* fragmentation = (TH1F*)nominator->Clone(histName.c_str());
//   fragmentation->Divide(denominator);
//   return fragmentation;
// }

// // Plots the fragmentation into a single hadron species for quark, gluon, and/or inclusive sample
// void plotSingleHadron(TFile* inFile, string hadron, int quarkOrGluon, string theory,
//                       double ptMin, double ptMax)
// {
//   TH2F* inclusiveJets; TH2F* gluonJets; TH2F* quarkJets;
//   TH1F* inclusiveFragmentation; TH1F* gluonFragmentation; TH1F* quarkFragmentation;
//   std::vector<TH1F*> histVector; std::vector<string> histNameVector;
//   std::vector<TF1*> funcVector; std::vector<string> funcNameVector;

//   string saveDir = "../Plots/SingleHadrons";
//   string saveName = TString::Format("%s", hadron.c_str()).Data();
//   string saveSuffix = "";

//   // Histogram settings
//   if (ptMin > -900) { saveSuffix = TString::Format("%s_pt%.0f_%.0f", saveSuffix.c_str(), ptMin, ptMax).Data(); }

//   // Names for legend
//   string formattedHadron = formatHadronName(hadron);
//   string quarkHistName = TString::Format("q #rightarrow %s", formattedHadron.c_str()).Data();
//   string gluonHistName = TString::Format("g #rightarrow %s", formattedHadron.c_str()).Data();
//   string inclusiveHistName = TString::Format("incl. #rightarrow %s", formattedHadron.c_str()).Data();

//   // Plot settings
//   string xTitle = "#it{z}", yTitle = "#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}}";
//   string histTitle = TString::Format("Fragmentation into %s", formattedHadron.c_str()).Data();
//   double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 1e-4, yMaxFrame = 1e4;
//   double xMinLegend = 0.6, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.9;
//   bool setLogY = true;
//   string setHistDrawOption = "hist", setFuncDrawOption = "";

//   switch(quarkOrGluon) {
//     case mQuark:
//       quarkFragmentation = loadFragmentation(inFile, hadron, mQuark, ptMin, ptMax);
//       histVector.push_back(quarkFragmentation);
//       histNameVector.push_back(quarkHistName);
//       saveSuffix = TString::Format("_q%s", saveSuffix.c_str()).Data();
//       break;
//     case mGluon:
//       gluonFragmentation = loadFragmentation(inFile, hadron, mGluon, ptMin, ptMax);
//       histVector.push_back(gluonFragmentation);
//       histNameVector.push_back(gluonHistName);
//       saveSuffix = TString::Format("_g%s", saveSuffix.c_str()).Data();
//       break;
//     case mQuarkGluon:
//       quarkFragmentation = loadFragmentation(inFile, hadron, mQuark, ptMin, ptMax);
//       histVector.push_back(quarkFragmentation);
//       histNameVector.push_back(quarkHistName);
//       gluonFragmentation = loadFragmentation(inFile, hadron, mGluon, ptMin, ptMax);
//       histVector.push_back(gluonFragmentation);
//       histNameVector.push_back(gluonHistName);
//       saveSuffix = TString::Format("_qg%s", saveSuffix.c_str()).Data();
//       break;
//     default: // mQuarkGluonInclusive
//       inclusiveFragmentation = loadFragmentation(inFile, hadron, mInclusive, ptMin, ptMax);
//       histVector.push_back(inclusiveFragmentation);
//       histNameVector.push_back(inclusiveHistName);
//       quarkFragmentation = loadFragmentation(inFile, hadron, mQuark, ptMin, ptMax);
//       histVector.push_back(quarkFragmentation);
//       histNameVector.push_back(quarkHistName);
//       gluonFragmentation = loadFragmentation(inFile, hadron, mGluon, ptMin, ptMax);
//       histVector.push_back(gluonFragmentation);
//       histNameVector.push_back(gluonHistName);
//   }

//   if (saveSuffix != "") { saveName = TString::Format("%s%s", saveName.c_str(), saveSuffix.c_str()).Data(); }
//   saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str()).Data();
//   plotNHists(histVector, histNameVector, funcVector, funcNameVector,
//              xTitle, yTitle, histTitle, saveName,
//              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
//              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
//              setLogY, setHistDrawOption, setFuncDrawOption);
// }

// void plotMultipleHadrons(TFile* inFile, std::vector<string> hadrons, int quarkOrGluon, string theory,
//                          double ptMin, double ptMax)
// {
//   std::vector<TH1F*> histVector; std::vector<string> histNameVector;
//   std::vector<TF1*> funcVector; std::vector<string> funcNameVector;

//   string saveDir = "../Plots/MultipleHadrons";
//   string saveName = "";
//   string saveSuffix = "";

//   // Plot settings
//   string xTitle = "#it{z}", yTitle = "#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}}";
//   string histTitle = "";
//   double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 1e-5, yMaxFrame = 1e3;
//   double xMinLegend = 0.6, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.9;
//   bool setLogY = true;
//   string setHistDrawOption = "hist", setFuncDrawOption = "";

//   if (quarkOrGluon == mQuark) {
//     histTitle = "Quark jet fragmentation";
//     saveSuffix = TString::Format("_q%s", saveSuffix.c_str()).Data();
//   }
//   else if (quarkOrGluon == mGluon) {
//     histTitle = "Gluon jet fragmentation";
//     saveSuffix = TString::Format("_g%s", saveSuffix.c_str()).Data();
//   }
//   else {
//     histTitle = "Jet fragmentation";
//   }

//   if (ptMin > -900) { saveSuffix = TString::Format("%s_pt%.0f_%.0f", saveSuffix.c_str(), ptMin, ptMax).Data(); }

//   for (auto hadron : hadrons) {
//     saveName = TString::Format("%s_%s", saveName.c_str(), hadron.c_str()).Data();
//     string histName = formatHadronName(hadron);
//     TH1F* fragmentation = loadFragmentation(inFile, hadron, quarkOrGluon, ptMin, ptMax);
//     histVector.push_back(fragmentation);
//     histNameVector.push_back(histName);
//   }

//   if (saveSuffix != "") { saveName = TString::Format("%s%s", saveName.c_str(), saveSuffix.c_str()).Data(); }
//   saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str()).Data();
//   plotNHists(histVector, histNameVector, funcVector, funcNameVector,
//              xTitle, yTitle, histTitle, saveName,
//              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
//              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
//              setLogY, setHistDrawOption, setFuncDrawOption);
// }

// void plotHadronRatios(TFile* inFile, std::vector<string> nominators, std::vector<string> denominators,
//                       std::vector<int> nominatorTypes, std::vector<int> denominatorTypes,
//                       string theory = "",
//                       double ptMin, double ptMax)
// {
//   // TODO: Add check to see if all input vectors are the same length
//   std::vector<TH1F*> histVector; std::vector<string> histNameVector;
//   std::vector<TF1*> funcVector; std::vector<string> funcNameVector;

//   string saveDir = "../Plots/HadronRatios";
//   string saveName = "";
//   string saveSuffix = "";

//   // Plot settings
//   string xTitle = "#it{z}", yTitle = "";
//   string histTitle = "Jet fragmentation ratio";
//   double xMinFrame = 0, xMaxFrame = 1, yMinFrame = 0e-2, yMaxFrame = 1e0;
//   double xMinLegend = 0.6, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.9;
//   bool setLogY = false;
//   string setHistDrawOption = "hist", setFuncDrawOption = "";

//   if (ptMin > -900) {
//     saveSuffix = TString::Format("%s_pt%.0f_%.0f", saveSuffix.c_str(), ptMin, ptMax).Data();
//     histTitle = TString::Format("%s (%.0f - %.0f GeV)", histTitle.c_str(), ptMin, ptMax).Data();
//   }

//   for (unsigned int i = 0; i < nominators.size(); i++) {
//     string nomName = nominators[i], denomName = denominators[i]; string histName;
//     int nomType = nominatorTypes[i], denomType = denominatorTypes[i];
//     string nomTypeName, denomTypeName;
//     if (nomType == mQuark) { nomTypeName = "q"; }
//     else { nomTypeName = "g"; }
//     if (denomTypeName == mQuark) { denomTypeName = "q"; }
//     else { denomTypeName = "g"; }

//     // Divide the histograms
//     TH1F* nominator = loadFragmentation(inFile, nomName, nomType, ptMin, ptMax);
//     TH1F* denominator = loadFragmentation(inFile, denomName, denomType, ptMin, ptMax);
//     TH1F* fragmentation = (TH1F*)nominator->Clone("fragmentation");
//     fragmentation->Divide(denominator);
//     fragmentation->SetLineWidth(3);

//     if (nomType == denomType) { // Same parton
//       saveName = TString::Format("%s_%s-%s", saveName.c_str(), nomName.c_str(), denomName.c_str()).Data();
//       histName = TString::Format("%s: %s/%s", nomTypeName.c_str(), formatHadronName(nomName).c_str(), formatHadronName(denomName).c_str()).Data();
//     } // Same parton
//     else { // Different parton
//       saveName = TString::Format("%s_%s%s-%s%s", saveName.c_str(), nomTypeName.c_str(), nomName.c_str(), denomTypeName.c_str(), denomName.c_str()).Data();
//       if (nomName == denomName) { // Same hadron
//         histName = TString::Format("%s: %s/%s", formatHadronName(nomName).c_str(), nomTypeName.c_str(), denomTypeName.c_str()).Data();
//       }
//       else { // Different hadron
//         histName = TString::Format("(%s/%s) / (%s/%s)", formatHadronName(nomName).c_str(), nomTypeName.c_str(), formatHadronName(denomName).c_str(), denomTypeName.c_str()).Data();
//       }
//     } // Different parton

//     histVector.push_back(fragmentation);
//     histNameVector.push_back(histName);
//     // delete nominator; delete denominator; delete fragmentation;
//   } // for i < nominators.size()

//   // if (quarkOrGluon == mQuarkGluon) { saveName = TString::Format("%s_q-g", saveName.c_str()).Data(); }
//   // else if (quarkOrGluon < 0 ) { saveName = TString::Format("%s_g-q", saveName.c_str()).Data(); }

//   if (saveSuffix != "") {
//     saveName = TString::Format("%s%s", saveName.c_str(), saveSuffix.c_str()).Data();
//   }
//   saveName = TString::Format("%s/%s", saveDir.c_str(), saveName.c_str()).Data();
//   plotNHists(histVector, histNameVector, funcVector, funcNameVector,
//              xTitle, yTitle, histTitle, saveName,
//              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
//              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
//              setLogY, setHistDrawOption, setFuncDrawOption);
// }

// void plotGluonQuarkJetsFraction(TFile* inFile)
// {
//   int nPtBins = 20;
//   double ptBins[21] = { 0., 10, 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200. };
//   TH1F* gluonFraction = new TH1F("gluonFraction", "Gluon Fraction", nPtBins, ptBins);
//   TH1F* quarkFraction = new TH1F("quarkFraction", "Quark Fraction", nPtBins, ptBins);
//   int nGluons = -1, nQuarks = -1, nJets = -1;
//   double fGluons = -1, fQuarks = -1;
//   double ptMin = -999, ptMax = -999;
//   double ptMean = -999;

//   std::vector<TH1F*> histVector; std::vector<string> histNameVector;
//   std::vector<TF1*> funcVector; std::vector<string> funcNameVector;
//   string saveDir = "../Plots/";
//   string saveName = "qgJetsFraction";
//   string saveSuffix = "";
//   // Plot settings
//   string xTitle = "#it{p}_{T}", yTitle = "";
//   string histTitle = "Relative quark and gluon fraction";
//   double xMinFrame = 0, xMaxFrame = ptBins[nPtBins], yMinFrame = 0e-2, yMaxFrame = 1e0;
//   double xMinLegend = 0.6, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.9;
//   bool setLogY = false;
//   string setHistDrawOption = "hist", setFuncDrawOption = "";

//   for (auto iPt = 0; iPt < nPtBins; iPt++) {
//     ptMin = ptBins[iPt];
//     ptMax = ptBins[iPt+1];
//     ptMean = (ptMin + ptMax)/2;
//     nGluons = getNJets(inFile, ptMin, ptMax, mGluon);
//     nQuarks = getNJets(inFile, ptMin, ptMax, mQuark);
//     nJets = getNJets(inFile, ptMin, ptMax, mInclusive);
//     if (nJets != 0) {
//       fGluons = ( 1. * nGluons) / (1. * nJets);
//       fQuarks = ( 1. * nQuarks) / (1. * nJets);
//     }
//     else {
//       fGluons = 0;
//       fQuarks = 0;
//     }
//     cout << "pT = (" << ptMin << ", " << ptMax << ")" << endl
//       << "nJets = " << nJets << " = " << nGluons << " + " << nQuarks << endl
//       << fGluons + fQuarks << " = " << fGluons << " + " << fQuarks << endl;
//     gluonFraction->Fill(ptMean, fGluons);
//     quarkFraction->Fill(ptMean, fQuarks);
//   }

//   histVector.push_back(gluonFraction);
//   histNameVector.push_back("g");
//   histVector.push_back(quarkFraction);
//   histNameVector.push_back("q");

//   saveName = TString::Format("%s%s", saveDir.c_str(), saveName.c_str()).Data();

//   plotNHists(histVector, histNameVector, funcVector, funcNameVector,
//              xTitle, yTitle, histTitle, saveName,
//              xMinFrame, xMaxFrame, yMinFrame, yMaxFrame,
//              xMinLegend, xMaxLegend, yMinLegend, yMaxLegend,
//              setLogY, setHistDrawOption, setFuncDrawOption);
// }
