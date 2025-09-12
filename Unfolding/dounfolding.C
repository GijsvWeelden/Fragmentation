#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "../plotting_macros/histUtils.C"

#ifndef DO_UNFOLDING
#define DO_UNFOLDING

namespace verbositylvls {
  enum Verbosity {kErrors, kWarnings, kInfo, kDebug, kDebugMax};
  bool is_valid(int v) { return (v >= kErrors && v <= kDebugMax); }
  string to_string(Verbosity v) {
    switch (v) {
      case kErrors:   return "kErrors";
      case kWarnings: return "kWarnings";
      case kInfo:     return "kInfo";
      case kDebug:    return "kDebug";
      case kDebugMax: return "kDebugMax";
      default:        return "Unknown";
    }
  }
  bool passVerbosityCheck(Verbosity level, Verbosity threshold) {
    return (is_valid(level) && level <= threshold);
  }
}

namespace rmutilities {
  const int axisJetRmPtJetRec = 0;
  const int axisJetRmPtJetGen = 1;
  const int ndimJetRm         = 2;

  const int axisv0PtRmPtJetGen = 0;
  const int axisv0PtRmV0GenPt  = 1;
  const int axisv0PtRmPtJetRec = 2;
  const int axisv0PtRmV0RecPt  = 3;
  const int nDimV0PtRm         = 4;

  const int axisv0ZRmPtJetRec  = 0;
  const int axisv0ZRmV0RecZ    = 1;
  const int axisv0ZRmPtJetzGen = 2;
  const int axisv0ZRmV0Gen     = 3;
  const int nDimV0ZRm          = 4;

  const string nameDirJetRm  = "jet-fragmentation/matching/jets/";
  const string nameJetRm     = nameDirJetRm + "matchDetJetPtPartJetPt";
  const string nameJetFake   = nameDirJetRm + "fakeJetPtEtaPhi";
  const string nameJetMiss   = nameDirJetRm + "missJetPtEtaPhi";

  const string nameDirV0Rm   = nameDirJetRm + "V0/";
  const string nameV0PtRm    = nameDirV0Rm + "partJetPtV0PtDetJetPtV0Pt";
  const string nameV0PtFake  = nameDirV0Rm + "fakeJetPtV0PtEtaPhi";
  const string nameV0PtMiss  = nameDirV0Rm + "missJetPtV0PtEtaPhi";
  const string nameV0ZRm     = nameDirV0Rm + "matchDetJetPtV0TrackProjPartJetPtV0TrackProj";
  const string nameV0ZFake   = nameDirV0Rm + "fakeJetPtV0TrackProj";
  const string nameV0ZMiss   = nameDirV0Rm + "missJetPtV0TrackProj";

  const string nameK0SPtRm   = nameDirV0Rm + "partJetPtK0SPtDetJetPtK0SPt";
  const string nameK0SPtFake = nameDirV0Rm + "fakeJetPtK0SPtEtaPhi";
  const string nameK0SPtMiss = nameDirV0Rm + "missJetPtK0SPtEtaPhi";
  const string nameK0SZRm    = nameDirV0Rm + "matchDetJetPtK0STrackProjPartJetPtK0STrackProj";
  const string nameK0SZFake  = nameDirV0Rm + "fakeJetPtK0STrackProj";
  const string nameK0SZMiss  = nameDirV0Rm + "missJetPtK0STrackProj";
}

struct InputSettings {
  private:
    string getNameFromVar(string prefix, string varstring, string suffix);
    string getNameFromPtJet(string prefix, double low, double high, string suffix);
    string getNameFromPtV0(string prefix, double low, double high, string suffix);
    string getNameFromZV0(string prefix, double low, double high, string suffix);
    bool setVariable(double a, double b, double &x, double &y);
  public:
    // Unfolding settings
    int doSmoothing = 0, nIterations = 1;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovToy;

    // General settings
    int train;
    verbositylvls::Verbosity verbosity = verbositylvls::kWarnings;
    string inputFileName, outputFileName, responseFileName;
    string hadron;
    string rmHistName;

    double ptjetminGen, ptjetmaxGen, ptjetminRec, ptjetmaxRec;
    double ptv0minGen, ptv0maxGen, ptv0minRec, ptv0maxRec;
    double zv0minGen, zv0maxGen, zv0minRec, zv0maxRec;
    double etamin, etamax;
    double binwidthptjet, binwidthptv0, binwidthzv0;

    bool makeplots = false, logplot = false, ratioplot = false;

    string getNameFromPtJetGen(string prefix, string suffix) { return getNameFromPtJet(prefix, ptjetminGen, ptjetmaxGen, suffix); }
    string getNameFromPtJetRec(string prefix, string suffix) { return getNameFromPtJet(prefix, ptjetminRec, ptjetmaxRec, suffix); }
    string getNameFromPtV0Gen(string prefix, string suffix) { return getNameFromPtV0(prefix, ptv0minGen, ptv0maxGen, suffix); }
    string getNameFromPtV0Rec(string prefix, string suffix) { return getNameFromPtV0(prefix, ptv0minRec, ptv0maxRec, suffix); }
    string getNameFromZV0Gen(string prefix, string suffix) { return getNameFromZV0(prefix, zv0minGen, zv0maxGen, suffix); }
    string getNameFromZV0Rec(string prefix, string suffix) { return getNameFromZV0(prefix, zv0minRec, zv0maxRec, suffix); }
    bool isVarInRange(double var, double min, double max) { return (var >= min && var < max); }
    bool isVarInRange(double var, array<double, 2> range) { return isVarInRange(var, range[0], range[1]); }
    bool printLog(string message, verbositylvls::Verbosity verbThreshold);
    string setInputFileNameFromTrain();
    bool setEta(double a, double b) { return setVariable(a, b, etamin, etamax); }
    bool setPtJetGen(double a, double b) { return setVariable(a, b, ptjetminGen, ptjetmaxGen); }
    bool setPtJetRec(double a, double b) { return setVariable(a, b, ptjetminRec, ptjetmaxRec); }
    bool setPtV0Gen(double a, double b) { return setVariable(a, b, ptv0minGen, ptv0maxGen); }
    bool setPtV0Rec(double a, double b) { return setVariable(a, b, ptv0minRec, ptv0maxRec); }
    bool setV0ZGen(double a, double b) { return setVariable(a, b, zv0minGen, zv0maxGen); }
    bool setV0ZRec(double a, double b) { return setVariable(a, b, zv0minRec, zv0maxRec); }
    template <typename T> bool writeOutputToFile(T* obj);
    template <typename T> bool writeOutputsToFile(vector<T*> obj);
};

string InputSettings::getNameFromVar(string prefix, string varstring, string suffix = "") {
  return prefix + "_" + varstring + suffix;
}

string InputSettings::getNameFromPtJet(string prefix, double low, double high, string suffix = "") {
  return getNameFromVar(prefix, TString::Format("ptjet%.f-%.f", low, high).Data(), suffix);
}

string InputSettings::getNameFromPtV0(string prefix, double low, double high, string suffix = "") {
  return getNameFromVar(prefix, TString::Format("ptv0%.1f-%.1f", low, high).Data(), suffix);
}

string InputSettings::getNameFromZV0(string prefix, double low, double high, string suffix = "") {
  return getNameFromVar(prefix, TString::Format("zv0%.3f-%.3f", low, high).Data(), suffix);
}

bool InputSettings::printLog(string message, verbositylvls::Verbosity messageVerbLevel) {
  if (!verbositylvls::passVerbosityCheck(messageVerbLevel, verbosity))
    return false;

  cout << message << endl;
  return true;
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  inputFileName = s;
  return s;
}

bool InputSettings::setVariable(double a, double b, double &x, double &y) {
  if (a > b) {
    printLog("InputSettings::setVariable() Error: min > max", verbositylvls::kErrors);
    return false;
  }
  x = a;
  y = b;
  return true;
}

template <typename T>
bool InputSettings::writeOutputToFile(T* obj) {
  if (!obj)
    return false;

  TFile* file = TFile::Open(outputFileName.c_str(), "UPDATE");
  obj->Write(obj->GetName(), TObject::kOverwrite);
  file->Close();
  return true;
}

template <typename T>
bool InputSettings::writeOutputsToFile(vector<T*> obj) {
  if (obj.empty())
    return false;

  for (auto obj : obj)
    writeOutputToFile(obj);

  return true;
}

// ----------------------------------------------------------

void CreateResponseJets(InputSettings& inputs) {
  inputs.printLog("Opening input file: " + inputs.inputFileName, verbositylvls::kDebug);
  TFile* file = TFile::Open(inputs.inputFileName.c_str());
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbositylvls::kErrors);
    return;
  }
  inputs.printLog("Retrieving histogram: " + inputs.rmHistName, verbositylvls::kDebug);
  TH2D* responseMatrix = (TH2D*)file->Get(inputs.rmHistName.c_str());
  if (!responseMatrix) {
    inputs.printLog("Error: could not find " + inputs.rmHistName + " in file " + inputs.inputFileName, verbositylvls::kErrors);
    return;
  }
  responseMatrix->SetName("responseMatrixJets");

  inputs.printLog("Creating histograms for jets.", verbositylvls::kDebug);

  array<int, 2> ptjetRecBins = getProjectionBins(responseMatrix->GetYaxis(), inputs.ptjetminRec, inputs.ptjetmaxRec);
  array<int, 2> ptjetGenBins = getProjectionBins(responseMatrix->GetYaxis(), inputs.ptjetminGen, inputs.ptjetmaxGen);
  array<double, 2> ptjetRecBinEdges = getProjectionEdges(responseMatrix->GetXaxis(), ptjetRecBins);
  array<double, 2> ptjetGenBinEdges = getProjectionEdges(responseMatrix->GetYaxis(), ptjetGenBins);

  int nbinsPtRec = (int)((ptjetRecBinEdges[1] - ptjetRecBinEdges[0]) / inputs.binwidthptjet);
  int nbinsPtGen = (int)((ptjetGenBinEdges[1] - ptjetGenBinEdges[0]) / inputs.binwidthptjet);

  TH1D* hRec = new TH1D("hRecJets", "Rec;p_{T,jet}", nbinsPtRec, ptjetRecBinEdges[0], ptjetRecBinEdges[1]);
  TH1D* hGen = new TH1D("hGenJets", "Gen;p_{T,jet}", nbinsPtGen, ptjetGenBinEdges[0], ptjetGenBinEdges[1]);

  TH1D* hMiss   = (TH1D*)hGen->Clone("hMissJets");
  TH1D* hKinEff = (TH1D*)hGen->Clone("hKinEffJets");
  TH1D* hFake   = (TH1D*)hGen->Clone("hFakeJets");

  // Create the RooUnfoldResponse and fill it
  RooUnfoldResponse* response = new RooUnfoldResponse("responseJets", "RMJets");
  response->Setup(hRec, hGen);
  inputs.printLog("Filling response.", verbositylvls::kDebug);

  for (int xBin = 0; xBin <= responseMatrix->GetNbinsX(); xBin++) {
    for (int yBin = 0; yBin <= responseMatrix->GetNbinsY(); yBin++) {
      double binContent = responseMatrix->GetBinContent(xBin, yBin);
      double ptjetRec = responseMatrix->GetXaxis()->GetBinCenter(xBin);
      double ptjetGen = responseMatrix->GetYaxis()->GetBinCenter(yBin);

      bool isAcceptedRec = inputs.isVarInRange(ptjetRec, ptjetRecBinEdges);
      //(ptjetRec >= ptjetRecBinEdges[0] && ptjetRec < ptjetRecBinEdges[1]);
      bool isAcceptedGen = inputs.isVarInRange(ptjetGen, ptjetGenBinEdges);
      //(ptjetGen >= ptjetGenBinEdges[0] && ptjetGen < ptjetGenBinEdges[1]);

      if (isAcceptedRec)
        hRec->Fill(ptjetRec, binContent);
      if (isAcceptedGen)
        hGen->Fill(ptjetGen, binContent);

      if (isAcceptedRec && isAcceptedGen) {
        response->Fill(ptjetRec, ptjetGen, binContent);
        hKinEff->Fill(ptjetGen, binContent);
      } else if (!isAcceptedRec && isAcceptedGen) {
        response->Miss(ptjetGen, binContent);
        hMiss->Fill(ptjetGen, binContent);
      } else if (isAcceptedRec && !isAcceptedGen) {
        response->Fake(ptjetRec, binContent);
        hFake->Fill(ptjetRec, binContent);
      }
    }
  }

  inputs.printLog("Saving objects to output file: " + inputs.outputFileName, verbositylvls::kDebug);

  hKinEff->Divide(hGen);
  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(responseMatrix);
  inputs.writeOutputsToFile(std::vector<TH1D*>({hRec, hGen, hMiss, hFake, hKinEff}));
}

void CreateResponseV0Pt(InputSettings& inputs) {
  inputs.printLog("Opening input file: " + inputs.inputFileName, verbositylvls::kDebug);
  TFile* file = TFile::Open(inputs.inputFileName.c_str());
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbositylvls::kErrors);
    return;
  }
  inputs.printLog("Retrieving histogram: " + inputs.rmHistName, verbositylvls::kDebug);
  THnSparseD* responseMatrix = (THnSparseD*)file->Get(inputs.rmHistName.c_str());
  if (!responseMatrix) {
    inputs.printLog("Error: could not find " + inputs.rmHistName + " in file " + inputs.inputFileName, verbositylvls::kErrors);
    return;
  }
  responseMatrix->SetName("responseMatrixV0Pt");

  inputs.printLog("Creating histograms for V0s in jets.", verbositylvls::kDebug);

  array<int, 2> ptjetRecBins = getProjectionBins(responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetRec), inputs.ptjetminRec, inputs.ptjetmaxRec);
  array<int, 2> ptjetGenBins = getProjectionBins(responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetGen), inputs.ptjetminGen, inputs.ptjetmaxGen);
  array<int, 2> ptv0RecBins  = getProjectionBins(responseMatrix->GetAxis(rmutilities::axisv0PtRmV0RecPt), inputs.ptv0minRec, inputs.ptv0maxRec);
  array<int, 2> ptv0GenBins  = getProjectionBins(responseMatrix->GetAxis(rmutilities::axisv0PtRmV0GenPt), inputs.ptv0minGen, inputs.ptv0maxGen);
  array<double, 2> ptjetRecBinEdges = getProjectionEdges(responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetRec), ptjetRecBins);
  array<double, 2> ptjetGenBinEdges = getProjectionEdges(responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetGen), ptjetGenBins);
  array<double, 2> ptv0RecBinEdges  = getProjectionEdges(responseMatrix->GetAxis(rmutilities::axisv0PtRmV0RecPt), ptv0RecBins);
  array<double, 2> ptv0GenBinEdges  = getProjectionEdges(responseMatrix->GetAxis(rmutilities::axisv0PtRmV0GenPt), ptv0GenBins);

  int nbinsPtJetRec = (int)((ptjetRecBinEdges[1] - ptjetRecBinEdges[0]) / inputs.binwidthptjet);
  int nbinsPtJetGen = (int)((ptjetGenBinEdges[1] - ptjetGenBinEdges[0]) / inputs.binwidthptjet);
  int nbinsPtV0Rec  = (int)((ptv0RecBinEdges[1] - ptv0RecBinEdges[0]) / inputs.binwidthptv0);
  int nbinsPtV0Gen  = (int)((ptv0GenBinEdges[1] - ptv0GenBinEdges[0]) / inputs.binwidthptv0);

  TH2D* hRec = new TH2D("hRecV0Pt", "Rec;p_{T,V0};p_{T,jet}", nbinsPtV0Rec, ptv0RecBinEdges[0], ptv0RecBinEdges[1], nbinsPtJetRec, ptjetRecBinEdges[0], ptjetRecBinEdges[1]);
  TH2D* hGen = new TH2D("hGenV0Pt", "Gen;p_{T,V0};p_{T,jet}", nbinsPtV0Gen, ptv0GenBinEdges[0], ptv0GenBinEdges[1], nbinsPtJetGen, ptjetGenBinEdges[0], ptjetGenBinEdges[1]);

  TH2D* hMiss   = (TH2D*)hGen->Clone("hMissV0Pt");
  TH2D* hKinEff = (TH2D*)hGen->Clone("hKinEffV0Pt");
  TH2D* hFake   = (TH2D*)hGen->Clone("hFakeV0Pt");

  // Create the RooUnfoldResponse and fill it
  RooUnfoldResponse* response = new RooUnfoldResponse("responseV0Pt", "RMV0Pt");
  response->Setup(hRec, hGen);
  inputs.printLog("Filling response.", verbositylvls::kDebug);

  int* coord = new int[responseMatrix->GetNdimensions()]; //Carries the bin coordinates
  for (int iBin = 0; iBin <= responseMatrix->GetNbins(); iBin++) {
    double binContent = responseMatrix->GetBinContent(iBin, coord);
    double ptjetRec  = responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetRec)->GetBinCenter(coord[rmutilities::axisv0PtRmPtJetRec]);
    double ptjetGen  = responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetGen)->GetBinCenter(coord[rmutilities::axisv0PtRmPtJetGen]);
    double ptv0Rec   = responseMatrix->GetAxis(rmutilities::axisv0PtRmV0RecPt)->GetBinCenter(coord[rmutilities::axisv0PtRmV0RecPt]);
    double ptv0Gen   = responseMatrix->GetAxis(rmutilities::axisv0PtRmV0GenPt)->GetBinCenter(coord[rmutilities::axisv0PtRmV0GenPt]);

    bool isAcceptedRec = inputs.isVarInRange(ptjetRec, ptjetRecBinEdges) && inputs.isVarInRange(ptv0Rec, ptv0RecBinEdges);
    bool isAcceptedGen = inputs.isVarInRange(ptjetGen, ptjetGenBinEdges) && inputs.isVarInRange(ptv0Gen, ptv0GenBinEdges);

    if (isAcceptedRec)
      hRec->Fill(ptv0Rec, ptjetRec, binContent);
    if (isAcceptedGen)
      hGen->Fill(ptv0Gen, ptjetGen, binContent);

    if (isAcceptedRec && isAcceptedGen) {
      response->Fill(ptv0Rec, ptjetRec, ptv0Gen, ptjetGen, binContent);
      hKinEff->Fill(ptv0Gen, ptjetGen, binContent);
    } else if (!isAcceptedRec && isAcceptedGen) {
      response->Miss(ptv0Gen, ptjetGen, binContent);
      hMiss->Fill(ptv0Gen, ptjetGen, binContent);
    } else if (isAcceptedRec && !isAcceptedGen) {
      response->Fake(ptv0Rec, ptjetRec, binContent);
      hFake->Fill(ptv0Rec, ptjetRec, binContent);
    }
  }

  inputs.printLog("Saving objects to output file: " + inputs.outputFileName, verbositylvls::kDebug);

  hKinEff->Divide(hGen);
  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(responseMatrix);
  inputs.writeOutputsToFile(std::vector<TH2D*>({hRec, hGen, hMiss, hFake, hKinEff}));
}

// void CreateResponseV0Z(InputSettings& inputs) {
//   inputs.printLog("Opening input file: " + inputs.inputFileName, verbositylvls::kDebug);
//   TFile* file = TFile::Open(inputs.inputFileName.c_str());
//   if (!file) {
//     inputs.printLog("Error: could not open file " + inputs.inputFileName, verbositylvls::kErrors);
//     return;
//   }
//   inputs.printLog("Retrieving histogram: " + inputs.rmHistName, verbositylvls::kDebug);
//   THnSparseD* responseMatrix = (THnSparseD*)file->Get(inputs.rmHistName.c_str());
//   if (!responseMatrix) {
//     inputs.printLog("Error: could not find " + inputs.rmHistName + " in file " + inputs.inputFileName, verbositylvls::kErrors);
//     return;
//   }
//   responseMatrix->SetName("responseMatrixV0Z");

//   inputs.printLog("Creating histograms for V0s in jets.", verbositylvls::kDebug);

//   array<int, 2> ptjetRecBins = getProjectionBins(responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetRec), inputs.ptjetminRec, inputs.ptjetmaxRec);
//   array<int, 2> ptjetGenBins = getProjectionBins(responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetGen), inputs.ptjetminGen, inputs.ptjetmaxGen);
//   array<int, 2> ptv0RecBins  = getProjectionBins(responseMatrix->GetAxis(rmutilities::axisv0PtRmV0RecPt), inputs.ptv0minRec, inputs.ptv0maxRec);
//   array<int, 2> ptv0GenBins  = getProjectionBins(responseMatrix->GetAxis(rmutilities::axisv0PtRmV0GenPt), inputs.ptv0minGen, inputs.ptv0maxGen);
//   array<double, 2> ptjetRecBinEdges = getProjectionEdges(responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetRec), ptjetRecBins);
//   array<double, 2> ptjetGenBinEdges = getProjectionEdges(responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetGen), ptjetGenBins);
//   array<double, 2> ptv0RecBinEdges  = getProjectionEdges(responseMatrix->GetAxis(rmutilities::axisv0PtRmV0RecPt), ptv0RecBins);
//   array<double, 2> ptv0GenBinEdges  = getProjectionEdges(responseMatrix->GetAxis(rmutilities::axisv0PtRmV0GenPt), ptv0GenBins);

//   int nbinsPtJetRec = (int)((ptjetRecBinEdges[1] - ptjetRecBinEdges[0]) / inputs.binwidthptjet);
//   int nbinsPtJetGen = (int)((ptjetGenBinEdges[1] - ptjetGenBinEdges[0]) / inputs.binwidthptjet);
//   int nbinsPtV0Rec  = (int)((ptv0RecBinEdges[1] - ptv0RecBinEdges[0]) / inputs.binwidthptv0);
//   int nbinsPtV0Gen  = (int)((ptv0GenBinEdges[1] - ptv0GenBinEdges[0]) / inputs.binwidthptv0);

//   TH2D* hRec = new TH2D("hRecV0Pt", "Rec;p_{T,V0};p_{T,jet}", nbinsPtV0Rec, ptv0RecBinEdges[0], ptv0RecBinEdges[1], nbinsPtJetRec, ptjetRecBinEdges[0], ptjetRecBinEdges[1]);
//   TH2D* hGen = new TH2D("hGenV0Pt", "Gen;p_{T,V0};p_{T,jet}", nbinsPtV0Gen, ptv0GenBinEdges[0], ptv0GenBinEdges[1], nbinsPtJetGen, ptjetGenBinEdges[0], ptjetGenBinEdges[1]);

//   TH2D* hMiss   = (TH2D*)hGen->Clone("hMissV0Pt");
//   TH2D* hKinEff = (TH2D*)hGen->Clone("hKinEffV0Pt");
//   TH2D* hFake   = (TH2D*)hGen->Clone("hFakeV0Pt");

//   // Create the RooUnfoldResponse and fill it
//   RooUnfoldResponse* response = new RooUnfoldResponse("responseV0Pt", "RMV0Pt");
//   response->Setup(hRec, hGen);
//   inputs.printLog("Filling response.", verbositylvls::kDebug);

//   int* coord = new int[responseMatrix->GetNdimensions()]; //Carries the bin coordinates
//   for (int iBin = 0; iBin <= responseMatrix->GetNbins(); iBin++) {
//     double binContent = responseMatrix->GetBinContent(iBin, coord);
//     double ptjetRec  = responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetRec)->GetBinCenter(coord[rmutilities::axisv0PtRmPtJetRec]);
//     double ptjetGen  = responseMatrix->GetAxis(rmutilities::axisv0PtRmPtJetGen)->GetBinCenter(coord[rmutilities::axisv0PtRmPtJetGen]);
//     double ptv0Rec   = responseMatrix->GetAxis(rmutilities::axisv0PtRmV0RecPt)->GetBinCenter(coord[rmutilities::axisv0PtRmV0RecPt]);
//     double ptv0Gen   = responseMatrix->GetAxis(rmutilities::axisv0PtRmV0GenPt)->GetBinCenter(coord[rmutilities::axisv0PtRmV0GenPt]);

//     bool isAcceptedRec = inputs.isVarInRange(ptjetRec, ptjetRecBinEdges) && inputs.isVarInRange(ptv0Rec, ptv0RecBinEdges);
//     bool isAcceptedGen = inputs.isVarInRange(ptjetGen, ptjetGenBinEdges) && inputs.isVarInRange(ptv0Gen, ptv0GenBinEdges);

//     if (isAcceptedRec)
//       hRec->Fill(ptv0Rec, ptjetRec, binContent);
//     if (isAcceptedGen)
//       hGen->Fill(ptv0Gen, ptjetGen, binContent);

//     if (isAcceptedRec && isAcceptedGen) {
//       response->Fill(ptv0Rec, ptjetRec, ptv0Gen, ptjetGen, binContent);
//       hKinEff->Fill(ptv0Gen, ptjetGen, binContent);
//     } else if (!isAcceptedRec && isAcceptedGen) {
//       response->Miss(ptv0Gen, ptjetGen, binContent);
//       hMiss->Fill(ptv0Gen, ptjetGen, binContent);
//     } else if (isAcceptedRec && !isAcceptedGen) {
//       response->Fake(ptv0Rec, ptjetRec, binContent);
//       hFake->Fill(ptv0Rec, ptjetRec, binContent);
//     }
//   }

//   inputs.printLog("Saving objects to output file: " + inputs.outputFileName, verbositylvls::kDebug);

//   hKinEff->Divide(hGen);
//   inputs.writeOutputToFile(response);
//   inputs.writeOutputToFile(responseMatrix);
//   inputs.writeOutputsToFile(std::vector<TH2D*>({hRec, hGen, hMiss, hFake, hKinEff}));
// }

// ----------------------------------------------------------

void UnfoldingJets(InputSettings& inputs) {
  TFile* responseFile = TFile::Open(inputs.responseFileName.c_str(), "READ");
  TH1D* trainingRec = (TH1D*)responseFile->Get("hRecJets");
  TH1D* trainingFake = (TH1D*)responseFile->Get("hFakeJets");
  RooUnfoldResponse* response = (RooUnfoldResponse*)responseFile->Get("responseJets");
  inputs.printLog("Rec hist is of size " + to_string(trainingRec->GetNbinsX()), verbositylvls::kDebug);

  // TH1D* testFake; // Get this from test file!

  string ruBayesName  = "ruBayes" + to_string(inputs.nIterations);
  string ruBayesTitle = ruBayesName;
  RooUnfoldBayes ruBayes(response, trainingRec, inputs.nIterations, inputs.doSmoothing, ruBayesName.c_str(), ruBayesTitle.c_str());

  string unfoldedName = "unfoldedJets" + to_string(inputs.nIterations);
  TH1D* unfolded = (TH1D*)ruBayes.Hreco(inputs.errorTreatment);
  unfolded->SetName(unfoldedName.c_str());
  unfolded->SetTitle(unfoldedName.c_str());

  string refoldedName = "refoldedJets" + to_string(inputs.nIterations);
  TH1D* refolded = (TH1D*)response->ApplyToTruth(unfolded, refoldedName.c_str());
  // refolded->Add(testFake);
  inputs.writeOutputToFile(&ruBayes);
  inputs.writeOutputsToFile(std::vector<TH1D*>{unfolded, refolded});

  // FIXME: Can't get cov matrix to work
  // Covariance matrix
  // string covMatrixName = "covMatrixJets" + to_string(inputs.nIterations);
  // TH1D tmp(ruBayes.Ereco(inputs.errorTreatment));
  // TH1D* covMatrix = (TH1D*)tmp.Clone(covMatrixName.c_str());
  // inputs.printLog("Covariance matrix is of size " + to_string(covMatrix->GetNbinsX()), verbositylvls::kDebug);

  // Pearson coefficients
  // string pearsonName = "pearsonJets" + to_string(inputs.nIterations);
  // TH2D* pearson = (TH2D*)covMatrix->Clone(pearsonName.c_str());
  // pearson->Reset();

  // for (int xCovBin = 1; xCovBin <= covMatrix->GetNbinsX(); xCovBin++) {
  //   for (int yCovBin = 1; yCovBin <= covMatrix->GetNbinsY(); yCovBin++) {
  //     double cov = covMatrix->GetBinContent(xCovBin, yCovBin);
  //     double sigmaX = sqrt(covMatrix->GetBinContent(xCovBin, xCovBin));
  //     double sigmaY = sqrt(covMatrix->GetBinContent(yCovBin, yCovBin));
  //     double pearsonCoeff = 0.;
  //     if (sigmaX > 0. && sigmaY > 0.)
  //       pearsonCoeff = cov / (sigmaX * sigmaY);
  //     pearson->SetBinContent(xCovBin, yCovBin, pearsonCoeff);
  //   }
  // }
}

// void ClosureTestJets(InputSettings& inputs) {
  // TFile* responseFile = TFile::Open(inputs.responseFileName.c_str(), "READ");
  // if (!responseFile) {
  //   inputs.printLog("Error: could not open file " + inputs.responseFileName, verbositylvls::kErrors);
  //   return;
  // }
  // RooUnfoldResponse* response = (RooUnfoldResponse*)responseFile->Get("responseJets");
  // TH2D* trainingRM = (TH2D*)responseFile->Get("responseMatrixJets");
  // TH1D* trainingRec = (TH1D*)responseFile->Get("hRecJets");
  // TH1D* trainingGen = (TH1D*)responseFile->Get("hGenJets");
  // TH1D* trainingFake = (TH1D*)responseFile->Get("hFakeJets");
  // TH1D* trainingMiss = (TH1D*)responseFile->Get("hMissJets");
  // TH1D* trainingKinEff = (TH1D*)responseFile->Get("hKinEffJets");
  // if (!response)
  //   inputs.printLog("Error: could not find responseJets in file " + inputs.responseFileName, verbositylvls::kErrors);
  // if (!trainingRM)
  //   inputs.printLog("Error: could not find responseMatrixJets in file " + inputs.responseFileName, verbositylvls::kErrors);
  // if (!trainingRec)
  //   inputs.printLog("Error: could not find hRecJets in file " + inputs.responseFileName, verbositylvls::kErrors);
  // if (!trainingGen)
  //   inputs.printLog("Error: could not find hGenJets in file " + inputs.responseFileName, verbositylvls::kErrors);
  // if (!trainingFake)
  //   inputs.printLog("Error: could not find hFakeJets in file " + inputs.responseFileName, verbositylvls::kErrors);
  // if (!trainingMiss)
  //   inputs.printLog("Error: could not find hMissJets in file " + inputs.responseFileName, verbositylvls::kErrors);
  // if (!trainingKinEff)
  //   inputs.printLog("Error: could not find hKinEffJets in file " + inputs.responseFileName, verbositylvls::kErrors);

  // if (!response || !trainingRM || !trainingRec || !trainingGen || !trainingFake || !trainingMiss || !trainingKinEff)
  //   return;



  // // Data that the response is tested on. For trivial closure testFile = trainingFile
  // TFile* testFile = TFile::Open(inputs.inputFileName.c_str(), "READ");
  // if (!testFile) {
  //   inputs.printLog("Error: could not open file " + inputs.inputFileName, verbositylvls::kErrors);
  //   return;
  // }
  // TH2D* testRM = (TH2D*)testFile->Get(inputs.rmHistName.c_str());
  // if (!testRM) {
  //   inputs.printLog("Error: could not find " + inputs.rmHistName + " in file " + inputs.inputFileName, verbositylvls::kErrors);
  //   return;
  // }
  // TH1D* testRec = (TH1D*)trainingRec->Clone("testRecJets");
  // testRec->Reset();
  // TH1D* testGen = (TH1D*)trainingGen->Clone("testGenJets");
  // testGen->Reset();
  // TH1D* testFake = (TH1D*)trainingRec->Clone("testFakeJets");
  // testFake->Reset();
  // TH1D* testMiss = (TH1D*)trainingGen->Clone("testMissJets");
  // testMiss->Reset();
  // TH1D* testKinEff = (TH1D*)trainingGen->Clone("testKinEffJets");
  // testKinEff->Reset();
// }
// ----------------------------------------------------------
//
// Interface
//
// ----------------------------------------------------------

void createresponses() {
  InputSettings x; x.verbosity = verbositylvls::kDebug;
  x.train = 468659;
  x.inputFileName = to_string(x.train) + ".root";
  x.outputFileName = "RooUnfoldResponse" + to_string(x.train) + ".root";
  x.rmHistName = rmutilities::nameJetRm;
  x.binwidthptjet = 5.;
  x.setPtJetRec(10., 100.);
  x.setPtJetGen(10., 100.);

  CreateResponseJets(x);

  x.rmHistName = rmutilities::nameV0PtRm;
  x.binwidthptv0 = 1.;
  x.setPtV0Rec(1., x.ptjetmaxRec);
  x.setPtV0Gen(1., x.ptjetmaxGen);
  CreateResponseV0Pt(x);
}

void doclosuretest() {
  InputSettings x; x.verbosity = verbositylvls::kDebug;
  x.inputFileName = "468659.root";
  x.responseFileName = "RooUnfoldResponse.root";
  x.outputFileName = "ClosureTest.root";
  x.nIterations = 3;
  UnfoldingJets(x);
}

# endif