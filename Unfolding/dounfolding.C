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
  bool passVerbosityCheck(Verbosity v, Verbosity threshold) {
    return (is_valid(v) && v <= threshold);
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

  const string nameDirRm     = "jet-fragmentation/matching/jets/";
  const string nameJetRm     = nameDirRm + "matchDetJetPtPartJetPt";
  const string nameJetFake   = nameDirRm + "fakeJetPtEtaPhi";
  const string nameJetMiss   = nameDirRm + "missJetPtEtaPhi";
  const string nameV0PtRm    = nameDirRm + "V0/partJetPtV0PtDetJetPtV0Pt";
  const string nameV0PtFake  = nameDirRm + "V0/fakeJetPtV0PtEtaPhi";
  const string nameV0PtMiss  = nameDirRm + "V0/missJetPtV0PtEtaPhi";
  const string nameV0ZRm     = nameDirRm + "V0/matchDetJetPtV0TrackProjPartJetPtV0TrackProj";
  const string nameV0ZFake   = nameDirRm + "V0/fakeJetPtV0TrackProj";
  const string nameV0ZMiss   = nameDirRm + "V0/missJetPtV0TrackProj";
  const string nameK0SPtRm   = nameDirRm + "V0/partJetPtK0SPtDetJetPtK0SPt";
  const string nameK0SPtFake = nameDirRm + "V0/fakeJetPtK0SPtEtaPhi";
  const string nameK0SPtMiss = nameDirRm + "V0/missJetPtK0SPtEtaPhi";
  const string nameK0SZRm    = nameDirRm + "V0/matchDetJetPtK0STrackProjPartJetPtK0STrackProj";
  const string nameK0SZFake  = nameDirRm + "V0/fakeJetPtK0STrackProj";
  const string nameK0SZMiss  = nameDirRm + "V0/missJetPtK0STrackProj";
}

struct InputSettings {
  private:
    string getNameFromVar(string prefix, string varstring, string suffix);
    string getNameFromPtJet(string prefix, double low, double high, string suffix);
    string getNameFromPtV0(string prefix, double low, double high, string suffix);
    string getNameFromZV0(string prefix, double low, double high, string suffix);
  public:
    int train;
    verbositylvls::Verbosity verbosity = verbositylvls::kWarnings;
    string inputFileName, outputFileName;
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
    bool printLog(string message, verbositylvls::Verbosity verbThreshold);
    string setInputFileNameFromTrain();
    bool setEta(double a, double b);
    bool setPtJetGen(double a, double b);
    bool setPtJetRec(double a, double b);
    bool setPtV0Gen(double a, double b);
    bool setPtV0Rec(double a, double b);
    bool setV0ZGen(double a, double b);
    bool setV0ZRec(double a, double b);
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

bool InputSettings::printLog(string message, verbositylvls::Verbosity verbThreshold) {
  if (!verbositylvls::passVerbosityCheck(verbosity, verbThreshold))
    return false;

  cout << message << endl;
  return true;
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  inputFileName = s;
  return s;
}

bool InputSettings::setEta(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setEta() Error: min > max", verbositylvls::kErrors);
    return false;
  }
  etamin = a;
  etamax = b;
  return true;
}

bool InputSettings::setPtJetGen(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setPtJetGen() Error: min > max", verbositylvls::kErrors);
    return false;
  }
  ptjetminGen = a;
  ptjetmaxGen = b;
  return true;
}

bool InputSettings::setPtJetRec(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setPtJetRec() Error: min > max", verbositylvls::kErrors);
    return false;
  }
  ptjetminRec = a;
  ptjetmaxRec = b;
  return true;
}

bool InputSettings::setPtV0Gen(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setPtV0Gen() Error: min > max", verbositylvls::kErrors);
    return false;
  }
  ptv0minGen = a;
  ptv0maxGen = b;
  return true;
}

bool InputSettings::setPtV0Rec(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setPtV0Rec() Error: min > max", verbositylvls::kErrors);
    return false;
  }
  ptv0minRec = a;
  ptv0maxRec = b;
  return true;
}

bool InputSettings::setV0ZGen(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setV0ZGen() Error: min > max", verbositylvls::kErrors);
    return false;
  }
  zv0minGen = a;
  zv0maxGen = b;
  return true;
}

bool InputSettings::setV0ZRec(double a, double b) {
  if (a > b) {
    printLog("InputSettings::setV0ZRec() Error: min > max", verbositylvls::kErrors);
    return false;
  }
  zv0minRec = a;
  zv0maxRec = b;
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
  TFile* file = TFile::Open(inputs.inputFileName.c_str());
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbositylvls::kErrors);
    return;
  }
  TH2D* responseMatrix = (TH2D*)file->Get(inputs.rmHistName.c_str());
  if (!responseMatrix) {
    inputs.printLog("Error: could not find " + inputs.rmHistName + " in file " + inputs.inputFileName, verbositylvls::kErrors);
    return;
  }
  responseMatrix->SetName("responseMatrixJets");

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

  for (int xBin = 0; xBin <= responseMatrix->GetNbinsX(); xBin++) {
    for (int yBin = 0; yBin <= responseMatrix->GetNbinsY(); yBin++) {
      double binContent = responseMatrix->GetBinContent(xBin, yBin);
      double ptjetRec = responseMatrix->GetXaxis()->GetBinCenter(xBin);
      double ptjetGen = responseMatrix->GetYaxis()->GetBinCenter(yBin);

      bool isAcceptedRec = (ptjetRec >= ptjetRecBinEdges[0] && ptjetRec < ptjetRecBinEdges[1]);
      bool isAcceptedGen = (ptjetGen >= ptjetGenBinEdges[0] && ptjetGen < ptjetGenBinEdges[1]);

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

  hKinEff->Divide(hGen);
  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(responseMatrix);
  inputs.writeOutputsToFile(std::vector<TH1D*>({hRec, hGen, hMiss, hFake, hKinEff}));
}

void Run() {
  InputSettings x; x.verbosity = verbositylvls::kDebug;
  x.inputFileName = "496208.root";
  x.outputFileName = "RooUnfoldResponse.root";
  x.rmHistName = rmutilities::nameJetRm;
  x.binwidthptjet = 5.;
  x.setPtJetRec(10., 100.);
  x.setPtJetGen(10., 100.);

  CreateResponseJets(x);
  // Do Closure test
}

# endif