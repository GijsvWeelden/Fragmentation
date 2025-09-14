
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "../plotting_macros/histUtils.C"

#ifndef DO_UNFOLDING
#define DO_UNFOLDING

namespace verbosityutilities {
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
  bool printLog(string message, Verbosity messageVerbLevel, Verbosity threshold) {
    if (!passVerbosityCheck(messageVerbLevel, threshold))
      return false;

    cout << message << endl;
    return true;
  }
}

namespace rmutilities {
  // For loading from train output
  namespace analysis {
    const int axisJetRmPtJetRec = 0;
    const int axisJetRmPtJetGen = 1;
    const int ndimJetRm         = 2;

    const int axisv0PtRmPtJetGen = 0;
    const int axisv0PtRmV0Gen    = 1;
    const int axisv0PtRmPtJetRec = 2;
    const int axisv0PtRmV0Rec    = 3;
    const int nDimV0PtRm         = 4;

    const int axisv0ZRmPtJetRec  = 0;
    const int axisv0ZRmV0Rec     = 1;
    const int axisv0ZRmPtJetGen  = 2;
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

  // For saving/loading from response file
  namespace unfolding {
    const string nameRooUnfoldBayesJets = "ruBayesJets";
    const string nameUnfoldedJets       = "unfoldedJets";
    const string nameRefoldedJets       = "refoldedJets";
    const string nameCovMatrixJets      = "covMatrixJets";
    const string namePearsonJets        = "pearsonJets";

    const string nameRooUnfoldBayesV0Pt = "ruBayesV0Pt";
    const string nameUnfoldedV0Pt       = "unfoldedV0Pt";
    const string nameRefoldedV0Pt       = "refoldedV0Pt";
    const string nameCovMatrixV0Pt      = "covMatrixV0Pt";
    const string namePearsonV0Pt        = "pearsonV0Pt";

    // For setting up training and testing namespaces
    const string prefixTraining   = "training";
    const string prefixTest       = "test";
    const string nameResponseJets = "responseJets";
    const string nameRmJets       = "responseMatrixJets";
    const string nameRecJets      = "RecJets";
    const string nameGenJets      = "GenJets";
    const string nameMissJets     = "MissJets";
    const string nameFakeJets     = "FakeJets";
    const string nameKinEffJets   = "KinEffJets";

    const string nameResponseV0Pt = "responseV0Pt";
    const string nameRmV0Pt       = "responseMatrixV0Pt";
    const string nameRecV0Pt      = "RecV0Pt";
    const string nameGenV0Pt      = "GenV0Pt";
    const string nameMissV0Pt     = "MissV0Pt";
    const string nameFakeV0Pt     = "FakeV0Pt";
    const string nameKinEffV0Pt   = "KinEffV0Pt";
  }
  namespace training {
    const string nameResponseJets = unfolding::nameResponseJets;
    const string nameRmJets       = unfolding::nameRmJets;
    const string nameRecJets      = unfolding::prefixTraining + unfolding::nameRecJets;
    const string nameGenJets      = unfolding::prefixTraining + unfolding::nameGenJets;
    const string nameMissJets     = unfolding::prefixTraining + unfolding::nameMissJets;
    const string nameFakeJets     = unfolding::prefixTraining + unfolding::nameFakeJets;
    const string nameKinEffJets   = unfolding::prefixTraining + unfolding::nameKinEffJets;

    const string nameResponseV0Pt = unfolding::nameResponseV0Pt;
    const string nameRmV0Pt       = unfolding::nameRmV0Pt;
    const string nameRecV0Pt      = unfolding::prefixTraining + unfolding::nameRecV0Pt;
    const string nameGenV0Pt      = unfolding::prefixTraining + unfolding::nameGenV0Pt;
    const string nameMissV0Pt     = unfolding::prefixTraining + unfolding::nameMissV0Pt;
    const string nameFakeV0Pt     = unfolding::prefixTraining + unfolding::nameFakeV0Pt;
    const string nameKinEffV0Pt   = unfolding::prefixTraining + unfolding::nameKinEffV0Pt;
  }
  namespace testing {
    const string nameResponseJets = unfolding::nameResponseJets;
    const string nameRmJets       = unfolding::nameRmJets;
    const string nameRecJets      = unfolding::prefixTest + unfolding::nameRecJets;
    const string nameGenJets      = unfolding::prefixTest + unfolding::nameGenJets;
    const string nameMissJets     = unfolding::prefixTest + unfolding::nameMissJets;
    const string nameFakeJets     = unfolding::prefixTest + unfolding::nameFakeJets;
    const string nameKinEffJets   = unfolding::prefixTest + unfolding::nameKinEffJets;

    const string nameResponseV0Pt = unfolding::nameResponseV0Pt;
    const string nameRmV0Pt       = unfolding::nameRmV0Pt;
    const string nameRecV0Pt      = unfolding::prefixTest + unfolding::nameRecV0Pt;
    const string nameGenV0Pt      = unfolding::prefixTest + unfolding::nameGenV0Pt;
    const string nameMissV0Pt     = unfolding::prefixTest + unfolding::nameMissV0Pt;
    const string nameFakeV0Pt     = unfolding::prefixTest + unfolding::nameFakeV0Pt;
    const string nameKinEffV0Pt   = unfolding::prefixTest + unfolding::nameKinEffV0Pt;
  }
}

struct InputSettings {
  private:
    string getNameFromVar(string prefix, string varstring, string suffix) {
      return prefix + "_" + varstring + suffix;
    }
    string getNameFromPtJet(string prefix, double low, double high, string suffix) {
      return getNameFromVar(prefix, TString::Format("ptjet%.f-%.f", low, high).Data(), suffix);
    }
    string getNameFromPtV0(string prefix, double low, double high, string suffix) {
      return getNameFromVar(prefix, TString::Format("ptv0%.1f-%.1f", low, high).Data(), suffix);
    }
    string getNameFromZV0(string prefix, double low, double high, string suffix) {
      return getNameFromVar(prefix, TString::Format("zv0%.3f-%.3f", low, high).Data(), suffix);
    }
    template <typename T> bool setVariable(T a, T b, T &x, T &y);
  public:
    // Unfolding settings
    int doSmoothing = 0, nIterations = 1;
    int minIteration = 1, maxIteration = 1;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovToy;

    // General settings
    int train;
    verbosityutilities::Verbosity verbosity = verbosityutilities::kWarnings;
    string inputFileName, outputFileName, responseFileName;
    string hadron;
    string rmHistName;

    double ptjetminGen, ptjetmaxGen, ptjetminRec, ptjetmaxRec;
    double ptv0minGen, ptv0maxGen, ptv0minRec, ptv0maxRec;
    double zv0minGen, zv0maxGen, zv0minRec, zv0maxRec;
    double etamin, etamax;
    double binwidthptjet, binwidthptv0, binwidthzv0;

    bool makeplots = false, logplot = false, ratioplot = false;
    bool doTrivialClosureTest = true;

    // Methods using private methods
    string getNameFromPtJetGen(string prefix, string suffix) {
      return getNameFromPtJet(prefix, ptjetminGen, ptjetmaxGen, suffix);
    }
    string getNameFromPtJetRec(string prefix, string suffix) {
      return getNameFromPtJet(prefix, ptjetminRec, ptjetmaxRec, suffix);
    }
    string getNameFromPtV0Gen(string prefix, string suffix) {
      return getNameFromPtV0(prefix, ptv0minGen, ptv0maxGen, suffix);
    }
    string getNameFromPtV0Rec(string prefix, string suffix) {
      return getNameFromPtV0(prefix, ptv0minRec, ptv0maxRec, suffix);
    }
    string getNameFromZV0Gen(string prefix, string suffix) {
      return getNameFromZV0(prefix, zv0minGen, zv0maxGen, suffix);
    }
    string getNameFromZV0Rec(string prefix, string suffix) {
      return getNameFromZV0(prefix, zv0minRec, zv0maxRec, suffix);
    }
    bool setEta(double a, double b) {
      return setVariable(a, b, etamin, etamax);
    }
    bool setIterations(int a, int b) {
      return setVariable(a, b, minIteration, maxIteration);
    }
    bool setPtJetGen(double a, double b) {
      return setVariable(a, b, ptjetminGen, ptjetmaxGen);
    }
    bool setPtJetRec(double a, double b) {
      return setVariable(a, b, ptjetminRec, ptjetmaxRec);
    }
    bool setPtV0Gen(double a, double b) {
      return setVariable(a, b, ptv0minGen, ptv0maxGen);
    }
    bool setPtV0Rec(double a, double b) {
      return setVariable(a, b, ptv0minRec, ptv0maxRec);
    }
    bool setV0ZGen(double a, double b) {
      return setVariable(a, b, zv0minGen, zv0maxGen);
    }
    bool setV0ZRec(double a, double b) {
      return setVariable(a, b, zv0minRec, zv0maxRec);
    }

    // Utilities
    bool isHistInRange(TH1* hist, double min, double max);
    bool isHistConsistentWithZero(TH1* hist, double threshold = 1e-10) {
      return isHistInRange(hist, -threshold, threshold);
    }
    bool isHistConsistentWithOne(TH1* hist, double threshold = 1e-10) {
      return isHistInRange(hist, 1. - threshold, 1. + threshold);
    }
    template <typename T> bool isVarInRange(T var, T min, T max, double epsilon = 1e-5) {
      return ((var >= min - epsilon) && (var < max + epsilon));
    }
    bool isVarInRange(double var, array<double, 2> range) {
      return isVarInRange(var, range[0], range[1]);
    }
    bool isVarConsistentWithZero(double var, double threshold = 1e-10) {
      return isVarInRange(var, -threshold, threshold);
    }
    bool printLog(string message, verbosityutilities::Verbosity messageVerbLevel) {
      return verbosityutilities::printLog(message, messageVerbLevel, verbosity);
    }
    string setInputFileNameFromTrain();
    template <typename T> bool writeOutputToFile(T* obj);
    template <typename T> bool writeOutputsToFile(vector<T*> obj);
};

bool InputSettings::isHistInRange(TH1* hist, double min, double max) {
  if (!hist) {
    printLog("InputSettings::isHistInRange() Error: histogram is null", verbosityutilities::kErrors);
    return false;
  }
  bool histInRange = true;
  for (int iBin = 1; iBin < hist->GetNbinsX(); iBin++) {
    double binContent = hist->GetBinContent(iBin);
    if (!isVarInRange(binContent, min, max))
      histInRange = false;
  }
  return histInRange;
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
  inputFileName = s;
  return s;
}

template <typename T>
bool InputSettings::setVariable(T a, T b, T &x, T &y) {
  if (a > b) {
    printLog("InputSettings::setVariable() Error: min > max", verbosityutilities::kErrors);
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
bool InputSettings::writeOutputsToFile(vector<T*> objs) {
  if (objs.empty())
    return false;

  for (auto obj : objs)
    writeOutputToFile(obj);

  return true;
}

// ----------------------------------------------------------

// NB: V0 Pt RM is filled (gen, rec), but the jet RM and V0 Z RM are filled (rec, gen)!
// NB: The RooUnfoldResponse is filled (rec, gen) always!
void FillFromRmJets(InputSettings& inputs, TH2D* responseMatrix, TH1D* hRec, TH1D* hGen, TH1D* hMiss, TH1D* hKinEff, TH1D* hFake, RooUnfoldResponse* response = nullptr) {
  // Fill the distributions from the response matrix
  // If the RooUnfoldResponse is given, fill it as well
  array<double, 2> ptjetRecBinEdges = {hRec->GetXaxis()->GetXmin(), hRec->GetXaxis()->GetXmax()};
  array<double, 2> ptjetGenBinEdges = {hGen->GetXaxis()->GetXmin(), hGen->GetXaxis()->GetXmax()};

  for (int xBin = 0; xBin <= responseMatrix->GetNbinsX(); xBin++) {
    for (int yBin = 0; yBin <= responseMatrix->GetNbinsY(); yBin++) {
      double binContent = responseMatrix->GetBinContent(xBin, yBin);
      double ptjetRec = responseMatrix->GetXaxis()->GetBinCenter(xBin);
      double ptjetGen = responseMatrix->GetYaxis()->GetBinCenter(yBin);

      bool isAcceptedRec = inputs.isVarInRange(ptjetRec, ptjetRecBinEdges);
      bool isAcceptedGen = inputs.isVarInRange(ptjetGen, ptjetGenBinEdges);

      if (isAcceptedRec)
        hRec->Fill(ptjetRec, binContent);
      if (isAcceptedGen)
        hGen->Fill(ptjetGen, binContent);

      if (isAcceptedRec && isAcceptedGen) {
        hKinEff->Fill(ptjetGen, binContent);
        if (response) response->Fill(ptjetRec, ptjetGen, binContent);
      } else if (!isAcceptedRec && isAcceptedGen) {
        hMiss->Fill(ptjetGen, binContent);
        if (response) response->Miss(ptjetGen, binContent);
      } else if (isAcceptedRec && !isAcceptedGen) {
        hFake->Fill(ptjetRec, binContent);
        if (response) response->Fake(ptjetRec, binContent);
      }
    }
  }
  hKinEff->Divide(hGen);
}

// NB: V0 Pt RM is filled (gen, rec), but the jet RM and V0 Z RM are filled (rec, gen)!
// NB: The RooUnfoldResponse is filled (rec, gen) always!
void FillFromRmV0Pt(InputSettings& inputs, THnSparseD* responseMatrix, TH2D* hRec, TH2D* hGen, TH2D* hMiss, TH2D* hKinEff, TH2D* hFake, RooUnfoldResponse* response = nullptr) {
  // Fill the distributions from the response matrix
  // If the RooUnfoldResponse is given, fill it as well
  array<double, 2> ptv0GenBinEdges  = {hGen->GetXaxis()->GetXmin(), hGen->GetXaxis()->GetXmax()};
  array<double, 2> ptjetGenBinEdges = {hGen->GetYaxis()->GetXmin(), hGen->GetYaxis()->GetXmax()};
  array<double, 2> ptv0RecBinEdges  = {hRec->GetXaxis()->GetXmin(), hRec->GetXaxis()->GetXmax()};
  array<double, 2> ptjetRecBinEdges = {hRec->GetYaxis()->GetXmin(), hRec->GetYaxis()->GetXmax()};

  int* coord = new int[responseMatrix->GetNdimensions()]; //Carries the bin coordinates
  for (int iBin = 0; iBin <= responseMatrix->GetNbins(); iBin++) {
    double binContent = responseMatrix->GetBinContent(iBin, coord);
    double ptjetGen  = responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmPtJetGen)->GetBinCenter(coord[rmutilities::analysis::axisv0PtRmPtJetGen]);
    double ptv0Gen   = responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmV0Gen)->GetBinCenter(coord[rmutilities::analysis::axisv0PtRmV0Gen]);
    double ptjetRec  = responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmPtJetRec)->GetBinCenter(coord[rmutilities::analysis::axisv0PtRmPtJetRec]);
    double ptv0Rec   = responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmV0Rec)->GetBinCenter(coord[rmutilities::analysis::axisv0PtRmV0Rec]);

    bool isAcceptedGen = inputs.isVarInRange(ptjetGen, ptjetGenBinEdges) && inputs.isVarInRange(ptv0Gen, ptv0GenBinEdges);
    bool isAcceptedRec = inputs.isVarInRange(ptjetRec, ptjetRecBinEdges) && inputs.isVarInRange(ptv0Rec, ptv0RecBinEdges);

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
  hKinEff->Divide(hGen);
}

// NB: V0 Pt RM is filled (gen, rec), but the jet RM and V0 Z RM are filled (rec, gen)!
// NB: The RooUnfoldResponse is filled (rec, gen) always!
void FillFromRmV0Z(InputSettings& inputs, THnSparseD* responseMatrix, TH2D* hRec, TH2D* hGen, TH2D* hMiss, TH2D* hKinEff, TH2D* hFake, RooUnfoldResponse* response = nullptr) {
  // Fill the distributions from the response matrix
  // If the RooUnfoldResponse is given, fill it as well
  array<double, 2> zv0GenBinEdges  = {hGen->GetXaxis()->GetXmin(), hGen->GetXaxis()->GetXmax()};
  array<double, 2> ptjetGenBinEdges = {hGen->GetYaxis()->GetXmin(), hGen->GetYaxis()->GetXmax()};
  array<double, 2> zv0RecBinEdges  = {hRec->GetXaxis()->GetXmin(), hRec->GetXaxis()->GetXmax()};
  array<double, 2> ptjetRecBinEdges = {hRec->GetYaxis()->GetXmin(), hRec->GetYaxis()->GetXmax()};

  int* coord = new int[responseMatrix->GetNdimensions()]; //Carries the bin coordinates
  for (int iBin = 0; iBin <= responseMatrix->GetNbins(); iBin++) {
    double binContent = responseMatrix->GetBinContent(iBin, coord);
    double ptjetGen  = responseMatrix->GetAxis(rmutilities::analysis::axisv0ZRmPtJetGen)->GetBinCenter(coord[rmutilities::analysis::axisv0ZRmPtJetGen]);
    double zv0Gen   = responseMatrix->GetAxis(rmutilities::analysis::axisv0ZRmV0Gen)->GetBinCenter(coord[rmutilities::analysis::axisv0ZRmV0Gen]);
    double ptjetRec  = responseMatrix->GetAxis(rmutilities::analysis::axisv0ZRmPtJetRec)->GetBinCenter(coord[rmutilities::analysis::axisv0ZRmPtJetRec]);
    double zv0Rec   = responseMatrix->GetAxis(rmutilities::analysis::axisv0ZRmV0Rec)->GetBinCenter(coord[rmutilities::analysis::axisv0ZRmV0Rec]);

    bool isAcceptedGen = inputs.isVarInRange(ptjetGen, ptjetGenBinEdges) && inputs.isVarInRange(zv0Gen, zv0GenBinEdges);
    bool isAcceptedRec = inputs.isVarInRange(ptjetRec, ptjetRecBinEdges) && inputs.isVarInRange(zv0Rec, zv0RecBinEdges);

    if (isAcceptedRec)
      hRec->Fill(zv0Rec, ptjetRec, binContent);
    if (isAcceptedGen)
      hGen->Fill(zv0Gen, ptjetGen, binContent);

    if (isAcceptedRec && isAcceptedGen) {
      response->Fill(zv0Rec, ptjetRec, zv0Gen, ptjetGen, binContent);
      hKinEff->Fill(zv0Gen, ptjetGen, binContent);
    } else if (!isAcceptedRec && isAcceptedGen) {
      response->Miss(zv0Gen, ptjetGen, binContent);
      hMiss->Fill(zv0Gen, ptjetGen, binContent);
    } else if (isAcceptedRec && !isAcceptedGen) {
      response->Fake(zv0Rec, ptjetRec, binContent);
      hFake->Fill(zv0Rec, ptjetRec, binContent);
    }
  }
  hKinEff->Divide(hGen);
}

// ----------------------------------------------------------

void FillPearsonJets(InputSettings& inputs, TH2D* covMatrix, TH2D* pearson) {
  pearson->Reset();
  for (int xCovBin = 1; xCovBin <= covMatrix->GetNbinsX(); xCovBin++) {
    for (int yCovBin = 1; yCovBin <= covMatrix->GetNbinsY(); yCovBin++) {
      double cov = covMatrix->GetBinContent(xCovBin, yCovBin);
      double sigmaX = sqrt(covMatrix->GetBinContent(xCovBin, xCovBin));
      double sigmaY = sqrt(covMatrix->GetBinContent(yCovBin, yCovBin));
      double pearsonCoeff = 0.;
      if (std::isnan(cov) || std::isnan(sigmaX) || std::isnan(sigmaY))
        continue;
      if (inputs.isVarConsistentWithZero(sigmaX) || inputs.isVarConsistentWithZero(sigmaY))
        continue;

      pearsonCoeff = cov / std::sqrt(sigmaX * sigmaY);
      pearson->SetBinContent(xCovBin, yCovBin, pearsonCoeff);
    }
  }
}

// ----------------------------------------------------------

void CreateResponseJets(InputSettings& inputs) {
  inputs.printLog("Opening input file: " + inputs.inputFileName, verbosityutilities::kDebug);
  TFile* file = TFile::Open(inputs.inputFileName.c_str());
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }
  inputs.printLog("Retrieving histogram: " + inputs.rmHistName, verbosityutilities::kDebug);
  TH2D* responseMatrix = (TH2D*)file->Get(inputs.rmHistName.c_str());
  if (!responseMatrix) {
    inputs.printLog("Error: could not find " + inputs.rmHistName + " in file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }
  responseMatrix->SetName(rmutilities::unfolding::nameRmJets.c_str());
  responseMatrix->SetTitle(rmutilities::unfolding::nameRmJets.c_str());

  inputs.printLog("Creating histograms for jets.", verbosityutilities::kDebug);
  array<int, 2> ptjetRecBins = getProjectionBins(responseMatrix->GetYaxis(), inputs.ptjetminRec, inputs.ptjetmaxRec);
  array<int, 2> ptjetGenBins = getProjectionBins(responseMatrix->GetYaxis(), inputs.ptjetminGen, inputs.ptjetmaxGen);
  array<double, 2> ptjetRecBinEdges = getProjectionEdges(responseMatrix->GetXaxis(), ptjetRecBins);
  array<double, 2> ptjetGenBinEdges = getProjectionEdges(responseMatrix->GetYaxis(), ptjetGenBins);

  int nbinsPtRec = (int)((ptjetRecBinEdges[1] - ptjetRecBinEdges[0]) / inputs.binwidthptjet);
  int nbinsPtGen = (int)((ptjetGenBinEdges[1] - ptjetGenBinEdges[0]) / inputs.binwidthptjet);

  TH1D* hRec = new TH1D(rmutilities::unfolding::nameRecJets.c_str(), "Rec;p_{T,jet}", nbinsPtRec, ptjetRecBinEdges[0], ptjetRecBinEdges[1]);
  TH1D* hGen = new TH1D(rmutilities::unfolding::nameGenJets.c_str(), "Gen;p_{T,jet}", nbinsPtGen, ptjetGenBinEdges[0], ptjetGenBinEdges[1]);

  TH1D* hMiss   = (TH1D*)hGen->Clone(rmutilities::unfolding::nameMissJets.c_str());
  TH1D* hKinEff = (TH1D*)hGen->Clone(rmutilities::unfolding::nameKinEffJets.c_str());
  TH1D* hFake   = (TH1D*)hRec->Clone(rmutilities::unfolding::nameFakeJets.c_str());

  inputs.printLog("Filling response.", verbosityutilities::kDebug);
  RooUnfoldResponse* response = new RooUnfoldResponse(rmutilities::unfolding::nameResponseJets.c_str(), rmutilities::unfolding::nameResponseJets.c_str());
  response->Setup(hRec, hGen);
  FillFromRmJets(inputs, responseMatrix, hRec, hGen, hMiss, hKinEff, hFake, response);

  inputs.printLog("Saving objects to output file: " + inputs.outputFileName, verbosityutilities::kDebug);
  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(responseMatrix);
  inputs.writeOutputsToFile(std::vector<TH1D*>({hRec, hGen, hMiss, hFake, hKinEff}));
}

// NB: V0 Pt RM is filled (gen, rec), but the jet RM and V0 Z RM are filled (rec, gen)!
// NB: The RooUnfoldResponse is filled (rec, gen) always!
void CreateResponseV0Pt(InputSettings& inputs) {
  inputs.printLog("Opening input file: " + inputs.inputFileName, verbosityutilities::kDebug);
  TFile* file = TFile::Open(inputs.inputFileName.c_str());
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }
  inputs.printLog("Retrieving histogram: " + inputs.rmHistName, verbosityutilities::kDebug);
  THnSparseD* responseMatrix = (THnSparseD*)file->Get(inputs.rmHistName.c_str());
  if (!responseMatrix) {
    inputs.printLog("Error: could not find " + inputs.rmHistName + " in file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }
  responseMatrix->SetName(rmutilities::unfolding::nameRmV0Pt.c_str());
  responseMatrix->SetTitle(rmutilities::unfolding::nameRmV0Pt.c_str());

  inputs.printLog("Creating histograms for V0 Pt in jets.", verbosityutilities::kDebug);
  array<int, 2> ptjetRecBins = getProjectionBins(responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmPtJetRec), inputs.ptjetminRec, inputs.ptjetmaxRec);
  array<int, 2> ptjetGenBins = getProjectionBins(responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmPtJetGen), inputs.ptjetminGen, inputs.ptjetmaxGen);
  array<int, 2> ptv0RecBins  = getProjectionBins(responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmV0Rec), inputs.ptv0minRec, inputs.ptv0maxRec);
  array<int, 2> ptv0GenBins  = getProjectionBins(responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmV0Gen), inputs.ptv0minGen, inputs.ptv0maxGen);
  array<double, 2> ptjetRecBinEdges = getProjectionEdges(responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmPtJetRec), ptjetRecBins);
  array<double, 2> ptjetGenBinEdges = getProjectionEdges(responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmPtJetGen), ptjetGenBins);
  array<double, 2> ptv0RecBinEdges  = getProjectionEdges(responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmV0Rec), ptv0RecBins);
  array<double, 2> ptv0GenBinEdges  = getProjectionEdges(responseMatrix->GetAxis(rmutilities::analysis::axisv0PtRmV0Gen), ptv0GenBins);

  int nbinsPtJetRec = (int)((ptjetRecBinEdges[1] - ptjetRecBinEdges[0]) / inputs.binwidthptjet);
  int nbinsPtJetGen = (int)((ptjetGenBinEdges[1] - ptjetGenBinEdges[0]) / inputs.binwidthptjet);
  int nbinsPtV0Rec  = (int)((ptv0RecBinEdges[1] - ptv0RecBinEdges[0]) / inputs.binwidthptv0);
  int nbinsPtV0Gen  = (int)((ptv0GenBinEdges[1] - ptv0GenBinEdges[0]) / inputs.binwidthptv0);

  TH2D* hRec = new TH2D(rmutilities::unfolding::nameRecV0Pt.c_str(), "Rec;p_{T,V0};p_{T,jet}", nbinsPtV0Rec, ptv0RecBinEdges[0], ptv0RecBinEdges[1], nbinsPtJetRec, ptjetRecBinEdges[0], ptjetRecBinEdges[1]);
  TH2D* hGen = new TH2D(rmutilities::unfolding::nameGenV0Pt.c_str(), "Gen;p_{T,V0};p_{T,jet}", nbinsPtV0Gen, ptv0GenBinEdges[0], ptv0GenBinEdges[1], nbinsPtJetGen, ptjetGenBinEdges[0], ptjetGenBinEdges[1]);

  TH2D* hMiss   = (TH2D*)hGen->Clone(rmutilities::unfolding::nameMissV0Pt.c_str());
  TH2D* hKinEff = (TH2D*)hGen->Clone(rmutilities::unfolding::nameKinEffV0Pt.c_str());
  TH2D* hFake   = (TH2D*)hGen->Clone(rmutilities::unfolding::nameFakeV0Pt.c_str());

  inputs.printLog("Filling response.", verbosityutilities::kDebug);
  RooUnfoldResponse* response = new RooUnfoldResponse(rmutilities::unfolding::nameResponseV0Pt.c_str(), rmutilities::unfolding::nameResponseV0Pt.c_str());
  response->Setup(hRec, hGen);
  FillFromRmV0Pt(inputs, responseMatrix, hRec, hGen, hMiss, hKinEff, hFake, response);

  inputs.printLog("Saving objects to output file: " + inputs.outputFileName, verbosityutilities::kDebug);
  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(responseMatrix);
  inputs.writeOutputsToFile(std::vector<TH2D*>({hRec, hGen, hMiss, hFake, hKinEff}));
}

// ----------------------------------------------------------

void DoUnfoldingJets(InputSettings& inputs, int nIterations = -1) {
  if (nIterations <= 0)
    nIterations = inputs.nIterations;

  inputs.printLog("Doing unfolding with " + to_string(nIterations) + " iterations.", verbosityutilities::kInfo);
  inputs.printLog("Getting response from file: " + inputs.responseFileName + "\nGetting test distributions from file: " + inputs.inputFileName, verbosityutilities::kDebug);
  TFile* responseFile = TFile::Open(inputs.responseFileName.c_str(), "UPDATE");
  if (!responseFile) {
    inputs.printLog("Error: could not open file " + inputs.responseFileName, verbosityutilities::kErrors);
    return;
  }

  inputs.printLog("Retrieving response and histograms from training.", verbosityutilities::kDebug);
  TH1D* trainingRec = (TH1D*)responseFile->Get(rmutilities::unfolding::nameRecJets.c_str());
  TH1D* trainingGen = (TH1D*)responseFile->Get(rmutilities::unfolding::nameGenJets.c_str());
  TH1D* trainingFake = (TH1D*)responseFile->Get(rmutilities::unfolding::nameFakeJets.c_str());
  TH1D* trainingMiss = (TH1D*)responseFile->Get(rmutilities::unfolding::nameMissJets.c_str());
  TH1D* trainingKinEff = (TH1D*)responseFile->Get(rmutilities::unfolding::nameKinEffJets.c_str());
  RooUnfoldResponse* response = (RooUnfoldResponse*)responseFile->Get(rmutilities::unfolding::nameResponseJets.c_str());
  if (!trainingRec)
  inputs.printLog("Error: could not find " + rmutilities::unfolding::nameRecJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingGen)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameGenJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingFake)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameFakeJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingMiss)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameMissJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingKinEff)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameKinEffJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!response)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameResponseJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingRec || !trainingGen || !trainingFake || !trainingMiss || !trainingKinEff || !response)
    return;

  // Change the names to distinguish these from the test histograms when writing to file
  trainingRec->SetName(rmutilities::training::nameRecJets.c_str());
  trainingGen->SetName(rmutilities::training::nameGenJets.c_str());
  trainingFake->SetName(rmutilities::training::nameFakeJets.c_str());
  trainingMiss->SetName(rmutilities::training::nameMissJets.c_str());
  trainingKinEff->SetName(rmutilities::training::nameKinEffJets.c_str());

  inputs.printLog("Retrieving test histograms from file: " + inputs.inputFileName, verbosityutilities::kDebug);
  // Create the test histograms. These must have the same binning as the training histograms
  // For trivial closure test, they are copies, otherwise they are filled with independent data
  TH1D* testRec = (TH1D*)trainingRec->Clone(rmutilities::testing::nameRecJets.c_str());
  TH1D* testGen = (TH1D*)trainingGen->Clone(rmutilities::testing::nameGenJets.c_str());
  TH1D* testFake = (TH1D*)trainingFake->Clone(rmutilities::testing::nameFakeJets.c_str());
  TH1D* testMiss = (TH1D*)trainingMiss->Clone(rmutilities::testing::nameMissJets.c_str());
  TH1D* testKinEff = (TH1D*)trainingKinEff->Clone(rmutilities::testing::nameKinEffJets.c_str());

  if (inputs.doTrivialClosureTest) {
    inputs.printLog("Trivial closure test: using training histograms as test histograms.", verbosityutilities::kInfo);
  } else {
    inputs.printLog("Statistically independent closure test: getting test histograms from test response matrix.", verbosityutilities::kInfo);
    testRec->Reset();
    testGen->Reset();
    testFake->Reset();
    testMiss->Reset();
    testKinEff->Reset();

    TFile* testFile = TFile::Open(inputs.inputFileName.c_str(), "READ");
    if (!testFile) {
      inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
      return;
    }
    TH2D* testResponseMatrix = (TH2D*)testFile->Get(inputs.rmHistName.c_str());
    FillFromRmJets(inputs, testResponseMatrix, testRec, testGen, testMiss, testKinEff, testFake);
    testFile->Close();
  }

  inputs.printLog("Creating RooUnfoldBayes object and unfolding.", verbosityutilities::kDebug);
  string ruBayesName  = rmutilities::unfolding::nameRooUnfoldBayesJets + to_string(nIterations);
  string ruBayesTitle = ruBayesName;
  RooUnfoldBayes ruBayes(response, trainingRec, nIterations, inputs.doSmoothing, ruBayesName.c_str(), ruBayesTitle.c_str());

  string unfoldedName = rmutilities::unfolding::nameUnfoldedJets + to_string(nIterations);
  TH1D* unfolded = (TH1D*)ruBayes.Hreco(inputs.errorTreatment);
  unfolded->SetName(unfoldedName.c_str());
  unfolded->SetTitle(unfoldedName.c_str());

  string refoldedName = rmutilities::unfolding::nameRefoldedJets + to_string(nIterations);
  TH1D* refolded = (TH1D*)response->ApplyToTruth(unfolded, refoldedName.c_str());
  refolded->Add(testFake);

  inputs.printLog("Calculating covariance matrix and Pearson coefficients.", verbosityutilities::kDebug);
  string covMatrixName = rmutilities::unfolding::nameCovMatrixJets + to_string(nIterations);
  TH2D tmp(ruBayes.Ereco(inputs.errorTreatment));
  TH2D* covMatrix = (TH2D*)tmp.Clone(covMatrixName.c_str());
  covMatrix->SetTitle(covMatrixName.c_str());

  inputs.printLog("Rec hist is of size " + to_string(trainingRec->GetNbinsX()) + ", covariance matrix is of size " + to_string(covMatrix->GetNbinsX()) + " x " + to_string(covMatrix->GetNbinsY()), verbosityutilities::kDebug);

  string pearsonName = rmutilities::unfolding::namePearsonJets + to_string(nIterations);
  TH2D* pearson = (TH2D*)covMatrix->Clone(pearsonName.c_str());
  pearson->SetTitle(pearsonName.c_str());
  FillPearsonJets(inputs, covMatrix, pearson);

  inputs.printLog("Unfolding done. Writing outputs to file " + inputs.outputFileName, verbosityutilities::kInfo);
  inputs.writeOutputToFile(response);
  inputs.writeOutputsToFile(std::vector<TH1D*>{trainingRec, trainingGen, trainingFake, trainingMiss, trainingKinEff});
  inputs.writeOutputToFile(&ruBayes);
  inputs.writeOutputsToFile(std::vector<TH1D*>{testRec, testGen, testFake, testMiss, testKinEff, unfolded, refolded});
  inputs.writeOutputsToFile(std::vector<TH2D*>{covMatrix, pearson});
}

void DoUnfoldingV0Pt(InputSettings& inputs, int nIterations = -1) {
  if (nIterations <= 0)
    nIterations = inputs.nIterations;

  inputs.printLog("Doing unfolding with " + to_string(nIterations) + " iterations.", verbosityutilities::kInfo);
  inputs.printLog("Getting response from file: " + inputs.responseFileName + "\nGetting test distributions from file: " + inputs.inputFileName, verbosityutilities::kDebug);

  TFile* responseFile = TFile::Open(inputs.responseFileName.c_str(), "READ");
  if (!responseFile) {
    inputs.printLog("Error: could not open file " + inputs.responseFileName, verbosityutilities::kErrors);
    return;
  }

  TH1D* trainingRec = (TH1D*)responseFile->Get(rmutilities::unfolding::nameRecJets.c_str());
  TH1D* trainingFake = (TH1D*)responseFile->Get(rmutilities::unfolding::nameFakeJets.c_str());
  RooUnfoldResponse* response = (RooUnfoldResponse*)responseFile->Get(rmutilities::unfolding::nameResponseJets.c_str());
  if (!response)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameResponseJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingRec)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameRecJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingFake)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameFakeJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  trainingRec->SetName("trainingRecJets");
  trainingFake->SetName("trainingFakeJets");

  TH1D* testFake;
  if (inputs.responseFileName == inputs.inputFileName) // Trivial closure test: testing = training
    testFake = (TH1D*)trainingFake->Clone("testFakeJets");
  else {
    TFile* testFile = TFile::Open(inputs.inputFileName.c_str(), "READ");
    if (!testFile) {
      inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
      return;
    }
    testFake = (TH1D*)testFile->Get(rmutilities::unfolding::nameFakeJets.c_str());
  }
  if (!testFake)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameFakeJets + " in file " + inputs.inputFileName, verbosityutilities::kErrors);

  if (!response || !trainingRec || !trainingFake || !testFake)
    return;

  string ruBayesName  = rmutilities::unfolding::nameRooUnfoldBayesV0Pt + to_string(nIterations);
  string ruBayesTitle = ruBayesName;
  RooUnfoldBayes ruBayes(response, trainingRec, nIterations, inputs.doSmoothing, ruBayesName.c_str(), ruBayesTitle.c_str());

  string unfoldedName = rmutilities::unfolding::nameUnfoldedV0Pt + to_string(nIterations);
  TH1D* unfolded = (TH1D*)ruBayes.Hreco(inputs.errorTreatment);
  unfolded->SetName(unfoldedName.c_str());
  unfolded->SetTitle(unfoldedName.c_str());

  string refoldedName = rmutilities::unfolding::nameRefoldedV0Pt + to_string(nIterations);
  TH1D* refolded = (TH1D*)response->ApplyToTruth(unfolded, refoldedName.c_str());
  refolded->Add(testFake);

  // Covariance matrix
  string covMatrixName = rmutilities::unfolding::nameCovMatrixV0Pt + to_string(nIterations);
  TH2D tmp(ruBayes.Ereco(inputs.errorTreatment));
  TH2D* covMatrix = (TH2D*)tmp.Clone(covMatrixName.c_str());
  covMatrix->SetTitle(covMatrixName.c_str());

  inputs.printLog("Rec hist is of size " + to_string(trainingRec->GetNbinsX()) + ", covariance matrix is of size " + to_string(covMatrix->GetNbinsX()) + " x " + to_string(covMatrix->GetNbinsY()), verbosityutilities::kDebug);

  // Pearson coefficients
  string pearsonName = rmutilities::unfolding::namePearsonV0Pt + to_string(nIterations);
  TH2D* pearson = (TH2D*)covMatrix->Clone(pearsonName.c_str());
  pearson->SetTitle(pearsonName.c_str());
  pearson->Reset();

  for (int xCovBin = 1; xCovBin <= covMatrix->GetNbinsX(); xCovBin++) {
    for (int yCovBin = 1; yCovBin <= covMatrix->GetNbinsY(); yCovBin++) {
      double cov = covMatrix->GetBinContent(xCovBin, yCovBin);
      double sigmaX = sqrt(covMatrix->GetBinContent(xCovBin, xCovBin));
      double sigmaY = sqrt(covMatrix->GetBinContent(yCovBin, yCovBin));
      double pearsonCoeff = 0.;
      if (std::isnan(cov) || std::isnan(sigmaX) || std::isnan(sigmaY))
        continue;
      if (inputs.isVarConsistentWithZero(sigmaX) || inputs.isVarConsistentWithZero(sigmaY))
        continue;

      pearsonCoeff = cov / std::sqrt(sigmaX * sigmaY);
      pearson->SetBinContent(xCovBin, yCovBin, pearsonCoeff);
    }
  }

  inputs.writeOutputToFile(response);
  inputs.writeOutputsToFile(std::vector<TH1D*>{trainingRec, trainingFake, testFake});
  inputs.writeOutputToFile(&ruBayes);
  inputs.writeOutputsToFile(std::vector<TH1D*>{unfolded, refolded});
  inputs.writeOutputsToFile(std::vector<TH2D*>{covMatrix, pearson});
}

// ----------------------------------------------------------

array<TH1D*, 3> MakeHistsClosureTestsJet(const TH1D* hRooUnfold, const TH1D* hAnalysis, string nRooUnfold = "", string nAnalysis = "") {
  TH1D* hDiff = (TH1D*)hRooUnfold->Clone("hDiff");
  hDiff->Add(hAnalysis, -1);
  hDiff->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
  if (!nRooUnfold.empty() && !nAnalysis.empty())
    hDiff->SetName((nRooUnfold + " - " + nAnalysis).c_str());

  TH1D* hRelDiff = (TH1D*)hRooUnfold->Clone("hRelDiff");
  hRelDiff->Add(hAnalysis, -1);
  hRelDiff->Divide(hAnalysis);
  hRelDiff->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
  hRelDiff->GetXaxis()->SetTitle(("#frac{" + nRooUnfold + " - " + nAnalysis + "}{" + nAnalysis + "}").c_str());

  TH1D* hRatio = (TH1D*)hRooUnfold->Clone("hRatio");
  hRatio->Divide(hAnalysis);
  hRatio->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
  hRatio->SetTitle(("#frac{" + nRooUnfold + "}{" + nAnalysis + "}").c_str());

  return {hDiff, hRelDiff, hRatio};
}

string TrivialClosureTestSummary(InputSettings& inputs, array<TH1D*, 3> hists, bool isUnfolded) {
  bool tctFailed = false;
  string s;
  if (!inputs.isHistConsistentWithZero(hists[0])) {
    tctFailed = true;
    s.append("\n"); s.append(hists[0]->GetTitle()); s.append(" distribution is NOT consistent with zero!");
  }
  if (!inputs.isHistConsistentWithZero(hists[1])) {
    tctFailed = true;
    s.append("\n"); s.append(hists[1]->GetTitle()); s.append(" distribution is NOT consistent with zero!");
  }
  if (!inputs.isHistConsistentWithOne(hists[2])) {
    tctFailed = true;
    s.append("\n"); s.append(hists[2]->GetTitle()); s.append(" distribution is NOT consistent with one!");
  }

  string result = "Trivial closure test: ";
  if (tctFailed) {
    result += "failed!" + s;
  } else {
    if (isUnfolded)
      result += "success! Unfolded distribution is consistent with generated distribution.";
    else
      result += "success! Refolded distribution is consistent with reconstructed distribution.";
  }
  return result;
}

void MakePlotsJets(InputSettings& inputs, int nIteration = -1, bool saveFigs = true, bool drawText = false) {
  if (nIteration <= 0)
    nIteration = inputs.nIterations;

  inputs.printLog("Making plots for closure test for iteration = " + to_string(nIteration), verbosityutilities::kInfo);

  TFile* file = TFile::Open(inputs.inputFileName.c_str(), "READ");
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }

  TH1D* unfolded = (TH1D*)file->Get((rmutilities::unfolding::nameUnfoldedJets + to_string(nIteration)).c_str());
  TH1D* refolded = (TH1D*)file->Get((rmutilities::unfolding::nameRefoldedJets + to_string(nIteration)).c_str());
  TH1D* testGen  = (TH1D*)file->Get(rmutilities::testing::nameGenJets.c_str());
  TH1D* testRec  = (TH1D*)file->Get(rmutilities::testing::nameRecJets.c_str());

  array<TH1D*, 3> hUnfGen = MakeHistsClosureTestsJet(unfolded, testGen);
  TH1D* unfoldedMinusGen        = hUnfGen[0];
  TH1D* unfoldedMinusGenOverGen = hUnfGen[1];
  TH1D* unfoldedOverGen         = hUnfGen[2];

  array<TH1D*, 3> hRefRec = MakeHistsClosureTestsJet(refolded, testRec);
  TH1D* refoldedMinusRec        = hRefRec[0];
  TH1D* refoldedMinusRecOverRec = hRefRec[1];
  TH1D* refoldedOverRec         = hRefRec[2];

  string drawOption = drawText ? "hist text" : "";
  string canvasName = "cUnfolded_" + to_string(nIteration);
  TCanvas* cUnfolded = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
  cUnfolded->Divide(3, 1);
  for (int iPad = 1; iPad <= 3; iPad++) {
    cUnfolded->cd(iPad);
    hUnfGen[iPad - 1]->Draw(drawOption.c_str());
  }
  if (saveFigs)
    cUnfolded->SaveAs(("unfolded-jets-iteration" + to_string(nIteration) + ".pdf").c_str());

  canvasName = "cRefolded_" + to_string(nIteration);
  TCanvas* cRefolded = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
  cRefolded->Divide(3, 1);
  for (int iPad = 1; iPad <= 3; iPad++) {
    cRefolded->cd(iPad);
    hRefRec[iPad - 1]->Draw(drawOption.c_str());
  }
  if (saveFigs)
    cRefolded->SaveAs(("refolded-jets-iteration" + to_string(nIteration) + ".pdf").c_str());

  if (inputs.doTrivialClosureTest) {
    const bool isUnfolded = true;
    inputs.printLog(TrivialClosureTestSummary(inputs, hUnfGen, isUnfolded), verbosityutilities::kInfo);
    inputs.printLog(TrivialClosureTestSummary(inputs, hRefRec, !isUnfolded), verbosityutilities::kInfo);
  }
}

// ----------------------------------------------------------
//
// Interface
//
// ----------------------------------------------------------

// Trivial Closure Test
void doclosuretest() {
  gROOT->SetBatch(kTRUE);
  InputSettings x; x.verbosity = verbosityutilities::kDebug;
  x.train = 468659;
  x.inputFileName = to_string(x.train) + ".root";
  x.outputFileName = "RooUnfoldResponse_" + to_string(x.train) + ".root";
  x.responseFileName = x.outputFileName;
  x.doTrivialClosureTest = true;
  x.minIteration = 3;
  x.maxIteration = 5;

  // Create the response from the input file name
  x.rmHistName = rmutilities::analysis::nameJetRm;
  x.binwidthptjet = 5.;
  array<double, 2> ptjetRecRangeForUnfolding = {10., 100.};
  array<double, 2> ptjetGenRangeForUnfolding = {10., 100.};
  x.setPtJetRec(ptjetRecRangeForUnfolding[0], ptjetRecRangeForUnfolding[1]);
  x.setPtJetGen(ptjetGenRangeForUnfolding[0], ptjetGenRangeForUnfolding[1]);
  CreateResponseJets(x);

  x.rmHistName = rmutilities::analysis::nameV0PtRm;
  x.binwidthptv0 = 1.;
  x.setPtV0Rec(1., x.ptjetmaxRec); // Should the max be 40 here? Avoid empty bins!
  x.setPtV0Gen(1., x.ptjetmaxGen);
  CreateResponseV0Pt(x);

  // Do unfolding with the given response, save the results to the output file
  // For trivial closure test, input = training. Otherwise, specify a different input file
  if (!x.doTrivialClosureTest)
    x.inputFileName  = to_string(x.train) + ".root";

  x.outputFileName = "ClosureTest.root";
  for (int iIteration = x.minIteration; iIteration <= x.maxIteration; iIteration++) {
    DoUnfoldingJets(x, iIteration);
    DoUnfoldingV0Pt(x, iIteration);
  }

  // Get the unfolded and refolded distributions and compare them to the generated and reconstructed ones from the test file
  x.inputFileName = x.outputFileName;
  for (int iIteration = x.minIteration; iIteration <= x.maxIteration; iIteration++) {
    MakePlotsJets(x, iIteration, true, false);
  }
}

# endif