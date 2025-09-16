
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
} // namespace verbosityutilities

namespace plotutilities {
  enum PlotType {kJet, kV0Pt, kV0Z};
  bool is_valid(int x) { return (x >= kJet && x <= kV0Z); }
  string to_string(PlotType x) {
    switch (x) {
      case kJet:  return "jet";
      case kV0Pt: return "v0Pt";
      case kV0Z:  return "v0Z";
      default:    return "Unknown";
    }
  }
} // namespace plotutilities

namespace rmutilities {
  // Strings for loading from train output
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

    const string nameRooUnfoldBayesV0Z  = "ruBayesV0Z";
    const string nameUnfoldedV0Z        = "unfoldedV0Z";
    const string nameRefoldedV0Z        = "refoldedV0Z";
    const string nameCovMatrixV0Z       = "covMatrixV0Z";
    const string namePearsonV0Z         = "pearsonV0Z";

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

    const string nameResponseV0Z  = "responseV0Z";
    const string nameRmV0Z        = "responseMatrixV0Z";
    const string nameRecV0Z       = "RecV0Z";
    const string nameGenV0Z       = "GenV0Z";
    const string nameMissV0Z      = "MissV0Z";
    const string nameFakeV0Z      = "FakeV0Z";
    const string nameKinEffV0Z    = "KinEffV0Z";
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

    const string nameResponseV0Z  = unfolding::nameResponseV0Z;
    const string nameRmV0Z        = unfolding::nameRmV0Z;
    const string nameRecV0Z       = unfolding::prefixTraining + unfolding::nameRecV0Z;
    const string nameGenV0Z       = unfolding::prefixTraining + unfolding::nameGenV0Z;
    const string nameMissV0Z      = unfolding::prefixTraining + unfolding::nameMissV0Z;
    const string nameFakeV0Z      = unfolding::prefixTraining + unfolding::nameFakeV0Z;
    const string nameKinEffV0Z    = unfolding::prefixTraining + unfolding::nameKinEffV0Z;
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

    const string nameResponseV0Z  = unfolding::nameResponseV0Z;
    const string nameRmV0Z        = unfolding::nameRmV0Z;
    const string nameRecV0Z       = unfolding::prefixTest + unfolding::nameRecV0Z;
    const string nameGenV0Z       = unfolding::prefixTest + unfolding::nameGenV0Z;
    const string nameMissV0Z      = unfolding::prefixTest + unfolding::nameMissV0Z;
    const string nameFakeV0Z      = unfolding::prefixTest + unfolding::nameFakeV0Z;
    const string nameKinEffV0Z    = unfolding::prefixTest + unfolding::nameKinEffV0Z;
  }
} // namespace rmutilities

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

    verbosityutilities::Verbosity verbosity = verbosityutilities::kInfo;
    plotutilities::PlotType plottype;
  public:
    // Unfolding settings
    int doSmoothing = 0, nIterations = 1;
    int minIteration = 1, maxIteration = 1;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovToy;

    // General settings
    int train;
    string inputFileName, outputFileName, responseFileName;
    string hadron;

    double ptjetminGen, ptjetmaxGen, ptjetminRec, ptjetmaxRec;
    double ptv0minGen, ptv0maxGen, ptv0minRec, ptv0maxRec;
    double zv0minGen, zv0maxGen, zv0minRec, zv0maxRec;
    double etamin, etamax;
    double binwidthptjet, binwidthptv0, binwidthzv0;
    double ptjetminProjection, ptjetmaxProjection;

    bool saveFigs = true, drawText = false;
    bool doTrivialClosureTest = true;

    // Using private methods
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
    bool setPtJetProjection(double a, double b) {
      return setVariable(a, b, ptjetminProjection, ptjetmaxProjection);
    }
    bool setPtV0Gen(double a, double b) {
      return setVariable(a, b, ptv0minGen, ptv0maxGen);
    }
    bool setPtV0Rec(double a, double b) {
      return setVariable(a, b, ptv0minRec, ptv0maxRec);
    }
    bool setZV0Gen(double a, double b) {
      return setVariable(a, b, zv0minGen, zv0maxGen);
    }
    bool setZV0Rec(double a, double b) {
      return setVariable(a, b, zv0minRec, zv0maxRec);
    }

    // Getters for private variables
    plotutilities::PlotType getPlotType() { return plottype; }
    string getRmHistName() {
      switch (plottype) {
        case plotutilities::kJet:  return rmutilities::analysis::nameJetRm;
        case plotutilities::kV0Pt: return rmutilities::analysis::nameV0PtRm;
        case plotutilities::kV0Z:  return rmutilities::analysis::nameV0ZRm;
        default: return "";
      }
    }
    verbosityutilities::Verbosity getVerbosity() { return verbosity; }
    // Setters for private variables
    bool setPlotType(plotutilities::PlotType p);
    bool setVerbosity(verbosityutilities::Verbosity v);

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
    bool passVerbosityCheck(verbosityutilities::Verbosity messageVerbLevel) {
      return verbosityutilities::passVerbosityCheck(messageVerbLevel, verbosity);
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

bool InputSettings::setPlotType(plotutilities::PlotType p) {
  if (!plotutilities::is_valid(p)) {
    printLog("InputSettings::setPlotType() Error: invalid plot type", verbosityutilities::kErrors);
    return false;
  }
  plottype = p;
  return true;
}

bool InputSettings::setVerbosity(verbosityutilities::Verbosity v) {
  if (!verbosityutilities::is_valid(v)) {
    printLog("InputSettings::setVerbosity() Error: invalid verbosity level", verbosityutilities::kErrors);
    return false;
  }
  verbosity = v;
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

void FillPearsonFromCovMatrix(InputSettings& inputs, TH2D* covMatrix, TH2D* pearson) {
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

array<TH1D*, 5> CreateDistributionHistsJet(InputSettings& inputs, TH2D* responseMatrix) {
  inputs.printLog("Creating distribution histograms for jets.", verbosityutilities::kDebug);
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
  return std::array<TH1D*, 5>{hRec, hGen, hMiss, hKinEff, hFake};
}

array<TH2D*, 5> CreateDistributionHistsV0(InputSettings& inputs, THnSparseD* responseMatrix) {
  inputs.printLog("Creating histograms for V0s.", verbosityutilities::kDebug);
  int axisPtJetRec, axisPtJetGen, axisV0Rec, axisV0Gen;
  double v0minRec, v0maxRec, v0minGen, v0maxGen, binwidthv0;
  string nRec, nGen, nMiss, nKinEff, nFake, nAxes;

  if (inputs.getPlotType() == plotutilities::PlotType::kV0Pt) {
    axisPtJetRec = rmutilities::analysis::axisv0PtRmPtJetRec;
    axisPtJetGen = rmutilities::analysis::axisv0PtRmPtJetGen;
    axisV0Rec    = rmutilities::analysis::axisv0PtRmV0Rec;
    axisV0Gen    = rmutilities::analysis::axisv0PtRmV0Gen;
    v0minRec     = inputs.ptv0minRec;
    v0maxRec     = inputs.ptv0maxRec;
    v0minGen     = inputs.ptv0minGen;
    v0maxGen     = inputs.ptv0maxGen;
    binwidthv0   = inputs.binwidthptv0;
    nRec         = rmutilities::unfolding::nameRecV0Pt;
    nGen         = rmutilities::unfolding::nameGenV0Pt;
    nMiss        = rmutilities::unfolding::nameMissV0Pt;
    nKinEff      = rmutilities::unfolding::nameKinEffV0Pt;
    nFake        = rmutilities::unfolding::nameFakeV0Pt;
    nAxes        = "p_{T,V0},p_{T,jet}";
  } else if (inputs.getPlotType() == plotutilities::PlotType::kV0Z) {
    axisPtJetRec = rmutilities::analysis::axisv0ZRmPtJetRec;
    axisPtJetGen = rmutilities::analysis::axisv0ZRmPtJetGen;
    axisV0Rec    = rmutilities::analysis::axisv0ZRmV0Rec;
    axisV0Gen    = rmutilities::analysis::axisv0ZRmV0Gen;
    v0minRec     = inputs.zv0minRec;
    v0maxRec     = inputs.zv0maxRec;
    v0minGen     = inputs.zv0minGen;
    v0maxGen     = inputs.zv0maxGen;
    binwidthv0   = inputs.binwidthzv0;
    nRec         = rmutilities::unfolding::nameRecV0Z;
    nGen         = rmutilities::unfolding::nameGenV0Z;
    nMiss        = rmutilities::unfolding::nameMissV0Z;
    nKinEff      = rmutilities::unfolding::nameKinEffV0Z;
    nFake        = rmutilities::unfolding::nameFakeV0Z;
    nAxes        = "z_{T,V0},p_{T,jet}";
  } else {
    inputs.printLog("CreateDistributionHistsV0() Error: invalid plot type " + plotutilities::to_string(inputs.getPlotType()), verbosityutilities::kErrors);
    return std::array<TH2D*, 5>{nullptr, nullptr, nullptr, nullptr, nullptr};
  }

  array<int, 2> ptjetRecBins = getProjectionBins(responseMatrix->GetAxis(axisPtJetRec), inputs.ptjetminRec, inputs.ptjetmaxRec);
  array<int, 2> ptjetGenBins = getProjectionBins(responseMatrix->GetAxis(axisPtJetGen), inputs.ptjetminGen, inputs.ptjetmaxGen);
  array<int, 2> v0RecBins    = getProjectionBins(responseMatrix->GetAxis(axisV0Rec), v0minRec, v0maxRec);
  array<int, 2> v0GenBins    = getProjectionBins(responseMatrix->GetAxis(axisV0Gen), v0minGen, v0maxGen);

  array<double, 2> ptjetRecBinEdges = getProjectionEdges(responseMatrix->GetAxis(axisPtJetRec), ptjetRecBins);
  array<double, 2> ptjetGenBinEdges = getProjectionEdges(responseMatrix->GetAxis(axisPtJetGen), ptjetGenBins);
  array<double, 2> v0RecBinEdges    = getProjectionEdges(responseMatrix->GetAxis(axisV0Rec), v0RecBins);
  array<double, 2> v0GenBinEdges    = getProjectionEdges(responseMatrix->GetAxis(axisV0Gen), v0GenBins);

  int nbinsPtJetRec = (int)((ptjetRecBinEdges[1] - ptjetRecBinEdges[0]) / inputs.binwidthptjet);
  int nbinsPtJetGen = (int)((ptjetGenBinEdges[1] - ptjetGenBinEdges[0]) / inputs.binwidthptjet);
  int nbinsv0Rec    = (int)((v0RecBinEdges[1] - v0RecBinEdges[0]) / binwidthv0);
  int nbinsv0Gen    = (int)((v0GenBinEdges[1] - v0GenBinEdges[0]) / binwidthv0);

  TH2D* hRec = new TH2D(nRec.c_str(), ("Rec;" + nAxes).c_str(), nbinsv0Rec, v0RecBinEdges[0], v0RecBinEdges[1], nbinsPtJetRec, ptjetRecBinEdges[0], ptjetRecBinEdges[1]);
  TH2D* hGen = new TH2D(nGen.c_str(), ("Gen;" + nAxes).c_str(), nbinsv0Gen, v0GenBinEdges[0], v0GenBinEdges[1], nbinsPtJetGen, ptjetGenBinEdges[0], ptjetGenBinEdges[1]);

  TH2D* hMiss   = (TH2D*)hGen->Clone(nMiss.c_str());
  TH2D* hKinEff = (TH2D*)hGen->Clone(nKinEff.c_str());
  TH2D* hFake   = (TH2D*)hRec->Clone(nFake.c_str());
  return std::array<TH2D*, 5>{hRec, hGen, hMiss, hKinEff, hFake};
}

template <typename T>
T LoadResponseMatrix(InputSettings& inputs) {
  inputs.printLog("Opening input file: " + inputs.inputFileName, verbosityutilities::kDebug);
  TFile* file = TFile::Open(inputs.inputFileName.c_str());
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return nullptr;
  }
  inputs.printLog("Retrieving histogram: " + inputs.getRmHistName(), verbosityutilities::kDebug);
  T responseMatrix = (T)file->Get(inputs.getRmHistName().c_str());
  if (!responseMatrix) {
    inputs.printLog("Error: could not find " + inputs.getRmHistName() + " in file " + inputs.inputFileName, verbosityutilities::kErrors);
    return nullptr;
  }
  return responseMatrix;
}

void CreateResponseJets(InputSettings& inputs) {
  inputs.printLog("Creating response for jets.", verbosityutilities::kInfo);
  inputs.setPlotType(plotutilities::PlotType::kJet);

  TH2D* responseMatrix = LoadResponseMatrix<TH2D*>(inputs);
  responseMatrix->SetName(rmutilities::unfolding::nameRmJets.c_str());
  responseMatrix->SetTitle(rmutilities::unfolding::nameRmJets.c_str());

  array<TH1D*, 5> hists = CreateDistributionHistsJet(inputs, responseMatrix);
  TH1D* hRec    = hists[0];
  TH1D* hGen    = hists[1];
  TH1D* hMiss   = hists[2];
  TH1D* hKinEff = hists[3];
  TH1D* hFake   = hists[4];

  inputs.printLog("Filling response.", verbosityutilities::kDebug);
  RooUnfoldResponse* response = new RooUnfoldResponse(rmutilities::unfolding::nameResponseJets.c_str(), rmutilities::unfolding::nameResponseJets.c_str());
  response->Setup(hRec, hGen);
  FillFromRmJets(inputs, responseMatrix, hRec, hGen, hMiss, hKinEff, hFake, response);

  inputs.printLog("Saving objects to output file: " + inputs.outputFileName, verbosityutilities::kDebug);
  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(responseMatrix);
  inputs.writeOutputsToFile(std::vector<TH1D*>({hRec, hGen, hMiss, hFake, hKinEff}));
  return;
}

void CreateResponseV0(InputSettings& inputs) {
  inputs.printLog("Creating response for " + plotutilities::to_string(inputs.getPlotType()), verbosityutilities::kInfo);

  string nameRm = rmutilities::unfolding::nameRmV0Pt;
  string nameResponse = rmutilities::unfolding::nameResponseV0Pt;
  if (inputs.getPlotType() == plotutilities::PlotType::kV0Z) {
    nameRm = rmutilities::unfolding::nameRmV0Z;
    nameResponse = rmutilities::unfolding::nameResponseV0Z;
  }

  THnSparseD* responseMatrix = LoadResponseMatrix<THnSparseD*>(inputs);
  responseMatrix->SetName(nameRm.c_str());
  responseMatrix->SetTitle(nameRm.c_str());

  array<TH2D*, 5> hists = CreateDistributionHistsV0(inputs, responseMatrix);
  TH2D* hRec    = hists[0];
  TH2D* hGen    = hists[1];
  TH2D* hMiss   = hists[2];
  TH2D* hKinEff = hists[3];
  TH2D* hFake   = hists[4];

  inputs.printLog("Filling response.", verbosityutilities::kDebug);
  RooUnfoldResponse* response = new RooUnfoldResponse(nameResponse.c_str(), nameResponse.c_str());
  response->Setup(hRec, hGen);
  if (inputs.getPlotType() == plotutilities::PlotType::kV0Pt)
    FillFromRmV0Pt(inputs, responseMatrix, hRec, hGen, hMiss, hKinEff, hFake, response);
  else if (inputs.getPlotType() == plotutilities::PlotType::kV0Z)
    FillFromRmV0Z(inputs, responseMatrix, hRec, hGen, hMiss, hKinEff, hFake, response);

  inputs.printLog("Saving objects to output file: " + inputs.outputFileName, verbosityutilities::kDebug);
  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(responseMatrix);
  inputs.writeOutputsToFile(std::vector<TH2D*>({hRec, hGen, hMiss, hKinEff, hFake}));
}

void CreateResponseV0Pt(InputSettings& inputs) {
  inputs.setPlotType(plotutilities::PlotType::kV0Pt);
  CreateResponseV0(inputs);
}

void CreateResponseV0Z(InputSettings& inputs) {
  inputs.setPlotType(plotutilities::PlotType::kV0Z);
  CreateResponseV0(inputs);
}

// ----------------------------------------------------------

void DoUnfoldingJets(InputSettings& inputs, int nIterations) {
  inputs.printLog("Doing unfolding for jets with " + to_string(nIterations) + " iterations.", verbosityutilities::kInfo);
  inputs.setPlotType(plotutilities::PlotType::kJet);

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
    TH2D* testResponseMatrix = (TH2D*)testFile->Get(inputs.getRmHistName().c_str());
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
  FillPearsonFromCovMatrix(inputs, covMatrix, pearson);

  inputs.printLog("Unfolding done. Writing outputs to file " + inputs.outputFileName, verbosityutilities::kInfo);
  inputs.writeOutputToFile(response);
  inputs.writeOutputsToFile(std::vector<TH1D*>{trainingRec, trainingGen, trainingFake, trainingMiss, trainingKinEff});
  inputs.writeOutputToFile(&ruBayes);
  inputs.writeOutputsToFile(std::vector<TH1D*>{testRec, testGen, testFake, testMiss, testKinEff, unfolded, refolded});
  inputs.writeOutputsToFile(std::vector<TH2D*>{covMatrix, pearson});
}

void DoUnfoldingV0(InputSettings& inputs, int nIterations) {
  inputs.printLog("Doing unfolding for " + plotutilities::to_string(inputs.getPlotType()) + " with " + to_string(nIterations) + " iterations.", verbosityutilities::kInfo);

  string nameRec, nameGen, nameFake, nameMiss, nameKinEff, nameResponse, ruBayesName, unfoldedName, refoldedName, covMatrixName, pearsonName;
  if (inputs.getPlotType() == plotutilities::PlotType::kV0Pt) {
    nameRec       = rmutilities::unfolding::nameRecV0Pt;
    nameGen       = rmutilities::unfolding::nameGenV0Pt;
    nameFake      = rmutilities::unfolding::nameFakeV0Pt;
    nameMiss      = rmutilities::unfolding::nameMissV0Pt;
    nameKinEff    = rmutilities::unfolding::nameKinEffV0Pt;
    nameResponse  = rmutilities::unfolding::nameResponseV0Pt;
    ruBayesName   = rmutilities::unfolding::nameRooUnfoldBayesV0Pt + to_string(nIterations);
    unfoldedName  = rmutilities::unfolding::nameUnfoldedV0Pt + to_string(nIterations);
    refoldedName  = rmutilities::unfolding::nameRefoldedV0Pt + to_string(nIterations);
    covMatrixName = rmutilities::unfolding::nameCovMatrixV0Pt + to_string(nIterations);
    pearsonName   = rmutilities::unfolding::namePearsonV0Pt + to_string(nIterations);
  } else if (inputs.getPlotType() == plotutilities::PlotType::kV0Z) {
    nameRec       = rmutilities::unfolding::nameRecV0Z;
    nameGen       = rmutilities::unfolding::nameGenV0Z;
    nameFake      = rmutilities::unfolding::nameFakeV0Z;
    nameMiss      = rmutilities::unfolding::nameMissV0Z;
    nameKinEff    = rmutilities::unfolding::nameKinEffV0Z;
    nameResponse  = rmutilities::unfolding::nameResponseV0Z;
    ruBayesName   = rmutilities::unfolding::nameRooUnfoldBayesV0Z + to_string(nIterations);
    unfoldedName  = rmutilities::unfolding::nameUnfoldedV0Z + to_string(nIterations);
    refoldedName  = rmutilities::unfolding::nameRefoldedV0Z + to_string(nIterations);
    covMatrixName = rmutilities::unfolding::nameCovMatrixV0Z + to_string(nIterations);
    pearsonName   = rmutilities::unfolding::namePearsonV0Z + to_string(nIterations);
  } else {
    inputs.printLog("DoUnfoldingV0() Error: invalid plot type " + plotutilities::to_string(inputs.getPlotType()), verbosityutilities::kErrors);
    return;
  }

  inputs.printLog("Getting response from file: " + inputs.responseFileName + "\nGetting test distributions from file: " + inputs.inputFileName, verbosityutilities::kDebug);
  TFile* responseFile = TFile::Open(inputs.responseFileName.c_str(), "READ");
  if (!responseFile) {
    inputs.printLog("Error: could not open file " + inputs.responseFileName, verbosityutilities::kErrors);
    return;
  }

  inputs.printLog("Retrieving response and histograms from training.", verbosityutilities::kDebug);
  TH2D* trainingRec = (TH2D*)responseFile->Get(nameRec.c_str());
  TH2D* trainingGen = (TH2D*)responseFile->Get(nameGen.c_str());
  TH2D* trainingFake = (TH2D*)responseFile->Get(nameFake.c_str());
  TH2D* trainingMiss = (TH2D*)responseFile->Get(nameMiss.c_str());
  TH2D* trainingKinEff = (TH2D*)responseFile->Get(nameKinEff.c_str());
  RooUnfoldResponse* response = (RooUnfoldResponse*)responseFile->Get(nameResponse.c_str());

  if (!trainingRec)
  inputs.printLog("Error: could not find " + nameRec + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingGen)
    inputs.printLog("Error: could not find " + nameGen + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingFake)
    inputs.printLog("Error: could not find " + nameFake + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingMiss)
    inputs.printLog("Error: could not find " + nameMiss + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingKinEff)
    inputs.printLog("Error: could not find " + nameKinEff + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!response)
    inputs.printLog("Error: could not find " + nameResponse + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingRec || !trainingGen || !trainingFake || !trainingMiss || !trainingKinEff || !response)
    return;

  // Change the names to distinguish these from the test histograms when writing to file
  string prefTraining = rmutilities::unfolding::prefixTraining;
  trainingRec->SetName((prefTraining + nameRec).c_str());
  trainingGen->SetName((prefTraining + nameGen).c_str());
  trainingFake->SetName((prefTraining + nameFake).c_str());
  trainingMiss->SetName((prefTraining + nameMiss).c_str());
  trainingKinEff->SetName((prefTraining + nameKinEff).c_str());

  // Create the test histograms. These must have the same binning as the training histograms
  // For trivial closure test, they are copies, otherwise they are filled with independent data
  string prefTest = rmutilities::unfolding::prefixTest;
  TH2D* testRec = (TH2D*)trainingRec->Clone((prefTest + nameRec).c_str());
  TH2D* testGen = (TH2D*)trainingGen->Clone((prefTest + nameGen).c_str());
  TH2D* testFake = (TH2D*)trainingFake->Clone((prefTest + nameFake).c_str());
  TH2D* testMiss = (TH2D*)trainingMiss->Clone((prefTest + nameMiss).c_str());
  TH2D* testKinEff = (TH2D*)trainingKinEff->Clone((prefTest + nameKinEff).c_str());

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
    THnSparseD* testResponseMatrix = LoadResponseMatrix<THnSparseD*>(inputs);
    if (inputs.getPlotType() == plotutilities::PlotType::kV0Pt)
      FillFromRmV0Pt(inputs, testResponseMatrix, testRec, testGen, testMiss, testKinEff, testFake);
    else if (inputs.getPlotType() == plotutilities::PlotType::kV0Z)
      FillFromRmV0Z(inputs, testResponseMatrix, testRec, testGen, testMiss, testKinEff, testFake);
    testFile->Close();
  }

  string ruBayesTitle = ruBayesName;
  RooUnfoldBayes ruBayes(response, trainingRec, nIterations, inputs.doSmoothing, ruBayesName.c_str(), ruBayesTitle.c_str());

  TH2D* unfolded = (TH2D*)ruBayes.Hreco(inputs.errorTreatment);
  unfolded->SetName(unfoldedName.c_str());
  unfolded->SetTitle(unfoldedName.c_str());

  TH2D* refolded = (TH2D*)response->ApplyToTruth(unfolded, refoldedName.c_str());
  refolded->Add(testFake);

  TH2D tmp(ruBayes.Ereco(inputs.errorTreatment));
  TH2D* covMatrix = (TH2D*)tmp.Clone(covMatrixName.c_str());
  covMatrix->SetTitle(covMatrixName.c_str());

  inputs.printLog("Rec hist is of size " + to_string(trainingRec->GetNbinsX()) + " x " + to_string(trainingRec->GetNbinsY()) + ", covariance matrix is of size " + to_string(covMatrix->GetNbinsX()) + " x " + to_string(covMatrix->GetNbinsY()), verbosityutilities::kDebug);

  TH2D* pearson = (TH2D*)covMatrix->Clone(pearsonName.c_str());
  pearson->SetTitle(pearsonName.c_str());
  FillPearsonFromCovMatrix(inputs, covMatrix, pearson);

  inputs.writeOutputToFile(response);
  inputs.writeOutputsToFile(std::vector<TH2D*>{trainingRec, trainingGen, trainingFake, trainingMiss, trainingKinEff});
  inputs.writeOutputToFile(&ruBayes);
  inputs.writeOutputsToFile(std::vector<TH2D*>{unfolded, refolded});
  inputs.writeOutputsToFile(std::vector<TH2D*>{testRec, testGen, testFake, testMiss, testKinEff});
  inputs.writeOutputsToFile(std::vector<TH2D*>{covMatrix, pearson});
}

void DoUnfoldingV0Pt(InputSettings& inputs, int nIterations) {
  inputs.setPlotType(plotutilities::PlotType::kV0Pt);
  DoUnfoldingV0(inputs, nIterations);
}

void DoUnfoldingV0Z(InputSettings& inputs, int nIterations) {
  inputs.setPlotType(plotutilities::PlotType::kV0Z);
  DoUnfoldingV0(inputs, nIterations);
}

// ----------------------------------------------------------

string TrivialClosureTestSummary(InputSettings& inputs, array<TH1D*, 3> hists, bool isUnfolded) {
  bool tctFailed = false;
  string s;
  if (!inputs.isHistConsistentWithZero(hists[0])) {
    tctFailed = true;
    s.append("\n"); s.append(hists[0]->GetTitle()); s.append(" distribution is NOT consistent with zero! (within 1e-10)");
  }
  if (!inputs.isHistConsistentWithZero(hists[1])) {
    tctFailed = true;
    s.append("\n"); s.append(hists[1]->GetTitle()); s.append(" distribution is NOT consistent with zero! (within 1e-10)");
  }
  if (!inputs.isHistConsistentWithOne(hists[2])) {
    tctFailed = true;
    s.append("\n"); s.append(hists[2]->GetTitle()); s.append(" distribution is NOT consistent with one! (within 1e-10)");
  }

  string result = "Trivial closure test: ";
  if (tctFailed) {
    result += "failed!" + s;
  } else {
    if (isUnfolded)
      result += "success! Unfolded distribution is consistent with generated distribution. (within 1e-10)";
    else
      result += "success! Refolded distribution is consistent with reconstructed distribution. (within 1e-10)";
  }
  return result;
}

void MakeClosureTestPlots(InputSettings& inputs, array<TH1D*, 3> hists, int nIteration, bool isUnfolded) {
  inputs.printLog(TString::Format("Making plots for %s closure test for iteration = %d, ptjet = %.f-%.f", plotutilities::to_string(inputs.getPlotType()).c_str(), nIteration, inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);

  string drawOption = inputs.drawText ? "hist text" : "";
  string canvasName = isUnfolded ? "unfolded" : "refolded";
  canvasName += plotutilities::to_string(inputs.getPlotType()) + "-iteration" + to_string(nIteration) + ".pdf";
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
  canvas->Divide(hists.size(), 1);
  for (int iPad = 1; iPad <= hists.size(); iPad++) {
    canvas->cd(iPad);
    hists[iPad - 1]->Draw(drawOption.c_str());
  }
  if (inputs.saveFigs)
    canvas->SaveAs(canvas->GetName());
}

template <typename T>
array<TH1D*, 3> MakeClosureTestHists(const T* hRooUnfold, const T* hAnalysis, double ptjetmin, double ptjetmax, bool isUnfolded, plotutilities::PlotType plot) {
  string nRooUnfold = isUnfolded ? "Unfolded" : "Refolded";
  string nAnalysis = isUnfolded ? "Generated" : "Reconstructed";
  string nObservable;
  TH1D* hRoo;
  TH1D* hAna;
  array<int, 2> ptjetBins;

  switch (plot) {
    case plotutilities::PlotType::kJet:
      nObservable = "#it{p}_{T,jet} (GeV/c)";
      hRoo = (TH1D*)hRooUnfold->Clone("hRoo");
      hAna = (TH1D*)hAnalysis->Clone("hAna");
      ptjetBins = getProjectionBins(hRooUnfold->GetXaxis(), ptjetmin, ptjetmax);
      hRoo->GetXaxis()->SetRange(ptjetBins[0], ptjetBins[1]);
      hAna->GetXaxis()->SetRange(ptjetBins[0], ptjetBins[1]);
      break;
    case plotutilities::PlotType::kV0Pt:
      nObservable = "#it{p}_{T,V0} (GeV/c)";
      ptjetBins = getProjectionBins(hRooUnfold->GetYaxis(), ptjetmin, ptjetmax);
      hRoo = (TH1D*)hRooUnfold->ProjectionX("hRoo", ptjetBins[0], ptjetBins[1]);
      hAna = (TH1D*)hAnalysis->ProjectionX("hAna", ptjetBins[0], ptjetBins[1]);
      break;
    case plotutilities::PlotType::kV0Z:
      nObservable = "#it{z}_{V0}";
      ptjetBins = getProjectionBins(hRooUnfold->GetYaxis(), ptjetmin, ptjetmax);
      hRoo = (TH1D*)hRooUnfold->ProjectionX("hDiff", ptjetBins[0], ptjetBins[1]);
      hAna = (TH1D*)hAnalysis->ProjectionX("hAnalysisProj", ptjetBins[0], ptjetBins[1]);
      break;
    default:
      verbosityutilities::printLog("MakeClosureTestHists: Unknown plot type!", verbosityutilities::kErrors, verbosityutilities::kErrors);
      break;
  }

  TH1D* hDiff = (TH1D*)hRoo->Clone("hDiff");
  hDiff->Add(hAna, -1);
  hDiff->GetXaxis()->SetTitle(nObservable.c_str());
  hDiff->SetName((nRooUnfold + " - " + nAnalysis).c_str());

  TH1D* hRelDiff = (TH1D*)hRoo->Clone("hRelDiff");
  hRelDiff->Add(hAna, -1);
  hRelDiff->Divide(hAna);
  hRelDiff->GetXaxis()->SetTitle(nObservable.c_str());
  hRelDiff->SetTitle(("#frac{" + nRooUnfold + " - " + nAnalysis + "}{" + nAnalysis + "}").c_str());

  TH1D* hRatio = (TH1D*)hRoo->Clone("hRatio");
  hRatio->Divide(hAna);
  hRatio->GetXaxis()->SetTitle(nObservable.c_str());
  hRatio->SetTitle(("#frac{" + nRooUnfold + "}{" + nAnalysis + "}").c_str());

  return {hDiff, hRelDiff, hRatio};
}

array<TH1D*, 3> MakeClosureTestHistsJet(const TH1D* hRooUnfold, const TH1D* hAnalysis, bool isUnfolded) {
  string nRooUnfold = isUnfolded ? "Unfolded" : "Refolded";
  string nAnalysis = isUnfolded ? "Generated" : "Reconstructed";

  TH1D* hDiff = (TH1D*)hRooUnfold->Clone("hDiff");
  hDiff->Add(hAnalysis, -1);
  hDiff->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
  hDiff->SetTitle((nRooUnfold + " - " + nAnalysis).c_str());

  TH1D* hRelDiff = (TH1D*)hRooUnfold->Clone("hRelDiff");
  hRelDiff->Add(hAnalysis, -1);
  hRelDiff->Divide(hAnalysis);
  hRelDiff->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
  hRelDiff->SetTitle(("#frac{" + nRooUnfold + " - " + nAnalysis + "}{" + nAnalysis + "}").c_str());

  TH1D* hRatio = (TH1D*)hRooUnfold->Clone("hRatio");
  hRatio->Divide(hAnalysis);
  hRatio->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
  hRatio->SetTitle(("#frac{" + nRooUnfold + "}{" + nAnalysis + "}").c_str());

  return {hDiff, hRelDiff, hRatio};
}

array<TH1D*, 3> MakeClosureTestHistsV0(InputSettings& inputs, const TH2D* hRooUnfold, const TH2D* hAnalysis, bool isUnfolded) {
  string nRooUnfold = isUnfolded ? "Unfolded" : "Refolded";
  string nAnalysis = isUnfolded ? "Generated" : "Reconstructed";

  string nObservable = "#it{p}_{T,V0} (GeV/c)";
  if (inputs.getPlotType() == plotutilities::PlotType::kV0Z)
    nObservable = "#it{z}_{V0}";

  // The binning of the two histograms is the same by construction
  array<int, 2> ptjetBins = getProjectionBins(hRooUnfold->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1D* hRooUnfoldProj = (TH1D*)hRooUnfold->ProjectionX("hRooUnfoldProj", ptjetBins[0], ptjetBins[1]);
  TH1D* hAnalysisProj  = (TH1D*)hAnalysis->ProjectionX("hAnalysisProj", ptjetBins[0], ptjetBins[1]);

  TH1D* hDiff = (TH1D*)hRooUnfoldProj->Clone("hDiff");
  hDiff->Add(hAnalysisProj, -1);
  hDiff->GetXaxis()->SetTitle(nObservable.c_str());
  hDiff->SetName((nRooUnfold + " - " + nAnalysis).c_str());

  TH1D* hRelDiff = (TH1D*)hRooUnfoldProj->Clone("hRelDiff");
  hRelDiff->Add(hAnalysisProj, -1);
  hRelDiff->Divide(hAnalysisProj);
  hRelDiff->GetXaxis()->SetTitle(nObservable.c_str());
  hRelDiff->SetTitle(("#frac{" + nRooUnfold + " - " + nAnalysis + "}{" + nAnalysis + "}").c_str());

  TH1D* hRatio = (TH1D*)hRooUnfoldProj->Clone("hRatio");
  hRatio->Divide(hAnalysisProj);
  hRatio->GetXaxis()->SetTitle(nObservable.c_str());
  hRatio->SetTitle(("#frac{" + nRooUnfold + "}{" + nAnalysis + "}").c_str());

  return {hDiff, hRelDiff, hRatio};
}

void DoClosureTestJets(InputSettings& inputs, int nIteration) {
  inputs.setPlotType(plotutilities::PlotType::kJet);

  TFile* file = TFile::Open(inputs.inputFileName.c_str(), "READ");
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }

  TH1D* unfolded = (TH1D*)file->Get((rmutilities::unfolding::nameUnfoldedJets + to_string(nIteration)).c_str());
  TH1D* refolded = (TH1D*)file->Get((rmutilities::unfolding::nameRefoldedJets + to_string(nIteration)).c_str());
  TH1D* testGen  = (TH1D*)file->Get(rmutilities::testing::nameGenJets.c_str());
  TH1D* testRec  = (TH1D*)file->Get(rmutilities::testing::nameRecJets.c_str());

  const bool isUnfolded = true;
  array<TH1D*, 3> hUnfGen = MakeClosureTestHistsJet(unfolded, testGen, isUnfolded);
  array<TH1D*, 3> hRefRec = MakeClosureTestHistsJet(refolded, testRec, !isUnfolded);
  MakeClosureTestPlots(inputs, hUnfGen, nIteration, isUnfolded);
  MakeClosureTestPlots(inputs, hRefRec, nIteration, !isUnfolded);
  if (inputs.doTrivialClosureTest) {
    inputs.printLog(TrivialClosureTestSummary(inputs, hUnfGen, isUnfolded), verbosityutilities::kInfo);
    inputs.printLog(TrivialClosureTestSummary(inputs, hRefRec, !isUnfolded), verbosityutilities::kInfo);
  }
}

void DoClosureTestV0(InputSettings& inputs, int nIteration) {
  TFile* file = TFile::Open(inputs.inputFileName.c_str(), "READ");
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }
  string nameUnfolded, nameRefolded, nameGen, nameRec;
  switch (inputs.getPlotType()) {
    case plotutilities::PlotType::kV0Pt:
      nameUnfolded = rmutilities::unfolding::nameUnfoldedV0Pt + to_string(nIteration);
      nameRefolded = rmutilities::unfolding::nameRefoldedV0Pt + to_string(nIteration);
      nameGen      = rmutilities::testing::nameGenV0Pt;
      nameRec      = rmutilities::testing::nameRecV0Pt;
    case plotutilities::PlotType::kV0Z:
      nameUnfolded = rmutilities::unfolding::nameUnfoldedV0Z + to_string(nIteration);
      nameRefolded = rmutilities::unfolding::nameRefoldedV0Z + to_string(nIteration);
      nameGen      = rmutilities::testing::nameGenV0Z;
      nameRec      = rmutilities::testing::nameRecV0Z;
      break;
    default:
      inputs.printLog("DoClosureTestV0() Error: invalid plot type " + plotutilities::to_string(inputs.getPlotType()), verbosityutilities::kErrors);
      return;
  }
  TH2D* unfolded = (TH2D*)file->Get(nameUnfolded.c_str());
  TH2D* refolded = (TH2D*)file->Get(nameRefolded.c_str());
  TH2D* testGen  = (TH2D*)file->Get(nameGen.c_str());
  TH2D* testRec  = (TH2D*)file->Get(nameRec.c_str());
  if (!unfolded)
    inputs.printLog("Error: could not find " + nameUnfolded + " in file " + inputs.inputFileName, verbosityutilities::kErrors);
  if (!refolded)
    inputs.printLog("Error: could not find " + nameRefolded + " in file " + inputs.inputFileName, verbosityutilities::kErrors);
  if (!testGen)
    inputs.printLog("Error: could not find " + nameGen + " in file " + inputs.inputFileName, verbosityutilities::kErrors);
  if (!testRec)
    inputs.printLog("Error: could not find " + nameRec + " in file " + inputs.inputFileName, verbosityutilities::kErrors);
  if (!unfolded || !refolded || !testGen || !testRec)
    return;

  const bool isUnfolded = true;
  array<TH1D*, 3> hUnfGen = MakeClosureTestHistsV0(inputs, unfolded, testGen, isUnfolded);
  array<TH1D*, 3> hRefRec = MakeClosureTestHistsV0(inputs, refolded, testRec, !isUnfolded);
  MakeClosureTestPlots(inputs, hUnfGen, nIteration, isUnfolded);
  MakeClosureTestPlots(inputs, hRefRec, nIteration, !isUnfolded);
  if (inputs.doTrivialClosureTest) {
    inputs.printLog(TrivialClosureTestSummary(inputs, hUnfGen, isUnfolded), verbosityutilities::kInfo);
    inputs.printLog(TrivialClosureTestSummary(inputs, hRefRec, !isUnfolded), verbosityutilities::kInfo);
  }
}

void DoClosureTestV0Pt(InputSettings& inputs, int nIteration) {
  inputs.setPlotType(plotutilities::PlotType::kV0Pt);
  DoClosureTestV0(inputs, nIteration);
}

void DoClosureTestV0Z(InputSettings& inputs, int nIteration) {
  inputs.setPlotType(plotutilities::PlotType::kV0Z);
  DoClosureTestV0(inputs, nIteration);
}

// ----------------------------------------------------------
//
// Interface
//
// ----------------------------------------------------------

InputSettings setup() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.train = 468659;
  x.inputFileName = to_string(x.train) + ".root";
  x.outputFileName = "RooUnfoldResponse_" + to_string(x.train) + ".root";
  x.responseFileName = x.outputFileName;
  return x;
}

void createresponses() {
  InputSettings x = setup();
  x.binwidthptjet = 5.;
  x.binwidthptv0 = 1.;
  x.binwidthzv0 = 0.1;
  array<double, 2> ptjetRecRangeForUnfolding = {10., 100.};
  array<double, 2> ptjetGenRangeForUnfolding = {10., 100.};
  array<double, 2> ptv0RecRangeForUnfolding = {1., 40.};
  array<double, 2> ptv0GenRangeForUnfolding = {1., 40.};
  array<double, 2> zv0RecRangeForUnfolding = {0., 1.};
  array<double, 2> zv0GenRangeForUnfolding = {0., 1.};

  x.setPtJetRec(ptjetRecRangeForUnfolding[0], ptjetRecRangeForUnfolding[1]);
  x.setPtJetGen(ptjetGenRangeForUnfolding[0], ptjetGenRangeForUnfolding[1]);
  x.setPtV0Rec(ptv0RecRangeForUnfolding[0], ptv0RecRangeForUnfolding[1]);
  x.setPtV0Gen(ptv0GenRangeForUnfolding[0], ptv0GenRangeForUnfolding[1]);
  x.setZV0Rec(zv0RecRangeForUnfolding[0], zv0RecRangeForUnfolding[1]);
  x.setZV0Gen(zv0GenRangeForUnfolding[0], zv0GenRangeForUnfolding[1]);

  x.setPlotType(plotutilities::PlotType::kJet);
  CreateResponseJets(x);
  CreateResponseV0Pt(x);
  CreateResponseV0Z(x);
}

void unfolding(string inputFileName = "") {
  InputSettings x = setup();
  x.minIteration = 3;
  x.maxIteration = 5;
  x.inputFileName = x.responseFileName;
  x.outputFileName = "ClosureTest.root";

  if (inputFileName.empty())
    x.doTrivialClosureTest = true;
  if (!x.doTrivialClosureTest)
    x.inputFileName = inputFileName;

  for (int iIteration = x.minIteration; iIteration <= x.maxIteration; iIteration++) {
    DoUnfoldingJets(x, iIteration);
    DoUnfoldingV0Pt(x, iIteration);
    DoUnfoldingV0Z(x, iIteration);
  }
}

void checkclosure() {
  gROOT->SetBatch(kTRUE);
  InputSettings x = setup();
  x.doTrivialClosureTest = true;
  x.minIteration = 3;
  x.maxIteration = 5;
  x.inputFileName = "ClosureTest.root";

  for (int iIteration = x.minIteration; iIteration <= x.maxIteration; iIteration++) {
    DoClosureTestJets(x, iIteration);
    x.setPtJetProjection(10., 20.);
    DoClosureTestV0Pt(x, iIteration);
    DoClosureTestV0Z(x, iIteration);
    x.setPtJetProjection(20., 30.);
    DoClosureTestV0Pt(x, iIteration);
    DoClosureTestV0Z(x, iIteration);
  }
}
// Trivial Closure Test
void doclosuretest() {
  createresponses();
  unfolding();
  checkclosure();
}

# endif