
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "../plotting_macros/histUtils.C"
#include "../plotting_macros/plotUtils.C"

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

namespace unfoldingutilities {
  enum VariableType {kJet, kK0SPt, kK0SZ};
  bool is_valid(int x) { return (x >= kJet && x <= kK0SZ); }
  string to_string(VariableType x) {
    switch (x) {
      case kJet:   return "jet";
      case kK0SPt: return "K0Spt";
      case kK0SZ:  return "K0Sz";
      default:     return "Unknown";
    }
  }
} // namespace unfoldingutilities

namespace rmutilities {
  // Strings for loading from train output
  namespace analysis {
    const int axisJetRmPtJetRec = 0;
    const int axisJetRmPtJetGen = 1;
    const int ndimJetRm         = 2;

    const int axisK0SRmPtJetRec = 0;
    const int axisK0SRmK0SRec   = 1;
    const int axisK0SRmPtJetGen = 2;
    const int axisK0SRmK0SGen   = 3;
    const int nDimK0SRm         = 4;

    const string nameDirJetRm  = "jet-fragmentation/matching/jets/";
    const string nameJetRm     = nameDirJetRm + "matchDetJetPtPartJetPt";
    const string nameJetFake   = nameDirJetRm + "fakeJetPtEtaPhi";
    const string nameJetMiss   = nameDirJetRm + "missJetPtEtaPhi";

    const string nameDirV0Rm   = nameDirJetRm + "V0/";
    const string nameK0SPtRm   = nameDirV0Rm + "partJetPtK0SPtDetJetPtK0SPtRightCollision";
    const string nameK0SPtFake = nameDirV0Rm + "fakeJetPtK0SPtEtaPhi";
    const string nameK0SPtMiss = nameDirV0Rm + "missJetPtK0SPtEtaPhi";
    const string nameK0SZRm    = nameDirV0Rm + "partJetPtK0STrackProjDetJetPtK0STrackProjRightCollision";
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

    const string nameRooUnfoldBayesK0SPt = "ruBayesK0SPt";
    const string nameUnfoldedK0SPt       = "unfoldedK0SPt";
    const string nameRefoldedK0SPt       = "refoldedK0SPt";
    const string nameCovMatrixK0SPt      = "covMatrixK0SPt";
    const string namePearsonK0SPt        = "pearsonK0SPt";

    const string nameRooUnfoldBayesK0SZ  = "ruBayesK0SZ";
    const string nameUnfoldedK0SZ        = "unfoldedK0SZ";
    const string nameRefoldedK0SZ        = "refoldedK0SZ";
    const string nameCovMatrixK0SZ       = "covMatrixK0SZ";
    const string namePearsonK0SZ         = "pearsonK0SZ";

    // For setting up training and testing namespaces
    const string prefixTraining    = "training";
    const string prefixTest        = "test";
    const string nameResponseJets  = "responseJets";
    const string nameRmJets        = "responseMatrixJets";
    const string nameRecJets       = "RecJets";
    const string nameGenJets       = "GenJets";
    const string nameMissJets      = "MissJets";
    const string nameFakeJets      = "FakeJets";
    const string nameKinEffJets    = "KinEffJets";
    const string nameGenAndRecJets = "GenAndRecJets";

    const string nameResponseK0SPt  = "responseK0SPt";
    const string nameRmK0SPt        = "responseMatrixK0SPt";
    const string nameRecK0SPt       = "RecK0SPt";
    const string nameGenK0SPt       = "GenK0SPt";
    const string nameMissK0SPt      = "MissK0SPt";
    const string nameFakeK0SPt      = "FakeK0SPt";
    const string nameKinEffK0SPt    = "KinEffK0SPt";
    const string nameGenAndRecK0SPt = "GenAndRecK0SPt";

    const string nameResponseK0SZ  = "responseK0SZ";
    const string nameRmK0SZ        = "responseMatrixK0SZ";
    const string nameRecK0SZ       = "RecK0SZ";
    const string nameGenK0SZ       = "GenK0SZ";
    const string nameMissK0SZ      = "MissK0SZ";
    const string nameFakeK0SZ      = "FakeK0SZ";
    const string nameKinEffK0SZ    = "KinEffK0SZ";
    const string nameGenAndRecK0SZ = "GenAndRecK0SZ";
  }
  namespace training {
    const string nameResponseJets  = unfolding::prefixTraining + unfolding::nameResponseJets;
    const string nameRmJets        = unfolding::prefixTraining + unfolding::nameRmJets;
    const string nameRecJets       = unfolding::prefixTraining + unfolding::nameRecJets;
    const string nameGenJets       = unfolding::prefixTraining + unfolding::nameGenJets;
    const string nameMissJets      = unfolding::prefixTraining + unfolding::nameMissJets;
    const string nameFakeJets      = unfolding::prefixTraining + unfolding::nameFakeJets;
    const string nameKinEffJets    = unfolding::prefixTraining + unfolding::nameKinEffJets;
    const string nameGenAndRecJets = unfolding::prefixTraining + unfolding::nameGenAndRecJets;

    const string nameResponseK0SPt  = unfolding::prefixTraining + unfolding::nameResponseK0SPt;
    const string nameRmK0SPt        = unfolding::prefixTraining + unfolding::nameRmK0SPt;
    const string nameRecK0SPt       = unfolding::prefixTraining + unfolding::nameRecK0SPt;
    const string nameGenK0SPt       = unfolding::prefixTraining + unfolding::nameGenK0SPt;
    const string nameMissK0SPt      = unfolding::prefixTraining + unfolding::nameMissK0SPt;
    const string nameFakeK0SPt      = unfolding::prefixTraining + unfolding::nameFakeK0SPt;
    const string nameKinEffK0SPt    = unfolding::prefixTraining + unfolding::nameKinEffK0SPt;
    const string nameGenAndRecK0SPt = unfolding::prefixTraining + unfolding::nameGenAndRecK0SPt;

    const string nameResponseK0SZ  = unfolding::prefixTraining + unfolding::nameResponseK0SZ;
    const string nameRmK0SZ        = unfolding::prefixTraining + unfolding::nameRmK0SZ;
    const string nameRecK0SZ       = unfolding::prefixTraining + unfolding::nameRecK0SZ;
    const string nameGenK0SZ       = unfolding::prefixTraining + unfolding::nameGenK0SZ;
    const string nameMissK0SZ      = unfolding::prefixTraining + unfolding::nameMissK0SZ;
    const string nameFakeK0SZ      = unfolding::prefixTraining + unfolding::nameFakeK0SZ;
    const string nameKinEffK0SZ    = unfolding::prefixTraining + unfolding::nameKinEffK0SZ;
    const string nameGenAndRecK0SZ = unfolding::prefixTraining + unfolding::nameGenAndRecK0SZ;
  }
  namespace testing {
    const string nameResponseJets  = unfolding::prefixTest + unfolding::nameResponseJets;
    const string nameRmJets        = unfolding::prefixTest + unfolding::nameRmJets;
    const string nameRecJets       = unfolding::prefixTest + unfolding::nameRecJets;
    const string nameGenJets       = unfolding::prefixTest + unfolding::nameGenJets;
    const string nameMissJets      = unfolding::prefixTest + unfolding::nameMissJets;
    const string nameFakeJets      = unfolding::prefixTest + unfolding::nameFakeJets;
    const string nameKinEffJets    = unfolding::prefixTest + unfolding::nameKinEffJets;
    const string nameGenAndRecJets = unfolding::prefixTest + unfolding::nameGenAndRecJets;

    const string nameResponseK0SPt  = unfolding::prefixTest + unfolding::nameResponseK0SPt;
    const string nameRmK0SPt        = unfolding::prefixTest + unfolding::nameRmK0SPt;
    const string nameRecK0SPt       = unfolding::prefixTest + unfolding::nameRecK0SPt;
    const string nameGenK0SPt       = unfolding::prefixTest + unfolding::nameGenK0SPt;
    const string nameMissK0SPt      = unfolding::prefixTest + unfolding::nameMissK0SPt;
    const string nameFakeK0SPt      = unfolding::prefixTest + unfolding::nameFakeK0SPt;
    const string nameKinEffK0SPt    = unfolding::prefixTest + unfolding::nameKinEffK0SPt;
    const string nameGenAndRecK0SPt = unfolding::prefixTest + unfolding::nameGenAndRecK0SPt;

    const string nameResponseK0SZ  = unfolding::prefixTest + unfolding::nameResponseK0SZ;
    const string nameRmK0SZ        = unfolding::prefixTest + unfolding::nameRmK0SZ;
    const string nameRecK0SZ       = unfolding::prefixTest + unfolding::nameRecK0SZ;
    const string nameGenK0SZ       = unfolding::prefixTest + unfolding::nameGenK0SZ;
    const string nameMissK0SZ      = unfolding::prefixTest + unfolding::nameMissK0SZ;
    const string nameFakeK0SZ      = unfolding::prefixTest + unfolding::nameFakeK0SZ;
    const string nameKinEffK0SZ    = unfolding::prefixTest + unfolding::nameKinEffK0SZ;
    const string nameGenAndRecK0SZ = unfolding::prefixTest + unfolding::nameGenAndRecK0SZ;
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
    template <typename T> void setTemplateHist(T* tmplt, T* hist);
    void setTemplateHistJet(const int nbins, const double* binedges, const bool generatorLevel);
    void setTemplateHistJet(const double min, const double max, const double binwidth, const bool generatorLevel);
    void setTemplateHistV0Pt(const int nbinsx, const double* binedgesx, const int nbinsy, const double* binedgesy, const bool generatorLevel);
    void setTemplateHistV0Pt(const double minx, const double maxx, const double binwidthx, const double miny, const double maxy, const double binwidthy, const bool generatorLevel);
    void setTemplateHistV0Z(const int nbinsx, const double* binedgesx, const int nbinsy, const double* binedgesy, const bool generatorLevel);
    void setTemplateHistV0Z(const double minx, const double maxx, const double binwidthx, const double miny, const double maxy, const double binwidthy, const bool generatorLevel);
    template <typename T> bool setVariable(T a, T b, T &x, T &y);

    verbosityutilities::Verbosity verbosity = verbosityutilities::kInfo;
    unfoldingutilities::VariableType vartype;

    TH1D* templatePtJetGen = nullptr;
    TH1D* templatePtJetRec = nullptr;
    TH2D* templatePtV0Gen  = nullptr;
    TH2D* templatePtV0Rec  = nullptr;
    TH2D* templateZV0Gen   = nullptr;
    TH2D* templateZV0Rec   = nullptr;
  public:
    // Unfolding settings
    int doSmoothing = 0, nIterations = 1;
    int minIteration = 1, maxIteration = 1;
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovToy;

    // General settings
    int train;
    string inputFileName, outputFileName, responseFileName;
    string hadron;

    double ptjetminGen = -1., ptjetmaxGen = -1., ptjetminRec = -1., ptjetmaxRec = -1.;
    double ptv0minGen = -1., ptv0maxGen = -1., ptv0minRec = -1., ptv0maxRec = -1.;
    double zv0minGen = -1., zv0maxGen = -1., zv0minRec = -1., zv0maxRec = -1.;
    double etamin = -1., etamax = -1.;
    double binwidthptjet = -1., binwidthptv0 = -1., binwidthzv0 = -1.;
    double ptjetminProjection = -1., ptjetmaxProjection = -1.;

    bool saveFigs = true, drawText = false;
    bool doTrivialClosureTest = true;

    // Using private methods
    string getNameFromPtJetGen(string prefix, string suffix) {
      return getNameFromPtJet(prefix, ptjetminGen, ptjetmaxGen, suffix);
    }
    string getNameFromPtJetRec(string prefix, string suffix) {
      return getNameFromPtJet(prefix, ptjetminRec, ptjetmaxRec, suffix);
    }
    string getNameFromPtJetProjection(string prefix, string suffix) {
      return getNameFromPtJet(prefix, ptjetminProjection, ptjetmaxProjection, suffix);
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
    unfoldingutilities::VariableType getVariableType() { return vartype; }
    string getRmHistName() {
      switch (vartype) {
        case unfoldingutilities::kJet:   return rmutilities::analysis::nameJetRm;
        case unfoldingutilities::kK0SPt: return rmutilities::analysis::nameK0SPtRm;
        case unfoldingutilities::kK0SZ:  return rmutilities::analysis::nameK0SZRm;
        default: return "";
      }
    }
    verbosityutilities::Verbosity getVerbosity() { return verbosity; }
    TH1D* getTemplateHistJetGen()  { return templatePtJetGen; }
    TH1D* getTemplateHistJetRec()  { return templatePtJetRec; }
    TH2D* getTemplateHistV0PtGen() { return templatePtV0Gen; }
    TH2D* getTemplateHistV0PtRec() { return templatePtV0Rec; }
    TH2D* getTemplateHistV0ZGen()  { return templateZV0Gen; }
    TH2D* getTemplateHistV0ZRec()  { return templateZV0Rec; }
    TH2D* getTemplateHistV0Gen()   {
      if (vartype == unfoldingutilities::kK0SPt)
        return getTemplateHistV0PtGen();
      else if (vartype == unfoldingutilities::kK0SZ)
        return getTemplateHistV0ZGen();
      else
        return nullptr;
    }
    TH2D* getTemplateHistV0Rec()   {
      if (vartype == unfoldingutilities::kK0SPt)
        return getTemplateHistV0PtRec();
      else if (vartype == unfoldingutilities::kK0SZ)
        return getTemplateHistV0ZRec();
      else
        return nullptr;
    }

    // Setters for private variables
    bool setVariableType(unfoldingutilities::VariableType p);
    bool setVerbosity(verbosityutilities::Verbosity v);
    void setTemplateHistJetGen(TH1D* hist)  { setTemplateHist(templatePtJetGen, hist); }
    void setTemplateHistJetRec(TH1D* hist)  { setTemplateHist(templatePtJetRec, hist); }
    void setTemplateHistV0PtGen(TH2D* hist) { setTemplateHist(templatePtV0Gen, hist); }
    void setTemplateHistV0PtRec(TH2D* hist) { setTemplateHist(templatePtV0Rec, hist); }
    void setTemplateHistV0ZGen(TH2D* hist)  { setTemplateHist(templateZV0Gen, hist); }
    void setTemplateHistV0ZRec(TH2D* hist)  { setTemplateHist(templateZV0Rec, hist); }

    void setTemplateHistJetGen(const int nbinsx, const double* binedges);
    void setTemplateHistJetRec(const int nbinsx, const double* binedges);
    void setTemplateHistJetGen(const double min, const double max, const double binwidth);
    void setTemplateHistJetRec(const double min, const double max, const double binwidth);

    void setTemplateHistV0PtGen(const int nbinsx, const double* binedgesx, const int nbinsy, const double* binedgesy) { setTemplateHistV0Pt(nbinsx, binedgesx, nbinsy, binedgesy, true); }
    void setTemplateHistV0PtRec(const int nbinsx, const double* binedgesx, const int nbinsy, const double* binedgesy) { setTemplateHistV0Pt(nbinsx, binedgesx, nbinsy, binedgesy, false); }
    void setTemplateHistV0PtGen(const double minx, const double maxx, const double binwidthx, const double miny, const double maxy, const double binwidthy) { setTemplateHistV0Pt(minx, maxx, binwidthx, miny, maxy, binwidthy, true); }
    void setTemplateHistV0PtRec(const double minx, const double maxx, const double binwidthx, const double miny, const double maxy, const double binwidthy) { setTemplateHistV0Pt(minx, maxx, binwidthx, miny, maxy, binwidthy, false); }

    void setTemplateHistV0ZGen(const int nbinsx, const double* binedgesx, const int nbinsy, const double* binedgesy) { setTemplateHistV0Z(nbinsx, binedgesx, nbinsy, binedgesy, true); }
    void setTemplateHistV0ZRec(const int nbinsx, const double* binedgesx, const int nbinsy, const double* binedgesy) { setTemplateHistV0Z(nbinsx, binedgesx, nbinsy, binedgesy, false); }
    void setTemplateHistV0ZGen(const double minx, const double maxx, const double binwidthx, const double miny, const double maxy, const double binwidthy) { setTemplateHistV0Z(minx, maxx, binwidthx, miny, maxy, binwidthy, true); }
    void setTemplateHistV0ZRec(const double minx, const double maxx, const double binwidthx, const double miny, const double maxy, const double binwidthy) { setTemplateHistV0Z(minx, maxx, binwidthx, miny, maxy, binwidthy, false); }

    // Utilities
    void autoTemplateHists();
    bool isHistInRange(TH1D* hist, double min, double max);
    bool isHistConsistentWithZero(TH1D* hist, double threshold = 1e-10) {
      return isHistInRange(hist, -threshold, threshold);
    }
    bool isHistConsistentWithOne(TH1D* hist, double threshold = 1e-10) {
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

void InputSettings::autoTemplateHists() {
  if (ptjetminGen < 0. || ptjetmaxGen < 0. || ptjetminRec < 0. || ptjetmaxRec < 0. || binwidthptjet < 0.) {
    ptjetminGen = 5.;
    ptjetmaxGen = 80.;
    ptjetminRec = 10.;
    ptjetmaxRec = 60.;
    binwidthptjet = 5.;
  }
  const int nbinsPtJetGen = (int)((ptjetmaxGen - ptjetminGen) / binwidthptjet);
  const int nbinsPtJetRec = (int)((ptjetmaxRec - ptjetminRec) / binwidthptjet);
  templatePtJetGen = new TH1D("templateJetPtGen", ";#it{p}_{T,jet} (GeV/c)", nbinsPtJetGen, ptjetminGen, ptjetmaxGen);
  templatePtJetRec = new TH1D("templateJetPtRec", ";#it{p}_{T,jet} (GeV/c)", nbinsPtJetRec, ptjetminRec, ptjetmaxRec);

  if (zv0minGen < 0. || zv0maxGen < 0. || zv0minRec < 0. || zv0maxRec < 0. || binwidthzv0 < 0.) {
    zv0minGen = 1e-3;
    zv0maxGen = 1. + 1e-3;
    zv0minRec = 1e-3;
    zv0maxRec = 1. + 1e-3;
    binwidthzv0 = 0.1;
  }
  const int nbinszv0Gen   = (int)((zv0maxGen - zv0minGen) / binwidthzv0);
  const int nbinszv0Rec   = (int)((zv0maxRec - zv0minRec) / binwidthzv0);
  templateZV0Gen   = new TH2D("templateV0ZGen", ";#it{z}_{V0};#it{p}_{T,jet} (GeV/c)", nbinszv0Gen, zv0minGen, zv0maxGen, nbinsPtJetGen, ptjetminGen, ptjetmaxGen);
  templateZV0Rec   = new TH2D("templateV0ZRec", ";#it{z}_{V0};#it{p}_{T,jet} (GeV/c)", nbinszv0Rec, zv0minRec, zv0maxRec, nbinsPtJetRec, ptjetminRec, ptjetmaxRec);

  if (ptv0minGen < 0. || ptv0maxGen < 0. || binwidthptv0 < 0.) {
    ptv0minGen = 1.;
    ptv0maxGen = 40.;
    binwidthptv0 = 1.;
  }
  const int nbinsptv0Gen  = (int)((ptv0maxGen - ptv0minGen) / binwidthptv0);
  templatePtV0Gen  = new TH2D("templateV0PtGen", ";#it{p}_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/c)", nbinsptv0Gen, ptv0minGen, ptv0maxGen, nbinsPtJetGen, ptjetminGen, ptjetmaxGen);

  if (ptv0minRec < 0. || ptv0maxRec < 0.) {
    const int nbinsptv0Rec  = 10;
    double ptv0RecBinEdges[nbinsptv0Rec + 1] = {1., 2., 3., 4., 5., 10., 15., 20., 25., 30., 40.};
    templatePtV0Rec  = new TH2D("templateV0PtRec", ";#it{p}_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/c)", nbinsptv0Rec, ptv0RecBinEdges, nbinsPtJetRec, ptjetminRec, ptjetmaxRec);
  } else {
    const int nbinsptv0Rec  = (int)((ptv0maxRec - ptv0minRec) / binwidthptv0);
    templatePtV0Rec  = new TH2D("templateV0PtRec", ";#it{p}_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/c)", nbinsptv0Rec, ptv0minRec, ptv0maxRec, nbinsPtJetRec, ptjetminRec, ptjetmaxRec);
  }

}

bool InputSettings::isHistInRange(TH1D* hist, double min, double max) {
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
void InputSettings::setTemplateHist(T* tmplt, T* hist) {
  if (!hist) {
    printLog("InputSettings::setTemplateHist: input histogram is null!", verbosityutilities::kErrors);
    return;
  }
  string name = hist->GetName();
  hist->SetName(TString::Format("%s_copy", name.c_str()).Data());
  printLog(TString::Format("%s: %d", name.c_str(), hist->GetNbinsX()).Data(), verbosityutilities::kDebug);
  tmplt = (T*)hist->Clone(name.c_str());
}

void InputSettings::setTemplateHistJet(const int nbins, const double* binedges, const bool generatorLevel) {
  if (generatorLevel) {
    TH1D* hist = new TH1D("templateJetPtGen", ";#it{p}_{T,jet} (GeV/c)", nbins, binedges);
    setTemplateHistJetGen(hist);
  } else {
    TH1D* hist = new TH1D("templateJetPtRec", ";#it{p}_{T,jet} (GeV/c)", nbins, binedges);
    setTemplateHistJetRec(hist);
  }
}

void InputSettings::setTemplateHistJet(const double min, const double max, const double binwidth, const bool generatorLevel) {
  int nbins = (int)((max - min) / binwidth);
  double edges[nbins + 1];
  for (int i = 0; i <= nbins; i++)
    edges[i] = min + i * binwidth;

  setTemplateHistJet(nbins, edges, generatorLevel);
}

void InputSettings::setTemplateHistV0Pt(const int nbinsx, const double* binedgesx, const int nbinsy, const double* binedgesy, const bool generatorLevel) {
  if (generatorLevel) {
    TH2D* hist = new TH2D("templateV0PtGen", ";#it{p}_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/c)", nbinsx, binedgesx, nbinsy, binedgesy);
    setTemplateHistV0PtGen(hist);
  } else {
    TH2D* hist = new TH2D("templateV0PtRec", ";#it{p}_{T,V0} (GeV/c);#it{p}_{T,jet} (GeV/c)", nbinsx, binedgesx, nbinsy, binedgesy);
    setTemplateHistV0PtRec(hist);
  }
}

void InputSettings::setTemplateHistV0Pt(const double xmin, const double xmax, const double binwidthx, const double ymin, const double ymax, const double binwidthy, const bool generatorLevel) {
  int nbinsx = (int)((xmax - xmin) / binwidthx);
  double edgesx[nbinsx + 1];
  for (int i = 0; i <= nbinsx; i++)
    edgesx[i] = xmin + i * binwidthx;

  int nbinsy = (int)((ymax - ymin) / binwidthy);
  double edgesy[nbinsy + 1];
  for (int i = 0; i <= nbinsy; i++)
    edgesy[i] = ymin + i * binwidthy;

  setTemplateHistV0Pt(nbinsx, edgesx, nbinsy, edgesy, generatorLevel);
}

void InputSettings::setTemplateHistV0Z(const int nbinsx, const double* binedgesx, const int nbinsy, const double* binedgesy, const bool generatorLevel) {
  if (generatorLevel) {
    TH2D* hist = new TH2D("templateV0ZGen", ";#it{z}_{V0};#it{p}_{T,jet} (GeV/c)", nbinsx, binedgesx, nbinsy, binedgesy);
    setTemplateHistV0ZGen(hist);
  } else {
    TH2D* hist = new TH2D("templateV0ZRec", ";#it{z}_{V0};#it{p}_{T,jet} (GeV/c)", nbinsx, binedgesx, nbinsy, binedgesy);
    setTemplateHistV0ZRec(hist);
  }
}

void InputSettings::setTemplateHistV0Z(const double minx, const double maxx, const double binwidthx, const double miny, const double maxy, const double binwidthy, const bool generatorLevel) {
  int nbinsx = (int)((maxx - minx) / binwidthx);
  double edgesx[nbinsx + 1];
  for (int i = 0; i <= nbinsx; i++)
    edgesx[i] = minx + i * binwidthx;

  int nbinsy = (int)((maxy - miny) / binwidthy);
  double edgesy[nbinsy + 1];
  for (int i = 0; i <= nbinsy; i++)
    edgesy[i] = miny + i * binwidthy;

  setTemplateHistV0Z(nbinsx, edgesx, nbinsy, edgesy, generatorLevel);
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

bool InputSettings::setVariableType(unfoldingutilities::VariableType p) {
  if (!unfoldingutilities::is_valid(p)) {
    printLog("InputSettings::setVariableType() Error: invalid plot type", verbosityutilities::kErrors);
    return false;
  }
  vartype = p;
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

// NB: K0S RM is filled (gen, rec), but the jet RM (rec, gen)!
// NB: The RooUnfoldResponse is filled (rec, gen) always!
void FillFromRmJets(InputSettings& inputs, TH2D* responseMatrix, TH1D* hRec, TH1D* hGen, TH1D* hMiss, TH1D* hGenAndRec, TH1D* hFake, RooUnfoldResponse* response = nullptr) {
  // Fill the distributions from the response matrix
  // If the RooUnfoldResponse is given, fill it as well
  array<double, 2> ptjetRecBinEdges = {hRec->GetXaxis()->GetXmin(), hRec->GetXaxis()->GetXmax()};
  array<double, 2> ptjetGenBinEdges = {hGen->GetXaxis()->GetXmin(), hGen->GetXaxis()->GetXmax()};

  for (int xBin = 0; xBin <= 1 + responseMatrix->GetNbinsX(); xBin++) {
    for (int yBin = 0; yBin <= 1 + responseMatrix->GetNbinsY(); yBin++) {
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
        hGenAndRec->Fill(ptjetGen, binContent);
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
}

void FillFromRmK0S(InputSettings& inputs, THnSparseD* responseMatrix, TH2D* hRec, TH2D* hGen, TH2D* hMiss, TH2D* hGenAndRec, TH2D* hFake, RooUnfoldResponse* response = nullptr) {
  // Fill the distributions from the response matrix
  // If the RooUnfoldResponse is given, fill it as well
  array<double, 2> k0sGenBinEdges   = {hGen->GetXaxis()->GetXmin(), hGen->GetXaxis()->GetXmax()};
  array<double, 2> ptjetGenBinEdges = {hGen->GetYaxis()->GetXmin(), hGen->GetYaxis()->GetXmax()};
  array<double, 2> k0sRecBinEdges   = {hRec->GetXaxis()->GetXmin(), hRec->GetXaxis()->GetXmax()};
  array<double, 2> ptjetRecBinEdges = {hRec->GetYaxis()->GetXmin(), hRec->GetYaxis()->GetXmax()};

  int* coord = new int[responseMatrix->GetNdimensions()]; //Carries the bin coordinates
  for (int iBin = 0; iBin <= responseMatrix->GetNbins(); iBin++) { // FIXME: How to get all overflow bins?
    double binContent = responseMatrix->GetBinContent(iBin, coord);
    double ptjetGen   = responseMatrix->GetAxis(rmutilities::analysis::axisK0SRmPtJetGen)->GetBinCenter(coord[rmutilities::analysis::axisK0SRmPtJetGen]);
    double ptK0SGen   = responseMatrix->GetAxis(rmutilities::analysis::axisK0SRmK0SGen)->GetBinCenter(coord[rmutilities::analysis::axisK0SRmK0SGen]);
    double ptjetRec   = responseMatrix->GetAxis(rmutilities::analysis::axisK0SRmPtJetRec)->GetBinCenter(coord[rmutilities::analysis::axisK0SRmPtJetRec]);
    double ptK0SRec   = responseMatrix->GetAxis(rmutilities::analysis::axisK0SRmK0SRec)->GetBinCenter(coord[rmutilities::analysis::axisK0SRmK0SRec]);

    bool isAcceptedGen = inputs.isVarInRange(ptjetGen, ptjetGenBinEdges) && inputs.isVarInRange(ptK0SGen, k0sGenBinEdges);
    bool isAcceptedRec = inputs.isVarInRange(ptjetRec, ptjetRecBinEdges) && inputs.isVarInRange(ptK0SRec, k0sRecBinEdges);

    if (isAcceptedRec)
      hRec->Fill(ptK0SRec, ptjetRec, binContent);
    if (isAcceptedGen)
      hGen->Fill(ptK0SGen, ptjetGen, binContent);

    if (isAcceptedRec && isAcceptedGen) {
      hGenAndRec->Fill(ptK0SGen, ptjetGen, binContent);
      if (response) response->Fill(ptK0SRec, ptjetRec, ptK0SGen, ptjetGen, binContent);
    } else if (!isAcceptedRec && isAcceptedGen) {
      hMiss->Fill(ptK0SGen, ptjetGen, binContent);
      if (response) response->Miss(ptK0SGen, ptjetGen, binContent);
    } else if (isAcceptedRec && !isAcceptedGen) {
      hFake->Fill(ptK0SRec, ptjetRec, binContent);
      if (response) response->Fake(ptK0SRec, ptjetRec, binContent);
    }
  }
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
  TH1D* hRec = (TH1D*)inputs.getTemplateHistJetRec()->Clone(rmutilities::unfolding::nameRecJets.c_str());
  TH1D* hGen = (TH1D*)inputs.getTemplateHistJetGen()->Clone(rmutilities::unfolding::nameGenJets.c_str());
  hRec->SetTitle("Rec");
  hGen->SetTitle("Gen");
  TH1D* hMiss      = (TH1D*)hGen->Clone(rmutilities::unfolding::nameMissJets.c_str());
  TH1D* hGenAndRec = (TH1D*)hGen->Clone(rmutilities::unfolding::nameGenAndRecJets.c_str());
  TH1D* hFake      = (TH1D*)hRec->Clone(rmutilities::unfolding::nameFakeJets.c_str());
  return std::array<TH1D*, 5>{hRec, hGen, hMiss, hGenAndRec, hFake};
}

array<TH2D*, 5> CreateDistributionHistsK0S(InputSettings& inputs, THnSparseD* responseMatrix) {
  inputs.printLog("Creating histograms for K0S.", verbosityutilities::kDebug);
  string nRec, nGen, nMiss, nGenAndRec, nFake, nAxes;

  if (inputs.getVariableType() == unfoldingutilities::VariableType::kK0SPt) {
    nRec        = rmutilities::unfolding::nameRecK0SPt;
    nGen        = rmutilities::unfolding::nameGenK0SPt;
    nMiss       = rmutilities::unfolding::nameMissK0SPt;
    nGenAndRec  = rmutilities::unfolding::nameGenAndRecK0SPt;
    nFake       = rmutilities::unfolding::nameFakeK0SPt;
    nAxes       = "p_{T,K^{0}_{S}};p_{T,jet}";
  } else if (inputs.getVariableType() == unfoldingutilities::VariableType::kK0SZ) {
    nRec        = rmutilities::unfolding::nameRecK0SZ;
    nGen        = rmutilities::unfolding::nameGenK0SZ;
    nMiss       = rmutilities::unfolding::nameMissK0SZ;
    nGenAndRec  = rmutilities::unfolding::nameGenAndRecK0SZ;
    nFake       = rmutilities::unfolding::nameFakeK0SZ;
    nAxes       = "z_{K^{0}_{S}};p_{T,jet}";
  } else {
    inputs.printLog("CreateDistributionHistsK0S() Error: invalid plot type " + unfoldingutilities::to_string(inputs.getVariableType()), verbosityutilities::kErrors);
    return std::array<TH2D*, 5>{nullptr, nullptr, nullptr, nullptr, nullptr};
  }

  TH2D* hRec = (TH2D*)inputs.getTemplateHistV0Rec()->Clone(nRec.c_str());
  TH2D* hGen = (TH2D*)inputs.getTemplateHistV0Gen()->Clone(nGen.c_str());
  hRec->SetTitle(("Rec;" + nAxes).c_str());
  hGen->SetTitle(("Gen;" + nAxes).c_str());
  TH2D* hMiss      = (TH2D*)hGen->Clone(nMiss.c_str());
  TH2D* hGenAndRec = (TH2D*)hGen->Clone(nGenAndRec.c_str());
  TH2D* hFake      = (TH2D*)hRec->Clone(nFake.c_str());
  return std::array<TH2D*, 5>{hRec, hGen, hMiss, hGenAndRec, hFake};
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
  inputs.setVariableType(unfoldingutilities::VariableType::kJet);

  TH2D* responseMatrix = LoadResponseMatrix<TH2D*>(inputs);
  responseMatrix->SetName(rmutilities::unfolding::nameRmJets.c_str());
  responseMatrix->SetTitle(rmutilities::unfolding::nameRmJets.c_str());

  array<TH1D*, 5> hists = CreateDistributionHistsJet(inputs, responseMatrix);
  TH1D* hRec       = hists[0];
  TH1D* hGen       = hists[1];
  TH1D* hMiss      = hists[2];
  TH1D* hGenAndRec = hists[3];
  TH1D* hFake      = hists[4];

  inputs.printLog("Filling response.", verbosityutilities::kDebug);
  RooUnfoldResponse* response = new RooUnfoldResponse(rmutilities::unfolding::nameResponseJets.c_str(), rmutilities::unfolding::nameResponseJets.c_str());
  response->Setup(hRec, hGen);
  FillFromRmJets(inputs, responseMatrix, hRec, hGen, hMiss, hGenAndRec, hFake, response);

  TH1D* hKinEff = (TH1D*)hGenAndRec->Clone(rmutilities::unfolding::nameKinEffJets.c_str());
  hKinEff->Divide(hGen);

  inputs.printLog("Saving objects to output file: " + inputs.outputFileName, verbosityutilities::kDebug);
  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(responseMatrix);
  inputs.writeOutputsToFile(std::vector<TH1D*>({hRec, hGen, hMiss, hFake, hGenAndRec, hKinEff}));
  return;
}

void CreateResponseV0(InputSettings& inputs) {
  inputs.printLog("Creating response for " + unfoldingutilities::to_string(inputs.getVariableType()), verbosityutilities::kInfo);

  string nameRm = rmutilities::unfolding::nameRmK0SPt;
  string nameResponse = rmutilities::unfolding::nameResponseK0SPt;
  string nameKinEff = rmutilities::unfolding::nameKinEffK0SPt;
  if (inputs.getVariableType() == unfoldingutilities::VariableType::kK0SZ) {
    nameRm = rmutilities::unfolding::nameRmK0SZ;
    nameResponse = rmutilities::unfolding::nameResponseK0SZ;
    nameKinEff = rmutilities::unfolding::nameKinEffK0SZ;
  }

  THnSparseD* responseMatrix = LoadResponseMatrix<THnSparseD*>(inputs);
  responseMatrix->SetName(nameRm.c_str());
  responseMatrix->SetTitle(nameRm.c_str());

  array<TH2D*, 5> hists = CreateDistributionHistsK0S(inputs, responseMatrix);
  TH2D* hRec       = hists[0];
  TH2D* hGen       = hists[1];
  TH2D* hMiss      = hists[2];
  TH2D* hGenAndRec = hists[3];
  TH2D* hFake      = hists[4];

  inputs.printLog("Filling response.", verbosityutilities::kDebug);
  RooUnfoldResponse* response = new RooUnfoldResponse(nameResponse.c_str(), nameResponse.c_str());
  response->Setup(hRec, hGen);
  FillFromRmK0S(inputs, responseMatrix, hRec, hGen, hMiss, hGenAndRec, hFake, response);

  TH2D* hKinEff = (TH2D*)hGenAndRec->Clone(nameKinEff.c_str());
  hKinEff->Divide(hGen);

  inputs.printLog("Saving objects to output file: " + inputs.outputFileName, verbosityutilities::kDebug);
  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(responseMatrix);
  inputs.writeOutputsToFile(std::vector<TH2D*>({hRec, hGen, hMiss, hGenAndRec, hFake, hKinEff}));
}

void CreateResponseV0Pt(InputSettings& inputs) {
  inputs.setVariableType(unfoldingutilities::VariableType::kK0SPt);
  CreateResponseV0(inputs);
}

void CreateResponseV0Z(InputSettings& inputs) {
  inputs.setVariableType(unfoldingutilities::VariableType::kK0SZ);
  CreateResponseV0(inputs);
}

// ----------------------------------------------------------

void DoUnfoldingJets(InputSettings& inputs, int nIterations) {
  inputs.printLog("Doing unfolding for jets with " + to_string(nIterations) + " iterations.", verbosityutilities::kInfo);
  inputs.setVariableType(unfoldingutilities::VariableType::kJet);

  inputs.printLog("Getting response from file: " + inputs.responseFileName + "\nGetting test distributions from file: " + inputs.inputFileName, verbosityutilities::kDebug);
  TFile* responseFile = TFile::Open(inputs.responseFileName.c_str(), "UPDATE");
  if (!responseFile) {
    inputs.printLog("Error: could not open file " + inputs.responseFileName, verbosityutilities::kErrors);
    return;
  }

  inputs.printLog("Retrieving response and histograms from training.", verbosityutilities::kDebug);
  TH2D* trainingRm = (TH2D*)responseFile->Get(rmutilities::unfolding::nameRmJets.c_str());
  TH1D* trainingRec = (TH1D*)responseFile->Get(rmutilities::unfolding::nameRecJets.c_str());
  TH1D* trainingGen = (TH1D*)responseFile->Get(rmutilities::unfolding::nameGenJets.c_str());
  TH1D* trainingFake = (TH1D*)responseFile->Get(rmutilities::unfolding::nameFakeJets.c_str());
  TH1D* trainingMiss = (TH1D*)responseFile->Get(rmutilities::unfolding::nameMissJets.c_str());
  TH1D* trainingKinEff = (TH1D*)responseFile->Get(rmutilities::unfolding::nameKinEffJets.c_str());
  TH1D* trainingGenAndRec = (TH1D*)responseFile->Get(rmutilities::unfolding::nameGenAndRecJets.c_str());
  RooUnfoldResponse* response = (RooUnfoldResponse*)responseFile->Get(rmutilities::unfolding::nameResponseJets.c_str());
  if (!trainingRm)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameRmJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
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
  if (!trainingGenAndRec)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameGenAndRecJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!response)
    inputs.printLog("Error: could not find " + rmutilities::unfolding::nameResponseJets + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingRm || !trainingRec || !trainingGen || !trainingFake || !trainingMiss || !trainingKinEff || !trainingGenAndRec || !response)
    return;

  // Change the names to distinguish these from the test histograms when writing to file
  trainingRm->SetName(rmutilities::training::nameRmJets.c_str());
  trainingRec->SetName(rmutilities::training::nameRecJets.c_str());
  trainingGen->SetName(rmutilities::training::nameGenJets.c_str());
  trainingFake->SetName(rmutilities::training::nameFakeJets.c_str());
  trainingMiss->SetName(rmutilities::training::nameMissJets.c_str());
  trainingKinEff->SetName(rmutilities::training::nameKinEffJets.c_str());
  trainingGenAndRec->SetName(rmutilities::training::nameGenAndRecJets.c_str());

  // Create the test histograms. These must have the same binning as the training histograms
  // For trivial closure test, they are copies, otherwise they are filled with independent data
  TH2D* testRm = (TH2D*)trainingRm->Clone(rmutilities::testing::nameRmJets.c_str());
  TH1D* testRec = (TH1D*)trainingRec->Clone(rmutilities::testing::nameRecJets.c_str());
  TH1D* testGen = (TH1D*)trainingGen->Clone(rmutilities::testing::nameGenJets.c_str());
  TH1D* testFake = (TH1D*)trainingFake->Clone(rmutilities::testing::nameFakeJets.c_str());
  TH1D* testMiss = (TH1D*)trainingMiss->Clone(rmutilities::testing::nameMissJets.c_str());
  TH1D* testGenAndRec = (TH1D*)trainingGenAndRec->Clone(rmutilities::testing::nameGenAndRecJets.c_str());

  if (inputs.doTrivialClosureTest) {
    inputs.printLog("Trivial closure test: using training histograms as test histograms.", verbosityutilities::kInfo);
  } else {
    inputs.printLog("Statistically independent closure test: getting test histograms from test response matrix.", verbosityutilities::kInfo);
    testRec->Reset();
    testGen->Reset();
    testFake->Reset();
    testMiss->Reset();
    testGenAndRec->Reset();

    TFile* testFile = TFile::Open(inputs.inputFileName.c_str(), "READ");
    if (!testFile) {
      inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
      return;
    }
    testRm = (TH2D*)testFile->Get(inputs.getRmHistName().c_str());
    FillFromRmJets(inputs, testRm, testRec, testGen, testMiss, testGenAndRec, testFake);
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
  inputs.writeOutputToFile(trainingRm);
  inputs.writeOutputsToFile(std::vector<TH1D*>{trainingRec, trainingGen, trainingFake, trainingMiss, trainingGenAndRec, trainingKinEff});
  inputs.writeOutputToFile(&ruBayes);

  TH1D* testKinEff = (TH1D*)testGenAndRec->Clone(rmutilities::testing::nameKinEffJets.c_str());
  testKinEff->Divide(testGen);
  inputs.writeOutputToFile(testRm);
  inputs.writeOutputsToFile(std::vector<TH1D*>{testRec, testGen, testFake, testMiss, testGenAndRec, testKinEff, unfolded, refolded});
  inputs.writeOutputsToFile(std::vector<TH2D*>{covMatrix, pearson});
}

void DoUnfoldingV0(InputSettings& inputs, int nIterations) {
  inputs.printLog("Doing unfolding for " + unfoldingutilities::to_string(inputs.getVariableType()) + " with " + to_string(nIterations) + " iterations.", verbosityutilities::kInfo);

  string nameRm, nameRec, nameGen, nameFake, nameMiss, nameGenAndRec, nameKinEff, nameResponse, ruBayesName, unfoldedName, refoldedName, covMatrixName, pearsonName;
  if (inputs.getVariableType() == unfoldingutilities::VariableType::kK0SPt) {
    nameRm        = rmutilities::unfolding::nameRmK0SPt;
    nameRec       = rmutilities::unfolding::nameRecK0SPt;
    nameGen       = rmutilities::unfolding::nameGenK0SPt;
    nameFake      = rmutilities::unfolding::nameFakeK0SPt;
    nameMiss      = rmutilities::unfolding::nameMissK0SPt;
    nameGenAndRec = rmutilities::unfolding::nameGenAndRecK0SPt;
    nameKinEff    = rmutilities::unfolding::nameKinEffK0SPt;
    nameResponse  = rmutilities::unfolding::nameResponseK0SPt;
    ruBayesName   = rmutilities::unfolding::nameRooUnfoldBayesK0SPt + to_string(nIterations);
    unfoldedName  = rmutilities::unfolding::nameUnfoldedK0SPt + to_string(nIterations);
    refoldedName  = rmutilities::unfolding::nameRefoldedK0SPt + to_string(nIterations);
    covMatrixName = rmutilities::unfolding::nameCovMatrixK0SPt + to_string(nIterations);
    pearsonName   = rmutilities::unfolding::namePearsonK0SPt + to_string(nIterations);
  } else if (inputs.getVariableType() == unfoldingutilities::VariableType::kK0SZ) {
    nameRm        = rmutilities::unfolding::nameRmK0SZ;
    nameRec       = rmutilities::unfolding::nameRecK0SZ;
    nameGen       = rmutilities::unfolding::nameGenK0SZ;
    nameFake      = rmutilities::unfolding::nameFakeK0SZ;
    nameMiss      = rmutilities::unfolding::nameMissK0SZ;
    nameGenAndRec = rmutilities::unfolding::nameGenAndRecK0SZ;
    nameKinEff    = rmutilities::unfolding::nameKinEffK0SZ;
    nameResponse  = rmutilities::unfolding::nameResponseK0SZ;
    ruBayesName   = rmutilities::unfolding::nameRooUnfoldBayesK0SZ + to_string(nIterations);
    unfoldedName  = rmutilities::unfolding::nameUnfoldedK0SZ + to_string(nIterations);
    refoldedName  = rmutilities::unfolding::nameRefoldedK0SZ + to_string(nIterations);
    covMatrixName = rmutilities::unfolding::nameCovMatrixK0SZ + to_string(nIterations);
    pearsonName   = rmutilities::unfolding::namePearsonK0SZ + to_string(nIterations);
  } else {
    inputs.printLog("DoUnfoldingV0() Error: invalid plot type " + unfoldingutilities::to_string(inputs.getVariableType()), verbosityutilities::kErrors);
    return;
  }

  inputs.printLog("Getting response from file: " + inputs.responseFileName + "\nGetting test distributions from file: " + inputs.inputFileName, verbosityutilities::kDebug);
  TFile* responseFile = TFile::Open(inputs.responseFileName.c_str(), "READ");
  if (!responseFile) {
    inputs.printLog("Error: could not open file " + inputs.responseFileName, verbosityutilities::kErrors);
    return;
  }

  inputs.printLog("Retrieving response and histograms from training.", verbosityutilities::kDebug);
  THnSparseD* trainingRm = (THnSparseD*)responseFile->Get(nameRm.c_str());
  TH2D* trainingRec = (TH2D*)responseFile->Get(nameRec.c_str());
  TH2D* trainingGen = (TH2D*)responseFile->Get(nameGen.c_str());
  TH2D* trainingFake = (TH2D*)responseFile->Get(nameFake.c_str());
  TH2D* trainingMiss = (TH2D*)responseFile->Get(nameMiss.c_str());
  TH2D* trainingKinEff = (TH2D*)responseFile->Get(nameKinEff.c_str());
  TH2D* trainingGenAndRec = (TH2D*)responseFile->Get(nameGenAndRec.c_str());
  RooUnfoldResponse* response = (RooUnfoldResponse*)responseFile->Get(nameResponse.c_str());

  if (!trainingRm)
    inputs.printLog("Error: could not find " + nameRm + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
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
  if (!trainingGenAndRec)
    inputs.printLog("Error: could not find " + nameGenAndRec + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!response)
    inputs.printLog("Error: could not find " + nameResponse + " in file " + inputs.responseFileName, verbosityutilities::kErrors);
  if (!trainingRm || !trainingRec || !trainingGen || !trainingFake || !trainingMiss || !trainingKinEff || !trainingGenAndRec || !response)
    return;

  // Change the names to distinguish these from the test histograms when writing to file
  string prefTraining = rmutilities::unfolding::prefixTraining;
  trainingRm->SetName((prefTraining + nameRm).c_str());
  trainingRec->SetName((prefTraining + nameRec).c_str());
  trainingGen->SetName((prefTraining + nameGen).c_str());
  trainingFake->SetName((prefTraining + nameFake).c_str());
  trainingMiss->SetName((prefTraining + nameMiss).c_str());
  trainingKinEff->SetName((prefTraining + nameKinEff).c_str());
  trainingGenAndRec->SetName((prefTraining + nameGenAndRec).c_str());

  // Create the test histograms. These must have the same binning as the training histograms
  // For trivial closure test, they are copies, otherwise they are filled with independent data
  string prefTest = rmutilities::unfolding::prefixTest;
  THnSparseD* testRm = (THnSparseD*)trainingRm->Clone((prefTest + nameRm).c_str());
  TH2D* testRec = (TH2D*)trainingRec->Clone((prefTest + nameRec).c_str());
  TH2D* testGen = (TH2D*)trainingGen->Clone((prefTest + nameGen).c_str());
  TH2D* testFake = (TH2D*)trainingFake->Clone((prefTest + nameFake).c_str());
  TH2D* testMiss = (TH2D*)trainingMiss->Clone((prefTest + nameMiss).c_str());
  TH2D* testGenAndRec = (TH2D*)trainingGenAndRec->Clone((prefTest + nameGenAndRec).c_str());

  if (inputs.doTrivialClosureTest) {
    inputs.printLog("Trivial closure test: using training histograms as test histograms.", verbosityutilities::kInfo);
  } else {
    inputs.printLog("Statistically independent closure test: getting test histograms from test response matrix.", verbosityutilities::kInfo);
    testRec->Reset();
    testGen->Reset();
    testFake->Reset();
    testMiss->Reset();
    testGenAndRec->Reset();

    TFile* testFile = TFile::Open(inputs.inputFileName.c_str(), "READ");
    if (!testFile) {
      inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
      return;
    }
    testRm = LoadResponseMatrix<THnSparseD*>(inputs);
    FillFromRmK0S(inputs, testRm, testRec, testGen, testMiss, testGenAndRec, testFake);
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

  TH2D* testKinEff = (TH2D*)testGenAndRec->Clone((prefTest + nameKinEff).c_str());
  testKinEff->Divide(testGen);

  inputs.writeOutputToFile(response);
  inputs.writeOutputToFile(trainingRm);
  inputs.writeOutputsToFile(std::vector<TH2D*>{trainingRec, trainingGen, trainingFake, trainingMiss, trainingKinEff, trainingGenAndRec});
  inputs.writeOutputToFile(&ruBayes);
  inputs.writeOutputsToFile(std::vector<TH2D*>{unfolded, refolded});
  inputs.writeOutputToFile(testRm);
  inputs.writeOutputsToFile(std::vector<TH2D*>{testRec, testGen, testFake, testMiss, testKinEff, testGenAndRec});
  inputs.writeOutputsToFile(std::vector<TH2D*>{covMatrix, pearson});
}

void DoUnfoldingV0Pt(InputSettings& inputs, int nIterations) {
  inputs.setVariableType(unfoldingutilities::VariableType::kK0SPt);
  DoUnfoldingV0(inputs, nIterations);
}

void DoUnfoldingV0Z(InputSettings& inputs, int nIterations) {
  inputs.setVariableType(unfoldingutilities::VariableType::kK0SZ);
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
  inputs.printLog(TString::Format("Making plots for %s closure test for iteration = %d, ptjet = %.f-%.f", unfoldingutilities::to_string(inputs.getVariableType()).c_str(), nIteration, inputs.ptjetminProjection, inputs.ptjetmaxProjection).Data(), verbosityutilities::kInfo);

  string drawOption = inputs.drawText ? "hist text" : "hist";
  string canvasName = isUnfolded ? "unfolded" : "refolded";
  canvasName += unfoldingutilities::to_string(inputs.getVariableType()) + "-iteration" + to_string(nIteration) + ".pdf";
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
array<TH1D*, 3> MakeClosureTestHists(const T* hRooUnfold, const T* hAnalysis, double ptjetmin, double ptjetmax, bool isUnfolded, unfoldingutilities::VariableType plot) {
  string nRooUnfold = isUnfolded ? "Unfolded" : "Refolded";
  string nAnalysis = isUnfolded ? "Generated" : "Reconstructed";
  string nObservable;
  TH1D* hRoo;
  TH1D* hAna;
  array<int, 2> ptjetBins;

  switch (plot) {
    case unfoldingutilities::VariableType::kJet:
      nObservable = "#it{p}_{T,jet} (GeV/c)";
      hRoo = (TH1D*)hRooUnfold->Clone("hRoo");
      hAna = (TH1D*)hAnalysis->Clone("hAna");
      ptjetBins = histutils::getProjectionBins(hRooUnfold->GetXaxis(), ptjetmin, ptjetmax);
      hRoo->GetXaxis()->SetRange(ptjetBins[0], ptjetBins[1]);
      hAna->GetXaxis()->SetRange(ptjetBins[0], ptjetBins[1]);
      break;
    case unfoldingutilities::VariableType::kK0SPt:
      nObservable = "#it{p}_{T,K^{0}_{S}} (GeV/c)";
      ptjetBins = histutils::getProjectionBins(hRooUnfold->GetYaxis(), ptjetmin, ptjetmax);
      hRoo = (TH1D*)hRooUnfold->ProjectionX("hRoo", ptjetBins[0], ptjetBins[1]);
      hAna = (TH1D*)hAnalysis->ProjectionX("hAna", ptjetBins[0], ptjetBins[1]);
      break;
    case unfoldingutilities::VariableType::kK0SZ:
      nObservable = "#it{z}_{K^{0}_{S}}";
      ptjetBins = histutils::getProjectionBins(hRooUnfold->GetYaxis(), ptjetmin, ptjetmax);
      hRoo = (TH1D*)hRooUnfold->ProjectionX("hRoo", ptjetBins[0], ptjetBins[1]);
      hAna = (TH1D*)hAnalysis->ProjectionX("hAna", ptjetBins[0], ptjetBins[1]);
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

  // Self-normalise the distributions, as the training and test distributions may have different integrals
  // RooUnfold does not change the total integral of the distribution
  TH1D* hRoo = (TH1D*)hRooUnfold->Clone("hRoo");
  TH1D* hAna = (TH1D*)hAnalysis->Clone("hAna");

  hRoo->Scale(1. / hRoo->Integral(), "width");
  hAna->Scale(1. / hAna->Integral(), "width");

  TH1D* hDiff = (TH1D*)hRoo->Clone("hDiff");
  hDiff->Add(hAna, -1);
  hDiff->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
  hDiff->SetTitle((nRooUnfold + " - " + nAnalysis).c_str());

  TH1D* hRelDiff = (TH1D*)hRoo->Clone("hRelDiff");
  hRelDiff->Add(hAna, -1);
  hRelDiff->Divide(hAna);
  hRelDiff->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
  hRelDiff->SetTitle(("#frac{" + nRooUnfold + " - " + nAnalysis + "}{" + nAnalysis + "}").c_str());

  TH1D* hRatio = (TH1D*)hRoo->Clone("hRatio");
  hRatio->Divide(hAna);
  hRatio->GetXaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
  hRatio->SetTitle(("#frac{" + nRooUnfold + "}{" + nAnalysis + "}").c_str());

  return {hDiff, hRelDiff, hRatio};
}

array<TH1D*, 3> MakeClosureTestHistsV0(InputSettings& inputs, const TH2D* hRooUnfold, const TH2D* hAnalysis, bool isUnfolded) {
  string nRooUnfold = isUnfolded ? "Unfolded" : "Refolded";
  string nAnalysis = isUnfolded ? "Generated" : "Reconstructed";

  string nObservable = "#it{p}_{T,K^{0}_{S}} (GeV/c)";
  if (inputs.getVariableType() == unfoldingutilities::VariableType::kK0SZ)
    nObservable = "#it{z}_{K^{0}_{S}}";

  // The binning of the two histograms is the same by construction
  array<int, 2> ptjetBins = histutils::getProjectionBins(hRooUnfold->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
  TH1D* hRooUnfoldProj = (TH1D*)hRooUnfold->ProjectionX("hRooUnfoldProj", ptjetBins[0], ptjetBins[1]);
  TH1D* hAnalysisProj  = (TH1D*)hAnalysis->ProjectionX("hAnalysisProj", ptjetBins[0], ptjetBins[1]);

  // Self-normalise the distributions, as the training and test distributions may have different integrals
  // RooUnfold does not change the total integral of the distribution
  hRooUnfoldProj->Scale(1. / hRooUnfoldProj->Integral(), "width");
  hAnalysisProj->Scale(1. / hAnalysisProj->Integral(), "width");

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
  inputs.setVariableType(unfoldingutilities::VariableType::kJet);

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
  switch (inputs.getVariableType()) {
    case unfoldingutilities::VariableType::kK0SPt:
      inputs.printLog("Doing closure test for K0S pt.", verbosityutilities::kInfo);
      nameUnfolded = rmutilities::unfolding::nameUnfoldedK0SPt + to_string(nIteration);
      nameRefolded = rmutilities::unfolding::nameRefoldedK0SPt + to_string(nIteration);
      nameGen      = rmutilities::testing::nameGenK0SPt;
      nameRec      = rmutilities::testing::nameRecK0SPt;
      break;
      case unfoldingutilities::VariableType::kK0SZ:
      inputs.printLog("Doing closure test for K0S z.", verbosityutilities::kInfo);
      nameUnfolded = rmutilities::unfolding::nameUnfoldedK0SZ + to_string(nIteration);
      nameRefolded = rmutilities::unfolding::nameRefoldedK0SZ + to_string(nIteration);
      nameGen      = rmutilities::testing::nameGenK0SZ;
      nameRec      = rmutilities::testing::nameRecK0SZ;
      break;
    default:
      inputs.printLog("DoClosureTestV0() Error: invalid plot type " + unfoldingutilities::to_string(inputs.getVariableType()), verbosityutilities::kErrors);
      return;
  }

  inputs.printLog("Retrieving histograms " + nameUnfolded + ", " + nameRefolded + ", " + nameGen + ", " + nameRec + " from file: " + inputs.inputFileName, verbosityutilities::kDebug);
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
  if (inputs.passVerbosityCheck(verbosityutilities::kDebug)) {
    inputs.printLog("Unfolded: (nx, ny) = (" + to_string(unfolded->GetNbinsX()) + ", " + to_string(unfolded->GetNbinsY()) + ")", verbosityutilities::kDebug);
    inputs.printLog("Refolded: (nx, ny) = (" + to_string(refolded->GetNbinsX()) + ", " + to_string(refolded->GetNbinsY()) + ")", verbosityutilities::kDebug);
    inputs.printLog("Gen: (nx, ny) = (" + to_string(testGen->GetNbinsX()) + ", " + to_string(testGen->GetNbinsY()) + ")", verbosityutilities::kDebug);
    inputs.printLog("Rec: (nx, ny) = (" + to_string(testRec->GetNbinsX()) + ", " + to_string(testRec->GetNbinsY()) + ")", verbosityutilities::kDebug);

    for (const auto& h : hUnfGen)
      inputs.printLog(TString::Format("%s: nx = %d", h->GetName(), h->GetNbinsX()).Data(), verbosityutilities::kDebug);
    for (const auto& h : hRefRec)
      inputs.printLog(TString::Format("%s: nx = %d", h->GetName(), h->GetNbinsX()).Data(), verbosityutilities::kDebug);
  }
  MakeClosureTestPlots(inputs, hUnfGen, nIteration, isUnfolded);
  MakeClosureTestPlots(inputs, hRefRec, nIteration, !isUnfolded);
  if (inputs.doTrivialClosureTest) {
    inputs.printLog(TrivialClosureTestSummary(inputs, hUnfGen, isUnfolded), verbosityutilities::kInfo);
    inputs.printLog(TrivialClosureTestSummary(inputs, hRefRec, !isUnfolded), verbosityutilities::kInfo);
  }
}

void DoClosureTestV0Pt(InputSettings& inputs, int nIteration) {
  inputs.setVariableType(unfoldingutilities::VariableType::kK0SPt);
  DoClosureTestV0(inputs, nIteration);
}

void DoClosureTestV0Z(InputSettings& inputs, int nIteration) {
  inputs.setVariableType(unfoldingutilities::VariableType::kK0SZ);
  DoClosureTestV0(inputs, nIteration);
}

// ----------------------------------------------------------

void MakePlotsTrainingAndTest(InputSettings& inputs, TCanvas* canvas, int nPadsX, int nPadsY, vector<TH1D*> trainingHists, vector<TH1D*> testHists, vector<string> titles, bool logplot) {
  for (int iPad = 1; iPad <= nPadsX * nPadsY; iPad++) {
    if (iPad > trainingHists.size()) {
      inputs.printLog("MakePlotsTrainingAndTest() Warning: not enough histograms to fill all pads", verbosityutilities::kWarnings);
      break;
    }

    canvas->cd(iPad);
    gPad->SetLogy(logplot);
    TH1D* training = trainingHists[iPad - 1];
    TH1D* test     = testHists[iPad - 1];
    plotutils::setStyle(training, 0, -1, 1);
    plotutils::setStyle(test, 1, -1, 1);
    TH1F* frame = plotutils::DrawFrame(training);
    frame->SetMaximum(histutils::roundToNextPowerOfTen(frame->GetMaximum()));
    frame->SetMinimum(histutils::roundToPrevPowerOfTen(frame->GetMinimum()));
    if (frame->GetMinimum() < 1)
      frame->SetMinimum(0.1);

    if ((iPad % nPadsX) == 0) {
      // Kinematic efficiency should not have log plot
      gPad->SetLogy(false);
      frame->SetMinimum(0);
      frame->SetMaximum(1.1);
    }
    frame->SetStats(0);
    frame->SetTitle(titles[(iPad - 1) % nPadsX].c_str());
    frame->Draw();
    training->Draw("same");
    test->Draw("same");
  }
  if (inputs.saveFigs)
    canvas->SaveAs(canvas->GetName());
}

void CompareTrainingAndTestJets(InputSettings& inputs) {
  inputs.setVariableType(unfoldingutilities::VariableType::kJet);
  inputs.printLog("Comparing training and test distributions for jets.", verbosityutilities::kInfo);

  TFile* file = TFile::Open(inputs.inputFileName.c_str(), "READ");
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }

  string canvasName = "distComparison_" + unfoldingutilities::to_string(inputs.getVariableType()) + ".pdf";
  TCanvas* compCanvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1600, 600);

  vector<TH1D*> trainingHists;
  vector<TH1D*> testHists;
  vector<string> trainingNames = {
    rmutilities::training::nameGenJets,
    rmutilities::training::nameRecJets,
    rmutilities::training::nameFakeJets,
    rmutilities::training::nameMissJets,
    rmutilities::training::nameKinEffJets
  };
  vector<string> testNames = {
    rmutilities::testing::nameGenJets,
    rmutilities::testing::nameRecJets,
    rmutilities::testing::nameFakeJets,
    rmutilities::testing::nameMissJets,
    rmutilities::testing::nameKinEffJets
  };
  vector<string> titles = {
    "Gen", "Rec", "Fake", "Miss", "KinEff"
  };

  for (const auto& name : trainingNames) {
    TH1D* h = (TH1D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestJets() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    trainingHists.push_back(h);
  }
  for (const auto& name : testNames) {
    TH1D* h = (TH1D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestJets() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    testHists.push_back(h);
  }
  if (trainingHists.size() != testHists.size()) {
    inputs.printLog("CompareTrainingAndTestJets() Error: loaded different number of test and training histograms", verbosityutilities::kErrors);
    return;
  }

  int nPadsX = 5, nPadsY = 1;
  compCanvas->Divide(nPadsX, nPadsY);
  MakePlotsTrainingAndTest(inputs, compCanvas, nPadsX, nPadsY, trainingHists, testHists, titles, true);

  // Now take the ratio
  vector<TH1D*> trainingRatios = trainingHists;
  vector<TH1D*> testRatios = testHists;
  canvasName = "distRatio_" + unfoldingutilities::to_string(inputs.getVariableType()) + ".pdf";
  TCanvas* ratioCanvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1600, 600);
  ratioCanvas->Divide(nPadsX, nPadsY);
  for (int i = 0; i < trainingHists.size(); i++) {
    testRatios[i] = (TH1D*)histutils::divideWithProtection(testHists[i], trainingHists[i]);
    trainingRatios[i] = (TH1D*)histutils::divideWithProtection(trainingHists[i], trainingHists[i]);
  }
  MakePlotsTrainingAndTest(inputs, ratioCanvas, nPadsX, nPadsY, trainingRatios, testRatios, titles, false);
}

void CompareTrainingAndTestV0Pt(InputSettings& inputs) {
  inputs.setVariableType(unfoldingutilities::VariableType::kK0SPt);
  inputs.printLog("Comparing training and test distributions for K0S pt.", verbosityutilities::kInfo);

  TFile* file = TFile::Open(inputs.inputFileName.c_str(), "READ");
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }

  string canvasName = "distComparison_" + unfoldingutilities::to_string(inputs.getVariableType());
  canvasName = inputs.getNameFromPtJetProjection(canvasName, ".pdf");
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1600, 1200);

  vector<string> titles = { "Gen", "Rec", "Fake", "Miss", "KinEff" };
  vector<TH1D*> trainingHists;
  vector<TH1D*> testHists;

  vector<string> trainingNames = {
    rmutilities::training::nameGenK0SPt,
    rmutilities::training::nameRecK0SPt,
    rmutilities::training::nameFakeK0SPt,
    rmutilities::training::nameMissK0SPt,
    rmutilities::training::nameGenAndRecK0SPt
  };
  vector<string> testNames = {
    rmutilities::testing::nameGenK0SPt,
    rmutilities::testing::nameRecK0SPt,
    rmutilities::testing::nameFakeK0SPt,
    rmutilities::testing::nameMissK0SPt,
    rmutilities::testing::nameGenAndRecK0SPt
  };

  inputs.printLog("Loading training histograms", verbosityutilities::kDebug);
  for (const auto& name : trainingNames) {
    TH2D* h = (TH2D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestV0() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    string hName = h->GetName();
    array<int, 2> ptjetBins = histutils::getProjectionBins(h->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
    TH1D* hv0 = h->ProjectionX((hName + "_ptv0").c_str(), ptjetBins[0], ptjetBins[1]);
    trainingHists.push_back(hv0);
  }
  for (const auto& name : trainingNames) {
    TH2D* h = (TH2D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestV0() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    string hName = h->GetName();
    TH1D* hjet = h->ProjectionY((hName + "_ptjet").c_str(), 1, h->GetNbinsX());
    trainingHists.push_back(hjet);
  }

  inputs.printLog("Loading test histograms", verbosityutilities::kDebug);
  for (const auto& name : testNames) {
    TH2D* h = (TH2D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestV0() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    string hName = h->GetName();
    array<int, 2> ptjetBins = histutils::getProjectionBins(h->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
    TH1D* hv0 = h->ProjectionX((hName + "_ptv0").c_str(), ptjetBins[0], ptjetBins[1]);
    testHists.push_back(hv0);
  }
  for (const auto& name : testNames) {
    TH2D* h = (TH2D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestV0() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    string hName = h->GetName();
    TH1D* hjet = h->ProjectionY((hName + "_ptjet").c_str(), 1, h->GetNbinsX());
    testHists.push_back(hjet);
  }

  if (trainingHists.size() != testHists.size()) {
    inputs.printLog("CompareTrainingAndTestV0() Error: loaded different number of test and training histograms", verbosityutilities::kErrors);
    return;
  }

  inputs.printLog("Calculating the kinematic efficiency", verbosityutilities::kInfo);
  TH1D* trainingKinEff_ptv0 = (TH1D*)trainingHists[4]->Clone("trainingKinEff_ptv0");
  trainingKinEff_ptv0->Divide(trainingHists[0]);
  trainingHists[4] = trainingKinEff_ptv0;
  TH1D* trainingKinEff_ptjet = (TH1D*)trainingHists[9]->Clone("trainingKinEff_ptjet");
  trainingKinEff_ptjet->Divide(trainingHists[5]);
  trainingHists[9] = trainingKinEff_ptjet;

  TH1D* testKinEff_ptv0 = (TH1D*)testHists[4]->Clone("testKinEff_ptv0");
  testKinEff_ptv0->Divide(testHists[0]);
  testHists[4] = testKinEff_ptv0;
  TH1D* testKinEff_ptjet = (TH1D*)testHists[9]->Clone("testKinEff_ptjet");
  testKinEff_ptjet->Divide(testHists[5]);
  testHists[9] = testKinEff_ptjet;

  inputs.printLog("Drawing histograms", verbosityutilities::kDebug);
  int nPadsX = 5, nPadsY = 2;
  canvas->Divide(nPadsX, nPadsY);
  MakePlotsTrainingAndTest(inputs, canvas, nPadsX, nPadsY, trainingHists, testHists, titles, true);

  // Now take the ratio
  vector<TH1D*> trainingRatios = trainingHists;
  vector<TH1D*> testRatios = testHists;
  canvasName = "distRatio_" + unfoldingutilities::to_string(inputs.getVariableType());
  canvasName = inputs.getNameFromPtJetProjection(canvasName, ".pdf");
  TCanvas* ratioCanvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1600, 1200);
  ratioCanvas->Divide(nPadsX, nPadsY);
  for (int i = 0; i < trainingHists.size(); i++) {
    testRatios[i] = (TH1D*)histutils::divideWithProtection(testHists[i], trainingHists[i]);
    trainingRatios[i] = (TH1D*)histutils::divideWithProtection(trainingHists[i], trainingHists[i]);
  }
  MakePlotsTrainingAndTest(inputs, ratioCanvas, nPadsX, nPadsY, trainingRatios, testRatios, titles, false);
}

void CompareTrainingAndTestV0Z(InputSettings& inputs) {
  inputs.setVariableType(unfoldingutilities::VariableType::kK0SZ);
  inputs.printLog("Comparing training and test distributions for K0S z.", verbosityutilities::kInfo);

  TFile* file = TFile::Open(inputs.inputFileName.c_str(), "READ");
  if (!file) {
    inputs.printLog("Error: could not open file " + inputs.inputFileName, verbosityutilities::kErrors);
    return;
  }

  string canvasName = "distComparison_" + unfoldingutilities::to_string(inputs.getVariableType());
  canvasName = inputs.getNameFromPtJetProjection(canvasName, ".pdf");
  TCanvas* canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1600, 1200);

  vector<string> titles = { "Gen", "Rec", "Fake", "Miss", "KinEff" };
  vector<TH1D*> trainingHists;
  vector<TH1D*> testHists;

  vector<string> trainingNames = {
    rmutilities::training::nameGenK0SZ,
    rmutilities::training::nameRecK0SZ,
    rmutilities::training::nameFakeK0SZ,
    rmutilities::training::nameMissK0SZ,
    rmutilities::training::nameGenAndRecK0SZ
  };
  vector<string> testNames = {
    rmutilities::testing::nameGenK0SZ,
    rmutilities::testing::nameRecK0SZ,
    rmutilities::testing::nameFakeK0SZ,
    rmutilities::testing::nameMissK0SZ,
    rmutilities::testing::nameGenAndRecK0SZ
  };

  inputs.printLog("Loading training histograms", verbosityutilities::kDebug);
  for (const auto& name : trainingNames) {
    TH2D* h = (TH2D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestV0() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    string hName = h->GetName();
    array<int, 2> ptjetBins = histutils::getProjectionBins(h->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
    TH1D* hv0 = h->ProjectionX((hName + "_zv0").c_str(), ptjetBins[0], ptjetBins[1]);
    trainingHists.push_back(hv0);
  }
  for (const auto& name : trainingNames) {
    TH2D* h = (TH2D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestV0() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    string hName = h->GetName();
    TH1D* hjet = h->ProjectionY((hName + "_ptjet").c_str(), 1, h->GetNbinsX());
    trainingHists.push_back(hjet);
  }

  inputs.printLog("Loading test histograms", verbosityutilities::kDebug);
  for (const auto& name : testNames) {
    TH2D* h = (TH2D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestV0() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    string hName = h->GetName();
    array<int, 2> ptjetBins = histutils::getProjectionBins(h->GetYaxis(), inputs.ptjetminProjection, inputs.ptjetmaxProjection);
    TH1D* hv0 = h->ProjectionX((hName + "_zv0").c_str(), ptjetBins[0], ptjetBins[1]);
    testHists.push_back(hv0);
  }
  for (const auto& name : testNames) {
    TH2D* h = (TH2D*)file->Get(name.c_str());
    if (!h) {
      inputs.printLog("CompareTrainingAndTestV0() Warning: could not find " + name + " in file " + inputs.inputFileName, verbosityutilities::kWarnings);
      continue;
    }
    string hName = h->GetName();
    TH1D* hjet = h->ProjectionY((hName + "_ptjet").c_str(), 1, h->GetNbinsX());
    testHists.push_back(hjet);
  }

  if (trainingHists.size() != testHists.size()) {
    inputs.printLog("CompareTrainingAndTestV0() Error: loaded different number of test and training histograms", verbosityutilities::kErrors);
    return;
  }

  inputs.printLog("Calculating the kinematic efficiency", verbosityutilities::kInfo);
  TH1D* trainingKinEff_zv0 = (TH1D*)trainingHists[4]->Clone("trainingKinEff_zv0");
  trainingKinEff_zv0->Divide(trainingHists[0]);
  trainingHists[4] = trainingKinEff_zv0;
  TH1D* trainingKinEff_ptjet = (TH1D*)trainingHists[9]->Clone("trainingKinEff_ptjet");
  trainingKinEff_ptjet->Divide(trainingHists[5]);
  trainingHists[9] = trainingKinEff_ptjet;

  TH1D* testKinEff_zv0 = (TH1D*)testHists[4]->Clone("testKinEff_zv0");
  testKinEff_zv0->Divide(testHists[0]);
  testHists[4] = testKinEff_zv0;
  TH1D* testKinEff_ptjet = (TH1D*)testHists[9]->Clone("testKinEff_ptjet");
  testKinEff_ptjet->Divide(testHists[5]);
  testHists[9] = testKinEff_ptjet;

  inputs.printLog("Drawing histograms", verbosityutilities::kDebug);
  int nPadsX = 5, nPadsY = 2;
  canvas->Divide(nPadsX, nPadsY);
  MakePlotsTrainingAndTest(inputs, canvas, nPadsX, nPadsY, trainingHists, testHists, titles, true);

  // Now take the ratio
  vector<TH1D*> trainingRatios = trainingHists;
  vector<TH1D*> testRatios = testHists;
  canvasName = "distRatio_" + unfoldingutilities::to_string(inputs.getVariableType());
  canvasName = inputs.getNameFromPtJetProjection(canvasName, ".pdf");
  TCanvas* ratioCanvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1600, 1200);
  ratioCanvas->Divide(nPadsX, nPadsY);
  for (int i = 0; i < trainingHists.size(); i++) {
    testRatios[i] = (TH1D*)histutils::divideWithProtection(testHists[i], trainingHists[i]);
    trainingRatios[i] = (TH1D*)histutils::divideWithProtection(trainingHists[i], trainingHists[i]);
  }
  MakePlotsTrainingAndTest(inputs, ratioCanvas, nPadsX, nPadsY, trainingRatios, testRatios, titles, false);
}

// ----------------------------------------------------------
//
// Interface
//
// ----------------------------------------------------------

void createresponses(InputSettings& x) {
  CreateResponseJets(x);
  CreateResponseV0Pt(x);
  CreateResponseV0Z(x);
}

void unfolding(InputSettings& x) {
  for (int iIteration = x.minIteration; iIteration <= x.maxIteration; iIteration++) {
    DoUnfoldingJets(x, iIteration);
    DoUnfoldingV0Pt(x, iIteration);
    DoUnfoldingV0Z(x, iIteration);
  }
}

void checkclosure(InputSettings& x) {
  gROOT->SetBatch(kTRUE);
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

void checktrainingtest(InputSettings& x) {
  gROOT->SetBatch(kTRUE);
  CompareTrainingAndTestJets(x);
  x.setPtJetProjection(10., 20.);
  CompareTrainingAndTestV0Pt(x);
  CompareTrainingAndTestV0Z(x);
  x.setPtJetProjection(20., 30.);
  CompareTrainingAndTestV0Pt(x);
  CompareTrainingAndTestV0Z(x);
}

enum ActionType { kCreateResponses, kDoUnfolding, kCheckClosure, kCompareDists, kDoAll };
void doclosuretest(ActionType action, InputSettings& inputs, array<string, 4> fileNames) {
  string trainingFileName = fileNames[0];
  string testFileName     = fileNames[1];
  string responseFileName = fileNames[2];
  string closureFileName  = fileNames[3];

  if (inputs.minIteration < 0 || inputs.maxIteration < 0) {
    if (action == kDoUnfolding || action == kCheckClosure || action == kDoAll) {
      inputs.printLog("doclosuretest(): minIteration and maxIteration must be set to do unfolding or check closure", verbosityutilities::kErrors);
      return;
    }
  }

  if (action == kCreateResponses || action == kDoAll) {
    inputs.inputFileName  = trainingFileName;
    inputs.outputFileName = responseFileName;
    createresponses(inputs);
  }
  if (action == kDoUnfolding || action == kDoAll) {
    inputs.inputFileName    = testFileName;
    inputs.responseFileName = responseFileName;
    inputs.outputFileName   = closureFileName;
    unfolding(inputs);
  }
  if (action == kCheckClosure || action == kDoAll) {
    inputs.inputFileName = closureFileName;
    checkclosure(inputs);
  }
  if (action == kCompareDists || action == kDoAll) {
    inputs.inputFileName = closureFileName;
    checktrainingtest(inputs);
  }
}

void closuretest520161(ActionType action, bool doTrivialClosureTest, int minIteration = -1, int maxIteration = -1) {
  int train = 520161;
  string trainingFileName, testFileName, responseFileName, closureFileName;
  if (doTrivialClosureTest) {
    trainingFileName = to_string(train) + ".root";
    testFileName     = trainingFileName;
  } else {
    trainingFileName = to_string(train) + "_40.root";
    testFileName     = to_string(train) + "_60.root";
  }
  responseFileName = "RooUnfoldResponse_" + trainingFileName;
  closureFileName  = "ClosureTest_" + trainingFileName;

  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);
  x.setPtV0Rec(1., 40);
  x.autoTemplateHists();
  x.doTrivialClosureTest = doTrivialClosureTest;
  x.minIteration = minIteration;
  x.maxIteration = maxIteration;
  doclosuretest(action, x, array<string, 4>{trainingFileName, testFileName, responseFileName, closureFileName});
}

# endif
