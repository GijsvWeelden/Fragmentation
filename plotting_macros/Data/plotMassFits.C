#include <vector>
#include <iostream>
#include <typeinfo>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TString.h"
#include "TLegend.h"

#include "../histUtils.C"
#include "../plotUtils.C"

#ifndef __PLOTMASSFITS_H__
#define __PLOTMASSFITS_H__

double chisqInRange(TH1* hist, TF1* fit, int minBin, int maxBin) {
  double chisq = 0;
  for (int i = minBin; i <= maxBin; i++) {
    double bc = hist->GetBinContent(i);
    double be = hist->GetBinError(i);
    double fx = fit->Eval(hist->GetBinCenter(i));

    if (bc < 1e-50) // Protect against zero division
      continue;

    chisq += (fx - bc) * (fx - bc) / (be * be);
  }
  return chisq;
}

double chisqInRange(TH1* hist, TF1* fit, double xmin, double xmax) {
  array<int, 2> bins = histutils::getProjectionBins(hist->GetXaxis(), xmin, xmax);
  return chisqInRange(hist, fit, bins[0], bins[1]);
}

double ndfInRange(TF1* fit, int minBin, int maxBin) {
  return (1 + maxBin - minBin - fit->GetNpar());
}

double ndfInRange(TH1* data, TF1* fit, double xmin, double xmax) {
  array<int, 2> bins = histutils::getProjectionBins(data->GetXaxis(), xmin, xmax);
  return ndfInRange(fit, bins[0], bins[1]);
}

// -------------------------------------------------------------------------------------------------
//
// Signal Fractions
//
// -------------------------------------------------------------------------------------------------

array<double, 2> GaussGaussSigRegionFromSteps(double n, TF1* f) {
  double ampNarrow = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigmaNarrow = f->GetParameter(2);
  double ampWide = f->GetParameter(3);
  double sigmaWide = f->GetParameter(4);

  double Sigma = (ampNarrow * sigmaNarrow + ampWide * sigmaWide) / (ampNarrow + ampWide);
  return {mu - n * Sigma, mu + n * Sigma};
}

array<double, 2> GaussGaussSigRegionFromFrac(double s, TF1* f) {
  double ampNarrow = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigmaNarrow = f->GetParameter(2);
  double ampWide = f->GetParameter(3);
  double sigmaWide = f->GetParameter(4);

  double Sigma = (ampNarrow * sigmaNarrow + ampWide * sigmaWide) / (ampNarrow + ampWide);
  double n = s * (ampNarrow * sigmaNarrow + ampWide * sigmaWide) / (ampNarrow + ampWide);

  return {mu - n * Sigma, mu + n * Sigma};
}

double GaussGaussSigFrac(double n, TF1* f) {
  double ampNarrow = f->GetParameter(0);
  double sigmaNarrow = f->GetParameter(2);
  double ampWide = f->GetParameter(3);
  double sigmaWide = f->GetParameter(4);

  double Sigma = (ampNarrow * sigmaNarrow + ampWide * sigmaWide) / (ampNarrow + ampWide);

  double x = n * Sigma / (sqrt(2) * sigmaNarrow);
  double y = n * Sigma / (sqrt(2) * sigmaWide);

  x = TMath::Erf(x);
  x *= ampNarrow * sigmaNarrow;

  y = TMath::Erf(y);
  y *= ampWide * sigmaWide;

  double s = (x + y) / (ampNarrow * sigmaNarrow + ampWide * sigmaWide);
  return s;
}

array<double, 2> GaussGaussExpSigRegion(double n, TF1* f) {
  double ampNarrow = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigmaNarrow = f->GetParameter(2);
  double lambda = f->GetParameter(3); // Crossover point
  double ampWide = f->GetParameter(4);
  double sigmaWide = f->GetParameter(5);

  double tau = sigmaNarrow / lambda;

  // Make double Gaussian from f
  TF1* g = new TF1("g", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4])", f->GetXmin(), f->GetXmax());
  g->SetParameters(ampNarrow, mu, sigmaNarrow, ampWide, sigmaWide);
  array<double, 2> gaussRegion = GaussGaussSigRegionFromSteps(n, g);

  double leftSide = gaussRegion[0];
  double rightSide = gaussRegion[1];
  if (n < lambda) rightSide = mu + n * sigmaNarrow;
  if (n > lambda) rightSide = mu + lambda * sigmaNarrow + (n - lambda) * tau;

  // TODO: Think about how to define the weird signal region here
  return {leftSide, rightSide};
}

double GaussGaussExpSigFrac(double n, TF1* f, bool leftSide) {
  double ampNarrow = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigmaNarrow = f->GetParameter(2);
  double lambda = f->GetParameter(3); // Crossover point
  double ampWide = f->GetParameter(4);
  double sigmaWide = f->GetParameter(5);

  if (leftSide) {
    // Make double Gaussian from f
    TF1* g = new TF1("g", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4])", f->GetXmin(), f->GetXmax());
    g->SetParameters(ampNarrow, mu, sigmaNarrow, ampWide, sigmaWide);
    return GaussGaussSigFrac(n, f);
  }

  double x = (lambda + n / lambda - 1) * sigmaNarrow / (sqrt(2) * sigmaWide);
  x = TMath::Erf(x);
  x *= ampWide * sigmaWide;

  double y = lambda / sqrt(2);
  y = TMath::Erf(y);
  y *= ampNarrow * sigmaNarrow;
  y *= sqrt(TMath::PiOver2());

  double z = exp(- lambda*lambda / 2);
  z *= ampNarrow * sigmaNarrow / lambda;

  double s = (x + y + z * (1 - exp(lambda - n))) / (ampNarrow * sigmaNarrow + y + z);
  return s;
}

array<double, 2> ExpGaussExpSigRegion(double n, TF1* f) {
  double amp = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  double lambdaL = f->GetParameter(3); // Crossover point left
  double lambdaR = f->GetParameter(4); // Crossover point right

  double tauL = sigma / lambdaL;
  double tauR = sigma / lambdaR;

  double leftSide = mu - n * sigma;
  if (n > lambdaL) leftSide = mu - lambdaL * sigma - (n - lambdaL) * tauL;

  double rightSide = mu + n * sigma;
  if (n > lambdaR) rightSide = mu + lambdaR * sigma + (n - lambdaR) * tauR;

  return {leftSide, rightSide};
}

double ExpGaussExpSigFrac(double n, TF1* f, bool leftSide) {
  double amp = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  double lambdaL = f->GetParameter(3); // Crossover point left
  double lambdaR = f->GetParameter(4); // Crossover point right

  double lambda = (leftSide) ? lambdaL : lambdaR;

  double x = lambda / sqrt(2);
  x = TMath::Erf(x);
  x *= lambda * sqrt(TMath::PiOver2());

  double y = exp(- lambda*lambda / 2);

  double s = (x + y * (1 - exp(lambda - n))) / (x + y);
  return s;
}

// -------------------------------------------------------------------------------------------------
//
// Struct for input settings
//
// -------------------------------------------------------------------------------------------------

struct InputSettings{
  private:
  public:
    int train = 0;
    int rebinNumber = -1;
    int fitType = -1;
    int verbosity = 1; // kWarnings
    string hadron = "";
    string fitName = "";
    string histName = "";
    string inputFileName = "";
    string outputFileName = "";
    double ptmin = -1e3, ptmax = -1e3, lowpt = -1e3, highpt = -1e3;
    double fitmin = -1e3, fitmax = -1e3; // Fit range
    double textSize = 0.04;
    bool logplot = false;
    bool normaliseData = false;
    bool fixMu = false;
    bool printDataSet = false;
    bool drawLegend = true;
    bool drawLatex = true;
    vector<vector<double>> ptBinEdges = {};

    double massWindowMin = -1., massWindowMax = -1.;
    double polInitx0 = -1., polInitx1 = -1., polInitx2 = -1.;
    double signalRegionMin = -1., signalRegionMax = -1.;
    double nSigma = -1., nSigmaL = -1., nSigmaR = -1.;

    enum FitType {kPol1BreitWigner, kPol2BreitWigner, kPol1BreitWignerXex, kPol2BreitWignerXex, kPol1ExpGausExp, kPol2ExpGausExp, kPol1GausExp, kPol2GausExp, kPol1GausGaus, kPol2GausGaus, kPol1GausGausExp, kPol2GausGausExp, kPol1GausGausXex, kPol2GausGausXex, kPol1GausXex, kPol2GausXex, kPol1Voigt, kPol2Voigt};

    enum Verbosity {kErrors, kWarnings, kInfo, kDebug, kDebugMax};

    double getMass();
    string getFitExpression();
    string getSaveNameFromPt(string prefix, string suffix);
    int inputIssue(string home, string obj);
    string printLog(string message, int verbThreshold);
    string setFitName();
    string setFitName(int fit);
    int setFitType();
    int setFitType(string fit);
    void setFitX(double a, double b);
    void setFitXFromHadron();
    string setInputFileNameFromTrain();
    string setInputFileNameFromFit();
    void setMassWindow(double a, double b);
    void setMassWindowDiff(double a, double b);
    void setMassWindowDiffFromHadron();
    void setPt(double a, double b);
    vector<vector<double>> setPtBinEdgesFromHadron();
    vector<vector<double>> setPtBinEdgesSorted(vector<vector<double>> x);
    void setPolInitX(double x0, double x1, double x2);
    void setPolInitXFromHadron();
    void setSignalRegion(double a, double b);
    template <typename T> int writeOutputToFile(T* obj);
};

double InputSettings::getMass() {
  if (this->hadron == "K0S")
    return histutils::MassK0S;
  if (this->hadron == "Lambda" || this->hadron == "AntiLambda")
    return histutils::MassLambda;

  return -1.;
}

string InputSettings::getFitExpression() {
  string s = "";
  switch (this->fitType) {
    case kPol1BreitWigner:
      s = "breitwigner(x,[0],[1],[2]) + [3]+[4]*x";
      break;
    case kPol2BreitWigner:
      s = "breitwigner(x,[0],[1],[2]) + [3]+[4]*x+[5]*x*x";
      break;
    case kPol1BreitWignerXex:
      s = "breitwigner(x,[0],[1],[2]) + max(0., [3] * (x-[4]) * TMath::Exp(-(x-[4])/[5])) + [6]+[7]*x";
      break;
    case kPol2BreitWignerXex:
      s = "breitwigner(x,[0],[1],[2]) + max(0., [3] * (x-[4]) * TMath::Exp(-(x-[4])/[5])) + [6]+[7]*x+[8]*x*x";
      break;
    case kPol1ExpGausExp:
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x > ([1]-[2]*[3])) * (x < ([1]+[2]*[4])) + [0]*TMath::Exp((x - [1] + [2]*[3]/2.)/([2]/[3])) * (x <= [1]-[2]*[3]) + [0]*TMath::Exp(-(x - [1] - [2]*[4]/2.)/([2]/[4])) * (x >= [1]+[2]*[4]) + [5]+[6]*x";
      break;
    case kPol2ExpGausExp:
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x > ([1]-[2]*[3])) * (x < ([1]+[2]*[4])) + [0]*TMath::Exp((x - [1] + [2]*[3]/2.)/(-[2]/[3])) * (x <= [1]-[2]*[3]) + [0]*TMath::Exp(-(x - [1] - [2]*[4]/2.)/([2]/[4])) * (x >= [1]+[2]*[4]) + [5]+[6]*x+[7]*x*x";
      break;
    case kPol1GausExp:
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - [1] - [2]*[3]/2.)/([2]/[3])) * (x >= [1]+[2]*[3]) + [4]+[5]*x";
      break;
      case kPol2GausExp:
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - [1] - [2]*[3]/2.)/([2]/[3])) * (x >= [1]+[2]*[3]) + [4]+[5]*x+[6]*x*x";
      break;
    case kPol1GausGaus:
      if (this->lowpt > 20. - 1e-5 && this->lowpt < 25. - 1e-5 && this->highpt > 20. - 1e-5 && this->highpt < 25. + 1e-5) {
        s = "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4]) + max(0., [5]+[6]*x)";
      } else {
        s = "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4]) + [5]+[6]*x";
      }
      break;
    case kPol2GausGaus:
      s = "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4]) + [5]+[6]*x+[7]*x*x";
      break;
    case kPol1GausGausExp:
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - [1] - [2]*[3]/2.)/([2]/[3])) * (x >= [1]+[2]*[3]) + [4]*TMath::Gaus(x,[1],[5]) + [6]+[7]*x";
      break;
      case kPol2GausGausExp:
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - [1] - [2]*[3]/2.)/([2]/[3])) * (x >= [1]+[2]*[3]) + [4]*TMath::Gaus(x,[1],[5]) + [6]+[7]*x+[8]*x*x";
      break;
    case kPol1GausGausXex:
      s = "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4]) + max(0., [5] * (x-[6]) * TMath::Exp(-(x-[6])/[7])) + [8]+[9]*x";
      break;
    case kPol2GausGausXex:
      s = "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4]) + max(0., [5] * (x-[6]) * TMath::Exp(-(x-[6])/[7])) + [8]+[9]*x+[10]*x*x";
      break;
    case kPol1GausXex:
      s = "[0]*TMath::Gaus(x,[1],[2]) + max(0., [3] * (x-[4]) * TMath::Exp(-(x-[4])/[5])) + [6]+[7]*x";
      break;
    case kPol2GausXex:
      s = "[0]*TMath::Gaus(x,[1],[2]) + max(0., [3] * (x-[4]) * TMath::Exp(-(x-[4])/[5])) + [6]+[7]*x+[8]*x*x";
      break;
    case kPol1Voigt:
      s = "[0]*TMath::Voigt(x-[1],[2],[3]) + [4]+[5]*x";
      break;
    case kPol2Voigt:
      s = "[0]*TMath::Voigt(x-[1],[2],[3]) + [4]+[5]*x+[6]*x*x";
      break;
    default:
      string s = "InputSettings::getFitExpression() Error: requested unknown function " + this->fitName;
      printLog(s, kErrors);
      break;
  }

  return s;
}

string InputSettings::getSaveNameFromPt(string prefix, string suffix = "") {
  string s = TString::Format("%s_pt%.1f-%.1f%s", prefix.c_str(), this->lowpt, this->highpt, suffix.c_str()).Data();
  return s;
}

int InputSettings::inputIssue(string home, string obj) {
  if (obj == "hadron") {
    if (this->hadron == "K0S" || this->hadron == "Lambda" || this->hadron == "AntiLambda") return 0;
  } else {
    string s = "InputSettings::inputIssue() Error: requested unknown object " + obj;
    printLog(s, kErrors);
    return 1;
  }

  string s = "InputSettings::" + home + "() Error: " + obj + " not initialised";
  printLog(s, kErrors);
  return 2;
}

string InputSettings::printLog(string message, int verbThreshold) {
  if (this->verbosity < verbThreshold)
    return "";

  cout << message << endl;
  return message;
}

string InputSettings::setFitName() {
  switch (this->fitType) {
    case kPol1BreitWigner:
      this->fitName = "pol1BreitWigner";
      break;
    case kPol2BreitWigner:
      this->fitName = "pol2BreitWigner";
      break;
    case kPol1BreitWignerXex:
      this->fitName = "pol1BreitWignerXex";
      break;
    case kPol2BreitWignerXex:
      this->fitName = "pol2BreitWignerXex";
      break;
    case kPol1ExpGausExp:
      this->fitName = "pol1ExpGausExp";
      break;
    case kPol2ExpGausExp:
      this->fitName = "pol2ExpGausExp";
      break;
    case kPol1GausExp:
      this->fitName = "pol1GausExp";
      break;
    case kPol2GausExp:
      this->fitName = "pol2GausExp";
      break;
    case kPol1GausGaus:
      this->fitName = "pol1GausGaus";
      break;
    case kPol2GausGaus:
      this->fitName = "pol2GausGaus";
      break;
    case kPol1GausGausExp:
      this->fitName = "pol1GausGausExp";
      break;
    case kPol2GausGausExp:
      this->fitName = "pol2GausGausExp";
      break;
    case kPol1GausGausXex:
      this->fitName = "pol1GausGausXex";
      break;
    case kPol2GausGausXex:
      this->fitName = "pol2GausGausXex";
      break;
    case kPol1GausXex:
      this->fitName = "pol1GausXex";
      break;
    case kPol2GausXex:
      this->fitName = "pol2GausXex";
      break;
    case kPol1Voigt:
      this->fitName = "pol1Voigt";
      break;
    case kPol2Voigt:
      this->fitName = "pol2Voigt";
      break;
    default:
      string s = "InputSettings::setFitName() Error: requested unknown function";
      printLog(s, kErrors);
  }
  return this->fitName;
}

string InputSettings::setFitName(int fit) {
  this->fitType = fit;
  return setFitName();
}

int InputSettings::setFitType() {
  if (this->fitName == "pol1BreitWigner") {
    this->fitType = kPol1BreitWigner;
  } else if (this->fitName == "pol2BreitWigner") {
    this->fitType = kPol2BreitWigner;
  } else if (this->fitName == "pol1BreitWignerXex") {
    this->fitType = kPol1BreitWignerXex;
  } else if (this->fitName == "pol2BreitWignerXex") {
    this->fitType = kPol2BreitWignerXex;
  } else if (this->fitName == "pol1ExpGausExp") {
    this->fitType = kPol1ExpGausExp;
  } else if (this->fitName == "pol2ExpGausExp") {
    this->fitType = kPol2ExpGausExp;
  } else if (this->fitName == "pol1GausExp") {
    this->fitType = kPol1GausExp;
  } else if (this->fitName == "pol2GausExp") {
    this->fitType = kPol2GausExp;
  } else if (this->fitName == "pol1GausGaus") {
    this->fitType = kPol1GausGaus;
  } else if (this->fitName == "pol2GausGaus") {
    this->fitType = kPol2GausGaus;
  } else if (this->fitName == "pol1GausGausExp") {
    this->fitType = kPol1GausGausExp;
  } else if (this->fitName == "pol2GausGausExp") {
    this->fitType = kPol2GausGausExp;
  } else if (this->fitName == "pol1GausGausXex") {
    this->fitType = kPol1GausGausXex;
  } else if (this->fitName == "pol2GausGausXex") {
    this->fitType = kPol2GausGausXex;
  } else if (this->fitName == "pol1GausXex") {
    this->fitType = kPol1GausXex;
  } else if (this->fitName == "pol2GausXex") {
    this->fitType = kPol2GausXex;
  } else if (this->fitName == "pol1Voigt") {
    this->fitType = kPol1Voigt;
  } else if (this->fitName == "pol2Voigt") {
    this->fitType = kPol2Voigt;
  } else {
    string s = "InputSettings::getFitType() Error: requested unknown function " + this->fitName;
    printLog(s, kErrors);
  }

  return this->fitType;
}

int InputSettings::setFitType(string fit) {
  this->fitName = fit;
  return setFitType();
}

string InputSettings::setInputFileNameFromTrain() {
  string s = "~/cernbox/TrainOutput/" + to_string(this->train) + "/AnalysisResults.root";
  this->inputFileName = s;
  return s;
}

string InputSettings::setInputFileNameFromFit() {
  string s = "./" + this->hadron + "_" + this->fitName + ".root";
  this->inputFileName = s;
  return s;
}

void InputSettings::setFitX(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setFitX() Error: fitmin > fitmax";
    printLog(s, kErrors);
    return;
  }
  this->fitmin = a;
  this->fitmax = b;
}

void InputSettings::setFitXFromHadron() {
  if (this->hadron == "K0S") {
    this->setFitX(0.45, 0.55);
  } else
    this->setFitX(-1., -1.);
}

void InputSettings::setMassWindow(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setMassWindow() Error: massWindowMin > massWindowMax";
    printLog(s, kErrors);
    return;
  }
  this->massWindowMin = a;
  this->massWindowMax = b;
}

void InputSettings::setMassWindowDiff(double a, double b) {
  double mass = this->getMass();
  this->setMassWindow(mass - a, mass + b);
}

void InputSettings::setMassWindowDiffFromHadron() {
  if (this->hadron == "K0S") {
    this->setMassWindowDiff(2e-2, 2e-2);
  } else {
    this->setMassWindowDiff(5e-3, 5e-3);
  }
}

void InputSettings::setPolInitX(double x0, double x1 = -1, double x2 = -2) {
  // vector<double> x = {x0, x1, x2};
  // std::sort(x.begin(), x.end());
  // x.erase( remove_if( x.begin(), x.end(), []( int i ){ return i < 0; } ), x.end() ); // Remove negative values

  // if (x.size() > 0) this->polInitx0 = x[0];
  // if (x.size() > 1) this->polInitx1 = x[1];
  // if (x.size() > 2) this->polInitx2 = x[2];

  this->polInitx0 = x0;
  this->polInitx1 = x1;
  this->polInitx2 = x2;
}

void InputSettings::setPolInitXFromHadron() {
  if (this->hadron == "K0S") {
    this->setPolInitX(0.45, 0.55, 0.57);
  } else {
    this->setPolInitX(1.09, 1.1, 1.13);
  }
}

void InputSettings::setPt(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setPt() Error: ptmin > ptmax";
    printLog(s, kErrors);
    return;
  }
  this->ptmin = a;
  this->ptmax = b;
  this->lowpt = a;
  this->highpt = b;
}

vector<vector<double>> InputSettings::setPtBinEdgesFromHadron() {
  vector<vector<double>> x;
  if (inputIssue("setPtBinEdgesFromHadron", "hadron"))
    return x;

  x.push_back({0., 1.});
  x.push_back({1., 2.});
  x.push_back({2., 3.});
  x.push_back({3., 4.});
  x.push_back({4., 5.});
  x.push_back({5., 10.});
  x.push_back({10., 15.});
  x.push_back({15., 20.});
  if (this->hadron == "K0S") {
    x.push_back({20., 25.});
    x.push_back({25., 30.});
    x.push_back({30., 40.});
  } else {
    x.push_back({20., 30.});
    x.push_back({30., 40.});
  }
  x = this->setPtBinEdgesSorted(x);
  return x;
}

vector<vector<double>> InputSettings::setPtBinEdgesSorted(vector<vector<double>> x) {
  sort(x.begin(), x.end());
  this->ptBinEdges = x;
  return x;
}

void InputSettings::setSignalRegion(double a, double b) {
  if (a > b) {
    string s = "InputSettings::setSignalRegion() Error: signalRegionMin > signalRegionMax";
    printLog(s, kErrors);
    return;
  }
  this->signalRegionMin = a;
  this->signalRegionMax = b;
}

template <typename T>
int InputSettings::writeOutputToFile(T* obj) {
  if (!obj)
    return 1;

  TFile* file = TFile::Open(this->outputFileName.c_str(), "UPDATE");
  obj->Write(obj->GetName(), TObject::kOverwrite);
  file->Close();
  return 0;
}

// -------------------------------------------------------------------------------------------------
//
// Struct for fitting mass spectra
//
// -------------------------------------------------------------------------------------------------

struct MassFitter {
  private:
  public:
    InputSettings* inputs;
    TH1* data = nullptr;
    TF1* fit = nullptr;
    vector<TF1*> fitParts = {};
    TH1* fitResults = nullptr;
    TH1* fitParams = nullptr;
    TH1* pull = nullptr;
    TH1* residual = nullptr;

    vector<double> signal = {};
    vector<double> background = {};
    vector<double> signalPlusBackground = {};

    MassFitter() { inputs = new InputSettings(); }
    MassFitter(InputSettings& x) { inputs = &x; }

    // double calcSignalRegionSize();
    void calcSigBkg();
    void doFitting(); // Full execution of workflow
    void fixFitInPost();
    TF1* loadFitFunction();
    TF1* loadSavedFitFunction(); // Load fit function from file
    TH1* loadSavedMassHist(); // Load 1D histogram from file
    TH1* loadMassHist(); // Load 1D hist from data THn
    void setFitInitialValues();
    vector<TF1*> loadSavedFitParts();
    vector<TF1*> loadFitParts();
    TH1* loadFitParams();
    TH1* loadFitResults();
    TH1* loadResidualHist();
    TH1* loadPullHist();
    string setHistNameFromTrain();
    int setRebinNumberFromHadronAndPt();
    void setSignalRegionFromSigma();
    void writeOutputsToFile();
};

void MassFitter::calcSigBkg() {
  if (!this->data || !this->fit || this->fitParts.size() == 0)
    return;

  if (this->signal.size() > 0 || this->background.size() > 0 || this->signalPlusBackground.size() > 0) {
    string s = "MassFitter::calcSigBkg() Warning: signal, background and/or signalPlusBackground already set, clearing the vectors";
    this->inputs->printLog(s, InputSettings::kWarnings);
    this->signal = {};
    this->background = {};
    this->signalPlusBackground = {};
  }

  // The number of the function is 1 + its index in fitParts
  // I.e. fitPart _f1 is at index 0
  vector<int> bkgFitNumbers = {};
  if (this->inputs->fitType == InputSettings::kPol1BreitWigner || this->inputs->fitType == InputSettings::kPol2BreitWigner) {
    bkgFitNumbers.push_back(2);
  } else if (this->inputs->fitType == InputSettings::kPol1BreitWignerXex || this->inputs->fitType == InputSettings::kPol2BreitWignerXex) {
    bkgFitNumbers.push_back(3);
  } else if (this->inputs->fitType == InputSettings::kPol1ExpGausExp || this->inputs->fitType == InputSettings::kPol2ExpGausExp) {
    bkgFitNumbers.push_back(4);
  } else if (this->inputs->fitType == InputSettings::kPol1GausExp || this->inputs->fitType == InputSettings::kPol2GausExp){
    bkgFitNumbers.push_back(3);
  } else if (this->inputs->fitType == InputSettings::kPol1GausGaus || this->inputs->fitType == InputSettings::kPol2GausGaus){
    bkgFitNumbers.push_back(3);
  } else if (this->inputs->fitType == InputSettings::kPol1GausGausExp || this->inputs->fitType == InputSettings::kPol2GausGausExp){
    bkgFitNumbers.push_back(4);
  } else if (this->inputs->fitType == InputSettings::kPol1GausGausXex || this->inputs->fitType == InputSettings::kPol2GausGausXex){
    bkgFitNumbers.push_back(4);
  } else if (this->inputs->fitType == InputSettings::kPol1GausXex || this->inputs->fitType == InputSettings::kPol2GausXex){
    bkgFitNumbers.push_back(3);
  } else if (this->inputs->fitType == InputSettings::kPol1Voigt || this->inputs->fitType == InputSettings::kPol2Voigt){
    bkgFitNumbers.push_back(2);
  } else {
    string s = "MassFitter::calcSigBkg() Error: do not know how to set parameters for this function";
    this->inputs->printLog(s, InputSettings::kErrors);
    return;
  }

  if (this->inputs->signalRegionMin < 0 || this->inputs->signalRegionMax < 0) {
    string s = "MassFitter::calcSigBkg() Automatically setting signal region to mu +/- n sigma";
    this->inputs->printLog(s, InputSettings::kInfo);
    this->setSignalRegionFromSigma();
  }

  // Load the signal and background fit functions
  vector<TF1*> bkgFits = {}; vector<TF1*> sigFits = {};
  for (int iBkg = 0; iBkg < this->fitParts.size(); iBkg++) {
    if(std::find(bkgFitNumbers.begin(), bkgFitNumbers.end(), iBkg+1) != bkgFitNumbers.end()) // iBkg in bkgFitNumbers
      bkgFits.push_back(this->fitParts[iBkg]);
    else // iBkg not in bkgFitNumbers
      sigFits.push_back(this->fitParts[iBkg]);
  }

  // Calculate signal and background with data hist and bkg fit function
  array<int, 2> sigBins = histutils::getProjectionBins(this->data->GetXaxis(), this->inputs->signalRegionMin, this->inputs->signalRegionMax);
  double xmin = this->data->GetXaxis()->GetBinLowEdge(sigBins[0]);
  double xmax = this->data->GetXaxis()->GetBinUpEdge(sigBins[1]);

  double signalPlusBackground = this->data->Integral(sigBins[0], sigBins[1], "width");
  double background = 0;
  for (auto bkgFit : bkgFits) {
    bkgFit->SetRange(xmin, xmax);
    background += bkgFit->Integral(xmin, xmax);
  }
  double signal = signalPlusBackground - background;
  this->signalPlusBackground.push_back(signalPlusBackground);
  this->background.push_back(background);
  this->signal.push_back(signal);

  // Calculate signal with fit function
  signalPlusBackground = 0, signal = 0;
  for (auto sigFit : sigFits) {
    sigFit->SetRange(xmin, xmax);
    signalPlusBackground += sigFit->Integral(xmin, xmax);
  }
  signal = signalPlusBackground - background;
  this->signalPlusBackground.push_back(signalPlusBackground);
  this->background.push_back(background);
  this->signal.push_back(signal);
}

void MassFitter::doFitting() {
  this->loadMassHist();
  this->loadFitFunction();
  this->setFitInitialValues();
  this->data->Fit(this->fit, "R");
  this->fixFitInPost(); // Applies any post-fit fixes, like swapping gaussians
  this->loadFitParts();
  this->loadFitParams();
  this->loadFitResults();
  this->loadResidualHist();
  this->loadPullHist();
  this->writeOutputsToFile();
}

void MassFitter::fixFitInPost() {
  if (this->inputs->hadron == "K0S" && this->inputs->fitType == InputSettings::kPol1GausGaus) {
    bool gausflip = false;
    if (this->inputs->highpt < 4 + 1e-3)
      gausflip = true;

    if (gausflip) {
      double A = this->fit->GetParameter(3);
      double sigma = this->fit->GetParameter(4);
      double B = this->fit->GetParameter(0);
      double rho = this->fit->GetParameter(2);
      this->fit->SetParameter(0, A);
      this->fit->SetParameter(2, sigma);
      this->fit->SetParameter(3, B);
      this->fit->SetParameter(4, rho);
    }
  }
  if (this->inputs->hadron == "K0S" && this->inputs->fitType == InputSettings::kPol2GausGaus) {
    bool gausflip = false;
    if (this->inputs->highpt > 1. + 1e-3 && this->inputs->highpt < 30. + 1e-3)
      gausflip = true;

    if (gausflip) {
      double A = this->fit->GetParameter(3);
      double sigma = this->fit->GetParameter(4);
      double B = this->fit->GetParameter(0);
      double rho = this->fit->GetParameter(2);
      this->fit->SetParameter(0, A);
      this->fit->SetParameter(2, sigma);
      this->fit->SetParameter(3, B);
      this->fit->SetParameter(4, rho);
    }
  }
}

TF1* MassFitter::loadFitFunction() {
  string saveName = this->inputs->getSaveNameFromPt("fit");
  string expression = this->inputs->getFitExpression();
  TF1* f = new TF1(saveName.c_str(), expression.c_str(), this->inputs->fitmin, this->inputs->fitmax);
  this->fit = (TF1*)f->Clone();

  if (inputs->verbosity >= InputSettings::kDebug)
    f->Print();

  return f;
}

TF1* MassFitter::loadSavedFitFunction() {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    string s = "Could not open file " + inputs->inputFileName;
    inputs->printLog(s, InputSettings::kErrors);
    return nullptr;
  }

  string fitName = inputs->getSaveNameFromPt("fit");
  TF1* f = (TF1*)file->Get(fitName.c_str());
  if (!f) {
    string s = "Could not find fit function " + fitName + " in file " + inputs->inputFileName;
    inputs->printLog(s, InputSettings::kErrors);
    return nullptr;
  }

  this->fit = (TF1*)f->Clone();
  return f;
}

TH1* MassFitter::loadSavedMassHist() {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    string s = "Could not open file " + inputs->inputFileName;
    inputs->printLog(s, InputSettings::kErrors);
    return nullptr;
  }

  string histName = inputs->histName;
  if (histName == "") {
    histName = inputs->getSaveNameFromPt("data");
  }

  TH1* hist = (TH1*)file->Get(histName.c_str());
  if (!hist) {
    string s = "Could not find histogram " + inputs->histName + " in file " + inputs->inputFileName;
    inputs->printLog(s, InputSettings::kErrors);
    return nullptr;
  }

  this->data = (TH1*)hist->Clone();
  return hist;
}

TH1* MassFitter::loadMassHist() {
  TFile* f = TFile::Open(this->inputs->inputFileName.c_str(), "READ");
  if (!f) {
    string s = "Could not open file " + this->inputs->inputFileName;
    this->inputs->printLog(s, InputSettings::kErrors);
    return nullptr;
  }

  if (this->inputs->histName == "") {
    this->setHistNameFromTrain();
    string s = "MassFitter::loadMassHist(): histName not set, setting to: " + this->inputs->histName;
    this->inputs->printLog(s, InputSettings::kInfo);
  }
  THnSparse* hist = (THnSparse*)f->Get(this->inputs->histName.c_str());
  if (!hist) {
    string s = "Could not find histogram " + this->inputs->histName + " in file " + this->inputs->inputFileName;
    this->inputs->printLog(s, InputSettings::kErrors);
    return nullptr;
  }

  int projectionAxis = 1;
  if (this->inputs->hadron == "Lambda")
    projectionAxis = 2;
  if (this->inputs->hadron == "AntiLambda")
    projectionAxis = 3;

  array<int, 2> bins = histutils::getProjectionBins(hist->GetAxis(0), this->inputs->ptmin, this->inputs->ptmax);
  int minBin = bins[0], maxBin = bins[1];
  this->inputs->lowpt = hist->GetAxis(0)->GetBinLowEdge(minBin);
  this->inputs->highpt = hist->GetAxis(0)->GetBinUpEdge(maxBin);

  hist->GetAxis(0)->SetRange(minBin, maxBin);
  TH1* h = hist->Projection(projectionAxis);

  if (this->inputs->rebinNumber < 0)
    this->setRebinNumberFromHadronAndPt();

  if (this->inputs->rebinNumber > 1)
    h->Rebin(this->inputs->rebinNumber);

  if (this->inputs->normaliseData) {
    bins = histutils::getProjectionBins(h->GetXaxis(), this->inputs->massWindowMin, this->inputs->massWindowMax);
    minBin = bins[0], maxBin = bins[1];
    // h->Scale(1./getScaleInRange(h, minBin, maxBin, false, false));
    double histScalse = histutils::getUpperBoundInRange(h, minBin, maxBin, false);
    h->Scale(1./histScalse);
  }

  string hName = this->inputs->getSaveNameFromPt("data");
  this->data = (TH1*)h->Clone(hName.c_str());

  if (inputs->verbosity >= InputSettings::kDebug)
    h->Print();

  return h;
}

void MassFitter::setFitInitialValues() {
  array<int, 2> massBins = histutils::getProjectionBins(this->data->GetXaxis(), this->inputs->massWindowMin, this->inputs->massWindowMax);
  int minBin = massBins[0], maxBin = massBins[1];
  int extremeBin = histutils::getUpperBoundBinInRange(this->data, minBin, maxBin, false);
  // int extremeBin = getExtremeBinInRange(this->data, minBin, maxBin, false, false);

  int b0 = this->data->GetXaxis()->FindBin(this->inputs->polInitx0);
  int b1 = this->data->GetXaxis()->FindBin(this->inputs->polInitx1);
  int b2 = this->data->GetXaxis()->FindBin(this->inputs->polInitx2);
  double x0 = this->data->GetBinLowEdge(b0);
  double x1 = this->data->GetBinLowEdge(b1);
  double x2 = this->data->GetBinLowEdge(b2);
  double y0 = this->data->GetBinContent(b0);
  double y1 = this->data->GetBinContent(b1);
  double y2 = this->data->GetBinContent(b2);

  double A, mu, sigma; // Signal Gaussian
  double Gamma; // Signal Breit-Wigner
  double lambda; // Gaus->Exp crossover in sigma
  double B, rho; // Second Gaussian
  double C, nu, tau; // Xex
  double a, b, c; // Polynomial background

  A = this->data->GetBinContent(extremeBin);
  mu = this->data->GetXaxis()->GetBinCenter(extremeBin);
  B = 1e-2 * A;
  C = 1e-2 * A;

  if (this->inputs->hadron == "K0S") {
    Gamma = 1e-2;
    sigma = 1e-2;
    lambda = 2.;
    rho = 2 * sigma;
    nu = 0.5;
    tau = 6e-3;
  } else {
    Gamma = 1e-3;
    sigma = 1e-3;
    lambda = 2.;
    rho = 2 * sigma;
    nu = 1.2;
    tau = 6e-3;
  }

  if (this->inputs->fitType % 2 == 0) { // pol1
    b = (y1-y0)/(x1-x0);
    a = y0 - b*x0;
  } else { // pol2
    c = ( (y2-y1)/(x2-x1) - (y1-y0)/(x1-x0) ) / (x2-x0);
    b = (y1-y0)/(x1-x0) - c*(x1+x0);
    a = y0 - b*x0 - c*x0*x0;
  }

  switch (this->inputs->fitType) {
    case InputSettings::kPol1BreitWigner:
      this->fit->SetParameters(A, mu, Gamma, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0, 1);
      this->fit->SetParLimits(4, -1, 1);
      break;
    case InputSettings::kPol2BreitWigner:
      this->fit->SetParameters(A, mu, Gamma, a, b, c);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0, 1);
      this->fit->SetParLimits(4, -1, 1);
      // this->fit->SetParLimits(5, -1, 1);
      break;
    case InputSettings::kPol1BreitWignerXex:
      this->fit->SetParameters(A, mu, Gamma, C, nu, tau, a, b);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 1.); // C
      this->fit->SetParLimits(4, 0.5, 0.6); // nu
      // this->fit->SetParLimits(5, 0., 1.); // tau
      this->fit->SetParLimits(6, 0., 1.); // a
      this->fit->SetParLimits(7, -1, 1); // b
      break;
    case InputSettings::kPol2BreitWignerXex:
      this->fit->SetParameters(A, mu, Gamma, C, nu, tau, a, b);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 1.); // C
      this->fit->SetParLimits(4, 0.5, 0.6); // nu
      // this->fit->SetParLimits(5, 0., 1.); // tau
      this->fit->SetParLimits(6, 0., 1.); // a
      this->fit->SetParLimits(7, -1, 1); // b
      break;
    case InputSettings::kPol1ExpGausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, lambda, a, b);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 10.); // lambda_L
      this->fit->SetParLimits(4, 0., 10.); // lambda_R
      this->fit->SetParLimits(5, 0., 1.); // a
      this->fit->SetParLimits(6, -1., 1.); // b
      break;
    case InputSettings::kPol2ExpGausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, lambda, a, b, c);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 10.); // lambda_L
      this->fit->SetParLimits(4, 0., 10.); // lambda_R
      this->fit->SetParLimits(5, 0., 1.); // a
      this->fit->SetParLimits(6, -1., 1.); // b
      // this->fit->SetParLimits(7, -1., 0); // c
      break;
    case InputSettings::kPol1GausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 10.);
      this->fit->SetParLimits(4, 0, 1);
      this->fit->SetParLimits(5, -1, 1);
      break;
    case InputSettings::kPol2GausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, a, b, c);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 10.); // lambda
      this->fit->SetParLimits(4, 0, 1); // a
      this->fit->SetParLimits(5, -1, 1); // b
      this->fit->SetParLimits(6, -1, 0); //c
      break;
    case InputSettings::kPol1GausGaus:
      this->fit->SetParameters(A, mu, sigma, B, rho, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 1.);
      this->fit->SetParLimits(4, 0., 1.);
      this->fit->SetParLimits(5, 0., 1.);
      this->fit->SetParLimits(6, -1, 1);
      break;
      case InputSettings::kPol2GausGaus:
      this->fit->SetParameters(A, mu, sigma, B, rho, a, b, c);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 1.);
      this->fit->SetParLimits(4, 0., 1.);
      this->fit->SetParLimits(5, 0., 1.);
      this->fit->SetParLimits(6, -1., 1.);
      break;
    case InputSettings::kPol1GausGausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, B, rho, a, b);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 10.); // lambda
      this->fit->SetParLimits(4, 0., 1.); // B
      this->fit->SetParLimits(5, 0., 1.); // rho
      this->fit->SetParLimits(6, 0., 1.); // a
      this->fit->SetParLimits(7, -1., 1.); // b
      break;
    case InputSettings::kPol2GausGausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, B, rho, a, b, c);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 10.); // lambda
      this->fit->SetParLimits(4, 0., 1.); // B
      this->fit->SetParLimits(5, 0., 1.); // rho
      this->fit->SetParLimits(6, 0., 1.); // a
      this->fit->SetParLimits(7, -1., 1.); // b
      // this->fit->SetParLimits(8, -1., 0.); // c
      break;
    case InputSettings::kPol1GausGausXex:
      this->fit->SetParameters(A, mu, sigma, B, rho, C, nu, tau, a, b);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 1.); // B
      this->fit->SetParLimits(4, 2e-2, 1.); // rho
      this->fit->SetParLimits(5, 0., 100.); // C
      this->fit->SetParLimits(6, 0.5, 0.6); // nu
      this->fit->SetParLimits(7, 0., 1e-2); // tau
      this->fit->SetParLimits(8, 0., 1.); // a
      this->fit->SetParLimits(9, -1., 1.); // b
      break;
    case InputSettings::kPol2GausGausXex:
      this->fit->SetParameters(A, mu, sigma, B, rho, C, nu, tau, a, b);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 1.); // B
      this->fit->SetParLimits(4, 0., 1.); // rho
      this->fit->SetParLimits(5, 0., 100.); // C
      this->fit->SetParLimits(6, 0.5, 0.6); // nu
      this->fit->SetParLimits(7, 0., 1e-2); // tau
      this->fit->SetParLimits(8, 0., 1.); // a
      this->fit->SetParLimits(9, -1., 1.); // b
      break;
    case InputSettings::kPol1GausXex:
      this->fit->SetParameters(A, mu, sigma, C, nu, tau, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 100.); // C
      this->fit->SetParLimits(4, 0.5, 0.6); // nu
      this->fit->SetParLimits(5, 0., 1e-2); // tau
      this->fit->SetParLimits(6, 0., 1.); // a
      this->fit->SetParLimits(7, -1., 1.); // b
      break;
    case InputSettings::kPol2GausXex:
      this->fit->SetParameters(A, mu, sigma, C, nu, tau, a, b, c);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 100.); // C
      this->fit->SetParLimits(4, 0.5, 0.6); // nu
      this->fit->SetParLimits(5, 0., 1e-2); // tau
      this->fit->SetParLimits(6, 0., 1.); // a
      this->fit->SetParLimits(7, -1., 1.); // b
      break;
    case InputSettings::kPol1Voigt:
      this->fit->SetParameters(A, mu, sigma, Gamma, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 1.);
      this->fit->SetParLimits(4, 0., 1.);
      this->fit->SetParLimits(5, -1., 1.);
      break;
    case InputSettings::kPol2Voigt:
      this->fit->SetParameters(A, mu, sigma, Gamma, a, b, c);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs->fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs->massWindowMin, this->inputs->massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 0.005);
      this->fit->SetParLimits(3, 0., 0.005);
      this->fit->SetParLimits(4, 0., 1.);
      this->fit->SetParLimits(5, -1., 1.);
      break;
    default:
      string s = "MassFitter::setFitInitialValues() Error: do not have initial values for function " + this->inputs->fitName;
      this->inputs->printLog(s, InputSettings::kErrors);
      return;
  }

  for (int iP = 0; iP < this->fit->GetNpar(); iP++) {
    bool parChanged = false;
    double parameter = this->fit->GetParameter(iP);
    double parMin = -1e6, parMax = 1e6;
    this->fit->GetParLimits(iP, parMin, parMax);
    if (parameter < parMin) {
      parameter = parMin;
      parChanged = true;
    }
    if (parameter > parMax) {
      parameter = parMax;
      parChanged = true;
    }
    if (parChanged) {
      string s = "MassFitter::setFitInitialValues() Warning: parameter " + to_string(iP) + " outside limits, changed to " + to_string(parameter);
      this->inputs->printLog(s, InputSettings::kWarnings);
      this->fit->SetParameter(iP, parameter);
    }
  }
}

vector<TF1*> MassFitter::loadSavedFitParts() {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    string s = "Could not open file " + inputs->inputFileName;
    inputs->printLog(s, InputSettings::kErrors);
    return {};
  }
  inputs->printLog("MassFitter::loadSavedFitParts() opening file " + inputs->inputFileName, InputSettings::kDebug);

  vector<TF1*> parts = {};
  int maxParts = 10; // To prevent infinite loop
  for (int i = 1; i < maxParts; i++) {
    inputs->printLog("Getting fit part " + inputs->getSaveNameFromPt("fit", "_f" + to_string(i)), InputSettings::kDebugMax);
    string fName = inputs->getSaveNameFromPt("fit", "_f" + to_string(i));
    TF1* f = (TF1*)file->Get(fName.c_str());
    if (!f) break;

    parts.push_back(f);
  }

  inputs->printLog("MassFitter::loadSavedFitParts() loaded " + to_string(parts.size()) + " fit parts", InputSettings::kDebug);

  fitParts = parts;
  return parts;
}

vector<TF1*> MassFitter::loadFitParts() {
  if (!this->fit)
    return {};

  string fName = this->fit->GetName();
  vector<TF1*> parts;

  if (this->inputs->fitType == InputSettings::kPol1BreitWigner || this->inputs->fitType == InputSettings::kPol2BreitWigner) {
    double A, mu, Gamma;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    Gamma = this->fit->GetParameter(2);
    a = this->fit->GetParameter(3);
    b = this->fit->GetParameter(4);
    if (this->inputs->fitType == InputSettings::kPol2BreitWigner) c = this->fit->GetParameter(5);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "breitwigner(x,[0],[1],[2])", this->inputs->fitmin, this->inputs->fitmax);
    f1->SetParameters(A, mu, Gamma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    if (this->inputs->fitType == InputSettings::kPol1BreitWigner) {
      TF1* f2 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs->fitmin, this->inputs->fitmax);
      f2->SetParameters(a, b);
      parts.push_back(f2);
    } else {
      TF1* f2 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs->fitmin, this->inputs->fitmax);
      f2->SetParameters(a, b, c);
      parts.push_back(f2);
    }
  } else if (this->inputs->fitType == InputSettings::kPol1BreitWignerXex || this->inputs->fitType == InputSettings::kPol2BreitWignerXex) {
    double A, mu, Gamma;
    double C, nu, tau;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    Gamma = this->fit->GetParameter(2);
    C = this->fit->GetParameter(3);
    nu = this->fit->GetParameter(4);
    tau = this->fit->GetParameter(5);
    a = this->fit->GetParameter(6);
    b = this->fit->GetParameter(7);
    if (this->inputs->fitType == InputSettings::kPol2BreitWignerXex) c = this->fit->GetParameter(8);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "breitwigner(x,[0],[1],[2])", this->inputs->fitmin, this->inputs->fitmax);
    f1->SetParameters(A, mu, Gamma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Exp(-[1]*(x-[2])*(x-[2]))", this->inputs->fitmin, this->inputs->fitmax);
    f2->SetParameters(C, nu, tau);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    if (this->inputs->fitType == InputSettings::kPol1BreitWignerXex) {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs->fitmin, this->inputs->fitmax);
      f3->SetParameters(a, b);
      parts.push_back(f3);
    } else {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs->fitmin, this->inputs->fitmax);
      f3->SetParameters(a, b, c);
      parts.push_back(f3);
    }
  } else if (this->inputs->fitType == InputSettings::kPol1ExpGausExp || this->inputs->fitType == InputSettings::kPol2ExpGausExp) {
    double A, mu, sigma, lambda_L, lambda_R;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    sigma = this->fit->GetParameter(2);
    lambda_L = this->fit->GetParameter(3);
    lambda_R = this->fit->GetParameter(4);
    a = this->fit->GetParameter(5);
    b = this->fit->GetParameter(6);
    if (this->inputs->fitType == InputSettings::kPol2ExpGausExp) c = this->fit->GetParameter(7);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2]) * (x > ([1]-[2]*[3])) * (x < ([1]+[2]*[4]))", mu - sigma*lambda_L, mu+ sigma*lambda_R);
    f1->SetParameters(A, mu, sigma, lambda_L, lambda_R);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Exp((x - [1] + [2]*[3]/2.)/([2]/[3])) * (x <= [1]-[2]*[3])", this->inputs->fitmin, mu - sigma*lambda_L);
    f2->SetParameters(A, mu, sigma, lambda_L);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    TF1* f3 = new TF1(s.c_str(), "[0]*TMath::Exp(-(x - [1] - [2]*[3]/2.)/([2]/[3])) * (x >= [1]+[2]*[3])", mu + sigma*lambda_R, this->inputs->fitmax);
    f3->SetParameters(A, mu, sigma, lambda_R);
    parts.push_back(f3);

    s = TString::Format("%s_%s", fName.c_str(), "f4").Data();
    if (this->inputs->fitType == InputSettings::kPol1ExpGausExp) {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs->fitmin, this->inputs->fitmax);
      f4->SetParameters(a, b);
      parts.push_back(f4);
    } else {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs->fitmin, this->inputs->fitmax);
      f4->SetParameters(a, b, c);
      parts.push_back(f4);
    }
  } else if (this->inputs->fitType == InputSettings::kPol1GausExp || this->inputs->fitType == InputSettings::kPol2GausExp) {
    double A, mu, sigma, lambda;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    sigma = this->fit->GetParameter(2);
    lambda = this->fit->GetParameter(3);
    a = this->fit->GetParameter(4);
    b = this->fit->GetParameter(5);
    if (this->inputs->fitType == InputSettings::kPol2GausExp) c = this->fit->GetParameter(6);

    // TODO: Change this to 2 functions: G->e and pol1/2
    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3]))", this->inputs->fitmin, mu+sigma*lambda);
    f1->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Exp(-(x - [1] - [2]*[3]/2.)/([2]/[3])) * (x >= [1]+[2]*[3])", mu+sigma*lambda, this->inputs->fitmax);
    f2->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    if (this->inputs->fitType == InputSettings::kPol1GausExp) {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs->fitmin, this->inputs->fitmax);
      f3->SetParameters(a, b);
      parts.push_back(f3);
    } else {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs->fitmin, this->inputs->fitmax);
      f3->SetParameters(a, b, c);
      parts.push_back(f3);
    }
  } else if (this->inputs->fitType == InputSettings::kPol1GausGaus || this->inputs->fitType == InputSettings::kPol2GausGaus) {
    double A, mu, sigma;
    double B, lambda;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    sigma = this->fit->GetParameter(2);
    B = this->fit->GetParameter(3);
    lambda = this->fit->GetParameter(4);
    a = this->fit->GetParameter(5);
    b = this->fit->GetParameter(6);
    if (this->inputs->fitType == InputSettings::kPol2GausGaus) c = this->fit->GetParameter(7);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs->fitmin, this->inputs->fitmax);
    f1->SetParameters(A, mu, sigma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs->fitmin, this->inputs->fitmax);
    f2->SetParameters(B, mu, lambda);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    if (this->inputs->fitType == InputSettings::kPol1GausGaus) {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs->fitmin, this->inputs->fitmax);
      f3->SetParameters(a, b);
      parts.push_back(f3);
    } else {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs->fitmin, this->inputs->fitmax);
      f3->SetParameters(a, b, c);
      parts.push_back(f3);
    }
  } else if (this->inputs->fitType == InputSettings::kPol1GausGausExp || this->inputs->fitType == InputSettings::kPol2GausGausExp) {
    double A, mu, sigma, lambda;
    double B, rho;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    sigma = this->fit->GetParameter(2);
    lambda = this->fit->GetParameter(3);
    B = this->fit->GetParameter(4);
    rho = this->fit->GetParameter(5);
    a = this->fit->GetParameter(6);
    b = this->fit->GetParameter(7);
    if (this->inputs->fitType == InputSettings::kPol2GausExp) c = this->fit->GetParameter(8);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3]))", this->inputs->fitmin, mu+sigma*lambda);
    f1->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Exp(-(x - [1] - [2]*[3]/2.)/([2]/[3])) * (x >= [1]+[2]*[3])", mu+sigma*lambda, this->inputs->fitmax);
    f2->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    TF1* f3 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs->fitmin, this->inputs->fitmax);
    f3->SetParameters(B, mu, rho);
    parts.push_back(f3);

    s = TString::Format("%s_%s", fName.c_str(), "f4").Data();
    if (this->inputs->fitType == InputSettings::kPol1GausGausExp) {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs->fitmin, this->inputs->fitmax);
      f4->SetParameters(a, b);
      parts.push_back(f4);
    } else {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs->fitmin, this->inputs->fitmax);
      f4->SetParameters(a, b, c);
      parts.push_back(f4);
    }
  } else if (this->inputs->fitType == InputSettings::kPol1GausGausXex || this->inputs->fitType == InputSettings::kPol2GausGausXex) {
    double A, mu, sigma;
    double B, rho;
    double C, nu, tau;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    sigma = this->fit->GetParameter(2);
    B = this->fit->GetParameter(3);
    rho = this->fit->GetParameter(4);
    C = this->fit->GetParameter(5);
    nu = this->fit->GetParameter(6);
    tau = this->fit->GetParameter(7);
    a = this->fit->GetParameter(8);
    b = this->fit->GetParameter(9);
    if (this->inputs->fitType == InputSettings::kPol2GausGausXex) c = this->fit->GetParameter(10);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs->fitmin, this->inputs->fitmax);
    f1->SetParameters(A, mu, sigma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs->fitmin, this->inputs->fitmax);
    f2->SetParameters(B, mu, rho);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    TF1* f3 = new TF1(s.c_str(), "max(0., [0]*(x-[1])*TMath::Exp(-(x-[1])/[2]))", this->inputs->fitmin, this->inputs->fitmax);
    f3->SetParameters(C, nu, tau);
    parts.push_back(f3);

    s = TString::Format("%s_%s", fName.c_str(), "f4").Data();
    if (this->inputs->fitType == InputSettings::kPol1GausGausXex) {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs->fitmin, this->inputs->fitmax);
      f4->SetParameters(a, b);
      parts.push_back(f4);
    } else {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs->fitmin, this->inputs->fitmax);
      f4->SetParameters(a, b, c);
      parts.push_back(f4);
    }
  } else if (this->inputs->fitType == InputSettings::kPol1GausXex || this->inputs->fitType == InputSettings::kPol2GausXex) {
    double A, mu, sigma;
    double B, nu, tau;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    sigma = this->fit->GetParameter(2);
    B = this->fit->GetParameter(3);
    nu = this->fit->GetParameter(4);
    tau = this->fit->GetParameter(5);
    a = this->fit->GetParameter(6);
    b = this->fit->GetParameter(7);
    if (this->inputs->fitType == InputSettings::kPol2GausXex) c = this->fit->GetParameter(8);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs->fitmin, this->inputs->fitmax);
    f1->SetParameters(A, mu, sigma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "max(0., [0]*(x-[1])*TMath::Exp(-(x-[1])/[2]))", this->inputs->fitmin, this->inputs->fitmax);
    f2->SetParameters(B, nu, tau);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    if (this->inputs->fitType == InputSettings::kPol1GausXex) {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs->fitmin, this->inputs->fitmax);
      f3->SetParameters(a, b);
      parts.push_back(f3);
    } else {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs->fitmin, this->inputs->fitmax);
      f3->SetParameters(a, b, c);
      parts.push_back(f3);
    }
  } else if (this->inputs->fitType == InputSettings::kPol1Voigt || this->inputs->fitType == InputSettings::kPol2Voigt) {
    double A, mu, sigma, Gamma;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    sigma = this->fit->GetParameter(2);
    Gamma = this->fit->GetParameter(3);
    a = this->fit->GetParameter(4);
    b = this->fit->GetParameter(5);
    if (this->inputs->fitType == InputSettings::kPol2Voigt) c = this->fit->GetParameter(6);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Voigt(x-[1],[2],[3])", this->inputs->fitmin, this->inputs->fitmax);
    f1->SetParameters(A, mu, sigma, Gamma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    if (this->inputs->fitType == InputSettings::kPol1Voigt) {
      TF1* f2 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs->fitmin, this->inputs->fitmax);
      f2->SetParameters(a, b);
      parts.push_back(f2);
    } else {
      TF1* f2 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs->fitmin, this->inputs->fitmax);
      f2->SetParameters(a, b, c);
      parts.push_back(f2);
    }
  } else {
    string s = "InputSettings::setFitParts() Error: do not know the parts of this function";
    this->inputs->printLog(s, InputSettings::kErrors);
  }
  this->fitParts = parts;
  return parts;
}

TH1* MassFitter::loadFitParams() {
  if (!this->fit)
    return nullptr;

  string saveName = this->inputs->getSaveNameFromPt("fitParams");

  int nPars = this->fit->GetNpar();
  TH1* h = new TH1D(saveName.c_str(), "fitParams", nPars, -0.5, 1.*nPars - 0.5);
  h->SetTitle("Fit parameters");
  h->GetXaxis()->SetTitle("");
  h->GetYaxis()->SetTitle("Value");
  for (int i = 0; i < nPars; i++) {
    h->GetXaxis()->SetBinLabel(i + 1, this->fit->GetParName(i));
    h->SetBinContent(i + 1, this->fit->GetParameter(i));
  }
  this->fitParams = (TH1*)h->Clone();
  return h;
}

// TODO: This should maybe be moved to SignalFinder and done separately of the fits
TH1* MassFitter::loadFitResults() {
  if (!this->fit)
    return nullptr;

  string saveName = this->inputs->getSaveNameFromPt("fitResults");

  int nBins = 15;
  TH1* h = new TH1D(saveName.c_str(), "fitResults", nBins, -0.5, 1.*nBins - 0.5);
  h->SetTitle("Fit results");
  h->GetXaxis()->SetTitle("");
  h->GetXaxis()->SetBinLabel(1, "#chi^{2} (fit)");
  h->GetXaxis()->SetBinLabel(2, "NDF (fit)");
  h->GetXaxis()->SetBinLabel(3, "#chi^{2}/NDF (fit)");
  h->GetXaxis()->SetBinLabel(4, "#chi^{2} (all)");
  h->GetXaxis()->SetBinLabel(5, "NDF (all)");
  h->GetXaxis()->SetBinLabel(6, "#chi^{2}/NDF (all)");
  h->GetXaxis()->SetBinLabel(7, "Sig+Bkg [hist]");
  h->GetXaxis()->SetBinLabel(8, "Sig [hist]");
  h->GetXaxis()->SetBinLabel(9, "Bkg [hist]");
  h->GetXaxis()->SetBinLabel(10, "Sig/(Sig+Bkg) [hist]");
  h->GetXaxis()->SetBinLabel(11, "Sig/Bkg [hist]");
  h->GetXaxis()->SetBinLabel(12, "Sig/sqrt(Sig+Bkg) [hist]");
  h->GetXaxis()->SetBinLabel(13, "Sig+Bkg [fit]");
  h->GetXaxis()->SetBinLabel(14, "Sig [fit]");
  h->GetXaxis()->SetBinLabel(15, "Bkg [fit]");

  double chisqFit = chisqInRange(this->data, this->fit, this->inputs->fitmin, this->inputs->fitmax);
  double ndfFit = ndfInRange(this->data, this->fit, this->inputs->fitmin, this->inputs->fitmax);
  double chisqAll = chisqInRange(this->data, this->fit, this->data->GetXaxis()->GetXmin(), this->data->GetXaxis()->GetXmax());
  double ndfAll = ndfInRange(this->data, this->fit, this->data->GetXaxis()->GetXmin(), this->data->GetXaxis()->GetXmax());
  this->calcSigBkg();

  h->SetBinContent(1, chisqFit);
  h->SetBinContent(2, ndfFit);
  h->SetBinContent(3, chisqFit / ndfFit);
  h->SetBinContent(4, chisqAll);
  h->SetBinContent(5, ndfAll);
  h->SetBinContent(6, chisqAll / ndfAll);
  h->SetBinContent(7, this->signalPlusBackground[0]);
  h->SetBinContent(8, this->signal[0]);
  h->SetBinContent(9, this->background[0]);
  h->SetBinContent(10, this->signal[0] / this->signalPlusBackground[0]);
  h->SetBinContent(11, this->signal[0] / this->background[0]);
  h->SetBinContent(12, this->signal[0] / sqrt(this->signalPlusBackground[0]));
  h->SetBinContent(13, this->signalPlusBackground[1]);
  h->SetBinContent(14, this->signal[1]);
  h->SetBinContent(15, this->background[1]);

  this->fitResults = (TH1*)h->Clone(saveName.c_str());
  return h;
}

TH1* MassFitter::loadResidualHist() {
  if (!this->data || !this->fit)
    return nullptr;

  string saveName = this->inputs->getSaveNameFromPt("residual");

  TH1* h = (TH1*)this->data->Clone(saveName.c_str());
  TF1* f = (TF1*)this->fit->Clone("function");
  f->SetRange(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax()); // Make sure function is defined everywhere
  h->Add(f, -1.);
  this->residual = (TH1*)h->Clone();
  return h;
}

TH1* MassFitter::loadPullHist() {
  if (!this->residual)
    this->loadResidualHist();

  string saveName = this->inputs->getSaveNameFromPt("pull");

  TH1* h = (TH1*)this->residual->Clone(saveName.c_str());
  h->Reset();
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double bc = this->residual->GetBinContent(i);
    double be = this->data->GetBinError(i);
    double p = bc / be;
    h->SetBinContent(i, p);
  }
  this->pull = (TH1*)h->Clone();
  return h;
}

string MassFitter::setHistNameFromTrain() {
  if (inputs->inputIssue("setHistNameFromTrain", "hadron"))
    return "";

  string s = "jet-fragmentation";
  if (inputs->train == 287744) s += "_id12406";
  if (inputs->train == 419996) s += "_id28293";
  if (inputs->train == 426828) s += "_id28293";

  s += "/data/V0/V0PtMass";

  inputs->histName = s;
  return s;
}

int MassFitter::setRebinNumberFromHadronAndPt() {
  int r = -1;

  if (this->inputs->hadron == "K0S" && this->inputs->lowpt > 25 - 1e3)
    r = 2;

  this->inputs->rebinNumber = r;
  return r;
}

void MassFitter::setSignalRegionFromSigma() {
  if (this->fitParts.size() == 0) {
    string s = "MassFitter::setSignalRegionFromSigma() Error: fit parts not set";
    this->inputs->printLog(s, InputSettings::kErrors);
    return;
  }

  TF1* f = this->fitParts[0];
  double mu = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  double nSigma = this->inputs->nSigma;
  this->inputs->setSignalRegion(mu - nSigma * sigma, mu + nSigma * sigma);
}

void MassFitter::writeOutputsToFile() {
  this->inputs->writeOutputToFile(this->data);
  this->inputs->writeOutputToFile(this->fit);
  for (auto f : this->fitParts) {
    this->inputs->writeOutputToFile(f);
  }
  this->inputs->writeOutputToFile(this->fitResults);
  this->inputs->writeOutputToFile(this->fitParams);
  this->inputs->writeOutputToFile(this->pull);
  this->inputs->writeOutputToFile(this->residual);
}

// -------------------------------------------------------------------------------------------------
//
// Struct for summarising fit information into TH2
//
// -------------------------------------------------------------------------------------------------

struct FitSummariser {
  private:
  public:
    InputSettings* inputs;

    FitSummariser() { inputs = new InputSettings(); }
    FitSummariser(InputSettings& x) { inputs = &x; }

    void summariseFitInfo();
    void summariseFitInfo(string inputHistName);
};

void FitSummariser::summariseFitInfo() {
  this->summariseFitInfo(this->inputs->histName);
}

// Puts fit info from all pt bins into a single TH2 histogram
void FitSummariser::summariseFitInfo(string inputHistName) {
  int nPtBins = this->inputs->ptBinEdges.size();
  int nFitInfos = -1;
  if (inputHistName == "fitParams") {
    TF1* f = new TF1("f", this->inputs->getFitExpression().c_str(), 0., 1.);
    nFitInfos = f->GetNpar();
  } else if (inputHistName == "fitResults") {
    nFitInfos = 12;
  } else {
    string s = "summariseFitInfo: unknown hist name " + inputHistName;
    this->inputs->printLog(s, InputSettings::kErrors);
    return;
  }

  TH2* fitInfo = new TH2D(inputHistName.c_str(), inputHistName.c_str(), nPtBins, -0.5, nPtBins - 0.5, nFitInfos, -0.5, nFitInfos - 0.5);
  for (int ix = 0; ix < nPtBins; ix++) {
    fitInfo->GetXaxis()->SetBinLabel(ix + 1, TString::Format("%.1f < #it{p}_{T, V0} < %.1f", this->inputs->ptBinEdges[ix][0], this->inputs->ptBinEdges[ix][1]).Data());
  }
  bool binLabelsSet = false;

  TFile* file = TFile::Open(this->inputs->inputFileName.c_str(), "READ");
  if (!file) {
    string s = "Could not open file " + this->inputs->inputFileName;
    this->inputs->printLog(s, InputSettings::kErrors);
    return;
  }

  for (int iPt = 0; iPt < nPtBins; iPt++) {
    this->inputs->setPt(this->inputs->ptBinEdges[iPt][0], this->inputs->ptBinEdges[iPt][1]);
    string histName = this->inputs->getSaveNameFromPt(inputHistName);
    TH1* hist = (TH1*)file->Get(histName.c_str());
    if (!hist) {
      string s = "Could not find histogram " + histName + " in file " + this->inputs->inputFileName;
      this->inputs->printLog(s, InputSettings::kErrors);
      continue;
    }

    for (int iInfo = 0; iInfo < nFitInfos; iInfo++) {
      double parValue = hist->GetBinContent(iInfo + 1);
      fitInfo->SetBinContent(iPt + 1, iInfo + 1, parValue);

      if (binLabelsSet)
        continue;

      fitInfo->GetYaxis()->SetBinLabel(iInfo + 1, hist->GetXaxis()->GetBinLabel(iInfo + 1));
    }
    binLabelsSet = true;
  }

  file->Close();
  this->inputs->writeOutputToFile(fitInfo);
}

// -------------------------------------------------------------------------------------------------
//
// Struct for determining signal region and fraction
//
// -------------------------------------------------------------------------------------------------

struct SignalFinder {
  private:
  public:
    // Use pointers so inputs = mf->inputs
    InputSettings* inputs;
    MassFitter* mf;

    int nBins = 100;
    double xMin = 0, xMax = 10;
    TH1* hSigFrac = nullptr;
    TH1* hSigFracL = nullptr;
    TH1* hSigFracR = nullptr;

    double signalFraction = -1, signalFractionLeft = -1, signalFractionRight = -1;
    double hSig, hBkg, hSigPlusBkg; // Signal from histogram
    double fSig, fBkg, fSigPlusBkg; // Signal from fit

    vector<int> bkgFits = {};

    SignalFinder() { inputs = new InputSettings(); mf = new MassFitter(*inputs); }
    SignalFinder(InputSettings& x) { inputs = &x; mf = new MassFitter(x); }
    SignalFinder(MassFitter& x) { mf = &x; inputs = x.inputs; }

    // Calculate the signal and background
    void calcSigBkg();
    void setBkgFits();
    void getSignalFraction(double n);
    array<double, 2> getSignalRegion(bool useSigma, bool set);

    // Defining the signal region and finding the expected fraction of total signal
    double getGG_Sigma(TF1* f); // Characteristic length

    // Signal region
    array<double, 2> GG_SigRegionFromSteps(double n, TF1* f = nullptr);
    array<double, 2> GG_SigRegionFromFrac(double s, TF1* f = nullptr);
    array<double, 2> GGE_SigRegionFromSteps(double n, TF1* f = nullptr);
    array<double, 2> GGE_SigRegionFromFrac(double s, TF1* f = nullptr);
    array<double, 2> EGE_SigRegionFromSteps(double n, TF1* f = nullptr);
    array<double, 2> EGE_SigRegionFromFrac(double s, TF1* f = nullptr);

    // Signal fraction
    double GG_SigFrac(double n, TF1* f = nullptr);
    double GGE_SigFrac(double n, TF1* f = nullptr, bool leftSide = false);
    double EGE_SigFrac(double n, TF1* f = nullptr, bool leftSide = false);

    // Make and load signal fraction histograms
    TH1* GG_makeSigFracHist(TF1* f);
    TH1* GG_loadSigFracHist();
    TH1* GGE_makeSigFracHist(TF1* f, bool leftSide);
    TH1* GGE_loadSigFracHist(bool leftSide);
    TH1* EGE_makeSigFracHist(TF1* f, bool leftSide);
    TH1* EGE_loadSigFracHist(bool leftSide);

    // Select data within signal region
    TH1* makeSigRegionHist();
};

void SignalFinder::calcSigBkg() {
  if (mf->fitParts.empty()) {
    inputs->printLog("SignalFinder::calcSigBkg: fit parts not set, loading saved fit parts", InputSettings::kInfo);
    mf->loadSavedFitParts();
  }

  array<int, 2> sigBins = histutils::getProjectionBins(mf->data->GetXaxis(), inputs->signalRegionMin, inputs->signalRegionMax);
  double xmin = mf->data->GetXaxis()->GetBinLowEdge(sigBins[0]);
  double xmax = mf->data->GetXaxis()->GetBinUpEdge(sigBins[1]);

  // TF1 integral is equivalent to TH1 integral with "width" option
  hSigPlusBkg = mf->data->Integral(sigBins[0], sigBins[1], "width");
  fSigPlusBkg = mf->fit->Integral(xmin, xmax);

  double background = 0;
  for (int i : bkgFits) {
    TF1* f = (TF1*)mf->fitParts[i - 1]->Clone(("f" + to_string(i)).c_str());
    f->SetRange(xmin, xmax);
    double bkg = f->Integral(xmin, xmax);
    inputs->printLog(TString::Format("f%d: %f", i, bkg).Data(), InputSettings::kDebugMax);
    background += bkg;
  }
  hBkg = background;
  fBkg = background;

  hSig = hSigPlusBkg - hBkg;
  fSig = fSigPlusBkg - fBkg;

  string s = TString::Format("SignalFinder::calcSigBkg: \nhSigPlusBkg = %f, hSig = %f, hBkg = %f,\nfSigPlusBkg = %f, fSig = %f, fBkg = %f", hSigPlusBkg, hSig, hBkg, fSigPlusBkg, fSig, fBkg).Data();
  inputs->printLog(s, InputSettings::kDebug);
}

void SignalFinder::setBkgFits() {
  bkgFits.clear();
  switch (inputs->fitType) {
    case InputSettings::kPol1GausGaus:
      bkgFits.push_back(3);
      break;
    case InputSettings::kPol1GausGausExp:
      bkgFits.push_back(4);
      break;
    case InputSettings::kPol1ExpGausExp:
      bkgFits.push_back(4);
      break;
    default:
      inputs->printLog("SignalFinder::setBkgFits: can't set bkg fits for fit type " + to_string(inputs->fitType), InputSettings::kErrors);
      return;
  }

  string s = "SignalFinder::setBkgFits: bkg fits for function " + inputs->fitName + ": ";
  for (auto i : bkgFits) {
    s += to_string(i) + ", ";
  }
  inputs->printLog(s, InputSettings::kDebug);
}

array<double, 2> SignalFinder::getSignalRegion(bool useSigma, bool set = true) {
  array<double, 2> signalRegion = {-1., -1.};
  if (!mf->data) {
    inputs->printLog("SignalFinder::getSignalRegion() no data hist loaded!", InputSettings::kErrors);
    return signalRegion;
  }
  if (!mf->fit) {
    inputs->printLog("SignalFinder::getSignalRegion() no fit loaded!", InputSettings::kErrors);
    return signalRegion;
  }

  switch (inputs->fitType) {
    case InputSettings::kPol1GausGaus:
      if (useSigma)
        signalRegion = GG_SigRegionFromSteps(inputs->nSigma);
      else {
        GG_loadSigFracHist();
        signalRegion = GG_SigRegionFromFrac(signalFraction);
      }
      break;
    case InputSettings::kPol1GausGausExp:
      if (useSigma)
        signalRegion = GGE_SigRegionFromSteps(inputs->nSigma);
      else {
        GGE_loadSigFracHist(true);
        GGE_loadSigFracHist(false);
        signalRegion = GGE_SigRegionFromFrac(signalFraction);
      }
      break;
    case InputSettings::kPol1ExpGausExp:
      if (useSigma)
        signalRegion = EGE_SigRegionFromSteps(inputs->nSigma);
      else {
        EGE_loadSigFracHist(true);
        EGE_loadSigFracHist(false);
        signalRegion = EGE_SigRegionFromFrac(signalFraction);
      }
      break;
    default:
      inputs->printLog("SignalFinder::getSignalRegion: fit type " + to_string(inputs->fitType) + " not supported", InputSettings::kErrors);
      return signalRegion;
  }
  if (set)
    inputs->setSignalRegion(signalRegion[0], signalRegion[1]);
  return signalRegion;
}

void SignalFinder::getSignalFraction(double n) {
  const bool left = true;
  switch (inputs->fitType) {
    case InputSettings::kPol1GausGaus:
      signalFraction = GG_SigFrac(n);
      break;
    case InputSettings::kPol1GausGausExp:
      signalFractionLeft = GGE_SigFrac(n, nullptr, left);
      signalFractionRight = GGE_SigFrac(n, nullptr, !left);
      break;
    case InputSettings::kPol1ExpGausExp:
      signalFractionLeft = EGE_SigFrac(n, nullptr, left);
      signalFractionRight = EGE_SigFrac(n, nullptr, !left);
      break;
    default:
      inputs->printLog("SignalFinder::getSignalFraction(): no sigfrac for fit " + inputs->fitName, InputSettings::kErrors);
  }
}

// For specific fit types
double SignalFinder::getGG_Sigma(TF1* f = nullptr) {
  if (!f) f = mf->fit;
  double ampNarrow = f->GetParameter(0);
  double sigmaNarrow = f->GetParameter(2);
  double ampWide = f->GetParameter(3);
  double sigmaWide = f->GetParameter(4);

  double Sigma = (ampNarrow * sigmaNarrow + ampWide * sigmaWide) / (ampNarrow + ampWide);
  return Sigma;
}

array<double, 2> SignalFinder::GG_SigRegionFromSteps(double n, TF1* f = nullptr) {
  if (!f) f = mf->fit;
  double mu = f->GetParameter(1);
  double Sigma = getGG_Sigma(f);
  return {mu - n * Sigma, mu + n * Sigma};
}

array<double, 2> SignalFinder::GG_SigRegionFromFrac(double s, TF1* f = nullptr) {
  if (!f) f = mf->fit;

  int bin = hSigFrac->FindFirstBinAbove(s);

  if (s > hSigFrac->GetBinContent(hSigFrac->GetNbinsX())) {
    string war = "SignalFinder::GG_SigRegionFromFrac() Warning: s = " + to_string(s) + " is larger than the maximum signal fraction in the histogram. Using the maximum value instead.";
    inputs->printLog(war, InputSettings::kWarnings);
    bin = hSigFrac->GetNbinsX();
  }

  double n = hSigFrac->GetXaxis()->GetBinCenter(bin);
  double Sigma = getGG_Sigma(f);
  double mu = f->GetParameter(1);

  inputs->nSigma = n;
  signalFraction = hSigFrac->GetBinContent(bin);
  return {mu - n * Sigma, mu + n * Sigma};
}

array<double, 2> SignalFinder::GGE_SigRegionFromSteps(double n, TF1* f = nullptr) {
  if (!f) f = mf->fit;

  double ampNarrow = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigmaNarrow = f->GetParameter(2);
  double lambda = f->GetParameter(3); // Crossover point
  double ampWide = f->GetParameter(4);
  double sigmaWide = f->GetParameter(5);

  TF1* g = new TF1("g", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4])", mf->inputs->fitmin, mf->inputs->fitmax);
  g->SetParameters(ampNarrow, mu, sigmaNarrow, ampWide, sigmaWide);
  double Sigma = getGG_Sigma(g);

  double tau = sigmaNarrow / lambda;

  double leftSide = mu - n * Sigma;

  double rightSide = mu + n * sigmaNarrow;
  if (n > lambda) rightSide = mu + lambda * sigmaNarrow + (n - lambda) * tau;

  return {leftSide, rightSide};
}

array<double, 2> SignalFinder::GGE_SigRegionFromFrac(double s, TF1* f = nullptr) {
  if (!f) f = mf->fit;

  int binL = hSigFracL->FindFirstBinAbove(s);
  int binR = hSigFracR->FindFirstBinAbove(s);

  if (s > hSigFracL->GetBinContent(hSigFracL->GetNbinsX())) {
    string war = "SignalFinder::GGE_SigRegionFromFrac() Warning: s = " + to_string(s) + " is larger than the maximum signal fraction in the histogram. Using the maximum value instead.";
    inputs->printLog(war, InputSettings::kWarnings);
    binL = hSigFracL->GetNbinsX();
  }
  if (s > hSigFracR->GetBinContent(hSigFracR->GetNbinsX())) {
    string war = "SignalFinder::GGE_SigRegionFromFrac() Warning: s = " + to_string(s) + " is larger than the maximum signal fraction in the histogram. Using the maximum value instead.";
    inputs->printLog(war, InputSettings::kWarnings);
    binR = hSigFracR->GetNbinsX();
  }

  double nL = hSigFracL->GetXaxis()->GetBinCenter(binL);
  double nR = hSigFracR->GetXaxis()->GetBinCenter(binR);

  double ampNarrow = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigmaNarrow = f->GetParameter(2);
  double lambda = f->GetParameter(3); // Crossover point
  double ampWide = f->GetParameter(4);
  double sigmaWide = f->GetParameter(5);
  double tau = sigmaNarrow / lambda;

  TF1* g = new TF1("g", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4])", mf->inputs->fitmin, mf->inputs->fitmax);
  g->SetParameters(ampNarrow, mu, sigmaNarrow, ampWide, sigmaWide);
  double Sigma = getGG_Sigma(g);

  inputs->nSigmaL = nL;
  inputs->nSigmaR = nR;
  signalFractionLeft = hSigFracL->GetBinContent(binL);
  signalFractionRight = hSigFracR->GetBinContent(binR);

  double leftSide = mu - nL * Sigma;
  double rightSide = mu + nR * sigmaNarrow;
  if (nR > lambda) rightSide = mu + lambda * sigmaNarrow + (nR - lambda) * tau;

  return {leftSide, rightSide};
}

array<double, 2> SignalFinder::EGE_SigRegionFromSteps(double n, TF1* f = nullptr) {
  if (!f) f = mf->fit;

  double amp = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  double lambdaL = f->GetParameter(3); // Crossover point left
  double lambdaR = f->GetParameter(4); // Crossover point right

  double tauL = sigma / lambdaL;
  double tauR = sigma / lambdaR;

  double leftSide = mu - n * sigma;
  if (n > lambdaL) leftSide = mu - lambdaL * sigma - (n - lambdaL) * tauL;

  double rightSide = mu + n * sigma;
  if (n > lambdaR) rightSide = mu + lambdaR * sigma + (n - lambdaR) * tauR;

  return {leftSide, rightSide};
}

array<double, 2> SignalFinder::EGE_SigRegionFromFrac(double s, TF1* f = nullptr) {
  if (!f) f = mf->fit;

  int binL = hSigFracL->FindFirstBinAbove(s);
  int binR = hSigFracR->FindFirstBinAbove(s);

  if (s > hSigFracL->GetBinContent(hSigFracL->GetNbinsX())) {
    string war = "SignalFinder::EGE_SigRegionFromFrac() Warning: s = " + to_string(s) + " is larger than the maximum signal fraction in the histogram. Using the maximum value instead.";
    inputs->printLog(war, InputSettings::kWarnings);
    binL = hSigFracL->GetNbinsX();
  }
  if (s > hSigFracR->GetBinContent(hSigFracR->GetNbinsX())) {
    string war = "SignalFinder::EGE_SigRegionFromFrac() Warning: s = " + to_string(s) + " is larger than the maximum signal fraction in the histogram. Using the maximum value instead.";
    inputs->printLog(war, InputSettings::kWarnings);
    binR = hSigFracR->GetNbinsX();
  }

  double nL = hSigFracL->GetXaxis()->GetBinCenter(binL);
  double nR = hSigFracR->GetXaxis()->GetBinCenter(binR);

  double mu = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  double lambdaL = f->GetParameter(3); // Crossover point left
  double lambdaR = f->GetParameter(4); // Crossover point right

  inputs->nSigmaL = nL;
  inputs->nSigmaR = nR;
  signalFractionLeft = hSigFracL->GetBinContent(binL);
  signalFractionRight = hSigFracR->GetBinContent(binR);

  double tauL = sigma / lambdaL;
  double tauR = sigma / lambdaR;

  double leftSide = mu - nL * sigma;
  if (nL > lambdaL) leftSide = mu - lambdaL * sigma - (nL - lambdaL) * tauL;

  double rightSide = mu + nR * sigma;
  if (nR > lambdaR) rightSide = mu + lambdaR * sigma + (nR - lambdaR) * tauR;

  return {leftSide, rightSide};
}

double SignalFinder::GG_SigFrac(double n, TF1* f = nullptr) {
  if (!f) f = mf->fit;

  double ampNarrow = f->GetParameter(0);
  double sigmaNarrow = f->GetParameter(2);
  double ampWide = f->GetParameter(3);
  double sigmaWide = f->GetParameter(4);

  inputs->printLog(TString::Format("GG_SigFrac: n = %f, \nampNarrow = %f, sigmaNarrow = %f, ampWide = %f, sigmaWide = %f", n, ampNarrow, sigmaNarrow, ampWide, sigmaWide).Data(), InputSettings::kDebug);

  double Sigma = getGG_Sigma(f);

  double x = n * Sigma / (sqrt(2) * sigmaNarrow);
  double y = n * Sigma / (sqrt(2) * sigmaWide);

  x = TMath::Erf(x);
  x *= ampNarrow * sigmaNarrow;

  y = TMath::Erf(y);
  y *= ampWide * sigmaWide;

  double s = (x + y) / (ampNarrow * sigmaNarrow + ampWide * sigmaWide);
  inputs->printLog(TString::Format("GG_SigFrac: s = %f", s).Data(), InputSettings::kDebug);
  return s;
}

double SignalFinder::GGE_SigFrac(double n, TF1* f = nullptr, bool leftSide = false) {
  if (!f) f = mf->fit;

  double ampNarrow = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigmaNarrow = f->GetParameter(2);
  double lambda = f->GetParameter(3); // Crossover point
  double ampWide = f->GetParameter(4);
  double sigmaWide = f->GetParameter(5);

  inputs->printLog(TString::Format("GGE_SigFrac: n = %f, \nampNarrow = %f, mu = %f, sigmaNarrow = %f, lambda = %f, \nampWide = %f, sigmaWide = %f", n, ampNarrow, mu, sigmaNarrow, lambda, ampWide, sigmaWide).Data(), InputSettings::kDebug);

  if (leftSide) {
    // Make double Gaussian from f
    TF1* g = new TF1("g", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4])", f->GetXmin(), f->GetXmax());
    g->SetParameters(ampNarrow, mu, sigmaNarrow, ampWide, sigmaWide);
    return GG_SigFrac(n, g);
  }

  double N = ampNarrow * sigmaNarrow * TMath::Sqrt(TMath::PiOver2());
  double W = ampWide * sigmaWide * TMath::Sqrt(TMath::PiOver2());
  double E = ampNarrow * sigmaNarrow * exp(-lambda*lambda / 2) / lambda;

  double s;
  if (n > lambda) {
    double x = W * TMath::Erf((lambda + n / lambda - 1) * sigmaNarrow / (sqrt(2) * sigmaWide));
    double y = N * TMath::Erf(lambda / sqrt(2));
    double z = E * (1 - exp(lambda - n));
    s = x + y + z;
  } else {
    TF1* g = new TF1("g", "[0]*TMath::Gaus(x,[1],[2]) + [3]*TMath::Gaus(x,[1],[4])", f->GetXmin(), f->GetXmax());
    g->SetParameters(ampNarrow, mu, sigmaNarrow, ampWide, sigmaWide);
    double Sigma = getGG_Sigma(g);

    double x = W * TMath::Erf((n * Sigma) / (sqrt(2) * sigmaWide));
    double y = N * TMath::Erf((n * Sigma) / (sqrt(2) * sigmaNarrow));
    s = x + y;
  }

  s /= (W + N * TMath::Erf(lambda / sqrt(2)) + E);
  inputs->printLog(TString::Format("GGE_SigFrac: s_R = %f", s).Data(), InputSettings::kDebug);
  return s;
}

double SignalFinder::EGE_SigFrac(double n, TF1* f = nullptr, bool leftSide = false) {
  if (!f) f = mf->fit;

  double amp = f->GetParameter(0);
  double mu = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  double lambdaL = f->GetParameter(3); // Crossover point left
  double lambdaR = f->GetParameter(4); // Crossover point right

  inputs->printLog(TString::Format("EGE_SigFrac: n = %f, \namp = %f, mu = %f, sigma = %f, lambdaL = %f, lambdaR = %f", n, amp, mu, sigma, lambdaL, lambdaR).Data(), InputSettings::kDebug);

  double lambda = (leftSide) ? lambdaL : lambdaR;

  double x = lambda * sqrt(TMath::PiOver2()) * TMath::Erf(lambda / TMath::Sqrt2());

  double y = exp(- lambda*lambda / 2);

  double s = x;
  if (n > lambda) {
    s += y * (1 - exp(lambda - n));
  }
  s /= (x + y);
  inputs->printLog(TString::Format("EGE_SigFrac: s_%s = %f", (leftSide) ? "L" : "R", s).Data(), InputSettings::kDebug);
  return s;
}

TH1* SignalFinder::GG_makeSigFracHist(TF1* f = nullptr) {
  if (!f) f = mf->fit;

  int nx = nBins;
  double xmin = xMin, xmax = xMax;
  TH1* hist = new TH1D("hist", "Signal Fraction;#it{n}#Sigma", nx, xmin, xmax);

  for (int i = 0; i < nx; i++) {
    double n = hist->GetXaxis()->GetBinCenter(i + 1);
    double s = GG_SigFrac(n, f);
    hist->SetBinContent(i + 1, s);
  }
  string histName = inputs->getSaveNameFromPt("signalFraction");
  hSigFrac = (TH1*)hist->Clone(histName.c_str());
  return hist;
}

TH1* SignalFinder::GG_loadSigFracHist() {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    string err = "SignalFinder::GG_loadSigFracHist() Error: could not open file " + inputs->inputFileName;
    inputs->printLog(err, InputSettings::kErrors);
    return nullptr;
  }

  string s = inputs->histName;
  if (s == "")
    s = inputs->getSaveNameFromPt("signalFraction");

  TH1* hist = (TH1*)file->Get(s.c_str());
  if (!hist) {
    string err = "SignalFinder::GG_loadSigFracHist() Error: could not find histogram " + s + " in file " + inputs->inputFileName;
    inputs->printLog(err, InputSettings::kErrors);
    return nullptr;
  }
  hSigFrac = (TH1*)hist->Clone();
  return hist;
}

TH1* SignalFinder::GGE_makeSigFracHist(TF1* f = nullptr, bool leftSide = false) {
  if (!f) f = mf->fit;

  int nx = nBins;
  double xmin = xMin, xmax = xMax;
  string histTitle = TString::Format("Signal Fraction, %s;#it{n}#Sigma", (leftSide) ? "left side" : "right side").Data();
  TH1* hist = new TH1D("hist", histTitle.c_str(), nx, xmin, xmax);

  for (int i = 0; i < nx; i++) {
    double n = hist->GetXaxis()->GetBinCenter(i + 1);
    double s = GGE_SigFrac(n, f, leftSide);
    hist->SetBinContent(i + 1, s);
  }
  string histName = TString::Format("signalFraction%s", (leftSide) ? "Left" : "Right").Data();
  histName = inputs->getSaveNameFromPt(histName.c_str());

  if (leftSide)
    hSigFracL = (TH1*)hist->Clone(histName.c_str());
  else
    hSigFracR = (TH1*)hist->Clone(histName.c_str());

  return hist;
}

TH1* SignalFinder::GGE_loadSigFracHist(bool leftSide = false) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    string err = "SignalFinder::GGE_loadSigFracHist() Error: could not open file " + inputs->inputFileName;
    inputs->printLog(err, InputSettings::kErrors);
    return nullptr;
  }

  string s = inputs->histName;
  if (s == "") {
    s = TString::Format("signalFraction%s", (leftSide) ? "Left" : "Right").Data();
    s = inputs->getSaveNameFromPt(s.c_str());
  }

  TH1* hist = (TH1*)file->Get(s.c_str());
  if (!hist) {
    string err = "SignalFinder::GGE_loadSigFracHist() Error: could not find histogram " + s + " in file " + inputs->inputFileName;
    inputs->printLog(err, InputSettings::kErrors);
    return nullptr;
  }

  if (leftSide)
    hSigFracL = (TH1*)hist->Clone();
  else
    hSigFracR = (TH1*)hist->Clone();

  return hist;
}

TH1* SignalFinder::EGE_makeSigFracHist(TF1* f = nullptr, bool leftSide = false) {
  if (!f) f = mf->fit;

  int nx = nBins;
  double xmin = xMin, xmax = xMax;
  string histTitle = TString::Format("Signal Fraction, %s;#it{n}#Sigma", (leftSide) ? "left side" : "right side").Data();
  string histName = TString::Format("signalFraction%s", (leftSide) ? "Left" : "Right").Data();
  histName = inputs->getSaveNameFromPt(histName.c_str());
  TH1* hist = new TH1D(histName.c_str(), histTitle.c_str(), nx, xmin, xmax);

  for (int i = 0; i < nx; i++) {
    double n = hist->GetXaxis()->GetBinCenter(i + 1);
    double s = EGE_SigFrac(n, f, leftSide);
    hist->SetBinContent(i + 1, s);
  }

  if (leftSide)
    hSigFracL = (TH1*)hist->Clone();
  else
    hSigFracR = (TH1*)hist->Clone();

  return hist;
}

TH1* SignalFinder::EGE_loadSigFracHist(bool leftSide = false) {
  TFile* file = TFile::Open(inputs->inputFileName.c_str(), "READ");
  if (!file) {
    string err = "SignalFinder::EGE_loadSigFracHist() Error: could not open file " + inputs->inputFileName;
    inputs->printLog(err, InputSettings::kErrors);
    return nullptr;
  }

  string s = inputs->histName;
  if (s == "") {
    s = TString::Format("signalFraction%s", (leftSide) ? "Left" : "Right").Data();
    s = inputs->getSaveNameFromPt(s.c_str());
  }

  TH1* hist = (TH1*)file->Get(s.c_str());
  if (!hist) {
    string err = "SignalFinder::EGE_loadSigFracHist() Error: could not find histogram " + s + " in file " + inputs->inputFileName;
    inputs->printLog(err, InputSettings::kErrors);
    return nullptr;
  }

  if (leftSide)
    hSigFracL = (TH1*)hist->Clone();
  else
    hSigFracR = (TH1*)hist->Clone();

  return hist;
}

TH1* SignalFinder::makeSigRegionHist() {
  array<int, 2> signalRegionBins = histutils::getProjectionBins(mf->data->GetXaxis(), inputs->signalRegionMin, inputs->signalRegionMax);
  return histutils::makeHistSubset(mf->data, signalRegionBins[0], signalRegionBins[1], this->inputs->getSaveNameFromPt("signalRegion"));
}

// -------------------------------------------------------------------------------------------------
//
// Struct for plotting
//
// -------------------------------------------------------------------------------------------------

struct FitPlotter {
  private:
  public:
    InputSettings* inputs;
    MassFitter* mf;

    TCanvas* canvas = nullptr;
    TH1F* frame = nullptr;
    TLegend* legend = nullptr;
    vector<TObject*> plottingObjects = {};

    double expectedSignalFraction = -1;
    double expectedSigFracLeft = -1;
    double expectedSigFracRight = -1;

    FitPlotter() { inputs = new InputSettings(); mf = new MassFitter(*inputs); }
    FitPlotter(InputSettings& x) { inputs = &x; mf = new MassFitter(x); }
    FitPlotter(MassFitter& x) { mf = &x; inputs = x.inputs; }

    void addLatex(double x, double y, string text, double textSize = 0.04);
    void createLegend(double x1, double x2, double y1, double y2, string title);
    void fillLegendWithFitParts();
    TH1* loadSignalRegionHist();
    void plotFitInfo();
    int plotFitInfo(int infoToPlot);
    void plotFitParts();
    void plotSignalRegion();

    void autoCanvas();
    void autoFrame();
};

void FitPlotter::addLatex(double x, double y, string text, double textSize) {
  TLatex* latex = plotutils::CreateLatex(x, y, text, textSize);
  this->plottingObjects.push_back(latex);
}

void FitPlotter::createLegend(double x1, double x2, double y1, double y2, string title = "") {
  this->legend = plotutils::CreateLegend(x1, x2, y1, y2, title, inputs->textSize);
}

void FitPlotter::fillLegendWithFitParts() {
  if (!this->legend) {
    this->inputs->printLog("FitPlotter::fillLegendWithFitParts() Error: Legend not set. Aborting", InputSettings::kErrors);
    return;
  }

  if (this->mf->fitParts.empty()) {
    this->inputs->printLog("FitPlotter::fillLegendWithFitParts() Error: fit parts not set. Aborting", InputSettings::kErrors);
    return;
  }

  legend->AddEntry(this->mf->fit, "Total Fit", "l");

  vector<string> partNames;
  switch (this->inputs->fitType) {
    case InputSettings::kPol1GausGaus:
      partNames = {"Signal (Narrow Gaussian)", "Signal (Wide Gaussian)", "Background"};
      break;
    case InputSettings::kPol1GausGausExp:
      partNames = {"Signal (Narrow Gaussian)", "Signal (Wide Gaussian)", "Signal (Exp tail)", "Background"};
      break;
    case InputSettings::kPol1ExpGausExp:
      partNames = {"Signal (Gaussian)", "Signal (Exp tail)", "Signal (Exp tail)", "Background"};
      break;
    default:
      this->inputs->printLog("FitPlotter::fillLegendWithFitParts() Error: fit type " + to_string(this->inputs->fitType) + " not supported. Aborting", InputSettings::kErrors);
      return;
  }

  for (size_t i = 0; i < this->mf->fitParts.size(); i++) {
    this->legend->AddEntry(mf->fitParts[i], partNames[i].c_str(), "l");
  }
}

TH1* FitPlotter::loadSignalRegionHist() {
  array<int, 2> signalRegionBins = histutils::getProjectionBins(this->mf->data->GetXaxis(), this->mf->inputs->signalRegionMin, this->mf->inputs->signalRegionMax);
  TH1* hist = histutils::makeHistSubset(this->mf->data, signalRegionBins[0], signalRegionBins[1], this->inputs->getSaveNameFromPt("signalRegion"));
  return hist;
}

void FitPlotter::plotFitInfo() {
  int maxInfos = 20; // To prevent infinite loop
  for (int i = 0; i < maxInfos; i++) {
    if (this->plotFitInfo(i) != 0)
      break;
  }
}

int FitPlotter::plotFitInfo(int infoToPlot) {
  inputs->printLog("FitPlotter::plotFitInfo(): plotting info " + to_string(infoToPlot) + " from histogram " + this->inputs->histName, InputSettings::kInfo);

  TFile* file = TFile::Open(this->inputs->inputFileName.c_str(), "READ");

  TH2* fitInfo = (TH2*)file->Get(this->inputs->histName.c_str());
  if (!fitInfo) {
    string s = "FitPlotter::plotFitInfo() Error: Could not find histogram " + this->inputs->histName + " in file " + this->inputs->inputFileName;
    this->inputs->printLog(s, InputSettings::kErrors);
    return 1;
  }

  if (infoToPlot < 0 || infoToPlot >= fitInfo->GetNbinsY()) {
    string s = "FitPlotter::plotFitInfo() Info: Invalid column " + to_string(infoToPlot) + " selected for plotting " + this->inputs->histName + "\nValid range is 0 to " + to_string(fitInfo->GetNbinsY() - 1);
    this->inputs->printLog(s, InputSettings::kInfo);
    return 2;
  }

  TH1* hist = (TH1*)fitInfo->ProjectionX("hist", infoToPlot+1, infoToPlot+1);
  plotutils::setStyle(hist, 0);
  if (!this->canvas) {
    string s = "FitPlotter::plotFitIfo(): Canvas not set, creating automatically";
    this->inputs->printLog(s, InputSettings::kInfo);
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetLogy(this->inputs->logplot);
    this->canvas = canvas;
  }

  this->canvas->cd();
  string histTitle;
  double yMinFrame, yMaxFrame;
  if (this->frame) {
    histTitle = this->frame->GetTitle();
    yMinFrame = this->frame->GetMinimum();
    yMaxFrame = this->frame->GetMaximum();
  } else {
    histTitle = TString::Format("%s, %s, ", histutils::getDataSet(this->inputs->train).c_str(), this->inputs->fitName.c_str()).Data();
    if (this->inputs->histName == "fitParams")
      histTitle += "par " + to_string(infoToPlot);
    if (this->inputs->histName == "fitResults")
      histTitle += fitInfo->GetYaxis()->GetBinLabel(infoToPlot + 1);

    yMinFrame = histutils::getLowerBound(hist, false);
    yMinFrame *= (yMinFrame < 0) ? 1.1 : 0.9;
    yMaxFrame = 1.1 * histutils::getUpperBound(hist, false);
  }

  hist->SetTitle(histTitle.c_str());
  hist->SetStats(0);
  hist->SetMinimum(yMinFrame);
  hist->SetMaximum(yMaxFrame);

  string outputFileName = this->inputs->outputFileName;
  if (outputFileName == "") {
    outputFileName = this->inputs->hadron + "_" + this->inputs->fitName + "_" + this->inputs->histName + to_string(infoToPlot) + ".pdf";
    string s = "FitPlotter::plotFitInfo(): outputFileName not set, using: " + outputFileName;
    this->inputs->printLog(s, InputSettings::kInfo);
  }
  hist->Draw();
  hist->Draw("same text45");
  canvas->SaveAs(outputFileName.c_str());
  return 0;
}

void FitPlotter::plotFitParts() {
  if (this->inputs->outputFileName == "") {
    this->inputs->outputFileName = this->inputs->getSaveNameFromPt(this->inputs->hadron + "_" + this->inputs->fitName, ".pdf");
    string s = "FitPlotter::plotFitParts(): outputFileName not set, using: " + this->inputs->outputFileName;
    this->inputs->printLog(s, InputSettings::kInfo);
  }

  if (!this->canvas) {
    string s = "FitPlotter::plotFitParts(): Canvas not set, creating automatically";
    this->inputs->printLog(s, InputSettings::kInfo);
    this->autoCanvas();
  }

  this->canvas->cd();
  if (!this->frame) {
    string s = "FitPlotter::plotFitParts(): Frame not set, creating automatically";
    this->inputs->printLog(s, InputSettings::kInfo);
    this->autoFrame();
  }

  this->frame->Draw();
  if (this->legend) this->legend->Draw("same");

  for (int i = 0; i < this->mf->fitParts.size(); i++) {
    TF1* f = this->mf->fitParts[i];
    plotutils::setStyle(f, i + 2);
    f->SetRange(this->mf->data->GetXaxis()->GetXmin(), this->mf->data->GetXaxis()->GetXmax());
    f->Draw("same c");
  }
  // Only do this if auto setup is requested
  plotutils::setStyle(this->mf->data, 0);
  plotutils::setStyle(this->mf->fit, 1);
  this->mf->fit->SetRange(this->mf->data->GetXaxis()->GetXmin(), this->mf->data->GetXaxis()->GetXmax());

  this->mf->data->Draw("same");
  this->mf->fit->Draw("same c");

  if (this->legend)
    legend->Draw();
  for (auto obj : this->plottingObjects) {
    obj->Draw("same");
  }
  this->canvas->SaveAs(this->inputs->outputFileName.c_str());
}

void FitPlotter::plotSignalRegion() {
  if (this->inputs->outputFileName == "") {
    this->inputs->outputFileName = this->inputs->getSaveNameFromPt(this->inputs->hadron + "_" + this->inputs->fitName, "_SR.pdf");
    string s = "FitPlotter::plotSignalRegion(): outputFileName not set, using: " + this->inputs->outputFileName;
    this->inputs->printLog(s, InputSettings::kInfo);
  }

  if (!this->canvas) {
    string s = "FitPlotter::plotSignalRegion(): Canvas not set, creating automatically";
    this->inputs->printLog(s, InputSettings::kInfo);
    this->autoCanvas();
  }

  this->canvas->cd();
  if (!this->frame) {
    string s = "FitPlotter::plotSignalRegion(): Frame not set, creating automatically";
    this->inputs->printLog(s, InputSettings::kInfo);
    this->autoFrame();
  }

  plotutils::setStyle(this->mf->data, 0);
  TH1* hSR = this->loadSignalRegionHist();
  hSR->SetFillColorAlpha(kGreen, 0.3);
  hSR->Print();

  this->frame->Draw();
  if (this->legend) {
    this->legend->AddEntry(this->mf->data, "Data", "p");
    this->legend->AddEntry(hSR, "Signal region", "f");
    string regionText = TString::Format("(%.3f, %.3f)", this->mf->inputs->signalRegionMin, this->mf->inputs->signalRegionMax).Data();
    legend->AddEntry((TObject*)0, regionText.c_str(), "");
    if (this->expectedSignalFraction > 0) {
      legend->AddEntry((TObject*)0, TString::Format("Signal: %.2f%%", 1e2*this->expectedSignalFraction).Data(), "");
    }
    this->legend->Draw("same");
  }

  hSR->Draw("bars same");
  this->mf->data->Draw("same");
  this->canvas->SaveAs(this->mf->inputs->outputFileName.c_str());
}

void FitPlotter::autoCanvas() {
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetLogy(this->inputs->logplot);
  this->canvas = canvas;
}

void FitPlotter::autoFrame() {
  inputs->printLog("FitPlotter::autoFrame(): Creating frame automatically", InputSettings::kInfo);

  double xMinFrame = this->mf->data->GetXaxis()->GetXmin();
  double xMaxFrame = this->mf->data->GetXaxis()->GetXmax();
  double yMinFrame = 0.;
  double yMaxFrame = 1.1 * histutils::getUpperBound(this->mf->data, false);
  string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", histutils::formatHadronDaughters(this->inputs->hadron).c_str()).Data();
  string yTitle = "#it{N}_{V0} (a. u.)";
  // string histTitle = TString::Format("%s, (%.1f < #it{p}_{T, V0} < %.1f)", histutils::getDataSet(this->inputs->train).c_str(), this->inputs->lowpt, this->inputs->highpt).Data();

  TH1F* frame = plotutils::DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->SetTitle("");
  frame->GetYaxis()->SetTitleOffset(1.);

  plotutils::SetPadMargins(0.15, 0.05, 0.15, 0.05);
  this->frame = frame;
}

// -------------------------------------------------------------------------------------------------
//
// Interface
//
// -------------------------------------------------------------------------------------------------

void printhelp() {
  cout << "Usage of plotMassFits\n"
       << "A) Fit invariant mass spectrum with MassFitter\n"
       << "B) Compile fit information with FitSummariser\n"
       << "C) Plot fit parts with FitPlotter\n"
       << "D) Plot fit information with FitPlotter\n"
       << endl;

  cout << "A) Fit mass:\n"
       << "1) Create InputSettings with hadron, fit type, input/output file name, pt range, signal region, fit region, and initial values for peak position and polynomial bkg\n"
       << "2) Do fitting with MassFitter\n" << endl;

  cout << "B) Summarise fit:\n"
       << "1) Create InputSettings with hadron, fit type, input/output file name\n"
       << "2) FitSummariser will summarise the TH1s outputted by MassFitter into a single TH2\n" << endl;

  cout << "C) Plot fit parts:\n"
       << "1) FitPlotter plots fit parts directly from MassFitter\n"
       << "1) Create InputSettings with hadron, fit type, input/output file name\n"
       << "2) FitPlotter will plot the fit parts from the histograms in the input file\n" << endl;

  cout << "D) Plot fit information:\n"
       << "1) Create InputSettings with hadron, fit type, input/output file name\n"
       << "2) FitPlotter plots fit information from the TH2s in the file\n"
       << endl;

  cout << "Hadrons implemented: K0S, Lambda, AntiLambda\n"
       << "Fit types implemented: \n"
       << "pol1BreitWigner, pol2BreitWigner\n"
       << "pol1GausExp, pol2GausExp\n"
       << "pol1GausGaus, pol2GausGaus\n"
       << "pol1GausGausExp, pol2GausGausExp\n"
       << "pol1GausXex, pol2GausXex\n"
       << "pol1Voigt, pol2Voigt\n" << endl;
}

void printFitExpression(string f) {
  InputSettings x;
  x.setFitType(f);
  cout << f << ": " << x.getFitExpression() << endl;
}

// Plot the fit parts by reading the file
void plotFitParts(InputSettings& x) {
  FitPlotter p(x);
  p.mf->loadSavedMassHist();
  p.mf->loadSavedFitFunction();
  p.mf->loadSavedFitParts();
  p.plotFitParts();
}

// Plot the fit parts by using the fitter
void plotFitParts(MassFitter& m) {
  FitPlotter p(m);
  p.inputs->outputFileName = p.inputs->getSaveNameFromPt(p.inputs->hadron + "_" + p.inputs->fitName, ".pdf");
  p.plotFitParts();
}

// Plot the fit information by reading the file
void plotFitInfo(InputSettings& x) {
  FitPlotter p(x);
  p.inputs->histName = "fitParams";
  p.plotFitInfo();
  p.inputs->histName = "fitResults";
  p.plotFitInfo();
}

// Plot the fit information by reading the file, automatic setup
void plotFitInfo() {
  InputSettings x;
  x.hadron = "K0S";
  x.train = 252064;
  x.setFitType("pol1ExpGausExp");
  x.setInputFileNameFromFit();
  x.nSigma = 3.;
  x.outputFileName = "";

  plotFitInfo(x);
}

// Summarise the fit information by reading the file
void summariseFitInfo(InputSettings& x) {
  FitSummariser f(x);
  f.summariseFitInfo("fitParams");
  f.summariseFitInfo("fitResults");
}

// Summarise the fit information by reading the file, automatic setup
void summariseFitInfo() {
  InputSettings x;
  x.hadron = "K0S";
  x.train = 252064;
  x.setFitType("pol1ExpGausExp");
  x.setInputFileNameFromFit();
  x.setPtBinEdgesFromHadron();
  x.nSigma = 3.;
  x.outputFileName = x.inputFileName;

  summariseFitInfo(x);
}

// -------------------------------------------------------------------------------------------------
//
// Fit the mass spectrum and plot the fit parts
//
// -------------------------------------------------------------------------------------------------

// This runs the entire workflow in order for given pt bins
void fitMassAndPlotPartsAllBins(InputSettings& x) {
  vector<vector<double>> ptBinEdges = x.ptBinEdges;

  // For each pt bin, do fit, and plot fit parts
  for (int iPt = 0; iPt < ptBinEdges.size(); iPt++) {
    x.setPt(ptBinEdges[iPt][0], ptBinEdges[iPt][1]);
    MassFitter m(x);
    m.doFitting();

    FitPlotter p(m);
    p.inputs->outputFileName = x.getSaveNameFromPt(x.hadron + "_" + x.fitName, ".pdf");
    p.plotFitParts();
    x.outputFileName = x.hadron + "_" + x.fitName + ".root";
  }

  // Compile the information into summary histograms
  x.setInputFileNameFromFit();
  x.outputFileName = x.inputFileName;
  FitSummariser f(x);
  f.summariseFitInfo("fitParams");
  f.summariseFitInfo("fitResults");

  // Plot the evolution of fit parameters and results with pt
  x.outputFileName = "";
  FitPlotter p(x);
  p.inputs->histName = "fitParams";
  p.plotFitInfo();
  p.inputs->histName = "fitResults";
  p.plotFitInfo();
}

// Setup for fitMassAndPlotPartsAllBins(x)
void fitMassAndPlotPartsAllBins() {
  gROOT->SetBatch();
  InputSettings x;
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();
  x.outputFileName = x.hadron + "_" + x.fitName + ".root";

  x.setPtBinEdgesFromHadron();
  fitMassAndPlotPartsAllBins(x);
}

// Perform fit in a single pt bin
void fitMassAndPlotPartsSingleBin() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(9., 10.);
  // x.setPt(30., 40.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;
  x.drawLegend = true;

  // x.setFitX(1.1, 1.13);
  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();
  // x.outputFileName = x.hadron + "_" + x.fitName + ".root";

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  histutils::printParLimits(m.fit);

  m.data->Fit(m.fit, "R");
  histutils::printParLimits(m.fit);
  // m.fixFitInPost(); // Applies `any post-fit fixes, like swapping gaussians
  m.loadFitParts();
  m.loadFitParams();
  // m.loadFitResults();
  m.loadResidualHist();
  m.loadPullHist();

  // m.writeOutputsToFile();

  FitPlotter p(m);
  if (p.inputs->drawLegend) {
    p.createLegend(0.55, 0.75, 0.7, 0.9);
    p.fillLegendWithFitParts();
  }
  if (p.inputs->drawLatex) {
    string lTitle = TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data();
    // Do Latex Stuff
  }

  // p.inputs->outputFileName = p.inputs->getSaveNameFromPt(p.inputs->hadron + "_" + p.inputs->fitName, ".pdf");
  p.plotFitParts();
}

// -------------------------------------------------------------------------------------------------
//
// Given a fit, calculate what the signal region should be
//
// -------------------------------------------------------------------------------------------------

// Creates the signal fraction histograms for the GG fit, as used in calcSignalRegion_GG
void saveSignalFractionHists_GG() {
  InputSettings x;
  x.hadron = "K0S";
  x.train = 252064;
  x.setFitType("pol1GausGaus");
  x.inputFileName = x.hadron + "_" + x.fitName + "/" + x.hadron + "_" + x.fitName + ".root";
  x.setPtBinEdgesFromHadron();

  for (int iPt = 0; iPt < x.ptBinEdges.size(); iPt++) {
    cout << "Processing pt bin " << iPt << ": (" << x.ptBinEdges[iPt][0] << ", " << x.ptBinEdges[iPt][1] << ")" << endl;
    x.setPt(x.ptBinEdges[iPt][0], x.ptBinEdges[iPt][1]);

    SignalFinder sf(x);
    sf.mf->loadSavedMassHist();
    sf.mf->loadSavedFitFunction();
    sf.nBins = 1e5;
    sf.GG_makeSigFracHist(sf.mf->fit);
    sf.inputs->outputFileName = sf.inputs->inputFileName;
    sf.inputs->writeOutputToFile(sf.hSigFrac);
  }
}

// Creates the signal fraction histograms for the EGE fit, as used in calcSignalRegion_EGE
void saveSignalFractionHists_EGE() {
  InputSettings x;
  x.hadron = "K0S";
  x.train = 252064;
  x.setFitType("pol1ExpGausExp");
  x.inputFileName = x.hadron + "_" + x.fitName + "/" + x.hadron + "_" + x.fitName + ".root";
  x.setPtBinEdgesFromHadron();
  x.verbosity = InputSettings::kInfo;

  const bool leftSide = true;

  for (int iPt = 0; iPt < x.ptBinEdges.size(); iPt++) {
    cout << "Processing pt bin " << iPt << ": (" << x.ptBinEdges[iPt][0] << ", " << x.ptBinEdges[iPt][1] << ")" << endl;
    x.setPt(x.ptBinEdges[iPt][0], x.ptBinEdges[iPt][1]);

    SignalFinder sf(x);
    sf.mf->loadSavedMassHist();
    sf.mf->loadSavedFitFunction();
    sf.nBins = 1e5;
    sf.inputs->outputFileName = sf.inputs->inputFileName;
    sf.EGE_makeSigFracHist(sf.mf->fit, leftSide);
    sf.EGE_makeSigFracHist(sf.mf->fit, !leftSide);

    sf.inputs->writeOutputToFile(sf.hSigFracL);
    sf.inputs->writeOutputToFile(sf.hSigFracR);
  }
}

// Creates the signal fraction histograms for the GGE fit, as used in calcSignalRegion_GGE
void saveSignalFractionHists_GGE() {
  InputSettings x;
  x.hadron = "K0S";
  x.train = 252064;
  x.setFitType("pol1GausGausExp");
  x.inputFileName = x.hadron + "_" + x.fitName + "/" + x.hadron + "_" + x.fitName + ".root";
  x.setPtBinEdgesFromHadron();
  x.verbosity = InputSettings::kInfo;

  const bool leftSide = true;

  for (int iPt = 0; iPt < x.ptBinEdges.size(); iPt++) {
    cout << "Processing pt bin " << iPt << ": (" << x.ptBinEdges[iPt][0] << ", " << x.ptBinEdges[iPt][1] << ")" << endl;
    x.setPt(x.ptBinEdges[iPt][0], x.ptBinEdges[iPt][1]);

    SignalFinder sf(x);
    sf.mf->loadSavedMassHist();
    sf.mf->loadSavedFitFunction();
    sf.nBins = 1e5;
    sf.inputs->outputFileName = sf.inputs->inputFileName;
    sf.GGE_makeSigFracHist(sf.mf->fit, leftSide);
    sf.GGE_makeSigFracHist(sf.mf->fit, !leftSide);

    sf.inputs->writeOutputToFile(sf.hSigFracL);
    sf.inputs->writeOutputToFile(sf.hSigFracR);
  }
}

// Calculates the signal region for the double gaussian, for 3 Sigma, s = 90%, 95%, and 99% signal fraction
// Uses saved signal fraction histogram
void calcSignalRegion_GG() {
  InputSettings x;
  x.hadron = "K0S";
  x.train = 252064;
  x.setFitType("pol1GausGaus");
  x.inputFileName = to_string(x.train) + "/MassFits/" + x.hadron + "_" + x.fitName + "_fixedMu/" + x.hadron + "_" + x.fitName + ".root";
  // x.inputFileName = x.hadron + "_" + x.fitName + "/" + x.hadron + "_" + x.fitName + ".root";

  // GG
  x.ptBinEdges = { {0., 1.}, {1., 2.}, {2., 3.}, {3., 4.}, {4., 5.}, {5., 10.}};

  for (int iPt = 0; iPt < x.ptBinEdges.size(); iPt++) {
    cout << "Processing pt bin " << iPt << ": (" << x.ptBinEdges[iPt][0] << ", " << x.ptBinEdges[iPt][1] << ")" << endl;
    x.setPt(x.ptBinEdges[iPt][0], x.ptBinEdges[iPt][1]);

    SignalFinder sf(x);
    sf.mf->loadSavedMassHist();
    sf.mf->loadSavedFitFunction();
    sf.GG_loadSigFracHist();

    // n = 3
    sf.inputs->nSigma = 3;
    sf.signalFraction = sf.GG_SigFrac(sf.inputs->nSigma, sf.mf->fit);
    array<double, 2> sigRegion = sf.GG_SigRegionFromSteps(sf.inputs->nSigma, sf.mf->fit);

    string coutput;
    // coutput = TString::Format("n = 3: \nSignal region: (%f, %f) \nSignal fraction: %f%%", sigRegion[0], sigRegion[1], sf.signalFraction * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(s = %.2f \\) \\%%}", sigRegion[0], sigRegion[1], sf.signalFraction * 100).Data();
    cout << coutput << endl;

    // TH1* hSR_n3 = sf.makeSigRegionHist();
    // hSR_n3->SetName("signalRegion_n3");

    // s = 0.90
    double desiredSignalFraction = 0.90;
    sf.signalFraction = desiredSignalFraction;
    sigRegion = sf.GG_SigRegionFromFrac(desiredSignalFraction, sf.mf->fit);

    // coutput = TString::Format("s = %.f%%: \nnSigma = %f \nSignal region: (%f, %f) \nSignal fraction: %f%%", desiredSignalFraction * 100, sf.inputs->nSigma, sigRegion[0], sigRegion[1], sf.signalFraction * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n = %.2f \\)}", sigRegion[0], sigRegion[1], sf.inputs->nSigma).Data();
    cout << coutput << endl;

    // TH1* hSR_s90 = sf.makeSigRegionHist();
    // hSR_s90->SetName("signalRegion_s90");

    // s = 0.95
    desiredSignalFraction = 0.95;
    sf.signalFraction = desiredSignalFraction;
    sigRegion = sf.GG_SigRegionFromFrac(desiredSignalFraction, sf.mf->fit);

    // coutput = TString::Format("s = %.f%%: \nnSigma = %f \nSignal region: (%f, %f) \nSignal fraction: %f%%", desiredSignalFraction * 100, sf.inputs->nSigma, sigRegion[0], sigRegion[1], sf.signalFraction * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n = %.2f \\)}", sigRegion[0], sigRegion[1], sf.inputs->nSigma).Data();
    cout << coutput << endl;

    // TH1* hSR_s95 = sf.makeSigRegionHist();
    // hSR_s95->SetName("signalRegion_s95");

    // s = 0.99
    desiredSignalFraction = 0.99;
    sf.signalFraction = desiredSignalFraction;
    sigRegion = sf.GG_SigRegionFromFrac(desiredSignalFraction, sf.mf->fit);

    // coutput = TString::Format("s = %.f%%: \nnSigma = %f \nSignal region: (%f, %f) \nSignal fraction: %f%%", desiredSignalFraction * 100, sf.inputs->nSigma, sigRegion[0], sigRegion[1], sf.signalFraction * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n = %.2f \\)}", sigRegion[0], sigRegion[1], sf.inputs->nSigma).Data();
    cout << coutput << endl;

    // TH1* hSR_s99 = sf.makeSigRegionHist();
    // hSR_s99->SetName("signalRegion_s99");
  }
}

// Calculates the signal region for the expgausexp, for 3 Sigma, s = 90%, 95%, and 99% signal fraction
// Uses saved signal fraction histogram
void calcSignalRegion_EGE() {
  InputSettings x;
  x.hadron = "K0S";
  x.train = 252064;
  x.setFitType("pol1ExpGausExp");
  x.inputFileName = x.hadron + "_" + x.fitName + "/" + x.hadron + "_" + x.fitName + ".root";

  // GGE & EGE
  x.ptBinEdges = { {10., 15.}, {15., 20.}, {20., 25}, {25., 30}, {30., 40.}};

  const bool leftSide = true;

  for (int iPt = 0; iPt < x.ptBinEdges.size(); iPt++) {
    cout << "Processing pt bin " << iPt << ": (" << x.ptBinEdges[iPt][0] << ", " << x.ptBinEdges[iPt][1] << ")" << endl;
    x.setPt(x.ptBinEdges[iPt][0], x.ptBinEdges[iPt][1]);

    SignalFinder sf(x);
    sf.mf->loadSavedMassHist();
    sf.mf->loadSavedFitFunction();
    sf.EGE_loadSigFracHist(leftSide);
    sf.EGE_loadSigFracHist(!leftSide);

    // n = 3
    sf.inputs->nSigma = 3;
    sf.signalFractionLeft = sf.EGE_SigFrac(sf.inputs->nSigma, sf.mf->fit, leftSide);
    sf.signalFractionRight = sf.EGE_SigFrac(sf.inputs->nSigma, sf.mf->fit, !leftSide);
    array<double, 2> sigRegion = sf.EGE_SigRegionFromSteps(sf.inputs->nSigma, sf.mf->fit);

    string coutput;
    coutput = TString::Format("n = 3: \nSignal region: (%f, %f) \nSignal fraction (left, right): %f%%, %f%%", sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(s_L = %.2f \\) \\%%, \\(s_R = %.2f \\) \\%%}", sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    cout << coutput << endl;

    // TH1* hSR_n3 = sf.makeSigRegionHist();
    // hSR_n3->SetName("signalRegion_n3");

    // s = 0.90
    double desiredSignalFraction = 0.90;
    sf.signalFractionLeft = desiredSignalFraction;
    sf.signalFractionRight = desiredSignalFraction;
    sigRegion = sf.EGE_SigRegionFromFrac(desiredSignalFraction, sf.mf->fit);
    // coutput = TString::Format("s = %.f%%: \nnSigma = (%f, %f) \nSignal region: (%f, %f) \nSignal fraction: (%f%%, %f%%)", desiredSignalFraction * 100, sf.inputs->nSigmaL, sf.inputs->nSigmaR, sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n_L = %.2f \\), \\(n_R = %.2f \\)}", sigRegion[0], sigRegion[1], sf.inputs->nSigmaL, sf.inputs->nSigmaR).Data();
    cout << coutput << endl;

    // TH1* hSR_s90 = sf.makeSigRegionHist();
    // hSR_s90->SetName("signalRegion_s90");

    // s = 0.95
    desiredSignalFraction = 0.95;
    sigRegion = sf.EGE_SigRegionFromFrac(desiredSignalFraction, sf.mf->fit);
    // coutput = TString::Format("s = %.f%%: \nnSigma = (%f, %f) \nSignal region: (%f, %f) \nSignal fraction: (%f%%, %f%%)", desiredSignalFraction * 100, sf.inputs->nSigmaL, sf.inputs->nSigmaR, sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n_L = %.2f \\), \\(n_R = %.2f \\)}", sigRegion[0], sigRegion[1], sf.inputs->nSigmaL, sf.inputs->nSigmaR).Data();
    cout << coutput << endl;

    // TH1* hSR_s95 = sf.makeSigRegionHist();
    // hSR_s95->SetName("signalRegion_s95");

    // s = 0.99
    desiredSignalFraction = 0.99;
    sigRegion = sf.EGE_SigRegionFromFrac(desiredSignalFraction, sf.mf->fit);
    // coutput = TString::Format("s = %.f%%: \nnSigma = (%f, %f) \nSignal region: (%f, %f) \nSignal fraction: (%f%%, %f%%)", desiredSignalFraction * 100, sf.inputs->nSigmaL, sf.inputs->nSigmaR, sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n_L = %.2f \\), \\(n_R = %.2f \\)}", sigRegion[0], sigRegion[1], sf.inputs->nSigmaL, sf.inputs->nSigmaR).Data();
    cout << coutput << endl;

    // TH1* hSR_s99 = sf.makeSigRegionHist();
    // hSR_s99->SetName("signalRegion_s99");
  }
}

// Calculates the signal region for the expgausexp, for 3 Sigma, s = 90%, 95%, and 99% signal fraction
// Uses saved signal fraction histogram
void calcSignalRegion_GGE() {
  InputSettings x; //x.verbosity = InputSettings::kDebug;
  x.hadron = "K0S";
  x.train = 252064;
  x.setFitType("pol1GausGausExp");
  x.inputFileName = to_string(x.train) + "/MassFits/" + x.hadron + "_" + x.fitName + "_fixedMu/" + x.hadron + "_" + x.fitName + ".root";
  // x.inputFileName = x.hadron + "_" + x.fitName + "/" + x.hadron + "_" + x.fitName + ".root";

  // GGE
  x.ptBinEdges = { {10., 15.}, {15., 20.}, {20., 25}, {25., 30}, {30., 40.}};

  const bool leftSide = true;

  for (int iPt = 0; iPt < x.ptBinEdges.size(); iPt++) {
    cout << "Processing pt bin " << iPt << ": (" << x.ptBinEdges[iPt][0] << ", " << x.ptBinEdges[iPt][1] << ")" << endl;
    x.setPt(x.ptBinEdges[iPt][0], x.ptBinEdges[iPt][1]);

    SignalFinder sf(x);
    sf.mf->loadSavedMassHist();
    sf.mf->loadSavedFitFunction();
    sf.GGE_loadSigFracHist(leftSide);
    sf.GGE_loadSigFracHist(!leftSide);

    // n = 3
    sf.inputs->nSigma = 3;
    sf.signalFractionLeft = sf.GGE_SigFrac(sf.inputs->nSigma, sf.mf->fit, leftSide);
    sf.signalFractionRight = sf.GGE_SigFrac(sf.inputs->nSigma, sf.mf->fit, !leftSide);
    array<double, 2> sigRegion = sf.GGE_SigRegionFromSteps(sf.inputs->nSigma, sf.mf->fit);

    string coutput;
    // coutput = TString::Format("n = 3: \nSignal region: (%f, %f) \nSignal fraction (left, right): %f%%, %f%%", sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(s_L = %.2f \\) \\%%, \\(s_R = %.2f \\) \\%%}", sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    cout << coutput << endl;

    // TH1* hSR_n3 = sf.makeSigRegionHist();
    // hSR_n3->SetName("signalRegion_n3");

    // s = 0.90
    double desiredSignalFraction = 0.90;
    sf.signalFractionLeft = desiredSignalFraction;
    sf.signalFractionRight = desiredSignalFraction;
    sigRegion = sf.GGE_SigRegionFromFrac(desiredSignalFraction, sf.mf->fit);
    // coutput = TString::Format("s = %.f%%: \nnSigma = (%f, %f) \nSignal region: (%f, %f) \nSignal fraction: (%f%%, %f%%)", desiredSignalFraction * 100, sf.inputs->nSigmaL, sf.inputs->nSigmaR, sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n_L = %.2f \\), \\(n_R = %.2f \\)}", sigRegion[0], sigRegion[1], sf.inputs->nSigmaL, sf.inputs->nSigmaR).Data();
    cout << coutput << endl;

    // TH1* hSR_s90 = sf.makeSigRegionHist();
    // hSR_s90->SetName("signalRegion_s90");

    // s = 0.95
    desiredSignalFraction = 0.95;
    sigRegion = sf.GGE_SigRegionFromFrac(desiredSignalFraction, sf.mf->fit);
    // coutput = TString::Format("s = %.f%%: \nnSigma = (%f, %f) \nSignal region: (%f, %f) \nSignal fraction: (%f%%, %f%%)", desiredSignalFraction * 100, sf.inputs->nSigmaL, sf.inputs->nSigmaR, sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n_L = %.2f \\), \\(n_R = %.2f \\)}", sigRegion[0], sigRegion[1], sf.inputs->nSigmaL, sf.inputs->nSigmaR).Data();
    cout << coutput << endl;

    // TH1* hSR_s95 = sf.makeSigRegionHist();
    // hSR_s95->SetName("signalRegion_s95");

    // s = 0.99
    desiredSignalFraction = 0.99;
    sigRegion = sf.GGE_SigRegionFromFrac(desiredSignalFraction, sf.mf->fit);
    // coutput = TString::Format("s = %.f%%: \nnSigma = (%f, %f) \nSignal region: (%f, %f) \nSignal fraction: (%f%%, %f%%)", desiredSignalFraction * 100, sf.inputs->nSigmaL, sf.inputs->nSigmaR, sigRegion[0], sigRegion[1], sf.signalFractionLeft * 100, sf.signalFractionRight * 100).Data();
    coutput = TString::Format("& \\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n_L = %.2f \\), \\(n_R = %.2f \\)}", sigRegion[0], sigRegion[1], sf.inputs->nSigmaL, sf.inputs->nSigmaR).Data();
    cout << coutput << endl;

    // TH1* hSR_s99 = sf.makeSigRegionHist();
    // hSR_s99->SetName("signalRegion_s99");
  }
}

// -------------------------------------------------------------------------------------------------
//
// Calculate the purity in the signal region, given a fit
//
// -------------------------------------------------------------------------------------------------

void calcPurity(double nSigma, double desiredSignalFraction, bool printLatex = false) {
  if (nSigma * desiredSignalFraction > 0) {
    cout << "Please set either nSigma or desiredSignalFraction, not both. Aborting" << endl;
    return;
  }

  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.hadron = "K0S";
  x.train = 252064;
  x.setFitType("pol1ExpGausExp");
  x.inputFileName = to_string(x.train) + "/MassFits/" + x.hadron + "_" + x.fitName + "_fixedMu/" + x.hadron + "_" + x.fitName + ".root";

  bool useSigma = false;
  const bool leftSide = true;
  SignalFinder sf(x);
  sf.setBkgFits();

  x.nSigma = nSigma;
  sf.signalFraction = desiredSignalFraction;
  if (x.nSigma > 0)
    useSigma = true;

  if (x.fitType == InputSettings::kPol1GausGaus)
    x.ptBinEdges = { {0., 1.}, {1., 2.}, {2., 3.}, {3., 4.}, {4., 5.}, {5., 10.}};
  else
    x.ptBinEdges = { {10., 15.}, {15., 20.}, {20., 25.}, {25., 30.}, {30., 40.}};

  // In a given pt bin:
  for (int iPt = 0; iPt < x.ptBinEdges.size(); iPt++) {
    x.setPt(x.ptBinEdges[iPt][0], x.ptBinEdges[iPt][1]);
    sf.mf->loadSavedMassHist();
    sf.mf->loadSavedFitFunction();
    sf.getSignalRegion(useSigma);
    sf.calcSigBkg();
    if (useSigma)
      sf.getSignalFraction(x.nSigma);

    double hPurity = sf.hSig / sf.hSigPlusBkg;
    double fPurity = sf.fSig / sf.fSigPlusBkg;

    x.printLog(TString::Format("Pt: %.f - %.f", x.lowpt, x.highpt).Data(), InputSettings::kInfo);
    string s;

    switch (x.fitType) {
      case InputSettings::kPol1GausGaus:
        if (printLatex && useSigma)
          s = TString::Format("\\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(s = %.2f \\) \\%%, \\newline \\(p = %.2f \\) \\%%}", x.signalRegionMin, x.signalRegionMax, sf.signalFraction * 100., hPurity * 100.).Data();
        if (printLatex && !useSigma)
          s = TString::Format("\\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n = %.2f \\), \\newline \\(p = %.2f \\) \\%%}", x.signalRegionMin, x.signalRegionMax, x.nSigma, hPurity * 100.).Data();

        if (!printLatex && useSigma)
          s = TString::Format("Signal region: %f, %f \nExpected signal fraction: %f \nPurity: %f (h), %f (f) \nh-f: %g, f/h: %g, (h-f)/h: %g", x.signalRegionMin, x.signalRegionMax, sf.signalFraction, hPurity, fPurity, fPurity - hPurity, hPurity / fPurity, (hPurity - fPurity) / hPurity).Data();

        if (!printLatex && !useSigma)
          s = TString::Format("Signal region: %f, %f \nExpected n: %f \nPurity: %f (h), %f (f) \nh-f: %g, f/h: %g, (h-f)/h: %g", x.signalRegionMin, x.signalRegionMax, x.nSigma, hPurity, fPurity, fPurity - hPurity, hPurity / fPurity, (hPurity - fPurity) / hPurity).Data();
        break;
      case InputSettings::kPol1ExpGausExp: // Same as kPol1GausGausExp
      case InputSettings::kPol1GausGausExp:
        if (printLatex && useSigma)
          s = TString::Format("\\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(s_L = %.2f \\) \\%%, \\newline \\(s_R = %.2f \\) \\%%, \\newline \\(p = %.2f \\) \\%%}", x.signalRegionMin, x.signalRegionMax, sf.signalFractionLeft * 100., sf.signalFractionRight * 100., hPurity * 100.).Data();
        if (printLatex && !useSigma)
          s = TString::Format("\\RaggedRight{\\( (%.3f, %.3f) \\) \\newline \\(n_L = %.2f \\), \\newline \\(n_R = %.2f \\), \\newline \\(p = %.2f \\) \\%%}", x.signalRegionMin, x.signalRegionMax, x.nSigmaL, x.nSigmaR, hPurity * 100.).Data();

        if (!printLatex && useSigma)
          s = TString::Format("Signal region: %f, %f \nExpected signal fraction: %f (L), %f (R) \nPurity: %f (h), %f (f) \nh-f: %g, f/h: %g, (h-f)/h: %g", x.signalRegionMin, x.signalRegionMax, sf.signalFractionLeft, sf.signalFractionRight, hPurity, fPurity, fPurity - hPurity, hPurity / fPurity, (hPurity - fPurity) / hPurity).Data();

        if (!printLatex && !useSigma)
          s = TString::Format("Signal region: %f, %f \nExpected n: %f \nPurity: %f (h), %f (f) \nh-f: %g, f/h: %g, (h-f)/h: %g", x.signalRegionMin, x.signalRegionMax, x.nSigma, hPurity, fPurity, fPurity - hPurity, hPurity / fPurity, (hPurity - fPurity) / hPurity).Data();
        break;
      default:
        s = TString::Format("Unknown fit type: %s", x.fitName.c_str()).Data();
    }

    x.printLog(s, InputSettings::kInfo);
  }
}

#endif
