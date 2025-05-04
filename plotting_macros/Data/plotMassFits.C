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

#ifndef __PLOTMASSFITS_H__
#define __PLOTMASSFITS_H__

// Return hist scale in a range of the histogram
double getHistScaleInRange(TH1* h, int minBin, int maxBin, bool doError, bool doMin) {
  double scale = doMin ? 1e12 : -1e12;
  for (int i = minBin; i <= maxBin; i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    double s = bc;
    if (doError)
      doMin ? s -= be : s += be;

    if (doMin)
      scale = min(scale, s);
    else
      scale = max(scale, s);
  }
  return scale;
}

int getExtremeBinInRange(TH1* h, int minBin, int maxBin, bool doError, bool doMin) {
  int BIN = -1;
  double scale = doMin ? 1e12 : -1e12;
  for (int i = minBin; i <= maxBin; i++) {
    double bc = h->GetBinContent(i);
    double be = h->GetBinError(i);
    double s = bc;
    if (doError)
      doMin ? s -= be : s += be;

    if (doMin && s < scale) {
      scale = s;
      BIN = i;
    }
    else if (s > scale) {
      scale = s;
      BIN = i;
    }
  }
  return BIN;
}

double chisqInRange(TH1* hist, TF1* fit, int minBin, int maxBin) {
  double chisq = 0;
  for (int i = minBin; i <= maxBin; i++) {
    double bc = hist->GetBinContent(i);
    double be = hist->GetBinError(i);
    double fx = fit->Eval(hist->GetBinCenter(i));
    chisq += (fx - bc) * (fx - bc) / (be * be);
  }
  return chisq;
}

double chisqInRange(TH1* hist, TF1* fit, double xmin, double xmax) {
  array<int, 2> bins = getProjectionBins(hist->GetXaxis(), xmin, xmax);
  return chisqInRange(hist, fit, bins[0], bins[1]);
}

double ndfInRange(TF1* fit, int minBin, int maxBin) {
  return (1 + maxBin - minBin - fit->GetNpar());
}

double ndfInRange(TH1* data, TF1* fit, double xmin, double xmax) {
  array<int, 2> bins = getProjectionBins(data->GetXaxis(), xmin, xmax);
  return ndfInRange(fit, bins[0], bins[1]);
}

double roundToNextPowerOfTen(double x) {
  return std::pow(10., std::ceil(std::log10(x)));
}
double roundToPrevPowerOfTen(double x) {
  return std::pow(10., std::floor(std::log10(x)));
}

// Prints the parameter limits of a function
void printParLimits(TF1* f) {
  for (int i = 0; i < f->GetNpar(); i++) {
    double min, max;
    f->GetParLimits(i, min, max);
    cout << f->GetName() << " (" << i << " = " << f->GetParName(i) << ") " << f->GetParameter(i) << " (" << min << ", " << max << ")" << endl;
  }
  cout << endl;
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
    string hadron = "";
    string fitName = "";
    string histName = "";
    string inputFileName = "";
    string outputFileName = "";
    double ptmin = -1e3, ptmax = -1e3, lowpt = -1e3, highpt = -1e3;
    double fitmin = -1e3, fitmax = -1e3; // Fit range
    bool logplot = false;
    bool normaliseData = false;
    bool fixMu = false;
    vector<vector<double>> ptBinEdges = {};

    double massWindowMin = -1., massWindowMax = -1.;
    double polInitx0 = -1., polInitx1 = -1., polInitx2 = -1.;
    double sigalRegionMin = -1., sigalRegionMax = -1.;
    double nSigma = -1.;

    enum FitType {kPol1BreitWigner, kPol2BreitWigner, kPol1GausExp, kPol2GausExp, kPol1GausGaus, kPol2GausGaus, kPol1GausGausExp, kPol2GausGausExp, kPol1GausGausXex, kPol2GausGausXex, kPol1GausXex, kPol2GausXex, kPol1Voigt, kPol2Voigt};

    double getMass();
    string getFitExpression();
    string getSaveNameFromPt(string prefix, string suffix);
    int inputIssue(string home, string obj);
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

string InputSettings::setFitName() {
  switch (this->fitType) {
    case kPol1BreitWigner:
      this->fitName = "pol1BreitWigner";
      break;
    case kPol2BreitWigner:
      this->fitName = "pol2BreitWigner";
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
      cout << "InputSettings::setFitName() Error: requested unknown function" << endl;
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
    cout << "InputSettings::getFitType() Error: requested unknown function " << this->fitName << endl;
  }

  return this->fitType;
}

int InputSettings::setFitType(string fit) {
  this->fitName = fit;
  return setFitType();
}

double InputSettings::getMass() {
  if (this->hadron == "K0S")
    return MassK0S;
  if (this->hadron == "Lambda" || this->hadron == "AntiLambda")
    return MassLambda;

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
    case kPol1GausExp:
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - ([1]+[2]/[3])/2.)/([2]/[3])) * (x > [1]+[2]*[3]) + [4]+[5]*x";
      break;
    case kPol2GausExp:
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - ([1]+[2]/[3])/2.)/([2]/[3])) * (x > [1]+[2]*[3]) + [4]+[5]*x+[6]*x*x";
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
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - ([1]+[2]/[3])/2.)/([2]/[3])) * (x > [1]+[2]*[3]) + [4]*TMath::Gaus(x,[1],[5]) + [6]+[7]*x";
      break;
    case kPol2GausGausExp:
      s = "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3])) + [0]*TMath::Exp(-(x - ([1]+[2]/[3])/2.)/([2]/[3])) * (x > [1]+[2]*[3]) + [4]*TMath::Gaus(x,[1],[5]) + [6]+[7]*x+[8]*x*x";
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
      cout << "InputSettings::getFitExpression() Error: requested unknown function " << this->fitName << endl;
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
    cout << "InputSettings::inputIssue() Error: requested unknown object " << obj << endl;
    return 1;
  }

  cout << "InputSettings::" << home << "() Error: " << obj << " not initialised" << endl;
  return 2;
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
    cout << "InputSettings::setFitX() Error: fitmin > fitmax" << endl;
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
    cout << "InputSettings::setMassWindow() Error: massWindowMin > massWindowMax" << endl;
    return;
  }
  this->massWindowMin = a;
  this->massWindowMax = b;
}

void InputSettings::setMassWindowDiff(double a, double b) {
  if (a > b) {
    cout << "InputSettings::setMassWindowDiff() Error: massWindowMin > massWindowMax" << endl;
    return;
  }
  double mass = this->getMass();
  this->setMassWindow(mass - a, mass + b);
}

void InputSettings::setMassWindowDiffFromHadron() {
  if (this->hadron == "K0S") {
    this->setMassWindowDiff(2e-2, 2e-2);
  } else {
    this->setMassWindowDiff(-1., -1.);
  }
}

void InputSettings::setPolInitX(double x0, double x1 = -1, double x2 = -2) {
  if (x0 > x1 || x0 > x2 || x1 > x2) {
    cout << "InputSettings::setPolInitX() Error: should have polInitx0 < polInitx1 < polInitx2" << endl;
    return;
  }
  this->polInitx0 = x0;
  this->polInitx1 = x1;
  this->polInitx2 = x2;
}

void InputSettings::setPolInitXFromHadron() {
  if (this->hadron == "K0S") {
    this->setPolInitX(0.45, 0.55, 0.57);
  } else {
    this->setPolInitX(-1., -1., -1.);
  }
}

void InputSettings::setPt(double a, double b) {
  if (a > b) {
    cout << "InputSettings::setPt() Error: ptmin > ptmax" << endl;
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
    cout << "InputSettings::setSignalRegion() Error: sigalRegionMin > sigalRegionMax" << endl;
    return;
  }
  this->sigalRegionMin = a;
  this->sigalRegionMax = b;
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
    InputSettings inputs;
    TH1* data = nullptr;
    TF1* fit = nullptr;
    vector<TF1*> fitParts = {};
    TH1* fitResults = nullptr;
    TH1* fitParams = nullptr;
    TH1* pull = nullptr;
    TH1* residual = nullptr;

    double signal = -1.;
    double background = -1.;
    double signalPlusBackground = -1.;

    MassFitter() { InputSettings(); }
    MassFitter(InputSettings& x) { this->inputs = x; }

    void calcSigBkg();
    void doFitting(); // Full execution of workflow
    void fixFitInPost();
    TF1* loadFitFunction();
    TH1* loadMassHist();
    void setFitInitialValues();
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

  int iBkg = 0;
  if (this->inputs.fitType == InputSettings::kPol1BreitWigner || this->inputs.fitType == InputSettings::kPol2BreitWigner) {
    iBkg = 2;
  } else if (this->inputs.fitType == InputSettings::kPol1GausExp || this->inputs.fitType == InputSettings::kPol2GausExp){
    iBkg = 3;
  } else if (this->inputs.fitType == InputSettings::kPol1GausGaus || this->inputs.fitType == InputSettings::kPol2GausGaus){
    iBkg = 3;
  } else if (this->inputs.fitType == InputSettings::kPol1GausGausExp || this->inputs.fitType == InputSettings::kPol2GausGausExp){
    iBkg = 4;
  } else if (this->inputs.fitType == InputSettings::kPol1GausGausXex || this->inputs.fitType == InputSettings::kPol2GausGausXex){
    iBkg = 4;
  } else if (this->inputs.fitType == InputSettings::kPol1GausXex || this->inputs.fitType == InputSettings::kPol2GausXex){
    iBkg = 3;
  } else if (this->inputs.fitType == InputSettings::kPol1Voigt || this->inputs.fitType == InputSettings::kPol2Voigt){
    iBkg = 2;
  } else {
    cout << "MassFitter::calcSigBkg() Error: do not know how to set parameters for this function" << endl;
    return;
  }

  if (this->inputs.sigalRegionMin < 0 || this->inputs.sigalRegionMax < 0) {
    cout << "MassFitter::calcSigBkg() Automatically setting signal region to mu +/- n sigma" << endl;
    this->setSignalRegionFromSigma();
  }

  array<int, 2> sigBins = getProjectionBins(this->data->GetXaxis(), this->inputs.sigalRegionMin, this->inputs.sigalRegionMax);
  double xmin = this->data->GetXaxis()->GetBinLowEdge(sigBins[0]);
  double xmax = this->data->GetXaxis()->GetBinUpEdge(sigBins[1]);

  TF1* bkgFit = this->fitParts[iBkg - 1];
  bkgFit->SetRange(xmin, xmax);
  this->background = bkgFit->Integral(xmin, xmax);
  this->signalPlusBackground = this->data->Integral(sigBins[0], sigBins[1], "width");
  this->signal = this->signalPlusBackground - this->background;
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
  if (this->inputs.hadron == "K0S" && this->inputs.fitType == InputSettings::kPol1GausGaus) {
    bool gausflip = false;
    if (this->inputs.highpt < 4 + 1e-3)
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
  if (this->inputs.hadron == "K0S" && this->inputs.fitType == InputSettings::kPol2GausGaus) {
    bool gausflip = false;
    if (this->inputs.highpt > 1. + 1e-3 && this->inputs.highpt < 30. + 1e-3)
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
  string saveName = this->inputs.getSaveNameFromPt("fit");
  string expression = this->inputs.getFitExpression();
  TF1* f = new TF1(saveName.c_str(), expression.c_str(), this->inputs.fitmin, this->inputs.fitmax);
  this->fit = (TF1*)f->Clone();
  return f;
}

TH1* MassFitter::loadMassHist() {
  TFile* f = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!f) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return nullptr;
  }

  if (this->inputs.histName == "") {
    this->setHistNameFromTrain();
    cout << "MassFitter::loadMassHist(): histName not set, setting to: " << this->inputs.histName << endl;
  }
  THnSparse* hist = (THnSparse*)f->Get(this->inputs.histName.c_str());
  if (!hist) {
    cout << "Could not find histogram " << this->inputs.histName << " in file " << this->inputs.inputFileName << endl;
    return nullptr;
  }

  int projectionAxis = 1;
  if (this->inputs.hadron == "Lambda")
    projectionAxis = 2;
  if (this->inputs.hadron == "AntiLambda")
    projectionAxis = 3;

  array<int, 2> bins = getProjectionBins(hist->GetAxis(0), this->inputs.ptmin, this->inputs.ptmax);
  int minBin = bins[0], maxBin = bins[1];
  this->inputs.lowpt = hist->GetAxis(0)->GetBinLowEdge(minBin);
  this->inputs.highpt = hist->GetAxis(0)->GetBinUpEdge(maxBin);

  hist->GetAxis(0)->SetRange(minBin, maxBin);
  TH1* h = hist->Projection(projectionAxis);

  if (this->inputs.rebinNumber < 0)
    this->setRebinNumberFromHadronAndPt();

  if (this->inputs.rebinNumber > 1)
    h->Rebin(this->inputs.rebinNumber);

  if (this->inputs.normaliseData) {
    bins = getProjectionBins(h->GetXaxis(), this->inputs.massWindowMin, this->inputs.massWindowMax);
    minBin = bins[0], maxBin = bins[1];
    h->Scale(1./getHistScaleInRange(h, minBin, maxBin, false, false));
  }

  string hName = this->inputs.getSaveNameFromPt("data");
  this->data = (TH1*)h->Clone(hName.c_str());
  return h;
}

void MassFitter::setFitInitialValues() {
  array<int, 2> massBins = getProjectionBins(this->data->GetXaxis(), this->inputs.massWindowMin, this->inputs.massWindowMax);
  int minBin = massBins[0], maxBin = massBins[1];
  int extremeBin = getExtremeBinInRange(this->data, minBin, maxBin, false, false);

  int b0 = this->data->GetXaxis()->FindBin(this->inputs.polInitx0);
  int b1 = this->data->GetXaxis()->FindBin(this->inputs.polInitx1);
  int b2 = this->data->GetXaxis()->FindBin(this->inputs.polInitx2);
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

  if (this->inputs.hadron == "K0S") {
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

  if (this->inputs.fitType % 2 == 0) {
    b = (y1-y0)/(x1-x0);
    a = y0 - b*x0;
  } else {
    c = ( (y2-y1)/(x2-x1) - (y1-y0)/(x1-x0) ) / (x2-x0);
    b = (y1-y0)/(x1-x0) - c*(x1+x0);
    a = y0 - b*x0 - c*x0*x0;
  }

  switch (this->inputs.fitType) {
    case InputSettings::kPol1BreitWigner:
      this->fit->SetParameters(A, mu, Gamma, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0, 1);
      this->fit->SetParLimits(4, -1, 1);
      break;
    case InputSettings::kPol2BreitWigner:
      this->fit->SetParameters(A, mu, Gamma, a, b, c);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0, 1);
      this->fit->SetParLimits(4, -1, 1);
      // this->fit->SetParLimits(5, -1, 1);
      break;
    case InputSettings::kPol1GausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 10.);
      this->fit->SetParLimits(4, 0, 1);
      this->fit->SetParLimits(5, -1, 1);
      break;
    case InputSettings::kPol2GausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, a, b, c);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 10.);
      this->fit->SetParLimits(4, 0, 1);
      this->fit->SetParLimits(5, -1, 1);
      break;
    case InputSettings::kPol1GausGaus:
      this->fit->SetParameters(A, mu, sigma, B, rho, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
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
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 1.);
      this->fit->SetParLimits(4, 0., 1.);
      this->fit->SetParLimits(5, 0., 1.);
      this->fit->SetParLimits(6, -1., 1.);
      break;
    case InputSettings::kPol1GausGausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, B, rho, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 10.);
      this->fit->SetParLimits(4, 0., 1.);
      this->fit->SetParLimits(5, 0., 1.);
      this->fit->SetParLimits(6, 0., 1.);
      this->fit->SetParLimits(7, -1., 1.);
      break;
    case InputSettings::kPol2GausGausExp:
      this->fit->SetParameters(A, mu, sigma, lambda, B, rho, a, b, c);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 10.);
      this->fit->SetParLimits(4, 0., 1.);
      this->fit->SetParLimits(5, 0., 1.);
      this->fit->SetParLimits(6, 0., 1.);
      this->fit->SetParLimits(7, -1., 1.);
      break;
    case InputSettings::kPol1GausGausXex:
      this->fit->SetParameters(A, mu, sigma, B, rho, C, nu, tau, a, b);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 1.); // B
      this->fit->SetParLimits(4, 2e-2, 1.); // rho
      this->fit->SetParLimits(5, 0., 1.); // C
      this->fit->SetParLimits(6, 0.5, 0.6); // nu
      // this->fit->SetParLimits(7, -1., 1.); // tau
      this->fit->SetParLimits(8, 0., 1.); // a
      this->fit->SetParLimits(9, -1., 1.); // b
      break;
    case InputSettings::kPol2GausGausXex:
      this->fit->SetParameters(A, mu, sigma, B, rho, C, nu, tau, a, b);
      this->fit->SetParLimits(0, 0., 1.); // A
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      this->fit->SetParLimits(3, 0., 1.); // B
      this->fit->SetParLimits(4, 0., 1.); // rho
      // this->fit->SetParLimits(5, 0., 1.); // C
      this->fit->SetParLimits(6, 0.5, 0.6); // nu
      // this->fit->SetParLimits(7, -1., 1.); // tau
      this->fit->SetParLimits(8, 0., 1.); // a
      this->fit->SetParLimits(9, -1., 1.); // b
      break;
    case InputSettings::kPol1GausXex:
      this->fit->SetParameters(A, mu, sigma, C, nu, tau, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.); // sigma
      // this->fit->SetParLimits(3, 0., 1.); // C
      this->fit->SetParLimits(4, 0.5, 0.6); // nu
      // this->fit->SetParLimits(5, 0., 10.); // tau
      this->fit->SetParLimits(6, 0., 1.); // a
      this->fit->SetParLimits(7, -1., 1.); // b
      break;
    case InputSettings::kPol2GausXex:
      this->fit->SetParameters(A, mu, sigma, C, nu, tau, a, b, c);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      // this->fit->SetParLimits(3, 0., 1.);
      this->fit->SetParLimits(4, 0.5, 0.6);
      // this->fit->SetParLimits(5, 0., 10.);
      this->fit->SetParLimits(6, 0., 1.);
      this->fit->SetParLimits(7, -1., 1.);
      break;
    case InputSettings::kPol1Voigt:
      this->fit->SetParameters(A, mu, sigma, Gamma, a, b);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 1.);
      this->fit->SetParLimits(3, 0., 1.);
      this->fit->SetParLimits(4, 0., 1.);
      this->fit->SetParLimits(5, -1., 1.);
      break;
    case InputSettings::kPol2Voigt:
      this->fit->SetParameters(A, mu, sigma, Gamma, a, b, c);
      this->fit->SetParLimits(0, 0., 1.);
      if (this->inputs.fixMu) {
        this->fit->SetParLimits(1, this->data->GetXaxis()->GetBinLowEdge(extremeBin), this->data->GetXaxis()->GetBinUpEdge(extremeBin));
      } else {
        this->fit->SetParLimits(1, this->inputs.massWindowMin, this->inputs.massWindowMax);
      }
      this->fit->SetParLimits(2, 0., 0.005);
      this->fit->SetParLimits(3, 0., 0.005);
      this->fit->SetParLimits(4, 0., 1.);
      this->fit->SetParLimits(5, -1., 1.);
      break;
    default:
      cout << "MassFitter::setFitInitialValues() Error: do not have initial values for function " << this->inputs.fitName << endl;
      return;
  }
}

vector<TF1*> MassFitter::loadFitParts() {
  if (!this->fit)
    return {};

  string fName = this->fit->GetName();
  vector<TF1*> parts;

  if (this->inputs.fitType == InputSettings::kPol1BreitWigner || this->inputs.fitType == InputSettings::kPol2BreitWigner) {
    double A, mu, Gamma;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    Gamma = this->fit->GetParameter(2);
    a = this->fit->GetParameter(3);
    b = this->fit->GetParameter(4);
    if (this->inputs.fitType == InputSettings::kPol2BreitWigner) c = this->fit->GetParameter(5);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "breitwigner(x,[0],[1],[2])", this->inputs.fitmin, this->inputs.fitmax);
    f1->SetParameters(A, mu, Gamma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    if (this->inputs.fitType == InputSettings::kPol1BreitWigner) {
      TF1* f2 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs.fitmin, this->inputs.fitmax);
      f2->SetParameters(a, b);
      parts.push_back(f2);
    } else {
      TF1* f2 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs.fitmin, this->inputs.fitmax);
      f2->SetParameters(a, b, c);
      parts.push_back(f2);
    }
  } else if (this->inputs.fitType == InputSettings::kPol1GausExp || this->inputs.fitType == InputSettings::kPol2GausExp) {
    double A, mu, sigma, lambda;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    sigma = this->fit->GetParameter(2);
    lambda = this->fit->GetParameter(3);
    a = this->fit->GetParameter(4);
    b = this->fit->GetParameter(5);
    if (this->inputs.fitType == InputSettings::kPol2GausExp) c = this->fit->GetParameter(6);

    // TODO: Change this to 2 functions: G->e and pol1/2
    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3]))", this->inputs.fitmin, mu+sigma*lambda);
    f1->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Exp(-(x - ([1]+[2]/[3])/2.)/([2]/[3])) * (x > [1]+[2]*[3])", mu+sigma*lambda, this->inputs.fitmax);
    f2->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    if (this->inputs.fitType == InputSettings::kPol1GausExp) {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs.fitmin, this->inputs.fitmax);
      f3->SetParameters(a, b);
      parts.push_back(f3);
    } else {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs.fitmin, this->inputs.fitmax);
      f3->SetParameters(a, b, c);
      parts.push_back(f3);
    }
  } else if (this->inputs.fitType == InputSettings::kPol1GausGaus || this->inputs.fitType == InputSettings::kPol2GausGaus) {
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
    if (this->inputs.fitType == InputSettings::kPol2GausGaus) c = this->fit->GetParameter(7);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs.fitmin, this->inputs.fitmax);
    f1->SetParameters(A, mu, sigma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs.fitmin, this->inputs.fitmax);
    f2->SetParameters(B, mu, lambda);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    if (this->inputs.fitType == InputSettings::kPol1GausGaus) {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs.fitmin, this->inputs.fitmax);
      f3->SetParameters(a, b);
      parts.push_back(f3);
    } else {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs.fitmin, this->inputs.fitmax);
      f3->SetParameters(a, b, c);
      parts.push_back(f3);
    }
  } else if (this->inputs.fitType == InputSettings::kPol1GausGausExp || this->inputs.fitType == InputSettings::kPol2GausGausExp) {
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
    if (this->inputs.fitType == InputSettings::kPol2GausExp) c = this->fit->GetParameter(8);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2]) * (x < ([1]+[2]*[3]))", this->inputs.fitmin, mu+sigma*lambda);
    f1->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Exp(-(x - ([1]+[2]/[3])/2.)/([2]/[3])) * (x > [1]+[2]*[3])", mu+sigma*lambda, this->inputs.fitmax);
    f2->SetParameters(A, mu, sigma, lambda);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    TF1* f3 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs.fitmin, this->inputs.fitmax);
    f3->SetParameters(B, mu, rho);
    parts.push_back(f3);

    s = TString::Format("%s_%s", fName.c_str(), "f4").Data();
    if (this->inputs.fitType == InputSettings::kPol1GausGausExp) {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs.fitmin, this->inputs.fitmax);
      f4->SetParameters(a, b);
      parts.push_back(f4);
    } else {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs.fitmin, this->inputs.fitmax);
      f4->SetParameters(a, b, c);
      parts.push_back(f4);
    }
  } else if (this->inputs.fitType == InputSettings::kPol1GausGausXex || this->inputs.fitType == InputSettings::kPol2GausGausXex) {
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
    if (this->inputs.fitType == InputSettings::kPol2GausGausXex) c = this->fit->GetParameter(10);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs.fitmin, this->inputs.fitmax);
    f1->SetParameters(A, mu, sigma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs.fitmin, this->inputs.fitmax);
    f2->SetParameters(B, mu, rho);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    TF1* f3 = new TF1(s.c_str(), "max(0., [0]*(x-[1])*TMath::Exp(-(x-[1])/[2]))", this->inputs.fitmin, this->inputs.fitmax);
    f3->SetParameters(C, nu, tau);
    parts.push_back(f3);

    s = TString::Format("%s_%s", fName.c_str(), "f4").Data();
    if (this->inputs.fitType == InputSettings::kPol1GausGausXex) {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs.fitmin, this->inputs.fitmax);
      f4->SetParameters(a, b);
      parts.push_back(f4);
    } else {
      TF1* f4 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs.fitmin, this->inputs.fitmax);
      f4->SetParameters(a, b, c);
      parts.push_back(f4);
    }
  } else if (this->inputs.fitType == InputSettings::kPol1GausXex || this->inputs.fitType == InputSettings::kPol2GausXex) {
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
    if (this->inputs.fitType == InputSettings::kPol2GausXex) c = this->fit->GetParameter(8);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Gaus(x,[1],[2])", this->inputs.fitmin, this->inputs.fitmax);
    f1->SetParameters(A, mu, sigma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    TF1* f2 = new TF1(s.c_str(), "max(0., [0]*(x-[1])*TMath::Exp(-(x-[1])/[2]))", this->inputs.fitmin, this->inputs.fitmax);
    f2->SetParameters(B, nu, tau);
    parts.push_back(f2);

    s = TString::Format("%s_%s", fName.c_str(), "f3").Data();
    if (this->inputs.fitType == InputSettings::kPol1GausXex) {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs.fitmin, this->inputs.fitmax);
      f3->SetParameters(a, b);
      parts.push_back(f3);
    } else {
      TF1* f3 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs.fitmin, this->inputs.fitmax);
      f3->SetParameters(a, b, c);
      parts.push_back(f3);
    }
  } else if (this->inputs.fitType == InputSettings::kPol1Voigt || this->inputs.fitType == InputSettings::kPol2Voigt) {
    double A, mu, sigma, Gamma;
    double a, b, c;
    A = this->fit->GetParameter(0);
    mu = this->fit->GetParameter(1);
    sigma = this->fit->GetParameter(2);
    Gamma = this->fit->GetParameter(3);
    a = this->fit->GetParameter(4);
    b = this->fit->GetParameter(5);
    if (this->inputs.fitType == InputSettings::kPol2Voigt) c = this->fit->GetParameter(6);

    string s = TString::Format("%s_%s", fName.c_str(), "f1").Data();
    TF1* f1 = new TF1(s.c_str(), "[0]*TMath::Voigt(x-[1],[2],[3])", this->inputs.fitmin, this->inputs.fitmax);
    f1->SetParameters(A, mu, sigma, Gamma);
    parts.push_back(f1);

    s = TString::Format("%s_%s", fName.c_str(), "f2").Data();
    if (this->inputs.fitType == InputSettings::kPol1Voigt) {
      TF1* f2 = new TF1(s.c_str(), "[0]+[1]*x", this->inputs.fitmin, this->inputs.fitmax);
      f2->SetParameters(a, b);
      parts.push_back(f2);
    } else {
      TF1* f2 = new TF1(s.c_str(), "[0]+[1]*x+[2]*x*x", this->inputs.fitmin, this->inputs.fitmax);
      f2->SetParameters(a, b, c);
      parts.push_back(f2);
    }
  } else {
    cout << "InputSettings::setFitParts() Error: do not know the parts of this function" << endl;
  }
  this->fitParts = parts;
  return parts;
}

TH1* MassFitter::loadFitParams() {
  if (!this->fit)
    return nullptr;

  string saveName = this->inputs.getSaveNameFromPt("fitParams");

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

TH1* MassFitter::loadFitResults() {
  if (!this->fit)
    return nullptr;

  string saveName = this->inputs.getSaveNameFromPt("fitResults");

  int nBins = 12;
  TH1* h = new TH1D(saveName.c_str(), "fitResults", nBins, -0.5, 1.*nBins - 0.5);
  h->SetTitle("Fit results");
  h->GetXaxis()->SetTitle("");
  h->GetXaxis()->SetBinLabel(1, "#chi^{2} (fit)");
  h->GetXaxis()->SetBinLabel(2, "NDF (fit)");
  h->GetXaxis()->SetBinLabel(3, "#chi^{2}/NDF (fit)");
  h->GetXaxis()->SetBinLabel(4, "#chi^{2} (all)");
  h->GetXaxis()->SetBinLabel(5, "NDF (all)");
  h->GetXaxis()->SetBinLabel(6, "#chi^{2}/NDF (all)");
  h->GetXaxis()->SetBinLabel(7, "Sig+Bkg");
  h->GetXaxis()->SetBinLabel(8, "Sig");
  h->GetXaxis()->SetBinLabel(9, "Bkg");
  h->GetXaxis()->SetBinLabel(10, "Sig/(Sig+Bkg)");
  h->GetXaxis()->SetBinLabel(11, "Sig/Bkg");
  h->GetXaxis()->SetBinLabel(12, "Sig/sqrt(Sig+Bkg)");

  double chisqFit = chisqInRange(this->data, this->fit, this->inputs.fitmin, this->inputs.fitmax);
  double ndfFit = ndfInRange(this->data, this->fit, this->inputs.fitmin, this->inputs.fitmax);
  double chisqAll = chisqInRange(this->data, this->fit, this->data->GetXaxis()->GetXmin(), this->data->GetXaxis()->GetXmax());
  double ndfAll = ndfInRange(this->data, this->fit, this->data->GetXaxis()->GetXmin(), this->data->GetXaxis()->GetXmax());
  this->calcSigBkg();

  h->SetBinContent(1, chisqFit);
  h->SetBinContent(2, ndfFit);
  h->SetBinContent(3, chisqFit / ndfFit);
  h->SetBinContent(4, chisqAll);
  h->SetBinContent(5, ndfAll);
  h->SetBinContent(6, chisqAll / ndfAll);
  h->SetBinContent(7, this->signalPlusBackground);
  h->SetBinContent(8, this->signal);
  h->SetBinContent(9, this->background);
  h->SetBinContent(10, this->signal / this->signalPlusBackground);
  h->SetBinContent(11, this->signal / this->background);
  h->SetBinContent(12, this->signal / sqrt(this->signalPlusBackground));

  this->fitResults = (TH1*)h->Clone(saveName.c_str());
  return h;
}

TH1* MassFitter::loadResidualHist() {
  if (!this->data || !this->fit)
    return nullptr;

  string saveName = this->inputs.getSaveNameFromPt("residual");

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

  string saveName = this->inputs.getSaveNameFromPt("pull");

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
  if (this->inputs.inputIssue("setHistNameFromTrain", "hadron"))
    return "";

  string s = "jet-fragmentation";
  if (this->inputs.train == 287744) s += "_id12406";
  s += "/data/V0/V0PtMass";

  this->inputs.histName = s;
  return s;
}

int MassFitter::setRebinNumberFromHadronAndPt() {
  int r = -1;

  if (this->inputs.hadron == "K0S" && this->inputs.lowpt > 25 - 1e3)
    r = 2;

  this->inputs.rebinNumber = r;
  return r;
}

void MassFitter::setSignalRegionFromSigma() {
  if (this->fitParts.size() == 0) {
    cout << "MassFitter::setSignalRegionFromSigma() Error: fit parts not set" << endl;
    return;
  }

  TF1* f = this->fitParts[0];
  double mu = f->GetParameter(1);
  double sigma = f->GetParameter(2);
  double nSigma = this->inputs.nSigma;
  this->inputs.setSignalRegion(mu - nSigma * sigma, mu + nSigma * sigma);
}

void MassFitter::writeOutputsToFile() {
  this->inputs.writeOutputToFile(this->data);
  this->inputs.writeOutputToFile(this->fit);
  for (auto f : this->fitParts) {
    this->inputs.writeOutputToFile(f);
  }
  this->inputs.writeOutputToFile(this->fitResults);
  this->inputs.writeOutputToFile(this->fitParams);
  this->inputs.writeOutputToFile(this->pull);
  this->inputs.writeOutputToFile(this->residual);
}

// -------------------------------------------------------------------------------------------------
//
// Struct for summarising fit information into TH2
//
// -------------------------------------------------------------------------------------------------

struct FitSummariser {
  private:
  public:
    InputSettings inputs;

    FitSummariser() { InputSettings(); }
    FitSummariser(InputSettings& x) { this->inputs = x; }

    void summariseFitInfo();
    void summariseFitInfo(string inputHistName);
};

void FitSummariser::summariseFitInfo() {
  this->summariseFitInfo(this->inputs.histName);
}

void FitSummariser::summariseFitInfo(string inputHistName) {
  int nPtBins = this->inputs.ptBinEdges.size();
  int nFitInfos = -1;
  if (inputHistName == "fitParams") {
    TF1* f = new TF1("f", this->inputs.getFitExpression().c_str(), 0., 1.);
    nFitInfos = f->GetNpar();
  } else if (inputHistName == "fitResults") {
    nFitInfos = 12;
  } else {
    cout << "writeFitInfo: unknown hist name " << inputHistName << endl;
    return;
  }

  TH2* fitInfo = new TH2D(inputHistName.c_str(), inputHistName.c_str(), nPtBins, -0.5, nPtBins - 0.5, nFitInfos, -0.5, nFitInfos - 0.5);
  for (int ix = 0; ix < nPtBins; ix++) {
    fitInfo->GetXaxis()->SetBinLabel(ix + 1, TString::Format("%.1f < #it{p}_{T, V0} < %.1f", this->inputs.ptBinEdges[ix][0], this->inputs.ptBinEdges[ix][1]).Data());
  }
  bool binLabelsSet = false;

  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!file) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return;
  }

  for (int iPt = 0; iPt < nPtBins; iPt++) {
    this->inputs.setPt(this->inputs.ptBinEdges[iPt][0], this->inputs.ptBinEdges[iPt][1]);
    string histName = this->inputs.getSaveNameFromPt(inputHistName);
    TH1* hist = (TH1*)file->Get(histName.c_str());
    if (!hist) {
      cout << "Could not find histogram " << histName << " in file " << this->inputs.inputFileName << endl;
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
  this->inputs.writeOutputToFile(fitInfo);
}

// -------------------------------------------------------------------------------------------------
//
// Struct for plotting
//
// -------------------------------------------------------------------------------------------------

struct FitPlotter {
  private:
  public:
    InputSettings inputs;
    MassFitter massFitter;

    TCanvas* canvas = nullptr;
    TH1F* frame = nullptr;
    TLegend* legend = nullptr;

    FitPlotter() { InputSettings(); MassFitter(); }
    FitPlotter(InputSettings& x) { this->inputs = x; this->massFitter = MassFitter(x); }
    FitPlotter(MassFitter& x) { this->massFitter = x; this->inputs = x.inputs; }

    TH1* loadMassHist();
    vector<TF1*> loadFitParts();
    TF1* loadFit();
    void plotFitInfo();
    int plotFitInfo(int infoToPlot);
    void plotFitParts();
};

vector<TF1*> FitPlotter::loadFitParts() {
  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!file) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return {};
  }

  vector<TF1*> parts = {};
  int maxParts = 10; // To prevent infinite loop
  for (int i = 1; i < maxParts; i++) {
    string fName = this->inputs.getSaveNameFromPt("fit", "_f" + to_string(i));
    TF1* f = (TF1*)file->Get(fName.c_str());
    if (!f)
      break;

    parts.push_back(f);
  }

  this->massFitter.fitParts = parts;
  return parts;
}

TF1* FitPlotter::loadFit() {
  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!file) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return nullptr;
  }

  string fName = this->inputs.getSaveNameFromPt("fit");
  TF1* f = (TF1*)file->Get(fName.c_str());
  this->massFitter.fit = (TF1*)f->Clone();
  return f;
}

TH1* FitPlotter::loadMassHist() {
  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");
  if (!file) {
    cout << "Could not open file " << this->inputs.inputFileName << endl;
    return nullptr;
  }

  string histName = this->inputs.getSaveNameFromPt("data");
  TH1* hist = (TH1*)file->Get(histName.c_str());
  if (!hist) {
    cout << "Could not find histogram " << histName << " in file " << this->inputs.inputFileName << endl;
    return nullptr;
  }
  this->massFitter.data = (TH1*)hist->Clone(histName.c_str());
  return hist;
}

void FitPlotter::plotFitInfo() {
  int maxInfos = 20; // To prevent infinite loop
  for (int i = 0; i < maxInfos; i++) {
    if (this->plotFitInfo(i) != 0)
      break;
  }
}

int FitPlotter::plotFitInfo(int infoToPlot)
{
  TFile* file = TFile::Open(this->inputs.inputFileName.c_str(), "READ");

  TH2* fitInfo = (TH2*)file->Get(this->inputs.histName.c_str());
  if (!fitInfo) {
    cout << "FitPlotter::plotFitInfo() Error: Could not find histogram " << this->inputs.histName << " in file " << this->inputs.inputFileName << endl;
    return 1;
  }

  if (infoToPlot < 0 || infoToPlot >= fitInfo->GetNbinsY()) {
    cout << "FitPlotter::plotFitInfo() Error: Invalid column " << infoToPlot << " selected for plotting " << this->inputs.histName << endl;
    cout << "Valid range is 0 to " << fitInfo->GetNbinsY() - 1 << endl;
    return 2;
  }

  TH1* hist = (TH1*)fitInfo->ProjectionX("hist", infoToPlot+1, infoToPlot+1);
  setStyle(hist, 0);
  if (!this->canvas) {
    cout << "FitPlotter::plotFitIfo() Warning: Canvas not set, creating automatically" << endl;
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetLogy(this->inputs.logplot);
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
    histTitle = TString::Format("%s, %s, ", getDataSet(this->inputs.train).c_str(), this->inputs.fitName.c_str()).Data();
    if (this->inputs.histName == "fitParams")
      histTitle += "par " + to_string(infoToPlot);
    if (this->inputs.histName == "fitResults")
      histTitle += fitInfo->GetYaxis()->GetBinLabel(infoToPlot + 1);

    yMinFrame = getHistLowerBound(hist, false);
    yMinFrame *= (yMinFrame < 0) ? 1.1 : 0.9;
    yMaxFrame = 1.1 * getHistUpperBound(hist, false);
  }

  hist->SetTitle(histTitle.c_str());
  hist->SetStats(0);
  hist->SetMinimum(yMinFrame);
  hist->SetMaximum(yMaxFrame);

  string outputFileName = this->inputs.outputFileName;
  if (outputFileName == "") {
    outputFileName = this->inputs.hadron + "_" + this->inputs.fitName + "_" + this->inputs.histName + to_string(infoToPlot) + ".pdf";
    cout << "FitPlotter::plotFitInfo() Warning: outputFileName not set, using: " << outputFileName << endl;
  }
  hist->Draw();
  hist->Draw("same text45");
  canvas->SaveAs(outputFileName.c_str());
  return 0;
}

void FitPlotter::plotFitParts() {
  if (this->inputs.outputFileName == "") {
    this->inputs.outputFileName = this->inputs.getSaveNameFromPt(this->inputs.hadron + "_" + this->inputs.fitName, ".pdf");
    cout << "FitPlotter::plotFitParts() Warning: outputFileName not set, using: " << this->inputs.outputFileName << endl;
  }

  if (!this->canvas) {
    cout << "FitPlotter::plotFitParts() Warning: Canvas not set, creating automatically" << endl;
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetLogy(this->inputs.logplot);
    this->canvas = canvas;
  }

  this->canvas->cd();
  if (!this->frame) {
    cout << "FitPlotter::plotFitParts() Warning: Frame not set, creating automatically" << endl;
    double xMinFrame = this->massFitter.data->GetXaxis()->GetXmin();
    double xMaxFrame = this->massFitter.data->GetXaxis()->GetXmax();
    double yMinFrame = 0.;
    double yMaxFrame = 1.1 * getHistUpperBound(this->massFitter.data, false);
    string xTitle = TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(this->inputs.hadron).c_str()).Data();
    string yTitle = "";
    string histTitle = TString::Format("%s, (%.1f < #it{p}_{T, V0} < %.1f)", getDataSet(this->inputs.train).c_str(), this->inputs.lowpt, this->inputs.highpt).Data();

    TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
    frame->SetTitle(histTitle.c_str());
    this->frame = frame;
  }

  this->frame->Draw();
  if (this->legend) this->legend->Draw("same");

  for (int i = 0; i < this->massFitter.fitParts.size(); i++) {
    TF1* f = this->massFitter.fitParts[i];
    setStyle(f, i + 2);
    f->SetRange(this->massFitter.data->GetXaxis()->GetXmin(), this->massFitter.data->GetXaxis()->GetXmax());
    f->Draw("same");
  }
  // Only do this if auto setup is requested
  setStyle(this->massFitter.data, 0);
  setStyle(this->massFitter.fit, 1);
  this->massFitter.fit->SetRange(this->massFitter.data->GetXaxis()->GetXmin(), this->massFitter.data->GetXaxis()->GetXmax());

  this->massFitter.data->Draw("same");
  this->massFitter.fit->Draw("same");
  this->canvas->SaveAs(this->inputs.outputFileName.c_str());
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
  p.loadMassHist();
  p.loadFit();
  p.loadFitParts();
  p.plotFitParts();
}

// Plot the fit parts by using the fitter
void plotFitParts(MassFitter& m) {
  FitPlotter p(m);
  p.inputs.outputFileName = p.inputs.getSaveNameFromPt(p.inputs.hadron + "_" + p.inputs.fitName, ".pdf");
  p.plotFitParts();
}

// Plot the fit information by reading the file
void plotFitInfo(InputSettings& x) {
  FitPlotter p(x);
  p.inputs.histName = "fitParams";
  p.plotFitInfo();
  p.inputs.histName = "fitResults";
  p.plotFitInfo();
}

// Summarise the fit information by reading the file
void summariseFitInfo(InputSettings& x) {
  FitSummariser f(x);
  f.summariseFitInfo("fitParams");
  f.summariseFitInfo("fitResults");
}

// This runs the entire workflow in order
void quickRun(InputSettings& x) {
  vector<vector<double>> ptBinEdges = x.ptBinEdges;

  // For each pt bin, do fit, and plot fit parts
  for (int iPt = 0; iPt < ptBinEdges.size(); iPt++) {
    x.setPt(ptBinEdges[iPt][0], ptBinEdges[iPt][1]);
    MassFitter m(x);
    m.doFitting();

    FitPlotter p(m);
    p.inputs.outputFileName = x.getSaveNameFromPt(x.hadron + "_" + x.fitName, ".pdf");
    p.plotFitParts();
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
  p.inputs.histName = "fitParams";
  p.plotFitInfo();
  p.inputs.histName = "fitResults";
  p.plotFitInfo();
}

void quickRun() {
  InputSettings x;
  x.train = 252064;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausXex");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitXFromHadron();
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();
  x.outputFileName = x.hadron + "_" + x.fitName + ".root";

  x.setPtBinEdgesFromHadron();
  quickRun(x);
}

void singleRun() {
  InputSettings x;
  x.setPt(3., 4.);
  x.train = 252064;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausXex");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitXFromHadron();
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();
  x.outputFileName = x.hadron + "_" + x.fitName + ".root";

  MassFitter m(x);
  // m.doFitting();
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  printParLimits(m.fit);
  m.data->Fit(m.fit, "R");
  printParLimits(m.fit);
  // m.fixFitInPost(); // Applies `any post-fit fixes, like swapping gaussians
  m.loadFitParts();
  m.loadFitParams();
  m.loadFitResults();
  m.loadResidualHist();
  m.loadPullHist();
  m.writeOutputsToFile();

  FitPlotter p(m);
  p.inputs.outputFileName = x.getSaveNameFromPt(x.hadron + "_" + x.fitName, ".pdf");
  p.plotFitParts();
}

#endif
