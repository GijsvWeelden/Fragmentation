
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TString.h"
#include "TLegend.h"

#include "../../histUtils.C"
#include "../../plotUtils.C"

//
// Plots for V0s in the perpendicular cones
//

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

struct InputSettings{
  private:
    verbosityutilities::Verbosity verbosity = verbosityutilities::kInfo;

    string getNameFromVar(string prefix, string varstring, string suffix) {
      return prefix + varstring + suffix;
    }
    template <typename T> bool setVariable(T a, T b, T &x, T &y);
  public:
    // Analysis settings
    int train;
    string hadron = "";
    string histName = "";
    string inputFileName = "";
    string outputFileName = "";

    double ptmin = -1e3, ptmax = -1e3, zmin = -1e3, zmax = -1e3;
    double ptjetmin = -1e3, ptjetmax = -1e3;
    double etamin = -0.9, etamax = 0.9;
    vector<vector<double>> ptBinEdges = {};

    // Plotting settings
    double textSize = 0.04;
    bool logplot = false;
    bool ratioplot = false;

    // Setters and getters
    verbosityutilities::Verbosity getVerbosity() { return verbosity; }
    bool setVerbosity(verbosityutilities::Verbosity v);
    bool setEta(double a, double b) { return setVariable(a, b, etamin, etamax); }
    bool setJetPt(double a, double b) { return setVariable(a, b, ptjetmin, ptjetmax); }
    bool setPt(double a, double b) { return setVariable(a, b, ptmin, ptmax); }

    // Functions
    string getNameFromJetPt(string prefix, string suffix) { return getNameFromVar(prefix, TString::Format("ptjet%.f-%.f", ptjetmin, ptjetmax).Data(), suffix); }
    bool passVerbosityCheck(verbosityutilities::Verbosity v) { return verbosityutilities::passVerbosityCheck(v, verbosity); }
    bool printLog(string message, verbosityutilities::Verbosity messageVerbLevel) { return verbosityutilities::printLog(message, messageVerbLevel, verbosity); }
    string setInputFileNameFromTrain() {
      string s = "~/cernbox/TrainOutput/" + to_string(train) + "/AnalysisResults.root";
      inputFileName = s;
      return s;
    }
    string setHistName(string v0type) {
      if (v0type != "Fake" && v0type != "Matched") {
        printLog("InputSettings::setHistName() Error: Invalid V0 type: " + v0type, verbosityutilities::kErrors);
        return "";
      }
      histName = "jet-fragmentation/matching/PC/jetPtEta" + v0type + "V0Pt";
      return histName;
    }
};

bool InputSettings::setVerbosity(verbosityutilities::Verbosity v) {
  if (!verbosityutilities::is_valid(v))
    return false;

  verbosity = v;
  return true;
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

// ----------------------------------

double getNjets(InputSettings& inputs) {
  TFile* file = TFile::Open(inputs.inputFileName.c_str());
  if (!file) {
    inputs.printLog(TString::Format("getNjets() Error: Could not open file %s", inputs.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return -1.;
  }
  TH2D* jets = (TH2D*)file->Get("jet-fragmentation/matching/jets/matchDetJetPtPartJetPt");
  if (!jets) {
    inputs.printLog("getNjets() Error: Could not find histogram jet-fragmentation/matching/jets/matchDetJetPtPartJetPt in file " + inputs.inputFileName, verbosityutilities::kErrors);
    return -1.;
  }

  array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), inputs.ptjetmin, inputs.ptjetmax);
  double nJets = jets->Integral(jetptbins[0], jetptbins[1], 0, 1 + jets->GetYaxis()->GetNbins(), "width");
  return nJets;
}

TH3D* getPerpConeHist(InputSettings& inputs) {
  TFile* file = TFile::Open(inputs.inputFileName.c_str());
  if (!file) {
    inputs.printLog("getPerpConeHist() Error: Could not open file " +  inputs.inputFileName, verbosityutilities::kErrors);
    return nullptr;
  }
  TH3D* h3 = (TH3D*)file->Get(inputs.histName.c_str());
  if (!h3) {
    inputs.printLog("getPerpConeHist() Error: Could not find histogram " + inputs.histName + " in file " + inputs.inputFileName, verbosityutilities::kErrors);
    return nullptr;
  }
  return h3;
}

void plotPerpCone(InputSettings& inputs) {
  TH3D* h3 = getPerpConeHist(inputs);
  if (!h3)
    return;

  plotutils::Plotter p(inputs.outputFileName, inputs.logplot, inputs.textSize);

  // TODO: Loop over jet pt bins
  array<int, 2> jetptbins = getProjectionBins(h3->GetXaxis(), inputs.ptjetmin, inputs.ptjetmax);
  array<int, 2> etabins   = getProjectionBins(h3->GetYaxis(), inputs.etamin, inputs.etamax);
  TH1D* h1 = h3->ProjectionZ("v0pt", jetptbins[0], jetptbins[1], etabins[0], etabins[1]);
  TH1* v0s = rebinHist(h1, rebinnedV0PtHist(inputs.hadron, inputs.getNameFromJetPt("v0pt", "")));
  v0s->Scale(1./(2.*getNjets(inputs)), "width"); // 2 cones per jet
  p.addHistogram(v0s);
  p.setHistStyles();
  p.plot();
}

InputSettings setup(int train) {
  InputSettings x;
  x.setVerbosity(verbosityutilities::kDebug);
  x.train = train;
  x.setInputFileNameFromTrain();
  x.hadron = "K0S";
  return x;
}

void plotPerpCone() {
  InputSettings x = setup(499300);
  x.setHistName("Matched");
  plotPerpCone(x);
  // Set jet pt bins

  // x.setHistName("Fake");
}

// -------------------------------------------------------------------------------------------------
//
// Old code
//
// -------------------------------------------------------------------------------------------------

// void myPlotter(TCanvas* c, TH1F* frame, vector<TObject*> objList)
// {
//   c->cd();
//   frame->Draw();
//   for (auto obj : objList) {
//     obj->Draw("same");
//   }
//   c->SaveAs(c->GetName());
// }

// void plotMass(array<string, 4> inputStrings, double ptmin, double ptmax, bool isTrueV0)
// {
//   gStyle->SetNdivisions(505);
//   string fileName = inputStrings[0];
//   string histName = inputStrings[1];
//   string dataSet  = inputStrings[2];
//   string hadron   = inputStrings[3];
//   if (("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1) {
//     cout << "Invalid input! Hadron (" << hadron << ") should be K0S, Lambda0 or AntiLambda0" << endl;
//     return;
//   }

//   const int nDim           = 4;
//   const int ptAxis         = 0;
//   const int K0SAxis        = 1;
//   const int LambdaAxis     = 2;
//   const int AntiLambdaAxis = 3;

//   int projectionAxis = K0SAxis;
//   if ("Lambda0" == hadron) projectionAxis = LambdaAxis;
//   if ("AntiLambda0" == hadron) projectionAxis = AntiLambdaAxis;

//   TFile* f = TFile::Open(fileName.c_str());
//   THnSparse* thn = (THnSparse*)f->Get(histName.c_str());

//   array<int, 2> ptBins = getProjectionBins(thn->GetAxis(ptAxis), ptmin, ptmax);
//   thn->GetAxis(ptAxis)->SetRange(ptBins[0], ptBins[1]);
//   double lowpt  = thn->GetAxis(ptAxis)->GetBinLowEdge(ptBins[0]);
//   double highpt = thn->GetAxis(ptAxis)->GetBinUpEdge(ptBins[1]);
//   TH1D* mass = thn->Projection(projectionAxis);
//   mass->SetName(Form("mass_%s_%s_pt%.1f-%.1f", hadron.c_str(), histName.c_str(), lowpt, highpt));
//   setStyle(mass, 0);

//   TLine* massLine = new TLine(MassK0S, 0., MassK0S, 1.1*mass->GetMaximum());
//   massLine->SetLineColor(GetColor(1));
//   massLine->SetLineStyle(9);
//   massLine->SetLineWidth(3);

//   string canvasName = "PCV0" + (isTrueV0 ? "True" : "Fake");
//   canvasName += TString::Format("_m%s_pt%.1f-%.1f", hadron.c_str(), lowpt, highpt).Data();
//   canvasName += ".pdf";
//   TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 900, 900);

//   string xTitle = "#it{M}("; xTitle += formatHadronDaughters(hadron).c_str(); xTitle += ") (GeV/#it{c}^{2})";
//   string yTitle = "Counts";

//   double xMinFrame = mass->GetXaxis()->GetXmin(), xMaxFrame = mass->GetXaxis()->GetXmax(),
//          yMinFrame = 0., yMaxFrame = 1.1 * mass->GetBinContent(mass->GetMaximumBin());
//   string text = dataSet;
//   text += TString::Format(", #it{p}_{T} = %.1f-%.1f GeV/#it{c}", lowpt, highpt).Data();
//   TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle.c_str(), yTitle.c_str());
//   frame->SetTitle(text.c_str());

//   vector<TObject*> objList = {mass, massLine};
//   for (auto l : latex) { objList.push_back(l); }
//   myPlotter(c, frame, objList);
// }

// void plot24b1bMatched(string hadron, double ptmin, double ptmax)
// {
//   string fileName = "~/cernbox/TrainOutput/267466/AnalysisResults.root";
//   string histName = "jet-fragmentation/matching/PC/matchedV0PtMass";
//   string dataSet  = "LHC24b1b";
//   array<string, 4> inputStrings = {fileName, histName, dataSet, hadron};
//   plotMass(inputStrings, ptmin, ptmax, true);
// }
// void plot24b1bFake(string hadron, double ptmin, double ptmax)
// {
//   string fileName = "~/cernbox/TrainOutput/267466/AnalysisResults.root";
//   string histName = "jet-fragmentation/matching/PC/fakeV0PtMass";
//   string dataSet  = "LHC24b1b";
//   array<string, 4> inputStrings = {fileName, histName, dataSet, hadron};
//   TLatex* latex = new TLatex(0.15, 0.85, "Combinatorial V0");
//   plotMass(inputStrings, ptmin, ptmax, false);
// }