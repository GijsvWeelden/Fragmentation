
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

#include "../histUtils.C"
#include "../plotUtils.C"
#include "../myStrings.C"

// This script compares jets when we assign Lambdas the K0S mass or the Lambda mass:
// * Jet pt spectrum
// * K0S z spectrum
// * Lambda z spectrum

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

namespace analysisutilities {
  enum AnaType {kJetLasK, kJetSchemes, kZK0SLasK, kZK0SSchemesV0, kZK0SSchemesK0, kZLambdaLasK, kZLambdaSchemesV0, kZLambdaSchemesK0};
}

struct Names {
  public:
    string inputFileName;
    string histogramName;
    string legendEntry;
    string jetHistName;

    Names() {}
    Names(string file, string hist, string legend) : inputFileName(file), histogramName(hist), legendEntry(legend) {}
    Names(string file, string hist, string legend, string jetHist) : inputFileName(file), histogramName(hist), legendEntry(legend), jetHistName(jetHist) {}
};

struct InputSettings{
  private:
    verbosityutilities::Verbosity verbosity = verbosityutilities::kInfo;

    string getNameFromVar(string prefix, string varstring, string suffix) {
      return prefix + varstring + suffix;
    }
    template <typename T> bool setVariable(T a, T b, T &x, T &y);
  public:
    // Analysis settings
    vector<Names> names = {};
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
    void addNames(string file, string hist, string legend) { names.push_back({file, hist, legend}); }
    void addNames(string file, string hist, string legend, string jethist) { names.push_back({file, hist, legend, jethist}); }
    string getNameFromJetPt(string prefix, string suffix) { return getNameFromVar(prefix, TString::Format("ptjet%.f-%.f", ptjetmin, ptjetmax).Data(), suffix); }
    bool passVerbosityCheck(verbosityutilities::Verbosity v) { return verbosityutilities::passVerbosityCheck(v, verbosity); }
    bool printLog(string message, verbosityutilities::Verbosity messageVerbLevel) { return verbosityutilities::printLog(message, messageVerbLevel, verbosity); }
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

struct Plotter {
  private:
    string drawOption; // Private because it requires caution with spaces

  public:
    InputSettings* inputs;

    TCanvas* canvas = nullptr;
    TH1F* frame = nullptr;
    TLegend* legend = nullptr;
    vector<TH1*> hists = {};
    vector<TObject*> objects = {};

    Plotter() { inputs = new InputSettings(); }
    Plotter(InputSettings& x) { inputs = &x; }

    void addLatex(double x, double y, string s);
    void addLine(double x0, double y0, double x1, double y1, int styleNumber, int lineStyle, int lineWidth);
    void makeCanvas(string s, double x, double y);
    void makeFrame(string sx, string sy);
    void makeFrame(double x0, double x1, double y0, double y1, string sx, string sy);
    void makeLegend(double x0, double x1, double y0, double y1, string s);
    void makeRatios(int baseIndex = 0);
    void plot();
    void reset();
    void setHistStyles();

    string getMassString();
    string setDrawOption(string s);
};

void Plotter::addLatex(double x, double y, string s) {
  TLatex* l = plotutils::CreateLatex(x, y, s.c_str(), inputs->textSize);
  objects.push_back(l);
}
void Plotter::addLine(double x0, double y0, double x1, double y1, int styleNumber, int lineStyle = 9, int lineWidth = 3) {
  TLine* l = new TLine(x0, y0, x1, y1);
  plotutils::setStyle(l, styleNumber, lineStyle, lineWidth);
  objects.push_back(l);
}
void Plotter::makeCanvas(string s = "canvas", double x = 800, double y = 600) {
  canvas = new TCanvas(s.c_str(), s.c_str(), x, y);
  canvas->SetLogy(inputs->logplot);
}
void Plotter::makeFrame(string sx, string sy) {
  if (hists.empty()) {
    inputs->printLog("Plotter::makeFrame() hists vector is empty! Aborting", verbosityutilities::kErrors);
    return;
  }
  if (!canvas) makeCanvas();

  double xMinFrame = hists[0]->GetXaxis()->GetXmin();
  double xMaxFrame = hists[0]->GetXaxis()->GetXmax();
  double yMinFrame = getLowerBound(hists, 0, 0) * 0.9;
  double yMaxFrame = getUpperBound(hists, 0, 0) * 1.2;

  inputs->printLog(TString::Format("Plotter::makeFrame() x: %.2f, %.2f \ny: %.2f, %.2f", xMinFrame, xMaxFrame, yMinFrame, yMaxFrame).Data(), verbosityutilities::kDebug);

  if (inputs->logplot) {
    roundToNextPowerOfTen(yMaxFrame);
    if (yMinFrame < 1e-12) {
      int xBin = hists[0]->FindLastBinAbove(0.);
      yMinFrame = hists[0]->GetBinContent(xBin) * 0.5;
      yMaxFrame *= 1.75;
      roundToPrevPowerOfTen(yMinFrame);
    }
  }
  inputs->printLog(TString::Format("Plotter::makeFrame() x: %.2f, %.2f \ny: %.2f, %.2f", xMinFrame, xMaxFrame, yMinFrame, yMaxFrame).Data(), verbosityutilities::kDebug);
  makeFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, sx, sy);
}
void Plotter::makeFrame(double x0, double x1, double y0, double y1, string sx, string sy) {
  if (!canvas) makeCanvas();
  frame = plotutils::DrawFrame(x0, x1, y0, y1, sx, sy);
}
void Plotter::makeLegend(double x0, double x1, double y0, double y1, string s) {
  legend = plotutils::CreateLegend(x0, x1, y0, y1, s.c_str(), inputs->textSize);
}

void Plotter::makeRatios(int baseIndex = 0) {
  if (hists.empty()) {
    inputs->printLog("Plotter::makeRatios(): Hist vector is empty!", verbosityutilities::kErrors);
    return;
  }
  TH1* baseCopy = (TH1*)hists[baseIndex]->Clone("baseCopy");
  for (auto& h : hists) h = divideWithProtection(h, baseCopy);
}

void Plotter::plot() {
  if (!canvas) makeCanvas();
  if (!frame) makeFrame("x", "y");

  frame->Draw();
  if (legend) legend->Draw("same");
  for (auto o : objects)
    o->Draw("same");

  for (auto h : hists) {
    h->Draw(("same" + drawOption).c_str());
    if (inputs->passVerbosityCheck(verbosityutilities::kDebugMax))
      h->Print("all");
    else if (inputs->passVerbosityCheck(verbosityutilities::kDebug))
      h->Print();
  }
  if (inputs->outputFileName != "")
    canvas->SaveAs(inputs->outputFileName.c_str());
}

void Plotter::reset() {
  canvas = nullptr;
  frame = nullptr;
  legend = nullptr;
  objects.clear();
  hists.clear();
}

void Plotter::setHistStyles() {
  for (unsigned int i = 0; i < hists.size(); i++)
    plotutils::setStyle(hists[i], i);
}

string Plotter::getMassString() {
  return TString::Format("#it{M}(%s) (GeV/#it{c}^{2})", formatHadronDaughters(inputs->hadron).c_str()).Data();
}
string Plotter::setDrawOption(string s) {
  if (s.at(0) == ' ') { // Single quotes, because checking for a char
    drawOption = s;
  } else {
    drawOption = " " + s;
  }
  return drawOption;
}

double getNevts(TFile* f, bool withV0s = false)
{
  TH1D* hNEvents = (TH1D*)f->Get("hNEvts");
  return hNEvents->GetBinContent(1 + (int)withV0s);
}

double getNjets(TH1* h, int minBin, int maxBin)
{
  return h->Integral(minBin, maxBin);
}

void comparePt(vector<string> inNames, vector<string> histNames, vector<string> legendEntries, TLatex* additionalLatex = nullptr, bool doRatio = false)
{
  double minEta = -0.35, maxEta = 0.35;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 150., yMinFrame = 1e-1, yMaxFrame = 2.;
  double xMinLegend = 0.25, xMaxLegend = 0.65, yMinLegend = 0.25, yMaxLegend = 0.4;
  double xLatex = 0.25, yLatex = 0.81;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = "#it{p}_{T, jet}";
  yTitle = "normalised count";
  yTitle = "#it{N}_{jets}";
  drawoption = "same";
  saveName = "Pythia-V0JetClustering-jetpt";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = plotutils::CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  latexText = TString::Format("#splitline{ %s }{#splitline{ %s }{#splitline{ %s }{ %s }}}",
                                "Pythia simulation pp",
                                "#sqrt{#it{s}} = 13.6 TeV",
                                "Anti-k_{T} jets, #it{R} = 0.4",
                                TString::Format("|#it{#eta}_{jet}| < %.2f", maxEta).Data()
                              );
  TLatex* latex = plotutils::CreateLatex(xLatex, yLatex, latexText, textSize);

  double scale = 1.;
  for (int i = 0; i < histNames.size(); i++) {
    string inName = inNames[i];
    TFile *inFile = TFile::Open(TString::Format("%s", inName.c_str()).Data());
    string histName = histNames[i];
    TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
    array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), minEta, maxEta);
    TH1D* th1 = jets->ProjectionX("jetpt", etabins[0], etabins[1]);
    th1->SetName(TString::Format("jetpt_%s_%s", histName.c_str(), inName.c_str()).Data());
    plotutils::setStyle(th1, i);
    th1->Rebin(rebinNumber);
    th1->Sumw2();
    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);
    if (th1->GetBinContent(th1->GetMaximumBin()) > scale) {
      scale = th1->GetBinContent(th1->GetMaximumBin());
    }
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.9; yMaxFrame = 1.1;
    scale = 1.;
    drawoption = "same hist";
    yTitle = TString::Format("ratio w.r.t. %s", legendEntries[0].c_str()).Data();
    for (int i = 1; i < histVector.size(); i++) {
      histVector[i]->Divide(histVector[0]);
    }
    histVector[0]->Divide(histVector[0]);
  }

  TH1F* frame = plotutils::DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame * scale, xTitle, yTitle);
  frame->Draw();
  for (int i = 0; i < histVector.size(); i++) {
    TH1D* hist = histVector[i];
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");
  latex->Draw("same");
  if (additionalLatex) additionalLatex->Draw("same");

  saveName += (doRatio ? "Ratio" : "Comparison");
  saveName += ".pdf";
  canvas->SaveAs(saveName.c_str());
}

// /*
void compareZ(vector<string> inNames, vector<string> histNames, vector<string> njetsNames, vector<string> legendEntries, string hadron, double jetptmin = 10., double jetptmax = 1e3, TLatex* additionalLatex = nullptr, bool doRatio = false)
{
  double minEta = -0.35, maxEta = 0.35;
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("V0" == hadron) != 1) {
    cout << "Hadron " << hadron << " not recognised. Should be K0S, Lambda0 or V0" << endl;
    return;
  }

  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int jetetaAxis = 1;
  const int jetphiAxis = 2;
  const int zAxis      = 3;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, drawoption;
  double lowjetpt, highjetpt;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 1e-3, xMaxFrame = 1.+1e-3, yMinFrame = 1e-3, yMaxFrame = 1e3;
  double xMinLegend = 0.25, xMaxLegend = 0.55, yMinLegend = 0.2, yMaxLegend = 0.35;
  double xLatex = 0.25, yLatex = 0.82;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 10;
  xTitle = "#it{z}_{" + formatHadronName(hadron) + "}";
  yTitle = TString::Format("#frac{1}{#it{N}_{jets}} #frac{d #it{N}}{d #it{z}_{%s}}", formatHadronName(hadron).c_str()).Data();
  drawoption = "same";
  saveName = "Pythia-z";
  saveName += hadron;

  vector<double> njets = { 0., 0. };
  vector<double> integrals; vector<double> nevents;
  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TLegend* legend = plotutils::CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);

  for (int i = 0; i < histNames.size(); i++) {
    string inName = inNames[i];
    string histName = histNames[i];
    TFile *inFile = TFile::Open(inName.c_str());
    THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

    std::array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), jetptmin, jetptmax);
    thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
    lowjetpt = thn->GetAxis(jetptAxis)->GetBinLowEdge(jetptbins[0]);
    highjetpt = thn->GetAxis(jetptAxis)->GetBinUpEdge(jetptbins[1]);
    array<int, 2> etabins = getProjectionBins(thn->GetAxis(jetetaAxis), minEta, maxEta);
    thn->GetAxis(jetetaAxis)->SetRange(etabins[0], etabins[1]);
    TH1D* th1 = (TH1D*)thn->Projection(zAxis);
    th1->SetName(TString::Format("z_%s_%s", histName.c_str(), inName.c_str()).Data());
    th1->Rebin(rebinNumber);
    th1->Sumw2();
    plotutils::setStyle(th1, i);
    legend->AddEntry(th1, legendEntries[i].c_str());
    histVector.push_back(th1);

    TH3D* jetPtEtaPhi = (TH3D*)inFile->Get(njetsNames[i].c_str());
    jetptbins = getProjectionBins(jetPtEtaPhi->GetXaxis(), jetptmin, jetptmax);
    etabins = getProjectionBins(jetPtEtaPhi->GetYaxis(), minEta, maxEta);
    TH1D* jetPt = jetPtEtaPhi->ProjectionX("jetpt", etabins[0], etabins[1]);
    jetPt->SetName(TString::Format("jetpt_%s", histName.c_str()).Data());
    njets[i] = getNjets(jetPt, jetptbins[0], jetptbins[1]);
  }
  if (doRatio) {
    canvas->SetLogy(false);
    yMinFrame = 0.9; yMaxFrame = 1.1;
    drawoption = "same hist";
    yTitle = TString::Format("ratio w.r.t. %s", legendEntries[0].c_str()).Data();
    for (int i = 1; i < histVector.size(); i++) {
      histVector[i]->Divide(histVector[0]);
    }
    histVector[0]->Divide(histVector[0]);
  }

  TH1F* frame = plotutils::DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  frame->Draw();
  for (int i = 0; i < histVector.size(); i++) {
    TH1D* hist = histVector[i];
    if (!doRatio) {
      // hist->Scale(1./getNevts(inFile), "width");
      hist->Scale(1./njets[i], "width");
    }
    hist->Draw(drawoption.c_str());
  }
  legend->Draw("same");

  latexText = TString::Format("#splitline{ %s }{#splitline{ %s }{#splitline{ %s }{ %s }}}",
                                "Pythia simulation pp",
                                "#sqrt{#it{s}} = 13.6 TeV",
                                TString::Format("Anti-k_{T} jets, #it{R} = 0.4, |#it{#eta}_{jet}| < %.1f", maxEta).Data(),
                                TString::Format("%.0f < #it{p}_{T, jet} < %.0f GeV/#it{c}", lowjetpt, highjetpt).Data()
                              );
  TLatex* latex = plotutils::CreateLatex(xLatex, yLatex, latexText, textSize);
  latex->Draw("same");
  if (additionalLatex) additionalLatex->Draw("same");

  saveName += (doRatio ? "_Ratio" : "_Comparison");
  saveName = TString::Format("%s_jetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName += ".pdf";
  canvas->SaveAs(saveName.c_str());
}
// */

string getNameZPlot(InputSettings& inputs, analysisutilities::AnaType type) {
  string prefix;
  switch (type) {
    case analysisutilities::kZK0SLasK:
    case analysisutilities::kZLambdaLasK:
      prefix = "LasK";
      break;
    case analysisutilities::kZK0SSchemesV0:
    case analysisutilities::kZLambdaSchemesV0:
      prefix = "schemesV0";
      break;
    case analysisutilities::kZK0SSchemesK0:
    case analysisutilities::kZLambdaSchemesK0:
      prefix = "schemesK0";
      break;
    default:
      inputs.printLog("getNameZPlot() Error: Analysis type not recognised!", verbosityutilities::kErrors);
      return "";
  }
  prefix += "-z" + inputs.hadron + "_";
  string suffix = (inputs.ratioplot) ? "Ratio.pdf" : "Comparison.pdf";
  return inputs.getNameFromJetPt(prefix, suffix);
}

double getNevts(Names names, bool withV0s = false) {
  TFile* f = TFile::Open(names.inputFileName.c_str());
  TH1D* hNEvents = (TH1D*)f->Get("hNEvts");
  return hNEvents->GetBinContent(1 + (int)withV0s);
}

double getNjets(InputSettings& inputs, Names names) {
  TFile* f = TFile::Open(names.inputFileName.c_str());
  if (!f) {
    inputs.printLog(TString::Format("getNjets() Error: Could not open file %s", names.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return -1;
  }
  TH3* jets = (TH3*)f->Get(names.jetHistName.c_str());
  if (!jets) {
    inputs.printLog(TString::Format("getNjets() Error: Could not find histogram %s in file %s", names.jetHistName.c_str(), names.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return -1;
  }
  array<int, 2> ptbins  = getProjectionBins(jets->GetXaxis(), inputs.ptjetmin, inputs.ptjetmax);
  array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), inputs.etamin, inputs.etamax);
  TH1* jetpt = jets->ProjectionX("jetpt", etabins[0], etabins[1]);
  return jetpt->Integral(ptbins[0], ptbins[1]);
}

TH1* getJetPtHist(InputSettings& inputs, Names name, const TH1* rebinTemplate = nullptr) {
  TFile *inFile = TFile::Open(name.inputFileName.c_str());
  if (!inFile) {
    inputs.printLog(TString::Format("getJetPtHist() Error: Could not open file %s", name.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return nullptr;
  }
  TH3* jets = (TH3*)inFile->Get(name.histogramName.c_str());
  if (!jets) {
    inputs.printLog(TString::Format("getJetPtHist() Error: Could not find histogram %s in file %s", name.histogramName.c_str(), name.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return nullptr;
  }

  array<int, 2> etabins = getProjectionBins(jets->GetYaxis(), inputs.etamin, inputs.etamax);
  TH1* th1 = jets->ProjectionX("jetpt", etabins[0], etabins[1]);
  th1->SetName(TString::Format("jetpt_%s_%s", name.histogramName.c_str(), name.inputFileName.c_str()).Data());
  th1->Sumw2();
  if (rebinTemplate) {
    TH1* tmp = rebinHist(th1, rebinTemplate);
    th1 = tmp;
  }
  return th1;
}

TH1* getV0ZHist(InputSettings& inputs, Names name, const TH1* rebinTemplate = nullptr) {
  TFile *inFile = TFile::Open(name.inputFileName.c_str());
  if (!inFile) {
    inputs.printLog(TString::Format("getV0ZHist() Error: Could not open file %s", name.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return nullptr;
  }
  THnSparseD* thn = (THnSparseD*)inFile->Get(name.histogramName.c_str());
  if (!thn) {
    inputs.printLog(TString::Format("getV0ZHist() Error: Could not find histogram %s in file %s", name.histogramName.c_str(), name.inputFileName.c_str()).Data(), verbosityutilities::kErrors);
    return nullptr;
  }

  const int nDim       = 4;
  const int jetptAxis  = 0;
  const int jetetaAxis = 1;
  const int jetphiAxis = 2;
  const int zAxis      = 3;

  array<int, 2> jetptbins = getProjectionBins(thn->GetAxis(jetptAxis), inputs.ptjetmin, inputs.ptjetmax);
  array<int, 2> etabins   = getProjectionBins(thn->GetAxis(jetetaAxis), inputs.etamin, inputs.etamax);
  thn->GetAxis(jetptAxis)->SetRange(jetptbins[0], jetptbins[1]);
  thn->GetAxis(jetetaAxis)->SetRange(etabins[0], etabins[1]);

  TH1* th1 = (TH1*)thn->Projection(zAxis);
  th1->SetName(TString::Format("z_%s_%s", name.histogramName.c_str(), name.inputFileName.c_str()).Data());
  th1->Sumw2();
  if (rebinTemplate) {
    TH1* tmp = rebinHist(th1, rebinTemplate);
    th1 = tmp;
  }
  return th1;
}

void compareJetPtLasK(InputSettings& inputs) {
  TH1* jetPtTemplate = new TH1D("jetpttemplate", "jetpttemplate", 12, 0., 60.);
  Plotter p(inputs);
  for (int i = 0; i < inputs.names.size(); i++) {
    TH1* h = getJetPtHist(inputs, inputs.names[i], jetPtTemplate);
    h->Scale(1./getNevts(inputs.names[i]), "width");
    p.hists.push_back(h);
  }

  if (inputs.ratioplot) {
    p.addLatex(0.35, 0.85, mystrings::sThisThesis);
    p.addLatex(0.35, 0.80, mystrings::sPythia + " simulation pp");
    p.addLatex(0.35, 0.75, mystrings::sSqrtS);
    p.addLatex(0.35, 0.70, mystrings::sAntikt + " jets, #it{E}-scheme");
    p.addLatex(0.35, 0.65, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);

    string xTitle = mystrings::sPtJet + " (" + mystrings::sGevC + ")";
    string yTitle = "Ratio";
    p.makeFrame(0., 60., 0.9, 1.1, xTitle, yTitle);
    p.makeRatios(0);
    p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  } else {
    p.addLatex(0.50, 0.85, mystrings::sThisThesis);
    p.addLatex(0.50, 0.80, mystrings::sPythia + " simulation pp");
    p.addLatex(0.50, 0.75, mystrings::sSqrtS);
    p.addLatex(0.50, 0.70, mystrings::sAntikt + " jets, #it{E}-scheme");
    p.addLatex(0.50, 0.65, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);

    string xTitle = mystrings::sPtJet + " (" + mystrings::sGevC + ")";
    string yTitle = mystrings::sJetsPerXsec;
    p.makeFrame(0., 60., 1e-8, 1e-3, xTitle, yTitle);
    p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  }

  for (int i = 0; i < p.hists.size(); i++) {
    p.legend->AddEntry(p.hists[i], inputs.names[i].legendEntry.c_str());
  }

  p.setHistStyles();
  p.plot();
}

void compareJetPtLasK() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string fileName = filePath + "v0jetclustering-Escheme.root";
  x.addNames(fileName, "hV0Jet", "#Lambda mass = #it{m}(#Lambda)");
  x.addNames(fileName, "hK0Jet", "#Lambda mass = #it{m}(K^{0}_{S})");

  x.setEta(-0.35, 0.35);

  x.outputFileName = "LasK-jetptComparison.pdf";
  x.ratioplot = false;
  x.logplot = true;
  compareJetPtLasK(x);

  x.outputFileName = "LasK-jetptRatio.pdf";
  x.ratioplot = true;
  x.logplot = false;
  compareJetPtLasK(x);
}

void compareJetPtSchemes(InputSettings& inputs) {
  TH1* jetPtTemplate = new TH1D("jetpttemplate", "jetpttemplate", 12, 0., 60.);
  Plotter p(inputs);
  for (int i = 0; i < inputs.names.size(); i++) {
    TH1* h = getJetPtHist(inputs, inputs.names[i], jetPtTemplate);
    h->Scale(1./getNevts(inputs.names[i]), "width");
    p.hists.push_back(h);
  }

  if (inputs.ratioplot) {
    p.addLatex(0.35, 0.85, mystrings::sThisThesis);
    p.addLatex(0.35, 0.80, mystrings::sPythia + " simulation pp");
    p.addLatex(0.35, 0.75, mystrings::sSqrtS);
    p.addLatex(0.35, 0.70, mystrings::sAntikt + " jets");
    p.addLatex(0.35, 0.65, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);

    string xTitle = mystrings::sPtJet + " (" + mystrings::sGevC + ")";
    string yTitle = "Ratio";
    p.makeFrame(0., 60., 0.7, 1.4, xTitle, yTitle);
    p.makeRatios(0);
  } else {
    p.addLatex(0.45, 0.80, mystrings::sThisThesis);
    p.addLatex(0.45, 0.75, mystrings::sPythia + " simulation pp");
    p.addLatex(0.45, 0.70, mystrings::sSqrtS);
    p.addLatex(0.45, 0.65, mystrings::sAntikt + " jets");
    p.addLatex(0.45, 0.60, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);

    string xTitle = mystrings::sPtJet + " (" + mystrings::sGevC + ")";
    string yTitle = mystrings::sJetsPerXsec;
    p.makeFrame(0., 60., 1e-8, 1e-3, xTitle, yTitle);
  }

  p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  for (int i = 0; i < p.hists.size(); i++) {
    p.legend->AddEntry(p.hists[i], inputs.names[i].legendEntry.c_str());
  }

  p.setHistStyles();
  p.plot();
}

void compareJetPtSchemes() {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string eFileName = filePath + "v0jetclustering-Escheme.root";
  string ptFileName = filePath + "v0jetclustering-ptscheme.root";
  x.addNames(eFileName, "hK0Jet", "#it{E}-scheme");
  x.addNames(ptFileName, "hK0Jet", "#it{p}_{T}-scheme");

  x.outputFileName = "schemes-jetptComparison.pdf";
  x.ratioplot = false;
  x.logplot = true;
  compareJetPtSchemes(x);

  x.outputFileName = "schemes-jetptRatio.pdf";
  x.ratioplot = true;
  x.logplot = false;
  compareJetPtSchemes(x);

  x.names.clear();
  x.addNames(eFileName, "hV0Jet", "#it{E}-scheme");
  x.addNames(ptFileName, "hV0Jet", "#it{p}_{T}-scheme");

  x.outputFileName = "schemesV0-jetptComparison.pdf";
  x.ratioplot = false;
  x.logplot = true;
  compareJetPtSchemes(x);

  x.outputFileName = "schemesV0-jetptRatio.pdf";
  x.ratioplot = true;
  x.logplot = false;
  compareJetPtSchemes(x);
}

Plotter setupPlotterV0Z(InputSettings& x, analysisutilities::AnaType setting) {
  Plotter p(x);
  p.reset();
  switch (setting) {
    case analysisutilities::kZK0SLasK:
    case analysisutilities::kZLambdaLasK:
      if (x.ratioplot) {
        string sZ = mystrings::getZString(formatHadronName(x.hadron));
        string xTitle = sZ;
        string yTitle = "Ratio";
        p.makeFrame(0., 1., 0.5, 2.0, xTitle, yTitle);
        p.frame->GetYaxis()->SetTitleOffset(0.8);
        gPad->SetLeftMargin(0.15);

        p.addLatex(0.30, 0.82, mystrings::sThisThesis);
        p.addLatex(0.30, 0.77, mystrings::sPythia + " simulation pp");
        p.addLatex(0.30, 0.72, mystrings::sSqrtS);
        p.addLatex(0.30, 0.67, mystrings::sAntikt + " jets, #it{E}-scheme");
        p.addLatex(0.30, 0.62, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);
        p.addLatex(0.30, 0.57, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax));
      } else {
        string sZ = mystrings::getZString(formatHadronName(x.hadron));
        string xTitle = sZ;
        string yTitle = mystrings::getOneOverString(mystrings::sNjets) + " " + mystrings::getdYdXString(mystrings::addSubscript(mystrings::sSigma, formatHadronName(x.hadron)), sZ);
        p.makeFrame(0., 1., 1e-3, 100, xTitle, yTitle);
        p.frame->GetYaxis()->SetTitleOffset(0.8);
        gPad->SetLeftMargin(0.15);

        p.addLatex(0.55, 0.80, mystrings::sThisThesis);
        p.addLatex(0.55, 0.75, mystrings::sPythia + " simulation pp");
        p.addLatex(0.55, 0.70, mystrings::sSqrtS);
        p.addLatex(0.55, 0.65, mystrings::sAntikt + " jets, #it{E}-scheme");
        p.addLatex(0.55, 0.60, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);
        p.addLatex(0.55, 0.55, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax));
      }
      break;
    case analysisutilities::kZK0SSchemesV0:
    case analysisutilities::kZK0SSchemesK0:
      if (x.ratioplot) {
        string sJets = mystrings::sAntikt + " jets, #Lambda mass = #it{m}";
        if (setting == analysisutilities::kZK0SSchemesV0 || setting == analysisutilities::kZLambdaSchemesV0)
          sJets += "(#Lambda)";
        else if (setting == analysisutilities::kZK0SSchemesK0 || setting == analysisutilities::kZLambdaSchemesK0)
          sJets += "(K^{0}_{S})";

        string sZ = mystrings::getZString(formatHadronName(x.hadron));
        string xTitle = sZ;
        string yTitle = "Ratio";
        p.makeFrame(0., 1., 0.5, 2.0, xTitle, yTitle);
        p.frame->GetYaxis()->SetTitleOffset(0.8);
        gPad->SetLeftMargin(0.15);

        p.makeLegend(0.60, 0.85, 0.70, 0.85, "");
        p.addLatex(0.20, 0.85, mystrings::sThisThesis);
        p.addLatex(0.20, 0.80, mystrings::sPythia + " simulation pp");
        p.addLatex(0.20, 0.75, mystrings::sSqrtS);
        p.addLatex(0.20, 0.70, sJets.c_str());
        p.addLatex(0.20, 0.65, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);
        p.addLatex(0.20, 0.60, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax));
      } else {
        string sZ = mystrings::getZString(formatHadronName(x.hadron));
        string xTitle = sZ;
        string yTitle = mystrings::getOneOverString(mystrings::sNjets) + " " + mystrings::getdYdXString(mystrings::addSubscript(mystrings::sSigma, formatHadronName(x.hadron)), sZ);
        p.makeFrame(0., 1., 1e-3, 100, xTitle, yTitle);

        p.addLatex(0.55, 0.80, mystrings::sThisThesis);
        p.addLatex(0.55, 0.75, mystrings::sPythia + " simulation pp");
        p.addLatex(0.55, 0.70, mystrings::sSqrtS);

        if (setting == analysisutilities::kZK0SSchemesV0 || setting == analysisutilities::kZLambdaSchemesV0)
          p.addLatex(0.55, 0.65, mystrings::sAntikt + " jets, #Lambda mass = #it{m}(#Lambda)");
        if (setting == analysisutilities::kZK0SSchemesK0 || setting == analysisutilities::kZLambdaSchemesK0)
          p.addLatex(0.55, 0.65, mystrings::sAntikt + " jets, #Lambda mass = #it{m}(K^{0}_{S})");

        p.addLatex(0.55, 0.60, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);
        p.addLatex(0.55, 0.55, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax));
      }
      break;
    case analysisutilities::kZLambdaSchemesV0:
    case analysisutilities::kZLambdaSchemesK0:
      if (x.ratioplot) {
        string sJets = mystrings::sAntikt + " jets, #Lambda mass = #it{m}";
        if (setting == analysisutilities::kZK0SSchemesV0 || setting == analysisutilities::kZLambdaSchemesV0)
          sJets += "(#Lambda)";
        else if (setting == analysisutilities::kZK0SSchemesK0 || setting == analysisutilities::kZLambdaSchemesK0)
          sJets += "(K^{0}_{S})";

        string sZ = mystrings::getZString(formatHadronName(x.hadron));
        string xTitle = sZ;
        string yTitle = "Ratio";
        p.makeFrame(0., 1., 0.5, 2.0, xTitle, yTitle);
        p.frame->GetYaxis()->SetTitleOffset(0.8);
        gPad->SetLeftMargin(0.15);

        p.makeLegend(0.35, 0.50, 0.20, 0.35, "");
        p.addLatex(0.20, 0.85, mystrings::sThisThesis);
        p.addLatex(0.20, 0.80, mystrings::sPythia + " simulation pp");
        p.addLatex(0.20, 0.75, mystrings::sSqrtS);
        p.addLatex(0.20, 0.70, sJets.c_str());
        p.addLatex(0.20, 0.65, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);
        p.addLatex(0.20, 0.60, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax));
      } else {
        string sZ = mystrings::getZString(formatHadronName(x.hadron));
        string xTitle = sZ;
        string yTitle = mystrings::getOneOverString(mystrings::sNjets) + " " + mystrings::getdYdXString(mystrings::addSubscript(mystrings::sSigma, formatHadronName(x.hadron)), sZ);
        p.makeFrame(0., 1., 1e-3, 100, xTitle, yTitle);

        p.addLatex(0.55, 0.80, mystrings::sThisThesis);
        p.addLatex(0.55, 0.75, mystrings::sPythia + " simulation pp");
        p.addLatex(0.55, 0.70, mystrings::sSqrtS);

        if (setting == analysisutilities::kZK0SSchemesV0 || setting == analysisutilities::kZLambdaSchemesV0)
          p.addLatex(0.55, 0.65, mystrings::sAntikt + " jets, #Lambda mass = #it{m}(#Lambda)");
        if (setting == analysisutilities::kZK0SSchemesK0 || setting == analysisutilities::kZLambdaSchemesK0)
          p.addLatex(0.55, 0.65, mystrings::sAntikt + " jets, #Lambda mass = #it{m}(K^{0}_{S})");

        p.addLatex(0.55, 0.60, mystrings::sRadius + ", " + mystrings::sEtaJetRange035);
        p.addLatex(0.55, 0.55, mystrings::getPtJetRangeString(x.ptjetmin, x.ptjetmax));
      }
      break;
    default:
      x.printLog("setupPlotterV0Z() Error: Analysis setting " + to_string(setting) + " not recognised!", verbosityutilities::kErrors);
  }
  return p;
}

void compareV0Z(InputSettings& inputs, Plotter& p) {
  TH1* zTemplate = new TH1D("ztemplate", "ztemplate", 10, 1e-3, 1.+1e-3);
  const int jetptAxis  = 0;
  const int jetetaAxis = 1;
  const int jetphiAxis = 2;
  const int zAxis      = 3;
  const int nDim       = 4;

  for (int i = 0; i < inputs.names.size(); i++) {
    TH1* h = getV0ZHist(inputs, inputs.names[i], zTemplate);
    h->Scale(1./getNjets(inputs, inputs.names[i]), "width");
    p.hists.push_back(h);
  }

  if (inputs.ratioplot) {
    p.makeRatios(0);
    if (!p.legend)
      p.makeLegend(0.25, 0.40, 0.20, 0.35, "");
  } else {
    p.makeLegend(0.25, 0.50, 0.20, 0.35, "");
  }

  for (int i = 0; i < p.hists.size(); i++) {
    p.legend->AddEntry(p.hists[i], inputs.names[i].legendEntry.c_str());
  }

  p.setHistStyles();
  p.plot();
}

void compareV0ZLasK(string h = "K0S") {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string fileName = filePath + "v0jetclustering-Escheme.root";
  vector<array<double, 2>> ptjetbins = { {10., 20.}, {20., 30.}, {30., 40.} };

  analysisutilities::AnaType option;
  if (h == "K0S") {
    option = analysisutilities::kZK0SLasK;
    x.hadron = "K0S";
    x.addNames(fileName, "hzV0_K0S", "#Lambda mass = #it{m}(#Lambda)", "hV0Jet");
    x.addNames(fileName, "hzK0_K0S", "#Lambda mass = #it{m}(K^{0}_{S})", "hK0Jet");
  } else if (h == "Lambda") { // Also AntiLambda
    option = analysisutilities::kZLambdaLasK;
    x.hadron = "Lambda";
    x.addNames(fileName, "hzV0_Lambda0", "#Lambda mass = #it{m}(#Lambda)", "hV0Jet");
    x.addNames(fileName, "hzK0_Lambda0", "#Lambda mass = #it{m}(K^{0}_{S})", "hK0Jet");
  }
  for (int i = 0; i < ptjetbins.size(); i++) {
    x.setJetPt(ptjetbins[i][0], ptjetbins[i][1]);
    x.ratioplot = false;
    x.logplot = true;
    x.outputFileName = getNameZPlot(x, option);
    Plotter p = setupPlotterV0Z(x, option);
    compareV0Z(x, p);

    x.ratioplot = true;
    x.logplot = false;
    x.outputFileName = getNameZPlot(x, option);
    p.reset();
    p = setupPlotterV0Z(x, option);
    compareV0Z(x, p);
  }
}

void compareV0ZSchemes(analysisutilities::AnaType option) {
  InputSettings x; x.setVerbosity(verbosityutilities::kDebug);

  string filePath = "../../inputfiles/pythia/V0Study/";
  string efileName = filePath + "v0jetclustering-Escheme.root";
  string ptfileName = filePath + "v0jetclustering-ptscheme.root";
  vector<array<double, 2>> ptjetbins = { {10., 20.}, {20., 30.}, {30., 40.} };

  switch (option) {
    case analysisutilities::kZK0SSchemesV0:
      x.hadron = "K0S";
      x.addNames(efileName, "hzV0_K0S", "#it{E}-scheme", "hV0Jet");
      x.addNames(ptfileName, "hzV0_K0S", "#it{p}_{T}-scheme", "hV0Jet");
      break;
    case analysisutilities::kZK0SSchemesK0:
      x.hadron = "K0S";
      x.addNames(efileName, "hzK0_K0S", "#it{E}-scheme", "hK0Jet");
      x.addNames(ptfileName, "hzK0_K0S", "#it{p}_{T}-scheme", "hK0Jet");
      break;
    case analysisutilities::kZLambdaSchemesV0:
      x.hadron = "Lambda";
      x.addNames(efileName, "hzV0_Lambda0", "#it{E}-scheme", "hV0Jet");
      x.addNames(ptfileName, "hzV0_Lambda0", "#it{p}_{T}-scheme", "hV0Jet");
      break;
    case analysisutilities::kZLambdaSchemesK0:
      x.hadron = "Lambda";
      x.addNames(efileName, "hzK0_Lambda0", "#it{E}-scheme", "hK0Jet");
      x.addNames(ptfileName, "hzK0_Lambda0", "#it{p}_{T}-scheme", "hK0Jet");
      break;
    default:
      cout << "compareV0ZSchemes() Error: option " << option << " not recognised" << endl;
      return;
  }
  for (int i = 0; i < ptjetbins.size(); i++) {
    x.setJetPt(ptjetbins[i][0], ptjetbins[i][1]);
    x.ratioplot = false;
    x.logplot = true;
    x.outputFileName = getNameZPlot(x, option);
    Plotter p = setupPlotterV0Z(x, option);
    compareV0Z(x, p);

    x.ratioplot = true;
    x.logplot = false;
    x.outputFileName = getNameZPlot(x, option);
    p.reset();
    p = setupPlotterV0Z(x, option);
    compareV0Z(x, p);
  }
}

void LasK_Escheme(string hadron, bool doRatio = false, double jetptmin = 10., double jetptmax = 1e3)
{
  string inName = "../../inputfiles/pythia/V0Study/v0jetclustering-Escheme.root";
  vector<string> jHistNames = { "hV0Jet", "hK0Jet"};
  vector<string> zHistNames = { "hzV0_" + hadron, "hzK0_" + hadron};
  vector<string> legendEntries = { "#Lambda mass = #it{m}_{#Lambda}", "#Lambda mass = #it{m}_{K^{0}_{S}}"};

  TLatex* additionalLatex = plotutils::CreateLatex(0.65, 0.81, "#it{E} scheme", 0.04);

  if ("" == hadron) {
    comparePt({inName, inName}, jHistNames, legendEntries, additionalLatex, doRatio);
  }
  else {
    compareZ({inName, inName}, zHistNames, jHistNames, legendEntries, hadron, jetptmin, jetptmax, additionalLatex, doRatio);
  }
}

void LasK_Ptscheme(string hadron, bool doRatio = false, double jetptmin = 10., double jetptmax = 1e3)
{
  string inName = "../../inputfiles/pythia/V0Study/v0jetclustering-ptscheme.root";
  vector<string> jHistNames = { "hV0Jet", "hK0Jet"};
  vector<string> zHistNames = { "hzV0_" + hadron, "hzK0_" + hadron};
  vector<string> legendEntries = { "#Lambda mass = #it{m}_{#Lambda}", "#Lambda mass = #it{m}_{K^{0}_{S}}"};

  TLatex* additionalLatex = plotutils::CreateLatex(0.65, 0.81, "#it{p}_{T} scheme", 0.04);

  if ("" == hadron) {
    comparePt({inName, inName}, jHistNames, legendEntries, additionalLatex, doRatio);
  }
  else {
    compareZ({inName, inName}, zHistNames, jHistNames, legendEntries, hadron, jetptmin, jetptmax, additionalLatex, doRatio);
  }
}

void EvsPtScheme(string hadron, bool V0Jets, bool doRatio = false, double jetptmin = 10., double jetptmax = 1e3)
{
  string EName = "../../inputfiles/pythia/V0Study/v0jetclustering-Escheme.root";
  string ptName = "../../inputfiles/pythia/V0Study/v0jetclustering-ptscheme.root";

  string jHistName = "hK0Jet";
  string zHistName = "hzK0_" + hadron;
  string legendEntry = "#Lambda mass = #it{m}_{K^{0}_{S}}";
  if (V0Jets) {
    jHistName = "hV0Jet";
    zHistName = "hzV0_" + hadron;
    legendEntry = "#Lambda mass = #it{m}_{#Lambda}";
  }

  vector<string> jHistNames = {jHistName, jHistName};
  vector<string> zHistNames = {zHistName, zHistName};
  vector<string> legendEntries = {"#it{E} scheme", "#it{p}_{T} scheme"};

  TLatex* additionalLatex = plotutils::CreateLatex(0.6, 0.6, legendEntry.c_str(), 0.04);

  if ("" == hadron) {
    comparePt({EName, ptName}, jHistNames, legendEntries, additionalLatex, doRatio);
  }
  else {
    compareZ({EName, ptName}, zHistNames, jHistNames, legendEntries, hadron, jetptmin, jetptmax, additionalLatex, doRatio);
  }
}
