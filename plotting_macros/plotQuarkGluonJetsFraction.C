
#include "histUtils.C"
#include "plotUtils.C"
#include "myStrings.C"

enum JetType { mQuark = 0, mGluon, mInclusive };
struct InputSettings {
  private:
  string _inputFileName, _outputFileName;
  double _ptmin, _ptmax;

  public:
  string getInputFileName() { return _inputFileName; }
  string getOutputFileName() { return _outputFileName; }
  double getPtMin() { return _ptmin; }
  double getPtMax() { return _ptmax; }

  void setInputFileName(string name) { _inputFileName = name; }
  void setOutputFileName(string name) { _outputFileName = name; }
  void setPtMin(double ptmin) { _ptmin = ptmin; }
  void setPtMax(double ptmax) { _ptmax = ptmax; }
  void setPt(double ptmin, double ptmax) { setPtMin(ptmin); setPtMax(ptmax); }
};

TFile* GetFile(InputSettings& inputs) {
  TFile* file = TFile::Open(inputs.getInputFileName().c_str());
  if (!file) {
    cout << "Could not open file " << inputs.getInputFileName() << endl;
    return nullptr;
  }
  return file;
}

double GetNJets(InputSettings& inputs, JetType jt) {
  TFile* file = GetFile(inputs);
  if (!file)
    return -1.;

  TH2F* hNJetTypes = (TH2F*)file->Get("hNJetTypes");
  if (!hNJetTypes) {
    cout << "Could not find histogram hNJetTypes in file " << inputs.getInputFileName() << endl;
    return -1.;
  }

  array<int, 2> ptbins = histutils::getProjectionBins(hNJetTypes->GetXaxis(), inputs.getPtMin(), inputs.getPtMax());
  TH1F* hNJets = (TH1F*)hNJetTypes->ProjectionX("nJets", ptbins[0], ptbins[1]);

  double nGluons = hNJets->GetBinContent(1);
  double nQuarks = hNJets->GetBinContent(2);
  double nJets = nQuarks + nGluons;

  switch (jt) {
    case mQuark:     return nQuarks;
    case mGluon:     return nGluons;
    case mInclusive: return nJets;
  }
}

void plotQuarkGluonJetsFraction() {
  InputSettings inputs;
  inputs.setInputFileName("../inputfiles/pythia/PythiaResult_pthat20-80.root");
  inputs.setOutputFileName("qgJetsFraction.pdf");

  TH1F* jetTemplate = new TH1F("hJetPtTemplate", "hJetPtTemplate", 18, 10., 100.);

  TFile* inFile = GetFile(inputs);
  TH2F* hNJetTypes = (TH2F*)inFile->Get("hNJetTypes");
  TH1F* hgJets = (TH1F*)hNJetTypes->ProjectionY("hgJets", 1, 1);
  TH1F* hqJets = (TH1F*)hNJetTypes->ProjectionY("hqJets", 2, 2);
  TH1F* hnJets = (TH1F*)hNJetTypes->ProjectionY("hnJets", 1, 2);

  hgJets = (TH1F*)histutils::rebinHist(hgJets, jetTemplate);
  hqJets = (TH1F*)histutils::rebinHist(hqJets, jetTemplate);
  hnJets = (TH1F*)histutils::rebinHist(hnJets, jetTemplate);

  TH1F* hgFraction = (TH1F*)histutils::divideWithProtection(hgJets, hnJets);
  hgFraction->SetName("gluonFraction");
  TH1F* hqFraction = (TH1F*)histutils::divideWithProtection(hqJets, hnJets);
  hqFraction->SetName("quarkFraction");

  plotutils::Plotter p(inputs.getOutputFileName(), false, 0.04);
  p.makeFrame(jetTemplate->GetXaxis()->GetXmin(), jetTemplate->GetXaxis()->GetXmax(), 0., 1., mystrings::sPtJetWithUnits, "Fraction of jets");
  p.getFrame()->GetYaxis()->SetLabelOffset(0.02);
  p.setHists({hgFraction, hqFraction});
  p.setHistStyles();

  p.makeLegend(0.4, 0.6, 0.20, 0.30, "");
  p.addLegendEntry(hgFraction, "Gluon-initiated jets");
  p.addLegendEntry(hqFraction, "Quark-initiated jets");

  double xLatex = 0.4, yLatex = 0.875;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisPythiaSim);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sAntikt + " jets, " + mystrings::sRadius + ", " + mystrings::sEtaJetRange20);
  p.plot();
}