
#include "plotMassFits.C"
#include "../myStrings.C"

// This macro uses plotMassFits.C to plot mass fits for different pt bins

// -------------------------------------------------------------------------------------------------
//
// K0S GausGaus
//
// -------------------------------------------------------------------------------------------------

void fitK0SGausGaus1_2() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(1., 2.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.fixFitInPost(); // Applies `any post-fit fixes, like swapping gaussians
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.createLegend(0.55, 0.75, 0.3, 0.5);
  p.fillLegendWithFitParts();

  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sEtaV0Range075);
  p.addLatex(xLatex, yLatex - 0.15, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGaus2_3() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(2., 3.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.fixFitInPost(); // Applies `any post-fit fixes, like swapping gaussians
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.addLatex(0.55, 0.85, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());
  p.plotFitParts();
}

void fitK0SGausGaus3_4() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(3., 4.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.fixFitInPost(); // Applies `any post-fit fixes, like swapping gaussians
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.addLatex(0.55, 0.85, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());
  p.plotFitParts();
}

void fitK0SGausGaus4_5() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(4., 5.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.addLatex(0.55, 0.85, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());
  p.plotFitParts();
}

void fitK0SGausGaus5_10() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(5., 10.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.addLatex(0.55, 0.85, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());
  p.plotFitParts();
}

void fitK0SGausGaus10_15() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(10., 15.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.addLatex(0.55, 0.85, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());
  p.plotFitParts();
}

void fitK0SGausGaus15_20() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(15., 20.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.addLatex(0.55, 0.85, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());
  p.plotFitParts();
}

void fitK0SGausGaus20_25() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(20., 25.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.addLatex(0.55, 0.85, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());
  p.plotFitParts();
}

void fitK0SGausGaus25_30() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(25., 30.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.addLatex(0.55, 0.85, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());
  p.plotFitParts();
}

void fitK0SGausGaus30_40() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(30., 40.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGaus");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.addLatex(0.55, 0.85, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());
  p.plotFitParts();
}

void fitK0SGausGaus() {
  gROOT->SetBatch(kTRUE);
  fitK0SGausGaus1_2();
  fitK0SGausGaus2_3();
  fitK0SGausGaus3_4();
  fitK0SGausGaus4_5();
  fitK0SGausGaus5_10();
  fitK0SGausGaus10_15();
  fitK0SGausGaus15_20();
  fitK0SGausGaus20_25();
  fitK0SGausGaus25_30();
  fitK0SGausGaus30_40();
}

// -------------------------------------------------------------------------------------------------
//
// K0S GausGausExp
//
// -------------------------------------------------------------------------------------------------

void fitK0SGausGausExp1_2() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(1., 2.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  // Try to stay close to GausGaus fit parameters
  m.fit->SetParLimits(0, 0.7, 1.); // A
  m.fit->SetParLimits(4, 0., 0.15); // B
  m.fit->SetParLimits(5, 1e-2, 1.); // rho
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  p.createLegend(0.55, 0.75, 0.3, 0.5);
  p.fillLegendWithFitParts();

  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, mystrings::sThisThesisAliceData);
  p.addLatex(xLatex, yLatex - 0.05, mystrings::sSqrtS);
  p.addLatex(xLatex, yLatex - 0.10, mystrings::sEtaV0Range075);
  p.addLatex(xLatex, yLatex - 0.15, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp2_3() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(2., 3.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  // Try to stay close to GausGaus fit parameters
  m.fit->SetParLimits(0, 0.7, 1.); // A
  m.fit->SetParLimits(4, 0., 0.15); // B
  m.fit->SetParLimits(5, 1e-2, 1.); // rho
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp3_4() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(3., 4.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  // Try to stay close to GausGaus fit parameters
  m.fit->SetParLimits(0, 0.7, 1.); // A
  m.fit->SetParLimits(4, 0., 0.2); // B
  m.fit->SetParLimits(5, 1e-2, 1.); // rho
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);
  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp4_5() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(4., 5.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  // Try to stay close to GausGaus fit parameters
  m.fit->SetParLimits(0, 0.7, 1.); // A
  m.fit->SetParLimits(4, 0., 0.2); // B
  m.fit->SetParLimits(5, 1e-2, 1.); // rho
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);

  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp5_10() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(5., 10.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  // Try to stay close to GausGaus fit parameters
  m.fit->SetParLimits(0, 0.6, 1.); // A
  m.fit->SetParLimits(4, 0., 0.25); // B
  m.fit->SetParLimits(5, 1e-2, 1.); // rho
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);

  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp10_15() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(10., 15.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);

  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp15_20() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(15., 20.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);

  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp20_25() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(20., 25.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);

  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp25_30() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(25., 30.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);

  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp30_40() {
  InputSettings x; x.verbosity = InputSettings::kInfo;
  x.setPt(30., 40.);
  x.train = 529130;
  x.hadron = "K0S";
  x.setFitType("pol1GausGausExp");
  x.normaliseData = true;
  x.nSigma = 3.;
  x.fixMu = true;

  x.setFitX(0.42, 0.55);
  x.setPolInitXFromHadron();
  x.setMassWindowDiffFromHadron();
  x.setInputFileNameFromTrain();

  MassFitter m(x);
  m.loadMassHist();
  m.loadFitFunction();
  m.setFitInitialValues();
  m.data->Fit(m.fit, "R");
  m.loadFitParts();
  m.loadFitParams();
  m.loadResidualHist();
  m.loadPullHist();

  FitPlotter p(m);

  double xLatex = 0.55, yLatex = 0.80;
  p.addLatex(xLatex, yLatex, TString::Format("%.f < #it{p}_{T, V0} < %.f GeV/#it{c}", p.inputs->lowpt, p.inputs->highpt).Data());

  p.plotFitParts();
}

void fitK0SGausGausExp() {
  gROOT->SetBatch(kTRUE);
  fitK0SGausGausExp1_2();
  fitK0SGausGausExp2_3();
  fitK0SGausGausExp3_4();
  fitK0SGausGausExp4_5();
  fitK0SGausGausExp5_10();
  fitK0SGausGausExp10_15();
  fitK0SGausGausExp15_20();
  fitK0SGausGausExp20_25();
  fitK0SGausGausExp25_30();
  fitK0SGausGausExp30_40();
}

// -------------------------------------------------------------------------------------------------
//
// K0S ExpGausExp
//
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
//
// Lambda GausGaus
//
// -------------------------------------------------------------------------------------------------