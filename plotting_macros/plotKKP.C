
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

double kkpParameter(string hadron, bool quark, int parameter, double sHat)
{
  if (parameter < 0 || parameter > 4) { cout << "kkpParameter: invalid input!" << endl; }
  double outPar = -999.;
  auto f3 = [](double a, double b, double c, double d, double x)
  {
    return a + b*x + c*x*x + d*x*x*x;
  };

  // Only charged pion implemented at this point
  if (hadron == "pi") {
    if (quark) { // Currently only up quark implemented
      switch (parameter) {
        case 0: // N
          outPar = f3(0.44809, -0.13828, -0.06951, 0.01354, sHat);
          break;
        case 1: // alpha
          outPar = f3(-1.47598, -0.30498, -0.01863, -0.12529, sHat);
          break;
        case 2: // beta
          outPar = f3(0.91338, 0.64145, 0.07270, -0.16989, sHat);
          break;
        case 3: // gamma
          outPar = f3(0, 0.07396, -0.07757, 0, sHat);
          break;
        default:
          break;
      }
    }
    else { // Gluon
      switch (parameter) {
        case 0: // N
          outPar = f3(3.73331, -3.16946, -0.47683, 0.70270 , sHat);
          break;
        case 1: // alpha
          outPar = f3(-0.74159, -0.51377, -0.19705, -0.17917, sHat);
          break;
        case 2: // beta
          outPar = f3(2.33092, 2.03394, -0.50764, -0.08565, sHat);
          break;
        case 3: // gamma
          outPar = f3(0, 0.09466, -0.10222, 0, sHat);
          break;
        default:
          break;
      }
    }
  }
  else { cout << "kkpParameter: this hadron type is not implemented" << endl; }
  return outPar;
}

// Load a theory prediction for fragmentation in TF1* format
TF1* loadKKP(string hadron, bool quark, double Q)
{
  // KKP: https://arxiv.org/pdf/hep-ph/0011155.pdf
  double Lambda = 213e-3; // In GeV, at NLO
  // double LambdaSquared = Lambda * Lambda;
  // double Q0Squared = 2; // In GeV^2
  double Q0 = sqrt(2); // In GeV
  // double N = 22.2815, alpha = 0.12732, beta = 6.13697, gamma = 0;
  // double sHat = log( log(QSquared / LambdaSquared) / log (Q0Squared / LambdaSquared) );
  double sHat = log( log(Q / Lambda) / log (Q0 / Lambda) );

  TF1* theory = new TF1("kkp", "[0]*TMath::Power(x, [1])*TMath::Power(1 - x, [2])*(1 + [3]/x)", 0, 1);
  double kkpParameters[4] = { 0, 0, 0, 0 };
  for (int iPar = 0; iPar <= 4; iPar++) {
    kkpParameters[iPar] = kkpParameter(hadron, quark, iPar, sHat);
  }
  theory->SetParameters(kkpParameters[0], kkpParameters[1], kkpParameters[2], kkpParameters[3]);
  if (quark) { cout << "KKP D(q->"; }
  else { cout << "KKP D(g->"; }
  cout << hadron << ") = " << kkpParameters[0] << " * x^" << kkpParameters[1] << " (1 - x)^" << kkpParameters[2] << " (1 + " << kkpParameters[3] << "/x)" << endl;
  // theory->Draw();
  return theory;
}

// Formats the hadron name to look nice (Greek letters, sub- and superscripts)
string formatHadronName(string hadron)
{
  string had = hadron;
  if (hadron == "pi"){
    had = "#pi^{#pm}";
  }
  else if (hadron == "pi0"){
    had = "#pi^{0}";
  }
  else if (hadron == "K0L"){
    had = "K^{0}_{L}";
  }
  else if (hadron == "K0S"){
    had = "K^{0}_{S}";
  }
  else if (hadron == "K0"){
    had = "K^{0}";
  }
  else if (hadron == "K"){
    had = "K^{#pm}";
  }
  else if (hadron == "Lambda0"){
    had = "#Lambda^{0}";
  }
  return had;
}

void plotKKP(double Q, string hadron = "pi")
{
  bool isQuark = true;
  TF1* q = loadKKP(hadron, isQuark, Q);
  isQuark = false;
  TF1* g = loadKKP(hadron, isQuark, Q);

  // Settings for example plot
  q->SetLineWidth(3);
  q->GetHistogram()->SetMinimum(1e-9); // Set the minimum so the plot looks nice
  q->GetHistogram()->GetXaxis()->SetTitle("z");
  q->GetHistogram()->GetYaxis()->SetTitle(TString::Format("#it{D}_{u}^{%s} (z, Q^{2} = (%.0f GeV)^{2})", formatHadronName(hadron).c_str(), Q).Data());
  q->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);

  g->SetLineWidth(3);
  g->SetLineColor(kBlack);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetLogy();
  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(q, "q");
  legend->AddEntry(g, "g");

  q->Draw();
  g->Draw("same");
  legend->Draw("same");
}
