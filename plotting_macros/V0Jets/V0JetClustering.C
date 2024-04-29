
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

const double MassK0S = 0.497611;
const double MassLambda0 = 1.115683;

double rap(double* x, double* par)
{
  double E = sqrt(par[0] * par[0] + x[0] * x[0]);
  return 0.5 * TMath::Log((E/x[0] + x[1]) / (E/x[0] - x[1]));
}

double deltaRap(double x, double y)
{
  TF2* yK = new TF2("yK", "rap", 0., 5., -1., 1., 1);
  yK->SetParameter(0, MassK0S);
  TF2* yL = new TF2("yL", "rap", 0., 5., -1., 1., 1);
  yL->SetParameter(0, MassLambda0);
  return yL->Eval(x, y) - yK->Eval(x, y);
}

void plotE()
{
  TF1* E = new TF1("E", "sqrt([0]*[0] + x*x)", 0., 20.);
  E->SetParameter(0, MassK0S);
  setStyle(E, 0);
  E->Draw();
}

void ploty()
{
  TF2* yK = new TF2("yK", "rap", 0., 5., -1., 1., 1);
  yK->SetParameter(0, MassK0S);
  setStyle(yK, 0);
  TF2* yL = new TF2("yL", "rap", 0., 5., -1., 1., 1);
  yL->SetParameter(0, MassLambda0);
  setStyle(yL, 1);

  // yK->Draw("colz");
  // TF2* deltaY = new TF2("deltaY", "rap - rap", 0., 5., 0., 5., 2);
  TF2* deltaY = new TF2("deltaY", deltaRap, 0., 5., 0., 5., 0);
  // deltaY->SetParameter(0, MassK0S);
  // deltaY->SetParameter(1, MassLambda0);
  // TF2* deltaY = new TF2("deltaY", yL - yK, 0., 5., 0., 5.);
  deltaY->Draw("colz");
}