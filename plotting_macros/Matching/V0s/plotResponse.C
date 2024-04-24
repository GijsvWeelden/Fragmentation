
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

double getNjets(TFile* inFile, double jetptmin, double jetptmax)
{
  string histName = "jet-fragmentation/matching/jets/matchPartJetPtEtaPhi";
  TH3D* jets = (TH3D*)inFile->Get(histName.c_str());
  std::array<int, 2> jetptbins = getProjectionBins(jets->GetXaxis(), jetptmin, jetptmax);
  return jets->Integral(jetptbins[0], jetptbins[1], 0, jets->GetNbinsY() + 1, 0, jets->GetNbinsZ() + 1);
}

// -------------------------------------------------------------------------------------------------

void plotResponse(string inName = "", string hadron = "", double partjetptmin = 10., double partjetptmax = 1e6, bool detector = false, bool doZ = false)
{
  if ("" == inName) {
    cout << "Error: inName must be specified" << endl;
    return;
  }
  if ( ("K0S" == hadron) + ("Lambda0" == hadron) + ("AntiLambda0" == hadron) != 1 ) {
    cout << "Error: hadron must be either V0, K0S, Lambda0, or AntiLambda0" << endl;
    return;
  }
  const int nDim               = 7;
  const int partJetPtAxis      = 0;
  const int partV0PtAxis       = 1;
  const int detJetPtAxis       = 2;
  const int detV0PtAxis        = 3;
  const int K0SMassAxis        = 4;
  const int LambdaMassAxis     = 5;
  const int AntiLambdaMassAxis = 6;

  gStyle->SetNdivisions(505, "xy");
  string saveName, histName, histTitle, xTitle, yTitle, legendTitle, latexText, dataSet;
  double textSize = 0.04;
  double labelSize = 0.04;
  double titleSize = 0.04;
  bool setLogY = true;
  double xMinFrame = 0., xMaxFrame = 100., yMinFrame = 1e-9, yMaxFrame = 1.;
  double xMinLegend = 0.5, xMaxLegend = 0.9, yMinLegend = 0.6, yMaxLegend = 0.8;
  double xLatex = 0.35, yLatex = 0.8;
  int xCanvas = 900, yCanvas = 900;
  int rebinNumber = 5;
  xTitle = TString::Format("#it{p}_{T, %s}^{%s.} (GeV/#it{c})", formatHadronName(hadron).c_str(), detector ? "det" : "part").Data();
  if (doZ) {
    xTitle = TString::Format("#it{z}_{%s}^{%s.}", formatHadronName(hadron).c_str(), detector ? "det" : "part").Data();
    xMinFrame = 1e-3; xMaxFrame = 1.+1e-3;
  }
  yTitle = "#frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{p}_{T}}";
  if (doZ) { yTitle = "#frac{1}{#it{N}_{jets}} #frac{d#it{N}}{d#it{z}}"; }
  dataSet = "LHC24b1";

  std::vector<TH1D*> histVector;
  TCanvas* canvas = new TCanvas("Plot", "Plot", xCanvas, yCanvas);
  if (setLogY) { canvas->SetLogy(); }
  TH1F* frame = DrawFrame(xMinFrame, xMaxFrame, yMinFrame, yMaxFrame, xTitle, yTitle);
  TLegend* legend = CreateLegend(xMinLegend, xMaxLegend, yMinLegend, yMaxLegend, legendTitle, textSize);
  TLatex* latex;

  histName = TString::Format("partJetPt%sPtDetJetPt%sPtAllMasses", hadron.c_str(), hadron.c_str()).Data();
  if (doZ) { histName = TString::Format("partJetPt%sTrackProjDetJetPt%sTrackProjAllMasses", hadron.c_str(), hadron.c_str()); }
  histName = TString::Format("jet-fragmentation/matching/jets/V0/%s", histName.c_str());
  TFile *inFile = TFile::Open(TString::Format("./%s", inName.c_str()).Data());
  THnSparseD* thn = (THnSparseD*)inFile->Get(histName.c_str());

  std::array<int, 2> partjetptbins = getProjectionBins(thn->GetAxis(partJetPtAxis), partjetptmin, partjetptmax);
  thn->GetAxis(partJetPtAxis)->SetRange(partjetptbins[0], partjetptbins[1]);
  double lowjetpt = thn->GetAxis(partJetPtAxis)->GetBinLowEdge(partjetptbins[0]);
  double highjetpt = thn->GetAxis(partJetPtAxis)->GetBinUpEdge(partjetptbins[1]);

  TH1D* RMprojection = (TH1D*)thn->Projection(detector ? detV0PtAxis : partV0PtAxis);
  RMprojection->SetName("RMprojection");
  RMprojection->Scale(1./getNjets(inFile, lowjetpt, highjetpt), "width");
  setStyle(RMprojection, 0);
  histVector.push_back(RMprojection);

  latexText = TString::Format("#splitline{ %s }{ #it{p}_{T, jet}^{part.} = %.0f - %.0f (GeV/#it{c}) }", dataSet.c_str(), lowjetpt, highjetpt).Data();
  latex = CreateLatex(xLatex, yLatex, latexText, textSize);

  saveName = TString::Format("RM%s_%s%s", hadron.c_str(), detector ? "det" : "part" , doZ ? "Z" : "Pt").Data();
  saveName = TString::Format("%s_partjetpt%.0f-%.0f", saveName.c_str(), lowjetpt, highjetpt);
  saveName = TString::Format("%s.pdf", saveName.c_str());
  plotNHists(canvas, frame, histVector, legend, latex, saveName, "");
}
