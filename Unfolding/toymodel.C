#include <TROOT.h>
#include <TStyle.h>
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TParameter.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

using std::cout;
using std::endl;
using std::string;

// Using a 2D toy model, this macro illustrates how to perform unfolding with RooUnfold
// It shows a both trivial and statistically independent closure tests
// Run the full chain with: toymodel()
// There is also an official RooUnfold tutorial at the end of this file

bool AreHistsEqual(const TH1* histA, const TH1* histB, double tolerance = 1e-10) {
  vector<int> bins;
  vector<double> diffs;
  for (int i = 0; i < histA->GetNbinsX(); i++) {
    double binDiff = histA->GetBinContent(i) - histB->GetBinContent(i);
    if (abs(binDiff) > tolerance) {
      bins.push_back(i);
      diffs.push_back(binDiff);
    }
  }
  cout << histA->GetName() << " - " << histB->GetName();
  if (bins.size() == 0) {
    cout << ": All bins match within tolerance " << tolerance << "\n";
    return true;
  } else {
    for (int i = 0; i < bins.size(); i++) {
      cout << " Bin " << bins[i] << " = " << diffs[i] << " ";
    }
    cout << "\n";
    return false;
  }
}
bool AreHistsEqual(const TH2* histA, const TH2* histB, double tolerance = 1e-10) {
  vector<array<int, 2>> bins;
  vector<double> diffs;
  for (int i = 0; i < histA->GetNbinsX(); i++) {
    for (int j = 0; j < histA->GetNbinsY(); j++) {
      double binDiff = histA->GetBinContent(i, j) - histB->GetBinContent(i, j);
      if (abs(binDiff) > tolerance) {
        bins.push_back(std::array<int, 2>{i, j});
        diffs.push_back(binDiff);
      }
    }
  }
  cout << histA->GetName() << " - " << histB->GetName();
  if (bins.size() == 0) {
    cout << ": All bins match within tolerance " << tolerance << "\n";
    return true;
  } else {
    for (int i = 0; i < bins.size(); i++) {
      cout << " Bin (" << bins[i][0] << ", " << bins[i][1] << ") = " << diffs[i] << " ";
    }
    cout << "\n";
    return false;
  }
}

// Trivial closure test means we use the same data to build the Truth and Measured distributions as we use to build the response matrix
void createResponse() {
  // Truth
  double xmin = -6, xmax = 6;
  double ymin = -3, ymax = 3;
  int nbinsx = 12, nbinsy = 6;

  // Measured
  double amin = -4, amax = 4;
  double bmin = -1, bmax = 3;
  int nbinsa = 16, nbinsb = 8;

  TF1 *fxTruth = new TF1("fxTruth", "gaus(0)", xmin, xmax);
  fxTruth->SetParameters(1., 0., 1.);
  TF1 *fyTruth = new TF1("fyTruth","TMath::Landau(x,[0],[1],0)", ymin, ymax);
  fyTruth->SetParameters(0.2, 0.3);

  TF1 *resx = new TF1("resx", "gaus(0)", -10, 10);
  resx->SetParameters(1., 0., 1.);
  TF1 *resy = new TF1("resy", "gaus(0)", -10, 10);
  resy->SetParameters(1., 0., 1.5);

  TH2F* hTruth = new TH2F("hTruth", "hTruth", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  TH2F* hMeasured = new TH2F("hMeasured", "hMeasured", nbinsa, amin, amax, nbinsb, bmin, bmax);
  TH2F* hMissTrivial = (TH2F*)hTruth->Clone("hMissTrivial");
  TH2F* hFakeTrivial = (TH2F*)hMeasured->Clone("hFakeTrivial");
  TH2F* hMissIndependent = (TH2F*)hTruth->Clone("hMissIndependent");
  TH2F* hFakeIndependent = (TH2F*)hMeasured->Clone("hFakeIndependent");

  RooUnfoldResponse *responseTrivial = new RooUnfoldResponse("responseTrivial", "responseTrivial");
  responseTrivial->Setup(hMeasured, hTruth);
  RooUnfoldResponse *responseIndependent = new RooUnfoldResponse("responseIndependent", "responseIndependent");
  responseIndependent->Setup(hMeasured, hTruth);

  int ntoys = 1e6;
  // Fill truth and measured distributions and response for trivial closure
  for (int i = 0; i < ntoys; i++) {
    double x = fxTruth->GetRandom();
    double a = x * resx->GetRandom();
    double y = fyTruth->GetRandom();
    double b = y * resy->GetRandom();

    bool isTruth = (x > xmin && x < xmax) * (y > ymin && y < ymax);
    bool isMeasured = (a > amin && a < amax) * (b > bmin && b < bmax);

    if ( isTruth ) {
      hTruth->Fill(x, y);
    }
    if ( isMeasured ) {
      hMeasured->Fill(a, b);
    }

    // Trivial Closure
    if ( isTruth && isMeasured ) {
      responseTrivial->Fill(a, b, x, y);
    }
    else if ( isTruth && !isMeasured ) {
      responseTrivial->Miss(x, y);
      hMissTrivial->Fill(x, y);
    }
    else if ( !isTruth && isMeasured ) {
      responseTrivial->Fake(a, b);
      hFakeTrivial->Fill(a, b);
    }
  }
  // Fill response for statistically independent closure test
  for (int i = 0; i < ntoys; i++) {
    double x = fxTruth->GetRandom();
    double a = x * resx->GetRandom();
    double y = fyTruth->GetRandom();
    double b = y * resy->GetRandom();

    bool isTruth = (x > xmin && x < xmax) * (y > ymin && y < ymax);
    bool isMeasured = (a > amin && a < amax) * (b > bmin && b < bmax);

    if ( isTruth && isMeasured ) {
      responseIndependent->Fill(a, b, x, y);
    }
    else if ( isTruth && !isMeasured ) {
      responseIndependent->Miss(x, y);
      hMissIndependent->Fill(x, y);
    }
    else if ( !isTruth && isMeasured ) {
      responseIndependent->Fake(a, b);
      hFakeIndependent->Fill(a, b);
    }
  }
  TFile* outFile = new TFile("toy2D.root", "RECREATE");
  hMeasured->Write();
  hTruth->Write();

  responseTrivial->Write();
  hMissTrivial->Write();
  hFakeTrivial->Write();

  responseIndependent->Write();
  hMissIndependent->Write();
  hFakeIndependent->Write();

  outFile->Write();
  outFile->Close();
}
void doUnfolding(bool trivialClosure, int nIter) {
  int doSmoothing = 0;
  string testType = trivialClosure ? "Trivial" : "Independent";
  TFile* inFile = TFile::Open("toy2D.root", "UPDATE");
  TH2F* hMeasured = (TH2F*)inFile->Get("hMeasured");
  TH2F* hTruth = (TH2F*)inFile->Get("hTruth");

  TH2F* hFake = (TH2F*)inFile->Get(("hFake" + testType).c_str());
  RooUnfoldResponse* response = (RooUnfoldResponse*)inFile->Get(("response" + testType).c_str());
  RooUnfoldBayes ruBayes(response, hMeasured, nIter, doSmoothing, "unf", "unf");
  TH2F* hUnfolded = (TH2F*)ruBayes.Hreco(RooUnfold::kCovToy);
  hUnfolded->SetName(TString::Format("hUnfolded%sIter%d", testType.c_str(), nIter).Data());
  TH2F* hRefolded = (TH2F*)response->ApplyToTruth(hUnfolded, TString::Format("hRefolded%sIter%d", testType.c_str(), nIter).Data());
  hRefolded->Add(hFake);

  hUnfolded->Write();
  hRefolded->Write();
  inFile->Close();
}
void plotClosureTest(bool trivialClosure, int nIter, bool saveFigs = false) {
  TFile* inFile = TFile::Open("toy2D.root", "READ");

  string testType = trivialClosure ? "Trivial" : "Independent";
  RooUnfoldResponse* response = (RooUnfoldResponse*)inFile->Get(("response" + testType).c_str());

  // Measured distribution
  TH2F* hMeasured  = (TH2F*)inFile->Get("hMeasured"); hMeasured->SetStats(0);
  TH1F* hMeasuredX = (TH1F*)hMeasured->ProjectionX("hMeasuredX");
  TH1F* hMeasuredY = (TH1F*)hMeasured->ProjectionY("hMeasuredY");
  hMeasuredX->SetLineColor(kBlue); hMeasuredX->SetMarkerColor(kBlue);
  hMeasuredY->SetLineColor(kBlue); hMeasuredY->SetMarkerColor(kBlue);

  // Truth distribution
  TH2F* hTruth  = (TH2F*)inFile->Get("hTruth"); hTruth->SetStats(0);
  TH1F* hTruthX = (TH1F*)hTruth->ProjectionX("hTruthX");
  TH1F* hTruthY = (TH1F*)hTruth->ProjectionY("hTruthY");
  hTruthX->SetLineColor(kBlue); hTruthX->SetMarkerColor(kBlue);
  hTruthY->SetLineColor(kBlue); hTruthY->SetMarkerColor(kBlue);

  // Unfolded distribution
  TH2F* hUnfolded  = (TH2F*)inFile->Get(TString::Format("hUnfolded%sIter%d", testType.c_str(), nIter).Data());
  TH1F* hUnfoldedX = (TH1F*)hUnfolded->ProjectionX(TString::Format("hUnfoldedX%sIter%d", testType.c_str(), nIter).Data());
  TH1F* hUnfoldedY = (TH1F*)hUnfolded->ProjectionY(TString::Format("hUnfoldedY%sIter%d", testType.c_str(), nIter).Data());
  hUnfoldedX->SetLineColor(kRed); hUnfoldedX->SetMarkerColor(kRed);
  hUnfoldedX->SetMarkerStyle(4);
  hUnfoldedY->SetLineColor(kRed); hUnfoldedY->SetMarkerColor(kRed);
  hUnfoldedY->SetMarkerStyle(4);

  // Refolded distribution
  TH2F* hRefolded  = (TH2F*)inFile->Get(TString::Format("hRefolded%sIter%d", testType.c_str(), nIter).Data());
  TH1F* hRefoldedX = (TH1F*)hRefolded->ProjectionX(TString::Format("hRefoldedX%sIter%d", testType.c_str(), nIter).Data());
  TH1F* hRefoldedY = (TH1F*)hRefolded->ProjectionY(TString::Format("hRefoldedY%sIter%d", testType.c_str(), nIter).Data());
  hRefoldedX->SetLineColor(kRed); hRefoldedX->SetMarkerColor(kRed);
  hRefoldedX->SetMarkerStyle(4);
  hRefoldedY->SetLineColor(kRed); hRefoldedY->SetMarkerColor(kRed);
  hRefoldedY->SetMarkerStyle(4);

  // What the response considers measured
  TH2F* respMeasured  = (TH2F*)response->Hmeasured();
  TH1F* respMeasuredX = (TH1F*)respMeasured->ProjectionX("respMeasuredX");
  TH1F* respMeasuredY = (TH1F*)respMeasured->ProjectionY("respMeasuredY");
  respMeasuredX->SetMarkerColor(kRed); respMeasuredX->SetLineColor(kRed);
  respMeasuredX->SetMarkerStyle(4);
  respMeasuredY->SetMarkerColor(kRed); respMeasuredY->SetLineColor(kRed);
  respMeasuredY->SetMarkerStyle(4);

  // What the response considers truth
  TH2F* respTruth  = (TH2F*)response->Htruth();
  TH1F* respTruthX = (TH1F*)respTruth->ProjectionX("respTruthX");
  TH1F* respTruthY = (TH1F*)respTruth->ProjectionY("respTruthY");
  respTruthX->SetMarkerColor(kRed); respTruthX->SetLineColor(kRed);
  respTruthX->SetMarkerStyle(4);
  respTruthY->SetMarkerColor(kRed); respTruthY->SetLineColor(kRed);
  respTruthY->SetMarkerStyle(4);

  // Response applied to truth distribution
  TH2F* respAppliedtoTruth  = (TH2F*)response->ApplyToTruth();
  TH1F* respAppliedtoTruthX = (TH1F*)respAppliedtoTruth->ProjectionX("respAppliedtoTruthX");
  TH1F* respAppliedtoTruthY = (TH1F*)respAppliedtoTruth->ProjectionY("respAppliedtoTruthY");
  respAppliedtoTruthX->SetMarkerColor(kRed); respAppliedtoTruthX->SetLineColor(kRed);
  respAppliedtoTruthX->SetMarkerStyle(4);
  respAppliedtoTruthY->SetMarkerColor(kRed); respAppliedtoTruthY->SetLineColor(kRed);
  respAppliedtoTruthY->SetMarkerStyle(4);

  // Sanity check: do measured/truth distributions match what the response thinks they are?
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  c1->Divide(2,2);
  c1->cd(1);
  hTruthX->Draw("hist");
  respTruthX->Draw("same p");
  TLegend* l11 = new TLegend(0.15,0.7,0.5,0.9);
  l11->AddEntry(hTruthX, "Truth");
  l11->AddEntry(respTruthX, "Response Truth");
  l11->Draw("same");

  c1->cd(2);
  hTruthY->Draw("hist");
  respTruthY->Draw("same p");
  TLegend* l12 = new TLegend(0.15,0.7,0.5,0.9);
  l12->AddEntry(hTruthY, "Truth");
  l12->AddEntry(respTruthY, "Response Truth");
  l12->Draw("same");

  c1->cd(3);
  hMeasuredX->Draw("hist");
  respMeasuredX->Draw("same p");
  TLegend* l13 = new TLegend(0.15,0.7,0.5,0.9);
  l13->AddEntry(hMeasuredX, "Measured");
  l13->AddEntry(respMeasuredX, "Response Measured");
  l13->Draw("same");

  c1->cd(4);
  hMeasuredY->Draw("hist");
  respMeasuredY->Draw("same p");
  TLegend* l14 = new TLegend(0.15,0.7,0.5,0.9);
  l14->AddEntry(hMeasuredY, "Measured");
  l14->AddEntry(respMeasuredY, "Response Measured");
  l14->Draw("same");

  if (saveFigs)
    c1->SaveAs(TString::Format("checkResponse%sIter%d.pdf", testType.c_str(), nIter).Data());

  // Compare truth/unfolded and measured/refolded:
  // Should match exactly for trivial closure, should be close for independent closure
  TCanvas* c2 = new TCanvas("c2", "c2", 800, 800);
  c2->Divide(2,2);
  c2->cd(1);
  hTruthX->Draw("hist");
  hUnfoldedX->Draw("same p");
  TLegend* l21 = new TLegend(0.15,0.7,0.5,0.9);
  l21->AddEntry(hTruthX, "Truth");
  l21->AddEntry(hUnfoldedX, "Unfolded");
  l21->Draw("same");

  c2->cd(2);
  hTruthY->Draw("hist");
  hUnfoldedY->Draw("same p");
  TLegend* l22 = new TLegend(0.15,0.7,0.5,0.9);
  l22->AddEntry(hTruthY, "Truth");
  l22->AddEntry(hUnfoldedY, "Unfolded");
  l22->Draw("same");

  c2->cd(3);
  hMeasuredX->Draw("hist");
  hRefoldedX->Draw("same p");
  TLegend* l23 = new TLegend(0.15,0.7,0.5,0.9);
  l23->AddEntry(hMeasuredX, "Measured");
  l23->AddEntry(hRefoldedX, "Refolded");
  l23->Draw("same");

  c2->cd(4);
  hMeasuredY->Draw("hist");
  hRefoldedY->Draw("same p");
  TLegend* l24 = new TLegend(0.15,0.7,0.5,0.9);
  l24->AddEntry(hMeasuredY, "Measured");
  l24->AddEntry(hRefoldedY, "Refolded");
  l24->Draw("same");

  if (saveFigs)
    c2->SaveAs(TString::Format("checkClosure%sIter%d.pdf", testType.c_str(), nIter).Data());

  // Compare measured and refolded with response applied to truth
  // Shows the effect of (not) adding fakes
  TCanvas* c3 = new TCanvas("c3", "c3", 800, 800);
  c3->Divide(2,2);
  c3->cd(1);
  hMeasuredX->Draw("hist");
  hRefoldedX->Draw("same p");
  TLegend* l31 = new TLegend(0.15,0.7,0.5,0.9);
  l31->AddEntry(hMeasuredX, "Measured");
  l31->AddEntry(hRefoldedX, "Refolded");
  l31->Draw("same");

  c3->cd(2);
  hMeasuredY->Draw("hist");
  hRefoldedY->Draw("same p");
  TLegend* l32 = new TLegend(0.15,0.7,0.5,0.9);
  l32->AddEntry(hMeasuredY, "Measured");
  l32->AddEntry(hRefoldedY, "Refolded");
  l32->Draw("same");

  c3->cd(3);
  hMeasuredX->Draw("hist");
  respAppliedtoTruthX->Draw("same p");
  TLegend* l33 = new TLegend(0.15,0.7,0.5,0.9);
  l33->AddEntry(hMeasuredX, "Measured");
  l33->AddEntry(respAppliedtoTruthX, "Response Applied to Truth");
  l33->Draw("same");

  c3->cd(4);
  hMeasuredY->Draw("hist");
  respAppliedtoTruthY->Draw("same p");
  TLegend* l34 = new TLegend(0.15,0.7,0.5,0.9);
  l34->AddEntry(hMeasuredY, "Measured");
  l34->AddEntry(respAppliedtoTruthY, "Response Applied to Truth");
  l34->Draw("same");

  if (saveFigs)
    c3->SaveAs(TString::Format("checkApplyToTruth%sIter%d.pdf", testType.c_str(), nIter).Data());
}

void toymodel() {
  int nIter = 3;
  const bool doTrivialClosure = true;
  createResponse();
  doUnfolding(doTrivialClosure, nIter);
  plotClosureTest(doTrivialClosure, nIter);

  doUnfolding(!doTrivialClosure, nIter);
  plotClosureTest(!doTrivialClosure, nIter);
}

// -----------------------------------------------------------------
// Official RooUnfold tutorial
// -----------------------------------------------------------------

const Double_t cutdummy= -99999.0;
Double_t smear (Double_t xt){
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}

void tutorial()
{
  RooUnfoldResponse response (40, -10.0, 10.0);
  auto *f0 = new TH1F("f0","f0",40,-10,10);
  auto *g0 = new TH1F("g0","g0",40,-10,10);

  for (Int_t i= 0; i<100000; i++) {
    Double_t xt= gRandom->BreitWigner (0.3, 2.5);
    f0->Fill(xt);
    Double_t x= smear (xt);
    if (x!=cutdummy){
      g0->Fill(x);
      response.Fill (x, xt);
    }
    else{
        response.Miss (xt);
    }
  }

  auto* c = new TCanvas();
  f0->SetStats(0);
  f0->SetTitle("");
  f0->SetFillColor(7);
  f0->Draw();
  g0->SetFillColor(42);
  g0->Draw("same");
  auto* leg = new TLegend(.55,0.7,.9,.9);
  leg->AddEntry(f0,"True Distribution");
  leg->AddEntry(g0,"Predicted Measured");
  leg->Draw();
  c->Draw();
  c->SaveAs("true-response.png");

  auto* R = response.HresponseNoOverflow();
  auto* c1 = new TCanvas();
  R->SetStats(0);
  R->Draw("colz");
  c1->Draw();
  c1->SaveAs("response.png");

  TH1D* hTrue= new TH1D ("true", "Test Truth",    40, -10.0, 10.0);
  TH1D* hMeas= new TH1D ("meas", "Test Measured", 40, -10.0, 10.0);
  // Test with a Gaussian, mean 0 and width 2.
  for (Int_t i=0; i<10000; i++) {
    Double_t xt= gRandom->Gaus (0.0, 2.0), x= smear (xt);
    hTrue->Fill(xt);
    if (x!=cutdummy) hMeas->Fill(x);
  }
  RooUnfoldInvert   unfold (&response, hMeas);
  auto* hReco = unfold.Hreco();

  auto* c2 = new TCanvas();
  hReco->SetStats(0);
  hReco->SetLineColor(6);
  hReco->SetMarkerColor(6);
  hTrue->SetLineColor(2);
  hReco->Draw();
  hTrue->Draw("same");
  hMeas->Draw("same");
  auto* leg2 = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg2->AddEntry(hTrue, "True distribution","pl");
  leg2->AddEntry(hMeas, "Measured distribution", "pl");
  leg2->AddEntry(hReco, "Unfolded distribution");
  leg2->Draw();
  c2->Draw();

  RooUnfoldResponse double_response (40, -10.0, 10.0);
  auto *f1 = new TH1F("f1","f1",40,-10,10);
  auto *g1 = new TH1F("g1","g1",40,-10,10);
  for (Int_t i= 0; i<5000; i++) {
    Double_t xt1= gRandom->Gaus (2, 1.5);
    f1->Fill(xt1);
    Double_t x1= gRandom->Gaus (xt1, 1.);
    if (x1!=cutdummy){
      g1->Fill(x1);
      double_response.Fill (x1, xt1);
    }
    else{
      double_response.Miss (xt1);
    }
    Double_t xt2= gRandom->Gaus (-2, 1.5);
    f1->Fill(xt2);
    Double_t x2= gRandom->Gaus (xt2, 1.);
    if (x2!=cutdummy){
      g1->Fill(x2);
      double_response.Fill (x2, xt2);
    }
    else{
      double_response.Miss (xt2);
    }
  }

  auto* c3 = new TCanvas();
  f1->SetStats(0);
  f1->SetTitle("");
  f1->SetFillColor(7);
  f1->Draw();
  g1->SetFillColor(42);
  g1->Draw("same");
  auto* leg3 = new TLegend(.55,0.7,.9,.9);
  leg3->AddEntry(f1,"True Distribution");
  leg3->AddEntry(g1,"Predicted Measured");
  leg3->Draw();
  c3->Draw();
  c3->SaveAs("double-peak-response.png");

  RooUnfoldInvert   unfold_double (&double_response, hMeas);
  auto* hReco2 = unfold_double.Hreco();
  auto* c4 = new TCanvas();
  hReco2->SetStats(0);
  hReco2->SetLineColor(6);
  hReco2->SetMarkerColor(6);
  hTrue->SetLineColor(2);
  hReco2->Draw();
  hTrue->Draw("same");
  hMeas->Draw("same");
  auto* leg4 = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg4->AddEntry(hTrue, "True distribution","pl");
  leg4->AddEntry(hMeas, "Measured distribution", "pl");
  leg4->AddEntry(hReco2, "Unfolded distribution");
  leg4->Draw();
  c4->Draw();

  RooUnfoldInvert   unfold_bin (&response, hMeas);
  auto* hReco_bin = unfold_bin.Hreco();
  auto* c5 = new TCanvas();
  hReco_bin->SetStats(0);
  hTrue->SetLineColor(2);
  hReco_bin->Draw();
  hTrue->Draw("same");
  hMeas->Draw("same");
  auto* leg5 = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg5->AddEntry(hTrue, "True distribution","pl");
  leg5->AddEntry(hMeas, "Measured distribution", "pl");
  leg5->AddEntry(hReco_bin, "Unfolded distribution");
  leg5->Draw();
  c5->Draw();

  RooUnfoldBayes   unfold_bayes (&response, hMeas, 4);
  auto* hReco_bayes = unfold_bayes.Hreco();
  auto* c6 = new TCanvas();
  hReco_bayes->SetStats(0);
  hTrue->SetLineColor(2);
  hReco_bayes->Draw();
  hTrue->Draw("same");
  hMeas->Draw("same");
  auto* leg6 = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg6->AddEntry(hTrue, "True distribution","pl");
  leg6->AddEntry(hMeas, "Measured distribution", "pl");
  leg6->AddEntry(hReco_bayes, "Unfolded distribution");
  leg6->Draw();
  c6->Draw();
}
