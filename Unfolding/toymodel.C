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

void test()
{
  double xTmin = -6, xTmax = 6;
  double yTmin = 0., yTmax = 20;
  double xMmin = -4, xMmax = 4;
  double yMmin = 4., yMmax = 10;

  TF1 *fx = new TF1("fx", "gaus(0)", xTmin, xTmax);
  fx->SetParameters(1., 0., 1.);
  TF1 *fy = new TF1("fy", "TMath::Landau(x,[0],[1],0)", yTmin, yTmax);
  fy->SetParameters(0.2, 0.3);

  TF1 *fresx = new TF1("fresx", "gaus(0)", -10., 10.);
  fresx->SetParameters(1., 1., 0.3);
  TF1 *fresy = new TF1("fresy", "gaus(0)", -10., 10.);
  fresy->SetParameters(1., 1., 3.);

  int nBinsxT = 10;
  int nBinsyT = 10;
  int nBinsxM = nBinsxT;
  int nBinsyM = 8;
  TH2F *h2Truth =
    new TH2F("h2Truth", "h2Truth", nBinsxT, xTmin, xTmax, nBinsyT, yTmin, yTmax);
  TH2F *h2Measured =
    new TH2F("h2Measured", "h2Measured", nBinsxM, xMmin, xMmax, nBinsyM, yMmin, yMmax);
  TH2F *h2Prior =
    new TH2F("h2Prior", "h2Prior", nBinsxT, xTmin, xTmax, nBinsyT, yTmin, yTmax);
  TH2F *h2TInRM =
    new TH2F("h2TInRM", "h2TInRM", nBinsxT, xTmin, xTmax, nBinsyT, yTmin, yTmax);

  RooUnfoldResponse *fResponse = new RooUnfoldResponse("resp", "resp");
  fResponse->Setup(h2Measured,h2Truth);
  const int ntoys = 100000;

  for (int i = 0; i < ntoys; i++) {
    double xT = fx->GetRandom();
    double yT = fy->GetRandom();
    double xM = xT * fresx->GetRandom();
    double yM = yT * fresy->GetRandom();

    bool isTruth = false;
    isTruth = (xT >= xTmin) * (xT < xTmax) * (yT >= yTmin) * (yT < yTmax);
    bool isMeasured = false;
    isMeasured = (xM >= xMmin) * (xM < xMmax) * (yM >= yMmin) * (yM < yMmax);

    if (isTruth && isMeasured) {
      h2Truth->Fill(xT, yT);
      h2Measured->Fill(xM, yM);

      h2TInRM->Fill(xT, yT);
      fResponse->Fill(xM, yM, xT, yT);
    }
    else if (isTruth) {
      h2Truth->Fill(xT, yT);
      fResponse->Miss(xT, yT);
    }
    else if (isMeasured) {
      h2Measured->Fill(xM, yM);
      fResponse->Fake(xM, yM);
    }
  } // for ntoys

  int iterDef = 1;//3;
  int iterMin = 1;//3;
  int iterMax = 1;//3;
  int nIterTmp = iterMax - iterMin+1;
  const int nIter = nIterTmp;

  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
  RooUnfoldBayes unfold[nIter];
  TH2F *hReco[nIter];
  TH2F *hFolded[nIter];
  TMatrixD covmat[nIter];
  TH2D *hCovMat[nIter];
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {
    cout << "Unfolding loop " << iter << endl;
    unfold[iter-iterMin] = RooUnfoldBayes(fResponse, h2Measured, iter);
    hReco[iter-iterMin] = (TH2F*)unfold[iter-iterMin].Hreco(errorTreatment);
    hReco[iter-iterMin]->SetName(Form("hReco_Iter%d",iter));

    hFolded[iter-iterMin] = (TH2F*)fResponse->ApplyToTruth(hReco[iter-iterMin],Form("hFolded_Iter%d",iter));
    //hFolded[iter-iterMin]->SetName(Form("hFolded_Iter%d",iter));

    //Get covariance matrix
    // cout << "get covariance for iter " << iter << endl;
    // TH2D htmp(unfold[iter-iterMin].Ereco(errorTreatment));
    // hCovMat[iter-iterMin] = dynamic_cast<TH2D*>(htmp.Clone(Form("hCovMat_Iter%d",iter)));//TH2D::TH2D(covmat);
    // hCovMat[iter-iterMin]->SetName(Form("hCovMat_Iter%d",iter));
  } // Unfolding iterations loop

  cout << "After unfolding loop" << endl;

  TH1D *hRecoP[2];
  hRecoP[0] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin]->ProjectionX("hRecoP_x"));
  hRecoP[1] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin]->ProjectionY("hRecoP_y"));

  TH1D *hFoldedP[2];
  hFoldedP[0] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin]->ProjectionX("hFoldedP_x"));
  hFoldedP[1] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin]->ProjectionY("hFoldedP_y"));

  TH1D *hMeasured[2];
  hMeasured[0] = dynamic_cast<TH1D*>(h2Measured->ProjectionX("hMeasured_x"));
  hMeasured[1] = dynamic_cast<TH1D*>(h2Measured->ProjectionY("hMeasured_y"));

  TH1D *hTruth[2];
  hTruth[0] = dynamic_cast<TH1D*>(h2Truth->ProjectionX("hTruth_x"));
  hTruth[1] = dynamic_cast<TH1D*>(h2Truth->ProjectionY("hTruth_y"));

  for (int iDim = 0; iDim < 2; iDim++) {
    hRecoP[iDim]->SetMarkerColor(1);
    hRecoP[iDim]->SetLineColor(1);
    hRecoP[iDim]->SetMarkerStyle(2);
    // hRecoP[iDim]->SetLineStyle();
    hRecoP[iDim]->SetLineWidth(3);
    // hRecoP[iDim]->SetOptStat(0);

    hFoldedP[iDim]->SetMarkerColor(1);
    hFoldedP[iDim]->SetLineColor(1);
    hFoldedP[iDim]->SetMarkerStyle(2);
    // hFoldedP[iDim]->SetLineStyle();
    hFoldedP[iDim]->SetLineWidth(3);
    // hFoldedP[iDim]->SetOptStat(0);

    hMeasured[iDim]->SetMarkerColor(3);
    hMeasured[iDim]->SetLineColor(3);
    hMeasured[iDim]->SetMarkerStyle(5);
    // hMeasured[iDim]->SetLineStyle();
    hMeasured[iDim]->SetLineWidth(3);
    // hMeasured[iDim]->SetOptStat(0);

    hTruth[iDim]->SetMarkerColor(4);
    hTruth[iDim]->SetLineColor(4);
    hTruth[iDim]->SetMarkerStyle(25);
    // hTruth[iDim]->SetLineStyle();
    hTruth[iDim]->SetLineWidth(3);
    // hTruth[iDim]->SetOptStat(0);
  }

  TCanvas* myCanvas = new TCanvas("myCanvas", "myCanvas", 800, 800);
  myCanvas->Divide(2, 2);
  gStyle->SetOptStat(0);

  cout << "Plot 1" << endl;
  myCanvas->cd(1);
  gPad->SetLogy();
  TLegend *leg1 = new TLegend(0.15,0.6,0.5,0.9);
  leg1->AddEntry(hTruth[0], "truth");
  hTruth[0]->Draw();
  hTruth[0]->Draw("p same");
  hTruth[0]->Print();
  leg1->AddEntry(hRecoP[0], "unfolded");
  hRecoP[0]->Draw("same hist p");
  hRecoP[0]->Print();
  leg1->Draw("same");

  cout << "Plot 2" << endl;
  myCanvas->cd(2);
  gPad->SetLogy();
  TLegend *leg2 = new TLegend(0.15,0.6,0.5,0.9);
  leg2->AddEntry(hTruth[1], "truth");
  hTruth[1]->Draw();
  hTruth[1]->Draw("p same");
  hTruth[1]->Print();
  leg2->AddEntry(hRecoP[1], "unfolded");
  hRecoP[1]->Draw("same hist p");
  hRecoP[1]->Print();
  leg2->Draw("same");

  cout << "Plot 3" << endl;
  myCanvas->cd(3);
  gPad->SetLogy();
  TLegend *leg3 = new TLegend(0.15,0.6,0.5,0.9);
  leg3->AddEntry(hMeasured[0], "measured");
  hMeasured[0]->Draw();
  hMeasured[0]->Draw("p same");
  hMeasured[0]->Print();
  leg3->AddEntry(hFoldedP[0], "refolded");
  hFoldedP[0]->Draw("same hist p");
  hFoldedP[0]->Print();
  leg3->Draw("same");

  cout << "Plot 4" << endl;
  myCanvas->cd(4);
  gPad->SetLogy();
  TLegend *leg4 = new TLegend(0.15,0.6,0.5,0.9);
  leg4->AddEntry(hMeasured[0], "measured");
  leg4->AddEntry(hFoldedP[0], "refolded");
  hMeasured[1]->Draw();
  hMeasured[1]->Draw("p same");
  hMeasured[1]->Print();
  hFoldedP[1]->Draw("same hist p");
  hFoldedP[1]->Print();
  leg4->Draw("same");

  cout << "Saving" << endl;
  myCanvas->SaveAs("./spectra.pdf");

  // ------------------------------------------------------------------------------------

  TH1D* hTruthOverUnfolded[2];
  hTruthOverUnfolded[0] = dynamic_cast<TH1D*>(hTruth[0]->Clone("hTruthOverUnfolded_x"));
  hTruthOverUnfolded[0]->Divide(hRecoP[0]);
  hTruthOverUnfolded[1] = dynamic_cast<TH1D*>(hTruth[1]->Clone("hTruthOverUnfolded_y"));
  hTruthOverUnfolded[1]->Divide(hRecoP[1]);

  TH1D* hMeasuredOverRefolded[2];
  hMeasuredOverRefolded[0] = dynamic_cast<TH1D*>(hMeasured[0]->Clone("hMeasuredOverRefolded_x"));
  hMeasuredOverRefolded[0]->Divide(hFoldedP[0]);
  hMeasuredOverRefolded[1] = dynamic_cast<TH1D*>(hMeasured[1]->Clone("hMeasuredOverRefolded_y"));
  hMeasuredOverRefolded[1]->Divide(hFoldedP[1]);

  TCanvas* ratioCanvas = new TCanvas("ratioCanvas", "ratioCanvas", 800, 800);
  ratioCanvas->Divide(2, 2);
  gStyle->SetOptStat(0);

  cout << "Plot 1" << endl;
  ratioCanvas->cd(1);
  TLegend *rLeg1 = new TLegend(0.15,0.6,0.5,0.9);
  rLeg1->AddEntry(hTruthOverUnfolded[0], "truth/unfolded");
  hTruthOverUnfolded[0]->Draw();
  hTruthOverUnfolded[0]->Draw("p same");
  hTruthOverUnfolded[0]->Print("all");
  rLeg1->Draw("same");

  cout << "Plot 2" << endl;
  ratioCanvas->cd(2);
  TLegend *rLeg2 = new TLegend(0.15,0.6,0.5,0.9);
  rLeg2->AddEntry(hTruthOverUnfolded[1], "truth/unfolded");
  hTruthOverUnfolded[1]->Draw();
  hTruthOverUnfolded[1]->Draw("p same");
  hTruthOverUnfolded[1]->Print("all");
  rLeg2->Draw("same");

  cout << "Plot 3" << endl;
  ratioCanvas->cd(3);
  TLegend *rLeg3 = new TLegend(0.15,0.6,0.5,0.9);
  rLeg3->AddEntry(hMeasuredOverRefolded[0], "measured/refolded");
  hMeasuredOverRefolded[0]->Draw();
  hMeasuredOverRefolded[0]->Draw("p same");
  hMeasuredOverRefolded[0]->Print("all");
  rLeg3->Draw("same");

  cout << "Plot 4" << endl;
  ratioCanvas->cd(4);
  TLegend *rLeg4 = new TLegend(0.15,0.6,0.5,0.9);
  rLeg4->AddEntry(hMeasuredOverRefolded[1], "measured/refolded");
  hMeasuredOverRefolded[1]->Draw();
  hMeasuredOverRefolded[1]->Draw("p same");
  hMeasuredOverRefolded[1]->Print("all");
  rLeg4->Draw("same");

  cout << "Saving" << endl;
  ratioCanvas->SaveAs("./ratio.pdf");
}

void CompareHists(TH1* histA, TH1* histB)
{
  cout << histA->GetName() << " - " << histB->GetName() << endl;
  for (int i = 0; i < histA->GetNbinsX(); i++)
  {
    double binDiff = histA->GetBinContent(i) - histB->GetBinContent(i);
    if ( abs(binDiff) > 0.) cout << "Bin " << i << " = " << binDiff << " ";
  }
  cout << endl;
}
void CompareHists(TH2* histA, TH2* histB)
{
  cout << histA->GetName() << " - " << histB->GetName() << endl;
  for (int i = 0; i < histA->GetNbinsX(); i++)
  {
    for (int j = 0; j < histA->GetNbinsY(); j++) {
      double binDiff = histA->GetBinContent(i, j) - histB->GetBinContent(i, j);
      if ( abs(binDiff) > 0.) cout << "Bin (" << i << ", " << j << ") = " << binDiff << " ";
    }
  }
  cout << endl;
}

void toyResponse1D(bool trivialClosure = true)
{
  double xmin = -2, xmax = 2;
  double ymin = -1.5, ymax = 1.5;
  int nbinsx = 20, nbinsy = 40;
  TF1 *fx = new TF1("fx", "gaus(0)", xmin, xmax);
  fx->SetParameters(1., 0., 1.);
  TF1 *gx = new TF1("gx", "gaus(0)", xmin, xmax);
  gx->SetParameters(1., 1., .3);

  TH1F* hT = new TH1F("hT", "hT", nbinsx, xmin, xmax);
  hT->Sumw2();
  TH1F* hM = new TH1F("hM", "hM", nbinsy, ymin, ymax);
  hM->Sumw2();
  TH1F* hMiss = new TH1F("hMiss", "hMiss", nbinsx, xmin, xmax);
  hMiss->Sumw2();
  RooUnfoldResponse *response = new RooUnfoldResponse("response", "response");
  response->Setup(hM, hT);

  for (int i = 0; i < 10000; i++)
  {
    double x = fx->GetRandom();
    double y = x * gx->GetRandom();
    bool xGood = (x >= xmin && x < xmax);
    bool yGood = (y >= ymin && y < ymax);

    if ( xGood ) {
      hT->Fill(x);
    }
    if ( yGood ) {
      hM->Fill(y);
    }
    if (trivialClosure) {
      if ( xGood && yGood ) {
        response->Fill(y, x);
      }
      else if ( xGood && !yGood ) {
        response->Miss(x);
        hMiss->Fill(x);
      }
      else if ( !xGood && yGood ) {
        response->Fake(y);
      }
    }
  }
  if (!trivialClosure) {
    for (int i = 0; i < 10000; i++)
    {
      double x = fx->GetRandom();
      double y = x * gx->GetRandom();
      bool xGood = (x >= xmin && x < xmax);
      bool yGood = (y >= ymin && y < ymax);
      if ( xGood && yGood ) {
        response->Fill(y, x);
      }
      else if ( xGood && !yGood ) {
        response->Miss(x);
        hMiss->Fill(x);
      }
      else if ( !xGood && yGood ) {
        response->Fake(y);
      }
    }
  }

  TFile* outFile = new TFile("toy1D.root", "RECREATE");
  response->Write();
  hM->Write();
  hMiss->Write();
  hT->Write();
  outFile->Write();
  outFile->Close();
}

void toyResponse2D(bool trivialClosure = true)
{
  double xmin = -6, xmax = 6;
  double ymin = -3, ymax = 3;
  int nbinsx = 12, nbinsy = 6;

  double amin = -4, amax = 4;
  double bmin = -1, bmax = 3;
  int nbinsa = 16, nbinsb = 8;

  TF1 *fx = new TF1("fx", "gaus(0)", xmin, xmax);
  fx->SetParameters(1., 0., 1.);
  TF1 *fy = new TF1("fy","TMath::Landau(x,[0],[1],0)", ymin, ymax);
  fy->SetParameters(0.2, 0.3);

  TF1 *gx = new TF1("gx", "gaus(0)", -10, 10);
  gx->SetParameters(1., 0., 1.);
  TF1 *gy = new TF1("gy", "gaus(0)", -10, 10);
  gy->SetParameters(1., 0., 1.5);

  TH2F* hT = new TH2F("hT", "hT", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  TH2F* hMiss = new TH2F("hMiss", "hMiss", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  TH2F* hM = new TH2F("hM", "hM", nbinsa, amin, amax, nbinsb, bmin, bmax);
  TH2F* hFake = new TH2F("hFake", "hFake", nbinsa, amin, amax, nbinsb, bmin, bmax);
  RooUnfoldResponse *response = new RooUnfoldResponse("response", "response");
  response->Setup(hM, hT);

  int ntoys = 1e6;
  for (int i = 0; i < ntoys; i++)
  {
    double x = fx->GetRandom();
    double a = x * gx->GetRandom();
    double y = fy->GetRandom();
    double b = y * gy->GetRandom();

    bool isTruth = (x > xmin && x < xmax) * (y > ymin && y < ymax);
    bool isMeasured = (a > amin && a < amax) * (b > bmin && b < bmax);

    if ( isTruth ) {
      hT->Fill(x, y);
    }
    if ( isMeasured ) {
      hM->Fill(a, b);
    }

    if (trivialClosure) {
      if ( isTruth && isMeasured ) {
        response->Fill(a, b, x, y);
      }
      else if ( isTruth && !isMeasured ) {
        response->Miss(x, y);
        hMiss->Fill(x, y);
      }
      else if ( !isTruth && isMeasured ) {
        response->Fake(a, b);
        hFake->Fill(a, b);
      }
    }
  }
  if (!trivialClosure)
  {
    for (int i = 0; i < ntoys; i++)
    {
      double x = fx->GetRandom();
      double a = x * gx->GetRandom();
      double y = fy->GetRandom();
      double b = y * gy->GetRandom();

      bool isTruth = (x > xmin && x < xmax) * (y > ymin && y < ymax);
      bool isMeasured = (a > amin && a < amax) * (b > bmin && b < bmax);

      if ( xGood && yGood ) {
        response->Fill(y, x);
      }
      else if ( xGood && !yGood ) {
        response->Miss(x);
        hMiss->Fill(x);
      }
      else if ( !xGood && yGood ) {
        response->Fake(y);
      }
    }
  }
  TFile* outFile = new TFile("toy2D.root", "RECREATE");
  response->Write();
  hM->Write();
  hT->Write();
  hMiss->Write();
  hFake->Write();
  outFile->Write();
  // outFile->Close();
}

void toyCompare1D(string inFileName = "toy1D.root")
{
  TFile* inFile = TFile::Open(inFileName.c_str(), "READ");
  RooUnfoldResponse* response = (RooUnfoldResponse*)inFile->Get("response");
  TH1F* hM = (TH1F*)inFile->Get("hM");
  TH1F* hT = (TH1F*)inFile->Get("hT");

  TH1F* rM = (TH1F*)response->Hmeasured();
  rM->SetName("response-measured");
  TH1F* rT = (TH1F*)response->Htruth();
  rT->SetName("response-truth");
  TH1F* rA = (TH1F*)response->ApplyToTruth();
  rA->SetName("response-applied");

  hM->Print("all");
  rM->Print("all");
  rA->Print("all");

  hT->Print("all");
  rT->Print("all");
}

void toyCompare2D(string inFileName = "toy2D.root")
{
  TFile* inFile = TFile::Open(inFileName.c_str(), "READ");

  RooUnfoldResponse* response = (RooUnfoldResponse*)inFile->Get("response");

  TH2F* hM = (TH2F*)inFile->Get("hM");
  TH1F* hMx = (TH1F*)hM->ProjectionX("hMx");
  TH1F* hMy = (TH1F*)hM->ProjectionY("hMy");

  TH2F* hT = (TH2F*)inFile->Get("hT");
  TH1F* hTx = (TH1F*)hT->ProjectionX("hTx");
  TH1F* hTy = (TH1F*)hT->ProjectionY("hTy");

  TH2F* rM = (TH2F*)response->Hmeasured();
  TH1F* rMx = (TH1F*)rM->ProjectionX("rMx");
  TH1F* rMy = (TH1F*)rM->ProjectionY("rMy");

  TH2F* rT = (TH2F*)response->Htruth();
  TH1F* rTx = (TH1F*)rT->ProjectionX("rTx");
  TH1F* rTy = (TH1F*)rT->ProjectionY("rTy");

  TH2F* rA = (TH2F*)response->ApplyToTruth();
  TH1F* rAx = (TH1F*)rA->ProjectionX("rAx");
  TH1F* rAy = (TH1F*)rA->ProjectionY("rAy");

  CompareHists(hMx, rMx);
  CompareHists(hMy, rMy);
  CompareHists(hTx, rTx);
  CompareHists(hTy, rTy);
  CompareHists(hMx, rAx);
  CompareHists(hMy, rAy);
}

void toyClosureTest1D(int nIter = 3)
{
  int doSmoothing = 0;
  TFile* inFile = TFile::Open("toy1D.root", "READ");
  RooUnfoldResponse* response = (RooUnfoldResponse*)inFile->Get("response");
  TH1F* hM = (TH1F*)inFile->Get("hM"); // measured test
  TH1F* hT = (TH1F*)inFile->Get("hT"); // truth test

  RooUnfoldBayes rub(response, hM, nIter, doSmoothing, "unf", "unf");
  TH1F* unfolded = (TH1F*)rub.Hreco(RooUnfold::kCovToy);
  unfolded->SetName("unfolded");
  TH1F* refolded = (TH1F*)response->ApplyToTruth(unfolded, "refolded");

  // cout << "truth - unfolded" << endl;
  // for (int i = 0; i < hT->GetNbinsX(); i++)
  // {
  //   double binDiff = hT->GetBinContent(i) - unfolded->GetBinContent(i);
  //   double binErr = hT->GetBinError(i);
  //   if ( abs(binDiff) > binErr) cout << "Bin " << i << " = " << binDiff << " ";
  // }
  // cout << endl;

  // cout << "measured - refolded" << endl;
  // for (int i = 0; i < hM->GetNbinsX(); i++)
  // {
  //   double binDiff = hM->GetBinContent(i) - refolded->GetBinContent(i);
  //   double binErr = hM->GetBinError(i);
  //   if ( abs(binDiff) > binErr) cout << "Bin " << i << " = " << binDiff << " ";
  // }
  // cout << endl;

  TCanvas* myCanvas = new TCanvas("canvas", "canvas", 800, 400);
  myCanvas->Divide(2,1);
  myCanvas->cd(1);

  int lowBin = hT->GetXaxis()->FindBin(-0.5);
  int highBin = hT->GetXaxis()->FindBin(0.5);

  TH1F* truthOverUnfolded = (TH1F*)hT->Clone("truthOverUnfolded");
  truthOverUnfolded->Divide(unfolded);
  truthOverUnfolded->GetXaxis()->SetRange(lowBin, highBin);
  truthOverUnfolded->SetLineWidth(3);
  truthOverUnfolded->SetLineColor(1);
  truthOverUnfolded->SetMarkerColor(1);
  truthOverUnfolded->SetMarkerStyle(4);
  truthOverUnfolded->Draw();
  TLine* truthLine = new TLine(truthOverUnfolded->GetXaxis()->GetBinLowEdge(lowBin), 1., truthOverUnfolded->GetXaxis()->GetBinLowEdge(highBin + 1), 1.);
  truthLine->Draw("same");

  myCanvas->cd(2);
  lowBin = hM->GetXaxis()->FindBin(-0.5);
  highBin = hM->GetXaxis()->FindBin(0.5);

  TH1F* measuredOverRefolded = (TH1F*)hM->Clone("measuredOverRefolded");
  measuredOverRefolded->Divide(refolded);
  measuredOverRefolded->GetXaxis()->SetRange(lowBin, highBin);
  measuredOverRefolded->SetLineWidth(3);
  measuredOverRefolded->SetLineColor(2);
  measuredOverRefolded->SetMarkerColor(2);
  measuredOverRefolded->SetMarkerStyle(27);
  measuredOverRefolded->SetMarkerSize(2);
  measuredOverRefolded->Draw("p");
  TLine* measuredLine = new TLine(measuredOverRefolded->GetXaxis()->GetBinLowEdge(lowBin), 1., measuredOverRefolded->GetXaxis()->GetBinLowEdge(highBin + 1), 1.);
  measuredLine->Draw("same");

  myCanvas->SaveAs(TString::Format("./closureTest1D-%d.pdf", nIter).Data());
}
void toyClosureTest2D(int nIter = 3)
{
  int doSmoothing = 0;
  TFile* inFile = TFile::Open("toy2D.root", "READ");
  RooUnfoldResponse* response = (RooUnfoldResponse*)inFile->Get("response");
  TH2F* hM = (TH2F*)inFile->Get("hM"); // measured test
  // hM->SetName("measured");
  TH2F* hT = (TH2F*)inFile->Get("hT"); // truth test
  // hT->SetName("truth");

  RooUnfoldBayes rub(response, hM, nIter, doSmoothing, "unf", "unf");
  TH2F* unfolded = (TH2F*)rub.Hreco(RooUnfold::kCovToy);
  unfolded->SetName("unfolded");
  TH2F* refolded = (TH2F*)response->ApplyToTruth(unfolded, "refolded");
  TH2F* hFake = (TH2F*)inFile->Get("hFake");
  refolded->Add(hFake);

  cout << "truth - unfolded" << endl;
  for (int i = 0; i < hT->GetNbinsX(); i++)
  {
    for (int j = 0; j < hT->GetNbinsY(); j++) {
      double binDiff = hT->GetBinContent(i, j) - unfolded->GetBinContent(i, j);
      if ( abs(binDiff) > 0.) cout << "Bin (" << i << ", " << j << ") = " << binDiff << " ";
    }
  }
  cout << endl;

  cout << "measured - refolded" << endl;
  for (int i = 0; i < hM->GetNbinsX(); i++)
  {
    for (int j = 0; j < hM->GetNbinsY(); j++) {
      double binDiff = hM->GetBinContent(i, j) - refolded->GetBinContent(i, j);
      if ( abs(binDiff) > 0.) cout << "Bin (" << i << ", " << j << ") = " << binDiff << " ";
    }
  }
  cout << endl;

  // TH1F* truthOverUnfolded = (TH1F*)hT->Clone("truthOverUnfolded");
  // truthOverUnfolded->Scale(2);
  // truthOverUnfolded->Divide(unfolded);
  // truthOverUnfolded->SetLineWidth(3);
  // truthOverUnfolded->SetLineColor(1);
  // truthOverUnfolded->SetMarkerColor(1);
  // truthOverUnfolded->SetMarkerStyle(4);

  // TH1F* measuredOverRefolded = (TH1F*)hM->Clone("measuredOverRefolded");
  // measuredOverRefolded->Divide(refolded);
  // measuredOverRefolded->SetLineWidth(3);
  // measuredOverRefolded->SetLineColor(2);
  // measuredOverRefolded->SetMarkerColor(2);
  // measuredOverRefolded->SetMarkerStyle(27);
  // measuredOverRefolded->SetMarkerSize(2);
}

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

  // -----------------------------------------------------------------

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
