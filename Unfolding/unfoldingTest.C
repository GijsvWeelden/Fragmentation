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


void unfoldingTest() {

  //RooUnfold library: set path to RooUnfold as environment variable
  //https://gitlab.cern.ch/RooUnfold/RooUnfold
  //gSystem->Load("/Users/mverweij/soft/unfold/RooUnfold/build/libRooUnfold.dylib");

  gStyle->SetOptStat(0000);
  gStyle->SetPalette(56);

  double minx = -6.;
  double maxx = 6.;
  TF1 *f1 = new TF1("f1","gaus(0)",minx,maxx);
  f1->SetParameters(1.,0.,1.);

  double miny = 0.;
  double maxy = 20.;
  TF1 *f2 = new TF1("f2","TMath::Landau(x,[0],[1],0)",miny,maxy);
  f2->SetParameters(0.2,0.3);

  TF1 *fresx = new TF1("fresx","gaus(0)",-10.,10.);
  fresx->SetParameters(1.,1.,0.3);

  TF1 *fresy = new TF1("fresy","gaus(0)",-10.,10.);
  fresy->SetParameters(1.,1.,3.);

  TCanvas *c1 = new TCanvas("c1","c1",400,400);
  f1->Draw();
  c1->SaveAs("f1.png");

  TCanvas *c2 = new TCanvas("c2","c2",400,400);
  f2->Draw();
  c2->SaveAs("f2.png");

  TCanvas *cres = new TCanvas("cres","cres",400,400);
  fresx->Draw();
  fresy->SetLineColor(4);
  fresy->Draw("same");
  cres->SaveAs("fres.png");

  int nbinsx = 15;
  int nbinsy = 10;
  int nbinsy2 = 8;
  double miny2 = 4.; //truncating the measured axis
  double maxy2 = 20.;
  TH2F *h2T = new TH2F("h2T","h2T",nbinsx,minx,maxx,nbinsy,miny,maxy); //truth
  TH2F *h2M = new TH2F("h2M","h2M",nbinsx,minx,maxx,nbinsy2,miny2,maxy2); //measured
  TH2F *h2P = new TH2F("h2P","h2P",nbinsx,minx,maxx,nbinsy,miny,maxy); //prior
  TH2F *h2TInRM = new TH2F("h2TInRM","h2TInRM",nbinsx,minx,maxx,nbinsy,miny,maxy); //truth

  RooUnfoldResponse *fResponse = new RooUnfoldResponse("resp","resp");
  fResponse->Setup(h2M,h2T);

  const int ntoys = 100000;
  TParameter<Double_t> *par = new TParameter<Double_t>("nEvt",(double)ntoys);

  //---------------------------------------------------------
  //  generate truth and measured distribution
  //---------------------------------------------------------
  // for(int i = 0; i<ntoys; ++i) {

  //   double x = f1->GetRandom();
  //   double y = f2->GetRandom();
  //   if(x>=minx && x<maxx && y>=miny && y<maxy)
  //     h2T->Fill(x,y);

  //   double x2 = x*fresx->GetRandom();
  //   double y2 = y*fresy->GetRandom();
  //   if(x2>=minx && x2<maxx && y2>=miny2 && y2<maxy2)
  //     h2M->Fill(x2,y2);
  // }

  //---------------------------------------------------------
  //  build response matrix
  //---------------------------------------------------------
  for(int i = 0; i<ntoys; ++i) {

    double x = f1->GetRandom();
    double y = f2->GetRandom();
    if(x>=minx && x<maxx && y>=miny && y<maxy){
      h2P->Fill(x,y);
      h2T->Fill(x,y);
    }

    double x2 = x*fresx->GetRandom();
    double y2 = y*fresy->GetRandom();

    //fResponse->Fill(x2,y2,x,y);

    if(x2>=minx && x2<maxx && y2>=miny2 && y2<maxy2) {
      fResponse->Fill(x2,y2,x,y);
      h2TInRM->Fill(x,y);
      h2M->Fill(x2, y2);
    } else {
      fResponse->Miss(x,y);
      // std::cout << "not in reponse. x2=" << x2 << "  y2=" << y2 << "  x=" << x << "  y=" << y << std::endl;
      // std::cout << "Not in response. ";
      // if (x2 < minx) { std::cout << "x2 < minx: " << x2 << " < " << minx << " "; }
      // if (x2 >= maxx) { std::cout << "x2 >= maxx: " << x2 << " >= " << maxx << " "; }
      // if (y2 < miny2) { std::cout << "y2 < miny2: " << y2 << " < " << miny2 << " "; }
      // if (y2 >= maxy2) { std::cout << "y2 >= maxy2: " << y2 << " >= " << maxy2 << " "; }
      // std::cout << endl;
    }

  }

  //---------------------------------------------------------
  //  Kinematic efficiency
  //---------------------------------------------------------
  TH2F *hKinEff = dynamic_cast<TH2F*>(h2TInRM->Clone("hKinEff"));
  hKinEff->SetTitle("hKinEff");
  hKinEff->Divide(h2P);

  //---------------------------------------------------------
  //  Fold prior with response
  //---------------------------------------------------------
  //apply kinematic efficciency first (hKinEff)
  TH2F *h2PKinEff = dynamic_cast<TH2F*>(h2P->Clone("h2PKinEff"));
  h2PKinEff->SetTitle("h2PKinEff");
  for(int i = 1; i<h2PKinEff->GetNbinsX()+1; ++i) {
    for(int j = 1; j<h2PKinEff->GetNbinsY()+1; ++j) {
      double val = h2P->GetBinContent(i,j);// * hKinEff->GetBinContent(i,j); //uncomment this if you don't use Miss AND do not include overflows in response
      h2PKinEff->SetBinContent(i,j,val);
    }
  }
  TH2F *hPriorFolded = (TH2F*)fResponse->ApplyToTruth(h2PKinEff,"hPriorFolded");
  //hPriorFolded->SetName("hPriorFolded");

  TCanvas *c3 = new TCanvas("c3","c3",400,400);
  gPad->SetLogz();
  h2T->Draw("colz");
  c3->SaveAs("h2T.png");

  TCanvas *c4 = new TCanvas("c4","c4",400,400);
  gPad->SetLogz();
  h2M->Draw("colz");
  c4->SaveAs("h2M.png");

  //---------------------------------------------------------
  //  unfold
  //---------------------------------------------------------
  int iterDef = 1;//3;
  int iterMin = 1;//3;
  int iterMax = 1;//3;
  int nIterTmp = iterMax - iterMin+1;
  const int nIter = nIterTmp;

  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;//RooUnfold::kCovToy;//
  RooUnfoldBayes unfold[nIter];
  TH2F *hReco[nIter];
  TH2F *hFolded[nIter];
  TMatrixD covmat[nIter];
  TH2D *hCovMat[nIter];
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {
    unfold[iter-iterMin] = RooUnfoldBayes(fResponse, h2M, iter);
    hReco[iter-iterMin] = (TH2F*)unfold[iter-iterMin].Hreco(errorTreatment);
    hReco[iter-iterMin]->SetName(Form("hReco_Iter%d",iter));

    hFolded[iter-iterMin] = (TH2F*)fResponse->ApplyToTruth(hReco[iter-iterMin],Form("hFolded_Iter%d",iter));
    //hFolded[iter-iterMin]->SetName(Form("hFolded_Iter%d",iter));

    //Get covariance matrix
    std::cout << "get covariance for iter " << iter << std::endl;
    TH2D htmp(unfold[iter-iterMin].Ereco(errorTreatment));
    hCovMat[iter-iterMin] = dynamic_cast<TH2D*>(htmp.Clone(Form("hCovMat_Iter%d",iter)));//TH2D::TH2D(covmat);
    hCovMat[iter-iterMin]->SetName(Form("hCovMat_Iter%d",iter));
  }

  //---------------------------------------------------------
  //  Get covariance matrix of default iter and calculate Pearson coefficients
  //---------------------------------------------------------
  std::cout << "get covariance for iterDef " << iterDef << std::endl;
  TMatrixD covmatDef = unfold[iterDef-iterMin].Ereco(errorTreatment);
  TH2D htmp(covmatDef);
  TH2D *hCovMatDef = dynamic_cast<TH2D*>(htmp.Clone(Form("hCovMatDef%d",iterDef)));
  TH2D *hPearsonDef =  dynamic_cast<TH2D*>(hCovMatDef->Clone("hPearsonDef"));
  hPearsonDef->Reset();

  //attempt to calculate Pearson coefficients
  for(int ibx = 0; ibx<nbinsy; ++ibx) {   //loop over blocks in cov matrix
    for(int iby = 0; iby<nbinsy; ++iby) {   //loop over blocks in cov matrix
      //std::cout << "ibx, iby " << ibx << "  " << iby << std::endl;
      for(int ix = ibx*nbinsx; ix<(ibx*nbinsx+nbinsx); ++ix) {
        for(int iy = iby*nbinsx; iy<(iby*nbinsx+nbinsx); ++iy) {
          //  std::cout << "ix, iy " << ix << "  " << iy << std::endl;
          double pearson = hCovMatDef->GetBinContent(ix+1,iy+1)/TMath::Sqrt(hCovMatDef->GetBinContent(ix+1,ix+1)*hCovMatDef->GetBinContent(iy+1,iy+1));
          hPearsonDef->SetBinContent(ix,iy,pearson);
        }
      }
    }
  }

  //Unfolding matrix. MV: don't know what this is
  TH2D htmp2(unfold[iterDef-iterMin].UnfoldingMatrix());;
  TH2D *hUnfMat = dynamic_cast<TH2D*>(htmp2.Clone(Form("hUnfMatDef%d",iterDef)));



  //---------------------------------------------------------
  //  Plotting of default iter
  //---------------------------------------------------------
  TCanvas *c5 = new TCanvas("c5","c5",800,750);
  c5->Divide(2,2);
  c5->cd(1);
  gPad->SetLogz();
  hReco[iterDef-iterMin]->Draw("colz");

  TH1D *hRecoP[2];
  hRecoP[0] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin]->ProjectionX("hRecoP_x"));
  hRecoP[1] = dynamic_cast<TH1D*>(hReco[iterDef-iterMin]->ProjectionY("hRecoP_y"));

  TH1D *hFoldedP[2];
  hFoldedP[0] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin]->ProjectionX("hFoldedP_x"));
  hFoldedP[1] = dynamic_cast<TH1D*>(hFolded[iterDef-iterMin]->ProjectionY("hFoldedP_y"));

  TH1D *hSP[2];
  hSP[0] = dynamic_cast<TH1D*>(h2M->ProjectionX("hSP_x"));
  hSP[1] = dynamic_cast<TH1D*>(h2M->ProjectionY("hSP_y"));

  TH1D *hTP[2];
  hTP[0] = dynamic_cast<TH1D*>(h2T->ProjectionX("hTP_x"));
  hTP[1] = dynamic_cast<TH1D*>(h2T->ProjectionY("hTP_y"));

  TH1D *hPriorFoldedP[2];
  hPriorFoldedP[0] = dynamic_cast<TH1D*>(hPriorFolded->ProjectionX("hPriorFoldedP_x"));
  hPriorFoldedP[1] = dynamic_cast<TH1D*>(hPriorFolded->ProjectionY("hPriorFoldedP_y"));


  c5->cd(3);
  gPad->SetLogy();
  hTP[0]->SetLineColor(1);
  hTP[0]->SetLineWidth(3);
  hSP[0]->SetLineColor(2);
  hSP[0]->SetLineWidth(3);
  hTP[0]->Draw();       //truth
  hSP[0]->Draw("same"); //measured
  hRecoP[0]->SetMarkerStyle(24);
  hRecoP[0]->SetMarkerColor(4);
  hRecoP[0]->SetLineColor(4);
  hRecoP[0]->Draw("same"); //unfolded

  hFoldedP[0]->SetMarkerStyle(24);
  hFoldedP[0]->SetMarkerColor(kMagenta+2);
  hFoldedP[0]->SetLineColor(kMagenta+2);
  hFoldedP[0]->SetLineWidth(3);
  hFoldedP[0]->SetLineStyle(2);
  hFoldedP[0]->Draw("same"); //refolded

  //hPriorFolded->SetMarkerStyle(24);
  hPriorFoldedP[0]->SetMarkerColor(kGreen+2);
  hPriorFoldedP[0]->SetLineColor(kGreen+2);
  hPriorFoldedP[0]->SetLineStyle(3);
  hPriorFoldedP[0]->SetLineWidth(3);
  // hPriorFoldedP[0]->Draw("same"); //folded prior

  TLegend *leg53 = new TLegend(0.15,0.6,0.5,0.9);
  leg53->SetFillColor(10);
  leg53->SetBorderSize(0);
  leg53->SetFillStyle(0);
  leg53->SetTextSize(0.05);
  leg53->SetTextFont(42);
  leg53->AddEntry(hTP[0],"truth","l");
  leg53->AddEntry(hSP[0],"measured","l");
  leg53->AddEntry(hRecoP[0],"unfolded","lp");
  leg53->AddEntry(hFoldedP[0],"unfolded #times RM","l");
  leg53->AddEntry(hPriorFoldedP[0],"prior #times RM","l");
  leg53->Draw();

  c5->cd(2);
  gPad->SetLogy();
  hSP[0]->Draw();       //measured
  hPriorFoldedP[0]->Draw("same");
  hFoldedP[0]->Draw("same");

  c5->cd(4);
  gPad->SetLogy();
  hTP[1]->SetLineColor(1);
  hTP[1]->SetLineWidth(3);
  hSP[1]->SetLineColor(2);
  hSP[1]->SetLineWidth(3);
  hTP[1]->Draw();
  hSP[1]->Draw("same");
  hRecoP[1]->SetMarkerStyle(24);
  hRecoP[1]->SetMarkerColor(4);
  hRecoP[1]->SetLineColor(4);
  hRecoP[1]->Draw("same");

  hFoldedP[1]->SetMarkerStyle(24);
  hFoldedP[1]->SetMarkerColor(kMagenta+2);
  hFoldedP[1]->SetLineColor(kMagenta+2);
  hFoldedP[1]->SetLineWidth(3);
  hFoldedP[1]->SetLineStyle(2);
  hFoldedP[1]->Draw("same");

  //hPriorFolded->SetMarkerStyle(24);
  hPriorFoldedP[1]->SetMarkerColor(kGreen+2);
  hPriorFoldedP[1]->SetLineColor(kGreen+2);
  hPriorFoldedP[1]->SetLineWidth(3);
  hPriorFoldedP[1]->SetLineStyle(3);
  hPriorFoldedP[1]->Draw("same"); //folded prior

  c5->SaveAs("control.png");


  //divisions
  TH2F *hRatFolded[nIter];
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {
    hRatFolded[iter-iterMin] = dynamic_cast<TH2F*>(hFolded[iter-iterMin]->Clone(Form("hRatFolded_Iter%d",iter)));
    hRatFolded[iter-iterMin]->Divide(h2M);
  }

  TCanvas *c6 = new TCanvas("c6","c6",800,750);
  c6->Divide(2,2);
  c6->cd(1);
  gPad->SetLogz();
  h2M->Draw("colz");
  c6->cd(2);
  gPad->SetLogz();
  hFolded[iterDef-iterMin]->Draw("colz");
  c6->cd(3);
  hRatFolded[iterDef-iterMin]->Draw("colz");

  c6->SaveAs("refolded.png");

  TCanvas *c7 = new TCanvas("c7","c7",800,750);
  gPad->SetLogz();
  hCovMat[iterDef-iterMin]->Draw("colz");
  c7->SaveAs("covmat.png");

  TCanvas *c8 = new TCanvas("c8","c8",800,750);
  gPad->SetLogz();
  fResponse->Hresponse()->Draw("colz");
  c8->SaveAs("response.png");

  TCanvas *c9 = new TCanvas("c9","c9",800,750);
  //gPad->SetLogz();
  hPearsonDef->Draw("colz");
  c9->SaveAs("pearson.png");

  TCanvas *c10 = new TCanvas("c10","c10",800,750);
  c10->Divide(2,2);
  c10->cd(1);
  gPad->SetLogz();
  h2P->Draw("colz");
  c10->cd(2);
  gPad->SetLogz();
  h2TInRM->Draw("colz");
  c10->cd(3);
  hKinEff->Draw("colz");
  c10->cd(4);
  gPad->SetLogz();
  h2PKinEff->Draw("colz");

  c10->SaveAs("kineff.png");

  TFile *fout = new TFile(Form("UnfoldedDistributions.root"),"RECREATE");
  h2P->Write("fh2Prior");
  h2T->Write("fh2True");
  h2TInRM->Write();
  hKinEff->Write();
  h2M->Write("fh2Smear");
  par->Write();
  //fResponse->Write();
  fResponse->Hresponse()->Write("hResponse");
  hPriorFolded->Write();
  for(Int_t iter = iterMin; iter<=iterMax; iter++) {
    hReco[iter-iterMin]->Write();
    hFolded[iter-iterMin]->Write();
    hCovMat[iter-iterMin]->Write();
    //covmat[iter-iterMin].Write(Form("covmat%d",iter));
    //unfold[iter-iterMin].Write(Form("unfold%d",iter));
  }
  hCovMatDef->Write("hCovMatDef");
  hPearsonDef->Write();
  hUnfMat->Write("hUnfMat");

  fout->Write();
  fout->Close();

}

int main() {
  unfoldingTest();
  return 0;
}
