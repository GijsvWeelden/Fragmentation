
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

//-------------------------------------------------------------
//
// plot_jetprops.C normalises and plots histograms against
// each other and compares them by taking their ratio
// it can do this for either 2 or 3 histograms and takes the
// pp as a reference
//
//-------------------------------------------------------------

void charged_VS_full(TH2F *hF, TH2F *hC, string obs, double min_pt, double max_pt, string setting, string type, string hadron);

void plot_hadron_frag(void){
  double time = clock();
  gROOT->SetBatch();
  string inName = "new_2dhists_frag";
  string sNN = "5tev02";
  std::vector<double> ptBins = {40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> obs = {"frag"};//,"orth"};
  std::vector<string> hadrons = {"pi", "K", "p"};
  //std::vector<string> chHadrons = {"pi", "K", "p"};
  //std::vector<string> nHadrons = {"pi0", "K0L", "K0S", "K0", "Lambda0"};


  //new_2dhists_frag_5tev02_pp_hadron_charged.root
  TFile *inFile = TFile::Open(TString::Format("./%s_%s_pp_hadron.root", inName.c_str(), sNN.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  TList *chList = (TList*)inFile->Get("pp_charged");
  if (!chList){
    std::cout << "Error: charged not found! Aborting program." << std::endl;
    return;
  }
  TList *fList = (TList*)inFile->Get("pp_full");
  if (!fList){
    std::cout << "Error: full not found! Aborting program." << std::endl;
    return;
  }

  TH2F *hC, *hF;
  //TH1F *hCharged, *hFull;

  //TCanvas *c1 = new TCanvas("name", "title", 900, 600);
  //c1->cd();

  //hC->GetYaxis()->SetRangeUser(ptBins[0], ptBins.back());
  /*
  hCharged = (TH1F*)hC->ProjectionX();
  hCharged->Scale(1./hCharged->Integral()); // Use option "width" to divide by bin width for jT
  hCharged->SetStats(0);
  //hCharged->SetMarkerStyle(kCircle);
  //hCharged->SetLineColor(kOrange+7);
  //hCharged->SetMarkerColor(kOrange+7);

  hF->GetYaxis()->SetRangeUser(ptBins[0], ptBins.back());
  hFull = (TH1F*)hF->ProjectionX();
  hFull->Scale(1./hFull->Integral());
  hFull->SetStats(0);
  hFull->SetMarkerStyle(kCircle);
  hFull->SetLineColor(kOrange+7);
  hFull->SetMarkerColor(kOrange+7);

  hCharged->Draw();
  hFull->Draw("same");
*/

  // /*
  for (int iobs=0; iobs<obs.size(); iobs++){
    std::cout << "Plotting " << obs[iobs] << std::endl;
    hC = (TH2F*)chList->FindObject(TString::Format("pp_%s", obs[iobs].c_str()).Data());
    hF = (TH2F*)fList->FindObject(TString::Format("pp_%s", obs[iobs].c_str()).Data());
    charged_VS_full(hF, hC, obs[iobs], ptBins[0], ptBins.back(), "pp", "", "");
    for (int ipt=0; ipt<ptBins.size()-1; ++ipt){
      //charged_VS_full(hF, hC, obs[iobs], ptBins[ipt], ptBins.[ipt+1], "pp", "", "");
    }
    for (int ihad = 0; ihad < hadrons.size(); ihad++){
      hC = (TH2F*)chList->FindObject(TString::Format("pp_charged_%s_%s", hadrons[ihad].c_str(), obs[iobs].c_str()).Data());
      hF = (TH2F*)fList->FindObject(TString::Format("pp_full_%s_%s", hadrons[ihad].c_str(), obs[iobs].c_str()).Data());
      charged_VS_full(hF, hC, "frag", ptBins[0], ptBins.back(), "pp", "", hadrons[ihad]);
      for (int ipt = 0; ipt < ptBins.size()-1; ++ipt){
        //charged_VS_full(hF, hC, obs[iobs], ptBins[ipt], ptBins[ipt+1], "pp", "", "");
      }
    }
  }
  // */
  //hC = (TH2F*)chList->FindObject("pp_charged_pi_frag");
  //hF = (TH2F*)fList->FindObject("pp_full_pi_frag");

  /*
  std::cout << "Charged: " << hC->Integral() << std::endl
    << "Full: " << hF->Integral() << std::endl;

  double min_pt = ptBins[0];
  double max_pt = ptBins.back();
  hF->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hFull = (TH1F*)hF->ProjectionX();
  //hFull->Scale(1./hFull->Integral());
  hC->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hCharged = (TH1F*)hC->ProjectionX();
  //hCharged->Scale(2./hCharged->Integral());
  std::cout << "Charged: " << hC->Integral() << std::endl
    << "Full: " << hF->Integral() << std::endl;

  hCharged->SetStats(0);
  hCharged->SetMarkerStyle(kMultiply);
  hCharged->SetLineColor(kRed);
  hCharged->SetMarkerColor(kRed);
  */

  //TCanvas *cC = new TCanvas("charged", "charged", 400, 400);
  //cC->cd();
  //hFull->Draw();
  //TCanvas *cF = new TCanvas("full", "full", 400, 400);
  //cF->cd();
  //hCharged->Draw("same");
  //charged_VS_full(hF, hC, "frag", ptBins[0], ptBins.back(), "pp", "", "pi");


  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  //chFile->Close();
  //nFile->Close();
}

void charged_VS_full(TH2F *hF, TH2F *hC, string obs, double min_pt, double max_pt, string setting, string type, string hadron){
  double jetR = 0.4;
  // Project with pt cut
  hF->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hFull = (TH1F*)hF->ProjectionX();
  //hFull->Scale(1./hFull->Integral());
  hC->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hCharged = (TH1F*)hC->ProjectionX();
  //hCharged->Scale(1./hCharged->Integral());
  if (abs(hFull->Integral() -  1.0) > 0.01){
    cout << "WARNING: NORMALISATION PROBLEM!" << endl
      << TString::Format("%s for pt %.0f-%.0f GeV/c", obs.c_str(), min_pt, max_pt).Data() << endl
      << "Full = " << hFull->Integral() << " (" << hFull->Integral()-1.0  << ")" << endl;
  }
  if (abs(hCharged->Integral() - 1.0) > 0.01){
    cout << "WARNING: NORMALISATION PROBLEM!" << endl
      << TString::Format("%s for pt %.0f-%.0f GeV/c", obs.c_str(), min_pt, max_pt).Data() << endl
      << "Charged = " << hCharged->Integral() << " (" << hCharged->Integral()-1.0 << ")" << endl;
  }

  // Top plot settings
  hFull->SetStats(0);
  //hFull->SetMarkerStyle(kPlus);
  //hFull->SetLineColor(kBlue);
  //hFull->SetMarkerColor(kBlue);
  if (type == "all"){
    hFull->SetTitle(TString::Format("Jet %s, pt #in [%.0f,%.0f) GeV/c, R = %.1f %s",
                                    obs.c_str(), min_pt, max_pt, jetR, setting.c_str()).Data());
  }
  else{
    hFull->SetTitle(TString::Format("%s %s, pt #in [%.0f,%.0f) GeV/c, R = %.1f %s",
                                    hadron.c_str(), obs.c_str(), min_pt, max_pt, jetR, type.c_str()).Data());
  }
  hFull->GetYaxis()->SetTitle(TString::Format("#frac{dN}{d %s}",
                                              obs.c_str()).Data());
  hFull->GetYaxis()->SetTitleFont(43);
  hFull->GetYaxis()->SetTitleSize(20);
  hFull->GetYaxis()->SetTitleOffset(1.5);//2.5);
  hFull->GetYaxis()->SetLabelFont(43);
  hFull->GetYaxis()->SetLabelSize(14);

  hCharged->SetStats(0);
  hCharged->SetMarkerStyle(kMultiply);
  hCharged->SetLineColor(kRed);
  hCharged->SetMarkerColor(kRed);

  // Ratio plot settings
  TH1F *ratioFull = (TH1F*)hFull->Clone("Ratio Base/Base");
  ratioFull->Divide(hFull);
  ratioFull->SetTitle("");
  ratioFull->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());
  ratioFull->GetXaxis()->SetTitleFont(43);
  ratioFull->GetXaxis()->SetTitleSize(20);
  ratioFull->GetXaxis()->SetTitleOffset(2.5);
  ratioFull->GetXaxis()->SetLabelFont(43);
  ratioFull->GetXaxis()->SetLabelSize(14);
  ratioFull->GetYaxis()->SetTitle("Ratio");
  ratioFull->GetYaxis()->SetTitleFont(43);
  ratioFull->GetYaxis()->SetTitleSize(20);
  ratioFull->GetYaxis()->SetTitleOffset(1.5);
  ratioFull->GetYaxis()->SetLabelFont(43);
  ratioFull->GetYaxis()->SetLabelSize(14);

  TH1F *ratioCharged = (TH1F*)hCharged->Clone("Ratio X/Base");
  ratioCharged->Divide(hFull);

  int XmaxBinR = ratioFull->FindLastBinAbove(0,1);
  int XminBinR = ratioFull->FindFirstBinAbove(0,1);
  double Ymin = 0;
  double Ymax = max({hFull->GetMaximum(),
      hCharged->GetMaximum()
      });
  double YminR = min({ratioFull->GetMinimum(),
      ratioCharged->GetMinimum()
      });
  double YmaxR = max({ratioFull->GetMaximum(),
      ratioCharged->GetMaximum()
      });

  hFull->GetXaxis()->SetRange(0,XmaxBinR);
  ratioFull->GetXaxis()->SetRange(0,XmaxBinR);
  if (YmaxR > 3) YmaxR = 3;
  else YmaxR *= 1.1;
  if (obs == "frag"){
    if (hadron == "p") Ymin = 1e-9; //min({hFull->GetMinimum(), hCharged->GetMinimum()}); //1e-6;
    else Ymin = 1e-8;
  }
  else if (obs == "orth") Ymin = 1e-11;

  hFull->GetYaxis()->SetRangeUser(Ymin, 2*Ymax);
  ratioFull->GetYaxis()->SetRangeUser(0.9*YminR,YmaxR);

  TCanvas *c_obs = new TCanvas(TString::Format("c_%s_pt_%.0f_%.0f",
                                               obs.c_str(), min_pt, max_pt).Data(),
                               TString::Format("%s_pt_%.0f_%.0f",
                                               obs.c_str(), min_pt, max_pt).Data(),
                               900, 600); //600, 900);
  TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.15);
  pad1->Draw();
  pad1->SetLogy();
  if (obs == "orth") pad1->SetLogx();
  pad1->cd();

  hFull->Draw();
  if (hCharged->GetEntries() != 0) hCharged->Draw("same");

  auto legend = new TLegend(0.6,0.7,0.9,0.9); // Top right corner

  if (type == "all"){
    legend->AddEntry(hFull,"Full set");
    if (hCharged->GetEntries() != 0) legend->AddEntry(hCharged,"Leading");
  }
  else{
    legend->AddEntry(hFull,"Full jets");
    if (hCharged->GetEntries() != 0) legend->AddEntry(hCharged, "Charged jets");
  }
  legend->Draw();

  c_obs->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.15);
  pad2->Draw();
  if (obs == "orth") pad2->SetLogx();
  pad2->cd();

  ratioFull->Draw();
  if (ratioCharged->GetEntries() != 0) ratioCharged->Draw("same");
  c_obs->cd();

  c_obs->SaveAs(TString::Format("../plots/fragmentation/ChargedVsFull/%s_%s_%s_pt_%.0f_%.0f.pdf",
                                hadron.c_str(), setting.c_str(), obs.c_str(), min_pt, max_pt).Data());
}

