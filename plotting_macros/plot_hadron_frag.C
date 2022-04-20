
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
// Compares fragmentation spectra in charged and full jets
//
//-------------------------------------------------------------

void charged_VS_full(TH2F *hF, TH2F *hC, string obs, double min_pt, double max_pt, string setting, string type, string hadron);

void plot_hadron_frag(void){
  double time = clock();
  gROOT->SetBatch();
  string inName = "partial_2dhists_frag";
  string sNN = "5tev02";
  std::vector<double> ptBins = {40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> observables = {"frag", "orth"};
  std::vector<string> chHadrons = {"pi", "K", "p"};
  std::vector<string> nHadrons = {"pi0", "K0L", "K0S", "K0", "Lambda0"};

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
  for (auto obs : observables){
    std::cout << "Plotting " << obs << std::endl;
    hC = (TH2F*)chList->FindObject(TString::Format("pp_charged_%s", obs.c_str()).Data());
    hF = (TH2F*)fList->FindObject(TString::Format("pp_full_%s", obs.c_str()).Data());
    charged_VS_full(hF, hC, obs, ptBins[0], ptBins.back(), "pp", "", "");
    for (int ipt=0; ipt<ptBins.size()-1; ++ipt){
      charged_VS_full(hF, hC, obs, ptBins[ipt], ptBins[ipt+1], "pp", "", "");
    }
    for (auto hadron : chHadrons){
      hC = (TH2F*)chList->FindObject(TString::Format("pp_charged_%s_%s", hadron.c_str(), obs.c_str()).Data());
      hF = (TH2F*)fList->FindObject(TString::Format("pp_full_%s_%s", hadron.c_str(), obs.c_str()).Data());
      charged_VS_full(hF, hC, "frag", ptBins[0], ptBins.back(), "pp", "", hadron);
      for (int ipt = 0; ipt < ptBins.size()-1; ++ipt){
        charged_VS_full(hF, hC, obs, ptBins[ipt], ptBins[ipt+1], "pp", "", hadron);
      }
    }
  }

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void charged_VS_full(TH2F *hF, TH2F *hC, string obs, double min_pt, double max_pt, string setting, string type, string hadron){
  double jetR = 0.4;
  if (hadron == "") hadron = "jet";
  hF->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hFull = (TH1F*)hF->ProjectionX();
  hC->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hCharged = (TH1F*)hC->ProjectionX();
  if (obs == "orth"){
    hFull->Scale(1./hFull->Integral(),"width");
    hCharged->Scale(1./hCharged->Integral(),"width");
  }

  // Top plot settings
  hFull->SetStats(0);
  hFull->SetTitle(TString::Format("%s %s, pt #in [%.0f,%.0f) GeV/c, R = %.1f %s",
                                  hadron.c_str(), obs.c_str(), min_pt, max_pt, jetR, setting.c_str()).Data());
  hFull->GetYaxis()->SetTitle(TString::Format("#frac{dN}{d %s}", obs.c_str()).Data());
  hFull->GetYaxis()->SetTitleFont(43);
  hFull->GetYaxis()->SetTitleSize(20);
  hFull->GetYaxis()->SetTitleOffset(1.5);
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
    if (hadron == "p") Ymin = 1e-9; //min({hFull->GetMinimum(), hCharged->GetMinimum()});
    else Ymin = 1e-8;
  }
  else if (obs == "orth") Ymin = 1e-11;

  hFull->GetYaxis()->SetRangeUser(Ymin, 2*Ymax);
  ratioFull->GetYaxis()->SetRangeUser(0.9*YminR,YmaxR);

  TCanvas *c_obs = new TCanvas(TString::Format("c_%s_pt_%.0f_%.0f", obs.c_str(), min_pt, max_pt).Data(),
                               TString::Format("%s_pt_%.0f_%.0f", obs.c_str(), min_pt, max_pt).Data(),
                               900, 600);
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

  legend->AddEntry(hFull,"Full jets");
  if (hCharged->GetEntries() != 0) legend->AddEntry(hCharged, "Charged jets");
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
                                setting.c_str(), hadron.c_str(), obs.c_str(), min_pt, max_pt).Data());
  delete c_obs;
}

