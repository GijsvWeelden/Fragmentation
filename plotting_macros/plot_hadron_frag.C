
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

using std::cout;
using std::endl;

//-------------------------------------------------------------
//
// Compares fragmentation spectra in charged and full jets
// both for inclusive sample as for specific hadron species
//
//-------------------------------------------------------------

void charged_VS_full(TH2F *hF, TH2F *hC, string obs, double min_pt, double max_pt, string setting, string hadron);
void baryon_meson(TH2F *hB, string baryon, TH2F* hM, string meson, string obs, double min_pt, double max_pt, string setting);
void fit_spectra(TH1F* hFull, TH1F* hCharged, string obs, string hadron);

void plot_hadron_frag(void){
  double time = clock();
  // gROOT->SetBatch();
  string inName = "partial_2dhists_frag";
  string sNN = "5tev02";
  std::vector<double> ptBins = {40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> observables = {"frag"};//, "orth", "track_pt"};

  std::vector<string> chMesons = {"pi", "K"};
  std::vector<string> chBaryons = {"p"};
  std::vector<string> chHadrons(chMesons);
  chHadrons.insert(chHadrons.end(), chBaryons.begin(), chBaryons.end());

  std::vector<string> nMesons = {"pi0", "K0L", "K0S", "K0"};
  std::vector<string> nBaryons = {"Lambda0"};
  std::vector<string> nHadrons(nMesons);
  nHadrons.insert(nHadrons.end(), nBaryons.begin(), nBaryons.end());

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

  /*
  TH2F *hC, *hF;
  for (auto obs : observables){
    std::cout << "Plotting " << obs << std::endl;
    hC = (TH2F*)chList->FindObject(TString::Format("pp_charged_%s", obs.c_str()).Data());
    hF = (TH2F*)fList->FindObject(TString::Format("pp_full_%s", obs.c_str()).Data());
    charged_VS_full(hF, hC, obs, ptBins[0], ptBins.back(), "pp", "");
    for (int ipt=0; ipt<ptBins.size()-1; ++ipt){
      charged_VS_full(hF, hC, obs, ptBins[ipt], ptBins[ipt+1], "pp", "");
    }
    for (auto hadron : chHadrons){
      hC = (TH2F*)chList->FindObject(TString::Format("pp_charged_%s_%s", hadron.c_str(), obs.c_str()).Data());
      hF = (TH2F*)fList->FindObject(TString::Format("pp_full_%s_%s", hadron.c_str(), obs.c_str()).Data());
      charged_VS_full(hF, hC, obs, ptBins[0], ptBins.back(), "pp", hadron);
      for (int ipt = 0; ipt < ptBins.size()-1; ++ipt){
        charged_VS_full(hF, hC, obs, ptBins[ipt], ptBins[ipt+1], "pp", hadron);
      }
    }
  }
  */

  TH2F *hM, *hB;
  string obs = "track_pt";
  // for (auto meson : chMesons){
  // for (auto baryon : chBaryons){
  hM = (TH2F*)chList->FindObject(TString::Format("pp_charged_pi_%s", obs.c_str()).Data());
  hB = (TH2F*)chList->FindObject(TString::Format("pp_charged_p_%s", obs.c_str()).Data());
  baryon_meson(hB, "p", hM, "pi", obs, ptBins[0], ptBins.back(), "pp");
  for (int i = 0; i < ptBins.size()-1; i++){
    baryon_meson(hB, "p", hM, "pi", obs, ptBins[i], ptBins[i+1], "pp");
  }
  // }
  // }
  // Do same for neutral

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
}

void charged_VS_full(TH2F *hF, TH2F *hC, string obs, double min_pt, double max_pt, string setting, string hadron){
  double jetR = 0.4;
  if (hadron == "") hadron = "jet";
  hF->GetYaxis()->SetRangeUser((3/2)*min_pt, (3/2)*max_pt);
  TH1F *hFull = (TH1F*)hF->ProjectionX();
  hC->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hCharged = (TH1F*)hC->ProjectionX();
  if (obs == "orth"){
    hFull->Scale(1./hFull->Integral(),"width");
    hCharged->Scale(1./hCharged->Integral(),"width");
  }

  // fit_spectra(hFull, hCharged, "frag", "pi");

  // Top plot settings
  hFull->SetStats(0);
  hFull->SetTitle(TString::Format("%s %s, pt_{ch} #in [%.0f,%.0f) GeV/c, R = %.1f %s",
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
  double Ymin = min({hFull->GetMinimum(), hCharged->GetMinimum()});
  Ymin *= 0.5;
  double Ymax = max({hFull->GetMaximum(), hCharged->GetMaximum()});
  Ymax *= 2;
  if (Ymin == 0){
    double tmpFull = hFull->FindLastBinAbove(0,1);
    double tmpCharged = hCharged->FindLastBinAbove(0,1);
    Ymin = min({hFull->GetBinContent(tmpFull), hCharged->GetBinContent(tmpCharged)});
  }
  double YminR = min({ratioFull->GetMinimum(), ratioCharged->GetMinimum()});
  YminR *= 0.9;
  double YmaxR = max({ratioFull->GetMaximum(), ratioCharged->GetMaximum()});
  if (YmaxR > 3) YmaxR = 3;
  else YmaxR *= 1.1;

  hFull->GetXaxis()->SetRange(0,XmaxBinR);
  ratioFull->GetXaxis()->SetRange(0,XmaxBinR);

  hFull->GetYaxis()->SetRangeUser(Ymin, Ymax);
  ratioFull->GetYaxis()->SetRangeUser(YminR,YmaxR);

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
  if (obs == "orth"){
    delete legend;
    legend = new TLegend(0.4, 0, 0.7, 0.2); // Bottom middle
  }
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

  // c_obs->SaveAs(TString::Format("../plots/fragmentation/ChargedVsFull/%s_%s_%s_ptch_%.0f_%.0f.pdf",
                                // setting.c_str(), hadron.c_str(), obs.c_str(), min_pt, max_pt).Data());
  // delete c_obs;
}

void baryon_meson(TH2F *hB, string baryon, TH2F* hM, string meson, string obs, double min_pt, double max_pt, string setting){
  double jetR = 0.4;
  hB->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hBaryon = (TH1F*)hB->ProjectionX();
  hM->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hMeson = (TH1F*)hM->ProjectionX();
  // hBaryon->Scale(1./hBaryon->Integral());
  // hMeson->Scale(1./hMeson->Integral());

  // Top plot settings
  hBaryon->SetStats(0);
  hBaryon->SetTitle(TString::Format("%s/%s, pt_{ch} #in [%.0f,%.0f) GeV/c, R = %.1f %s",
                                    baryon.c_str(), meson.c_str(), min_pt, max_pt, jetR, setting.c_str()).Data());
  hBaryon->GetYaxis()->SetTitle(TString::Format("%s/%s", baryon.c_str(), meson.c_str()).Data());
  hBaryon->GetYaxis()->SetTitleFont(43);
  hBaryon->GetYaxis()->SetTitleSize(20);
  hBaryon->GetYaxis()->SetTitleOffset(1.5);
  hBaryon->GetYaxis()->SetLabelFont(43);
  hBaryon->GetYaxis()->SetLabelSize(14);

  hMeson->SetStats(0);
  hMeson->SetMarkerStyle(kMultiply);
  hMeson->SetLineColor(kRed);
  hMeson->SetMarkerColor(kRed);

  // Ratio plot settings
  TH1F *ratio = (TH1F*)hBaryon->Clone("Ratio Base/Base");
  ratio->Divide(hMeson);
  // ratio->SetTitle("");
  ratio->GetXaxis()->SetTitle(TString::Format("#it{p}_T").Data());
  ratio->GetXaxis()->SetTitleFont(43);
  ratio->GetXaxis()->SetTitleSize(20);
  // ratio->GetXaxis()->SetTitleOffset(2.5);
  ratio->GetXaxis()->SetLabelFont(43);
  ratio->GetXaxis()->SetLabelSize(14);
  ratio->GetYaxis()->SetTitle("Ratio");
  ratio->GetYaxis()->SetTitleFont(43);
  ratio->GetYaxis()->SetTitleSize(20);
  // ratio->GetYaxis()->SetTitleOffset(1.5);
  ratio->GetYaxis()->SetLabelFont(43);
  ratio->GetYaxis()->SetLabelSize(14);

  TH1F *ratioBaryon = (TH1F*)hBaryon->Clone("Ratio X/Base");
  ratioBaryon->Divide(hMeson);

  int XmaxBinR = ratio->FindLastBinAbove(0,1);
  // int XminBinR = ratio->FindFirstBinAbove(0,1);
  // double Ymin = min({hBaryon->GetMinimum(), hMeson->GetMinimum()});
  // Ymin *= 0.5;
  // double Ymax = max({hBaryon->GetMaximum(), hMeson->GetMaximum()});
  // Ymax *= 2;
  // if (Ymin == 0){
  //   double tmpBaryon = hBaryon->FindLastBinAbove(0,1);
  //   double tmpMeson = hMeson->FindLastBinAbove(0,1);
  //   Ymin = min({hBaryon->GetBinContent(tmpBaryon), hMeson->GetBinContent(tmpMeson)});
  // }
  double YminR = ratio->GetMinimum();
  // YminR *= 0.9;
  double YmaxR = ratio->GetMaximum();
  if (YmaxR > 0.8) YmaxR = 0.8;
  // else YmaxR *= 1.1;
  // */

  // hBaryon->GetXaxis()->SetRange(0,XmaxBinR);
  ratio->GetXaxis()->SetRange(0,XmaxBinR);

  // hBaryon->GetYaxis()->SetRangeUser(Ymin, Ymax);
  ratio->GetYaxis()->SetRangeUser(YminR,YmaxR);

  TCanvas *c_obs = new TCanvas(TString::Format("c_%s_pt_%.0f_%.0f", obs.c_str(), min_pt, max_pt).Data(),
                               TString::Format("%s_pt_%.0f_%.0f", obs.c_str(), min_pt, max_pt).Data(),
                               900, 600);
  // TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
  // pad1->SetBottomMargin(0);
  // pad1->SetLeftMargin(0.15);
  // pad1->Draw();
  // // pad1->SetLogy();
  // if (obs == "orth") pad1->SetLogx();
  // pad1->cd();

  // hBaryon->Draw();
  // if (hMeson->GetEntries() != 0) hMeson->Draw("same");

  auto legend = new TLegend(0.1,0.7,0.3,0.9); // Top right corner
  if (legend) std::cout << "legend" << std::endl;
  // if (obs == "orth"){
  //   delete legend;
  //   legend = new TLegend(0.4, 0, 0.7, 0.2); // Bottom middle
  // }
  // legend->AddEntry(hBaryon, baryon.c_str());
  // legend->AddEntry(hMeson, meson.c_str());
  legend->AddEntry(ratioBaryon, TString::Format("#it{p}_{T}^{ch jet} %.0f - %.0f", min_pt, max_pt).Data());


  // c_obs->cd();
  // TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
  // pad2->SetTopMargin(0);
  // pad2->SetBottomMargin(0.25);
  // pad2->SetLeftMargin(0.15);
  // pad2->Draw();
  // if (obs == "orth") pad2->SetLogx();
  // pad2->cd();

  ratio->Draw();
  legend->Draw("same");
  // if (ratioBaryon->GetEntries() != 0) ratioBaryon->Draw();
  // c_obs->cd();

  // c_obs->SaveAs(TString::Format("../plots/fragmentation/ChargedVsFull/%s_%s_%s_ptch_%.0f_%.0f.pdf",
                                // setting.c_str(), hadron.c_str(), obs.c_str(), min_pt, max_pt).Data());
  // delete c_obs;
}

void fit_spectra(TH1F* hFull, TH1F* hCharged, string obs, string hadron){
  TFitResultPtr fit1 =  hFull->Fit("expo", "S", "", 0.4, 0.8);
  auto params1 = fit1->Parameters();
  cout << "Params " << hadron << " " << obs << " (full): ";
  for (auto element : params1){
    cout << element << ", ";
  }
  cout << endl;

  TFitResultPtr fit2 =  hCharged->Fit("expo", "S", "", 0.4, 0.8);
  auto params2 = fit2->Parameters();
  cout << "Params " << hadron << " " << obs << " (charged): ";
  for (auto element : params2){
    cout << element << ", ";
  }
  cout << endl;
  return;
}
