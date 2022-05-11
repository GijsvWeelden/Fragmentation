
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

#include "./plotUtils.C"

using std::cout;
using std::endl;

//-------------------------------------------------------------
//
// * Compares fragmentation spectra in charged and full jets
// both for inclusive sample as for specific hadron species,
// and fits those spectra
// * Compares hadron spectra and shows their ratios for various
// jet pt bins
//
//-------------------------------------------------------------

void charged_VS_full(TH2F *hF, TH2F *hC, string obs, double min_pt, double max_pt, string setting, string hadron);
void spectra_ratio(TH2F *hN, string sNum, TH2F* hD, string sDen, std::vector<double> ptBins, string setting);
TH1F* make_ratio(TH2F* hNum, string sNum, TH2F* hDen, string sDen, double min_pt, double max_pt);
void fit_spectra(TH1F* hFull, TH1F* hCharged, string obs, string hadron);

void plot_hadron_frag(void){
  double time = clock();
  bool doChargedVSFull = true, doRatios = false;
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

  if (doChargedVSFull){
    TH2F *hC, *hF;
    for (auto obs : observables){
      std::cout << "Plotting " << obs << std::endl;
      hC = (TH2F*)chList->FindObject(TString::Format("pp_charged_%s", obs.c_str()).Data());
      hF = (TH2F*)fList->FindObject(TString::Format("pp_full_%s", obs.c_str()).Data());
      charged_VS_full(hF, hC, obs, ptBins[0], ptBins.back(), "pp", "");
      for (int ipt=0; ipt < ptBins.size()-1; ++ipt){
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
  }

  if (doRatios){
    TH2F *hPi, *hK, *hP, *hPi0, *hK0L, *hK0S, *hK0, *hLambda0;
    string obs = "track_pt";
    hPi = (TH2F*)chList->FindObject(TString::Format("pp_charged_pi_%s", obs.c_str()).Data());
    hK = (TH2F*)chList->FindObject(TString::Format("pp_charged_K_%s", obs.c_str()).Data());
    hP = (TH2F*)chList->FindObject(TString::Format("pp_charged_p_%s", obs.c_str()).Data());
    hPi0 = (TH2F*)fList->FindObject(TString::Format("pp_full_pi0_%s", obs.c_str()).Data());
    hK0L = (TH2F*)fList->FindObject(TString::Format("pp_full_K0L_%s", obs.c_str()).Data());
    hK0S = (TH2F*)fList->FindObject(TString::Format("pp_full_K0S_%s", obs.c_str()).Data());
    hK0 = (TH2F*)fList->FindObject(TString::Format("pp_full_K0_%s", obs.c_str()).Data());
    hLambda0 = (TH2F*)fList->FindObject(TString::Format("pp_full_Lambda0_%s", obs.c_str()).Data());

    spectra_ratio(hP, "p", hPi, "pi", ptBins, "pp");
    spectra_ratio(hK, "K", hPi, "pi", ptBins, "pp");
    spectra_ratio(hK0L, "K0L", hPi0, "pi0", ptBins, "pp");
    // TODO: Add others as needed
  }

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
    hFull->Scale(1./hFull->Integral(), "width");
    hCharged->Scale(1./hCharged->Integral(), "width");
  }

  // fit_spectra(hFull, hCharged, "frag", "pi");
  double fitmin = 0.4;
  double fitmax = 0.8;
  TFitResultPtr fitFull =  hFull->Fit("expo", "S", "", fitmin, fitmax);
  TFitResultPtr fitCharged = hCharged->Fit("expo", "S", "", fitmin, fitmax);
  // return;

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
                               1600, 1000);
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
  auto paramsFull = fitFull->Parameters();
  auto paramsCharged = fitCharged->Parameters();
  legend->AddEntry(hFull, TString::Format("Full: exp(%.1f+%.1f x)", paramsCharged[0], paramsCharged[1]).Data());
  if (hCharged->GetEntries() != 0) legend->AddEntry(hCharged, TString::Format("Charged: exp(%.1f+%.1f x)", paramsFull[0], paramsFull[1]).Data());
  // legend->AddEntry(fitFull, TString::Format("exp(%.1f+%.1fx)", paramsFull[0], paramsFull[1]).Data())
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

  c_obs->SaveAs(TString::Format("../plots/fragmentation/ChargedVsFull/%s_%s_%s_ptch_%.0f_%.0f.pdf",
                                setting.c_str(), hadron.c_str(), obs.c_str(), min_pt, max_pt).Data());
  delete c_obs;
}

void spectra_ratio(TH2F* hNum, string sNum, TH2F* hDen, string sDen, std::vector<double> ptBins, string setting){
  double jetR = 0.4;
  std::vector<TH1F*> hists;
  std::vector<TH1F*> ratios;
  TH1F *hRatio, *ptRatio;

  TCanvas *c_obs = new TCanvas(TString::Format("c_%s_%s_%s", sNum.c_str(), sDen.c_str(), setting.c_str()).Data(),
                               TString::Format("%s_%s_%s", sNum.c_str(), sDen.c_str(), setting.c_str()).Data(),
                               1200, 800);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 0.9);
  pad1->Draw();
  c_obs->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.4);
  pad2->Draw();
  pad1->cd();
  auto legend = CreateLegend(0.25, 0.45, 0.6, 0.9, "", 0.05);

  int ptNum = ptBins.size()-1;
  ptNum = 3;
  for (int i = 0; i < ptNum; i++){
    double min_pt = ptBins[i];
    double max_pt = ptBins[i+1];

    hRatio = make_ratio(hNum, sNum, hDen, sDen, min_pt, max_pt);
    hRatio->SetMarkerStyle(GetMarker(i));
    hRatio->SetMarkerColor(GetColor(i));
    hRatio->SetLineColor(GetColor(i));
    legend->AddEntry(hRatio, TString::Format("#it{p}_{T}^{ch jet} %.0f - %.0f", min_pt, max_pt).Data());
    hists.push_back(hRatio);

    ptRatio = (TH1F*)hRatio->Clone("bin/base");
    ptRatio->Divide(hists[0]);
    ratios.push_back(ptRatio);
  }

  string yTitle;
  if (sDen == "pi" || sDen == "pi0"){
    yTitle = TString::Format("%s/#%s (%s)", sNum.c_str(), sDen.c_str(), setting.c_str()).Data();
  }
  else{
    yTitle = TString::Format("%s/%s (%s)", sNum.c_str(), sDen.c_str(), setting.c_str()).Data();
  }
  double xmin = 0, xmax = 0, ymin = 0, ymax = 0, yminr = 0, ymaxr = 0;
  for (auto hist : hists){
    int lastBin = hist->FindLastBinAbove(0,1);
    xmax = max(xmax, hist->GetBinLowEdge(lastBin) + hist->GetBinWidth(lastBin));
  }

  /* If you want to manually set the range, do it here.
  xmin = 5;
  xmax = 15;
  // */

  for (auto hist : hists){
    hist->GetXaxis()->SetRangeUser(xmin, xmax);
    ymax = max(ymax, hist->GetMaximum());
  }
  for (auto ratio : ratios){
    ratio->GetXaxis()->SetRangeUser(xmin, xmax);
    ymaxr = max(ymax, ratio->GetMaximum());
  }
  if (ymax > 3) ymax = 3;
  if (ymaxr > 3) ymaxr = 3;

  TH1F* frame1 = DrawFrame(xmin, xmax, ymin, ymax, "", yTitle);
  frame1->Draw();
  for (auto hist : hists){
    hist->Draw("same");
  }
  legend->Draw("same");

  pad2->cd();
  TH1F* frame2 = DrawFrame(xmin, xmax, yminr, ymaxr, "#it{p}_{T}^{ had}", "Ratio");
  frame2->Draw("same");
  for (auto ratio : ratios){
    ratio->Draw("same");
  }
  c_obs->cd();
  c_obs->SaveAs(TString::Format("../plots/fragmentation/HadronRatios/%s_%s_%s_ratio_%.0f_%.0f.pdf",
                                setting.c_str(), sNum.c_str(), sDen.c_str(), xmin, xmax).Data());
  // delete c_obs;
}

TH1F* make_ratio(TH2F* hNum, string sNum, TH2F* hDen, string sDen, double min_pt, double max_pt){
  hNum->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hNumerator = (TH1F*)hNum->ProjectionX();
  hDen->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hDenominator = (TH1F*)hDen->ProjectionX();

  // Top plot settings
  hNumerator->SetStats(0);
  hNumerator->SetTitle(TString::Format("%s/%s, pt_{ch} #in [%.0f,%.0f) GeV/c",
                                       sNum.c_str(), sDen.c_str(), min_pt, max_pt).Data());
  hNumerator->GetYaxis()->SetTitle(TString::Format("%s/%s", sNum.c_str(), sDen.c_str()).Data());
  hDenominator->SetStats(0);

  // Ratio plot settings
  TH1F *ratio = (TH1F*)hNumerator->Clone("Ratio baryon/meson");
  ratio->Divide(hDenominator);
  ratio->GetXaxis()->SetTitle(TString::Format("#it{p}_{T}").Data());
  return ratio;
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

// TFitResultPtr* fit_spectrum(TH1F* hist, xmin, xmax){
//   TFitResultPtr fit = hist->Fit("expo", "S", "", xmin, xmax);
//   return fit;
// }
