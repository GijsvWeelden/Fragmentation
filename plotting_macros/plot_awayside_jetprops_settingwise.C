
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
// plot_wayside_jetprops.C compares the contribution of
// leading and wayside jets to the full jet spectrum for
// various jet properties
// it also shows the ratio of the full and leading spectra
// to the wayside spectrum
//
//-------------------------------------------------------------

void plot_histograms(TH1F *hA, TH1F *hL, TH1F *hF, string setting, string obs, double min_pt, double max_pt, bool compare_to_full){
  double jetR = 0.4;
  hA->Scale(1./hA->Integral());
  hL->Scale(1./hL->Integral());
  hF->Scale(1./hF->Integral());

  // Check normalisation within precision
  if (abs(hA->Integral() -  1.0) > 0.01
      ||
      abs(hL->Integral() - 1.0) > 0.01
      ||
      abs(hF->Integral() - 1.0) > 0.01){
    cout << "WARNING: NORMALISATION PROBLEM!" << endl
      << TString::Format("%s for pt %.0f-%.0f GeV/c", obs.c_str(), min_pt, max_pt).Data() << endl
      << "Wayside = " << hA->Integral() << " (" << hA->Integral()-1.0  << ")" << endl 
      << "Leading = " << hL->Integral() << " (" << hL->Integral()-1.0 << ")" << endl
      << "Full = " << hF->Integral() << " (" << hF->Integral()-1.0 << ")" << endl;
  }

  // Top plot settings
  hA->SetStats(0);
  hA->SetMarkerStyle(kPlus);
  hA->SetLineColor(kBlue);
  hA->SetMarkerColor(kBlue);
  hA->SetTitle(TString::Format("Jet %s, pt #in [%.0f,%.0f) GeV/c, R = %.1f (%s)",
        obs.c_str(), min_pt, max_pt, jetR, setting.c_str()).Data());
  //hA->GetXaxis()->SetRange(0,XmaxBin);
  hA->GetYaxis()->SetTitle(TString::Format("#frac{1}{N_{jets}} #frac{dN}{d %s}",
        obs.c_str()).Data());
  hA->GetYaxis()->SetTitleFont(43);
  hA->GetYaxis()->SetTitleSize(20);
  hA->GetYaxis()->SetTitleOffset(2.5);
  hA->GetYaxis()->SetLabelFont(43);
  hA->GetYaxis()->SetLabelSize(14);

  hL->SetStats(0);
  hL->SetMarkerStyle(kMultiply);
  hL->SetLineColor(kRed);
  hL->SetMarkerColor(kRed);

  hF->SetStats(0);
  hF->SetMarkerStyle(kStar);
  hF->SetLineColor(kMagenta);
  hF->SetMarkerColor(kMagenta);

  // Ratio plot settings
  TH1F *rA = (TH1F*)hA->Clone("Ratio A/A");
  rA->Divide(hA);
  rA->SetTitle("");

  rA->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());
  rA->GetXaxis()->SetTitleFont(43);
  rA->GetXaxis()->SetTitleSize(20);
  rA->GetXaxis()->SetTitleOffset(2.5);
  rA->GetXaxis()->SetLabelFont(43);
  rA->GetXaxis()->SetLabelSize(14);
  rA->GetYaxis()->SetTitle("Ratio");
  rA->GetYaxis()->SetTitleFont(43);
  rA->GetYaxis()->SetTitleSize(20);
  rA->GetYaxis()->SetTitleOffset(2.5);
  rA->GetYaxis()->SetLabelFont(43);
  rA->GetYaxis()->SetLabelSize(14);

  TH1F *rL = (TH1F*)hL->Clone("Ratio L/A");
  rL->Divide(hA);

  TH1F *rF = (TH1F*)hF->Clone("Ratio F/F");
  rF->Divide(hF);

  // Set plot ranges
  int XmaxBinR = rA->FindLastBinAbove(0,1);
  double Ymax = max({hA->GetMaximum(),
      hL->GetMaximum(),
      hF->GetMaximum()});
  double YminR = min({rA->GetMinimum(),
      rL->GetMinimum(),
      rF->GetMinimum()});
  double YmaxR = max({rA->GetMaximum(),
      rL->GetMaximum(),
      rF->GetMaximum()});

  if (obs == "zg"){
    hA->GetXaxis()->SetRangeUser(0.1,hA->GetBinLowEdge(XmaxBinR));
    rA->GetXaxis()->SetRangeUser(0.1,hA->GetBinLowEdge(XmaxBinR));
    rL->GetXaxis()->SetRangeUser(0.1,hA->GetBinLowEdge(XmaxBinR));
    rF->GetXaxis()->SetRangeUser(0.1,hA->GetBinLowEdge(XmaxBinR));
    YminR = min({rA->GetMinimum(),
        rL->GetMinimum(),
        rF->GetMinimum()});
    YmaxR = max({rA->GetMaximum(),
        rL->GetMaximum(),
        rF->GetMaximum()});
  }
  else if (obs == "mr2"){
    hA->GetXaxis()->SetRangeUser(0,0.3);
    rA->GetXaxis()->SetRangeUser(0,0.3);
  }
  else if (obs == "mz2"){
    hA->GetXaxis()->SetRangeUser(0,0.2);
    rA->GetXaxis()->SetRangeUser(0,0.2);
  }
  else{
    hA->GetXaxis()->SetRange(0,XmaxBinR);
    rA->GetXaxis()->SetRange(0,XmaxBinR);
  }
  hA->GetYaxis()->SetRangeUser(0,1.1*Ymax);
  YmaxR = min(YmaxR, 5.);
  rA->GetYaxis()->SetRangeUser(0.9*YminR,1.1*YmaxR);

  // Draw histograms
  TCanvas *c_obs = new TCanvas(TString::Format("c_%s_pt_%.0f_%.0f",
        obs.c_str(), min_pt, max_pt).Data(),
      TString::Format("%s_pt_%.0f_%.0f",
        obs.c_str(), min_pt, max_pt).Data(),600,900);

  TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.15);
  pad1->Draw();
  pad1->cd();

  hA->Draw();
  if (hL->GetEntries() != 0) hL->Draw("same");
  if (hF->GetEntries() != 0)  hF->Draw("same");

  auto legend = new TLegend(0.7,0.8,0.9,0.9); // Top right corner

  legend->AddEntry(hA,"Awayside");
  if (hL->GetEntries() != 0) legend->AddEntry(hL,"Leading");
  if (hF->GetEntries() != 0) legend->AddEntry(hF,"Full");
  legend->Draw();

  c_obs->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.15);
  pad2->Draw();
  pad2->cd();

  rA->Draw();
  if (rL->GetEntries() != 0) rL->Draw("same");
  //if (rF->GetEntries() != 0)  rF->Draw("same");
  c_obs->cd();
  if (compare_to_full) c_obs->SaveAs(TString::Format("../plots/awayside/%s/ALF/%s_%s_ALF_pt%.0f_%.0f.pdf",
        setting.c_str(), setting.c_str(), obs.c_str(), min_pt, max_pt).Data());
  else c_obs->SaveAs(TString::Format("../plots/awayside/%s/AL/%s_%s_AL_pt%.0f_%.0f.pdf",
        setting.c_str(), setting.c_str(), obs.c_str(), min_pt, max_pt).Data());


}

//-------------------------------------------------------------
//
// Main function
//
//-------------------------------------------------------------

void plot_awayside_jetprops_settingwise(void){
  double time = clock();
  bool compare_to_full = false;
  string setting = "";

  bool pp = true;
  bool AAnr = true;
  bool AAr = true;

  // File containing leading/wayside histograms
  string wayName = "awayside_2tev76_ppAAnrAAr";
  TFile *wayFile = TFile::Open(TString::Format("./%s.root",wayName.c_str()).Data());
  if(!wayFile){
    cout << "File " << wayFile << " not found. Aborting program." << endl;
    return;
  }

  // File containing full histograms
  string fullName = "compare_2tev76_ppAAnrAAr";
  TFile *fullFile = TFile::Open(TString::Format("./%s.root",fullName.c_str()).Data());
  if(!fullFile){
    cout << "File " << fullFile << " not found. Aborting program." << endl;
    return;
  }

  std::vector<double> ptBins = {0.,20.,40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> obs = {"nconst","zg","Rg","nSD","mass","mz2","mr","mr2","rz","r2z"};

  for (int iobs=0; iobs<obs.size(); iobs++){
    cout << "Plotting " << obs[iobs] << endl;

    string ppName = TString::Format("pp_%s", obs[iobs].c_str()).Data();
    string AAnrName = TString::Format("AA_norecoil_%s", obs[iobs].c_str()).Data();
    string AArName = TString::Format("AA_recoil_%s", obs[iobs].c_str()).Data();

    TList *ppList = (TList*)wayFile->Get(TString::Format("%s",
          ppName.c_str()).Data());
    if (!ppList) cout << TString::Format("WARNING: %s not found!",
        ppName.c_str()).Data() << endl;
    TList *ppList_F = (TList*)fullFile->Get(TString::Format("%s",
          ppName.c_str()).Data());
    if (!ppList_F) cout << TString::Format("WARNING: %s/%s not found!",
        fullName.c_str(), ppName.c_str()).Data() << endl;

    TList *AAnrList = (TList*)wayFile->Get(TString::Format("%s",
          AAnrName.c_str()).Data());
    if (!AAnrList) cout << TString::Format("WARNING: %s not found!",
        AAnrName.c_str()).Data() << endl;
    TList *AAnrList_F = (TList*)fullFile->Get(TString::Format("%s",
          AAnrName.c_str()).Data());
    if (!AAnrList_F) cout << TString::Format("WARNING: %s/%s not found!",
        fullName.c_str(), AAnrName.c_str()).Data() << endl;

    TList *AArList = (TList*)wayFile->Get(TString::Format("%s",
          AArName.c_str()).Data());
    if (!AArList) cout << TString::Format("WARNING: %s/%s not found!",
        wayName.c_str(), AArName.c_str()).Data() << endl;
    TList *AArList_F = (TList*)fullFile->Get(TString::Format("%s",
          AArName.c_str()).Data());
    if (!AArList_F) cout << TString::Format("WARNING: %s/%s not found!",
        fullName.c_str(), AArName.c_str()).Data() << endl;

    // Loop over pt bins
    for (int ipt=0; ipt<ptBins.size()-1; ipt++){
      //------------
      // pp
      //------------
      if (pp){
        setting = "pp";
        TH1F *hpp_F = new TH1F();
        TH1F *hpp_W = new TH1F();
        TH1F *hpp_L = new TH1F();

        if (ppList){
          hpp_W = (TH1F*)ppList->FindObject(TString::Format("pp_Wayside_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
          hpp_L = (TH1F*)ppList->FindObject(TString::Format("pp_Leading_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        if (ppList_F && compare_to_full){
          hpp_F = (TH1F*)ppList_F->FindObject(TString::Format("pp_h%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        plot_histograms(hpp_W, hpp_L, hpp_F, setting, obs[iobs], ptBins[ipt], ptBins[ipt+1], compare_to_full);
      }
      //------------
      // AA_norecoil
      //------------
      if (AAnr){
        setting = "AA_norecoil";
        TH1F *hAAnr_F = new TH1F();
        TH1F *hAAnr_W = new TH1F();
        TH1F *hAAnr_L = new TH1F();

        if (AAnrList){
          hAAnr_W = (TH1F*)AAnrList->FindObject(TString::Format("AA_norecoil_Wayside_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
          hAAnr_L = (TH1F*)AAnrList->FindObject(TString::Format("AA_norecoil_Leading_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        if (AAnrList_F && compare_to_full){
          hAAnr_F = (TH1F*)AAnrList_F->FindObject(TString::Format("AA_norecoil_h%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        plot_histograms(hAAnr_W, hAAnr_L, hAAnr_F, setting, obs[iobs], ptBins[ipt], ptBins[ipt+1], compare_to_full);
      }
      //------------
      // AA_recoil
      //------------
      if (AAr){
        setting = "AA_recoil";
        TH1F *hAAr_F = new TH1F();
        TH1F *hAAr_W = new TH1F();
        TH1F *hAAr_L = new TH1F();

        if (AArList){
          hAAr_W = (TH1F*)AArList->FindObject(TString::Format("AA_recoil_Wayside_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
          hAAr_L = (TH1F*)AArList->FindObject(TString::Format("AA_recoil_Leading_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        if (AArList_F && compare_to_full){
          hAAr_F = (TH1F*)AArList_F->FindObject(TString::Format("AA_recoil_h%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        plot_histograms(hAAr_W, hAAr_L, hAAr_F, setting, obs[iobs], ptBins[ipt], ptBins[ipt+1], compare_to_full);
      }
    }
  }
  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  // Do not close the file, it will delete the ratio histogram
  //inFile->Close();
}
