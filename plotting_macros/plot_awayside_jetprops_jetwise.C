
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
// plot_wayside_jetprops_jetwise.C compares the leading and
// wayside jets between JEWEL settings for various jet properties
//
//-------------------------------------------------------------

void plot_hists(TH1F *hpp, TH1F *hAAnr, TH1F *hAAr, string setting, string obs, double min_pt, double max_pt){
  double jetR = 0.4;
  hpp->Scale(1./hpp->Integral());
  hAAnr->Scale(1./hAAnr->Integral());
  hAAr->Scale(1./hAAr->Integral());

  // Check normalisation within precision
  if (abs(hpp->Integral() -  1.0) > 0.01
      ||
      abs(hAAnr->Integral() - 1.0) > 0.01
      ||
      abs(hAAr->Integral() - 1.0) > 0.01){
    cout << "WARNING: NORMALISATION PROBLEM!" << endl
      << TString::Format("%s, %s for pt %.0f-%.0f GeV/c",
          setting.c_str(), obs.c_str(), min_pt, max_pt).Data() << endl
      << "pp = " << hpp->Integral() << " (" << hpp->Integral()-1.0  << ")" << endl
      << "AA_norecoil = " << hAAnr->Integral() << " (" << hAAnr->Integral()-1.0 << ")" << endl
      << "AA_recoil = " << hAAr->Integral() << " (" << hAAr->Integral()-1.0 << ")" << endl;
  }

  // Top plot settings
  hpp->SetStats(0);
  hpp->SetMarkerStyle(kPlus);
  hpp->SetLineColor(kBlue);
  hpp->SetMarkerColor(kBlue);
  hpp->SetTitle(TString::Format("%s jet %s, pt #in [%.0f,%.0f) GeV/c, R = %.1f",
        setting.c_str(), obs.c_str(), min_pt, max_pt, jetR).Data()); // TODO: add setting in title
  //hpp->GetXaxis()->SetRange(0,XmaxBin);
  hpp->GetYaxis()->SetTitle(TString::Format("#frac{1}{N_{jets}} #frac{dN}{d %s}",
        obs.c_str()).Data());
  hpp->GetYaxis()->SetTitleFont(43);
  hpp->GetYaxis()->SetTitleSize(20);
  hpp->GetYaxis()->SetTitleOffset(2.5);
  hpp->GetYaxis()->SetLabelFont(43);
  hpp->GetYaxis()->SetLabelSize(14);

  hAAnr->SetStats(0);
  hAAnr->SetMarkerStyle(kMultiply);
  hAAnr->SetLineColor(kRed);
  hAAnr->SetMarkerColor(kRed);

  hAAr->SetStats(0);
  hAAr->SetMarkerStyle(kStar);
  hAAr->SetLineColor(kMagenta);
  hAAr->SetMarkerColor(kMagenta);

  // Ratio plot settings
  TH1F *rpp = (TH1F*)hpp->Clone("Ratio pp/pp");
  rpp->Divide(hpp);
  rpp->SetTitle("");

  rpp->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());
  rpp->GetXaxis()->SetTitleFont(43);
  rpp->GetXaxis()->SetTitleSize(20);
  rpp->GetXaxis()->SetTitleOffset(2.5);
  rpp->GetXaxis()->SetLabelFont(43);
  rpp->GetXaxis()->SetLabelSize(14);
  rpp->GetYaxis()->SetTitle("Ratio");
  rpp->GetYaxis()->SetTitleFont(43);
  rpp->GetYaxis()->SetTitleSize(20);
  rpp->GetYaxis()->SetTitleOffset(2.5);
  rpp->GetYaxis()->SetLabelFont(43);
  rpp->GetYaxis()->SetLabelSize(14);

  TH1F *rAAnr = (TH1F*)hAAnr->Clone("Ratio AAnr/pp");
  rAAnr->Divide(hpp);

  TH1F *rAAr = (TH1F*)hAAr->Clone("Ratio AAr/pp");
  rAAr->Divide(hpp);

  // Set plot ranges
  int XmaxBin = max({hpp->FindLastBinAbove(0,1),
      hAAnr->FindLastBinAbove(0,1),
      hAAr->FindLastBinAbove(0,1)});
  double Ymax = max({hpp->GetMaximum(),
      hAAnr->GetMaximum(),
      hAAr->GetMaximum()});
  double YminR = min({rpp->GetMinimum(),
      rAAnr->GetMinimum(),
      rAAr->GetMinimum()});
  double YmaxR = max({rpp->GetMaximum(),
      rAAnr->GetMaximum(),
      rAAr->GetMaximum()});

  if (obs == "zg"){
    hpp->GetXaxis()->SetRangeUser(0.1,hpp->GetBinLowEdge(XmaxBin));
    rpp->GetXaxis()->SetRangeUser(0.1,hpp->GetBinLowEdge(XmaxBin));
    rAAnr->GetXaxis()->SetRangeUser(0.1,hpp->GetBinLowEdge(XmaxBin));
    rAAr->GetXaxis()->SetRangeUser(0.1,hpp->GetBinLowEdge(XmaxBin));
    YminR = min({rpp->GetMinimum(),
        rAAnr->GetMinimum(),
        rAAr->GetMinimum()});
    YmaxR = max({rpp->GetMaximum(),
        rAAnr->GetMaximum(),
        rAAr->GetMaximum()});
  }
  else if (obs == "mr2"){
    hpp->GetXaxis()->SetRangeUser(0,0.3);
    rpp->GetXaxis()->SetRangeUser(0,0.3);
  }
  else if (obs == "mz2"){
    hpp->GetXaxis()->SetRangeUser(0,0.2);
    rpp->GetXaxis()->SetRangeUser(0,0.2);
  }
  else{
    hpp->GetXaxis()->SetRange(0,XmaxBin);
    rpp->GetXaxis()->SetRange(0,XmaxBin);
  }
  hpp->GetYaxis()->SetRangeUser(0,1.1*Ymax);
  YmaxR = min(YmaxR, 5.);
  rpp->GetYaxis()->SetRangeUser(0.9*YminR,1.1*YmaxR);

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

  hpp->Draw();
  if (hAAnr->GetEntries() != 0) hAAnr->Draw("same");
  if (hAAr->GetEntries() != 0)  hAAr->Draw("same");

  auto legend = new TLegend(0.7,0.8,0.9,0.9); // Top right corner

  legend->AddEntry(hpp,"pp");
  if (hAAnr->GetEntries() != 0) legend->AddEntry(hAAnr,"No recoil");
  if (hAAr->GetEntries() != 0) legend->AddEntry(hAAr,"Recoil");
  legend->Draw();

  c_obs->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.15);
  pad2->Draw();
  pad2->cd();

  rpp->Draw();
  if (rAAnr->GetEntries() != 0) rAAnr->Draw("same");
  if (rAAr->GetEntries() != 0)  rAAr->Draw("same");
  c_obs->cd();

  if (hAAnr->GetEntries() != 0 && hAAr->GetEntries() != 0){
    c_obs->SaveAs(TString::Format("../plots/awayside/%s/ppAAnrAAr/%s_%s_ppAAnrAAr_pt%.0f_%.0f.pdf",
          setting.c_str(), setting.c_str(), obs.c_str(), min_pt, max_pt).Data());
  }
  else if (hAAnr->GetEntries() != 0){
    c_obs->SaveAs(TString::Format("../plots/awayside/%s/ppAAnr/%s_%s_ppAAnr_pt%.0f_%.0f.pdf",
          setting.c_str(), setting.c_str(), obs.c_str(), min_pt, max_pt).Data());
  }
  else if (hAAr->GetEntries() != 0){
    c_obs->SaveAs(TString::Format("../plots/awayside/%s/ppAAr/%s_%s_ppAAr_pt%.0f_%.0f.pdf",
          setting.c_str(), setting.c_str(), obs.c_str(), min_pt, max_pt).Data());
  }
}

//-------------------------------------------------------------
//
// Main function
//
//-------------------------------------------------------------

void plot_awayside_jetprops_jetwise(void){
  double time = clock();
  string setting = "";

  bool leading = true;
  bool away = true;

  bool pp = true; // Don't change pp
  bool AAnr = true;
  bool AAr = true;

  // File containing leading/wayside histograms
  string wayName = "awayside_2tev76_ppAAnrAAr";
  TFile *wayFile = TFile::Open(TString::Format("./%s.root",wayName.c_str()).Data());
  if(!wayFile){
    cout << "File " << wayFile << " not found. Aborting program." << endl;
    return;
  }

  std::vector<double> ptBins = {0.,20.,40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> obs = {"nconst","zg","Rg","nSD","mass","mz2","mr","mr2","rz","r2z","t2t1","t3t2","t2dist","t3dist"};

  for (int iobs=0; iobs<obs.size(); iobs++){
    cout << "Plotting " << obs[iobs] << endl;

    string ppName = TString::Format("pp_%s", obs[iobs].c_str()).Data();
    string AAnrName = TString::Format("AA_norecoil_%s", obs[iobs].c_str()).Data();
    string AArName = TString::Format("AA_recoil_%s", obs[iobs].c_str()).Data();

    TList *ppList = (TList*)wayFile->Get(TString::Format("%s",
          ppName.c_str()).Data());
    if (!ppList) cout << TString::Format("WARNING: %s not found!",
        ppName.c_str()).Data() << endl;

    TList *AAnrList = (TList*)wayFile->Get(TString::Format("%s",
          AAnrName.c_str()).Data());
    if (!AAnrList) cout << TString::Format("WARNING: %s not found!",
        AAnrName.c_str()).Data() << endl;

    TList *AArList = (TList*)wayFile->Get(TString::Format("%s",
          AArName.c_str()).Data());
    if (!AArList) cout << TString::Format("WARNING: %s/%s not found!",
        wayName.c_str(), AArName.c_str()).Data() << endl;

    // Loop over pt bins
    for (int ipt=0; ipt<ptBins.size()-1; ipt++){
      //------------
      // Leading
      //------------
      if (leading){
        setting = "Leading";
        TH1F *ppL = new TH1F();
        TH1F *AAnrL = new TH1F();
        TH1F *AArL = new TH1F();

        if (ppList && pp){
          ppL = (TH1F*)ppList->FindObject(TString::Format("pp_Leading_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        if (AAnrList && AAnr){
          AAnrL = (TH1F*)AAnrList->FindObject(TString::Format("AA_norecoil_Leading_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        if (AArList && AAr){
          AArL = (TH1F*)AArList->FindObject(TString::Format("AA_recoil_Leading_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }

        plot_hists(ppL, AAnrL, AArL, setting, obs[iobs], ptBins[ipt], ptBins[ipt+1]);
      }
      //------------
      // Awayside
      //------------
      if (away){
        setting = "Awayside";
        TH1F *ppA = new TH1F();
        TH1F *AAnrA = new TH1F();
        TH1F *AArA = new TH1F();

        if (ppList && pp){
          ppA = (TH1F*)ppList->FindObject(TString::Format("pp_Wayside_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        if (AAnrList && AAnr){
          AAnrA = (TH1F*)AAnrList->FindObject(TString::Format("AA_norecoil_Wayside_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        if (AArList && AAr){
          AArA = (TH1F*)AArList->FindObject(TString::Format("AA_recoil_Wayside_%s_pt%.0f_%.0f",
                obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
        }
        plot_hists(ppA, AAnrA, AArA, setting, obs[iobs], ptBins[ipt], ptBins[ipt+1]);
      }
    }
  }
  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  // Do not close the file, it will delete the ratio histogram
  //inFile->Close();
}
