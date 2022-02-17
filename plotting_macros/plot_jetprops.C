
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

void plot_histograms(TH1F *h_pp, TH1F *h_AA_nr, TH1F *h_AA_r, string obs, double min_pt, double max_pt){
  double jetR = 0.4;
  h_pp->Scale(1./h_pp->Integral());
  h_AA_nr->Scale(1./h_AA_nr->Integral());
  h_AA_r->Scale(1./h_AA_r->Integral());

  // Check normalisation within precision
  if (abs(h_pp->Integral() -  1.0) > 0.01
      ||
      abs(h_AA_nr->Integral() - 1.0) > 0.01
      ||
      abs(h_AA_r->Integral() - 1.0) > 0.01){
    cout << "WARNING: NORMALISATION PROBLEM!" << endl
      << TString::Format("%s for pt %.0f-%.0f GeV/c", obs.c_str(), min_pt, max_pt).Data() << endl
      << "pp = " << h_pp->Integral() << " (" << h_pp->Integral()-1.0  << ")" << endl
      << "AA_norecoil = " << h_AA_nr->Integral() << " (" << h_AA_nr->Integral()-1.0 << ")" << endl
      << "AA_recoil = " << h_AA_r->Integral() << " (" << h_AA_r->Integral()-1.0 << ")" << endl;
  }

  // Top plot settings
  h_pp->SetStats(0);
  h_pp->SetMarkerStyle(kPlus);
  h_pp->SetLineColor(kBlue);
  h_pp->SetMarkerColor(kBlue);
  h_pp->SetTitle(TString::Format("Jet %s, pt #in [%.0f,%.0f) GeV/c, R = %.1f",
        obs.c_str(), min_pt, max_pt, jetR).Data());
  //h_pp->GetXaxis()->SetRange(0,XmaxBin);
  h_pp->GetYaxis()->SetTitle(TString::Format("#frac{1}{N_{jets}} #frac{dN}{d %s}",
        obs.c_str()).Data());
  h_pp->GetYaxis()->SetTitleFont(43);
  h_pp->GetYaxis()->SetTitleSize(20);
  h_pp->GetYaxis()->SetTitleOffset(2.5);
  h_pp->GetYaxis()->SetLabelFont(43);
  h_pp->GetYaxis()->SetLabelSize(14);

  h_AA_nr->SetStats(0);
  h_AA_nr->SetMarkerStyle(kMultiply);
  h_AA_nr->SetLineColor(kRed);
  h_AA_nr->SetMarkerColor(kRed);

  h_AA_r->SetStats(0);
  h_AA_r->SetMarkerStyle(kStar);
  h_AA_r->SetLineColor(kMagenta);
  h_AA_r->SetMarkerColor(kMagenta);

  // Ratio plot settings
  TH1F *ratio_pp = (TH1F*)h_pp->Clone("Ratio pp/pp");
  ratio_pp->Divide(h_pp);
  ratio_pp->SetTitle("");
  // TODO: Add units?
  ratio_pp->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());
  ratio_pp->GetXaxis()->SetTitleFont(43);
  ratio_pp->GetXaxis()->SetTitleSize(20);
  ratio_pp->GetXaxis()->SetTitleOffset(2.5);
  ratio_pp->GetXaxis()->SetLabelFont(43);
  ratio_pp->GetXaxis()->SetLabelSize(14);
  ratio_pp->GetYaxis()->SetTitle("Ratio");
  ratio_pp->GetYaxis()->SetTitleFont(43);
  ratio_pp->GetYaxis()->SetTitleSize(20);
  ratio_pp->GetYaxis()->SetTitleOffset(2.5);
  ratio_pp->GetYaxis()->SetLabelFont(43);
  ratio_pp->GetYaxis()->SetLabelSize(14);

  TH1F *ratio_AA_nr = (TH1F*)h_AA_nr->Clone("Ratio AA_norecoil/pp");
  ratio_AA_nr->Divide(h_pp);

  TH1F *ratio_AA_r = (TH1F*)h_AA_r->Clone("Ratio AA_recoil/pp");
  ratio_AA_r->Divide(h_pp);

  // Set plot ranges
  int XmaxBinR = ratio_pp->FindLastBinAbove(0,1);
  double Ymax = max({h_pp->GetMaximum(),
      h_AA_nr->GetMaximum(),
      h_AA_r->GetMaximum()});
  double YminR = min({ratio_pp->GetMinimum(),
      ratio_AA_nr->GetMinimum(),
      ratio_AA_r->GetMinimum()});
  double YmaxR = max({ratio_pp->GetMaximum(),
      ratio_AA_nr->GetMaximum(),
      ratio_AA_r->GetMaximum()});


  if (obs == "dphi"){
    h_pp->GetXaxis()->SetRangeUser(0.1,h_pp->GetBinLowEdge(XmaxBinR));
    ratio_pp->GetXaxis()->SetRangeUser(0.1,h_pp->GetBinLowEdge(XmaxBinR));
    h_AA_nr->GetXaxis()->SetRangeUser(0.1,h_pp->GetBinLowEdge(XmaxBinR));
    h_AA_r->GetXaxis()->SetRangeUser(0.1,h_pp->GetBinLowEdge(XmaxBinR));
    Ymax = max({h_pp->GetMaximum(),
        h_AA_nr->GetMaximum(),
        h_AA_r->GetMaximum()});
  }
  else if (obs == "zg"){
    h_pp->GetXaxis()->SetRangeUser(0.1,h_pp->GetBinLowEdge(XmaxBinR));
    ratio_pp->GetXaxis()->SetRangeUser(0.1,h_pp->GetBinLowEdge(XmaxBinR));
    ratio_AA_nr->GetXaxis()->SetRangeUser(0.1,h_pp->GetBinLowEdge(XmaxBinR));
    ratio_AA_r->GetXaxis()->SetRangeUser(0.1,h_pp->GetBinLowEdge(XmaxBinR));
    YminR = min({ratio_pp->GetMinimum(),
        //ratio_AA_nr->GetMinimum(),
        ratio_AA_r->GetMinimum()});
    YmaxR = max({ratio_pp->GetMaximum(),
        //ratio_AA_nr->GetMaximum(),
        ratio_AA_r->GetMaximum()});
  }
  else if (obs == "mr2"){
    h_pp->GetXaxis()->SetRangeUser(0,0.3);
    ratio_pp->GetXaxis()->SetRangeUser(0,0.3);
  }
  else if (obs == "mz2"){
    h_pp->GetXaxis()->SetRangeUser(0,0.2);
    ratio_pp->GetXaxis()->SetRangeUser(0,0.2);
  }
  else{
    h_pp->GetXaxis()->SetRange(0,XmaxBinR);
    ratio_pp->GetXaxis()->SetRange(0,XmaxBinR);
  }
  h_pp->GetYaxis()->SetRangeUser(0,1.1*Ymax);
  ratio_pp->GetYaxis()->SetRangeUser(0.9*YminR,1.1*YmaxR);

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

  h_pp->Draw();
  if (h_AA_nr->GetEntries() != 0) h_AA_nr->Draw("same");
  if (h_AA_r->GetEntries() != 0)  h_AA_r->Draw("same");

  // if (obs != "dphi"){
  auto legend = new TLegend(0.7,0.8,0.9,0.9); // Top right corner
  if (obs == "dphi"){
    delete legend;
    auto legend = new TLegend(0.15,0.8,0.35,0.9); // Top left corner
  }
  legend->AddEntry(h_pp,"pp");
  if (h_AA_nr->GetEntries() != 0) legend->AddEntry(h_AA_nr,"AA_norecoil");
  if (h_AA_r->GetEntries() != 0) legend->AddEntry(h_AA_r,"AA_recoil");
  legend->Draw();

  c_obs->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.15);
  pad2->Draw();
  pad2->cd();

  ratio_pp->Draw();
  if (ratio_AA_nr->GetEntries() != 0) ratio_AA_nr->Draw("same");
  if (ratio_AA_r->GetEntries() != 0)  ratio_AA_r->Draw("same");
  c_obs->cd();
  c_obs->SaveAs(TString::Format("../plots/%s_pt%.0f-%.0f.pdf",
        obs.c_str(), min_pt, max_pt).Data());
}

//-------------------------------------------------------------
//
// Main function
//
//-------------------------------------------------------------

void plot_jetprops(void){
  double time = clock();

  // TODO: make user input
  string inName = "compare_2tev76_ppAAnrAAr";
  TFile *inFile = TFile::Open(TString::Format("./%s.root",inName.c_str()).Data());
  if(!inFile){
    cout << "File " << inFile << " not found. Aborting program." << endl;
    return;
  }

  std::vector<double> ptBins = {0.,20.,40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> obs = {"dphi","nconst","zg","Rg","nSD","mass","mz2","mr","mr2","rz","r2z","ptD","t2t1","t3t2","t2dist"};
  //int ipt = 3;

  for (int iobs=0; iobs<obs.size(); iobs++){
    cout << "Plotting " << obs[iobs] << endl;

    string ppName = TString::Format("pp_%s", obs[iobs].c_str()).Data();
    string AA_nrName = TString::Format("AA_norecoil_%s", obs[iobs].c_str()).Data();
    string AA_rName = TString::Format("AA_recoil_%s", obs[iobs].c_str()).Data();

    TList *ppList = (TList*)inFile->Get(TString::Format("%s",
          ppName.c_str()).Data());
    if (!ppList) cout << TString::Format("WARNING: %s not found!",
        ppName.c_str()).Data() << endl;
    TList *AA_nrList = (TList*)inFile->Get(TString::Format("%s",
          AA_nrName.c_str()).Data());
    if (!AA_nrList) cout << TString::Format("WARNING: %s not found!",
        AA_nrName.c_str()).Data() << endl;
    TList *AA_rList = (TList*)inFile->Get(TString::Format("%s",
          AA_rName.c_str()).Data());
    if (!AA_rList) cout << TString::Format("WARNING: %s not found!",
        AA_rName.c_str()).Data() << endl;

    for (int ipt=0; ipt<ptBins.size()-1; ipt++){
      TH1F *hpp = new TH1F();
      if (ppList){
        hpp = (TH1F*)ppList->FindObject(TString::Format("pp_h%s_pt%.0f_%.0f",
              obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
      }
      TH1F *hAA_nr = new TH1F();
      if (AA_nrList){
        hAA_nr = (TH1F*)AA_nrList->FindObject(TString::Format("AA_norecoil_h%s_pt%.0f_%.0f",
              obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
      }
      TH1F *hAA_r = new TH1F();
      if (AA_rList){
        hAA_r = (TH1F*)AA_rList->FindObject(TString::Format("AA_recoil_h%s_pt%.0f_%.0f",
              obs[iobs].c_str(), ptBins[ipt], ptBins[ipt+1]).Data());
      }
      plot_histograms(hpp, hAA_nr, hAA_r, obs[iobs], ptBins[ipt], ptBins[ipt+1]);
    }
  }
  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  // Do not close the file, it will delete the ratio histogram
  //inFile->Close();
}
