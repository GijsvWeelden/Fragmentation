
#include <vector>
#include <iostream>

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


void plot_histograms(TH1F *h_pp, TH1F *h_AA_nr, string obs, double min_pt, double max_pt){
  double jetR = 0.4;
  h_pp->Scale(1./h_pp->Integral(), "width");
  h_AA_nr->Scale(1./h_AA_nr->Integral(), "width");

  // Check normalisation within precision
  if (abs(h_pp->Integral() -  1.0) > 0.01
      ||
      abs(h_AA_nr->Integral() - 1.0) > 0.01){
    cout << "WARNING: NORMALISATION PROBLEM!" << endl
      << "pp = " << h_pp->Integral() << " (" << h_pp->Integral()-1.0  << ")" << endl 
      << "AA_norecoil = " << h_AA_nr->Integral() << " (" << h_AA_nr->Integral()-1.0 << ")" << endl;
  }

  // Plot settings
  h_AA_nr->SetStats(0);
  h_AA_nr->SetMarkerStyle(kPlus);
  h_AA_nr->SetTitle(TString::Format("Jet %s, pt #in [%.0f,%.0f), R = %.1f",
        obs.c_str(), min_pt, max_pt, jetR).Data());
  //h_AA_nr->SetAxisRange(0,30,"X");
  // h_AA_nr y axis
  h_AA_nr->GetYaxis()->SetTitle(TString::Format("#frac{1}{N_{jets}} #frac{dN}{d %s}",
        obs.c_str()).Data());
  h_AA_nr->GetYaxis()->SetTitleFont(43);
  h_AA_nr->GetYaxis()->SetTitleSize(20);//0.05);
  h_AA_nr->GetYaxis()->SetTitleOffset(2.5);
  h_AA_nr->GetYaxis()->SetLabelFont(43);
  h_AA_nr->GetYaxis()->SetLabelSize(14);
  // h_pp y axis
  h_pp->SetStats(0);
  h_pp->SetMarkerStyle(kMultiply);
  h_pp->SetLineColor(kRed);
  h_pp->SetMarkerColor(kRed);

  // ratio
  TH1F *ratio = (TH1F*)h_AA_nr->Clone("Ratio");
  ratio->Divide(h_pp);
  ratio->SetStats(0);
  ratio->SetMarkerStyle(kDot);
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->SetTitle("");
  // ratio x axis
  // TODO: add units?
  ratio->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());
  ratio->GetXaxis()->SetTitleFont(43);
  ratio->GetXaxis()->SetTitleSize(20);
  ratio->GetXaxis()->SetTitleOffset(2.5);
  ratio->GetXaxis()->SetLabelFont(43);
  ratio->GetXaxis()->SetLabelSize(14);
  // ratio y axis
  ratio->GetYaxis()->SetTitle("PbPb/pp");
  ratio->GetYaxis()->SetTitleFont(43);
  ratio->GetYaxis()->SetTitleSize(20);
  ratio->GetYaxis()->SetTitleOffset(2.5);
  ratio->GetYaxis()->SetLabelFont(43);
  ratio->GetYaxis()->SetLabelSize(14);
  //ratio->SetAxisRange(0,1.5,"Y");

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
  h_AA_nr->Draw();
  h_pp->Draw("same");
  auto legend = new TLegend(0.8,0.8,0.9,0.9);
  legend->AddEntry(h_AA_nr,"AA_norecoil");
  legend->AddEntry(h_pp,"pp");
  legend->AddEntry(ratio,"ratio");
  legend->Draw();
  c_obs->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.15);
  pad2->Draw();
  pad2->cd();
  ratio->Draw();
  c_obs->cd();
}

//-------------------------------------------------------------
//
// Main function
//
//-------------------------------------------------------------

void plot_jetprops(void){
  double time = clock();

  // TODO: Add filename
  // TODO: make user input
  string inName = "test";
  TFile *inFile = TFile::Open(TString::Format("./%s.root",inName.c_str()).Data());
  if(!inFile){
    cout << "File " << inFile << " not found. Aborting programme." << endl;
    return;
  }

  std::vector<double> ptBins = {0.,20.,40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> obs = {"mass"}; //{"nconst","zg","Rg","nSD","mass","mz2","mr","mr2","rz","r2z"};
  //int ipt = 3;
  int iobs = 0;

  cout << "Plotting " << obs[iobs] << endl;
  TList *ppList = (TList*)inFile->Get(TString::Format("pp_%s",
        obs[iobs].c_str()).Data());
  if (!ppList) cout << TString::Format("WARNING: pp_%s not found!",
      obs[iobs].c_str()).Data() << endl; 
  TList *AA_nrList = (TList*)inFile->Get(TString::Format("AA_norecoil_%s",
        obs[iobs].c_str()).Data());
  if (!AA_nrList) cout << TString::Format("WARNING: AA_norecoil_%s not found!",
      obs[iobs].c_str()).Data() << endl; 
  TList *AA_rList = (TList*)inFile->Get(TString::Format("AA_recoil_%s",
        obs[iobs].c_str()).Data());
  if (!AA_rList) cout << TString::Format("WARNING: AA_recoil_%s not found!",
      obs[iobs].c_str()).Data() << endl;

  for (int ipt=0; ipt<ptBins.size()-1; ipt++){
    TH1F *hpp = new TH1F();
    if (ppList){
      hpp = (TH1F*)ppList->FindObject(TString::Format("pp_hmass_pt%.0f_%.0f",
            ptBins[ipt], ptBins[ipt+1]).Data());
    }
    TH1F *hAA_nr = new TH1F();
    if (AA_nrList){
      hAA_nr = (TH1F*)AA_nrList->FindObject(TString::Format("AA_norecoil_hmass_pt%.0f_%.0f",
            ptBins[ipt], ptBins[ipt+1]).Data());
    }
    TH1F *hAA_recoil = new TH1F();
    if (AA_rList){
      hAA_recoil = (TH1F*)AA_rList->FindObject(TString::Format("hAA_recoil_hmass_pt%.0f_%.0f",
            ptBins[ipt], ptBins[ipt+1]).Data());
    }
    plot_histograms(hpp, hAA_nr, obs[iobs], ptBins[ipt], ptBins[ipt+1]);
  }

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  // Do not close the file, it will delete the ratio histogram
  //inFile->Close();
}
