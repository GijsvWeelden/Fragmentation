
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

void plot_histograms(TH2F *hB, TH2F *h1, TH2F *h2, string obs, double min_pt, double max_pt, string setting, string type);

void plot_jetprops(void){
  double time = clock();
  string inName = "2dhists_2tev76_ppAAnrAAr";
  TFile *inFile = TFile::Open(TString::Format("./%s.root",inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  std::vector<double> ptBins = {0.,20.,40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> obs = {"dphi","nconst","zg","Rg","nSD","mass","mz2","mr","mr2","rz","r2z",
                             "ptD","t2t1","t3t2","t2dist","t3dist[0]","t3dist[1]","t3dist[2]"};

  TList *ppList = (TList*)inFile->Get("pp");
  if (!ppList){
    std::cout << "Error: pp not found! Is necessary for baseline. Aborting program." << std::endl;
    return;
  }
  TList *AA_nrList = (TList*)inFile->Get("AA_norecoil");
  if (!AA_nrList) std::cout << "Warning: AA_norecoil not found!" << std::endl;
  TList *AA_rList = (TList*)inFile->Get("AA_recoil");
  if (!AA_rList) cout << "Warning: AA_recoil not found!" << std::endl;

  TH2F *hpp, *hpp_L, *hpp_A, *hAA_nr, *hAA_nr_L, *hAA_nr_A, *hAA_r, *hAA_r_L, *hAA_r_A;

  for (int iobs=0; iobs<obs.size(); iobs++){
    std::cout << "Plotting " << obs[iobs] << std::endl;
    if (ppList){
      hpp = (TH2F*)ppList->FindObject(TString::Format("pp_%s",
                                                      obs[iobs].c_str()).Data());
      hpp_L = (TH2F*)ppList->FindObject(TString::Format("pp_%s_leading",
                                                        obs[iobs].c_str()).Data());
      hpp_A = (TH2F*)ppList->FindObject(TString::Format("pp_%s_awayside",
                                                        obs[iobs].c_str()).Data());
    }
    if (AA_nrList){
      hAA_nr = (TH2F*)AA_nrList->FindObject(TString::Format("AA_norecoil_%s",
                                                            obs[iobs].c_str()).Data());
      hAA_nr_L = (TH2F*)AA_nrList->FindObject(TString::Format("AA_norecoil_%s_leading",
                                                              obs[iobs].c_str()).Data());
      hAA_nr_A = (TH2F*)AA_nrList->FindObject(TString::Format("AA_norecoil_%s_awayside",
                                                              obs[iobs].c_str()).Data());
    }
    if (AA_rList){
      hAA_r = (TH2F*)AA_rList->FindObject(TString::Format("AA_recoil_%s",
                                                          obs[iobs].c_str()).Data());
      hAA_r_L = (TH2F*)AA_rList->FindObject(TString::Format("AA_recoil_%s_leading",
                                                            obs[iobs].c_str()).Data());
      hAA_r_A = (TH2F*)AA_rList->FindObject(TString::Format("AA_recoil_%s_awayside",
                                                            obs[iobs].c_str()).Data());
    }

    plot_histograms(hpp, hAA_nr, hAA_r, obs[iobs], ptBins[0], ptBins.back(), "", "");
    plot_histograms(hpp_L, hAA_nr_L, hAA_r_L, obs[iobs], ptBins[0], ptBins.back(), "", "leading");
    plot_histograms(hpp_A, hAA_nr_A, hAA_r_A, obs[iobs], ptBins[0], ptBins.back(), "", "away");
    plot_histograms(hpp, hpp_L, hpp_A, obs[iobs], ptBins[0], ptBins.back(), "pp", "all");
    plot_histograms(hAA_nr, hAA_nr_L, hAA_nr_A, obs[iobs], ptBins[0], ptBins.back(), "AA_nr", "all");
    plot_histograms(hAA_r, hAA_r_L, hAA_r_A, obs[iobs], ptBins[0], ptBins.back(), "AA_r", "all");
    for (int ipt=0; ipt<ptBins.size()-1; ipt++){
      std::cout << "pt bin: " << ptBins[ipt] << "-" << ptBins[ipt+1] << std::endl;
      plot_histograms(hpp, hAA_nr, hAA_r, obs[iobs], ptBins[ipt], ptBins[ipt+1], "", "");
      plot_histograms(hpp_L, hAA_nr_L, hAA_r_L, obs[iobs], ptBins[ipt], ptBins[ipt+1], "", "leading");
      plot_histograms(hpp_A, hAA_nr_A, hAA_r_A, obs[iobs], ptBins[ipt], ptBins[ipt+1], "", "away");
      plot_histograms(hpp, hpp_L, hpp_A, obs[iobs], ptBins[ipt], ptBins[ipt+1], "pp", "all");
      plot_histograms(hAA_nr, hAA_nr_L, hAA_nr_A, obs[iobs], ptBins[ipt], ptBins[ipt+1], "AA_nr", "all");
      plot_histograms(hAA_r, hAA_r_L, hAA_r_A, obs[iobs], ptBins[ipt], ptBins[ipt+1], "AA_r", "all");
    }
  }
  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  // Do not close the file, it will delete the ratio histogram
  //inFile->Close();
}

void plot_histograms(TH2F *hB, TH2F *h1, TH2F *h2, string obs, double min_pt, double max_pt, string setting, string type){
  double jetR = 0.4; // TODO: Should be taken from tree
  // Project with pt cut
  hB->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hBase = (TH1F*)hB->ProjectionX();
  hBase->Scale(1./hBase->Integral());
  h1->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hX = (TH1F*)h1->ProjectionX();
  hX->Scale(1./hX->Integral());
  h2->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hY = (TH1F*)h2->ProjectionX();
  hY->Scale(1./hY->Integral());
  // Check normalisation within precision
  if (abs(hBase->Integral() -  1.0) > 0.01
      ||
      abs(hX->Integral() - 1.0) > 0.01
      ||
      abs(hY->Integral() - 1.0) > 0.01){
    cout << "WARNING: NORMALISATION PROBLEM!" << endl
      << TString::Format("%s for pt %.0f-%.0f GeV/c", obs.c_str(), min_pt, max_pt).Data() << endl
      << "Base = " << hBase->Integral() << " (" << hBase->Integral()-1.0  << ")" << endl
      << "X = " << hX->Integral() << " (" << hX->Integral()-1.0 << ")" << endl
      << "Y = " << hY->Integral() << " (" << hY->Integral()-1.0 << ")" << endl;
  }

  // Top plot settings
  hBase->SetStats(0);
  hBase->SetMarkerStyle(kPlus);
  hBase->SetLineColor(kBlue);
  hBase->SetMarkerColor(kBlue);
  if (type == "all"){
    hBase->SetTitle(TString::Format("Jet %s, pt #in [%.0f,%.0f) GeV/c, R = %.1f %s",
                                    obs.c_str(), min_pt, max_pt, jetR, setting.c_str()).Data());
  }
  else{
    hBase->SetTitle(TString::Format("Jet %s, pt #in [%.0f,%.0f) GeV/c, R = %.1f %s",
                                    obs.c_str(), min_pt, max_pt, jetR, type.c_str()).Data());
  }
  hBase->GetYaxis()->SetTitle(TString::Format("#frac{1}{N_{jets}} #frac{dN}{d %s}",
                                              obs.c_str()).Data());
  hBase->GetYaxis()->SetTitleFont(43);
  hBase->GetYaxis()->SetTitleSize(20);
  hBase->GetYaxis()->SetTitleOffset(2.5);
  hBase->GetYaxis()->SetLabelFont(43);
  hBase->GetYaxis()->SetLabelSize(14);

  hX->SetStats(0);
  hX->SetMarkerStyle(kMultiply);
  hX->SetLineColor(kRed);
  hX->SetMarkerColor(kRed);

  hY->SetStats(0);
  hY->SetMarkerStyle(kStar);
  hY->SetLineColor(kMagenta);
  hY->SetMarkerColor(kMagenta);

  // Ratio plot settings
  TH1F *ratioBase = (TH1F*)hBase->Clone("Ratio Base/Base");
  ratioBase->Divide(hBase);
  ratioBase->SetTitle("");
  ratioBase->GetXaxis()->SetTitle(TString::Format("%s", obs.c_str()).Data());
  ratioBase->GetXaxis()->SetTitleFont(43);
  ratioBase->GetXaxis()->SetTitleSize(20);
  ratioBase->GetXaxis()->SetTitleOffset(2.5);
  ratioBase->GetXaxis()->SetLabelFont(43);
  ratioBase->GetXaxis()->SetLabelSize(14);
  ratioBase->GetYaxis()->SetTitle("Ratio");
  ratioBase->GetYaxis()->SetTitleFont(43);
  ratioBase->GetYaxis()->SetTitleSize(20);
  ratioBase->GetYaxis()->SetTitleOffset(2.5);
  ratioBase->GetYaxis()->SetLabelFont(43);
  ratioBase->GetYaxis()->SetLabelSize(14);

  TH1F *ratioX = (TH1F*)hX->Clone("Ratio X/Base");
  ratioX->Divide(hBase);

  TH1F *ratioY = (TH1F*)hY->Clone("Ratio Y/Base");
  ratioY->Divide(hBase);

  // Set plot ranges TODO: Change these to make more sense
  int XmaxBinR = ratioBase->FindLastBinAbove(0,1);
  int XminBinR = ratioBase->FindFirstBinAbove(0,1);
  double Ymax = max({hBase->GetMaximum(),
      hX->GetMaximum(),
      hY->GetMaximum()});
  double YminR = min({ratioBase->GetMinimum(),
      ratioX->GetMinimum(),
      ratioY->GetMinimum()});
  double YmaxR = max({ratioBase->GetMaximum(),
      ratioX->GetMaximum(),
      ratioY->GetMaximum()});

  if (obs == "dphi"){
    hBase->GetXaxis()->SetRangeUser(0.1,hBase->GetBinLowEdge(XmaxBinR));
    ratioBase->GetXaxis()->SetRangeUser(0.1,hBase->GetBinLowEdge(XmaxBinR));
    hX->GetXaxis()->SetRangeUser(0.1,hBase->GetBinLowEdge(XmaxBinR));
    hY->GetXaxis()->SetRangeUser(0.1,hBase->GetBinLowEdge(XmaxBinR));
    Ymax = max({hBase->GetMaximum(),
        hX->GetMaximum(),
        hY->GetMaximum()});
  }
  else if (obs == "zg"){
    hBase->GetXaxis()->SetRangeUser(0.1,hBase->GetBinLowEdge(XmaxBinR));
    ratioBase->GetXaxis()->SetRangeUser(0.1,hBase->GetBinLowEdge(XmaxBinR));
    ratioX->GetXaxis()->SetRangeUser(0.1,hBase->GetBinLowEdge(XmaxBinR));
    ratioY->GetXaxis()->SetRangeUser(0.1,hBase->GetBinLowEdge(XmaxBinR));
    YminR = min({ratioBase->GetMinimum(),
        ratioX->GetMinimum(),
        ratioY->GetMinimum()});
    YmaxR = max({ratioBase->GetMaximum(),
        ratioX->GetMaximum(),
        ratioY->GetMaximum()});
  }
  else if (obs == "mr2"){
    hBase->GetXaxis()->SetRangeUser(0,0.3);
    ratioBase->GetXaxis()->SetRangeUser(0,0.3);
  }
  else if (obs == "mz2"){
    hBase->GetXaxis()->SetRangeUser(0,0.2);
    ratioBase->GetXaxis()->SetRangeUser(0,0.2);
  }
  else if (obs == "ptD"){
    hBase->GetXaxis()->SetRangeUser(hBase->GetBinLowEdge(XminBinR),1.);
    ratioBase->GetXaxis()->SetRangeUser(hBase->GetBinLowEdge(XminBinR),1.);
  }
  else{
    hBase->GetXaxis()->SetRange(0,XmaxBinR);
    ratioBase->GetXaxis()->SetRange(0,XmaxBinR);
  }
  if (YmaxR > 3) YmaxR = 3;
  else YmaxR *= 1.1;
  hBase->GetYaxis()->SetRangeUser(0,1.1*Ymax);
  ratioBase->GetYaxis()->SetRangeUser(0.9*YminR,YmaxR);

  TCanvas *c_obs = new TCanvas(TString::Format("c_%s_pt_%.0f_%.0f",
                                               obs.c_str(), min_pt, max_pt).Data(),
                               TString::Format("%s_pt_%.0f_%.0f",
                                               obs.c_str(), min_pt, max_pt).Data(),600,900);
  TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
  pad1->SetBottomMargin(0);
  pad1->SetLeftMargin(0.15);
  pad1->Draw();
  pad1->cd();

  hBase->Draw();
  if (hX->GetEntries() != 0) hX->Draw("same");
  if (hY->GetEntries() != 0)  hY->Draw("same");
  if (type == "all"){
    hBase->Reset(); // Only use for ratio
  }

  auto legend = new TLegend(0.7,0.8,0.9,0.9); // Top right corner
  if (obs == "dphi" || obs == "t3t2" || obs == "t2t1"){
    delete legend;
    auto legend = new TLegend(0.15,0.8,0.35,0.9); // Top left corner
  }

  if (type == "all"){
    legend->AddEntry(hBase,"Full set");
    if (hX->GetEntries() != 0) legend->AddEntry(hX,"Leading");
    if (hY->GetEntries() != 0) legend->AddEntry(hY,"Awayside");
  }
  else{
    legend->AddEntry(hBase,"pp");
    if (hX->GetEntries() != 0) legend->AddEntry(hX,"AA_norecoil");
    if (hY->GetEntries() != 0) legend->AddEntry(hY,"AA_recoil");
  }
  legend->Draw();

  c_obs->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->SetLeftMargin(0.15);
  pad2->Draw();
  pad2->cd();

  ratioBase->Draw();
  if (ratioX->GetEntries() != 0) ratioX->Draw("same");
  if (ratioY->GetEntries() != 0)  ratioY->Draw("same");
  c_obs->cd();
  if (type == "all"){
    c_obs->SaveAs(TString::Format("../plots/%s_%s_pt%.0f-%.0f.pdf",
                                  setting.c_str(), obs.c_str(), min_pt, max_pt).Data());
  }
  else{
    c_obs->SaveAs(TString::Format("../plots/%s_%s_pt%.0f-%.0f.pdf",
                                  type.c_str(), obs.c_str(), min_pt, max_pt).Data());
  }
}
