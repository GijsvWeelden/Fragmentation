
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
// Plot fragmentation in pp, AAnr, AAr with and without
// background subtraction
//
//-------------------------------------------------------------

void plot_hists_ppAA(TH2F *hB, TH2F *h1, TH2F *h2, TH2F *h3, string obs, double min_pt, double max_pt, string setting, string type, string suff);
void plot_hists_ILC(TH2F *hB, TH2F *h1, TH2F *h2, string obs, double min_pt, double max_pt, string setting, string type, string suff); // Inclusive, Leading, Awayside

void plot_frag(void){
  double time = clock();
  gROOT->SetBatch();
  string inName = "2dhists_frag";
  string sNN = "5tev02";
  string jetType = "charged";
  TFile *inFile = TFile::Open(TString::Format("./%s_%s_%s.root", inName.c_str(), sNN.c_str(), jetType.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  std::vector<double> ptBins = {20.,40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> obs = {"frag"};//,"orth"};

  TList *ppList = (TList*)inFile->Get("pp");
  if (!ppList){
    std::cout << "Error: pp not found! Is necessary for baseline. Aborting program." << std::endl;
    return;
  }
  TList *AAnrList = (TList*)inFile->Get("AAnr");
  if (!AAnrList) std::cout << "Warning: AAnr not found!" << std::endl;
  TList *AAr_nobkgList = (TList*)inFile->Get("AAr_nobkg");
  if (!AAr_nobkgList) cout << "Warning: AAr_nobkg not found!" << std::endl;
  TList *AAr_bkgList = (TList*)inFile->Get("AAr_bkg");
  if (!AAr_bkgList) cout << "Warning: AAr_bkg not found!" << std::endl;

  //TH2F *hpp, *hpp_L, *hpp_A, *hAAnr, *hAAnr_L, *hAAnr_A, *hAAr, *hAAr_L, *hAAr_A;
  TH2F *hpp, *hAAnr, *hAAr_nobkg, *hAAr_bkg;

  for (int iobs=0; iobs<obs.size(); iobs++){
    std::cout << "Plotting " << obs[iobs] << std::endl;
    hpp = (TH2F*)ppList->FindObject(TString::Format("pp_%s", obs[iobs].c_str()).Data());
    //hpp_L = (TH2F*)ppList->FindObject(TString::Format("pp_%s_leading",
    //                                                  obs[iobs].c_str()).Data());
    //hpp_A = (TH2F*)ppList->FindObject(TString::Format("pp_%s_awayside",
    //                                                  obs[iobs].c_str()).Data());
    if (AAnrList){
      hAAnr = (TH2F*)AAnrList->FindObject(TString::Format("AAnr_%s", obs[iobs].c_str()).Data());
      //hAAnr_L = (TH2F*)AAnrList->FindObject(TString::Format("AA_norecoil_%s_leading",
      //                                                        obs[iobs].c_str()).Data());
      //hAAnr_A = (TH2F*)AAnrList->FindObject(TString::Format("AA_norecoil_%s_awayside",
      //                                                        obs[iobs].c_str()).Data());
    }
    if (AAr_nobkgList){
      hAAr_nobkg = (TH2F*)AAr_nobkgList->FindObject(TString::Format("AAr_nobkg_%s", obs[iobs].c_str()).Data());
      //hAAr_L = (TH2F*)AAr_nobkgList->FindObject(TString::Format("AAr_nobkg_%s_leading",
      //                                                      obs[iobs].c_str()).Data());
      //hAAr_A = (TH2F*)AAr_nobkgList->FindObject(TString::Format("AAr_nobkg_%s_awayside",
      //                                                      obs[iobs].c_str()).Data());
    }
    if (AAr_bkgList){
      hAAr_bkg = (TH2F*)AAr_bkgList->FindObject(TString::Format("AAr_bkg_%s", obs[iobs].c_str()).Data());
    }

    plot_hists_ppAA(hpp, hAAnr, hAAr_nobkg, hAAr_bkg, obs[iobs], ptBins[0], ptBins.back(), "", "", jetType);
    for (int ipt=0; ipt<ptBins.size()-1; ++ipt){
      plot_hists_ppAA(hpp, hAAnr, hAAr_nobkg, hAAr_bkg, obs[iobs], ptBins[ipt], ptBins[ipt+1], "", "", jetType);
    }

    /*
    plot_hists_ppAA(hpp, hAAnr, hAAr, obs[iobs], ptBins[0], ptBins.back(), "", "", suffix);
    plot_hists_ppAA(hpp_L, hAAnr_L, hAAr_L, obs[iobs], ptBins[0], ptBins.back(), "", "leading", suffix);
    plot_hists_ppAA(hpp_A, hAAnr_A, hAAr_A, obs[iobs], ptBins[0], ptBins.back(), "", "away", suffix);
    plot_hists_ppAA(hpp, hpp_L, hpp_A, obs[iobs], ptBins[0], ptBins.back(), "pp", "all", suffix);
    if (AAnrList) plot_hists_ppAA(hAAnr, hAAnr_L, hAAnr_A, obs[iobs], ptBins[0], ptBins.back(), "AAnr", "all", suffix);
    if (AAr_bkgList) plot_hists_ppAA(hAAr, hAAr_L, hAAr_A, obs[iobs], ptBins[0], ptBins.back(), "AAr_bkg", "all", suffix);
    for (int ipt=0; ipt<ptBins.size()-1; ipt++){
      std::cout << "pt bin: " << ptBins[ipt] << "-" << ptBins[ipt+1] << std::endl;
      plot_hists_ppAA(hpp, hAAnr, hAAr, obs[iobs], ptBins[ipt], ptBins[ipt+1], "", "", suffix);
      plot_hists_ppAA(hpp_L, hAAnr_L, hAAr_L, obs[iobs], ptBins[ipt], ptBins[ipt+1], "", "leading", suffix);
      plot_hists_ppAA(hpp_A, hAAnr_A, hAAr_A, obs[iobs], ptBins[ipt], ptBins[ipt+1], "", "away", suffix);
      plot_hists_ppAA(hpp, hpp_L, hpp_A, obs[iobs], ptBins[ipt], ptBins[ipt+1], "pp", "all", suffix);
      if (AAnrList) plot_hists_ppAA(hAAnr, hAAnr_L, hAAnr_A, obs[iobs], ptBins[ipt], ptBins[ipt+1], "AAnr", "all", suffix);
      if (AAr_bkgList) plot_hists_ppAA(hAAr, hAAr_L, hAAr_A, obs[iobs], ptBins[ipt], ptBins[ipt+1], "AAr_bkg", "all", suffix);
    }
    */
  }
  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  // Do not close the file, it will delete the ratio histogram
  inFile->Close();
}

void plot_hists_ppAA(TH2F *hB, TH2F *h1, TH2F *h2, TH2F *h3, string obs, double min_pt, double max_pt, string setting, string type, string suff){
  double jetR = 0.4;
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
  h3->GetYaxis()->SetRangeUser(min_pt, max_pt);
  TH1F *hZ = (TH1F*)h3->ProjectionX();
  hZ->Scale(1./hZ->Integral()); // Check normalisation within precision
  if (abs(hBase->Integral() -  1.0) > 0.01
      ||
      abs(hX->Integral() - 1.0) > 0.01
      ||
      abs(hY->Integral() - 1.0) > 0.01
      ||
      abs(hZ->Integral() - 1.0) > 0.01){
    cout << "WARNING: NORMALISATION PROBLEM!" << endl
      << TString::Format("%s for pt %.0f-%.0f GeV/c", obs.c_str(), min_pt, max_pt).Data() << endl
      << "Base = " << hBase->Integral() << " (" << hBase->Integral()-1.0  << ")" << endl
      << "X = " << hX->Integral() << " (" << hX->Integral()-1.0 << ")" << endl
      << "Y = " << hY->Integral() << " (" << hY->Integral()-1.0 << ")" << endl
      << "Z = " << hZ->Integral() << " (" << hZ->Integral()-1.0 << ")" << endl;
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
  hBase->GetYaxis()->SetTitleOffset(1.5);//2.5);
  hBase->GetYaxis()->SetLabelFont(43);
  hBase->GetYaxis()->SetLabelSize(14);
  //hBase->GetYaxis()->SetLabelOffset(0.001);

  hX->SetStats(0);
  hX->SetMarkerStyle(kMultiply);
  hX->SetLineColor(kGreen+3); //kRed);
  hX->SetMarkerColor(kGreen+3); //kRed);

  hY->SetStats(0);
  hY->SetMarkerStyle(kStar);
  hY->SetLineColor(kMagenta);
  hY->SetMarkerColor(kMagenta);

  hZ->SetStats(0);
  hZ->SetMarkerStyle(kCircle);
  hZ->SetLineColor(kOrange+7);
  hZ->SetMarkerColor(kOrange+7);

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
  ratioBase->GetYaxis()->SetTitleOffset(1.5);
  ratioBase->GetYaxis()->SetLabelFont(43);
  ratioBase->GetYaxis()->SetLabelSize(14);

  TH1F *ratioX = (TH1F*)hX->Clone("Ratio X/Base");
  ratioX->Divide(hBase);

  TH1F *ratioY = (TH1F*)hY->Clone("Ratio Y/Base");
  ratioY->Divide(hBase);

  TH1F *ratioZ = (TH1F*)hZ->Clone("Ratio Z/Base");
  ratioZ->Divide(hBase);

  // Set plot ranges TODO: Change these to make more sensee
  int XmaxBinR = ratioBase->FindLastBinAbove(0,1);
  int XminBinR = ratioBase->FindFirstBinAbove(0,1);
  double Ymin = 0;
  double Ymax = max({hBase->GetMaximum(),
      hX->GetMaximum(),
      hY->GetMaximum(),
      hZ->GetMaximum()});
  double YminR = min({ratioBase->GetMinimum(),
      ratioX->GetMinimum(),
      ratioY->GetMinimum(),
      ratioZ->GetMinimum()});
  double YmaxR = max({ratioBase->GetMaximum(),
      ratioX->GetMaximum(),
      ratioY->GetMaximum(),
      ratioZ->GetMaximum()});

  hBase->GetXaxis()->SetRange(0,XmaxBinR);
  ratioBase->GetXaxis()->SetRange(0,XmaxBinR);
  if (YmaxR > 3) YmaxR = 3;
  else YmaxR *= 1.1;
  if (obs == "frag") Ymin = 1e-6;
  else if (obs == "orth") Ymin = 1e-11;

  hBase->GetYaxis()->SetRangeUser(Ymin, 1);//1.1*Ymax);
  ratioBase->GetYaxis()->SetRangeUser(0.9*YminR,YmaxR);

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

  hBase->Draw();
  if (hX->GetEntries() != 0) hX->Draw("same");
  if (hY->GetEntries() != 0) hY->Draw("same");
  if (hZ->GetEntries() != 0) hZ->Draw("same");
  if (type == "all"){
    hBase->Reset(); // Only use for ratio
  }

  auto legend = new TLegend(0.6,0.7,0.9,0.9); // Top right corner

  if (type == "all"){
    legend->AddEntry(hBase,"Full set");
    if (hX->GetEntries() != 0) legend->AddEntry(hX,"Leading");
    if (hY->GetEntries() != 0) legend->AddEntry(hY,"Awayside");
  }
  else{
    legend->AddEntry(hBase,"pp");
    if (hX->GetEntries() != 0) legend->AddEntry(hX, "AA no recoil");
    if (hY->GetEntries() != 0) legend->AddEntry(hY, "AA recoil no bkg");
    if (hZ->GetEntries() != 0) legend->AddEntry(hZ, "AA recoil bkg");
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

  ratioBase->Draw();
  if (ratioX->GetEntries() != 0) ratioX->Draw("same");
  if (ratioY->GetEntries() != 0) ratioY->Draw("same");
  if (ratioZ->GetEntries() != 0) ratioZ->Draw("same");
  c_obs->cd();

  //c_obs->SaveAs(TString::Format("../plots/%s_pt_%.0f_%.0f.pdf", obs.c_str(), min_pt, max_pt).Data());
  c_obs->SaveAs("./myPlot.pdf");
}

