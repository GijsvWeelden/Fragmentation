
#include <vector>
#include <iostream>
#include <typeinfo>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TLegend.h"

#include "./plotUtils.C"

//-------------------------------------------------------------
//
// Plot pythia fragmentation
//
//-------------------------------------------------------------

string format_hadron_name(string hadron);
void prep_single_hadron(TH2F* hInclJet, TH2F* hgJet, TH2F* hqJet, TH2F* hNJetTypes,
                        TH1F* hIncl, TH1F* hg, TH1F* hq, TH1F* hNgq,
                        double min_pt, double max_pt
                        );
void plot_FFs_single(TH1F* hIncl, TH1F* hg, TH1F* hq, string hadron, double min_pt, double max_pt);
void plot_FFs_mult(TH1F* h0, TH1F* h1, TH1F* h2, TH1F* h3,
                   string hadron0, string hadron1, string hadron2, string hadron3,
                   double min_pt, double max_pt, string IGQ
                   );
void plot_matchDist(TH1F* dist, std::vector<double> matchDist);

void plot_frags(void)
{
  double time = clock();
  gROOT->SetBatch();
  bool plotSingleHadron = false, plotMultHadron = true, plotMatchDist = false;
  string inName = "PythiaResultJob12128342_66_pthat80_200";
  TFile *inFile = TFile::Open(TString::Format("../data/%s.root", inName.c_str()).Data());
  if(!inFile){
    std::cout << "File " << inFile << " not found. Aborting program." << std::endl;
    return;
  }
  std::vector<double> jetPtBins = {40.,60.,80.,100.,120.,160.,200.};
  std::vector<string> Hadrons = {"pi", "K", "p", "pi0", "K0L", "K0S", "K0", "Lambda0"};
	std::vector<string> Partons = {"g", "q"};

  TH2F* hNJetTypes = (TH2F*) inFile->Get("hNJetTypes");
  TH1F* hNgqJets = (TH1F*) hNJetTypes->ProjectionX();
  int nGluons_all = (int) hNgqJets->GetBinContent(1);
  int nQuarks_all = (int) hNgqJets->GetBinContent(2);
  int nJets = nGluons_all + nQuarks_all;

  if (plotSingleHadron){
    TH2F *hInclJet, *hgJet, *hqJet;
    TH1F *hIncl, *hg, *hq, *hNgq;

    for (int iHad = 0; iHad < Hadrons.size(); iHad++){
      std::cout << "Plotting " << Hadrons[iHad] << std::endl;
      hInclJet = (TH2F*) inFile->Get(TString::Format("hJetFrag_%s", Hadrons[iHad].c_str()).Data());
      hgJet = (TH2F*) inFile->Get(TString::Format("hgFrags_%s", Hadrons[iHad].c_str()).Data());
      hqJet = (TH2F*) inFile->Get(TString::Format("hqFrags_%s", Hadrons[iHad].c_str()).Data());

      TH1F* hIncl_all = (TH1F*) hInclJet->ProjectionX();
      TH1F* hg_all = (TH1F*) hgJet->ProjectionX();
      TH1F* hq_all = (TH1F*) hqJet->ProjectionX();

      hIncl_all->Scale(1./nJets);
      hg_all->Scale(1./nGluons_all);
      hq_all->Scale(1./nQuarks_all);
      plot_FFs_single(hIncl_all, hg_all, hq_all, Hadrons[iHad], jetPtBins.front(), jetPtBins.back());

      for (int ipt = 1; ipt < jetPtBins.size(); ipt++){
        double min_pt = jetPtBins[ipt - 1], max_pt = jetPtBins[ipt];
        hInclJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
        hgJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
        hqJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
        hNJetTypes->GetYaxis()->SetRangeUser(min_pt, max_pt);

        hIncl = (TH1F*) hInclJet->ProjectionX();
        hg = (TH1F*) hgJet->ProjectionX();
        hq = (TH1F*) hqJet->ProjectionX();
        hNgq = (TH1F*) hNJetTypes->ProjectionX();
        int nGluons_Pt = (int) hNgq->GetBinContent(1);
        int nQuarks_Pt = (int) hNgq->GetBinContent(2);
        int nJets_Pt = nGluons_Pt + nQuarks_Pt;

        hIncl_all->Scale(1./nJets_Pt);
        hg_all->Scale(1./nGluons_Pt);
        hq_all->Scale(1./nQuarks_Pt);
        // prep_single_hadron(hInclJet, hgJet, hqJet, hNJetTypes, hIncl, hg, hq, hNgq, min_pt, max_pt);
        // cout << "Prepped hadron" << endl;
        plot_FFs_single(hIncl, hg, hq, Hadrons[iHad], min_pt, max_pt);
      }
    }
  } // Single hadron plotting
  else cout << "Skipping single hadron fragmentation." << endl;
  if (plotMultHadron){
    TH2F* hgPiJet = (TH2F*) inFile->Get("hgFrags_pi");
    TH2F* hgKJet = (TH2F*) inFile->Get("hgFrags_K");
    TH2F* hgLambda0Jet = (TH2F*) inFile->Get("hgFrags_Lambda0");
    TH2F* hgpJet = (TH2F*) inFile->Get("hgFrags_p");

    TH1F* hgPi = (TH1F*) hgPiJet->ProjectionX();
    TH1F* hgK = (TH1F*) hgKJet->ProjectionX();
    TH1F* hgLambda0 = (TH1F*) hgLambda0Jet->ProjectionX();
    TH1F* hgp = (TH1F*) hgpJet->ProjectionX();

    hgPi->Scale(1./nGluons_all);
    hgK->Scale(1./nGluons_all);
    hgLambda0->Scale(1./nGluons_all);
    hgp->Scale(1./nGluons_all);

    plot_FFs_mult(hgPi, hgK, hgp, hgLambda0,
                  Hadrons[0], Hadrons[1], Hadrons[3], Hadrons[7],
                  jetPtBins.front(), jetPtBins.back(), "g");

    TH2F* hqPiJet = (TH2F*) inFile->Get("hqFrags_pi");
    TH2F* hqKJet = (TH2F*) inFile->Get("hqFrags_K");
    TH2F* hqLambda0Jet = (TH2F*) inFile->Get("hqFrags_Lambda0");
    TH2F* hqpJet = (TH2F*) inFile->Get("hqFrags_p");

    TH1F* hqPi = (TH1F*) hqPiJet->ProjectionX();
    TH1F* hqK = (TH1F*) hqKJet->ProjectionX();
    TH1F* hqLambda0 = (TH1F*) hqLambda0Jet->ProjectionX();
    TH1F* hqp = (TH1F*) hqpJet->ProjectionX();

    hqPi->Scale(1./nQuarks_all);
    hqK->Scale(1./nQuarks_all);
    hqLambda0->Scale(1./nQuarks_all);
    hqp->Scale(1./nQuarks_all);

    plot_FFs_mult(hqPi, hqK, hqp, hqLambda0,
                  Hadrons[0], Hadrons[1], Hadrons[3], Hadrons[7],
                  jetPtBins.front(), jetPtBins.back(), "q");

    for (int ipt = 1; ipt < jetPtBins.size(); ipt++){
      double min_pt = jetPtBins[ipt - 1], max_pt = jetPtBins[ipt];
      hNJetTypes->GetYaxis()->SetRangeUser(min_pt, max_pt);
      TH1F* hNgq = (TH1F*) hNJetTypes->ProjectionX();
      int nGluons_Pt = (int) hNgq->GetBinContent(1);
      int nQuarks_Pt = (int) hNgq->GetBinContent(2);
      int nJets_Pt = nGluons_Pt + nQuarks_Pt;

      hgPiJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
      hgKJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
      hgLambda0Jet->GetYaxis()->SetRangeUser(min_pt, max_pt);
      hgpJet->GetYaxis()->SetRangeUser(min_pt, max_pt);

      hgPi = (TH1F*) hgPiJet->ProjectionX();
      hgK = (TH1F*) hgKJet->ProjectionX();
      hgLambda0 = (TH1F*) hgLambda0Jet->ProjectionX();
      hgp = (TH1F*) hgpJet->ProjectionX();

      hgPi->Scale(1./nGluons_Pt);
      hgK->Scale(1./nGluons_Pt);
      hgLambda0->Scale(1./nGluons_Pt);
      hgp->Scale(1./nGluons_Pt);

      plot_FFs_mult(hgPi, hgK, hgp, hgLambda0,
                    Hadrons[0], Hadrons[1], Hadrons[3], Hadrons[7],
                    min_pt, max_pt, "g");

      hqPiJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
      hqKJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
      hqLambda0Jet->GetYaxis()->SetRangeUser(min_pt, max_pt);
      hqpJet->GetYaxis()->SetRangeUser(min_pt, max_pt);

      hqPi = (TH1F*) hqPiJet->ProjectionX();
      hqK = (TH1F*) hqKJet->ProjectionX();
      hqLambda0 = (TH1F*) hqLambda0Jet->ProjectionX();
      hqp = (TH1F*) hqpJet->ProjectionX();

      hqPi->Scale(1./nQuarks_Pt);
      hqK->Scale(1./nQuarks_Pt);
      hqLambda0->Scale(1./nQuarks_Pt);
      hqp->Scale(1./nQuarks_Pt);

      plot_FFs_mult(hqPi, hqK, hqp, hqLambda0,
                    Hadrons[0], Hadrons[1], Hadrons[3], Hadrons[7],
                    min_pt, max_pt, "q");
    }

    // hPiJet = (TH2F*) inFile->Get("hqFrags_pi");
    // hKJet = (TH2F*) inFile->Get("hqFrags_K");
    // hLambda0Jet = (TH2F*) inFile->Get("hqFrags_Lambda0");
    // hpJet = (TH2F*) inFile->Get("hqFrags_p");

    // hPi = (TH1F*) hPiJet->ProjectionX();
    // hK = (TH1F*) hKJet->ProjectionX();
    // hLambda0 = (TH1F*) hLambda0Jet->ProjectionX();
    // hp = (TH1F*) hpJet->ProjectionX();

    // plot_FFs_mult(hPi, hK, hp, hLambda0, Hadrons[0], Hadrons[1], Hadrons[3], Hadrons[7], jetPtBins.front(), jetPtBins.back());
  } // Multiple hadron plotting
  else cout << "Skipping multiple hadron fragmentation comparison." << endl;
  if (plotMatchDist){
    TH1F* hDeltaPartonJet = (TH1F*) inFile->Get("hDeltaPartonJet");
    std::vector<double> matchDist = {0.3, 0.5, 0.7};
    plot_matchDist(hDeltaPartonJet, matchDist);
  }

  time = (clock() - time)/CLOCKS_PER_SEC;
  cout << "Time taken: " << time << " seconds." << endl;
  // Do not close the file, it will delete the ratio histogram
  // inFile->Close();
}

string format_hadron_name(string hadron)
{
  string had = hadron;
  if (hadron == "pi"){
    had = "#pi^{#pm}";
  }
  else if (hadron == "pi0"){
    had = "#pi^{0}";
  }
  else if (hadron == "K0L"){
    had = "K^{0}_{L}";
  }
  else if (hadron == "K0S"){
    had = "K^{0}_{S}";
  }
  else if (hadron == "K0"){
    had = "K^{0}";
  }
  else if (hadron == "Lambda0"){
    had = "#Lambda^{0}";
  }
  return had;
}

void plot_matchDist(TH1F* dist, std::vector<double> matchDist)
{
  TCanvas *c_matchDist = new TCanvas("matchDist", "Distance between parton and jet", 1600, 1000);
  dist->Scale(1./dist->Integral());
  dist->SetLineColor(GetColor(2));
  dist->SetMarkerColor(GetColor(2));
  TH1F* distClone = (TH1F*) dist->Clone();
  double ymax = dist->GetBinContent(dist->GetMaximumBin());
  ymax *= 1.1;
  auto legend = CreateLegend(0.5, 0.95, 0.6, 0.9, "", 0.05);
  TH1F* frame = DrawFrame(0, 1., 0, ymax, "#Delta", "");

  distClone->GetXaxis()->SetRangeUser(0, matchDist[2]);
  double frac3 = distClone->Integral();
  TLine* L3 = new TLine(matchDist[2], 0, matchDist[2], distClone->GetBinContent(distClone->FindBin(matchDist[2])));
  L3->SetLineColor(GetColor(3));
  legend->AddEntry(L3, TString::Format("#Delta: %.1f, frac: %.1f", matchDist[2], frac3).Data());

  distClone->GetXaxis()->SetRangeUser(0, matchDist[1]);
  double frac2 = distClone->Integral();
  TLine* L2 = new TLine(matchDist[1], 0, matchDist[1], distClone->GetBinContent(distClone->FindBin(matchDist[1])));
  L2->SetLineColor(GetColor(1));
  legend->AddEntry(L2, TString::Format("#Delta: %.1f, frac: %.1f", matchDist[1], frac2).Data());

  distClone->GetXaxis()->SetRangeUser(0, matchDist[0]);
  double frac1 = distClone->Integral();
  TLine* L1 = new TLine(matchDist[0], 0, matchDist[0], distClone->GetBinContent(distClone->FindBin(matchDist[0])));
  L1->SetLineColor(GetColor(0));
  legend->AddEntry(L1, TString::Format("#Delta: %.1f, frac: %.1f", matchDist[0], frac1).Data());

  frame->Draw();
  dist->Draw("same hist");
  L1->Draw("same");
  L2->Draw("same");
  L3->Draw("same");
  legend->Draw("same");
  c_matchDist->SaveAs(TString::Format("../plots/matchDist.pdf"));
  return;
}

void prep_single_hadron(TH2F* hInclJet, TH2F* hgJet, TH2F* hqJet, TH2F* hNJetTypes,
                        TH1F* hIncl, TH1F* hg, TH1F* hq, TH1F* hNgq,
                        double min_pt, double max_pt)
{
  hInclJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
  hgJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
  hqJet->GetYaxis()->SetRangeUser(min_pt, max_pt);
  hNJetTypes->GetYaxis()->SetRangeUser(min_pt, max_pt);

  hIncl = (TH1F*) hInclJet->ProjectionX();
  hg = (TH1F*) hgJet->ProjectionX();
  hq = (TH1F*) hqJet->ProjectionX();
  hNgq = (TH1F*) hNJetTypes->ProjectionX();

  int nGluons_Pt = (int) hNgq->GetBinContent(1);
  int nQuarks_Pt = (int) hNgq->GetBinContent(2);
  int nJets_Pt = nGluons_Pt + nQuarks_Pt;

  hIncl->Scale(1./nJets_Pt);
  hg->Scale(1./nGluons_Pt);
  hq->Scale(1./nQuarks_Pt);
  return;
}

// Plot Fragmentation Function for single hadron species. Plot inclusive, gluon, and quark jets
void plot_FFs_single(TH1F* hIncl, TH1F* hg, TH1F* hq, string hadron, double min_pt, double max_pt)
{
  hIncl->SetStats(0);
  hIncl->SetLineColor(GetColor(0));
  hIncl->SetMarkerColor(GetColor(0));
  hIncl->SetMarkerStyle(GetMarker(0));

  hg->SetStats(1);
  hg->SetLineColor(GetColor(0));
  hg->SetMarkerColor(GetColor(1));
  hg->SetMarkerStyle(GetMarker(1));

  hq->SetStats(2);
  hq->SetLineColor(GetColor(0));
  hq->SetMarkerColor(GetColor(2));
  hq->SetMarkerStyle(GetMarker(2));

  auto legend = CreateLegend(0.75, 0.95, 0.6, 0.9, "", 0.05);
  legend->AddEntry(hIncl, "incl.");
  legend->AddEntry(hg, "g");
  legend->AddEntry(hq, "q");

  TCanvas *c_had = new TCanvas(TString::Format("c_%s_pt_%.0f_%.0f",
                                               hadron.c_str(), min_pt, max_pt).Data(),
                               TString::Format("%s_pt_%.0f_%.0f",
                                               hadron.c_str(), min_pt, max_pt).Data(),
                               1600, 1000);
  c_had->SetLogy();

  // cout << "Plot_FFs_single after setup" << endl;

  double xmin = 0., xmax = 1., ymin = 1e-5, ymax;
  ymax = max({hIncl->GetMaximum(), hg->GetMaximum(), hq->GetMaximum()});
  ymax *= 1.1;
  string s = format_hadron_name(hadron);
  string yTitle = TString::Format("#frac{1}{N_{jets}} #frac{dN_{%s}}{dz}", s.c_str()).Data();
  TH1F* frame1 = DrawFrame(xmin, xmax, ymin, ymax, "#it{z}", yTitle);

  // cout << "Plot_FFs_single constructed frame" << endl;

  frame1->Draw();
  hIncl->Draw("same");
  hg->Draw("same");
  hq->Draw("same");
  legend->Draw("same");
  // cout << "Plot_FFs_single after draw" << endl;
  c_had->SaveAs(TString::Format("../plots/SingleHadrons/%s_pt_%.0f_%.0f.pdf", hadron.c_str(), min_pt, max_pt).Data());
  return;
}

void plot_FFs_mult(TH1F* h0, TH1F* h1, TH1F* h2, TH1F* h3,
                   string hadron0, string hadron1, string hadron2, string hadron3,
                   double min_pt, double max_pt, string IGQ)
{
  auto legend = CreateLegend(0.75, 0.95, 0.6, 0.9, "", 0.05);

  h0->SetStats(0);
  h0->SetLineColor(GetColor(0));
  h0->SetMarkerColor(GetColor(0));
  h0->SetMarkerStyle(GetMarker(0));
  string had0 = format_hadron_name(hadron0);
  legend->AddEntry(h0, had0.c_str());

  h1->SetStats(1);
  h1->SetLineColor(GetColor(1));
  h1->SetMarkerColor(GetColor(1));
  h1->SetMarkerStyle(GetMarker(1));
  string had1 = format_hadron_name(hadron1);
  legend->AddEntry(h1, had1.c_str());

  h2->SetStats(2);
  h2->SetLineColor(GetColor(2));
  h2->SetMarkerColor(GetColor(2));
  h2->SetMarkerStyle(GetMarker(2));
  string had2 = format_hadron_name(hadron2);
  legend->AddEntry(h2, had2.c_str());

  h3->SetStats(3);
  h3->SetLineColor(GetColor(3));
  h3->SetMarkerColor(GetColor(3));
  h3->SetMarkerStyle(GetMarker(3));
  string had3 = format_hadron_name(hadron3);
  legend->AddEntry(h3, had3.c_str());

  TCanvas *c_mult = new TCanvas(TString::Format("c_mult_pt_%.0f_%.0f",
                                               min_pt, max_pt).Data(),
                               TString::Format("mult_pt_%.0f_%.0f",
                                               min_pt, max_pt).Data(),
                               1600, 1000);
  c_mult->SetLogy();

  double xmin = 0., xmax = 1., ymin = 1e-5, ymax;
  ymax = max({h0->GetMaximum(), h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum()});
  ymax *= 1.1;
  string yTitle = TString::Format("#frac{1}{N_{jets}^{%s}} #frac{dN_{h}}{dz}", IGQ.c_str()).Data();
  TH1F* frame1 = DrawFrame(xmin, xmax, ymin, ymax, "#it{z}", yTitle);

  // cout << "Plot_FFs_single constructed frame" << endl;

  frame1->Draw();
  h0->Draw("same");
  h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
  legend->Draw("same");
  // cout << "Plot_FFs_single after draw" << endl;
  string str;
  if (IGQ == "incl"){
    str = "inclusive";
  }
  else if (IGQ == "g"){
    str = "gluon";
  }
  else if (IGQ == "q"){
    str = "quark";
  }
  c_mult->SaveAs(TString::Format("../plots/MultipleHadrons/%s_%s_%s_%s_%s_pt_%.0f_%.0f.pdf",
                                str.c_str(), hadron0.c_str(), hadron1.c_str(), hadron2.c_str(), hadron3.c_str(),
                                min_pt, max_pt).Data());
  return;
}